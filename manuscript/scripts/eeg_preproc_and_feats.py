#!/usr/bin/env python3
"""
eeg_preproc_and_feats.py

Exact Python reconstruction of the MATLAB EEG preprocessing, FOOOF spectral
parameterization, and single-trial feature extraction pipeline from:
- Preproc_OrientWheel_Exo_ML.m
- VpixxEarlyTriggerFix.m
- gratton_emcp.m
- mark_Heyemove_ML.m & step_t.m
- reject_trials_inspect_Exo_ML.m
- rej_beh_trials_Exo_ML.m & Sep_byCond_OrientExo_ML.m
- fooof_OrientExo_ML.m & fooof_analysis.py
- Feat_Extract_OrientExo_ML.m & Conn_Feats_OrientExoML.m
- PrepData_DDTBOX.m & PrepConnData_DDTBOX.m
"""

import os
import re
import time
import warnings
from pathlib import Path
from typing import Dict, Any, List, Tuple

import numpy as np
import scipy.io as sio
import scipy.signal as sig
from scipy.special import lpmv
from joblib import Parallel, delayed
from numba import njit, prange

warnings.filterwarnings("ignore")

SCRIPTS_DIR = Path(__file__).resolve().parent
MANUSCRIPT_DIR = SCRIPTS_DIR.parent
PROJECT_ROOT = MANUSCRIPT_DIR.parent
EEG_DIR = PROJECT_ROOT / "data" / "eeg"
BEH_DIR = PROJECT_ROOT / "data" / "behavior"
CACHE_DIR = PROJECT_ROOT / "data" / "cache"
CED_FILE = PROJECT_ROOT / "repos" / "MVPA_CuedOrientTask_paper" / "Analysis_Code" / "EOG-electrode-locs-32_orientwheel.ced"

INCLUDED_PARTICIPANTS = [
    "004", "005", "007", "009", "010", "012", "013", "014", "015", "017",
    "018", "019", "021", "022", "024", "025", "026", "027", "029", "031",
    "032", "034", "035", "036", "037", "038", "039", "040",
]

BRAIN_ELEC_NAMES = [
    "Oz", "Pz", "Cz", "FCz", "Fz", "O1", "O2", "PO3", "PO4", "P7", "P8",
    "P5", "P6", "P3", "P4", "CP5", "CP6", "CP1", "CP2", "C3", "C4", "FC5",
    "FC6", "FC1", "FC2", "F7", "F8", "F3", "F4", "Fp1", "Fp2",
]

# Exact 1-based visual inspection trial rejections from reject_trials_inspect_Exo_ML.m (byTargets_ML_v1)
VISUAL_INSPECT_REJ_1BASED: Dict[str, List[int]] = {
    "004": [86, 164, 165, 167, 170, 179, 210, 218, 243, 306, 308, 310, 314, 315, 316, 318],
    "005": [242, 273],
    "007": [],
    "009": [311],
    "010": [1, 53, 129, 175, 176, 195, 198, 297],
    "012": [3, 64, 106, 153, 181, 238, 324],
    "013": [84, 267, 344],
    "014": [308, 334, 355],
    "015": [156],
    "017": [154, 211, 299] + list(range(301, 307)) + list(range(310, 339)),
    "018": [147, 183],
    "019": [250, 251],
    "021": [136],
    "022": [51, 114, 171, 199, 203, 214],
    "024": [],
    "025": [],
    "026": [],
    "027": [95, 127, 128, 242],
    "029": [2, 123, 142, 202, 228, 229, 263, 266, 267, 273, 277, 281, 299],
    "031": [223, 364, 368],
    "032": [],
    "034": [155, 183, 214, 281],
    "035": [278],
    "036": [48, 138, 199, 263, 292, 308, 344],
    "037": [46, 82],
    "038": [132, 287, 303, 321],
    "039": [],
    "040": [13, 77, 79, 93, 122, 125, 131, 142, 146, 151, 180, 199, 205, 206, 236, 249],
}


def load_electrode_xyz() -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[Dict[str, Any]]]:
    """Loads the 31 brain electrode coordinates from EOG-electrode-locs-32_orientwheel.ced."""
    lines = CED_FILE.read_text(encoding="utf-8").strip().splitlines()
    elecs = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        idx = int(parts[0].strip())
        label = parts[1].strip()
        theta = float(parts[2].strip())
        radius = float(parts[3].strip())
        x = float(parts[4].strip())
        y = float(parts[5].strip())
        z = float(parts[6].strip())
        elecs.append({
            "idx": idx, "label": label, "theta": theta, "radius": radius,
            "X": x, "Y": y, "Z": z,
        })
    # Brain electrodes are indices 2..32 (1-based) = 1..31 (0-based)
    brain = elecs[1:32]
    X = np.array([e["X"] for e in brain], dtype=np.float64)
    Y = np.array([e["Y"] for e in brain], dtype=np.float64)
    Z = np.array([e["Z"] for e in brain], dtype=np.float64)
    return X, Y, Z, brain


def compute_laplacian_matrices(X: np.ndarray, Y: np.ndarray, Z: np.ndarray, leg_order: int = 20, smoothing: float = 1e-5):
    """Precomputes the Perrin et al. (1989) surface Laplacian operator matrix L_op so surf_lap = L_op @ data."""
    numelectrodes = len(X)
    m = 4 if numelectrodes <= 100 else 3
    radii = np.sqrt(X**2 + Y**2 + Z**2)
    maxrad = np.max(radii)
    x = X / maxrad
    y = Y / maxrad
    z = Z / maxrad

    cosdist = np.zeros((numelectrodes, numelectrodes), dtype=np.float64)
    for i in range(numelectrodes):
        for j in range(i + 1, numelectrodes):
            cosdist[i, j] = 1.0 - (((x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2) / 2.0)
    cosdist = cosdist + cosdist.T + np.eye(numelectrodes)
    cosdist = np.clip(cosdist, -1.0, 1.0)

    legpoly = np.zeros((leg_order, numelectrodes, numelectrodes), dtype=np.float64)
    for ni in range(1, leg_order + 1):
        legpoly[ni - 1] = lpmv(0, ni, cosdist)

    n_vec = np.arange(1, leg_order + 1, dtype=np.float64)
    twoN1 = 2.0 * n_vec + 1.0
    gdenom = (n_vec * (n_vec + 1.0)) ** m
    hdenom = (n_vec * (n_vec + 1.0)) ** (m - 1)

    G = np.einsum("n,nij->ij", twoN1 / gdenom, legpoly) / (4.0 * np.pi)
    H = np.einsum("n,nij->ij", twoN1 / hdenom, legpoly) / (4.0 * np.pi)

    Gs = G + np.eye(numelectrodes) * smoothing
    Gs_inv = np.linalg.inv(Gs)
    GsinvS = np.sum(Gs_inv, axis=0, keepdims=True)  # (1, E)
    sum_GsinvS = np.sum(GsinvS)

    # For data shape (E, K):
    # dataGs = data.T @ Gs_inv  (K, E)
    # C = dataGs - (np.sum(dataGs, axis=1, keepdims=True) / sum_GsinvS) @ GsinvS
    # surf_lap = (C @ H.T).T = H @ C.T
    # Let's express C.T as M_C @ data:
    # dataGs.T = Gs_inv @ data
    # C.T = Gs_inv @ data - (GsinvS.T / sum_GsinvS) @ (np.ones((1, E)) @ Gs_inv @ data)
    ones_row = np.ones((1, numelectrodes), dtype=np.float64)
    M_C = Gs_inv - (GsinvS.T / sum_GsinvS) @ (ones_row @ Gs_inv)
    L_op = H @ M_C  # (E, E) linear operator!
    return L_op


def gratton_emcp_numpy(epochs: np.ndarray, trig_codes: np.ndarray, srate: float = 1000.0) -> np.ndarray:
    """
    Exact NumPy translation of gratton_emcp.m (Gratton, Coles & Donchin, 1983).
    epochs: (34, pnts, trials) in uV, where ch 32 is VEOG and ch 33 is HEOG (0-based).
    trig_codes: (trials,) array of target trigger codes (20 or 30).
    """
    nbchan, pnts, ntrials = epochs.shape
    temp_ocular = epochs[32:34, :, :].copy()

    # Prepend VEOG and HEOG as channels 0 and 1 of EEGdata (36, pnts, ntrials)
    EEGdata = np.concatenate([epochs[32:34, :, :], epochs], axis=0).astype(np.float64)

    blinkCriteria = 14.0
    msTen = int(round((10.0 / (1000.0 / srate)) + 0.5))
    thirdOfBlink = msTen * 7
    if thirdOfBlink % 2 == 0:
        thirdOfBlink += 1
    middleOfBlink = int(round((thirdOfBlink + 1) / 2.0))
    lengthOfWind = thirdOfBlink * 3.0
    initBlinkScan = thirdOfBlink + middleOfBlink  # 1-based
    endBlinkScan = pnts - initBlinkScan - thirdOfBlink + 1  # 1-based
    windVariance = 2.0

    # Fast cumulative sum scan for blink template covariance on VEOG (EEGdata[0])
    veog = EEGdata[0, :, :]  # (pnts, ntrials)
    cs = np.vstack([np.zeros((1, ntrials), dtype=np.float64), np.cumsum(veog, axis=0)])

    pts_1b = np.arange(initBlinkScan, endBlinkScan + 1)
    start1 = pts_1b - initBlinkScan + 1
    end1 = start1 + thirdOfBlink
    start2 = end1 + 1
    end2 = start1 + 2 * thirdOfBlink
    start3 = end2 + 1
    end3 = pts_1b + initBlinkScan - 1

    sum1 = cs[end1, :] - cs[start1 - 1, :]
    sum2 = cs[end2, :] - cs[start2 - 1, :]
    sum3 = cs[end3, :] - cs[start3 - 1, :]
    cov_scan = (-sum1 + 2.0 * sum2 - sum3) / lengthOfWind
    slope_scan = cov_scan / (windVariance * windVariance)

    markBlink = np.zeros((pnts, ntrials), dtype=bool)
    hits_pt_idx, hits_tr = np.where(np.abs(slope_scan) > blinkCriteria)
    for idx_pt, tr in zip(hits_pt_idx, hits_tr):
        pt_1b = initBlinkScan + idx_pt
        s_0b = max(0, pt_1b - middleOfBlink)
        e_0b = min(pnts, pt_1b + middleOfBlink - 1)
        markBlink[s_0b:e_0b, tr] = True

    # Subtract epoch mean from each waveform
    EEGdata_centered = EEGdata - np.mean(EEGdata, axis=1, keepdims=True)
    subdata = EEGdata_centered.copy()

    # Selection cards: bin 30 and bin 20
    for code in (30, 20):
        sel = np.where(trig_codes == code)[0]
        if len(sel) > 0:
            avg_bin = np.mean(subdata[:, :, sel], axis=2, keepdims=True)
            subdata[:, :, sel] -= avg_bin

    rblinkpts = int(np.sum(markBlink))
    if rblinkpts == 0:
        markBlink[:, :] = False

    PropBlink = np.zeros((nbchan + 2, 2), dtype=np.float64)
    PropSaccade = np.zeros((nbchan + 2, 2), dtype=np.float64)

    if rblinkpts > 0:
        b_vals = subdata[:, markBlink]  # (36, rblinkpts)
        b_vals = b_vals - np.mean(b_vals, axis=1, keepdims=True)
        xvecB = b_vals[0:2, :].T  # (rblinkpts, 2)
        XTXI_B = np.linalg.pinv(xvecB.T @ xvecB)
        PropBlink = (XTXI_B @ (xvecB.T @ b_vals.T)).T  # (36, 2)

    s_vals = subdata[:, ~markBlink]  # (36, rsaccpts)
    s_vals = s_vals - np.mean(s_vals, axis=1, keepdims=True)
    xvecS = s_vals[0:2, :].T
    XTXI_S = np.linalg.pinv(xvecS.T @ xvecS)
    PropSaccade = (XTXI_S @ (xvecS.T @ s_vals.T)).T  # (36, 2)

    # Reconstitute EEGdata with only epoch mean subtracted
    subdata_raw = EEGdata_centered
    out_epochs = epochs.copy()

    if rblinkpts > 0:
        EEGblinkmean = np.mean(subdata_raw[:, markBlink], axis=1)  # (36,)
        deltaBlink = EEGdata_centered[0:2, markBlink] - EEGblinkmean[0:2, None]  # (2, rblinkpts)
        adjB = PropBlink[2:, :] @ deltaBlink + EEGblinkmean[2:, None]  # (34, rblinkpts)
        out_epochs[:, markBlink] -= adjB

    EEGsaccademean = np.mean(subdata_raw[:, ~markBlink], axis=1)  # (36,)
    deltaSaccade = EEGdata_centered[0:2, ~markBlink] - EEGsaccademean[0:2, None]  # (2, rsaccpts)
    adjS = PropSaccade[2:, :] @ deltaSaccade + EEGsaccademean[2:, None]  # (34, rsaccpts)
    out_epochs[:, ~markBlink] -= adjS

    # Restore eye channels
    out_epochs[32:34, :, :] = temp_ocular
    return out_epochs


def parse_vmrk_events(vmrk_path: Path) -> List[Tuple[int, int]]:
    """
    Parses BrainVision .vmrk file and returns list of (code_int, latency_0based_shifted)
    after applying VpixxEarlyTriggerFix (+6 samples to all stimulus events after the first).
    """
    events = []
    text = vmrk_path.read_text(encoding="utf-8", errors="ignore")
    mk_idx = 0
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("Mk") and "=" in line:
            mk_idx += 1
            rhs = line.split("=", 1)[1]
            parts = rhs.split(",")
            if len(parts) >= 3 and parts[0] == "Stimulus":
                m = re.search(r"(\d+)", parts[1])
                if m:
                    code = int(m.group(1))
                    pos_1b = int(parts[2].strip())
                    # VpixxEarlyTriggerFix adds +6 samples to events 2..end
                    pos_0b = (pos_1b - 1) + 6
                    events.append((code, pos_0b))
    return events


def butter_bandpass_twopass(data: np.ndarray, lowcut: float, highcut: float, fs: float = 1000.0, order: int = 4) -> np.ndarray:
    """
    Exact Python equivalent of FieldTrip ft_preproc_bandpassfilter(..., filt_order, 'but', 'twopass', 'no'):
    4th-order Butterworth bandpass filter applied forward-backward (filtfilt).
    """
    nyq = 0.5 * fs
    low = max(1e-4, lowcut / nyq)
    high = min(0.999, highcut / nyq)
    sos = sig.butter(order, [low, high], btype="bandpass", output="sos")
    return sig.sosfiltfilt(sos, data, axis=-1)


def extract_cfc_features_fast(
    cond_data: np.ndarray,
    iaf_per_chan: np.ndarray,
    fs: float = 1000.0,
) -> Dict[str, np.ndarray]:
    """
    Vectorized extraction of the 6 channel-level features (lowFreq_amp, AAC, PACoz, PACtort, PACglm, PACplv)
    from cond_data of shape (31, 3800, n_trials) where times_out is -2400..+1399 ms.
    Pre-stimulus window [-1000, 0] ms corresponds to indices 1400..2400 (1001 samples) in cond_data,
    or 2000..3000 in 600-sample zero-padded data (length 5000).
    Returns DDTBOX-binned arrays:
    - lowFreq_amp: (50, 20, 31, n_trials)
    - AAC, PACoz, PACtort, PACglm, PACplv: (50, 15, 31, n_trials)
    """
    n_chan, n_pnts, n_trials = cond_data.shape
    pad = 600
    # Transpose to (31, n_trials, n_pnts) and zero-pad 600 samples on both ends along last axis
    padded = np.pad(np.transpose(cond_data, (0, 2, 1)), ((0, 0), (0, 0), (pad, pad)), mode="constant")
    # [-1000, 0] ms in padded array: original -2400 ms is index 600, so -1000 ms is index 600 + 1400 = 2000,
    # and 0 ms is index 600 + 2400 = 3000 (inclusive -> slice 2000:3001, 1001 samples).
    s_pre, e_pre = 2000, 3001
    Ntime = 1001
    winSize = 250
    n_steps = 750  # (1001 - 250) / 1 = 751 in 1:(Ntime-winSize), which has length 750 in MATLAB (1:750)!

    # High-frequency gamma [15, 60] Hz across all channels & trials at once
    filt_high = butter_bandpass_twopass(padded, 15.0, 60.0, fs=fs, order=4)
    analytic_high = sig.hilbert(filt_high, axis=-1)
    amp_high_all = np.abs(analytic_high[:, :, s_pre:e_pre])  # (31, n_trials, 1001)

    # Low-frequency filter per channel around iaf_per_chan[ch] +/- 1 Hz
    phase_low_all = np.empty((n_chan, n_trials, Ntime), dtype=np.float64)
    amp_low_all = np.empty((n_chan, n_trials, Ntime), dtype=np.float64)
    phase_ampH_all = np.empty((n_chan, n_trials, Ntime), dtype=np.float64)

    for ch in range(n_chan):
        lf = float(iaf_per_chan[ch])
        if not np.isfinite(lf) or lf < 8.0 or lf > 14.0:
            lf = 10.5
        filt_low = butter_bandpass_twopass(padded[ch], lf - 1.0, lf + 1.0, fs=fs, order=4)
        an_low = sig.hilbert(filt_low, axis=-1)[:, s_pre:e_pre]
        phase_low_all[ch] = np.angle(an_low)
        amp_low_all[ch] = np.abs(an_low)
        # For cohen_PLV: phase of Hilbert-transformed high-frequency amplitude envelope
        an_ampH = sig.hilbert(amp_high_all[ch], axis=-1)
        phase_ampH_all[ch] = np.angle(an_ampH)

    # 1. lowFreq_amp: first 1000 samples of amp_low_all -> (31, n_trials, 1000) -> (50, 20, 31, n_trials)
    lf_1000 = amp_low_all[:, :, :1000].reshape(n_chan, n_trials, 50, 20)
    feat_lowFreq_amp = np.transpose(lf_1000, (2, 3, 0, 1)).astype(np.float32)

    # Helper for cumulative window sums of length W=250 over steps 0..749 along last axis (length 1001)
    def win_sum(arr: np.ndarray) -> np.ndarray:
        cs = np.concatenate([np.zeros(arr.shape[:-1] + (1,), dtype=arr.dtype), np.cumsum(arr, axis=-1)], axis=-1)
        return cs[..., winSize : winSize + n_steps] - cs[..., :n_steps]

    # 2. PACoz (Ozkurt MVL): |sum(A_H * exp(i * phi_L))| / sqrt(W * sum(A_H^2))
    z_oz = amp_high_all * np.exp(1j * phase_low_all)
    num_oz = np.abs(win_sum(z_oz))
    den_oz = np.sqrt(winSize * np.maximum(win_sum(amp_high_all ** 2), 1e-15))
    pac_oz = (num_oz / den_oz).reshape(n_chan, n_trials, 50, 15)
    feat_PACoz = np.transpose(pac_oz, (2, 3, 0, 1)).astype(np.float32)

    # 3. PACplv (Cohen PLV): |mean(exp(i * (phi_L - phi_ampH)))|
    z_plv = np.exp(1j * (phase_low_all - phase_ampH_all))
    pac_plv = (np.abs(win_sum(z_plv)) / winSize).reshape(n_chan, n_trials, 50, 15)
    feat_PACplv = np.transpose(pac_plv, (2, 3, 0, 1)).astype(np.float32)

    # 4. PACtort (Tort KL Modulation Index with 18 phase bins)
    nbins = 18
    phase_2pi = np.mod(phase_low_all, 2.0 * np.pi)
    bin_idx = np.clip((phase_2pi / (2.0 * np.pi / nbins)).astype(np.int32), 0, nbins - 1)
    bin_means = np.empty((nbins, n_chan, n_trials, n_steps), dtype=np.float64)
    for b in range(nbins):
        mask_b = (bin_idx == b).astype(np.float64)
        cnt_b = win_sum(mask_b)
        sum_b = win_sum(amp_high_all * mask_b)
        bin_means[b] = np.where(cnt_b > 0, sum_b / np.maximum(cnt_b, 1.0), 1e-12)
    P_bins = bin_means / np.maximum(np.sum(bin_means, axis=0, keepdims=True), 1e-15)
    P_bins = np.maximum(P_bins, 1e-12)
    kl_mi = np.sum(P_bins * np.log(P_bins * nbins), axis=0) / np.log(nbins)
    pac_tort = kl_mi.reshape(n_chan, n_trials, 50, 15)
    feat_PACtort = np.transpose(pac_tort, (2, 3, 0, 1)).astype(np.float32)

    # 5. PACglm (GLM R^2 of Amp_H ~ cos(phi_L) + sin(phi_L) + 1)
    c = np.cos(phase_low_all)
    s = np.sin(phase_low_all)
    y = amp_high_all
    Sc = win_sum(c)
    Ss = win_sum(s)
    Scc = win_sum(c * c)
    Sss = win_sum(s * s)
    Scs = win_sum(c * s)
    Sy = win_sum(y)
    Syc = win_sum(y * c)
    Sys = win_sum(y * s)
    Syy = win_sum(y * y)

    # Center variables within each 250-sample window to eliminate intercept analytically
    yc = Syc - (Sy * Sc) / winSize
    ys = Sys - (Sy * Ss) / winSize
    cc = Scc - (Sc * Sc) / winSize + 1e-12
    ss = Sss - (Ss * Ss) / winSize + 1e-12
    cs_cov = Scs - (Sc * Ss) / winSize
    sst = np.maximum(Syy - (Sy * Sy) / winSize, 1e-12)

    det = np.maximum(cc * ss - cs_cov * cs_cov, 1e-12)
    b1 = (ss * yc - cs_cov * ys) / det
    b2 = (cc * ys - cs_cov * yc) / det
    ssr_explained = np.clip(b1 * yc + b2 * ys, 0.0, sst)
    pac_glm = (ssr_explained / sst).reshape(n_chan, n_trials, 50, 15)
    feat_PACglm = np.transpose(pac_glm, (2, 3, 0, 1)).astype(np.float32)

    # 6. AAC (Spearman rank correlation between amp_low and amp_high in each 250-ms window)
    aac_all = compute_aac_numba(amp_low_all, amp_high_all, winSize, n_steps)
    feat_AAC = np.transpose(aac_all.reshape(n_chan, n_trials, 50, 15), (2, 3, 0, 1))

    return {
        "lowFreq_amp": feat_lowFreq_amp,
        "AAC": feat_AAC,
        "PACoz": feat_PACoz,
        "PACtort": feat_PACtort,
        "PACglm": feat_PACglm,
        "PACplv": feat_PACplv,
    }


@njit(parallel=True, fastmath=True)
def compute_aac_numba(amp_low: np.ndarray, amp_high: np.ndarray, winSize: int = 250, n_steps: int = 750) -> np.ndarray:
    n_chan, n_trials, _ = amp_low.shape
    out = np.empty((n_chan, n_trials, n_steps), dtype=np.float32)
    denom = float(winSize * (winSize * winSize - 1))
    for ch in prange(n_chan):
        for tr in range(n_trials):
            xL = amp_low[ch, tr]
            xH = amp_high[ch, tr]
            rL = np.empty(winSize, dtype=np.int32)
            rH = np.empty(winSize, dtype=np.int32)
            for i in range(winSize):
                vL = xL[i]
                vH = xH[i]
                cL = 0
                cH = 0
                for j in range(winSize):
                    if xL[j] < vL or (xL[j] == vL and j < i):
                        cL += 1
                    if xH[j] < vH or (xH[j] == vH and j < i):
                        cH += 1
                rL[i] = cL
                rH[i] = cH
            d2 = 0
            for i in range(winSize):
                diff = rL[i] - rH[i]
                d2 += diff * diff
            out[ch, tr, 0] = 1.0 - (6.0 * d2) / denom

            for t in range(1, n_steps):
                newL = xL[t + winSize - 1]
                newH = xH[t + winSize - 1]
                rk_oldL = rL[0]
                rk_oldH = rH[0]
                rk_newL = 0
                rk_newH = 0
                d2 = 0
                for k in range(winSize - 1):
                    curL = rL[k + 1]
                    if curL > rk_oldL:
                        curL -= 1
                    vLk = xL[t + k]
                    if vLk <= newL:
                        rk_newL += 1
                    else:
                        curL += 1
                    rL[k] = curL

                    curH = rH[k + 1]
                    if curH > rk_oldH:
                        curH -= 1
                    vHk = xH[t + k]
                    if vHk <= newH:
                        rk_newH += 1
                    else:
                        curH += 1
                    rH[k] = curH
                    diff = curL - curH
                    d2 += diff * diff
                rL[winSize - 1] = rk_newL
                rH[winSize - 1] = rk_newH
                diff = rk_newL - rk_newH
                d2 += diff * diff
                out[ch, tr, t] = 1.0 - (6.0 * d2) / denom
    return out


@njit(parallel=True, fastmath=True)
def compute_fc_numba(an_raw: np.ndarray, an_lap: np.ndarray, winSize: int = 375, winStep: int = 5, n_steps: int = 120):
    n_chan, n_trials, _ = an_raw.shape
    plv_out = np.empty((n_steps, 465, n_trials), dtype=np.float32)
    corr_out = np.empty((n_steps, 465, n_trials), dtype=np.float32)
    icoh_out = np.empty((n_steps, 465, n_trials), dtype=np.float32)
    wpli_out = np.empty((n_steps, 465, n_trials), dtype=np.float32)

    u_lap = an_lap / np.maximum(np.abs(an_lap), 1e-15)
    amp_lap = np.abs(an_lap)

    norm_rank_denom = np.sqrt(winSize * (winSize * winSize - 1) / 12.0)
    mean_rank = (winSize - 1) / 2.0

    for ij in prange(n_steps):
        st = ij * winStep
        en = st + winSize
        r_norm = np.empty((n_chan, n_trials, winSize), dtype=np.float64)
        spec_diag = np.empty((n_chan, n_trials), dtype=np.float64)
        for ch in range(n_chan):
            for tr in range(n_trials):
                w = amp_lap[ch, tr, st:en]
                order = np.argsort(w)
                for rk in range(winSize):
                    r_norm[ch, tr, order[rk]] = (rk - mean_rank) / norm_rank_denom
                s_pow = 0.0
                for t in range(st, en):
                    z = an_raw[ch, tr, t]
                    s_pow += z.real * z.real + z.imag * z.imag
                spec_diag[ch, tr] = s_pow

        pair_idx = 0
        for c1 in range(n_chan - 1):
            for c2 in range(c1 + 1, n_chan):
                for tr in range(n_trials):
                    plv_re = 0.0
                    plv_im = 0.0
                    corr_sum = 0.0
                    cross_im = 0.0
                    wpli_imagsum = 0.0
                    wpli_imagsumW = 0.0
                    wpli_debias = 0.0
                    for k in range(winSize):
                        t = st + k
                        u1 = u_lap[c1, tr, t]
                        u2 = u_lap[c2, tr, t]
                        plv_re += u1.real * u2.real + u1.imag * u2.imag
                        plv_im += u1.imag * u2.real - u1.real * u2.imag
                        corr_sum += r_norm[c1, tr, k] * r_norm[c2, tr, k]
                        z1 = an_raw[c1, tr, t]
                        z2 = an_raw[c2, tr, t]
                        cdi = z1.imag * z2.real - z1.real * z2.imag
                        cross_im += cdi
                        wpli_imagsum += cdi
                        wpli_imagsumW += abs(cdi)
                        wpli_debias += cdi * cdi

                    plv_out[ij, pair_idx, tr] = np.sqrt(plv_re * plv_re + plv_im * plv_im) / winSize
                    corr_out[ij, pair_idx, tr] = corr_sum
                    denom_ic = np.sqrt(max(spec_diag[c1, tr] * spec_diag[c2, tr], 1e-30))
                    icoh_out[ij, pair_idx, tr] = abs(cross_im) / denom_ic
                    denom_dwpli = wpli_imagsumW * wpli_imagsumW - wpli_debias
                    wpli_out[ij, pair_idx, tr] = (
                        (wpli_imagsum * wpli_imagsum - wpli_debias) / denom_dwpli
                        if denom_dwpli > 1e-30
                        else 0.0
                    )
                pair_idx += 1

    return plv_out, corr_out, icoh_out, wpli_out


def extract_fc_features_fast(
    cond_data: np.ndarray,
    L_op: np.ndarray,
    fs: float = 1000.0,
) -> Dict[str, np.ndarray]:
    """
    Vectorized extraction of the 4 Functional Connectivity features (wPLI, corr, PLV, iCOH)
    across all 465 upper-triangular electrode pairs and 120 time steps (12 windows x 10 steps)
    in [-1000, 0] ms (winSize = 375 ms, winStep = 5 ms).
    Returns DDTBOX-binned arrays of shape (12, 10, 465, n_trials).
    """
    n_chan, n_pnts, n_trials = cond_data.shape
    lap_data = (L_op @ cond_data.reshape(n_chan, -1)).reshape(n_chan, n_pnts, n_trials)

    pad = 600
    padded_raw = np.pad(np.transpose(cond_data, (0, 2, 1)), ((0, 0), (0, 0), (pad, pad)), mode="constant")
    padded_lap = np.pad(np.transpose(lap_data, (0, 2, 1)), ((0, 0), (0, 0), (pad, pad)), mode="constant")

    s_pre, e_pre = 2000, 3001
    filt_raw = butter_bandpass_twopass(padded_raw, 8.0, 14.0, fs=fs, order=4)
    filt_lap = butter_bandpass_twopass(padded_lap, 8.0, 14.0, fs=fs, order=4)

    an_raw = sig.hilbert(filt_raw, axis=-1)[:, :, s_pre:e_pre]
    an_lap = sig.hilbert(filt_lap, axis=-1)[:, :, s_pre:e_pre]

    plv_out, corr_out, icoh_out, wpli_out = compute_fc_numba(an_raw, an_lap, 375, 5, 120)

    return {
        "PLV": plv_out.reshape(12, 10, 465, n_trials),
        "corr": corr_out.reshape(12, 10, 465, n_trials),
        "iCOH": icoh_out.reshape(12, 10, 465, n_trials),
        "wPLI": wpli_out.reshape(12, 10, 465, n_trials),
    }


@njit(parallel=True, fastmath=True)
def fooof_batch_numba(freqs: np.ndarray, psds: np.ndarray):
    n_chan, n_trials, n_all = psds.shape
    idx_list = []
    for k in range(n_all):
        if 2.0 <= freqs[k] <= 50.0:
            idx_list.append(k)
    K = len(idx_list)
    f_sub = np.empty(K, dtype=np.float64)
    log_f = np.empty(K, dtype=np.float64)
    for k in range(K):
        f_sub[k] = freqs[idx_list[k]]
        log_f[k] = np.log10(f_sub[k])

    iaf_out = np.full((n_chan, n_trials), np.nan, dtype=np.float64)
    ap_out = np.full((n_chan, n_trials, 2), np.nan, dtype=np.float64)
    r2_out = np.full((n_chan, n_trials), np.nan, dtype=np.float64)
    err_out = np.full((n_chan, n_trials), np.nan, dtype=np.float64)

    for ch in prange(n_chan):
        for tr in range(n_trials):
            L = np.empty(K, dtype=np.float64)
            for k in range(K):
                L[k] = np.log10(max(psds[ch, tr, idx_list[k]], 1e-30))

            slope_init = (L[K - 1] - L[0]) / (log_f[K - 1] - log_f[0])
            off_init = L[0] - slope_init * log_f[0]
            flat_init = np.empty(K, dtype=np.float64)
            for k in range(K):
                flat_init[k] = L[k] - (off_init + slope_init * log_f[k])
            med_flat = np.median(flat_init)
            sx = 0.0; sy = 0.0; sxx = 0.0; sxy = 0.0; cnt = 0
            for k in range(K):
                if flat_init[k] <= med_flat:
                    x = log_f[k]
                    y = L[k]
                    sx += x; sy += y; sxx += x * x; sxy += x * y; cnt += 1
            if cnt < 2:
                for k in range(K):
                    x = log_f[k]; y = L[k]
                    sx += x; sy += y; sxx += x * x; sxy += x * y; cnt += 1
            denom = cnt * sxx - sx * sx
            slope = (cnt * sxy - sx * sy) / denom
            offset = (sy - slope * sx) / cnt

            F = np.empty(K, dtype=np.float64)
            for k in range(K):
                F[k] = L[k] - (offset + slope * log_f[k])
            G_sum = np.zeros(K, dtype=np.float64)

            best_alpha_cf = np.nan
            best_alpha_h = -1.0
            for _ in range(7):
                std_f = np.std(F)
                thresh = max(0.22, 0.2 * std_f)
                max_idx = 0
                max_val = F[0]
                for k in range(1, K):
                    if F[k] > max_val:
                        max_val = F[k]
                        max_idx = k
                if max_val < thresh:
                    break
                cf = f_sub[max_idx]
                h = max_val
                if 0 < max_idx < K - 1:
                    y0 = F[max_idx - 1]
                    y1 = F[max_idx]
                    y2 = F[max_idx + 1]
                    d = 2.0 * (2.0 * y1 - y0 - y2)
                    if abs(d) > 1e-12:
                        delta = (y2 - y0) / d
                        if -1.0 < delta < 1.0:
                            df = f_sub[1] - f_sub[0]
                            cf = f_sub[max_idx] + delta * df
                half_h = 0.5 * h
                left_f = f_sub[0]
                for k in range(max_idx, -1, -1):
                    if F[k] <= half_h:
                        left_f = f_sub[k]
                        break
                right_f = f_sub[K - 1]
                for k in range(max_idx, K):
                    if F[k] <= half_h:
                        right_f = f_sub[k]
                        break
                fwhm = max(min(right_f - left_f, 5.0), 1.0)
                sd = fwhm / 2.0
                for k in range(K):
                    g = h * np.exp(-((f_sub[k] - cf) ** 2) / (2.0 * sd * sd))
                    F[k] -= g
                    G_sum[k] += g
                if 8.0 <= cf <= 14.0 and h > best_alpha_h:
                    best_alpha_h = h
                    best_alpha_cf = cf

            sx = 0.0; sy = 0.0; sxx = 0.0; sxy = 0.0
            for k in range(K):
                x = log_f[k]
                y = L[k] - G_sum[k]
                sx += x; sy += y; sxx += x * x; sxy += x * y
            denom = K * sxx - sx * sx
            slope_final = (K * sxy - sx * sy) / denom
            offset_final = (sy - slope_final * sx) / K
            exponent_final = -slope_final

            L_hat = np.empty(K, dtype=np.float64)
            mae = 0.0
            mean_L = 0.0
            mean_hat = 0.0
            for k in range(K):
                lh = offset_final + slope_final * log_f[k] + G_sum[k]
                L_hat[k] = lh
                mae += abs(L[k] - lh)
                mean_L += L[k]
                mean_hat += lh
            mae /= K
            mean_L /= K
            mean_hat /= K
            cov = 0.0; v1 = 0.0; v2 = 0.0
            for k in range(K):
                d1 = L[k] - mean_L
                d2 = L_hat[k] - mean_hat
                cov += d1 * d2
                v1 += d1 * d1
                v2 += d2 * d2
            r = cov / max(np.sqrt(v1 * v2), 1e-15)
            iaf_out[ch, tr] = best_alpha_cf
            ap_out[ch, tr, 0] = offset_final
            ap_out[ch, tr, 1] = exponent_final
            r2_out[ch, tr] = r * r
            err_out[ch, tr] = mae
    return iaf_out, ap_out, r2_out, err_out


def run_fooof_on_condition(cond_data: np.ndarray, fs: float = 1000.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Runs FOOOF spectral parameterization on pre-stimulus [-1000, 0] ms Welch PSDs (nperseg=500)
    across all 31 channels and trials for one condition.
    Returns (iaf_per_chan, aperiodic_medians_31x2, r2_medians_31, err_medians_31).
    """
    pre_data = np.transpose(cond_data[:, 1400:2401, :], (0, 2, 1))  # (31, n_trials, 1001)
    freqs, psds = sig.welch(pre_data, fs=fs, window="hamming", nperseg=500, noverlap=250, nfft=512, axis=-1)
    iaf_tr, ap_tr, r2_tr, err_tr = fooof_batch_numba(freqs, psds)
    iaf_per_chan = np.nanmedian(iaf_tr, axis=1)
    ap_med = np.nanmedian(ap_tr, axis=1)
    r2_med = np.nanmedian(r2_tr, axis=1)
    err_med = np.nanmedian(err_tr, axis=1)
    return iaf_per_chan, ap_med, r2_med, err_med


def process_single_subject(subj_id: str) -> Dict[str, Any]:
    """
    Runs the complete preprocessing, FOOOF, and 10-feature extraction for a single subject
    and caches the result in `data/cache/{subj_id}_feats.npz`.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / f"{subj_id}_feats.npz"
    if cache_file.exists():
        return {"subj_id": subj_id, "status": "cached", "cache_file": str(cache_file)}

    eeg_file = EEG_DIR / f"{subj_id}_Orient_Exo_Exp.eeg"
    vmrk_file = EEG_DIR / f"{subj_id}_Orient_Exo_Exp.vmrk"
    beh_file = BEH_DIR / f"{subj_id}_Orient_Exo.mat"
    if not eeg_file.exists() or not vmrk_file.exists() or not beh_file.exists():
        return {"subj_id": subj_id, "status": "missing_raw_file"}

    t0 = time.time()
    # 1. Load raw multiplexed INT_16 EEG (34 channels, 0.1 uV resolution)
    raw = np.fromfile(eeg_file, dtype="<i2")
    n_samples = raw.size // 34
    data = raw[: n_samples * 34].reshape(n_samples, 34).T.astype(np.float64) * 0.1

    # 2. Linked-mastoid re-reference for brain channels 1..31 (Ch0 is M2)
    data[1:32, :] -= 0.5 * data[0:1, :]

    # 3. Bandpass filter 0.1 - 50 Hz
    sos_bp = sig.butter(4, [0.1 / 500.0, 50.0 / 500.0], btype="bandpass", output="sos")
    data = sig.sosfiltfilt(sos_bp, data, axis=-1)

    # 4. Parse markers + VpixxEarlyTriggerFix (+6 samples) + epoch [-2400, +1400] ms around target triggers (20, 30)
    events = parse_vmrk_events(vmrk_file)
    target_events = [(code, pos) for (code, pos) in events if code in (20, 30)]

    epoch_list = []
    trig_list = []
    for code, pos in target_events:
        if pos - 2400 >= 0 and pos + 1400 <= n_samples:
            epoch_list.append(data[:, pos - 2400 : pos + 1400])
            trig_list.append(code)
    epochs = np.stack(epoch_list, axis=-1)  # (34, 3800, n_total)
    trigs = np.array(trig_list, dtype=np.int32)

    # 5. Remove first 20 practice trials (+ 009 truncation)
    epochs = epochs[:, :, 20:]
    trigs = trigs[20:]

    mat = sio.loadmat(str(beh_file), squeeze_me=True, struct_as_record=False)
    bdata = mat["data"]
    err_deg = np.asarray(bdata.errorDegrees, dtype=np.float64)[20:]
    valid = np.asarray(bdata.valid, dtype=np.int32)[20:]
    position = np.asarray(bdata.position, dtype=np.int32)[20:]

    if subj_id == "009":
        epochs = epochs[:, :, :343]
        trigs = trigs[:343]
        err_deg = err_deg[:343]
        valid = valid[:343]
        position = position[:343]

    # Ensure length matches minimum in case of trailing trigger
    n_min = min(epochs.shape[2], len(err_deg))
    epochs = epochs[:, :, :n_min]
    trigs = trigs[:n_min]
    err_deg = err_deg[:n_min]
    valid = valid[:n_min]
    position = position[:n_min]

    # 6. Pre-EMCP baseline subtraction (pop_rmbase(EEG, exp.epochbaseline) with exp.epochbaseline = [] in byTargets_ML_v1_Settings.mat) + Artifact Rejection 1: range [-1000, 1000] uV on channels 1..34
    epochs -= np.mean(epochs, axis=1, keepdims=True)
    bad1 = np.any((epochs < -1000.0) | (epochs > 1000.0), axis=(0, 1))
    keep1 = ~bad1
    epochs, trigs, err_deg, valid, position = (
        epochs[:, :, keep1], trigs[keep1], err_deg[keep1], valid[keep1], position[keep1]
    )

    # 7. Gratton EMCP Ocular Correction
    epochs = gratton_emcp_numpy(epochs, trigs, srate=1000.0)

    # 8. Post-EMCP baseline subtraction (pop_rmbase(EEG, exp.epochbaseline)) + Artifact Rejection 2: range [-500, 500] uV on channels 1..32 (0..31)
    epochs -= np.mean(epochs, axis=1, keepdims=True)
    bad2 = np.any((epochs[:32, :, :] < -500.0) | (epochs[:32, :, :] > 500.0), axis=(0, 1))
    keep2 = ~bad2
    epochs, trigs, err_deg, valid, position = (
        epochs[:, :, keep2], trigs[keep2], err_deg[keep2], valid[keep2], position[keep2]
    )

    # 9. Horizontal Eye Movement Rejection (mark_Heyemove_ML.m + step_t.m)
    # Window [-1600, 600] ms -> indices 800..3001 on HEOG (channel index 33)
    heog_win = epochs[33, 800:3001, :]  # (2201, n_trials)
    n_heog = heog_win.shape[0]
    bad3 = np.zeros(epochs.shape[2], dtype=bool)
    for j in range(0, n_heog - 100 + 1, 10):
        w1 = np.mean(heog_win[j : j + 50, :], axis=0)
        w2 = np.mean(heog_win[j + 50 : j + 100, :], axis=0)
        bad3 |= np.abs(w2 - w1) > 25.0
    keep3 = ~bad3
    epochs, trigs, err_deg, valid, position = (
        epochs[:, :, keep3], trigs[keep3], err_deg[keep3], valid[keep3], position[keep3]
    )

    # 10. Visual inspection trial rejection (reject_trials_inspect_Exo_ML.m)
    vis_rej_1b = VISUAL_INSPECT_REJ_1BASED.get(subj_id, [])
    bad4 = np.zeros(epochs.shape[2], dtype=bool)
    for idx_1b in vis_rej_1b:
        if 1 <= idx_1b <= len(bad4):
            bad4[idx_1b - 1] = True
    keep4 = ~bad4
    epochs, trigs, err_deg, valid, position = (
        epochs[:, :, keep4], trigs[keep4], err_deg[keep4], valid[keep4], position[keep4]
    )

    # Keep only the 31 brain electrodes (indices 1..31, 0-based)
    brain_epochs = epochs[1:32, :, :]  # (31, 3800, n_clean_trials)

    # Separate into 4 conditions:
    # Liv: position==0 & valid==0
    # Lv:  position==0 & valid==1
    # Riv: position==1 & valid==0
    # Rv:  position==1 & valid==1
    cond_masks = {
        "Liv": (position == 0) & (valid == 0),
        "Lv":  (position == 0) & (valid == 1),
        "Riv": (position == 1) & (valid == 0),
        "Rv":  (position == 1) & (valid == 1),
    }

    X, Y, Z, _ = load_electrode_xyz()
    L_op = compute_laplacian_matrices(X, Y, Z)

    save_dict: Dict[str, Any] = {
        "subj_id": subj_id,
        "err_deg": err_deg,
        "valid": valid,
        "position": position,
    }

    for cname, mask in cond_masks.items():
        c_epochs = brain_epochs[:, :, mask]
        c_svr_labels = np.abs(err_deg[mask]).astype(np.float64)
        save_dict[f"SVR_labels_{cname}"] = c_svr_labels

        # FOOOF spectral parameterization
        iaf_c, ap_c, r2_c, err_c = run_fooof_on_condition(c_epochs, fs=1000.0)
        save_dict[f"fooof_iaf_{cname}"] = iaf_c
        save_dict[f"fooof_ap_{cname}"] = ap_c
        save_dict[f"fooof_r2_{cname}"] = r2_c
        save_dict[f"fooof_err_{cname}"] = err_c

        # CFC & Alpha Amplitude features
        cfc_feats = extract_cfc_features_fast(c_epochs, iaf_c, fs=1000.0)
        for mname, arr in cfc_feats.items():
            save_dict[f"feat_{mname}_{cname}"] = arr

        # Functional Connectivity features
        fc_feats = extract_fc_features_fast(c_epochs, L_op, fs=1000.0)
        for mname, arr in fc_feats.items():
            save_dict[f"feat_{mname}_{cname}"] = arr

    np.savez(cache_file, **save_dict)
    dt = time.time() - t0
    print(f"[eeg_preproc_and_feats] Subject {subj_id}: {brain_epochs.shape[2]} clean trials, done in {dt:.1f}s", flush=True)
    return {"subj_id": subj_id, "status": "computed", "clean_trials": int(brain_epochs.shape[2]), "seconds": round(dt, 1)}


def main():
    avail = [
        s for s in INCLUDED_PARTICIPANTS
        if (EEG_DIR / f"{s}_Orient_Exo_Exp.eeg").exists()
        and (EEG_DIR / f"{s}_Orient_Exo_Exp.vmrk").exists()
    ]
    print(f"[eeg_preproc_and_feats] Processing {len(avail)} available subjects with Numba multi-core acceleration...", flush=True)
    results = []
    for idx, s in enumerate(avail, 1):
        r = process_single_subject(s)
        results.append(r)
        print(f"  [{idx}/{len(avail)}] {r}", flush=True)


if __name__ == "__main__":
    main()
