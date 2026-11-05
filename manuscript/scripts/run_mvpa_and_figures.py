#!/usr/bin/env python3
"""
run_mvpa_and_figures.py

Complete, zero-hotwiring Python reconstruction of:
1. FOOOF group-level F_max permutation statistics (fooof_OrientExo_ML.m)
2. Spatiotemporal DDTBOX LIBSVM Support Vector Classification (SVM, C=1.0) and
   Support Vector Regression (SVR, C=0.1, eps=0.1) across all 10 pre-stimulus
   oscillatory metrics (Alpha Amplitude, AAC, MVL Ozkurt, MI Tort, PAC GLM,
   PAC PLV, wPLI, Amplitude Correlation, PLV, iCOH) with 10-fold x 10-repeat
   cross-validation, shuffled-label empirical chance baselines, Haufe et al. (2014)
   activation pattern feature weights, and 10,000-permutation Phipson & Smyth (2010)
   t_max multiple-comparison correction (analyse_decoding_erp_mine.m,
   tmax_classifier_accuracies.m, mxt_perm1_Phipson_Smyth.m, FWanal_steps_*.m).
3. Publication-ready rendering of Figures 3, 4, 5, 6, and 7 (300 DPI PNG + vector PDF)
   adhering to Reviewer 2's formatting requirements and clean scientific style.
4. Automatic registration of all statistics and figures into `analysis_results.json`
   and `generated_stats.tex` via `pipeline_api.py`.
"""

import os
import sys
import time
import json
import warnings
from pathlib import Path
from typing import Dict, Any, List, Tuple

import numpy as np
import scipy.interpolate as interp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, Polygon
from numba import njit, prange

warnings.filterwarnings("ignore")

SCRIPTS_DIR = Path(__file__).resolve().parent
MANUSCRIPT_DIR = SCRIPTS_DIR.parent
PROJECT_ROOT = MANUSCRIPT_DIR.parent
CACHE_DIR = PROJECT_ROOT / "data" / "cache"
FIGURES_DIR = MANUSCRIPT_DIR / "figures"

sys.path.insert(0, str(SCRIPTS_DIR))
from eeg_preproc_and_feats import INCLUDED_PARTICIPANTS
from pipeline_api import OrientExoPipelineAPI


# ==============================================================================
# 1. REDBLUE COLORMAP (exact match to redblue.m in OrientTask_paper)
# ==============================================================================
def get_redblue_cmap(m: int = 256) -> LinearSegmentedColormap:
    """Exact match to redblue.m: blue (0,0,1) -> white (1,1,1) -> red (1,0,0)."""
    return LinearSegmentedColormap.from_list(
        "redblue_matlab",
        [(0.0, 0.0, 1.0), (1.0, 1.0, 1.0), (1.0, 0.0, 0.0)],
        N=m,
    )


# ==============================================================================
# 2. NUMBA-ACCELERATED LIBSVM SMO SOLVERS & DDTBOX DECODING + HAUFE WEIGHTS
# ==============================================================================
@njit(fastmath=True)
def solve_csvc_smo(K: np.ndarray, y: np.ndarray, C: float = 1.0, tol: float = 1e-2, max_iter: int = 600):
    """Fan et al. (2005) LIBSVM dual SMO solver for C-SVC (-s 0 -t 0 -c 1)."""
    n = len(y)
    alpha = np.zeros(n, dtype=np.float64)
    G = np.full(n, -1.0, dtype=np.float64)
    for _ in range(max_iter):
        i = -1
        G_max = -1e30
        G_min = 1e30
        for t in range(n):
            if (y[t] == 1 and alpha[t] < C) or (y[t] == -1 and alpha[t] > 0.0):
                if -y[t] * G[t] >= G_max:
                    i = t
                    G_max = -y[t] * G[t]
        j = -1
        obj_min = 1e30
        for t in range(n):
            if (y[t] == 1 and alpha[t] > 0.0) or (y[t] == -1 and alpha[t] < C):
                b = G_max + y[t] * G[t]
                if -y[t] * G[t] <= G_min:
                    G_min = -y[t] * G[t]
                if b > 0.0 and i != -1:
                    quad = K[i, i] + K[t, t] - 2.0 * K[i, t]
                    if quad <= 0.0:
                        quad = 1e-12
                    val = -(b * b) / quad
                    if val <= obj_min:
                        j = t
                        obj_min = val
        if G_max - G_min < tol or i == -1 or j == -1:
            break
        quad = K[i, i] + K[j, j] - 2.0 * K[i, j]
        if quad <= 0.0:
            quad = 1e-12
        old_ai = alpha[i]
        old_aj = alpha[j]
        if y[i] != y[j]:
            delta = (-G[i] - G[j]) / quad
            diff = old_ai - old_aj
            alpha[i] += delta
            alpha[j] += delta
            if diff > 0.0:
                if alpha[j] < 0.0:
                    alpha[j] = 0.0
                    alpha[i] = diff
                if alpha[i] > C:
                    alpha[i] = C
                    alpha[j] = C - diff
            else:
                if alpha[i] < 0.0:
                    alpha[i] = 0.0
                    alpha[j] = -diff
                if alpha[j] > C:
                    alpha[j] = C
                    alpha[i] = C + diff
        else:
            delta = (G[i] - G[j]) / quad
            sm = old_ai + old_aj
            alpha[i] -= delta
            alpha[j] += delta
            if sm > C:
                if alpha[i] > C:
                    alpha[i] = C
                    alpha[j] = sm - C
                if alpha[j] > C:
                    alpha[j] = C
                    alpha[i] = sm - C
            else:
                if alpha[j] < 0.0:
                    alpha[j] = 0.0
                    alpha[i] = sm
                if alpha[i] < 0.0:
                    alpha[i] = 0.0
                    alpha[j] = sm
        d_ai = alpha[i] - old_ai
        d_aj = alpha[j] - old_aj
        yi = y[i]
        yj = y[j]
        for t in range(n):
            G[t] += y[t] * (d_ai * yi * K[t, i] + d_aj * yj * K[t, j])

    nr_free = 0
    ub = 1e30
    lb = -1e30
    sum_free = 0.0
    for t in range(n):
        yG = y[t] * G[t]
        if alpha[t] == C:
            if y[t] == -1:
                ub = min(ub, yG)
            else:
                lb = max(lb, yG)
        elif alpha[t] == 0.0:
            if y[t] == 1:
                ub = min(ub, yG)
            else:
                lb = max(lb, yG)
        else:
            nr_free += 1
            sum_free += yG
    rho = sum_free / nr_free if nr_free > 0 else 0.5 * (ub + lb)
    return alpha * y, -rho


@njit(fastmath=True)
def solve_esvr_smo(K: np.ndarray, y_target: np.ndarray, C: float = 0.1, eps: float = 0.1, tol: float = 1e-2, max_iter: int = 600):
    """Fan et al. (2005) LIBSVM dual SMO solver for epsilon-SVR (-s 3 -t 0 -c 0.1 -p 0.1)."""
    n = len(y_target)
    n2 = 2 * n
    alpha = np.zeros(n2, dtype=np.float64)
    sign = np.empty(n2, dtype=np.int32)
    G = np.empty(n2, dtype=np.float64)
    for t in range(n):
        sign[t] = 1
        sign[t + n] = -1
        G[t] = eps - y_target[t]
        G[t + n] = eps + y_target[t]

    for _ in range(max_iter):
        i = -1
        G_max = -1e30
        G_min = 1e30
        for t in range(n2):
            if (sign[t] == 1 and alpha[t] < C) or (sign[t] == -1 and alpha[t] > 0.0):
                if -sign[t] * G[t] >= G_max:
                    i = t
                    G_max = -sign[t] * G[t]
        j = -1
        obj_min = 1e30
        ii = i % n if i >= 0 else 0
        for t in range(n2):
            if (sign[t] == 1 and alpha[t] > 0.0) or (sign[t] == -1 and alpha[t] < C):
                b = G_max + sign[t] * G[t]
                if -sign[t] * G[t] <= G_min:
                    G_min = -sign[t] * G[t]
                if b > 0.0 and i != -1:
                    tt = t % n
                    quad = K[ii, ii] + K[tt, tt] - 2.0 * K[ii, tt]
                    if quad <= 0.0:
                        quad = 1e-12
                    val = -(b * b) / quad
                    if val <= obj_min:
                        j = t
                        obj_min = val
        if G_max - G_min < tol or i == -1 or j == -1:
            break
        jj = j % n
        quad = K[ii, ii] + K[jj, jj] - 2.0 * K[ii, jj]
        if quad <= 0.0:
            quad = 1e-12
        old_ai = alpha[i]
        old_aj = alpha[j]
        si = sign[i]
        sj = sign[j]
        if si != sj:
            delta = (-G[i] - G[j]) / quad
            diff = old_ai - old_aj
            alpha[i] += delta
            alpha[j] += delta
            if diff > 0.0:
                if alpha[j] < 0.0:
                    alpha[j] = 0.0
                    alpha[i] = diff
                if alpha[i] > C:
                    alpha[i] = C
                    alpha[j] = C - diff
            else:
                if alpha[i] < 0.0:
                    alpha[i] = 0.0
                    alpha[j] = -diff
                if alpha[j] > C:
                    alpha[j] = C
                    alpha[i] = C + diff
        else:
            delta = (G[i] - G[j]) / quad
            sm = old_ai + old_aj
            alpha[i] -= delta
            alpha[j] += delta
            if sm > C:
                if alpha[i] > C:
                    alpha[i] = C
                    alpha[j] = sm - C
                if alpha[j] > C:
                    alpha[j] = C
                    alpha[i] = sm - C
            else:
                if alpha[j] < 0.0:
                    alpha[j] = 0.0
                    alpha[i] = sm
                if alpha[i] < 0.0:
                    alpha[i] = 0.0
                    alpha[j] = sm
        d_ai = alpha[i] - old_ai
        d_aj = alpha[j] - old_aj
        for t in range(n):
            kt = d_ai * si * K[t, ii] + d_aj * sj * K[t, jj]
            G[t] += kt
            G[t + n] -= kt

    nr_free = 0
    ub = 1e30
    lb = -1e30
    sum_free = 0.0
    for t in range(n2):
        yG = sign[t] * G[t]
        if alpha[t] == C:
            if sign[t] == -1:
                ub = min(ub, yG)
            else:
                lb = max(lb, yG)
        elif alpha[t] == 0.0:
            if sign[t] == 1:
                ub = min(ub, yG)
            else:
                lb = max(lb, yG)
        else:
            nr_free += 1
            sum_free += yG
    rho = sum_free / nr_free if nr_free > 0 else 0.5 * (ub + lb)
    beta = np.empty(n, dtype=np.float64)
    for t in range(n):
        beta[t] = alpha[t] - alpha[t + n]
    return beta, -rho


@njit(parallel=True, fastmath=True)
def ddtbox_svm_subject(X_c1: np.ndarray, X_c2: np.ndarray, seed: int = 42):
    """
    Runs 10-fold x 10-repeat spatiotemporal SVM classification (real + permuted labels)
    and Haufe-corrected absolute feature weight extraction for one subject.
    X_c1, X_c2: shape (n_steps, win_len, n_chan, n_trials)
    """
    n_steps, win_len, n_chan, n1 = X_c1.shape
    n2 = X_c2.shape[3]
    n_bal = min(n1, n2)
    fold_sz = max(1, n_bal // 10)
    n_use = fold_sz * 10
    n_train_c = n_use - fold_sz
    n_train = 2 * n_train_c
    n_test = 2 * fold_sz
    P = win_len * n_chan

    real_acc = np.zeros(n_steps, dtype=np.float64)
    perm_acc = np.zeros((n_steps, 10, 10), dtype=np.float64)
    fw_z = np.zeros((n_steps, n_chan), dtype=np.float64)

    for s in prange(n_steps):
        rng_state = np.uint64(seed + s * 10007)
        F_all = np.empty((n1 + n2, P), dtype=np.float64)
        for tr in range(n1):
            idx = 0
            for ch in range(n_chan):
                for w in range(win_len):
                    F_all[tr, idx] = X_c1[s, w, ch, tr]
                    idx += 1
        for tr in range(n2):
            idx = 0
            for ch in range(n_chan):
                for w in range(win_len):
                    F_all[n1 + tr, idx] = X_c2[s, w, ch, tr]
                    idx += 1

        for p in range(P):
            mn = np.min(F_all[:, p])
            mx = np.max(F_all[:, p])
            if mx > mn + 1e-12:
                F_all[:, p] = (F_all[:, p] - mn) / (mx - mn)
        K_all = F_all @ F_all.T

        sum_abs_fw = np.zeros(P, dtype=np.float64)
        acc_sum = 0.0

        for rep in range(10):
            p1 = np.arange(n1)
            p2 = np.arange(n2)
            for i in range(n1 - 1, 0, -1):
                rng_state = (rng_state * np.uint64(1664525) + np.uint64(1013904223)) & np.uint64(0xFFFFFFFF)
                j = int(rng_state % np.uint64(i + 1))
                p1[i], p1[j] = p1[j], p1[i]
            for i in range(n2 - 1, 0, -1):
                rng_state = (rng_state * np.uint64(1664525) + np.uint64(1013904223)) & np.uint64(0xFFFFFFFF)
                j = int(rng_state % np.uint64(i + 1))
                p2[i], p2[j] = p2[j], p2[i]
            u1 = p1[:n_use]
            u2 = p2[:n_use] + n1

            for cv in range(10):
                te_start = cv * fold_sz
                te_end = te_start + fold_sz
                tr_idx = np.empty(n_train, dtype=np.int32)
                te_idx = np.empty(n_test, dtype=np.int32)
                ti = 0
                vi = 0
                for k in range(n_use):
                    if te_start <= k < te_end:
                        te_idx[vi] = u1[k]
                        te_idx[vi + fold_sz] = u2[k]
                        vi += 1
                    else:
                        tr_idx[ti] = u1[k]
                        tr_idx[ti + n_train_c] = u2[k]
                        ti += 1

                K_tr = np.empty((n_train, n_train), dtype=np.float64)
                for i in range(n_train):
                    ii = tr_idx[i]
                    for j in range(n_train):
                        K_tr[i, j] = K_all[ii, tr_idx[j]]

                y_tr = np.empty(n_train, dtype=np.int32)
                y_tr[:n_train_c] = 1
                y_tr[n_train_c:] = -1
                y_te = np.empty(n_test, dtype=np.int32)
                y_te[:fold_sz] = 1
                y_te[fold_sz:] = -1

                coef, b = solve_csvc_smo(K_tr, y_tr, C=1.0)
                correct = 0
                for i in range(n_test):
                    ii = te_idx[i]
                    dec = b
                    for j in range(n_train):
                        dec += coef[j] * K_all[ii, tr_idx[j]]
                    pred = 1 if dec >= 0.0 else -1
                    if pred == y_te[i]:
                        correct += 1
                acc_sum += 100.0 * correct / n_test

                # Haufe activation pattern: a = cov(X_tr) @ w = X_tr^T @ ((K_tr @ coef - mean(K_tr @ coef)) / (n_train - 1))
                Kw = K_tr @ coef
                mean_Kw = np.mean(Kw)
                Kw_centered = (Kw - mean_Kw) / max(1, n_train - 1)
                for p in range(P):
                    ap = 0.0
                    for i in range(n_train):
                        ap += F_all[tr_idx[i], p] * Kw_centered[i]
                    sum_abs_fw[p] += abs(ap)

                # Permuted labels null
                y_perm = y_tr.copy()
                for i in range(n_train - 1, 0, -1):
                    rng_state = (rng_state * np.uint64(1664525) + np.uint64(1013904223)) & np.uint64(0xFFFFFFFF)
                    j = int(rng_state % np.uint64(i + 1))
                    y_perm[i], y_perm[j] = y_perm[j], y_perm[i]
                coef_p, b_p = solve_csvc_smo(K_tr, y_perm, C=1.0)
                correct_p = 0
                for i in range(n_test):
                    ii = te_idx[i]
                    dec = b_p
                    for j in range(n_train):
                        dec += coef_p[j] * K_all[ii, tr_idx[j]]
                    pred = 1 if dec >= 0.0 else -1
                    if pred == y_te[i]:
                        correct_p += 1
                perm_acc[s, cv, rep] = 100.0 * correct_p / n_test

        real_acc[s] = acc_sum / 100.0
        ch_abs = np.empty(n_chan, dtype=np.float64)
        for ch in range(n_chan):
            ch_abs[ch] = np.mean(sum_abs_fw[ch * win_len : (ch + 1) * win_len]) / 20.0
        m_ch = np.mean(ch_abs)
        sd_ch = np.std(ch_abs)
        if sd_ch > 1e-12:
            fw_z[s] = (ch_abs - m_ch) / sd_ch

    return real_acc, perm_acc, fw_z


@njit(parallel=True, fastmath=True)
def ddtbox_svr_subject(X_c: np.ndarray, y_raw: np.ndarray, seed: int = 1337):
    """
    Runs 10-fold x 10-repeat spatiotemporal epsilon-SVR regression (C=0.1, eps=0.1,
    real + permuted labels, Fisher-z transformed Pearson r) and Haufe-corrected
    absolute feature weight extraction for one subject and condition.
    X_c: shape (n_steps, win_len, n_chan, n_trials)
    y_raw: shape (n_trials,)
    """
    n_steps, win_len, n_chan, n_trials = X_c.shape
    fold_sz = max(2, n_trials // 10)
    n_use = fold_sz * 10
    n_train = n_use - fold_sz
    n_test = fold_sz
    P = win_len * n_chan

    # Normalize SVR target labels to [0, 1] (or z-score)
    y_min = np.min(y_raw)
    y_max = np.max(y_raw)
    if y_max > y_min + 1e-12:
        y_norm = (y_raw - y_min) / (y_max - y_min)
    else:
        y_norm = y_raw.copy()

    real_rz = np.zeros(n_steps, dtype=np.float64)
    perm_rz = np.zeros((n_steps, 10, 10), dtype=np.float64)
    fw_z = np.zeros((n_steps, n_chan), dtype=np.float64)

    for s in prange(n_steps):
        rng_state = np.uint64(seed + s * 10007)
        F_all = np.empty((n_trials, P), dtype=np.float64)
        for tr in range(n_trials):
            idx = 0
            for ch in range(n_chan):
                for w in range(win_len):
                    F_all[tr, idx] = X_c[s, w, ch, tr]
                    idx += 1

        for p in range(P):
            mn = np.min(F_all[:, p])
            mx = np.max(F_all[:, p])
            if mx > mn + 1e-12:
                F_all[:, p] = (F_all[:, p] - mn) / (mx - mn)
        K_all = F_all @ F_all.T

        sum_abs_fw = np.zeros(P, dtype=np.float64)
        rz_sum = 0.0

        for rep in range(10):
            p_ord = np.arange(n_trials)
            for i in range(n_trials - 1, 0, -1):
                rng_state = (rng_state * np.uint64(1664525) + np.uint64(1013904223)) & np.uint64(0xFFFFFFFF)
                j = int(rng_state % np.uint64(i + 1))
                p_ord[i], p_ord[j] = p_ord[j], p_ord[i]
            u_ord = p_ord[:n_use]

            for cv in range(10):
                te_start = cv * fold_sz
                te_end = te_start + fold_sz
                tr_idx = np.empty(n_train, dtype=np.int32)
                te_idx = np.empty(n_test, dtype=np.int32)
                ti = 0
                vi = 0
                for k in range(n_use):
                    if te_start <= k < te_end:
                        te_idx[vi] = u_ord[k]
                        vi += 1
                    else:
                        tr_idx[ti] = u_ord[k]
                        ti += 1

                K_tr = np.empty((n_train, n_train), dtype=np.float64)
                y_tr = np.empty(n_train, dtype=np.float64)
                for i in range(n_train):
                    ii = tr_idx[i]
                    y_tr[i] = y_norm[ii]
                    for j in range(n_train):
                        K_tr[i, j] = K_all[ii, tr_idx[j]]

                y_te = np.empty(n_test, dtype=np.float64)
                for i in range(n_test):
                    y_te[i] = y_norm[te_idx[i]]

                beta, b = solve_esvr_smo(K_tr, y_tr, C=0.1, eps=0.1)
                preds = np.empty(n_test, dtype=np.float64)
                for i in range(n_test):
                    ii = te_idx[i]
                    dec = b
                    for j in range(n_train):
                        dec += beta[j] * K_all[ii, tr_idx[j]]
                    preds[i] = dec

                # Fisher-z Pearson correlation between preds and y_te
                mp = np.mean(preds)
                my = np.mean(y_te)
                cov = 0.0
                vp = 0.0
                vy = 0.0
                for i in range(n_test):
                    dp = preds[i] - mp
                    dy = y_te[i] - my
                    cov += dp * dy
                    vp += dp * dp
                    vy += dy * dy
                denom = np.sqrt(vp * vy)
                r = cov / denom if denom > 1e-12 else 0.0
                r = min(max(r, -0.9999), 0.9999)
                rz_sum += 0.5 * np.log((1.0 + r) / (1.0 - r))

                # Haufe feature weights
                Kw = K_tr @ beta
                mean_Kw = np.mean(Kw)
                Kw_centered = (Kw - mean_Kw) / max(1, n_train - 1)
                for p in range(P):
                    ap = 0.0
                    for i in range(n_train):
                        ap += F_all[tr_idx[i], p] * Kw_centered[i]
                    sum_abs_fw[p] += abs(ap)

                # Permuted labels null
                y_perm = y_tr.copy()
                for i in range(n_train - 1, 0, -1):
                    rng_state = (rng_state * np.uint64(1664525) + np.uint64(1013904223)) & np.uint64(0xFFFFFFFF)
                    j = int(rng_state % np.uint64(i + 1))
                    y_perm[i], y_perm[j] = y_perm[j], y_perm[i]
                beta_p, b_p = solve_esvr_smo(K_tr, y_perm, C=0.1, eps=0.1)
                preds_p = np.empty(n_test, dtype=np.float64)
                for i in range(n_test):
                    ii = te_idx[i]
                    dec = b_p
                    for j in range(n_train):
                        dec += beta_p[j] * K_all[ii, tr_idx[j]]
                    preds_p[i] = dec
                mp_p = np.mean(preds_p)
                cov_p = 0.0
                vp_p = 0.0
                for i in range(n_test):
                    dp = preds_p[i] - mp_p
                    dy = y_te[i] - my
                    cov_p += dp * dy
                    vp_p += dp * dp
                denom_p = np.sqrt(vp_p * vy)
                rp = cov_p / denom_p if denom_p > 1e-12 else 0.0
                rp = min(max(rp, -0.9999), 0.9999)
                perm_rz[s, cv, rep] = 0.5 * np.log((1.0 + rp) / (1.0 - rp))

        real_rz[s] = rz_sum / 100.0
        ch_abs = np.empty(n_chan, dtype=np.float64)
        for ch in range(n_chan):
            ch_abs[ch] = np.mean(sum_abs_fw[ch * win_len : (ch + 1) * win_len]) / 20.0
        m_ch = np.mean(ch_abs)
        sd_ch = np.std(ch_abs)
        if sd_ch > 1e-12:
            fw_z[s] = (ch_abs - m_ch) / sd_ch

    return real_rz, perm_rz, fw_z


# ==============================================================================
# 3. GROUP PERMUTATION TESTS: mxt_perm1_Phipson_Smyth & F_max ANOVA
# ==============================================================================
def mxt_perm1_phipson_smyth(
    diff_data: np.ndarray,
    n_perm: int = 10000,
    alph: float = 0.05,
    tail: int = 1,
    seed: int = 2024,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Exact Python implementation of mxt_perm1_Phipson_Smyth.m.
    diff_data: shape (n_subs, n_pts) where each column is a time step or channel.
    Returns (corrected_p, corrected_h, t_orig, tmx_ptile).
    """
    n_subs, n_pts = diff_data.shape
    sqrt_nXnM1 = np.sqrt(n_subs * (n_subs - 1))

    # Observed t-scores
    sm_obs = np.sum(diff_data, axis=0)
    mn_obs = sm_obs / n_subs
    sm_sqrs_obs = np.sum(diff_data ** 2, axis=0) - (sm_obs ** 2) / n_subs
    stder_obs = np.maximum(np.sqrt(np.maximum(sm_sqrs_obs, 0.0)) / sqrt_nXnM1, 1e-12)
    t_orig = mn_obs / stder_obs

    rng = np.random.default_rng(seed)
    signs = rng.choice(np.array([-1.0, 1.0]), size=(n_perm, n_subs))
    sm_perm = signs @ diff_data  # (n_perm, n_pts)
    mn_perm = sm_perm / n_subs
    # Notice sum((d * s)^2) == sum(d^2) for s in {-1, +1}!
    sum_d2 = np.sum(diff_data ** 2, axis=0)[None, :]
    sm_sqrs_perm = np.maximum(sum_d2 - (sm_perm ** 2) / n_subs, 0.0)
    stder_perm = np.maximum(np.sqrt(sm_sqrs_perm) / sqrt_nXnM1, 1e-12)
    t_perm = mn_perm / stder_perm

    mxt_pos = np.max(np.abs(t_perm), axis=1)
    mxt = np.concatenate([mxt_pos, -mxt_pos])

    if tail == 1:
        tmx_ptile = float(np.percentile(mxt, 100.0 * (1.0 - alph)))
        b = np.sum(mxt[:, None] >= t_orig[None, :], axis=0)
        corrected_p = np.minimum((b + 1.0) / (n_perm + 1.0), 1.0)
    elif tail == 0:
        tmx_ptile = float(np.percentile(mxt, 100.0 * (1.0 - alph / 2.0)))
        b = np.sum(mxt[:, None] >= np.abs(t_orig)[None, :], axis=0)
        corrected_p = np.minimum(2.0 * (b + 1.0) / (n_perm + 1.0), 1.0)
    else:
        tmx_ptile = float(np.percentile(mxt, 100.0 * alph))
        b = np.sum(mxt[:, None] >= -t_orig[None, :], axis=0)
        corrected_p = np.minimum((b + 1.0) / (n_perm + 1.0), 1.0)

    corrected_h = (corrected_p < alph).astype(np.int32)
    return corrected_p, corrected_h, t_orig, tmx_ptile


def calc_fmax_2x2(
    data_4d: np.ndarray,
    effect: str,
    n_perm: int = 10000,
    alpha: float = 0.05 / 3.0,
    seed: int = 2024,
) -> Dict[str, Any]:
    """
    Exact Python implementation of calc_Fmax_ModelFit.m + perm_rbANOVA_OrientExo.m
    for a 2x2 repeated-measures design across 31 electrodes:
    data_4d: shape (31, 2_side, 2_cue, n_subs)
    effect: 'side' (dim 1), 'cue' (dim 2), or 'interaction' (dims [1, 2]).
    By reduce_data_ModelFit.m, each effect in a 2x2 design reduces to a 2-level
    one-way repeated-measures ANOVA across 31 electrodes!
    """
    n_elec, _, _, n_subs = data_4d.shape
    if effect == "side":
        reduced = np.mean(data_4d, axis=2)  # (31, 2, n_subs)
    elif effect == "cue":
        reduced = np.mean(data_4d, axis=1)  # (31, 2, n_subs)
    elif effect == "interaction":
        # Subtract across side (dim 1): shape (31, 2_cue, n_subs)
        reduced = data_4d[:, 1, :, :] - data_4d[:, 0, :, :]
    else:
        raise ValueError(f"Unknown effect: {effect}")

    # For a 2-condition repeated-measures one-way ANOVA, F = t^2 on the difference d = reduced[:, 1, :] - reduced[:, 0, :]!
    diff = (reduced[:, 1, :] - reduced[:, 0, :]).T  # (n_subs, 31)
    sqrt_nXnM1 = np.sqrt(n_subs * (n_subs - 1))

    sm_obs = np.sum(diff, axis=0)
    mn_obs = sm_obs / n_subs
    sm_sqrs_obs = np.maximum(np.sum(diff ** 2, axis=0) - (sm_obs ** 2) / n_subs, 1e-15)
    t_obs = mn_obs / (np.sqrt(sm_sqrs_obs) / sqrt_nXnM1)
    F_obs = t_obs ** 2  # (31,)

    rng = np.random.default_rng(seed)
    signs = rng.choice(np.array([-1.0, 1.0]), size=(n_perm, n_subs))
    signs[0, :] = 1.0  # first permutation is observed data (matching perm_rbANOVA_OrientExo.m line 90)
    sm_perm = signs @ diff
    mn_perm = sm_perm / n_subs
    sum_d2 = np.sum(diff ** 2, axis=0)[None, :]
    sm_sqrs_perm = np.maximum(sum_d2 - (sm_perm ** 2) / n_subs, 1e-15)
    t_perm = mn_perm / (np.sqrt(sm_sqrs_perm) / sqrt_nXnM1)
    F_dist = t_perm ** 2  # (n_perm, 31)

    Fmax_dist = np.sort(np.max(F_dist, axis=1))
    crit_idx = min(int(np.ceil((1.0 - alpha) * n_perm)) - 1, n_perm - 1)
    Fmax_crit = float(Fmax_dist[crit_idx])
    p_vals = np.mean(Fmax_dist[:, None] >= F_obs[None, :], axis=0)
    h_vals = (F_obs > Fmax_crit).astype(np.int32)

    max_elec_idx = int(np.argmax(F_obs))
    return {
        "F_obs": F_obs,
        "p_vals": p_vals,
        "h_vals": h_vals,
        "Fmax_crit": Fmax_crit,
        "max_F": float(F_obs[max_elec_idx]),
        "min_p": float(p_vals[max_elec_idx]),
        "n_sig": int(np.sum(h_vals)),
        "df": [1, n_subs - 1],
    }


def find_contiguous_runs(h_vec: np.ndarray) -> List[np.ndarray]:
    """Returns list of 0-based index arrays for each contiguous run of 1s in h_vec."""
    runs = []
    in_run = False
    start = 0
    for i, val in enumerate(h_vec):
        if val and not in_run:
            in_run = True
            start = i
        elif not val and in_run:
            in_run = False
            runs.append(np.arange(start, i))
    if in_run:
        runs.append(np.arange(start, len(h_vec)))
    return runs


# ==============================================================================
# 4. TOPOPLOT & CONNECTOME PLOTTING HELPERS (Figures 3, 4, 5, 6)
# ==============================================================================
def get_2d_topo_coords() -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Computes 2D polar-to-cartesian scalp projection (x_plot, y_plot) matching
    EEGLAB topoplot / topoplot_connect.m for the 31 brain electrodes.
    """
    ced_path = PROJECT_ROOT / "repos" / "MVPA_CuedOrientTask_paper" / "Analysis_Code" / "EOG-electrode-locs-32_orientwheel.ced"
    lines = ced_path.read_text(encoding="utf-8", errors="ignore").strip().splitlines()
    thetas = []
    radii = []
    labels = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= 4:
            lbl = parts[1].strip()
            if lbl in ("M2", "VEOG", "HEOG") or len(labels) >= 31:
                continue
            th_deg = float(parts[2].strip())
            rd = float(parts[3].strip())
            thetas.append(np.deg2rad(th_deg))
            radii.append(rd)
            labels.append(lbl)
    th = np.array(thetas)
    rd = np.array(radii)
    # In EEGLAB topoplot.m / topoplot_connect.m:
    # [x, y] = pol2cart(Th, Rd); and then plots (y, x) so nose (Th=0) is at +Y (top) and +X is right!
    x_cart = rd * np.cos(th)
    y_cart = rd * np.sin(th)
    plotrad = max(min(1.0, float(np.max(rd)) * 1.02), 0.5)
    squeezefac = 0.5 / plotrad
    x_plot = y_cart * squeezefac
    y_plot = x_cart * squeezefac
    return x_plot, y_plot, labels


def draw_eeglab_head(ax: plt.Axes, rmax: float = 0.5, lw: float = 1.2):
    """Draws EEGLAB cartoon head circle, nose, and ears on `ax`."""
    head_circle = Circle((0.0, 0.0), rmax, facecolor="none", edgecolor="black", linewidth=lw, zorder=4)
    ax.add_patch(head_circle)
    # Nose
    base = rmax - 0.0046
    basex = 0.18 * rmax
    tip = 1.15 * rmax
    tiphw = 0.04 * rmax
    tipr = 0.01 * rmax
    ax.plot(
        [basex, tiphw, 0.0, -tiphw, -basex],
        [base, tip - tipr, tip, tip - tipr, base],
        color="black", linewidth=lw, zorder=4, clip_on=False,
    )
    # Ears
    q = 0.04
    ear_x = np.array([0.492, 0.510, 0.518, 0.5299, 0.5419, 0.54, 0.547, 0.532, 0.510, 0.484])
    ear_y = np.array([q + 0.0555, q + 0.0775, q + 0.0783, q + 0.0746, q + 0.0555, -0.0055, -0.0932, -0.1313, -0.1384, -0.1199])
    ax.plot(ear_x, ear_y, color="black", linewidth=lw, zorder=4, clip_on=False)
    ax.plot(-ear_x, ear_y, color="black", linewidth=lw, zorder=4, clip_on=False)
    ax.set_xlim(-rmax * 1.22, rmax * 1.22)
    ax.set_ylim(-rmax * 1.15, rmax * 1.25)
    ax.set_aspect("equal")
    ax.axis("off")


def plot_scalp_topo(
    ax: plt.Axes,
    values_31: np.ndarray,
    sig_mask_31: np.ndarray,
    x_el: np.ndarray,
    y_el: np.ndarray,
    cmap: LinearSegmentedColormap,
    vmin: float = -1.0,
    vmax: float = 1.0,
    head_lw: float = 1.1,
    dot_size: float = 8.0,
    star_size: float = 28.0,
):
    """
    Renders a smooth interpolated EEGLAB-style scalp topography on `ax`,
    marking all 31 electrodes with black dots and statistically significant
    electrodes (`sig_mask_31 == 1`) with green stars (`#39FF14` / `#10B981`).
    """
    rmax = 0.5
    grid_n = 80
    gx = np.linspace(-rmax, rmax, grid_n)
    gy = np.linspace(-rmax, rmax, grid_n)
    GX, GY = np.meshgrid(gx, gy)

    # Add boundary ghost points on circle r=0.58 using nearest-neighbor weighted average for smooth extrapolation to head rim
    angles = np.linspace(0, 2 * np.pi, 36, endpoint=False)
    bx = 0.58 * np.cos(angles)
    by = 0.58 * np.sin(angles)
    bvals = np.empty_like(angles)
    for i in range(len(angles)):
        d2 = (x_el - bx[i]) ** 2 + (y_el - by[i]) ** 2
        w = 1.0 / np.maximum(d2, 1e-4) ** 2
        bvals[i] = np.sum(w * values_31) / np.sum(w)

    pts_x = np.concatenate([x_el, bx])
    pts_y = np.concatenate([y_el, by])
    pts_v = np.concatenate([values_31, bvals])

    rbf = interp.Rbf(pts_x, pts_y, pts_v, function="multiquadric", smooth=0.02)
    GZ = rbf(GX, GY)
    GZ[GX ** 2 + GY ** 2 > rmax ** 2] = np.nan

    im = ax.imshow(
        GZ,
        extent=(-rmax, rmax, -rmax, rmax),
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
        zorder=1,
    )
    ax.contour(
        GX, GY, GZ,
        levels=5,
        colors="black",
        linewidths=0.35,
        alpha=0.45,
        zorder=2,
    )
    draw_eeglab_head(ax, rmax=rmax, lw=head_lw)

    nonsig = ~sig_mask_31.astype(bool)
    ax.scatter(x_el[nonsig], y_el[nonsig], s=dot_size, c="black", zorder=5, linewidths=0)
    if np.any(sig_mask_31):
        sig_idx = sig_mask_31.astype(bool)
        ax.scatter(
            x_el[sig_idx], y_el[sig_idx],
            s=star_size, marker="*", c="#15C615", edgecolors="black", linewidths=0.35, zorder=6,
        )
    return im


def plot_connectome_topo(
    ax: plt.Axes,
    weights_465: np.ndarray,
    x_el: np.ndarray,
    y_el: np.ndarray,
    cmap: LinearSegmentedColormap,
    thresh: float = 0.30,
    vmin: float = -1.0,
    vmax: float = 1.0,
    head_lw: float = 1.1,
    line_lw: float = 1.1,
    dot_size: float = 8.0,
):
    """
    Exact Python implementation of topoplot_connect.m:
    Draws the EEGLAB cartoon head, 31 electrode dots, and colored lines connecting
    channel pairs with |z| >= thresh on [-1, 1] redblue scale.
    """
    rmax = 0.5
    draw_eeglab_head(ax, rmax=rmax, lw=head_lw)
    ax.scatter(x_el, y_el, s=dot_size, c="black", zorder=5, linewidths=0)

    tri_i, tri_j = np.triu_indices(31, k=1)
    keep = np.where(np.abs(weights_465) >= thresh)[0]
    if len(keep) > 0:
        # Sort by |weight| so strongest connections draw on top
        order = keep[np.argsort(np.abs(weights_465[keep]))]
        for idx in order:
            val = float(np.clip(weights_465[idx], vmin, vmax))
            norm_val = (val - vmin) / (vmax - vmin)
            col = cmap(norm_val)
            c1 = tri_i[idx]
            c2 = tri_j[idx]
            ax.plot(
                [x_el[c1], x_el[c2]],
                [y_el[c1], y_el[c2]],
                color=col,
                linewidth=line_lw,
                alpha=0.88,
                zorder=3,
            )


# ==============================================================================
# 5. RENDERING FIGURES 3, 4, 5, AND 6
# ==============================================================================
def render_figure3_or_5(
    fig_num: int,
    results_dict: Dict[str, Any],
    x_el: np.ndarray,
    y_el: np.ndarray,
) -> Tuple[Path, Path]:
    """
    Renders Figure 3 (SVM Classification, 6 rows x 2 cols) or
    Figure 5 (SVR Regression, 6 rows x 4 cols) with Panel A (time-courses + inset
    significant-window topographies) and Panel B (epoch-averaged Haufe feature weight
    scalp topographies with significance stars).
    """
    is_svm = (fig_num == 3)
    metrics = [
        ("AAC", "AAC"),
        ("lowFreq_amp", "Alpha Amplitude"),
        ("PACtort", "MI"),
        ("PACoz", "MVL_Özkurt"),
        ("PACglm", "PAC_GLM"),
        ("PACplv", "PAC_PLV"),
    ]
    if is_svm:
        cols = [("L", "Left Targets"), ("R", "Right Targets")]
        fig_w, fig_h = 11.5, 13.5
        out_stem = "figure3_svm_cfc"
    else:
        cols = [
            ("Liv", "Non-Informative Left"),
            ("Riv", "Non-Informative Right"),
            ("Lv", "Informative Left"),
            ("Rv", "Informative Right"),
        ]
        fig_w, fig_h = 14.5, 13.5
        out_stem = "figure5_svr_cfc"

    n_rows = len(metrics)
    n_cols = len(cols)
    cmap_rb = get_redblue_cmap(256)

    from matplotlib.transforms import blended_transform_factory

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300, facecolor="white")
    # Top 62% of figure: Panel A (6 rows x n_cols time courses)
    # Bottom 33% of figure: Panel B (6 rows x n_cols scalp topographies)
    gs_a = fig.add_gridspec(
        n_rows, n_cols,
        left=0.175, right=0.97, top=0.95, bottom=0.40,
        hspace=0.38, wspace=0.18 if is_svm else 0.22,
    )
    gs_b = fig.add_gridspec(
        n_rows, n_cols,
        left=0.22 if is_svm else 0.175, right=0.82 if is_svm else 0.91,
        top=0.33, bottom=0.03,
        hspace=0.15, wspace=0.08,
    )

    fig.text(0.012, 0.965, "A", fontsize=16, fontweight="bold", va="top", ha="left")
    fig.text(0.012, 0.350, "B", fontsize=16, fontweight="bold", va="top", ha="left")

    last_im = None
    for r_idx, (m_key, m_label) in enumerate(metrics):
        if m_key == "lowFreq_amp":
            t_ms = np.arange(0, 1000, 20, dtype=np.float64) - 1000.0
            step_w = 20.0
        else:
            t_ms = np.arange(125, 875, 15, dtype=np.float64) - 1000.0
            step_w = 15.0

        for c_idx, (c_key, c_title) in enumerate(cols):
            ax = fig.add_subplot(gs_a[r_idx, c_idx])
            key_prefix = f"{'svm' if is_svm else 'svr'}_{m_key}_{c_key}"
            mean_real = results_dict[f"{key_prefix}_mean_real"]
            se_real = results_dict[f"{key_prefix}_se_real"]
            mean_perm = results_dict[f"{key_prefix}_mean_perm"]
            se_perm = results_dict[f"{key_prefix}_se_perm"]
            h_steps = results_dict[f"{key_prefix}_h_steps"]
            fw_all_z = results_dict[f"{key_prefix}_fw_all_z"]  # (n_subs, 50, 31)

            # Green shaded significant intervals + inset topographies
            runs = find_contiguous_runs(h_steps)
            for run_idx, run in enumerate(runs):
                x0 = t_ms[run[0]] - 0.3 * step_w
                x1 = t_ms[run[-1]] + 0.3 * step_w
                ax.axvspan(x0, x1, color="#55FF44", alpha=0.55, lw=0, zorder=1)
                # Inset topography for significant window (up to 2 insets per panel to prevent crowding)
                if run_idx < 2:
                    avg_run_z = np.mean(fw_all_z[:, run, :], axis=1)  # (n_subs, 31)
                    mean_run_z = np.mean(avg_run_z, axis=0)
                    _, h_chan_run, _, _ = mxt_perm1_phipson_smyth(avg_run_z, n_perm=2000, alph=0.05, tail=1, seed=900 + r_idx * 10 + c_idx)
                    # Place inset in upper space near the significant window
                    norm_x = float(np.clip((np.mean(t_ms[run]) + 1000.0) / 1000.0, 0.12, 0.78))
                    if len(runs) == 2:
                        norm_x = 0.22 if run_idx == 0 else 0.68
                    ins = ax.inset_axes([norm_x - 0.13, 0.56, 0.26, 0.42])
                    plot_scalp_topo(
                        ins, mean_run_z, h_chan_run, x_el, y_el, cmap_rb,
                        vmin=-1.0, vmax=1.0, head_lw=0.7, dot_size=2.5, star_size=12.0,
                    )

            # Plot permutation chance and real decoding curves
            ax.errorbar(
                t_ms, mean_perm, yerr=se_perm,
                color="#7F7F7F", linestyle="--", linewidth=0.9,
                marker="^", markersize=2.2, elinewidth=0.5, capsize=1.2, zorder=2,
                label="Permutation chance" if (r_idx == 0 and c_idx == 0) else None,
            )
            ax.errorbar(
                t_ms, mean_real, yerr=se_real,
                color="black", linestyle="-", linewidth=1.0,
                marker="s", markersize=2.3, elinewidth=0.55, capsize=1.2, zorder=3,
                label="Classifier accuracy" if (is_svm and r_idx == 0 and c_idx == 0) else ("SVR Fisher-z r" if (not is_svm and r_idx == 0 and c_idx == 0) else None),
            )

            if is_svm:
                ax.axhline(50.0, color="#999999", linestyle=":", linewidth=0.6, zorder=1)
                ax.set_ylim(35.0, 72.0)
                ax.set_yticks([40, 50, 60, 70])
            else:
                ax.axhline(0.0, color="#999999", linestyle=":", linewidth=0.6, zorder=1)
                ax.set_ylim(-0.18, 0.22)
                ax.set_yticks([-0.1, 0.0, 0.1, 0.2])

            ax.set_xlim(-1000.0, 0.0)
            ax.set_xticks([-1000, -800, -600, -400, -200, 0])
            ax.tick_params(axis="both", labelsize=7.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if r_idx == 0:
                ax.set_title(c_title, fontsize=9.5, fontweight="bold", pad=6)
            if r_idx == n_rows - 1:
                ax.set_xlabel("Time-steps [ms]", fontsize=8.5, fontweight="bold")
            if c_idx == 0:
                row_trans_a = blended_transform_factory(fig.transFigure, ax.transAxes)
                ax.text(
                    0.132, 0.5, m_label,
                    transform=row_trans_a, fontsize=9.0, fontweight="bold",
                    va="center", ha="right",
                )

            # Panel B: Epoch-averaged Haufe feature weights topography
            ax_b = fig.add_subplot(gs_b[r_idx, c_idx])
            mean_fw_epoch = results_dict[f"{key_prefix}_fw_epoch_mean"]
            h_fw_epoch = results_dict[f"{key_prefix}_fw_epoch_h"]
            last_im = plot_scalp_topo(
                ax_b, mean_fw_epoch, h_fw_epoch, x_el, y_el, cmap_rb,
                vmin=-1.0, vmax=1.0, head_lw=0.95, dot_size=4.5, star_size=22.0,
            )
            if r_idx == 0:
                ax_b.set_title(c_title, fontsize=8.5, fontweight="bold", pad=3)
            if c_idx == 0:
                row_trans_b = blended_transform_factory(fig.transFigure, ax_b.transAxes)
                ax_b.text(
                    0.195 if is_svm else 0.150, 0.5, m_label,
                    transform=row_trans_b, fontsize=8.5, fontweight="bold",
                    va="center", ha="right",
                )

    # Shared Y-axis label for Panel A
    fig.text(
        0.020, 0.675,
        "Classification Accuracy [%]" if is_svm else "Fisher-transformed correlation coeff",
        fontsize=10.5, fontweight="bold", rotation=90, va="center", ha="center",
    )

    # Colorbar for Panel B
    cax = fig.add_axes([0.85 if is_svm else 0.93, 0.08, 0.014, 0.18])
    cb = fig.colorbar(last_im, cax=cax, ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
    cb.ax.tick_params(labelsize=7.5)
    cb.set_label("Z-standardised absolute feature weights", fontsize=8.0, fontweight="bold")

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIGURES_DIR / f"{out_stem}.png"
    pdf_path = FIGURES_DIR / f"{out_stem}.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


def render_figure4_or_6(
    fig_num: int,
    results_dict: Dict[str, Any],
    x_el: np.ndarray,
    y_el: np.ndarray,
) -> Tuple[Path, Path]:
    """
    Renders Figure 4 (SVM Functional Connectivity, 4 rows x 2 cols) or
    Figure 6 (SVR Functional Connectivity, 4 rows x 4 cols) with Panel A
    (time-courses + inset connectomes for significant steps, |z| >= 0.45) and
    Panel B (epoch-averaged connectomes, |z| >= 0.30).
    """
    is_svm = (fig_num == 4)
    metrics = [
        ("corr", "Correlation"),
        ("wPLI", "wPLI"),
        ("PLV", "PLV"),
        ("iCOH", "iCOH"),
    ]
    if is_svm:
        cols = [("L", "Left Targets"), ("R", "Right Targets")]
        fig_w, fig_h = 11.5, 10.5
        out_stem = "figure4_svm_fc"
    else:
        cols = [
            ("Liv", "Non-Informative Left"),
            ("Riv", "Non-Informative Right"),
            ("Lv", "Informative Left"),
            ("Rv", "Informative Right"),
        ]
        fig_w, fig_h = 14.5, 10.5
        out_stem = "figure6_svr_fc"

    n_rows = len(metrics)
    n_cols = len(cols)
    cmap_rb = get_redblue_cmap(64)

    from matplotlib.transforms import blended_transform_factory

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300, facecolor="white")
    gs_a = fig.add_gridspec(
        n_rows, n_cols,
        left=0.175, right=0.97, top=0.94, bottom=0.42,
        hspace=0.40, wspace=0.18 if is_svm else 0.22,
    )
    gs_b = fig.add_gridspec(
        n_rows, n_cols,
        left=0.22 if is_svm else 0.175, right=0.82 if is_svm else 0.91,
        top=0.34, bottom=0.03,
        hspace=0.16, wspace=0.08,
    )

    fig.text(0.012, 0.96, "A", fontsize=16, fontweight="bold", va="top", ha="left")
    fig.text(0.012, 0.36, "B", fontsize=16, fontweight="bold", va="top", ha="left")

    # 12 windows: 187.5 : 50 : 737.5 ms -> -812.5 to -262.5 ms relative to target
    t_ms = np.arange(187.5, 787.0, 50.0, dtype=np.float64) - 1000.0
    step_w = 50.0

    for r_idx, (m_key, m_label) in enumerate(metrics):
        for c_idx, (c_key, c_title) in enumerate(cols):
            ax = fig.add_subplot(gs_a[r_idx, c_idx])
            key_prefix = f"{'svm' if is_svm else 'svr'}_{m_key}_{c_key}"
            mean_real = results_dict[f"{key_prefix}_mean_real"]
            se_real = results_dict[f"{key_prefix}_se_real"]
            mean_perm = results_dict[f"{key_prefix}_mean_perm"]
            se_perm = results_dict[f"{key_prefix}_se_perm"]
            h_steps = results_dict[f"{key_prefix}_h_steps"]
            fw_all_z = results_dict[f"{key_prefix}_fw_all_z"]  # (n_subs, 12, 465)

            runs = find_contiguous_runs(h_steps)
            for run_idx, run in enumerate(runs):
                x0 = t_ms[run[0]] - 0.3 * step_w
                x1 = t_ms[run[-1]] + 0.3 * step_w
                ax.axvspan(x0, x1, color="#55FF44", alpha=0.55, lw=0, zorder=1)
                if run_idx < 2:
                    mean_run_z = np.mean(np.mean(fw_all_z[:, run, :], axis=1), axis=0)
                    norm_x = float(np.clip((np.mean(t_ms[run]) + 1000.0) / 1000.0, 0.14, 0.78))
                    if len(runs) == 2:
                        norm_x = 0.24 if run_idx == 0 else 0.68
                    ins = ax.inset_axes([norm_x - 0.13, 0.54, 0.26, 0.44])
                    plot_connectome_topo(
                        ins, mean_run_z, x_el, y_el, cmap_rb,
                        thresh=0.45, vmin=-1.0, vmax=1.0, head_lw=0.7, line_lw=0.65, dot_size=2.2,
                    )

            ax.errorbar(
                t_ms, mean_perm, yerr=se_perm,
                color="#7F7F7F", linestyle="--", linewidth=0.95,
                marker="^", markersize=2.6, elinewidth=0.55, capsize=1.4, zorder=2,
            )
            ax.errorbar(
                t_ms, mean_real, yerr=se_real,
                color="black", linestyle="-", linewidth=1.05,
                marker="s", markersize=2.7, elinewidth=0.6, capsize=1.4, zorder=3,
            )

            if is_svm:
                ax.axhline(50.0, color="#999999", linestyle=":", linewidth=0.6, zorder=1)
                ax.set_ylim(35.0, 72.0)
                ax.set_yticks([40, 50, 60, 70])
            else:
                ax.axhline(0.0, color="#999999", linestyle=":", linewidth=0.6, zorder=1)
                ax.set_ylim(-0.18, 0.22)
                ax.set_yticks([-0.1, 0.0, 0.1, 0.2])

            ax.set_xlim(-1000.0, 0.0)
            ax.set_xticks([-1000, -800, -600, -400, -200, 0])
            ax.tick_params(axis="both", labelsize=7.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if r_idx == 0:
                ax.set_title(c_title, fontsize=9.5, fontweight="bold", pad=6)
            if r_idx == n_rows - 1:
                ax.set_xlabel("Time-steps [ms]", fontsize=8.5, fontweight="bold")
            if c_idx == 0:
                row_trans_a = blended_transform_factory(fig.transFigure, ax.transAxes)
                ax.text(
                    0.132, 0.5, m_label,
                    transform=row_trans_a, fontsize=9.0, fontweight="bold",
                    va="center", ha="right",
                )

            # Panel B: Epoch-averaged connectome (|z| >= 0.30)
            ax_b = fig.add_subplot(gs_b[r_idx, c_idx])
            mean_fw_epoch = results_dict[f"{key_prefix}_fw_epoch_mean"]
            plot_connectome_topo(
                ax_b, mean_fw_epoch, x_el, y_el, cmap_rb,
                thresh=0.30, vmin=-1.0, vmax=1.0, head_lw=0.95, line_lw=0.85, dot_size=4.5,
            )
            if r_idx == 0:
                ax_b.set_title(c_title, fontsize=8.5, fontweight="bold", pad=3)
            if c_idx == 0:
                row_trans_b = blended_transform_factory(fig.transFigure, ax_b.transAxes)
                ax_b.text(
                    0.195 if is_svm else 0.150, 0.5, m_label,
                    transform=row_trans_b, fontsize=8.5, fontweight="bold",
                    va="center", ha="right",
                )

    fig.text(
        0.020, 0.68,
        "Classification Accuracy [%]" if is_svm else "Fisher-transformed correlation coeff",
        fontsize=10.5, fontweight="bold", rotation=90, va="center", ha="center",
    )

    sm = plt.cm.ScalarMappable(cmap=cmap_rb, norm=plt.Normalize(vmin=-1.0, vmax=1.0))
    sm.set_array([])
    cax = fig.add_axes([0.85 if is_svm else 0.93, 0.08, 0.014, 0.18])
    cb = fig.colorbar(sm, cax=cax, ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
    cb.ax.tick_params(labelsize=7.5)
    cb.set_label("Z-standardised absolute feature weights", fontsize=8.0, fontweight="bold")

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIGURES_DIR / f"{out_stem}.png"
    pdf_path = FIGURES_DIR / f"{out_stem}.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


# ==============================================================================
# 6. MAIN END-TO-END EXECUTION
# ==============================================================================
def main():
    t_start = time.time()
    cache_files = sorted(CACHE_DIR.glob("*_feats.npz"))
    subj_ids = [f.name.split("_")[0] for f in cache_files if f.name.split("_")[0] in INCLUDED_PARTICIPANTS]
    print(f"[run_mvpa_and_figures] Found {len(subj_ids)} cached EEG subjects: {subj_ids}", flush=True)

    # Load all subject caches into memory
    subj_data: Dict[str, Any] = {}
    for sid in subj_ids:
        subj_data[sid] = dict(np.load(CACHE_DIR / f"{sid}_feats.npz"))

    n_subs = len(subj_ids)

    # --------------------------------------------------------------------------
    # STEP A: Group FOOOF Spectral Parameterization Statistics (fooof_OrientExo_ML.m)
    # --------------------------------------------------------------------------
    print("[run_mvpa_and_figures] Running FOOOF 10,000-permutation F_max tests...", flush=True)
    # Shape: (31, 2_side [L=0, R=1], 2_cue [NI=0, I=1], n_subs)
    ap_offset_4d = np.full((31, 2, 2, n_subs), np.nan, dtype=np.float64)
    ap_exp_4d = np.full((31, 2, 2, n_subs), np.nan, dtype=np.float64)
    r2_4d = np.full((31, 2, 2, n_subs), np.nan, dtype=np.float64)
    err_4d = np.full((31, 2, 2, n_subs), np.nan, dtype=np.float64)
    iaf_4d = np.full((31, 2, 2, n_subs), np.nan, dtype=np.float64)

    cond_map = [("Liv", 0, 0), ("Lv", 0, 1), ("Riv", 1, 0), ("Rv", 1, 1)]
    for s_idx, sid in enumerate(subj_ids):
        d = subj_data[sid]
        for cname, side_i, cue_i in cond_map:
            ap = d[f"fooof_ap_{cname}"]
            ap_offset_4d[:, side_i, cue_i, s_idx] = ap[:, 0]
            ap_exp_4d[:, side_i, cue_i, s_idx] = ap[:, 1]
            r2_4d[:, side_i, cue_i, s_idx] = d[f"fooof_r2_{cname}"]
            err_4d[:, side_i, cue_i, s_idx] = d[f"fooof_err_{cname}"]
            iaf_4d[:, side_i, cue_i, s_idx] = d[f"fooof_iaf_{cname}"]

    # Replace any isolated NaNs with channel median across subjects
    for arr4 in (ap_offset_4d, ap_exp_4d, r2_4d, err_4d, iaf_4d):
        nan_mask = ~np.isfinite(arr4)
        if np.any(nan_mask):
            arr4[nan_mask] = np.nanmedian(arr4)

    off_cue = calc_fmax_2x2(ap_offset_4d, "cue", n_perm=10000, alpha=0.05 / 3.0)
    off_side = calc_fmax_2x2(ap_offset_4d, "side", n_perm=10000, alpha=0.05 / 3.0)
    off_int = calc_fmax_2x2(ap_offset_4d, "interaction", n_perm=10000, alpha=0.05 / 3.0)
    # Electrodes [7, 11] in 1-based MATLAB -> [6, 10] in 0-based Python (PO3, P7)
    po3_p7_by_subj = np.mean(ap_offset_4d[[6, 10], :, :, :], axis=1)  # (2_elec, 2_cue, n_subs)
    po3_p7_mean_across_elec = np.mean(po3_p7_by_subj, axis=0)  # (2_cue, n_subs)
    ni_offset_mean = float(np.mean(po3_p7_mean_across_elec[0, :]))
    ni_offset_sd = float(np.std(po3_p7_mean_across_elec[0, :], ddof=1))
    i_offset_mean = float(np.mean(po3_p7_mean_across_elec[1, :]))
    i_offset_sd = float(np.std(po3_p7_mean_across_elec[1, :], ddof=1))

    exp_cue = calc_fmax_2x2(ap_exp_4d, "cue", n_perm=10000, alpha=0.05 / 3.0)
    exp_side = calc_fmax_2x2(ap_exp_4d, "side", n_perm=10000, alpha=0.05 / 3.0)
    exp_int = calc_fmax_2x2(ap_exp_4d, "interaction", n_perm=10000, alpha=0.05 / 3.0)

    r2_cue = calc_fmax_2x2(r2_4d, "cue", n_perm=10000, alpha=0.05 / 3.0)
    r2_side = calc_fmax_2x2(r2_4d, "side", n_perm=10000, alpha=0.05 / 3.0)
    r2_int = calc_fmax_2x2(r2_4d, "interaction", n_perm=10000, alpha=0.05 / 3.0)

    err_cue = calc_fmax_2x2(err_4d, "cue", n_perm=10000, alpha=0.05 / 4.0)
    err_side = calc_fmax_2x2(err_4d, "side", n_perm=10000, alpha=0.05 / 4.0)
    err_int = calc_fmax_2x2(err_4d, "interaction", n_perm=10000, alpha=0.05 / 4.0)

    # Grand IAF mean, SD, min, max across participants (fooof_Avg_OrientExo_ML.m)
    iaf_subj_means = np.nanmean(iaf_4d, axis=(0, 1, 2))
    iaf_grand_mean = float(np.nanmean(iaf_subj_means))
    iaf_grand_sd = float(np.nanstd(iaf_subj_means, ddof=1))
    iaf_min = float(np.nanmin(iaf_subj_means))
    iaf_max = float(np.nanmax(iaf_subj_means))

    off_mean_crit = float(np.mean([off_cue["Fmax_crit"], off_side["Fmax_crit"], off_int["Fmax_crit"]]))
    exp_mean_crit = float(np.mean([exp_cue["Fmax_crit"], exp_side["Fmax_crit"], exp_int["Fmax_crit"]]))
    r2_mean_crit = float(np.mean([r2_cue["Fmax_crit"], r2_side["Fmax_crit"], r2_int["Fmax_crit"]]))
    err_mean_crit = float(np.mean([err_cue["Fmax_crit"], err_side["Fmax_crit"], err_int["Fmax_crit"]]))

    api = OrientExoPipelineAPI()
    api.handle_request({
        "action": "update_stats",
        "source_script": "scripts/run_mvpa_and_figures.py",
        "status": "computed_from_raw_data",
        "updates": {
            "fooof": {
                "offset_mean_crit_fmax": round(off_mean_crit, 2),
                "exponent_mean_crit_fmax": round(exp_mean_crit, 2),
                "r2_mean_crit_fmax": round(r2_mean_crit, 2),
                "error_mean_crit_fmax": round(err_mean_crit, 2),
                "iaf_grand_mean_hz": f"{iaf_grand_mean:.1f}",
                "iaf_grand_sd_hz": f"{iaf_grand_sd:.1f}",
                "iaf_range_hz": f"{iaf_min:.1f}--{iaf_max:.1f}",
                "aperiodic_offset": {
                    "cue_validity_fmax_crit": round(off_cue["Fmax_crit"], 2),
                    "target_side_fmax_crit": round(off_side["Fmax_crit"], 2),
                    "interaction_fmax_crit": round(off_int["Fmax_crit"], 2),
                    "non_informative_mean_log_power": round(ni_offset_mean, 2),
                    "non_informative_sd_log_power": round(ni_offset_sd, 2),
                    "informative_mean_log_power": round(i_offset_mean, 2),
                    "informative_sd_log_power": round(i_offset_sd, 2),
                },
            }
        },
    })
    print(f"[run_mvpa_and_figures] Updated FOOOF stats in analysis_results.json (IAF M={iaf_grand_mean:.2f} Hz).", flush=True)

    # --------------------------------------------------------------------------
    # STEP B: Spatiotemporal DDTBOX SVM & SVR Decoding + Haufe Feature Weights
    # --------------------------------------------------------------------------
    mvpa_cache_file = CACHE_DIR / "group_mvpa_results.npz"
    if mvpa_cache_file.exists():
        print(f"[run_mvpa_and_figures] Loading cached group MVPA results from {mvpa_cache_file}...", flush=True)
        results_dict = dict(np.load(mvpa_cache_file))
    else:
        results_dict = {}
        all_metrics = [
            "lowFreq_amp", "AAC", "PACoz", "PACtort", "PACglm", "PACplv",
            "wPLI", "corr", "PLV", "iCOH",
        ]
        # Warm up Numba JIT compilers on tiny synthetic arrays
        print("[run_mvpa_and_figures] Warming up Numba SVM and SVR JIT compilers...", flush=True)
        _dummy = np.random.default_rng(0).standard_normal((2, 5, 4, 20)).astype(np.float32)
        _ = ddtbox_svm_subject(_dummy, _dummy, seed=1)
        _ = ddtbox_svr_subject(_dummy, np.arange(20, dtype=np.float64), seed=1)

        for m_idx, m_key in enumerate(all_metrics, 1):
            t_m = time.time()
            print(f"[run_mvpa_and_figures] [{m_idx}/10] Running SVM & SVR across {n_subs} subjects for `{m_key}`...", flush=True)

            # 1. SVM Classification: Left (Liv vs Lv) and Right (Riv vs Rv)
            for side_key, c1_name, c2_name in [("L", "Liv", "Lv"), ("R", "Riv", "Rv")]:
                sample_arr = subj_data[subj_ids[0]][f"feat_{m_key}_{c1_name}"]
                n_steps, _, n_chan, _ = sample_arr.shape
                all_real = np.zeros((n_subs, n_steps), dtype=np.float64)
                all_perm_reps = np.zeros((n_subs, n_steps, 10, 10), dtype=np.float64)
                all_fw_z = np.zeros((n_subs, n_steps, n_chan), dtype=np.float64)

                for s_idx, sid in enumerate(subj_ids):
                    X1 = subj_data[sid][f"feat_{m_key}_{c1_name}"]
                    X2 = subj_data[sid][f"feat_{m_key}_{c2_name}"]
                    ra, pa, fz = ddtbox_svm_subject(X1, X2, seed=1000 + m_idx * 100 + s_idx)
                    all_real[s_idx] = ra
                    all_perm_reps[s_idx] = pa
                    all_fw_z[s_idx] = fz

                all_perm_mean = np.mean(all_perm_reps, axis=(2, 3))
                # Draw random repetition per subject as in tmax_classifier_accuracies.m (default: mean or drawn perm)
                diff_scores = all_real - all_perm_mean
                p_steps, h_steps, t_steps, tcrit_steps = mxt_perm1_phipson_smyth(
                    diff_scores, n_perm=10000, alph=0.05, tail=1, seed=5000 + m_idx * 10 + (0 if side_key == "L" else 1)
                )

                # Epoch-averaged feature weights (steps 1..n_steps)
                fw_epoch_subj = np.mean(all_fw_z, axis=1)  # (n_subs, n_chan)
                p_fw, h_fw, _, _ = mxt_perm1_phipson_smyth(
                    fw_epoch_subj, n_perm=10000, alph=0.05, tail=1, seed=6000 + m_idx * 10 + (0 if side_key == "L" else 1)
                )

                kp = f"svm_{m_key}_{side_key}"
                results_dict[f"{kp}_mean_real"] = np.mean(all_real, axis=0)
                results_dict[f"{kp}_se_real"] = np.std(all_real, axis=0, ddof=1) / np.sqrt(n_subs)
                results_dict[f"{kp}_mean_perm"] = np.mean(all_perm_mean, axis=0)
                results_dict[f"{kp}_se_perm"] = np.std(all_perm_mean, axis=0, ddof=1) / np.sqrt(n_subs)
                results_dict[f"{kp}_p_steps"] = p_steps
                results_dict[f"{kp}_h_steps"] = h_steps
                results_dict[f"{kp}_t_steps"] = t_steps
                results_dict[f"{kp}_fw_all_z"] = all_fw_z
                results_dict[f"{kp}_fw_epoch_mean"] = np.mean(fw_epoch_subj, axis=0)
                results_dict[f"{kp}_fw_epoch_h"] = h_fw
                results_dict[f"{kp}_fw_epoch_p"] = p_fw

            # 2. SVR Regression: Liv, Riv, Lv, Rv
            for cname in ["Liv", "Riv", "Lv", "Rv"]:
                sample_arr = subj_data[subj_ids[0]][f"feat_{m_key}_{cname}"]
                n_steps, _, n_chan, _ = sample_arr.shape
                all_real = np.zeros((n_subs, n_steps), dtype=np.float64)
                all_perm_reps = np.zeros((n_subs, n_steps, 10, 10), dtype=np.float64)
                all_fw_z = np.zeros((n_subs, n_steps, n_chan), dtype=np.float64)

                for s_idx, sid in enumerate(subj_ids):
                    Xc = subj_data[sid][f"feat_{m_key}_{cname}"]
                    yc = subj_data[sid][f"SVR_labels_{cname}"]
                    rr, pr, fz = ddtbox_svr_subject(Xc, yc, seed=3000 + m_idx * 100 + s_idx)
                    all_real[s_idx] = rr
                    all_perm_reps[s_idx] = pr
                    all_fw_z[s_idx] = fz

                all_perm_mean = np.mean(all_perm_reps, axis=(2, 3))
                diff_scores = all_real - all_perm_mean
                p_steps, h_steps, t_steps, _ = mxt_perm1_phipson_smyth(
                    diff_scores, n_perm=10000, alph=0.05, tail=1, seed=7000 + m_idx * 10
                )

                fw_epoch_subj = np.mean(all_fw_z, axis=1)
                p_fw, h_fw, _, _ = mxt_perm1_phipson_smyth(
                    fw_epoch_subj, n_perm=10000, alph=0.05, tail=1, seed=8000 + m_idx * 10
                )

                kp = f"svr_{m_key}_{cname}"
                results_dict[f"{kp}_mean_real"] = np.mean(all_real, axis=0)
                results_dict[f"{kp}_se_real"] = np.std(all_real, axis=0, ddof=1) / np.sqrt(n_subs)
                results_dict[f"{kp}_mean_perm"] = np.mean(all_perm_mean, axis=0)
                results_dict[f"{kp}_se_perm"] = np.std(all_perm_mean, axis=0, ddof=1) / np.sqrt(n_subs)
                results_dict[f"{kp}_p_steps"] = p_steps
                results_dict[f"{kp}_h_steps"] = h_steps
                results_dict[f"{kp}_t_steps"] = t_steps
                results_dict[f"{kp}_fw_all_z"] = all_fw_z
                results_dict[f"{kp}_fw_epoch_mean"] = np.mean(fw_epoch_subj, axis=0)
                results_dict[f"{kp}_fw_epoch_h"] = h_fw
                results_dict[f"{kp}_fw_epoch_p"] = p_fw

            print(f"  -> Completed `{m_key}` in {time.time() - t_m:.1f}s", flush=True)

        np.savez(mvpa_cache_file, **results_dict)
        print(f"[run_mvpa_and_figures] Saved group MVPA results cache to {mvpa_cache_file}", flush=True)

    # --------------------------------------------------------------------------
    # STEP C: Render Publication-Ready Figures 3, 4, 5, and 6
    # --------------------------------------------------------------------------
    x_el, y_el, _ = get_2d_topo_coords()

    print("[run_mvpa_and_figures] Rendering Figure 3 (SVM Alpha & CFC)...", flush=True)
    f3_png, _ = render_figure3_or_5(3, results_dict, x_el, y_el)
    api.handle_request({
        "action": "register_figure",
        "figure_id": "figure3",
        "file": f"figures/{f3_png.name}",
        "status": "generated_from_analysis",
        "source_script": "scripts/run_mvpa_and_figures.py",
    })

    print("[run_mvpa_and_figures] Rendering Figure 4 (SVM Functional Connectivity)...", flush=True)
    f4_png, _ = render_figure4_or_6(4, results_dict, x_el, y_el)
    api.handle_request({
        "action": "register_figure",
        "figure_id": "figure4",
        "file": f"figures/{f4_png.name}",
        "status": "generated_from_analysis",
        "source_script": "scripts/run_mvpa_and_figures.py",
    })

    print("[run_mvpa_and_figures] Rendering Figure 5 (SVR Alpha & CFC)...", flush=True)
    f5_png, _ = render_figure3_or_5(5, results_dict, x_el, y_el)
    api.handle_request({
        "action": "register_figure",
        "figure_id": "figure5",
        "file": f"figures/{f5_png.name}",
        "status": "generated_from_analysis",
        "source_script": "scripts/run_mvpa_and_figures.py",
    })

    print("[run_mvpa_and_figures] Rendering Figure 6 (SVR Functional Connectivity)...", flush=True)
    f6_png, _ = render_figure4_or_6(6, results_dict, x_el, y_el)
    api.handle_request({
        "action": "register_figure",
        "figure_id": "figure6",
        "file": f"figures/{f6_png.name}",
        "status": "generated_from_analysis",
        "source_script": "scripts/run_mvpa_and_figures.py",
    })

    # --------------------------------------------------------------------------
    # STEP D: Update MVPA & MVPA Summary Metrics in analysis_results.json + Render Fig 7
    # --------------------------------------------------------------------------
    # Order matching analysis_results.json["mvpa_summary_metrics"]["labels"]:
    # ["Alpha Amplitude", "CFC: MVL (Özkurt)", "CFC: AAC", "CFC: MI", "CFC: PAC (GLM)",
    #  "CFC: PAC (PLV)", "FC: wPLI", "FC: Correlation", "FC: PLV", "FC: iCOH"]
    summary_order = [
        ("Alpha Amplitude", "lowFreq_amp"),
        ("CFC: MVL (Özkurt)", "PACoz"),
        ("CFC: AAC", "AAC"),
        ("CFC: MI", "PACtort"),
        ("CFC: PAC (GLM)", "PACglm"),
        ("CFC: PAC (PLV)", "PACplv"),
        ("FC: wPLI", "wPLI"),
        ("FC: Correlation", "corr"),
        ("FC: PLV", "PLV"),
        ("FC: iCOH", "iCOH"),
    ]
    svm_left_peak = []
    svm_left_sem = []
    svm_right_peak = []
    svm_right_sem = []
    svm_left_pct_sig = []
    svm_right_pct_sig = []
    svm_left_sig_indices = []
    svm_right_sig_indices = []
    svr_noninf_peak = []
    svr_noninf_sem = []
    svr_inf_peak = []
    svr_inf_sem = []
    svr_noninf_sig_indices = []
    svr_inf_sig_indices = []

    for idx_m, (_, m_key) in enumerate(summary_order):
        # SVM Left
        m_l = results_dict[f"svm_{m_key}_L_mean_real"]
        se_l = results_dict[f"svm_{m_key}_L_se_real"]
        h_l = results_dict[f"svm_{m_key}_L_h_steps"]
        pk_l_idx = int(np.argmax(m_l))
        svm_left_peak.append(round(float(m_l[pk_l_idx]), 1))
        svm_left_sem.append(round(float(se_l[pk_l_idx]), 1))
        pct_l = round(100.0 * float(np.sum(h_l)) / len(h_l), 1)
        svm_left_pct_sig.append(pct_l)
        if pct_l > 0:
            svm_left_sig_indices.append(idx_m)

        # SVM Right
        m_r = results_dict[f"svm_{m_key}_R_mean_real"]
        se_r = results_dict[f"svm_{m_key}_R_se_real"]
        h_r = results_dict[f"svm_{m_key}_R_h_steps"]
        pk_r_idx = int(np.argmax(m_r))
        svm_right_peak.append(round(float(m_r[pk_r_idx]), 1))
        svm_right_sem.append(round(float(se_r[pk_r_idx]), 1))
        pct_r = round(100.0 * float(np.sum(h_r)) / len(h_r), 1)
        svm_right_pct_sig.append(pct_r)
        if pct_r > 0:
            svm_right_sig_indices.append(idx_m)

        # SVR Non-Informative (max across Liv and Riv)
        m_liv = results_dict[f"svr_{m_key}_Liv_mean_real"]
        se_liv = results_dict[f"svr_{m_key}_Liv_se_real"]
        h_liv = results_dict[f"svr_{m_key}_Liv_h_steps"]
        m_riv = results_dict[f"svr_{m_key}_Riv_mean_real"]
        se_riv = results_dict[f"svr_{m_key}_Riv_se_real"]
        h_riv = results_dict[f"svr_{m_key}_Riv_h_steps"]
        if np.max(m_liv) >= np.max(m_riv):
            pk_ni = int(np.argmax(m_liv))
            svr_noninf_peak.append(round(float(m_liv[pk_ni]), 3))
            svr_noninf_sem.append(round(float(se_liv[pk_ni]), 3))
        else:
            pk_ni = int(np.argmax(m_riv))
            svr_noninf_peak.append(round(float(m_riv[pk_ni]), 3))
            svr_noninf_sem.append(round(float(se_riv[pk_ni]), 3))
        if np.any(h_liv) or np.any(h_riv):
            svr_noninf_sig_indices.append(idx_m)

        # SVR Informative (max across Lv and Rv)
        m_lv = results_dict[f"svr_{m_key}_Lv_mean_real"]
        se_lv = results_dict[f"svr_{m_key}_Lv_se_real"]
        h_lv = results_dict[f"svr_{m_key}_Lv_h_steps"]
        m_rv = results_dict[f"svr_{m_key}_Rv_mean_real"]
        se_rv = results_dict[f"svr_{m_key}_Rv_se_real"]
        h_rv = results_dict[f"svr_{m_key}_Rv_h_steps"]
        if np.max(m_lv) >= np.max(m_rv):
            pk_i = int(np.argmax(m_lv))
            svr_inf_peak.append(round(float(m_lv[pk_i]), 3))
            svr_inf_sem.append(round(float(se_lv[pk_i]), 3))
        else:
            pk_i = int(np.argmax(m_rv))
            svr_inf_peak.append(round(float(m_rv[pk_i]), 3))
            svr_inf_sem.append(round(float(se_rv[pk_i]), 3))
        if np.any(h_lv) or np.any(h_rv):
            svr_inf_sig_indices.append(idx_m)

    def format_sig_windows(h_arr: np.ndarray, t_arr: np.ndarray) -> str:
        runs = find_contiguous_runs(h_arr)
        parts = []
        for r in runs:
            t0 = int(round(float(t_arr[r[0]])))
            t1 = int(round(float(t_arr[r[-1]])))
            if t0 == t1:
                parts.append(f"{t0} ms")
            else:
                parts.append(f"{t0} to {t1} ms")
        if not parts:
            return "none"
        if len(parts) == 1:
            return parts[0]
        if len(parts) == 2:
            return f"{parts[0]} and {parts[1]}"
        return ", ".join(parts[:-1]) + f", and {parts[-1]}"

    t_alpha = np.arange(0, 1000, 20, dtype=np.float64) - 1000.0
    alpha_pct_min = int(round(min(svm_left_pct_sig[0], svm_right_pct_sig[0])))
    alpha_pct_max = int(round(max(svm_left_pct_sig[0], svm_right_pct_sig[0])))
    mvl_pct_min = int(round(min(svm_left_pct_sig[1], svm_right_pct_sig[1])))
    mvl_pct_max = int(round(max(svm_left_pct_sig[1], svm_right_pct_sig[1])))

    api.handle_request({
        "action": "update_stats",
        "source_script": "scripts/run_mvpa_and_figures.py",
        "status": "computed_from_raw_data",
        "updates": {
            "mvpa": {
                "n_permutations": 10000,
                "alpha_level": 0.05,
                "svm_results": {
                    "alpha_pct_sig_range": f"{alpha_pct_min}\\%--{alpha_pct_max}\\%",
                    "alpha_pct_sig_min": alpha_pct_min,
                    "alpha_left_sig_windows": format_sig_windows(results_dict["svm_lowFreq_amp_L_h_steps"], t_alpha),
                    "alpha_right_sig_windows": format_sig_windows(results_dict["svm_lowFreq_amp_R_h_steps"], t_alpha),
                    "mvl_pct_sig_range": f"{mvl_pct_min}\\%--{mvl_pct_max}\\%",
                },
            },
            "mvpa_summary_metrics": {
                "labels": [lbl for (lbl, _) in summary_order],
                "svm_left_peak": svm_left_peak,
                "svm_left_sem": svm_left_sem,
                "svm_right_peak": svm_right_peak,
                "svm_right_sem": svm_right_sem,
                "svm_left_pct_sig": svm_left_pct_sig,
                "svm_right_pct_sig": svm_right_pct_sig,
                "svm_left_sig_indices": svm_left_sig_indices,
                "svm_right_sig_indices": svm_right_sig_indices,
                "svr_noninf_peak": svr_noninf_peak,
                "svr_noninf_sem": svr_noninf_sem,
                "svr_inf_peak": svr_inf_peak,
                "svr_inf_sem": svr_inf_sem,
                "svr_noninf_sig_indices": svr_noninf_sig_indices,
                "svr_inf_sig_indices": svr_inf_sig_indices,
            },
        },
    })
    api.handle_request({
        "action": "generate_summary_figure",
        "status": "generated_from_analysis",
        "source_script": "scripts/run_mvpa_and_figures.py",
    })

    print(f"[run_mvpa_and_figures] Done in {time.time() - t_start:.1f}s!", flush=True)


if __name__ == "__main__":
    main()
