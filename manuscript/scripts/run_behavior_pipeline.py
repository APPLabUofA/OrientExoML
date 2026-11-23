#!/usr/bin/env python3
"""
run_behavior_pipeline.py

Reconstructs the empirical behavioral analysis and Figure 2 directly from the
raw single-trial MATLAB behavioral data files (OSF: https://osf.io/u6fgm/).
Implements the 2x2 repeated-measures factorial permutation ANOVA (Fmax and tmax,
N = 10,000 permutations) following `perm_rbANOVA_OrientExo.m` from `CuedOrientTask_paper`.

Outputs:
1. Recomputed behavioral statistics (Fmax, crit_Fmax, p-values, mean errors, staircase luminance).
2. Publication-ready Figure 2 (300 DPI PNG + vector PDF) with clean scientific figure styling.
3. Automatically updates `analysis_results.json` and `generated_stats.tex` via `pipeline_api.py`.
"""

import sys
import json
from pathlib import Path
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPTS_DIR = Path(__file__).resolve().parent
MANUSCRIPT_DIR = SCRIPTS_DIR.parent
PROJECT_ROOT = MANUSCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data" / "behavior"
CACHE_DIR = PROJECT_ROOT / "data" / "cache"
FIGURES_DIR = MANUSCRIPT_DIR / "figures"

sys.path.insert(0, str(SCRIPTS_DIR))
from pipeline_api import OrientExoPipelineAPI
from eeg_preproc_and_feats import VISUAL_INSPECT_REJ_1BASED

# Included participants per AnalysisWrapper_OrientWheelExo_ML.m (N = 28)
INCLUDED_PARTICIPANTS = [
    "004", "005", "007", "009", "010",
    "012", "013", "014", "015", "017",
    "018", "019", "021", "022", "024",
    "025", "026", "027", "029", "031",
    "032", "034", "035", "036", "037",
    "038", "039", "040"
]

def get_arr(val):
    if hasattr(val, "item") and getattr(val, "ndim", 0) == 0:
        val = val.item()
    return np.asarray(val).ravel()

def load_behavior_data(restrict_to_eeg_trials=True):
    n_subs = len(INCLUDED_PARTICIPANTS)
    # Shape: (n_subs, cue_type: [0=info, 1=noninfo], side: [0=left, 1=right])
    abs_err = np.zeros((n_subs, 2, 2))
    resp_err = np.zeros((n_subs, 2, 2))
    staircase_grays = []

    for i, sub in enumerate(INCLUDED_PARTICIPANTS):
        fpath = DATA_DIR / f"{sub}_Orient_Exo.mat"
        if not fpath.exists():
            raise FileNotFoundError(f"Missing behavioral data file: {fpath}")
        mat = sio.loadmat(fpath, squeeze_me=True)
        data = mat["data"]
        prefs = mat["prefs"]
        tg = prefs["targ_gray"]
        staircase_grays.append(float(tg.item() if hasattr(tg, "item") else tg))

        cpath = CACHE_DIR / f"{sub}_feats.npz"
        if restrict_to_eeg_trials and cpath.exists():
            c_npz = np.load(cpath)
            # Note: data.errorDegrees in MATLAB is (Response - Target); negate to match (Target - Response) in main.tex & Fig 1B
            err = -np.asarray(c_npz["err_deg"], dtype=np.float64)
            val = np.asarray(c_npz["valid"], dtype=np.int32)
            pos = np.asarray(c_npz["position"], dtype=np.int32)
        else:
            # Remove first 20 practice trials
            err = -get_arr(data["errorDegrees"])[20:].astype(np.float64)
            val = get_arr(data["valid"])[20:].astype(np.int32)
            pos = get_arr(data["position"])[20:].astype(np.int32)

            # Participant 009 recording error removal per getBEHdata_OrientWheel_Exo.m
            if sub == "009":
                err = err[:343]
                val = val[:343]
                pos = pos[:343]

            if restrict_to_eeg_trials:
                vis_rej = VISUAL_INSPECT_REJ_1BASED.get(sub, [])
                keep = np.ones(len(err), dtype=bool)
                for idx_1b in vis_rej:
                    if 1 <= idx_1b <= len(keep):
                        keep[idx_1b - 1] = False
                err, val, pos = err[keep], val[keep], pos[keep]

        # Info Left: val == 1, pos == 0
        resp_err[i, 0, 0] = np.mean(err[(val == 1) & (pos == 0)])
        abs_err[i, 0, 0] = np.mean(np.abs(err[(val == 1) & (pos == 0)]))
        # Info Right: val == 1, pos == 1
        resp_err[i, 0, 1] = np.mean(err[(val == 1) & (pos == 1)])
        abs_err[i, 0, 1] = np.mean(np.abs(err[(val == 1) & (pos == 1)]))
        # Non-info Left: val == 0, pos == 0
        resp_err[i, 1, 0] = np.mean(err[(val == 0) & (pos == 0)])
        abs_err[i, 1, 0] = np.mean(np.abs(err[(val == 0) & (pos == 0)]))
        # Non-info Right: val == 0, pos == 1
        resp_err[i, 1, 1] = np.mean(err[(val == 0) & (pos == 1)])
        abs_err[i, 1, 1] = np.mean(np.abs(err[(val == 0) & (pos == 1)]))

    return abs_err, resp_err, np.array(staircase_grays)


def perm_rm_anova_2x2(data, n_perms=10000, seed=42):
    """
    Computes 2x2 repeated measures ANOVA with exact FMUT permutation-based Fmax testing
    matching `perm_rbANOVA_OrientExo.m` + `reduce_data_ModelFit.m` (Fields & Kuperberg, 2019).
    Collapses to 1-way 2-level difference vectors for each main effect and the 2x2 interaction,
    and performs within-subject sign-flipping (with permutation 0 = observed sample).
    """
    rng = np.random.RandomState(seed)
    N = data.shape[0]

    diff_cue = np.mean(data[:, 0, :] - data[:, 1, :], axis=1)
    diff_side = np.mean(data[:, :, 0] - data[:, :, 1], axis=1)
    diff_int = (data[:, 0, 0] - data[:, 0, 1]) - (data[:, 1, 0] - data[:, 1, 1])

    def run_1way_fmax(diff_vec, local_rng):
        signs = local_rng.choice([-1.0, 1.0], size=(n_perms, N))
        signs[0, :] = 1.0  # Permutation 0 is the observed data per FMUT
        perm_diffs = signs * diff_vec[None, :]
        means = np.mean(perm_diffs, axis=1)
        stds = np.std(perm_diffs, ddof=1, axis=1)
        t_vals = np.where(stds > 1e-15, means / (stds / np.sqrt(N)), 0.0)
        F_dist = t_vals ** 2
        F_obs = float(F_dist[0])
        crit_F = float(np.percentile(F_dist, 95))
        p_val = float(np.mean(F_dist >= F_obs))
        return {
            "fmax": round(F_obs, 2),
            "crit_fmax": round(crit_F, 2),
            "p": "< 0.001" if p_val < 0.001 else f"{p_val:.3f}",
        }

    return {
        "cue": run_1way_fmax(diff_cue, np.random.RandomState(seed + 1)),
        "side": run_1way_fmax(diff_side, np.random.RandomState(seed + 2)),
        "int": run_1way_fmax(diff_int, np.random.RandomState(seed + 3)),
    }


def pairwise_tmax_by_fdr(data, n_perms=10000, seed=42):
    """
    Executes two-tailed repeated-measures tmax permutation tests (Blair & Karniski, 1993;
    mxt_perm1) across the 4 pairwise simple-effect comparisons and applies the
    Benjamini & Yekutieli (2001) FDR correction at alpha = 0.05.
    """
    N = data.shape[0]
    pairs = {
        "left_info_vs_noninfo": data[:, 0, 0] - data[:, 1, 0],
        "right_info_vs_noninfo": data[:, 0, 1] - data[:, 1, 1],
        "info_left_vs_right": data[:, 0, 0] - data[:, 0, 1],
        "noninfo_left_vs_right": data[:, 1, 0] - data[:, 1, 1],
    }
    raw_p = {}
    t_obs_dict = {}
    for idx, (k, dvec) in enumerate(pairs.items()):
        rng = np.random.RandomState(seed + 100 + idx)
        signs = rng.choice([-1.0, 1.0], size=(n_perms, N))
        signs[0, :] = 1.0
        p_diffs = signs * dvec[None, :]
        m = np.mean(p_diffs, axis=1)
        s = np.std(p_diffs, ddof=1, axis=1)
        t_dist = np.where(s > 1e-15, m / (s / np.sqrt(N)), 0.0)
        t_obs = float(t_dist[0])
        t_obs_dict[k] = round(t_obs, 2)
        raw_p[k] = float(np.mean(np.abs(t_dist) >= abs(t_obs)))

    # Benjamini & Yekutieli (2001) step-up FDR adjustment across the m=4 pairwise tests
    m_tests = len(pairs)
    c_m = sum(1.0 / i for i in range(1, m_tests + 1))
    sorted_keys = sorted(raw_p.keys(), key=lambda k: raw_p[k])
    adj_p = {}
    running_min = 1.0
    for rank_rev in range(m_tests - 1, -1, -1):
        k = sorted_keys[rank_rev]
        rank = rank_rev + 1
        val = min(1.0, raw_p[k] * m_tests * c_m / rank)
        running_min = min(running_min, val)
        adj_p[k] = round(float(running_min), 4)

    out = {}
    for k in pairs:
        out[k] = {
            "t": t_obs_dict[k],
            "p_perm": round(raw_p[k], 4),
            "p_by_fdr": adj_p[k],
            "sig": bool(adj_p[k] < 0.05),
        }
    return out


def plot_publication_figure2(abs_err, resp_err, staircase_grays, pairwise_abs, out_path):
    """
    Renders high-resolution publication-quality Figure 2.
    Clean scientific style conforming to Cereb Cortex and Reviewer 2.
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1.1, 1.1, 0.8]})

    # Panel A: Absolute Response Error
    ax_a = axes[0]
    # Data means & SEMs
    means_abs = [
        [np.mean(abs_err[:, 0, 0]), np.mean(abs_err[:, 1, 0])], # Left: Info, Non-info
        [np.mean(abs_err[:, 0, 1]), np.mean(abs_err[:, 1, 1])], # Right: Info, Non-info
    ]
    sems_abs = [
        [np.std(abs_err[:, 0, 0], ddof=1)/np.sqrt(28), np.std(abs_err[:, 1, 0], ddof=1)/np.sqrt(28)],
        [np.std(abs_err[:, 0, 1], ddof=1)/np.sqrt(28), np.std(abs_err[:, 1, 1], ddof=1)/np.sqrt(28)],
    ]

    x = np.array([0, 1])
    width = 0.32
    ax_a.bar(x - width/2, [means_abs[0][0], means_abs[1][0]], width, yerr=[sems_abs[0][0], sems_abs[1][0]],
             capsize=3, color="#2b5c8f", label="Informative Cue", edgecolor="black", linewidth=0.8)
    ax_a.bar(x + width/2, [means_abs[0][1], means_abs[1][1]], width, yerr=[sems_abs[0][1], sems_abs[1][1]],
             capsize=3, color="#7fcdbb", label="Non-Informative Cue", edgecolor="black", linewidth=0.8)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(["Left Target", "Right Target"])
    ax_a.set_ylabel("Absolute Response Error (°)")
    ax_a.set_ylim(0, 52)
    ax_a.set_title("A. Absolute Orientation Error", fontweight="bold")
    ax_a.legend(frameon=False, loc="upper right")

    # Data-driven significance brackets from Benjamini-Yekutieli FDR-corrected pairwise tmax tests
    lbl_l = "*" if pairwise_abs["left_info_vs_noninfo"]["sig"] else "ns"
    fs_l = 13 if lbl_l == "*" else 9.5
    y_top_l = max(means_abs[0][0] + sems_abs[0][0], means_abs[0][1] + sems_abs[0][1]) + 2.2
    ax_a.plot([x[0] - width/2, x[0] - width/2, x[0] + width/2, x[0] + width/2],
              [y_top_l - 0.8, y_top_l, y_top_l, y_top_l - 0.8], color="black", linewidth=1.0)
    ax_a.text(x[0], y_top_l + 0.4, lbl_l, ha="center", va="bottom", fontsize=fs_l, fontweight="bold")

    lbl_r = "*" if pairwise_abs["right_info_vs_noninfo"]["sig"] else "ns"
    fs_r = 13 if lbl_r == "*" else 9.5
    y_top_r = max(means_abs[1][0] + sems_abs[1][0], means_abs[1][1] + sems_abs[1][1]) + 2.2
    ax_a.plot([x[1] - width/2, x[1] - width/2, x[1] + width/2, x[1] + width/2],
              [y_top_r - 0.8, y_top_r, y_top_r, y_top_r - 0.8], color="black", linewidth=1.0)
    ax_a.text(x[1], y_top_r + 0.5, lbl_r, ha="center", va="bottom", fontsize=fs_r, fontweight="bold")

    # Panel B: Signed Response Error
    ax_b = axes[1]
    means_sign = [
        [np.mean(resp_err[:, 0, 0]), np.mean(resp_err[:, 1, 0])],
        [np.mean(resp_err[:, 0, 1]), np.mean(resp_err[:, 1, 1])],
    ]
    sems_sign = [
        [np.std(resp_err[:, 0, 0], ddof=1)/np.sqrt(28), np.std(resp_err[:, 1, 0], ddof=1)/np.sqrt(28)],
        [np.std(resp_err[:, 0, 1], ddof=1)/np.sqrt(28), np.std(resp_err[:, 1, 1], ddof=1)/np.sqrt(28)],
    ]
    ax_b.bar(x - width/2, [means_sign[0][0], means_sign[1][0]], width, yerr=[sems_sign[0][0], sems_sign[1][0]],
             capsize=3, color="#2b5c8f", label="Informative Cue", edgecolor="black", linewidth=0.8)
    ax_b.bar(x + width/2, [means_sign[0][1], means_sign[1][1]], width, yerr=[sems_sign[0][1], sems_sign[1][1]],
             capsize=3, color="#7fcdbb", label="Non-Informative Cue", edgecolor="black", linewidth=0.8)
    ax_b.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(["Left Target", "Right Target"])
    ax_b.set_ylabel("Signed Response Error (°)")
    ax_b.set_ylim(-6, 6)
    ax_b.set_title("B. Directional Bias", fontweight="bold")

    # Panel C: Staircase Luminance Distribution
    ax_c = axes[2]
    ax_c.hist(staircase_grays, bins=10, color="#41b6c4", edgecolor="black", linewidth=0.8)
    ax_c.axvline(np.mean(staircase_grays), color="#d95f0e", linestyle="--", linewidth=1.5,
                label=f"Mean: {np.mean(staircase_grays):.1f}")
    ax_c.set_xlabel("Target Luminance (8-bit Gray)")
    ax_c.set_ylabel("Participant Count")
    ax_c.set_title("C. Staircase Threshold", fontweight="bold")
    ax_c.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[run_behavior_pipeline] Generated {out_path}")


def main():
    print("[run_behavior_pipeline] Loading EEG-artifact-rejected single-trial behavioral datasets for N=28...")
    abs_err, resp_err, staircase_grays = load_behavior_data(restrict_to_eeg_trials=True)

    print("[run_behavior_pipeline] Running FMUT permutation-based repeated measures ANOVA (10,000 perms)...")
    abs_anova = perm_rm_anova_2x2(abs_err, n_perms=10000)
    resp_anova = perm_rm_anova_2x2(resp_err, n_perms=10000)
    pairwise_abs = pairwise_tmax_by_fdr(abs_err, n_perms=10000)

    print("Absolute Error ANOVA:")
    print("  Cue: ", abs_anova["cue"])
    print("  Side:", abs_anova["side"])
    print("  Int: ", abs_anova["int"])
    print("  Pairwise tmax (BY-FDR):", pairwise_abs)

    print("Signed Error ANOVA:")
    print("  Cue: ", resp_anova["cue"])
    print("  Side:", resp_anova["side"])
    print("  Int: ", resp_anova["int"])

    fig2_path = FIGURES_DIR / "figure2_behavior.png"
    plot_publication_figure2(abs_err, resp_err, staircase_grays, pairwise_abs, fig2_path)

    right_floor = int(round(float(np.mean(abs_err[:, :, 1]))))

    # Push updates to pipeline API
    api = OrientExoPipelineAPI()
    update_payload = {
        "action": "update_stats",
        "source_script": "manuscript/scripts/run_behavior_pipeline.py",
        "updates": {
            "participants": {
                "n_recruited": 35,
                "n_excluded_eog": 7,
                "n_included": len(INCLUDED_PARTICIPANTS)
            },
            "task": {
                "staircase_gray_mean": round(float(np.mean(staircase_grays)), 1),
                "staircase_gray_sd": round(float(np.std(staircase_grays, ddof=1)), 1)
            },
            "behavior": {
                "df1": 1,
                "df2": len(INCLUDED_PARTICIPANTS) - 1,
                "right_side_floor_deg": right_floor,
                "pairwise_abs_error": pairwise_abs,
                "response_error": {
                    "target_side": {
                        "fmax": resp_anova["side"]["fmax"],
                        "crit_fmax": resp_anova["side"]["crit_fmax"],
                        "p": resp_anova["side"]["p"]
                    },
                    "cue_type": {
                        "fmax_str": "< 1" if resp_anova["cue"]["fmax"] < 1 else str(resp_anova["cue"]["fmax"]),
                        "p": "ns" if resp_anova["cue"]["p"] == "ns" or float(resp_anova["cue"]["p"]) > 0.05 else resp_anova["cue"]["p"]
                    },
                    "interaction": {
                        "fmax": resp_anova["int"]["fmax"],
                        "crit_fmax": resp_anova["int"]["crit_fmax"],
                        "p": resp_anova["int"]["p"]
                    }
                },
                "absolute_response_error": {
                    "target_side": {
                        "fmax": abs_anova["side"]["fmax"],
                        "crit_fmax": abs_anova["side"]["crit_fmax"],
                        "p": abs_anova["side"]["p"]
                    },
                    "cue_type": {
                        "fmax": abs_anova["cue"]["fmax"],
                        "crit_fmax": abs_anova["cue"]["crit_fmax"],
                        "p": abs_anova["cue"]["p"]
                    },
                    "interaction": {
                        "fmax": abs_anova["int"]["fmax"],
                        "crit_fmax": abs_anova["int"]["crit_fmax"],
                        "p": abs_anova["int"]["p"]
                    }
                }
            }
        }
    }
    api.handle_request(update_payload)

    # Register figure2 with provenance
    api.handle_request({
        "action": "register_figure",
        "source_script": "manuscript/scripts/run_behavior_pipeline.py",
        "figure_id": "figure2",
        "file": "figures/figure2_behavior.png",
        "status": "generated_from_analysis",
        "description": "Empirical Figure 2 generated from N=28 single-trial behavioral datasets (OSF u6fgm)"
    })

    print("[run_behavior_pipeline] Successfully updated behavior provenance and generated Figure 2!")

if __name__ == "__main__":
    main()

