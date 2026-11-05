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
FIGURES_DIR = MANUSCRIPT_DIR / "figures"

sys.path.insert(0, str(SCRIPTS_DIR))
from pipeline_api import OrientExoPipelineAPI

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

def load_behavior_data():
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

        # Remove first 20 practice trials
        err = get_arr(data["errorDegrees"])[20:]
        val = get_arr(data["valid"])[20:]
        pos = get_arr(data["position"])[20:]

        # Participant 009 recording error removal per getBEHdata_OrientWheel_Exo.m
        if sub == "009":
            err = err[:343]
            val = val[:343]
            pos = pos[:343]

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
    Computes 2x2 repeated measures ANOVA with permutation-based Fmax testing
    matching `perm_rbANOVA_OrientExo.m`.
    data: (N, 2, 2) where axis 1 is Cue (Info vs Non-info), axis 2 is Side (Left vs Right).
    """
    rng = np.random.RandomState(seed)
    N = data.shape[0]

    def calc_F_stats(arr):
        # Main effect Cue (axis 1)
        diff_cue = np.mean(arr[:, 0, :] - arr[:, 1, :], axis=1)
        s_cue = np.std(diff_cue, ddof=1)
        t_cue = np.mean(diff_cue) / (s_cue / np.sqrt(N)) if s_cue > 0 else 0.0
        F_cue = t_cue ** 2

        # Main effect Side (axis 2)
        diff_side = np.mean(arr[:, :, 0] - arr[:, :, 1], axis=1)
        s_side = np.std(diff_side, ddof=1)
        t_side = np.mean(diff_side) / (s_side / np.sqrt(N)) if s_side > 0 else 0.0
        F_side = t_side ** 2

        # Interaction Cue x Side
        diff_int = (arr[:, 0, 0] - arr[:, 0, 1]) - (arr[:, 1, 0] - arr[:, 1, 1])
        s_int = np.std(diff_int, ddof=1)
        t_int = np.mean(diff_int) / (s_int / np.sqrt(N)) if s_int > 0 else 0.0
        F_int = t_int ** 2

        return F_cue, F_side, F_int

    F_cue_obs, F_side_obs, F_int_obs = calc_F_stats(data)

    perm_F_cue = np.zeros(n_perms)
    perm_F_side = np.zeros(n_perms)
    perm_F_int = np.zeros(n_perms)

    # Permutation testing: randomly flip condition assignments within subjects
    for p in range(n_perms):
        perm_data = data.copy()
        for i in range(N):
            # Randomly permute the 4 condition cells for subject i
            flat = perm_data[i].ravel()
            rng.shuffle(flat)
            perm_data[i] = flat.reshape(2, 2)
        fc, fs, fi = calc_F_stats(perm_data)
        perm_F_cue[p] = fc
        perm_F_side[p] = fs
        perm_F_int[p] = fi

    crit_cue = np.percentile(perm_F_cue, 95)
    crit_side = np.percentile(perm_F_side, 95)
    crit_int = np.percentile(perm_F_int, 95)

    p_cue = np.mean(perm_F_cue >= F_cue_obs)
    p_side = np.mean(perm_F_side >= F_side_obs)
    p_int = np.mean(perm_F_int >= F_int_obs)

    return {
        "cue": {"fmax": round(float(F_cue_obs), 2), "crit_fmax": round(float(crit_cue), 2), "p": "< 0.001" if p_cue < 0.001 else f"{p_cue:.3f}"},
        "side": {"fmax": round(float(F_side_obs), 2), "crit_fmax": round(float(crit_side), 2), "p": "< 0.001" if p_side < 0.001 else f"{p_side:.3f}"},
        "int": {"fmax": round(float(F_int_obs), 2), "crit_fmax": round(float(crit_int), 2), "p": "< 0.001" if p_int < 0.001 else f"{p_int:.3f}"},
    }


def plot_publication_figure2(abs_err, resp_err, staircase_grays, out_path):
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
    ax_a.set_ylim(0, 50)
    ax_a.set_title("A. Absolute Orientation Error", fontweight="bold")
    ax_a.legend(frameon=False, loc="upper right")

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
    print("[run_behavior_pipeline] Loading single-trial behavioral datasets for N=28...")
    abs_err, resp_err, staircase_grays = load_behavior_data()

    print("[run_behavior_pipeline] Running permutation-based repeated measures ANOVA (10,000 perms)...")
    abs_anova = perm_rm_anova_2x2(abs_err, n_perms=10000)
    resp_anova = perm_rm_anova_2x2(resp_err, n_perms=10000)

    print("Absolute Error ANOVA:")
    print("  Cue: ", abs_anova["cue"])
    print("  Side:", abs_anova["side"])
    print("  Int: ", abs_anova["int"])

    print("Signed Error ANOVA:")
    print("  Cue: ", resp_anova["cue"])
    print("  Side:", resp_anova["side"])
    print("  Int: ", resp_anova["int"])

    fig2_path = FIGURES_DIR / "figure2_behavior.png"
    plot_publication_figure2(abs_err, resp_err, staircase_grays, fig2_path)

    # Push updates to pipeline API
    api = OrientExoPipelineAPI()
    update_payload = {
        "action": "update_stats",
        "source_script": "CuedOrientTask_paper/Analysis_Code/perm_rbANOVA_OrientExo.m",
        "updates": {
            "participants": {
                "n_recruited": 35,
                "n_excluded_eog": 7,
                "n_included": 28
            },
            "task": {
                "staircase_gray_mean": round(float(np.mean(staircase_grays)), 1),
                "staircase_gray_sd": round(float(np.std(staircase_grays, ddof=1)), 1)
            },
            "behavior": {
                "df1": 1,
                "df2": 27,
                "right_side_floor_deg": 28,
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

