#!/usr/bin/env python3
"""
regenerate_mvpa_figures.py

Template & API Generator for Figures 1–7 of OrientExoML (CerCor-2022-00354).
Can be called directly or via `pipeline_api.py` when updated MATLAB/Python
MVPA outputs (.mat / .json / .csv) are produced by the analysis scripts.

Addresses Reviewer 2's formatting requirements for Figures 3–6:
- High resolution (300 DPI PNG + vector PDF)
- Time axis ticks every 200 ms (-1000, -800, -600, -400, -200, 0 ms) in larger font (11 pt)
- Concise, readable abbreviations explained in figure legends:
  * "Accuracy (%)" instead of "Classification Accuracy (%)"
  * "SVR r_z" instead of "Fisher-Z Transformed Correlation Coefficient"
  * "AAC", "Alpha", "MI", "MVL_Özkurt", "PAC_GLM", "PAC_PLV"
  * "Corr", "wPLI", "PLV", "iCOH"
- Clean scientific figure styling (no boxed statistical text overlays, despined axes,
  non-overlapping legends).
"""

import argparse
import json
from pathlib import Path
from typing import Dict, Any, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent
MANUSCRIPT_DIR = SCRIPTS_DIR.parent
FIGURES_DIR = MANUSCRIPT_DIR / "figures"
JSON_PATH = MANUSCRIPT_DIR / "analysis_results.json"


def apply_publication_style():
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


def plot_mvpa_timecourse_grid(
    time_ms: np.ndarray,
    empirical_curves: Dict[str, np.ndarray],
    empirical_sems: Dict[str, np.ndarray],
    null_curves: Dict[str, np.ndarray],
    null_sems: Dict[str, np.ndarray],
    sig_masks: Dict[str, np.ndarray],
    row_labels: list,
    col_labels: list,
    ylabel_short: str,
    chance_val: float,
    ylim: tuple,
    out_path: Path,
):
    """
    Renders a high-readability multi-row MVPA time-course figure conforming to
    Reviewer 2's specifications (200 ms tick spacing, abbreviated labels, large font).
    """
    apply_publication_style()
    n_rows = len(row_labels)
    n_cols = len(col_labels)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(4.2 * n_cols, 2.1 * n_rows),
        sharex=True, sharey=True
    )
    if n_rows == 1:
        axes = np.expand_dims(axes, axis=0)
    if n_cols == 1:
        axes = np.expand_dims(axes, axis=1)

    xticks_200 = [-1000, -800, -600, -400, -200, 0]

    for r, r_label in enumerate(row_labels):
        for c, c_label in enumerate(col_labels):
            ax = axes[r, c]
            key = f"{r_label}__{c_label}"
            emp = empirical_curves[key]
            emp_sem = empirical_sems[key]
            nul = null_curves[key]
            nul_sem = null_sems[key]
            sig = sig_masks[key]

            # Shade significant time intervals in green
            for idx in np.where(sig)[0]:
                dt = (time_ms[1] - time_ms[0]) / 2.0 if len(time_ms) > 1 else 10.0
                ax.axvspan(time_ms[idx] - dt, time_ms[idx] + dt, color="#7fcdbb", alpha=0.55, lw=0)

            ax.axvline(0, color="#888888", linestyle="--", linewidth=1.0)
            ax.errorbar(time_ms, nul, yerr=nul_sem, color="#2b5c8f", lw=1.2, capsize=1.5, label="Chance")
            ax.errorbar(time_ms, emp, yerr=emp_sem, color="#111111", lw=1.4, capsize=1.5, label="Actual")

            ax.set_xlim(-1000, 10)
            ax.set_ylim(*ylim)
            ax.set_xticks(xticks_200)

            if r == 0:
                ax.set_title(c_label, fontweight="bold", pad=8)
            if r == n_rows - 1:
                ax.set_xlabel("Time (ms)")
            if c == 0:
                ax.set_ylabel(ylabel_short)
            if c == n_cols - 1:
                ax.text(
                    1.04, 0.5, r_label,
                    transform=ax.transAxes,
                    rotation=-90, va="center", ha="left",
                    fontweight="bold", fontsize=11.5
                )
            if r == 0 and c == 0:
                ax.legend(loc="upper right", frameon=False, ncol=2)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[regenerate_mvpa_figures] Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Regenerate high-readability MVPA figures from data payload")
    parser.add_argument("--data-json", type=str, help="Path to JSON file containing time-series arrays from MATLAB/Python")
    args = parser.parse_args()

    if args.data_json:
        with open(args.data_json, "r", encoding="utf-8") as f:
            payload = json.load(f)
        time_ms = np.array(payload["time_ms"])
        emp = {k: np.array(v) for k, v in payload["empirical"].items()}
        emp_sem = {k: np.array(v) for k, v in payload["empirical_sem"].items()}
        nul = {k: np.array(v) for k, v in payload["null"].items()}
        nul_sem = {k: np.array(v) for k, v in payload["null_sem"].items()}
        sig = {k: np.array(v, dtype=bool) for k, v in payload["significant"].items()}
        out_file = FIGURES_DIR / payload.get("out_filename", "figure_custom_mvpa.png")
        plot_mvpa_timecourse_grid(
            time_ms=time_ms,
            empirical_curves=emp,
            empirical_sems=emp_sem,
            null_curves=nul,
            null_sems=nul_sem,
            sig_masks=sig,
            row_labels=payload["row_labels"],
            col_labels=payload["col_labels"],
            ylabel_short=payload.get("ylabel_short", "Accuracy (%)"),
            chance_val=payload.get("chance_val", 50.0),
            ylim=tuple(payload.get("ylim", [35, 70])),
            out_path=out_file,
        )
    else:
        print("Usage: python3 regenerate_mvpa_figures.py --data-json <path_to_timeseries.json>")


if __name__ == "__main__":
    main()

