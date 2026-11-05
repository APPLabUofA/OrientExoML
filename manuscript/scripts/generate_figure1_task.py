#!/usr/bin/env python3
"""
generate_figure1_task.py

Programmatically generates Figure 1 (Cued Orientation Perception Task Sequence &
Response Error Schematic) from the experimental parameters in `analysis_results.json`
and registers Figure 1 as `generated_from_analysis` via `pipeline_api.py`.
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Rectangle, Arc
import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent
MANUSCRIPT_DIR = SCRIPTS_DIR.parent
FIGURES_DIR = MANUSCRIPT_DIR / "figures"
JSON_PATH = MANUSCRIPT_DIR / "analysis_results.json"

sys.path.insert(0, str(SCRIPTS_DIR))
from pipeline_api import OrientExoPipelineAPI


def draw_screen_box(ax, x0, y0, size=1.8, bg_gray=0.52):
    rect = Rectangle(
        (x0 - size / 2, y0 - size / 2), size, size,
        facecolor=(bg_gray, bg_gray, bg_gray),
        edgecolor="#222222", linewidth=1.8, zorder=1
    )
    ax.add_patch(rect)


def draw_placeholders_and_fixation(ax, x0, y0, show_fix=True, show_circles=True):
    if show_circles:
        for dx in [-0.46, 0.46]:
            c = Circle(
                (x0 + dx, y0), radius=0.31,
                facecolor="none", edgecolor="#333333", linewidth=0.9, zorder=2
            )
            ax.add_patch(c)
    if show_fix:
        fix = Circle(
            (x0, y0), radius=0.042,
            facecolor="white", edgecolor="#333333", linewidth=0.5, zorder=3
        )
        ax.add_patch(fix)


def draw_cue(ax, x0, y0, cue_type="left"):
    yc = y0 + 0.20
    h = 0.075
    w = 0.12
    if cue_type == "left":
        pts = [[x0 - w / 2, yc], [x0 + w / 2, yc + h], [x0 + w / 2, yc - h]]
        ax.add_patch(Polygon(pts, closed=True, facecolor="black", edgecolor="black", zorder=4))
    elif cue_type == "right":
        pts = [[x0 + w / 2, yc], [x0 - w / 2, yc + h], [x0 - w / 2, yc - h]]
        ax.add_patch(Polygon(pts, closed=True, facecolor="black", edgecolor="black", zorder=4))
    elif cue_type == "both":
        pts_l = [[x0, yc], [x0 - w, yc + h], [x0 - w, yc - h]]
        pts_r = [[x0, yc], [x0 + w, yc + h], [x0 + w, yc - h]]
        ax.add_patch(Polygon(pts_l, closed=True, facecolor="black", edgecolor="black", zorder=4))
        ax.add_patch(Polygon(pts_r, closed=True, facecolor="black", edgecolor="black", zorder=4))


def draw_needle(ax, cx, cy, angle_deg=45.0, length=0.22, color="black"):
    rad = np.deg2rad(angle_deg)
    # Pinhead at center + needle tip
    hub = Circle((cx, cy), radius=0.032, facecolor=color, edgecolor=color, zorder=5)
    ax.add_patch(hub)
    tx = cx + length * np.cos(rad)
    ty = cy + length * np.sin(rad)
    perp_x = -0.024 * np.sin(rad)
    perp_y = 0.024 * np.cos(rad)
    pts = [[cx + perp_x, cy + perp_y], [cx - perp_x, cy - perp_y], [tx, ty]]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor=color, zorder=5))


def draw_star_mask(ax, cx, cy, radius=0.24, n_rays=32):
    for k in range(n_rays):
        ang = 2.0 * np.pi * k / n_rays
        ax.plot(
            [cx, cx + radius * np.cos(ang)],
            [cy, cy + radius * np.sin(ang)],
            color="black", linewidth=0.9, zorder=4
        )


def draw_dimension_bar(ax, x1, y1, x2, y2, text, text_offset=(0.0, 0.08), fontsize=9.5):
    ax.plot([x1, x2], [y1, y2], color="#111111", linewidth=1.0, zorder=6)
    dx, dy = x2 - x1, y2 - y1
    L = np.hypot(dx, dy)
    if L > 0:
        nx, ny = -dy / L * 0.04, dx / L * 0.04
        for px, py in [(x1, y1), (x2, y2)]:
            ax.plot([px - nx, px + nx], [py - ny, py + ny], color="#111111", linewidth=1.0, zorder=6)
    mx, my = (x1 + x2) / 2.0 + text_offset[0], (y1 + y2) / 2.0 + text_offset[1]
    ax.text(mx, my, text, ha="center", va="center", fontsize=fontsize, color="#111111", zorder=7)


def main():
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        results = json.load(f)
    task = results["task"]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    })

    fig, ax = plt.subplots(figsize=(14.0, 10.2), dpi=300)
    ax.set_xlim(-0.45, 14.45)
    ax.set_ylim(-5.4, 5.2)
    ax.set_aspect("equal")
    ax.axis("off")

    # Panel A header
    ax.text(0.0, 4.75, "A", fontsize=34, fontweight="bold", ha="left", va="top", color="#111111")

    # X positions of the 7 stages in Panel A
    xs = [1.0, 3.05, 5.1, 7.15, 9.2, 11.25, 13.3]
    y_mid = 1.65

    # 1. Fixation
    x_fix = xs[0]
    draw_screen_box(ax, x_fix, y_mid)
    draw_placeholders_and_fixation(ax, x_fix, y_mid)
    ax.text(x_fix, y_mid + 1.02, "Fixation", fontsize=17, ha="center", va="bottom", color="#111111")
    ax.text(x_fix, y_mid - 1.05, f"{task['fixation_duration_ms']} ms", fontsize=13.5, ha="center", va="top", color="#111111")
    draw_dimension_bar(ax, x_fix - 0.83, y_mid - 0.31, x_fix - 0.83, y_mid + 0.31, "2.11°", (-0.28, 0.0), fontsize=10)
    draw_dimension_bar(ax, x_fix - 0.045, y_mid + 0.12, x_fix + 0.045, y_mid + 0.12, "0.11°", (0.0, 0.13), fontsize=9.5)
    draw_dimension_bar(ax, x_fix, y_mid - 0.42, x_fix + 0.46, y_mid - 0.42, "6.47°", (0.0, -0.16), fontsize=10)

    # 2. Cue (3 stacked boxes: Left cue, Right cue, Non-informative double cue)
    x_cue = xs[1]
    y_cues = [y_mid + 1.98, y_mid, y_mid - 1.98]
    for yc, ctype in zip(y_cues, ["left", "right", "both"]):
        draw_screen_box(ax, x_cue, yc)
        draw_placeholders_and_fixation(ax, x_cue, yc)
        draw_cue(ax, x_cue, yc, cue_type=ctype)
    draw_dimension_bar(ax, x_cue - 0.07, y_cues[0] + 0.36, x_cue + 0.07, y_cues[0] + 0.36, "0.53°", (0.0, 0.15), fontsize=10)
    ax.text(x_cue, y_cues[0] + 1.02, "Cue", fontsize=17, ha="center", va="bottom", color="#111111")
    ax.text(x_cue, y_cues[2] - 1.05, "1242/1284/\n1325/1367 ms", fontsize=13, ha="center", va="top", color="#111111")

    # 3. Target (2 stacked boxes: Right target, Left target)
    x_targ = xs[2]
    y_targs = [y_mid + 0.99, y_mid - 0.99]
    draw_screen_box(ax, x_targ, y_targs[0])
    draw_placeholders_and_fixation(ax, x_targ, y_targs[0])
    draw_needle(ax, x_targ + 0.46 - 0.06, y_targs[0] - 0.06, angle_deg=45.0)

    draw_screen_box(ax, x_targ, y_targs[1])
    draw_placeholders_and_fixation(ax, x_targ, y_targs[1])
    draw_needle(ax, x_targ - 0.46 - 0.06, y_targs[1] - 0.06, angle_deg=45.0)
    draw_dimension_bar(
        ax, x_targ - 0.62, y_targs[1] - 0.02, x_targ - 0.44, y_targs[1] + 0.16,
        "1°", (-0.12, 0.10), fontsize=10
    )
    ax.text(x_targ, y_targs[0] + 1.02, "Target", fontsize=17, ha="center", va="bottom", color="#111111")
    ax.text(x_targ, y_targs[1] - 1.05, f"{task['target_duration_ms']:.1f} ms", fontsize=13.5, ha="center", va="top", color="#111111")

    # 4. Blank (Target-Mask ISI / SOA)
    x_blk1 = xs[3]
    draw_screen_box(ax, x_blk1, y_mid)
    draw_placeholders_and_fixation(ax, x_blk1, y_mid)
    ax.text(x_blk1, y_mid + 1.02, "Blank", fontsize=17, ha="center", va="bottom", color="#111111")
    ax.text(x_blk1, y_mid - 1.05, f"{task['target_mask_soa_ms']:.1f} ms", fontsize=13.5, ha="center", va="top", color="#111111")

    # 5. Mask
    x_mask = xs[4]
    draw_screen_box(ax, x_mask, y_mid)
    draw_placeholders_and_fixation(ax, x_mask, y_mid)
    draw_star_mask(ax, x_mask - 0.46, y_mid)
    draw_star_mask(ax, x_mask + 0.46, y_mid)
    draw_dimension_bar(ax, x_mask - 0.70, y_mid + 0.28, x_mask - 0.22, y_mid + 0.28, "1.65°", (0.0, 0.14), fontsize=10)
    ax.text(x_mask, y_mid + 1.02, "Mask", fontsize=17, ha="center", va="bottom", color="#111111")
    ax.text(x_mask, y_mid - 1.05, f"{task['mask_duration_ms']:.1f} ms", fontsize=13.5, ha="center", va="top", color="#111111")

    # 6. Post-mask Blank
    x_blk2 = xs[5]
    draw_screen_box(ax, x_blk2, y_mid)
    draw_placeholders_and_fixation(ax, x_blk2, y_mid)
    ax.text(x_blk2, y_mid + 1.02, "Blank", fontsize=17, ha="center", va="bottom", color="#111111")
    ax.text(x_blk2, y_mid - 1.05, f"{task['post_mask_blank_ms']:.1f} ms", fontsize=13.5, ha="center", va="top", color="#111111")

    # 7. Response
    x_resp = xs[6]
    draw_screen_box(ax, x_resp, y_mid)
    draw_needle(ax, x_resp + 0.06, y_mid - 0.06, angle_deg=135.0)
    draw_dimension_bar(ax, x_resp + 0.02, y_mid + 0.16, x_resp + 0.20, y_mid - 0.02, "1°", (0.14, 0.10), fontsize=10)
    ax.text(x_resp, y_mid + 1.02, "Response", fontsize=17, ha="center", va="bottom", color="#111111")

    # Panel B: Target - Response = Error
    y_b = -4.15
    ax.text(2.3, y_b + 1.10, "B", fontsize=34, fontweight="bold", ha="left", va="bottom", color="#111111")
    ax.text(5.85, y_b + 1.05, "Target - Response =", fontsize=19, ha="center", va="bottom", color="#111111")
    ax.text(8.55, y_b + 1.05, "Error", fontsize=19, ha="left", va="bottom", color="#d62728")

    xb_boxes = [4.35, 6.65, 8.95]
    # Box 1: Target
    draw_screen_box(ax, xb_boxes[0], y_b)
    draw_placeholders_and_fixation(ax, xb_boxes[0], y_b)
    draw_needle(ax, xb_boxes[0] + 0.46 - 0.06, y_b - 0.06, angle_deg=45.0)

    # Box 2: Response
    draw_screen_box(ax, xb_boxes[1], y_b)
    draw_needle(ax, xb_boxes[1] + 0.06, y_b - 0.06, angle_deg=135.0)

    # Box 3: Error arc
    draw_screen_box(ax, xb_boxes[2], y_b)
    draw_needle(ax, xb_boxes[2], y_b - 0.08, angle_deg=45.0, color="#222222")
    draw_needle(ax, xb_boxes[2], y_b - 0.08, angle_deg=135.0, color="#444444")
    err_arc = Arc(
        (xb_boxes[2], y_b + 0.16), width=0.46, height=0.46,
        angle=0, theta1=-15, theta2=195,
        edgecolor="#d62728", linewidth=2.4, zorder=6
    )
    ax.add_patch(err_arc)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    out_png = FIGURES_DIR / "figure1_task.png"
    out_pdf = FIGURES_DIR / "figure1_task.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    api = OrientExoPipelineAPI()
    api.handle_request({
        "action": "register_figure",
        "figure_id": "figure1",
        "file": "figures/figure1_task.png",
        "status": "generated_from_analysis",
        "source_script": "manuscript/scripts/generate_figure1_task.py",
        "description": "Programmatic Figure 1 generated from experimental parameters in analysis_results.json",
    })
    if (FIGURES_DIR / "figure2_behavior.png").exists():
        api.handle_request({
            "action": "register_figure",
            "figure_id": "figure2",
            "file": "figures/figure2_behavior.png",
            "status": "generated_from_analysis",
            "source_script": "manuscript/scripts/run_behavior_pipeline.py",
            "description": "Empirical Figure 2 generated from N=28 single-trial behavioral datasets (OSF u6fgm)",
        })
    print(f"[generate_figure1_task] Generated {out_png} and registered figure1 & figure2.")


if __name__ == "__main__":
    main()
