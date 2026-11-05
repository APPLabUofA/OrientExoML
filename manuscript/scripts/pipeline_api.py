#!/usr/bin/env python3
r"""
pipeline_api.py

End-to-End Reproducible Manuscript & Figure Generation API for OrientExoML
(CerCor-2022-00354: "Connecting Covert Attention and Visual Perception to the
Spatiotemporal Dynamics of Alpha Band Activity, Cross-Frequency Coupling (CFC),
and Functional Connectivity using Multivariate Pattern Analysis (MVPA)").

Anti-Hotwiring & Provenance Features:
1. ZERO hardcoded stats in figure generators: Figure 7 reads its arrays directly
   from `analysis_results.json["mvpa_summary_metrics"]`.
2. Provenance tracking: Every statistical block and figure entry tracks
   `status` ("placeholder_*" vs "computed_from_raw_data" / "generated_from_analysis"),
   `source_script`, `updated_at`, and `sha256` file checksums.
3. Strict Submission Audit (`--audit` / `--strict`):
   - Scans `main.tex`, `cover_letter.tex`, and `response_to_reviewers.tex` for
     hotwired statistical expressions (e.g., literal `Fmax = ...`, `p < 0.0...`,
     or hardcoded `\includegraphics{figures/...}` bypassing `\figPath*` macros).
   - Reports which stats and figures are still placeholders vs. verified analysis outputs.
   - When `--strict` is passed, blocks final submission compilation if any
     placeholder stats, placeholder screenshots, or hotwired literals remain.
"""

import argparse
import datetime
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent
MANUSCRIPT_DIR = SCRIPTS_DIR.parent
PROJECT_ROOT = MANUSCRIPT_DIR.parent
JSON_PATH = MANUSCRIPT_DIR / "analysis_results.json"
FIGURES_DIR = MANUSCRIPT_DIR / "figures"
BIORXIV_PDF = PROJECT_ROOT / "recerebralcortexmanuscriptdecisioncercor20220035" / "BIORXIV-2022-493754v1-Sheldon.pdf"


def compute_sha256(file_path: Path) -> str:
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


class OrientExoPipelineAPI:
    """
    Programmatic API for analysis scripts (MATLAB / Python) and figure generators
    to update statistics, regenerate figures, audit provenance, and compile the
    LaTeX manuscript, cover letter, and reviewer response documents.
    """

    def __init__(self, json_path: Path = JSON_PATH):
        self.json_path = json_path
        self.data = self._load()

    def _load(self) -> Dict[str, Any]:
        if not self.json_path.exists():
            raise FileNotFoundError(f"Missing analysis results file: {self.json_path}")
        with open(self.json_path, "r", encoding="utf-8") as f:
            return json.load(f)

    def save(self) -> None:
        with open(self.json_path, "w", encoding="utf-8") as f:
            json.dump(self.data, f, indent=2)
        subprocess.run(["python3", str(SCRIPTS_DIR / "generate_latex_stats.py")], check=True)

    def handle_request(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a JSON API request from an analysis script or figure generator.

        Supported actions:
        - "update_stats": Merge nested dictionary keys into `analysis_results.json`,
          automatically record provenance (`status: "computed_from_raw_data"`,
          `source_script`, `updated_at`), and regenerate `generated_stats.tex`.
        - "register_figure": Update a figure's file path, compute its SHA256 checksum,
          and mark its status as `"generated_from_analysis"`.
        - "generate_summary_figure": Generate Figure 7 from `analysis_results.json["mvpa_summary_metrics"]`.
        - "extract_placeholders": Extract 300 DPI screenshots from existing bioRxiv PDF for Figures 1-6.
        - "audit": Run anti-hotwiring and provenance audit across JSON and LaTeX files.
        - "compile_all": Recompile `main.tex`, `cover_letter.tex`, and `response_to_reviewers.tex`.
        """
        action = payload.get("action", "")
        source_script = payload.get("source_script", "external_analysis_api")
        now_iso = datetime.datetime.now(datetime.timezone.utc).isoformat()

        if action == "update_stats":
            updates = payload.get("updates", {})
            self._deep_update(self.data, updates)
            # Mark updated blocks as computed_from_raw_data in provenance
            prov = self.data.setdefault("provenance", {})
            for section_key in updates.keys():
                if section_key in prov:
                    prov[section_key] = {
                        "status": payload.get("status", "computed_from_raw_data"),
                        "source_script": source_script,
                        "updated_at": now_iso,
                    }
            # If mvpa_summary_metrics was updated, automatically re-render Figure 7
            if "mvpa_summary_metrics" in updates:
                self.generate_figure7_mvpa_summary(
                    status=payload.get("status", "generated_from_analysis"),
                    source_script=source_script,
                )
            self._refresh_strict_readiness()
            self.save()
            return {
                "status": "ok",
                "action": action,
                "updated_sections": list(updates.keys()),
                "strict_mode_ready": self.data["metadata"]["strict_mode_ready"],
            }

        elif action == "register_figure":
            fig_id = payload["figure_id"]
            if fig_id not in self.data["figures"]:
                self.data["figures"][fig_id] = {"id": fig_id}
            rel_file = payload.get("file", self.data["figures"][fig_id].get("file"))
            abs_file = MANUSCRIPT_DIR / rel_file
            if not abs_file.exists():
                raise FileNotFoundError(f"Registered figure file does not exist: {abs_file}")
            sha = compute_sha256(abs_file)
            self.data["figures"][fig_id].update({
                "file": rel_file,
                "status": payload.get("status", "generated_from_analysis"),
                "source_script": source_script,
                "updated_at": now_iso,
                "sha256": sha,
                "description": payload.get("description", self.data["figures"][fig_id].get("description", "")),
            })
            self._refresh_strict_readiness()
            self.save()
            return {"status": "ok", "action": action, "figure_id": fig_id, "sha256": sha}

        elif action == "generate_summary_figure":
            out_path = self.generate_figure7_mvpa_summary(
                status=payload.get("status", self.data["provenance"]["mvpa_summary_metrics"]["status"]),
                source_script=source_script,
            )
            self.save()
            return {"status": "ok", "action": action, "file": str(out_path)}

        elif action == "extract_placeholders":
            extracted = self.extract_placeholder_screenshots()
            self.save()
            return {"status": "ok", "action": action, "files": extracted}

        elif action == "audit":
            report = self.audit_pipeline()
            return report

        elif action == "compile_all":
            strict = bool(payload.get("strict", False))
            self.save()
            pdfs = self.compile_latex_documents(strict=strict)
            return {"status": "ok", "action": action, "pdfs": pdfs}

        else:
            raise ValueError(f"Unknown API action: {action}")

    @staticmethod
    def _deep_update(target: Dict[str, Any], source: Dict[str, Any]) -> None:
        for k, v in source.items():
            if isinstance(v, dict) and isinstance(target.get(k), dict):
                OrientExoPipelineAPI._deep_update(target[k], v)
            else:
                target[k] = v

    def _refresh_strict_readiness(self) -> None:
        prov_ok = all(
            not str(v.get("status", "")).startswith("placeholder")
            for v in self.data.get("provenance", {}).values()
        )
        figs_ok = all(
            not str(v.get("status", "")).startswith("placeholder")
            for v in self.data.get("figures", {}).values()
        )
        self.data["metadata"]["strict_mode_ready"] = bool(prov_ok and figs_ok)

    def extract_placeholder_screenshots(self) -> List[str]:
        """
        Extracts 300 DPI placeholder screenshots from the existing bioRxiv PDF
        (Pages 50-55) into `manuscript/figures/` and records their SHA256 hashes
        and `placeholder_screenshot` status in `analysis_results.json`.
        """
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        mapping = {
            "figure1": (51, "figure1_task.png"),
            "figure2": (52, "figure2_behavior.png"),
            "figure3": (53, "figure3_svm_cfc.png"),
            "figure4": (54, "figure4_svm_fc.png"),
            "figure5": (55, "figure5_svr_cfc.png"),
            "figure6": (50, "figure6_svr_fc.png"),
        }
        extracted_files = []
        for fig_id, (page_num, fname) in mapping.items():
            out_file = FIGURES_DIR / fname
            existing_status = self.data.get("figures", {}).get(fig_id, {}).get("status", "")
            if existing_status == "generated_from_analysis" and out_file.exists():
                self.data["figures"][fig_id]["sha256"] = compute_sha256(out_file)
                extracted_files.append(str(out_file))
                continue
            if not out_file.exists() and BIORXIV_PDF.exists():
                cmd = [
                    "/opt/homebrew/bin/gs",
                    "-dNOPAUSE", "-dBATCH",
                    f"-dFirstPage={page_num}", f"-dLastPage={page_num}",
                    "-sDEVICE=png16m", "-r300",
                    f"-sOutputFile={out_file}",
                    str(BIORXIV_PDF)
                ]
                subprocess.run(cmd, check=True)
            if out_file.exists():
                self.data["figures"][fig_id]["sha256"] = compute_sha256(out_file)
                self.data["figures"][fig_id]["status"] = "placeholder_screenshot"
            extracted_files.append(str(out_file))
        return extracted_files

    def generate_figure7_mvpa_summary(
        self,
        status: str = "placeholder_from_manuscript",
        source_script: str = "pipeline_api.py::generate_figure7_mvpa_summary",
    ) -> Path:
        """
        Generates Figure 7 dynamically from `analysis_results.json["mvpa_summary_metrics"]`
        with ZERO hardcoded numerical arrays in Python code.
        """
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        out_png = FIGURES_DIR / "figure7_mvpa_summary.png"
        out_pdf = FIGURES_DIR / "figure7_mvpa_summary.pdf"

        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9.5,
        })

        sm = self.data["mvpa_summary_metrics"]
        metrics = sm["labels"]
        y_pos = np.arange(len(metrics))[::-1]

        svm_left_peak = np.array(sm["svm_left_peak"], dtype=float)
        svm_left_sem  = np.array(sm["svm_left_sem"], dtype=float)
        svm_right_peak = np.array(sm["svm_right_peak"], dtype=float)
        svm_right_sem  = np.array(sm["svm_right_sem"], dtype=float)

        svm_left_pct_sig  = np.array(sm["svm_left_pct_sig"], dtype=float)
        svm_right_pct_sig = np.array(sm["svm_right_pct_sig"], dtype=float)

        svr_noninf_peak = np.array(sm["svr_noninf_peak"], dtype=float)
        svr_noninf_sem  = np.array(sm["svr_noninf_sem"], dtype=float)
        svr_inf_peak    = np.array(sm["svr_inf_peak"], dtype=float)
        svr_inf_sem     = np.array(sm["svr_inf_sem"], dtype=float)

        fig, axes = plt.subplots(
            1, 3, figsize=(13.5, 5.2), sharey=True,
            gridspec_kw={"width_ratios": [1.25, 1.0, 1.25], "wspace": 0.14}
        )

        bar_h = 0.34
        c_left = "#1f77b4"
        c_right = "#e05a47"
        c_noninf = "#d94882"
        c_inf = "#1f77b4"

        # Panel A: Peak SVM Classification Accuracy (%)
        ax = axes[0]
        ax.axvline(50.0, color="#666666", linestyle="--", linewidth=1.0, zorder=1)
        ax.barh(y_pos + bar_h / 2, svm_left_peak, height=bar_h, xerr=svm_left_sem,
                color=c_left, alpha=0.9, capsize=2.5, label="Left Targets", zorder=3)
        ax.barh(y_pos - bar_h / 2, svm_right_peak, height=bar_h, xerr=svm_right_sem,
                color=c_right, alpha=0.9, capsize=2.5, label="Right Targets", zorder=3)

        for idx in sm.get("svm_left_sig_indices", []):
            ax.text(svm_left_peak[idx] + svm_left_sem[idx] + 0.35, y_pos[idx] + bar_h / 2,
                    "*", va="center", ha="left", fontsize=12, fontweight="bold", color=c_left)
        for idx in sm.get("svm_right_sig_indices", []):
            ax.text(svm_right_peak[idx] + svm_right_sem[idx] + 0.35, y_pos[idx] - bar_h / 2,
                    "*", va="center", ha="left", fontsize=12, fontweight="bold", color=c_right)

        svm_max_val = float(np.max(np.maximum(svm_left_peak + svm_left_sem, svm_right_peak + svm_right_sem)))
        ax.set_xlim(45, max(68.0, svm_max_val + 8.5))
        ax.set_yticks(y_pos)
        ax.set_yticklabels(metrics, fontweight="bold")
        ax.set_xlabel("Peak Classification Accuracy (%)")
        ax.set_title("SVM Cue Classification Peak", fontweight="bold", pad=10)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(loc="lower right", frameon=False)
        ax.text(-0.28, 1.06, "A", transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")

        # Panel B: % Significant Pre-Stimulus Time Windows (SVM)
        ax2 = axes[1]
        ax2.barh(y_pos + bar_h / 2, svm_left_pct_sig, height=bar_h,
                 color=c_left, alpha=0.9, label="Left Targets", zorder=3)
        ax2.barh(y_pos - bar_h / 2, svm_right_pct_sig, height=bar_h,
                 color=c_right, alpha=0.9, label="Right Targets", zorder=3)
        ax2.set_xlim(0, 100)
        ax2.set_xlabel("Significant Time Windows (%)")
        ax2.set_title("SVM Temporal Extent (p < 0.05)", fontweight="bold", pad=10)
        ax2.spines[["top", "right"]].set_visible(False)
        ax2.text(-0.08, 1.06, "B", transform=ax2.transAxes, fontsize=16, fontweight="bold", va="top")

        # Panel C: Peak SVR Prediction Accuracy (Fisher-Z r)
        ax3 = axes[2]
        ax3.axvline(0.0, color="#666666", linestyle="--", linewidth=1.0, zorder=1)
        ax3.barh(y_pos + bar_h / 2, svr_noninf_peak, height=bar_h, xerr=svr_noninf_sem,
                 color=c_noninf, alpha=0.9, capsize=2.5, label="Non-Informative Cue", zorder=3)
        ax3.barh(y_pos - bar_h / 2, svr_inf_peak, height=bar_h, xerr=svr_inf_sem,
                 color=c_inf, alpha=0.9, capsize=2.5, label="Informative Cue", zorder=3)

        for idx in sm.get("svr_noninf_sig_indices", []):
            ax3.text(svr_noninf_peak[idx] + svr_noninf_sem[idx] + 0.004, y_pos[idx] + bar_h / 2,
                     "*", va="center", ha="left", fontsize=12, fontweight="bold", color=c_noninf)
        for idx in sm.get("svr_inf_sig_indices", []):
            ax3.text(svr_inf_peak[idx] + svr_inf_sem[idx] + 0.004, y_pos[idx] - bar_h / 2,
                     "*", va="center", ha="left", fontsize=12, fontweight="bold", color=c_inf)

        svr_max_val = float(np.max(np.maximum(svr_noninf_peak + svr_noninf_sem, svr_inf_peak + svr_inf_sem)))
        ax3.set_xlim(-0.02, max(0.24, svr_max_val + 0.09))
        ax3.set_xlabel("Peak Fisher-Z Correlation ($r_z$)")
        ax3.set_title("SVR Error Prediction Peak", fontweight="bold", pad=10)
        ax3.spines[["top", "right"]].set_visible(False)
        ax3.legend(loc="lower right", frameon=False)
        ax3.text(-0.08, 1.06, "C", transform=ax3.transAxes, fontsize=16, fontweight="bold", va="top")

        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)

        self.data["figures"]["figure7"].update({
            "file": "figures/figure7_mvpa_summary.png",
            "status": status,
            "source_script": source_script,
            "sha256": compute_sha256(out_png),
        })
        print(f"[pipeline_api] Generated summary figure: {out_png} (status={status})")
        return out_png

    def stamp_figure2_panels(self) -> Path:
        """
        Addresses Reviewer 2's comment: 'Figure 2 should be divided into A and B panel'
        on the placeholder screenshot if not already stamped.
        """
        fig2_path = FIGURES_DIR / "figure2_behavior.png"
        if not fig2_path.exists():
            return fig2_path
        img = plt.imread(str(fig2_path))
        h, w, c = img.shape
        # Avoid double-padding if already stamped
        if h > 440:
            self.data["figures"]["figure2"]["sha256"] = compute_sha256(fig2_path)
            return fig2_path
        pad_top, pad_left = 45, 30
        padded = np.ones((h + pad_top, w + pad_left, c), dtype=img.dtype)
        padded[pad_top:, pad_left:, :] = img

        ph, pw = padded.shape[:2]
        fig, ax = plt.subplots(figsize=(pw / 300.0, ph / 300.0), dpi=300)
        ax.imshow(padded)
        ax.axis("off")
        ax.text(0.01, 0.98, "A", transform=ax.transAxes, fontsize=15, fontweight="bold", va="top", ha="left")
        ax.text(0.49, 0.98, "B", transform=ax.transAxes, fontsize=15, fontweight="bold", va="top", ha="left")
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig.savefig(fig2_path, dpi=300, bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)
        self.data["figures"]["figure2"]["sha256"] = compute_sha256(fig2_path)
        print(f"[pipeline_api] Stamped Panel A & B labels onto {fig2_path}")
        return fig2_path

    def audit_pipeline(self) -> Dict[str, Any]:
        r"""
        Anti-Hotwiring & Provenance Audit:
        1. Checks `main.tex`, `cover_letter.tex`, and `response_to_reviewers.tex`
           for any hardcoded figure paths (`\includegraphics{...}` not using `\figPath*`
           except cover-letter logos) or hardcoded statistical test literals (`F_{\max}(...) = [0-9]`).
        2. Checks `provenance` and `figures` in `analysis_results.json` to list
           which components are still `placeholder_*` vs `computed_from_raw_data` / `generated_from_analysis`.
        3. Verifies SHA256 checksums of all figure files on disk.
        """
        hotwired_violations: List[str] = []
        tex_files = ["main.tex", "cover_letter.tex", "response_to_reviewers.tex"]

        # Patterns that must ALWAYS use macros from generated_stats.tex, never raw numbers:
        forbidden_patterns = [
            (re.compile(r"F_\{\\max\}\([^)]*\)\s*=\s*[0-9]"), "Hardcoded F_max value (must use \\stat*Fmax macro)"),
            (re.compile(r"critical\s+F_\{\\max\}\s*=\s*\\pm\s*[0-9]"), "Hardcoded critical F_max value (must use \\stat*Crit macro)"),
            (re.compile(r"N\s*=\s*[0-9]+"), "Hardcoded sample size N (must use \\statNIncluded or \\statNRecruited)"),
            (re.compile(r"M\s*=\s*[0-9]+\.?[0-9]*"), "Hardcoded mean M (must use \\stat*Mean* macro)"),
            (re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{(?!\\figPath|figures/coverletter_)([^}]+)\}"), "Hardcoded figure path in \\includegraphics (must use \\figPath* macro)"),
        ]

        for tname in tex_files:
            tpath = MANUSCRIPT_DIR / tname
            if not tpath.exists():
                continue
            lines = tpath.read_text(encoding="utf-8").splitlines()
            for lineno, line in enumerate(lines, start=1):
                if line.lstrip().startswith("%"):
                    continue
                for regex, msg in forbidden_patterns:
                    if regex.search(line):
                        hotwired_violations.append(f"{tname}:{lineno}: {msg} -> {line.strip()}")

        placeholder_stats = [
            k for k, v in self.data.get("provenance", {}).items()
            if str(v.get("status", "")).startswith("placeholder")
        ]
        placeholder_figures = [
            k for k, v in self.data.get("figures", {}).items()
            if str(v.get("status", "")).startswith("placeholder")
        ]

        checksum_mismatches = []
        for fig_id, fmeta in self.data.get("figures", {}).items():
            fpath = MANUSCRIPT_DIR / fmeta["file"]
            if not fpath.exists():
                checksum_mismatches.append(f"{fig_id}: file missing ({fpath})")
            elif fmeta.get("sha256"):
                actual_sha = compute_sha256(fpath)
                if actual_sha != fmeta["sha256"]:
                    checksum_mismatches.append(
                        f"{fig_id}: SHA256 mismatch (disk file modified outside pipeline API)"
                    )

        self._refresh_strict_readiness()
        return {
            "hotwired_tex_violations": hotwired_violations,
            "placeholder_stat_blocks": placeholder_stats,
            "placeholder_figures": placeholder_figures,
            "figure_checksum_mismatches": checksum_mismatches,
            "strict_submission_ready": (
                len(hotwired_violations) == 0
                and len(placeholder_stats) == 0
                and len(placeholder_figures) == 0
                and len(checksum_mismatches) == 0
            ),
        }

    def compile_latex_documents(self, strict: bool = False) -> List[str]:
        """
        Audits the LaTeX files for hotwired numbers/figures, and compiles
        `main.tex`, `cover_letter.tex`, and `response_to_reviewers.tex`.
        If `strict=True`, raises an error if any placeholder stats/figures
        or hotwired LaTeX literals remain.
        """
        audit = self.audit_pipeline()
        if audit["hotwired_tex_violations"]:
            raise RuntimeError(
                "HOTWIRING DETECTED IN LATEX FILES:\n  "
                + "\n  ".join(audit["hotwired_tex_violations"])
            )
        if strict and not audit["strict_submission_ready"]:
            raise RuntimeError(
                "STRICT SUBMISSION GATE FAILED — Placeholder data or figures still active:\n"
                f"  Placeholder stat blocks: {audit['placeholder_stat_blocks']}\n"
                f"  Placeholder figures:     {audit['placeholder_figures']}\n"
                f"  Checksum mismatches:     {audit['figure_checksum_mismatches']}\n"
                "Run your analysis & figure scripts via `pipeline_api.py` to replace all placeholders before `--strict` build."
            )

        compiled = []
        for tex_name in ["main.tex", "cover_letter.tex", "response_to_reviewers.tex"]:
            tex_file = MANUSCRIPT_DIR / tex_name
            if tex_file.exists():
                cmd = [
                    "/Library/TeX/texbin/latexmk",
                    "-pdf",
                    "-interaction=nonstopmode",
                    "-halt-on-error",
                    tex_name
                ]
                subprocess.run(cmd, cwd=str(MANUSCRIPT_DIR), check=True)
                compiled.append(str(tex_file.with_suffix(".pdf")))
        return compiled


def main():
    parser = argparse.ArgumentParser(description="OrientExoML Reproducible Manuscript & Figure Pipeline API")
    parser.add_argument("--request-json", type=str, help="JSON string or path to JSON request file")
    parser.add_argument("--init-placeholders", action="store_true", help="Initialize placeholder figures and Figure 7")
    parser.add_argument("--audit", action="store_true", help="Audit LaTeX files for hotwiring and check provenance status")
    parser.add_argument("--compile", action="store_true", help="Generate stats macros and compile all LaTeX documents")
    parser.add_argument("--strict", action="store_true", help="Fail compilation if any placeholder stats or figures remain")
    args = parser.parse_args()

    api = OrientExoPipelineAPI()

    if args.init_placeholders:
        api.extract_placeholder_screenshots()
        if api.data.get("figures", {}).get("figure2", {}).get("status") != "generated_from_analysis":
            api.stamp_figure2_panels()
        fig7_status = (
            "generated_from_analysis"
            if api.data.get("provenance", {}).get("mvpa_summary_metrics", {}).get("status") == "computed_from_raw_data"
            else "placeholder_from_manuscript"
        )
        api.generate_figure7_mvpa_summary(status=fig7_status)
        api.save()

    if args.request_json:
        req_str = args.request_json
        if os.path.exists(req_str):
            with open(req_str, "r", encoding="utf-8") as f:
                payload = json.load(f)
        else:
            payload = json.loads(req_str)
        resp = api.handle_request(payload)
        print(json.dumps(resp, indent=2))

    if args.audit:
        report = api.audit_pipeline()
        print(json.dumps(report, indent=2))
        if report["hotwired_tex_violations"]:
            sys.exit(1)
        if args.strict and not report["strict_submission_ready"]:
            sys.exit(2)

    if args.compile:
        api.save()
        pdfs = api.compile_latex_documents(strict=args.strict)
        print("Compiled PDFs:", pdfs)


if __name__ == "__main__":
    main()
