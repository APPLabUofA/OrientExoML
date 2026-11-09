# OrientExoML End-to-End Reproducible Manuscript & Figure Pipeline (`CerCor-2022-00354`)

This directory contains the automated, anti-hotwired LaTeX and Python pipeline for the revised *Cerebral Cortex* resubmission:
> **"Connecting Covert Attention and Visual Perception to the Spatiotemporal Dynamics of Alpha Band Activity, Cross-Frequency Coupling (CFC), and Functional Connectivity using Multivariate Pattern Analysis (MVPA)"**  
> Sarah S. Sheldon, Alona Fyshe, and Kyle E. Mathewson

---

## 1. Anti-Hotwiring & Provenance Architecture

To guarantee that no statistical numbers, figure arrays, or file paths are ever hardcoded ("hotwired") in `.tex` files or plotting scripts:

1. **Zero Hardcoded Numbers in LaTeX (`generated_stats.tex`)**:
   Every empirical number, sample size, $F_{\max}$ / $t_{\max}$ statistic, FOOOF parameter, individual alpha frequency, MVPA decoding window, percentage of significant windows, and figure file path in [`main.tex`](main.tex), [`cover_letter.tex`](cover_letter.tex), and [`response_to_reviewers.tex`](response_to_reviewers.tex) is injected via a `\stat*` or `\figPath*` macro generated from [`analysis_results.json`](analysis_results.json) by [`scripts/generate_latex_stats.py`](scripts/generate_latex_stats.py).
2. **Zero Hardcoded Arrays in Figure Generators**:
   [`scripts/pipeline_api.py`](scripts/pipeline_api.py) (`generate_figure7_mvpa_summary`) and [`scripts/regenerate_mvpa_figures.py`](scripts/regenerate_mvpa_figures.py) contain zero hardcoded data arrays; they read all means, SEMs, and significance indices directly from [`analysis_results.json`](analysis_results.json) (`mvpa_summary_metrics`) or external JSON payloads exported by MATLAB/Python.
3. **Provenance & SHA256 Checksum Tracking**:
   [`analysis_results.json`](analysis_results.json) tracks a `provenance` block for every statistical category (`participants`, `behavior`, `fooof`, `mvpa`, `mvpa_summary_metrics`) and every figure (`figure1`--`figure7`), recording:
   - `status`: `"placeholder_from_manuscript"` / `"placeholder_screenshot"` vs. `"computed_from_raw_data"` / `"generated_from_analysis"`
   - `source_script`: which analysis script pushed the data
   - `updated_at`: UTC ISO-8601 timestamp
   - `sha256`: cryptographic hash of the generated figure file on disk
4. **Automated Anti-Hotwiring Linter & Strict Submission Gate (`--audit` and `--strict`)**:
   - Running `make audit` (or `python3 scripts/pipeline_api.py --audit`) scans all `.tex` files for any literal numbers in statistical expressions (`F_{\max}(...) = 12.3`, `N = 28`, `M = 10.4`) or direct `\includegraphics{...}` paths that bypass `\figPath*` macros, and verifies all figure SHA256 hashes.
   - Running `make submission-ready` (or `python3 scripts/pipeline_api.py --audit --strict --compile`) enforces a **Strict Submission Gate** that blocks final PDF compilation if *any* statistical block or figure is still marked as a `"placeholder_*"` or has a mismatched SHA256 hash.

---

## 2. Build & Audit Commands

```bash
cd manuscript

# Draft build (allows placeholder stats/screenshots, blocks hotwired LaTeX literals):
make all

# Run anti-hotwiring and provenance audit:
make audit

# Final submission gate (fails until all placeholders are replaced by real analysis outputs):
make submission-ready
```

---

## 3. API Request Payloads for Analysis Scripts & Figure Generators

When your MATLAB/Python analysis scripts finish running on raw data, have them call [`scripts/pipeline_api.py`](scripts/pipeline_api.py). Pushing via `update_stats` or `register_figure` automatically flips that component's provenance status from `"placeholder_*"` to `"computed_from_raw_data"` / `"generated_from_analysis"`, records the `source_script` and timestamp, updates `generated_stats.tex`, and (for `mvpa_summary_metrics`) automatically re-renders **Figure 7**.

### A. Push Real Analysis Statistics (`action: "update_stats"`)

```bash
python3 manuscript/scripts/pipeline_api.py --request-json '{
  "action": "update_stats",
  "source_script": "MVPA_CuedOrientTask_paper/run_fooof_and_behavior.m",
  "updates": {
    "fooof": {
      "iaf_grand_mean_hz": "10.42",
      "iaf_grand_sd_hz": "0.88",
      "iaf_range_hz": "9.2--11.8"
    },
    "behavior": {
      "absolute_response_error": {
        "interaction": { "fmax": 17.19, "crit_fmax": 4.06, "p": "< 0.001" }
      }
    }
  }
}' --compile
```

### B. Push Real MVPA Summary Arrays & Auto-Regenerate Figure 7

```bash
python3 manuscript/scripts/pipeline_api.py --request-json '{
  "action": "update_stats",
  "source_script": "MVPA_CuedOrientTask_paper/summarize_mvpa_metrics.m",
  "updates": {
    "mvpa_summary_metrics": {
      "svm_left_peak": [58.2, 54.1, 51.4, 53.0, 51.8, 51.9, 50.8, 50.9, 50.4, 51.0],
      "svm_right_peak": [58.8, 53.9, 51.2, 52.6, 52.1, 51.5, 52.3, 50.5, 51.0, 50.6]
    }
  }
}' --compile
```

### C. Register a Newly Generated Figure Replacing a Placeholder (`action: "register_figure"`)

```bash
python3 manuscript/scripts/pipeline_api.py --request-json '{
  "action": "register_figure",
  "source_script": "manuscript/scripts/regenerate_mvpa_figures.py",
  "figure_id": "figure3",
  "file": "figures/figure3_svm_cfc.png",
  "status": "generated_from_analysis",
  "description": "Updated 300 DPI SVM CFC figure with 200 ms tick spacing and abbreviated labels"
}' --compile
```
