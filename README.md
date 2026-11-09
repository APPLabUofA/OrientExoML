# OrientExoML: End-to-End Reproducible Pipeline for Covert Attention, Alpha/CFC/Connectivity MVPA, and Manuscript Generation

**Title:** *Connecting Covert Attention and Visual Perception to the Spatiotemporal Dynamics of Alpha Band Activity, Cross-Frequency Coupling (CFC), and Functional Connectivity using Multivariate Pattern Analysis (MVPA)*  
**Authors:** Sarah S. Sheldon, Alona Fyshe, & Kyle E. Mathewson  
**Affiliation:** Attention, Perception, and Performance Lab (APPLab), Department of Psychology, Department of Computing Science, & Neuroscience and Mental Health Institute, University of Alberta

---

## Overview & Companion Links

This repository contains the complete, anti-hotwired Python and LaTeX reproducibility pipeline that reconstructs all behavioral analyses, EEG spectral parameterization (FOOOF), single-trial Cross-Frequency Coupling (CFC) and Functional Connectivity (FC) feature extraction, spatiotemporal LIBSVM/SVR Multivariate Pattern Analysis (MVPA) with Haufe activation patterns, publication-ready Figures 1–7, dynamic LaTeX statistical macros, and compiled manuscript PDFs from the raw data.

- **Open Science Framework (Raw Behavioral & BrainVision EEG Data):** <https://osf.io/u6fgm/>
- **Cued Orientation Task & Behavioral Analysis Repository (MATLAB):** <https://github.com/APPLabUofA/CuedOrientTask_paper>
- **Original MVPA & DDTBOX Analysis Repository (MATLAB/Python):** <https://github.com/APPLabUofA/MVPA_CuedOrientTask_paper>
- **Baseline Orientation Perception Task Repository:** <https://github.com/APPLabUofA/OrientTask_paper>

---

## Repository Structure

- [`manuscript/`](manuscript/)
  - [`main.tex`](manuscript/main.tex) / [`main.pdf`](manuscript/main.pdf) — Full revised manuscript in *Cerebral Cortex* format with embedded tables and Figures 1–7.
  - [`cover_letter.tex`](manuscript/cover_letter.tex) / [`cover_letter.pdf`](manuscript/cover_letter.pdf) — Resubmission cover letter.
  - [`response_to_reviewers.tex`](manuscript/response_to_reviewers.tex) / [`response_to_reviewers.pdf`](manuscript/response_to_reviewers.pdf) — Point-by-point response to Editor-in-Chief and Reviewers 1 & 2.
  - [`analysis_results.json`](manuscript/analysis_results.json) — Single source of truth for all computed statistics, figure SHA256 hashes, and provenance tracking.
  - [`generated_stats.tex`](manuscript/generated_stats.tex) — Auto-generated LaTeX `\newcommand` macros injected into the manuscript and response documents (zero hardcoded numbers in `.tex`).
  - [`figures/`](manuscript/figures/) — Publication-ready 300 DPI PNG and vector PDF figures (`figure1_task` through `figure7_mvpa_summary`).
  - [`scripts/`](manuscript/scripts/)
    - [`download_osf_eeg.py`](manuscript/scripts/download_osf_eeg.py) — Downloads all raw behavioral `.mat` and BrainVision `.eeg/.vhdr/.vmrk` files from OSF (`https://osf.io/u6fgm/`).
    - [`reconstruct_behavior_and_fig1.py`](manuscript/scripts/reconstruct_behavior_and_fig1.py) — Reconstructs Figure 1 (task schematic) and Figure 2 + 10,000-permutation $2 \times 2$ $F_{\max}$ behavioral ANOVAs from `data/behavior/*.mat`.
    - [`eeg_preproc_and_feats.py`](manuscript/scripts/eeg_preproc_and_feats.py) — Numba-accelerated EEG preprocessing (FIR filtering, mastoid re-referencing, HEOG split-half step saccade & blink rejection), FOOOF $1/f$ + alpha peak parameterization, and single-trial extraction of Alpha Amplitude, 5 CFC metrics (`AAC`, `MVL_Özkurt`, `MI`, `PAC_GLM`, `PAC_PLV`), and 4 FC metrics (`wPLI`, `Correlation`, `PLV`, `iCOH`).
    - [`run_mvpa_and_figures.py`](manuscript/scripts/run_mvpa_and_figures.py) — Numba-accelerated LIBSVM C-SVC ($C=1$) & $\epsilon$-SVR ($C=0.1, \epsilon=0.1$) 10-fold $\times$ 10-repeat spatiotemporal decoding, associative Haufe feature weight transformation, 10,000-permutation Phipson & Smyth (2010) $t_{\max}$ testing, and rendering of Figures 3–7.
    - [`pipeline_api.py`](manuscript/scripts/pipeline_api.py) — Strict anti-hotwiring audit, LaTeX macro generator, and `latexmk` PDF compiler.
- [`data/`](data/)
  - `data/behavior/` — Raw single-trial MATLAB behavioral files (`004_Orient_Exo.mat` – `040_Orient_Exo.mat`).
  - `data/eeg/` — BrainVision header (`.vhdr`) and marker (`.vmrk`) files (raw `.eeg` binaries are downloaded automatically by `download_osf_eeg.py`).
  - `data/cache/group_mvpa_results.npz` — Precomputed group-level MVPA trajectories, $t_{\max}$ statistics, and Haufe feature weights across all 10 metrics.

---

## Reproducing the Pipeline

```bash
# 1. Re-run strict audit, regenerate LaTeX macros & Figure 7, and compile all PDFs (~10 seconds):
cd manuscript
make submission-ready

# 2. Regenerate Figures 1-7 and all stats from cached group MVPA + behavioral data (~15 seconds):
cd manuscript
make reconstruct

# 3. Full raw EEG preprocessing & feature extraction from scratch (requires raw .eeg files from OSF):
python3 manuscript/scripts/download_osf_eeg.py
python3 manuscript/scripts/eeg_preproc_and_feats.py
rm -f data/cache/group_mvpa_results.npz
python3 manuscript/scripts/run_mvpa_and_figures.py
```

