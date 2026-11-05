# End-to-End Reproducible Manuscript & Figure Makefile for OrientExoML (CerCor-2022-00354)

PYTHON := $(shell if [ -x ../.venv/bin/python ]; then echo ../.venv/bin/python; else echo python3; fi)
LATEXMK := /Library/TeX/texbin/latexmk

.PHONY: all reconstruct stats figures audit pdf submission-ready clean

all: reconstruct submission-ready

reconstruct:
	$(PYTHON) scripts/generate_figure1_task.py
	$(PYTHON) scripts/run_behavior_pipeline.py
	$(PYTHON) scripts/eeg_preproc_and_feats.py
	$(PYTHON) scripts/run_mvpa_and_figures.py

stats:
	$(PYTHON) scripts/generate_latex_stats.py

figures:
	$(PYTHON) scripts/pipeline_api.py --init-placeholders

audit:
	$(PYTHON) scripts/pipeline_api.py --audit

pdf: stats audit
	$(PYTHON) scripts/pipeline_api.py --compile

# Strict submission gate: fails if ANY stat block or figure is still a placeholder
# or if any LaTeX file has hotwired numbers/paths or unverified figure SHA256 hashes.
submission-ready:
	$(PYTHON) scripts/pipeline_api.py --audit --strict --compile

clean:
	$(LATEXMK) -c main.tex cover_letter.tex response_to_reviewers.tex
