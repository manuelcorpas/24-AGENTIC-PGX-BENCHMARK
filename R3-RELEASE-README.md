# Agentic pharmacogenomics benchmark, version 1.9.0

Version DOI: https://doi.org/10.5281/zenodo.22806813
Concept DOI: https://doi.org/10.5281/zenodo.20567742
Prior version: https://doi.org/10.5281/zenodo.22141182

Third-revision reproducibility data and code for “Authored skills improve accuracy in agentic pharmacogenomics but not end-to-end clinical validation”. Prepared 17 September 2026.

The ZIP is the complete public reproducibility package. Extract it and work from its REPOSITORY folder. It contains preserved historical data plus separately identified R3 corrections and derived results. No model calls are needed. Original references, prompts and outputs are preserved; this update does not silently replace them.

## Start here

- docs/R3-2026-09-16/ANALYSIS-DEFINITIONS.md: frozen scorers, denominator corrections, both wording rules, original and exclusion analyses, and 96.31% zero-credit arithmetic sensitivity.
- docs/R3-2026-09-17/README.md and DATA/publication-checks.json: all 110 screen means, 102 cases outside the screen, one additional case under a five-percentage-point relaxation, and the frozen scorer's uninformative 0.5 fallback (96.76% overall for authored-rule execution).
- docs/R3-2026-09-16/CANDIDATE-AUDIT.md: eight candidates, observed responses, cited guideline evidence and dispositions. This is candidate review, not independent clinical adjudication of all cases.
- docs/R3-2026-09-16/DATA/test-cases-corrected-v1.json: separately versioned RYR1 reference correction. Opaque inputs, evaluator references, mapping and 550 prompt-identity checks are in the same directory.
- docs/R3-2026-09-17/DATA/Table-S10.csv and Table-S11.csv: complete exclusion analysis and candidate-screen table. Values are fractions before display rounding.
- data/: all raw inputs, including the twelve historical files holding 9,557 registered re-parsed rows. The separately labelled canonical V4 Flash panel and eight-model freeze are also present. Historical DeepSeek files are retained for the parser audit; they do not replace V4 Flash in the current panel.
- RELEASE-MANIFEST.json: every packaged file with size and SHA-256 hash, source commit, and deliberate omissions of confidential document builders and their tests.

## Reproduce

Use Python 3.11 or later with numpy, scipy, pandas, matplotlib, statsmodels, Pillow, pytest, python-docx, lxml, requests and python-dotenv. From REPOSITORY:

    python docs/R3-2026-09-16/SOURCE/01-regenerate.py
    python docs/R3-2026-09-16/SOURCE/02-build-audit.py
    python code/99-publication-checks.py
    python code/79-reparse-normalisation.py
    python -m pytest tests/test_third_revision.py tests/test_scorer_freeze.py tests/test_rescore_matched.py tests/test_extracted_rules_eval.py tests/test_hierarchical_stats.py tests/test_matched_factorial.py tests/test_caller_truth_eval.py tests/test_figure_r1_layout.py tests/test_r3_supplementary_figures.py tests/test_r3_regeneration.py tests/test_r3_audit.py

The existing RAG-generation REML fit is not converged and is retained as a labelled diagnostic. The historical notes dated 16 September describe pre-publication staging; this version's public record and this README supersede their availability wording. The older root README describes the historical project; use this release guide for R3.

Manuscript DOCX/PDF files, response and cover letters, redlines and source embedding confidential correspondence are excluded. They belong to the separate author review bundle. Individual-level genotypes from the three human cohorts are not redistributed.
