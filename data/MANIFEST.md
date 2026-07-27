# Manifest: every manuscript number to the file and script that produce it

Added 2026-07-27 in response to a referee request. The objection was that `data/`
contained several similarly named five-cell reports with materially different
numbers and nothing identifying which was authoritative, so a reader auditing the
retrieval-executed cell could arrive at 0.8998 or 0.9652 depending on which file
they opened. Only one of those appears in the paper.

Superseded scorings now live in `data/SUPERSEDED/` with a README explaining why
each was retired. This file is the forward map.

## Main text

| Manuscript item | Data file | Script that regenerates it |
|---|---|---|
| Table 1, matched factorial | `v3_five_cell_live_report.txt` | `code/61-rescore-matched.py --input data/v3_five_cell_live.json` |
| Table 2, lethal errors by locus type | `v3_five_cell_live_rows.json` | `code/61-rescore-matched.py --rows` then the locus split in `code/70-figures-revision.py` |
| Table 3, clustered intervals and variance | `v3_five_cell_live_stats.txt` | `code/65-hierarchical-stats.py --input data/v3_five_cell_live_rows.json` |
| Table 4, four-cohort vocabulary coverage | `v3_ancestry_four_cohorts.txt` | `code/64-ancestry-matched-analysis.py` |
| Figure 1, study design | (schematic, no data) | `code/46-figure-study-design.py` |
| Figure 2, factorial result | `v3_five_cell_live.json` | `code/70-figures-revision.py` |
| Figure 3, error localisation | `v3_model_caller_eval.json` | `code/70-figures-revision.py` |
| Figure 4, between-model spread | `v3_five_cell_live.json` | `code/70-figures-revision.py` |
| Figure 5, ancestry coverage | `v3_ancestry_four_cohorts.json` | `code/70-figures-revision.py` |
| Figure 6, adversarial specification | `v3_adversarial*.json` | `code/31-figure2-adversarial-spec.py` |
| Figure 7, real-genome interpretation | real-genome-arm outputs | `real-genome-arm/scripts/05_score_report.py` |
| PharmCAT comparison | `v3_pharmcat_comparison_report.txt` | `code/63-pharmcat-comparator.py` |
| N4, drug-keyed retrieval | `v3_rag_drugkeyed_fullgrid.json` | `code/62-rag-drugkeyed-fullgrid.py` |
| Per-cell provenance (models, calls, spend) | `v3_cell_provenance.txt` | `code/69-cell-provenance.py` |
| GeT-RM caller concordance, decomposed | `v3_getrm_disagreement_classes.txt` | `real-genome-arm/scripts/14_getrm_disagreement_classes.py` |
| Extracted-rule diff and execution | `v3_extracted_rules_eval.txt` | `code/67-extract-rules-from-prose.py`, then `code/68-evaluate-extracted-rules.py` |

## Not published, do not cite

| File | What it is |
|---|---|
| `v3_matched_factorial_report.txt` | default output path of `61-rescore-matched.py`; the report of whatever was last scored |
| `SUPERSEDED/v3_five_cell_report.txt` | pre-parser-fix, retired `drug_match` metric |
| `SUPERSEDED/v3_five_cell_common6_report.txt` | six-model intermediate subset |
| `v3_five_cell_matched.json` | earlier merge retaining legacy skill rows under their original model labels; useful only as the comparator in `69-cell-provenance.py` |

## Raw evaluations

Raw per-call records are excluded from git by `.gitignore` and archived on Zenodo
(concept DOI 10.5281/zenodo.20567742). Every file listed above is regenerable
from them with the script named beside it.
