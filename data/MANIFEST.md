# Manifest: submitted results to authoritative inputs and generators

This is the forward map for the corrected Cell Genomics revision. Raw JSON
evaluation rows are excluded from git and must be downloaded from the corrected
Zenodo version into `data/`. The frozen file
`v3_input_normalisation_seven_model_freeze.json` records every source hash,
paired analysis unit, operational denominator and sample-cluster interval.

## Main manuscript

| Manuscript item | Authoritative input | Generator or validator |
|---|---|---|
| Table 1, matched five-cell comparison | `v3_five_cell_live.json` | `code/61-rescore-matched.py` |
| Table 2, lethal errors by locus type | `v3_lethal_case_level.json` | `code/71-lethal-case-level-stats.py` |
| Table 3, clustered intervals and variance | `v3_five_cell_live_rows.json` | `code/65-hierarchical-stats.py` |
| Table 4, four-cohort vocabulary coverage | `v3_ancestry_four_cohorts.json` | `code/64-ancestry-matched-analysis.py` |
| Figure 1, study design | schematic plus matched design constants | `code/46-figure-study-design.py` |
| Figure 2, matched five-cell result | `v3_five_cell_live.json` | `code/70-figures-revision.py` |
| Figure 3, conditional error localisation | `v3_model_caller_eval.json`, `v3_five_cell_live.json` | `code/70-figures-revision.py` |
| Figure 4, between-model spread | `v3_five_cell_live.json` | `code/70-figures-revision.py` |
| Figure 5, four-cohort coverage | `v3_ancestry_four_cohorts.json` | `code/70-figures-revision.py` |
| Figure 6, bidirectional corrupted contract | `v3_adversarial_scrambled.json`, `v3_adversarial_reverse.json` | `code/80-figure6-adversarial-contract.py` |
| Figure 7, real-genome interpretation | `v3_realgenome_preds_4cohorts.tsv`, `v3_agent_vs_deterministic.json` | `real-genome-arm/scripts/15_agent_vs_deterministic.py`, `code/82-figure7-real-genome.py` |
| Figure 8, input normalisation | `v3_input_normalisation_seven_model_freeze.json` | `code/83-freeze-seven-model-normalisation.py`, `code/78-figure-input-normalisation.py` |

## Supplementary and response items

| Item | Authoritative input | Generator or validator |
|---|---|---|
| Tables S8 and S9 | corrected input-normalisation JSON rows and seven-model freeze | `code/75-evaluate-input-normalisation.py`, `code/83-freeze-seven-model-normalisation.py`, submission builder |
| Table S7, crossed model | `v3_five_cell_live_rows.json` | `code/65-hierarchical-stats.py` |
| PharmCAT comparison | `v3_pharmcat_comparison.json` | `code/63-pharmcat-comparator.py` |
| Drug-keyed retrieval control | `v3_rag_drugkeyed_fullgrid.json` | `code/62-rag-drugkeyed-fullgrid.py` |
| Extraction experiment | `v3_extracted_rules.json` | `code/68-evaluate-extracted-rules.py` |
| Per-cell provenance | `v3_five_cell_live.json` | `code/69-cell-provenance.py` |
| C14 and C15 effects | stored raw input-normalisation responses | `code/79-reparse-normalisation.py`, `CORRECTIONS.md` |

## Release checks

Run the registered-number validator against the exact manuscript being sent,
run `code/81-validate-submission-package.py` against the complete submission
package, then rebuild and inspect Figures 3, 6, 7 and 8 from the deposited rows. A release is
not synchronized until the repository tag, corrected Zenodo version, manuscript,
supplement, response letter and figure files all refer to the same snapshot.

## Superseded files

Files in `data/SUPERSEDED/` and `v3_input_normalisation_definitions.txt` are
historical intermediates and must not be used for the corrected submission.
