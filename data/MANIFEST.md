# Manifest: submitted results to authoritative inputs and generators

This is the forward map for the Cell Genomics revision released under tag
`agentic-pgx-benchmark-v2.2`. Raw JSON evaluation rows are excluded from git by
design and must be downloaded from the Zenodo version cited in the manuscript
into `data/`; the tracked `.txt` summaries beside them carry the same numbers in
readable form. The frozen file `v3_input_normalisation_eight_model_freeze.json`
records every source hash, paired analysis unit, operational denominator and
sample-cluster interval for the input-normalisation arm, and is the sole input
to the input-normalisation figure and its operational table.

Display-item numbers below are those of the submitted manuscript, in which
seven display items appear in the main text and the remainder are supplementary.

## Main manuscript

| Manuscript item | Authoritative input | Generator or validator |
|---|---|---|
| Figure 1, study design | schematic plus matched design constants | `code/46-figure-study-design.py` |
| Table 1, matched five-cell comparison | `v3_five_cell_live.json`, `v3_five_cell_live_report.txt` | `code/61-rescore-matched.py` |
| Figure 2, accuracy and lethal-class errors across the five cells | `v3_five_cell_live.json` | `code/70-figures-revision.py` |
| Table 2, rule tables extracted from guideline prose | `v3_extracted_rules.json`, `v3_extracted_rules_eval.json` | `code/67-extract-rules-from-prose.py`, `code/68-evaluate-extracted-rules.py` |
| Table 3, four-cohort vocabulary coverage | `v3_ancestry_four_cohorts.json` | `code/64-ancestry-matched-analysis.py` |
| Figure 3, real-genome interpretation across four cohorts | `v3_realgenome_preds_4cohorts.tsv`, `v3_agent_vs_deterministic.json` | `real-genome-arm/scripts/15_agent_vs_deterministic.py`, `code/82-figure7-real-genome.py` |
| Figure 4, input normalisation with and without definitions | `v3_input_normalisation_eight_model_freeze.json` | `code/85-freeze-eight-model-normalisation.py`, `code/87-figure9-input-normalisation-eight.py` |

## Supplementary

| Supplementary item | Authoritative input | Generator or validator |
|---|---|---|
| Figure S1, the query issued in each of the five cells | runner dry-run output | `code/84-figure-five-cell-queries.py` |
| Figure S2, information without action | `v3_five_cell_live.json` | `code/86-figures-s2-s3-five-cell.py` |
| Figure S3, correctness by coincidence | `v3_five_cell_live.json` | `code/86-figures-s2-s3-five-cell.py` |
| Figure S4, residual error localised to the input call | `v3_model_caller_eval.json`, `v3_five_cell_live.json` | `code/70-figures-revision.py` |
| Figure S5, between-model spread | `v3_five_cell_live.json` | `code/70-figures-revision.py` |
| Figure S6, four-cohort coverage | `v3_ancestry_four_cohorts.json` | `code/70-figures-revision.py` |
| Figure S7, bidirectional corrupted contract | `v3_adversarial_scrambled.json`, `v3_adversarial_reverse.json` | `code/80-figure6-adversarial-contract.py` |
| Figure S8, conditional invariance under population framing | `v3_armA9_armBv2_POP.json` | `code/42-armAB-population-sweep.py` |
| Table S1, headline metrics under both scorers | `v3_five_cell_live_report.txt` | `code/61-rescore-matched.py`, `code/89-supp-headline-and-per-gene.py` |
| Table S2, lethal-class errors by locus type | `v3_lethal_case_level.json` | `code/71-lethal-case-level-stats.py` |
| Table S3, per-gene lethal-class errors | `v3_matched_scored_rows_all5.json` | `code/61-rescore-matched.py`, `code/89-supp-headline-and-per-gene.py` |
| Table S4, clustered intervals and variance components | `v3_five_cell_live_rows.json` | `code/65-hierarchical-stats.py` |
| Table S5, crossed random-effects model | `v3_five_cell_live_rows.json`, `v3_crossed_mixed_model.json` | `code/73-crossed-mixed-model.py` |
| Table S6, bidirectional adversarial specification test | `v3_adversarial_scrambled.json`, `v3_adversarial_reverse.json` | `code/14-adversarial-scrambled-spec.py`, `code/14b-adversarial-reverse-spec.py` |
| Table S7, retrieval indexing in the six-model rerun | `v3_rag_drugkeyed_fullgrid.json` | `code/62-rag-drugkeyed-fullgrid.py` |
| Table S8, input normalisation by gene | deposited raw normalisation rows, `v3_input_normalisation_eight_model_freeze.json` | `code/88-supp-normalisation-by-gene.py` |
| Table S9, input-normalisation operational outcomes | `v3_input_normalisation_eight_model_freeze.json` | `code/85-freeze-eight-model-normalisation.py` |
| Key resources table (STAR Methods, unnumbered) | model identifiers and corpus releases | maintained by hand; `MODEL-VERSIONS.md` |

## Other reported results

| Item | Authoritative input | Generator or validator |
|---|---|---|
| PharmCAT comparison | `v3_pharmcat_comparison.json` | `code/63-pharmcat-comparator.py` |
| GeT-RM caller validation and its decomposition | `v3_getrm_disagreement_classes.json`, `v3_getrm_caller_truth.json` | `code/77-normalisation-disagreement-classes.py` |
| Pooled input-normalisation figures without definitions | `v3_input_normalisation_eval.json` | `code/75-evaluate-input-normalisation.py` |
| Per-cell provenance and spend | `v3_cell_provenance.json` | `code/69-cell-provenance.py` |
| C14 and C15 parser effects | stored raw input-normalisation responses | `code/79-reparse-normalisation.py`, `CORRECTIONS.md` |

## Release checks

1. `python3 code/72-validate-v30-numbers.py --manuscript <final .docx>` must
   report zero failures. It recomputes every registered claim from the data and
   checks that the manuscript states it. A claim absent from that script is
   unverified, so the count it prints matters as much as the result.
2. `python3 code/88-supp-normalisation-by-gene.py` must exit 0; it refuses to
   write if its own totals disagree with the eight-model freeze.
3. `python3 code/81-validate-submission-package.py --package-dir <dir>` against
   the complete submission package.
4. Rebuild and inspect the figures from the deposited rows.

A release is not synchronised until the repository tag, the Zenodo version, the
manuscript, the supplement, the response letter and the figure files all refer
to the same snapshot. The tag cited by the manuscript is immutable: if the
content must change, cut a new tag rather than moving the cited one.

## Superseded files

Files in `data/SUPERSEDED/`, `v3_input_normalisation_definitions.txt` and
`v3_input_normalisation_seven_model_freeze.json` are historical intermediates.
The seven-model freeze remains deposited for provenance and must not be used
for any submitted number; the eight-model freeze supersedes it.
