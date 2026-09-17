# Claude Science on the agentic-PGx benchmark, September 2026

Everything needed to check the claims in the accompanying email, and to disagree with them.
Nothing here is a summary of a result you cannot recompute; every number in the email comes
out of one command below, run against the CSVs in this folder and the frozen scorer already in
this repository.

## Recompute in one command

```bash
python evaluations/2026-09-claude-science/verify.py
```

That rescores all six replication runs and the pilot with `code/61-rescore-matched.py` at this
commit, prints the per-cell means and spreads, reproduces the identifier-leak split, and exits
non-zero if anything differs from `results/expected.json`. `pytest evaluations/2026-09-claude-science/tests`
runs the same checks plus a hermetic one on synthetic answers.

## What was done

A blinded package of the 110 CPIC Level A cases (gene, drug, genotype only; identifiers hashed)
and the same CPIC guideline corpus used by the paper's RAG arm was handed to Claude Science
(Anthropic's desktop research agent, Opus 5 by the UI) with a written brief and no skill. Two
conditions, three replicates each, one fresh project per run so nothing carried over:

| cell | knowledge | runs |
|---|---|---|
| `api_corpus` | corpus plus live access to `api.cpicpgx.org` | `runs/api_corpus/rep1..3` |
| `corpus_only` | corpus only; every network request denied | `runs/corpus_only/rep1..3` |

A pilot run on 15 September used the original, semantic case identifiers (`cyp2d6_codeine_pm`)
and is kept in `pilot/` because it is the evidence for the identifier-leak analysis, not
because its accuracy figure is usable.

## The numbers, as reproduced by `verify.py`

| cell | phenotype accuracy | abstentions per run | "lethal-class" errors by the paper's wording rule |
|---|---|---|---|
| `api_corpus` | 0.942, sd 0.004 | 0, 5, 0 | 6, 6, 6 of 14 |
| `corpus_only` | 0.918, sd 0.021 | 2, 2, 2 (same two CYP4F2 cases every time) | 5, 5, 0 of 14 |
| paper, authored skill arms | 0.966 / 0.967 | | |

The pilot scored 0.973 and is withdrawn: on the 104 cases whose identifier encoded the
phenotype, the pilot scored 0.990 and the hashed rerun 0.962; on the 6 cases whose identifier
did not, both scored 0.667. `results/DATASET-ISSUE-case-id-leak.md` documents the defect in the
published dataset itself. The paper's own model arms never saw the identifiers.

## What the numbers do not say

- The recommendation column is scored by a wording rule that accepts only "avoid" or
  "contraindicated". Every one of the "lethal-class" failures above says "do not use", "is not
  recommended" or names an alternative. Clinical correctness has not been independently
  adjudicated in either direction. Do not read 0.57 as a clinical error rate.
- Three replicates per cell, one gene-drug domain, one agent product. No interval is claimed.
- The paper's arms are not a matched control: different model, runtime, snapshot and workflow.
- "No skill used" was inferred from the agent's own methods notes; its runtime inputs and
  built-in skills were not audited. `audit/BLINDING-AUDIT.json` lists every leakage opportunity
  found, including the corpus metadata.

## The finding that matters more than the accuracy

`audit/TARGETED-CHECKS.json` records three probes run through the pilot's own generated
pipeline: with no genotype supplied it fills reference alleles and returns Normal Metabolizer;
with one locus unreported it fills that locus as reference; and given "laboratory-confirmed
\*1/\*1, copy number 2; no duplication detected" it returns \*1/\*1xN, Ultrarapid Metabolizer,
because a substring match ignored the negation. None of that moves the accuracy score, because
the benchmark never asks those questions.

## Folder map

- `runs/` the six replication runs: answers CSV, the agent's methods note, its lookup tables and
  code. The agent's cached copies of the CPIC API (up to 115 MB per run) are omitted; the source
  is `api.cpicpgx.org`, and each `METHODS.md` names the tables and retrieval date.
- `keys/id_map.json` surrogate id to real case id. Publishing it un-blinds these runs for you;
  a future run needs fresh hashes (`scripts/00-build_package.py`).
- `package/` the exact blinded inputs handed over, with SHA-256 in `MANIFEST.json`.
- `pilot/` the 15 September run, semantic ids, plus its scored rows.
- `results/` `COMPARISON.md` (the working write-up, with every correction left in),
  `expected.json` (what `verify.py` must reproduce), the dataset issue.
- `audit/` an independent audit of the pilot session: locked scores recomputed, an exact offline
  replay of all 110 rows, the three input-validity probes, the blinding audit.
- `clawbio-skill-audit/` the classification of all 97 ClawBio skills against `PHASE1-RUBRIC.md`
  by five independent readers, with the citation check of every justification. Signal labels in
  `inventory.csv` are a file-level heuristic that the readers overruled on 42 of 97 rows; they
  are not the classification.
- `scripts/` the package builder, run builder and aggregator as used. They resolve the scorer
  through `~/dev/24-AGENTIC-PGX-BENCHMARK`; `verify.py` is the repo-relative entry point.

## Provenance

Benchmark and scorer: this repository at tag `agentic-pgx-benchmark-v2.10`, commit `eed8834`,
Zenodo 10.5281/zenodo.22141182. Runs: 15 to 16 September 2026, Claude Science beta, Opus 5.
Audit: 15 September 2026, independent of the session that produced the pilot.
