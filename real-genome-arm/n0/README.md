# N0: end-to-end executed pipeline (mapping concordance)

Reproduces the N0 result cited in the Cell Genomics plan of revision (CELL-GENOMICS-D-26-00551).

## What N0 measures (read this before quoting the number)

N0 runs the deployment architecture end to end on three real cohorts: a deterministic
caller (PyPGx 0.26.0, GRCh37) calls the diplotype, the versioned skill executes the CPIC
mapping in code, and the pipeline abstains on no-calls, uncertain-function and out-of-scope
diplotypes.

**The scored quantity is concordance of the executed skill's CPIC mapping with an
independent implementation of CPIC (PyPGx), on emitted in-scope diplotype states.** It is
NOT end-to-end correctness. The calling step has no independent ground truth here: PyPGx
both calls the diplotype and supplies the phenotype it is scored against. Substituting PyPGx
for the model therefore makes input-call error unobservable, not zero. External truth for the
calling step (GeT-RM consensus diplotypes on the 1000 Genomes / Coriell overlap) is a
separate, planned experiment; until it is run, residual risk sits in caller accuracy and is
not measured here.

## Result (mapping concordance among emitted in-scope states)

| Cohort | Emitted | Concordance | 95% CI (Clopper–Pearson) |
|--------|---------|-------------|--------------------------|
| European family | 9/9 | 100% | 66.4%–100% |
| Peruvian Genome Project | 12/12 | 100% | 73.5%–100% |
| Uganda Genome Resource | 17/17 | 100% | 80.5%–100% |
| All cohorts | 38/38 | 100% | 90.7%–100% |

Zero disagreements. Coverage of observed diplotype states is 22% / 16% / 6% (European family
/ Peru / Uganda); the balance is abstentions, dominated by diplotypes outside the validated
skill's vocabulary. The disparity across cohorts appears as coverage (abstention), not as
wrong answers.

## Files
- `n0_input_3cohorts.tsv` — aggregate per-cohort (gene, diplotype, phenotype, n_carriers) state
  table. Counts only; no individual-level genotypes (consistent with the Data availability
  statement). Cohorts: European family, Peru (PGP), Uganda (UGR).
- `n0_result.json` — output of the run below.

## Reproduce (deterministic, 0 API calls)
```
python3 real-genome-arm/scripts/07_executed_pipeline_n0.py real-genome-arm/n0/n0_input_3cohorts.tsv /tmp/n0_out.json
# compare /tmp/n0_out.json to n0_result.json
```
