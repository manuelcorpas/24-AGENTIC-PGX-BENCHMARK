# Real-genome arm: agentic pharmacogenomics on real genomes across three continents

This directory contains the complete, reproducible pipeline for the **real-genome
validation** of the agentic pharmacogenomics benchmark. It is the empirical
implementation of the agentic-genomics framework set out in the companion
Perspective (Corpas, Fatumo & Guio, 2026), and the first test of large-scale
agentic genomic interpretation on **real human genomes** rather than curated cases.

Where the main benchmark uses 110 canonical CPIC cases, this arm runs the agent on
star-allele diplotypes **actually observed in real individuals** from three
ancestrally distinct cohorts, called by a single deterministic caller (PyPGx) so the
cross-population comparison is fair:

| Cohort | Population | Individuals | Source |
|---|---|---|---|
| Uganda Genome Resource (UGR) | East African | 6,407 | controlled-access (DTA) |
| Peruvian Genome Project | Admixed Latin American (7 subpopulations) | 150 / 736 | published / cohort |
| Corpas family | European | family genome | personal genome |

The design follows the Perspective's architecture exactly: a **deterministic caller**
(PyPGx) handles the input step (variants → diplotype), and the **agent** performs the
interpretation step (diplotype → CPIC phenotype/recommendation). We score the agent
against the caller's CPIC phenotype and report three things the canonical benchmark
cannot show:

1. **Real-genome accuracy** vs the canonical 96.2% (how far curated accuracy transfers).
2. **Abstention** on uncertain/Indeterminate diplotypes (real genomes are full of them).
3. **Population/equity behaviour** on ancestry-specific real alleles.

## Pipeline (GRCh37 throughout: no liftover)

```
data (PLINK or VCF, GRCh37)
  └─ 01_extract_pgx_regions.sh   plink2: extract PGx gene regions -> bgzipped VCF
       └─ 02_call_diplotypes_pypgx.sh   PyPGx run-chip-pipeline per gene (GRCh37)
            └─ 03_aggregate_diplotypes.py   per-(gene,diplotype) table + phenotype
                 └─ 04_run_agent_realgenome.py   model panel: diplotype -> phenotype
                      └─ 05_score_report.py   accuracy / abstention / errors by population
```

### Quickstart (per cohort)

```bash
# 1. extract PGx regions from a GRCh37 PLINK set (UGR / Peru chip)
scripts/01_extract_pgx_regions.sh /path/to/ugr_data /work/ugr_pgx_grch37

# 2. call diplotypes with PyPGx (single consistent caller, GRCh37)
scripts/02_call_diplotypes_pypgx.sh /work/ugr_pgx_grch37.vcf.gz /work/ugr_calls GRCh37

# 3. aggregate observed real diplotypes
python3 scripts/03_aggregate_diplotypes.py /work/ugr_calls UGR /work/ugr_diplotypes.tsv

# 4. run the agent (model panel) on the real diplotypes  (needs API keys in env)
python3 scripts/04_run_agent_realgenome.py /work/ugr_diplotypes.tsv /work/ugr_predictions.tsv

# 5. score + decompose vs cohort CPIC phenotype
python3 scripts/05_score_report.py /work/ugr_predictions.tsv /work/ugr_report.txt
```

Concatenate the per-cohort `*_diplotypes.tsv` before step 4 to score the panel across
all three populations in one run.

## Step 7 (N0): end-to-end executed pipeline with abstention

`scripts/07_executed_pipeline_n0.py` runs the deployment architecture the paper proposes:
it takes the deterministic caller's diplotypes (the `*_diplotypes.tsv` from step 3),
executes the validated skill's CPIC mapping in code, and **abstains** on no-calls,
uncertain-function (ambiguous) and out-of-scope diplotypes. It reports, per cohort,
**coverage** and **accuracy among emitted in-scope answers**, scored against the caller's
own CPIC phenotype (PyPGx) with the manuscript's phenotype scorer (`code/10-rescore-v3.py`).
This measures, rather than predicts, the "100% correctness among emitted answers" claim;
any disagreement with PyPGx is recorded as a falsifier. No new model calls (deterministic).

```bash
# after step 3 has produced a per-cohort (or concatenated) diplotypes.tsv:
python3 scripts/07_executed_pipeline_n0.py /work/all_cohorts_diplotypes.tsv /work/n0_executed.json
# offline sanity check on a built-in fixture (no data needed):
python3 scripts/07_executed_pipeline_n0.py --demo
```

Deterministic logic is unit-tested in `tests/test_n0_executed_pipeline.py`.

## Reproducing exactly

- **Software and versions**: see [`SOFTWARE.md`](SOFTWARE.md).
- **Obtaining the cohort data** (controlled-access UGR, Peruvian, Corpas): see
  [`DATA-ACCESS.md`](DATA-ACCESS.md). **No genotype data is included in this repo.**
- **Python deps**: `pip install -r requirements.txt`.
- **Config**: `config/pgx_regions_grch37.bed` (extraction regions),
  `config/genes.txt` (genes called), `config/models.txt` (model panel).

## Data governance (important)

This is a public repository. The Uganda Genome Resource is **controlled-access** under
a data transfer agreement; the Peruvian and Corpas genomes are personal/cohort data.
**No individual-level genotypes, VCFs, PLINK files, caller outputs, or API keys are or
should ever be committed here** (enforced by `.gitignore`). Reviewers reproduce by
obtaining each dataset through its own access route (DATA-ACCESS.md) and running the
pipeline above. Only code, configuration, software versions, and instructions are
public.
