# Real-genome arm: agentic pharmacogenomics across four cohorts

This directory contains the complete, reproducible pipeline for the **real-genome
validation** of the agentic pharmacogenomics benchmark. It is the empirical
implementation of the agentic-genomics framework set out in the companion
Perspective (Corpas, Fatumo & Guio, 2026), and the first test of large-scale
agentic genomic interpretation on **real human genomes** rather than curated cases.

Where the main benchmark uses 110 canonical CPIC cases, this arm runs the agent on
star-allele diplotypes **actually observed in real individuals** from four
cohorts, called by a single deterministic caller (PyPGx) so the comparison uses
one input pipeline throughout:

| Cohort | Population | Individuals | Source |
|---|---|---|---|
| Uganda Genome Resource (UGR) | East African | 6,407 | controlled-access (DTA) |
| Peruvian Genome Project | Admixed Latin American (7 subpopulations) | 150 / 736 | published / cohort |
| Corpas family | Spanish European | 4 | personal whole-genome sequences |
| 1000 Genomes IBS | Iberian European | 93 | 1000 Genomes Project |

The design follows the Perspective's architecture exactly: a **deterministic caller**
(PyPGx) handles the input step (variants → diplotype), and the **agent** performs the
interpretation step (diplotype → CPIC phenotype/recommendation). We score the agent
against the caller's CPIC phenotype and report four mutually exclusive outcomes the canonical benchmark
cannot show:

1. **Correct**, matching the caller's CPIC phenotype.
2. **Abstained**, an explicit non-answer on a scorable caller state.
3. **Clinically wrong**, a usable phenotype that does not match the caller.
4. **Invalid output**, no unique usable `PHENOTYPE:` line, empty output or provider failure.

Caller states labelled `Indeterminate` have no reference phenotype and are
excluded. `Uncertain Susceptibility` is a valid CPIC phenotype and remains in
the scored denominator.

## Pipeline (GRCh37 throughout: no liftover)

```
data (PLINK or VCF, GRCh37)
  └─ 01_extract_pgx_regions.sh   plink2: extract PGx gene regions -> bgzipped VCF
       └─ 02_call_diplotypes_pypgx.sh   PyPGx run-chip-pipeline per gene (GRCh37)
            └─ 03_aggregate_diplotypes.py   per-(gene,diplotype) table + phenotype
                 └─ 21_run_agent_four_cohorts.py   model panel plus raw provider provenance
                      └─ 22_score_agent_realgenome.py   four-category scoring and model audit
                           └─ ../../code/82-figure7-real-genome.py   Figure 3
```

### Quickstart (per cohort)

```bash
# 1. extract PGx regions from a GRCh37 PLINK set (UGR / Peru chip)
scripts/01_extract_pgx_regions.sh /path/to/ugr_data /work/ugr_pgx_grch37

# 2. call diplotypes with PyPGx (single consistent caller, GRCh37)
scripts/02_call_diplotypes_pypgx.sh /work/ugr_pgx_grch37.vcf.gz /work/ugr_calls GRCh37

# 3. aggregate observed real diplotypes
python3 scripts/03_aggregate_diplotypes.py /work/ugr_calls UGR /work/ugr_diplotypes.tsv

# 4. concatenate the four cohort state files, then run the model panel
#    (needs provider API keys in the environment)
python3 scripts/21_run_agent_four_cohorts.py

# 5. derive four-category scores and the provider/model audit
python3 scripts/22_score_agent_realgenome.py

# 6. regenerate Figure 3
python3 ../code/82-figure7-real-genome.py
```

The runner checkpoints each complete provider response to JSONL and resumes only
missing or failed calls. It preserves the exact prompt and request, raw provider
JSON, returned model identifier, response text, stop reason, token usage,
response identifier, timestamps, latency, retry count and error body. The audit
requires 536 returned records per model before Figure 3 is generated.

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
