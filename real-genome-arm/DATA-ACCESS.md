# Obtaining the cohort data

No genotype data is included in this repository. Each cohort is obtained through its
own access route. The pipeline (see `README.md`) runs once you have the data locally.

## Uganda Genome Resource (UGR) — East African, n = 6,407 — CONTROLLED ACCESS

The UGR is controlled-access genomic data governed by a data transfer / access
agreement. It is **not** redistributable and is not in this repository. Access is
granted on application to the data custodians (MRC/UVRI and LSHTM Uganda Research
Unit / Uganda Medical Informatics Centre). Approved users receive a PLINK set
(`ugr_data.bed/.bim/.fam`, GRCh37, numeric chromosomes). The pipeline expects that
PLINK prefix as input to `01_extract_pgx_regions.sh`.

Governance: hold the data only on approved, access-controlled storage; do not make
it public or place it in any public repository or cloud share without authorisation.

## Peruvian Genome Project — Admixed Latin American (7 subpopulations) — n = 150 / 736

Published as the Peruvian Clinical Genome (Frontiers in Genetics, 2025,
DOI 10.3389/fgene.2025.1614021). Per-subpopulation PGx phenotype tables are in the
paper's supplementary materials (SUPPL-TABLE2 diplotype→phenotype by subpopulation;
SUPPL-TABLE3 per-individual FDA recommendations). Genotype data (chip VCF, GRCh37,
e.g. 736 individuals) is available from the project (INBIOMEDIC / UTEC, H. Guio).
Feed the chip VCF as a PLINK set or VCF to step 1/2; or use the published
diplotype→phenotype table directly as the step-3 output.

## Corpas family — European — family genome

The Corpas family genome is a publicly described personal-genomics dataset (the
"Corpasome" / Corpas family genome project). VCFs are available from the project
(see the Corpas family genome publications). Feed the family VCF to PyPGx (step 2).

## Reference data (fetched by the pipeline, not stored here)

- PyPGx phasing panels: `git clone https://github.com/sbslee/pypgx-bundle` (GRCh37).
- No GRCh38 reference or liftover chain is needed (GRCh37-native calling).
