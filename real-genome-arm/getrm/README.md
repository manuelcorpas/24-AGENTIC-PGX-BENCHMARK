# GeT-RM caller-truth set (revision item N5)

## Why this exists

N0 validated the executed skill's **mapping** step against PyPGx. It did not
validate the **calling** step, because PyPGx both called the diplotype and
supplied the phenotype it was scored against. Input-call error is therefore
unobservable in N0 rather than eliminated, and the plan of revision says so
explicitly to the editor.

GeT-RM supplies external truth for the calling step. With it, the claim "100%
correctness among emitted answers" becomes genuinely falsifiable: a call that
disagrees with the GeT-RM consensus is an end-to-end error attributable to
input interpretation.

## The source to use (updated 2026-07-27)

Use the **GeT-RM Consolidated PGx and HLA Table**. It is a single spreadsheet
carrying consensus genotypes for **363 samples across 34 genes/loci**, pooled
from nine GeT-RM PGx and HLA studies. This supersedes the earlier plan of
reconstructing truth from the 2010 and 2016 supplementary tables gene by gene.

Cite:

> Scheinfeldt L, Kusic D, Gaedigk A, Turner AJ, Moyer AM, Pratt VM, Kalman LV.
> New Resources to Identify Characterized DNA Reference Materials for
> Pharmacogenetic (PGx) and Human Leukocyte Antigen (HLA) Testing: The Genetic
> Testing Reference Material (GeT-RM) Program PGx Search Tool and GeT-RM
> Consolidated PGx and HLA Table. *J Mol Diagn* 2025;27(6):457-464.
> doi:10.1016/j.jmoldx.2025.02.008

Landing pages:

- Table (Excel): https://www.cdc.gov/lab-quality/php/get-rm/reference-materials.html
- Search tool: https://www.coriell.org/GetRM/PGxSearch

## Acquisition, and a correction

**This section previously said the download had to be manual because cdc.gov
returns HTTP 403 to non-browser clients. That was wrong.** cdc.gov 403s on
requests missing a full browser header set, not on automated clients as such.
With `Accept`, `Accept-Language`, `Sec-Fetch-*` and `Upgrade-Insecure-Requests`
it returns 200 for both the landing page and the `.xlsx` asset.

The original conclusion came from two failed fetches with a bare user-agent. Two
negative results from one method were generalised into a property of the server.
That is the same error class this paper is about, so it is recorded here rather
than quietly amended.

The download is therefore automated and reproducible:

```bash
bash real-genome-arm/scripts/10b_fetch_getrm_consolidated.sh
```

It pins the asset URL, refuses to proceed if CDC serves an error page instead of
a spreadsheet, and checks the SHA-256 against the copy retrieved 2026-07-27
(`c66174be…3be109`). A changed digest means CDC revised the table, which is
information rather than an error: re-ingest, re-evaluate, and re-record it.

The CDC asset path is dated and may move when the table is revised. If it 404s,
take the current link from the landing page above.

The Coriell PGx Search tool remains a form-driven page with no bulk export, so
the consolidated spreadsheet is the right source.

## Acquisition status, 2026-07-27

Truth set **acquired and ingested**:

    3,554 genotypes over 323 samples and 35 genes
    includes HLA-A and HLA-B, which carry four of the benchmark's
    lethal-class cases

Eleven columns yielded no star-allele call and are excluded, and the ingester
prints them rather than dropping them silently: the rsID-specific SNP columns
(`CYP2C Cluster … rs12777823`, two `GGCX` columns, two `VKORC1` columns, plain
`VKORC1`), `GSTT1` and `SLCO2B1` (reported in a parenthesised or `WT/WT`
vocabulary rather than star alleles), and three ENA accession columns.

**The evaluation itself is not yet run**, and this is the honest remaining gap.
It needs PyPGx calls for GeT-RM samples, and only **2 of the 323** GeT-RM
samples (HG01680, HG01697) are in the locally called IBS cohort. Concordance on
two samples would be an underpowered number, not a validation. Closing it
requires calling PyPGx on the GeT-RM samples, whose 1000 Genomes and ENA data
are not held locally. That is a download and compute job, not an analysis
decision.

## What the ingester refuses to do

It does not invent, impute or carry forward a genotype. Blank cells, `Not
tested`, `ND`, `N/A` and free-text cells are dropped rather than admitted, and a
table that yields no usable genotypes raises rather than writing an empty file.
The evaluator likewise refuses to run on an empty or malformed truth set instead
of reporting a vacuous 100% concordance.

`tests/test_getrm_ingest.py` covers the shape handling with **synthetic
fixtures**. Those fixtures are not a truth set and must never be used as one.

## Expected schema

`getrm_consensus.tsv`, tab separated:

    sample	gene	diplotype
    NA12878	CYP2C19	*1/*2

The evaluator is source-agnostic: any `(sample, gene, diplotype)` table works,
so the harness does not rot if the truth source changes.
