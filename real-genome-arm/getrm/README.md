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

## Acquisition is manual, and this is not an oversight

`www.cdc.gov` returns **HTTP 403 to non-browser clients** (verified 2026-07-27
with both a plain fetch and a browser user-agent), and the Coriell PGx Search
tool is a form-driven ASP.NET page with no bulk export. There is no stable
machine-readable endpoint to pin, so an automated fetch would be a scraper that
breaks silently and produces a truth set nobody can reconstruct.

The download is therefore a documented human step. Everything after it is
automated.

## Steps

1. Open the CDC page above in a browser and download the consolidated table.
2. Save it here, for example `real-genome-arm/getrm/GeT-RM_consolidated.xlsx`.
   A CSV export of the first sheet works equally well.
3. Convert it to the caller-truth schema:

   ```bash
   python real-genome-arm/scripts/11_ingest_getrm_consolidated.py \
       --input real-genome-arm/getrm/GeT-RM_consolidated.xlsx \
       --out   real-genome-arm/getrm/getrm_consensus.tsv
   # optional: restrict to the genes the benchmark calls
   #   --genes CYP2C19 CYP2D6 CYP2C9 TPMT NUDT15 DPYD SLCO1B1 CYP3A5 CYP2B6
   ```

   Reading `.xlsx` needs `openpyxl`; a CSV export needs nothing.

4. Evaluate:

   ```bash
   python real-genome-arm/scripts/08_caller_truth_eval.py \
       --truth  real-genome-arm/getrm/getrm_consensus.tsv \
       --called <pypgx calls tsv>
   ```

The join key is the Coriell identifier (`NA12878` style), which is also how the
1000 Genomes calls from `scripts/09_call_1000g_eur.sh` and
`scripts/10_call_1000g_ibs.sh` are named, so samples overlap without a mapping
file.

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
