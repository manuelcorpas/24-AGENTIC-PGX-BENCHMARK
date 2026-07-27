#!/usr/bin/env python3
"""
Convert the GeT-RM Consolidated PGx and HLA Table into the caller-truth schema.

WHY THIS EXISTS
The consolidated table is a WIDE sheet, one row per Coriell sample and one
column per gene, published by CDC GeT-RM as a single spreadsheet covering 363
samples and 34 genes/loci (Scheinfeldt et al., J Mol Diagn 2025;27:457-464,
doi:10.1016/j.jmoldx.2025.02.008). The caller-truth evaluator
(08_caller_truth_eval.py) wants a LONG (sample, gene, diplotype) table. This
script is the bridge, so the only human step left is the download itself.

WHY THE DOWNLOAD IS STILL MANUAL
www.cdc.gov returns HTTP 403 to non-browser clients (verified 2026-07-27), and
the Coriell PGx Search tool is a form-driven ASP.NET page with no bulk export.
An automated fetch would be a scraper that breaks silently and yields a truth
set nobody can reconstruct. See real-genome-arm/getrm/README.md.

WHAT IT REFUSES TO DO
It does not invent, impute or carry forward a genotype. Blank cells, "Not
tested", "ND", "N/A" and free text are dropped rather than admitted, and a
table that yields no usable genotypes is an error rather than an empty truth
set. A truth set that silently loses rows produces a concordance number that
looks fine and means nothing.

USAGE
    python real-genome-arm/scripts/11_ingest_getrm_consolidated.py \
        --input real-genome-arm/getrm/GeT-RM_consolidated.xlsx \
        --out   real-genome-arm/getrm/getrm_consensus.tsv
    # optionally restrict to the benchmark's genes:
    #   --genes CYP2C19 CYP2D6 TPMT NUDT15 DPYD
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
GETRM = HERE.parent / "getrm"

# Headings the sample column has appeared under across GeT-RM releases.
SAMPLE_HEADINGS = ("coriell", "sample", "cell line", "dna", "nigms", "nhgri")

# Cells that mean "no genotype here". Matched case-insensitively after strip.
NON_GENOTYPE = {"", "-", "--", "nd", "n/a", "na", "none", "not tested",
                "not determined", "not available", "unknown", "nt"}

# Columns that carry provenance or data-availability flags, not genotypes.
METADATA_HINTS = ("reference", "bam", "fastq", "1000 genomes", "note", "comment",
                  "source", "study", "method", "population", "sex", "ethnicity")

# A diplotype is two alleles round a slash, or a single named allele. Star
# alleles, HLA names and named haplotypes all pass; prose does not.
ALLELE = r"(?:\*[\w\.]+|[A-Za-z0-9][\w\*\.\:\-]*)"
DIPLOTYPE_RE = re.compile(rf"^{ALLELE}\s*/\s*{ALLELE}$")
SINGLE_RE = re.compile(rf"^{ALLELE}$")


def read_table(path: Path) -> list[dict]:
    """Read the consolidated table from .xlsx, .xls or a CSV/TSV export."""
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xlsm", ".xls"):
        try:
            from openpyxl import load_workbook
        except ImportError:
            raise SystemExit(
                "reading .xlsx needs openpyxl (pip install openpyxl), or export "
                "the sheet to CSV and pass that instead")
        wb = load_workbook(path, read_only=True, data_only=True)
        ws = wb[wb.sheetnames[0]]
        rows = [[("" if c is None else str(c)) for c in r]
                for r in ws.iter_rows(values_only=True)]
        wb.close()
        if not rows:
            raise SystemExit(f"{path} has no rows")
        header, body = rows[0], rows[1:]
        return [dict(zip(header, r)) for r in body]

    delim = "\t" if suffix in (".tsv", ".tab") else ","
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def find_sample_column(fieldnames) -> str:
    for name in fieldnames:
        if name and any(h in name.lower() for h in SAMPLE_HEADINGS):
            return name
    raise SystemExit(
        "could not identify the sample-identifier column in "
        f"{list(fieldnames)!r}. Rename it to 'Coriell #' or pass a CSV whose "
        "first column is the Coriell identifier.")


def is_metadata_column(name: str) -> bool:
    low = (name or "").lower()
    return any(h in low for h in METADATA_HINTS)


def clean_diplotype(value: str) -> str | None:
    """Return a normalised diplotype, or None if the cell is not a genotype."""
    v = (value or "").strip()
    if v.lower() in NON_GENOTYPE:
        return None
    # GeT-RM sometimes footnotes a call; drop trailing footnote markers.
    v = re.sub(r"\s*\[[^\]]*\]$", "", v).strip()
    v = re.sub(r"\s+", " ", v)
    if DIPLOTYPE_RE.match(v):
        a, b = (p.strip() for p in v.split("/", 1))
        return f"{a}/{b}"
    if SINGLE_RE.match(v) and ("*" in v or ":" in v):
        return v
    return None


def to_long(records: list[dict], genes: set[str] | None = None) -> list[dict]:
    if not records:
        return []
    fieldnames = list(records[0].keys())
    sample_col = find_sample_column(fieldnames)
    gene_cols = [c for c in fieldnames
                 if c and c != sample_col and not is_metadata_column(c)]

    out: list[dict] = []
    for rec in records:
        sample = (rec.get(sample_col) or "").strip()
        if not sample:
            continue
        for col in gene_cols:
            gene = col.strip()
            if genes is not None and gene not in genes:
                continue
            dip = clean_diplotype(rec.get(col, ""))
            if dip is None:
                continue
            out.append({"sample": sample, "gene": gene, "diplotype": dip})
    return out


def build(src: Path, out: Path | None = None, genes: set[str] | None = None) -> list[dict]:
    rows = to_long(read_table(src), genes)
    if not rows:
        raise SystemExit(
            f"{src} yielded no usable genotypes. Refusing to write an empty "
            "truth set: the evaluator would report a vacuous 100% concordance.")
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["sample", "gene", "diplotype"],
                               delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--input", type=Path, required=True,
                    help="GeT-RM consolidated table (.xlsx or a CSV export)")
    ap.add_argument("--out", type=Path, default=GETRM / "getrm_consensus.tsv")
    ap.add_argument("--genes", nargs="*", default=None,
                    help="restrict to these gene columns")
    args = ap.parse_args(argv)

    if not args.input.exists():
        sys.stderr.write(
            f"no input at {args.input}\n"
            "The consolidated table is a manual download; see "
            "real-genome-arm/getrm/README.md\n")
        return 1

    rows = build(args.input, args.out, set(args.genes) if args.genes else None)
    samples = len({r["sample"] for r in rows})
    genes = sorted({r["gene"] for r in rows})
    print(f"wrote {args.out}")
    print(f"  {len(rows)} genotypes over {samples} samples and {len(genes)} genes")
    print(f"  genes: {', '.join(genes)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
