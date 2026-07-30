#!/usr/bin/env python3
"""
Flatten PyPGx per-gene results archives into a (sample, gene, diplotype) table.

08_caller_truth_eval.py compares called diplotypes against an external truth set
and is deliberately source-agnostic about where the calls came from. PyPGx emits
one zipped results archive per gene, keyed by sample; this turns that into the
one table the evaluator reads.

Two behaviours matter more than the parsing, and both exist because the failure
mode this project keeps meeting is an evaluation harness that quietly reports a
cleaner number than the data supports:

1. A gene whose archive is missing or unreadable is REPORTED, never skipped.
   Silently dropping it would shrink the denominator and raise concordance.
2. Indeterminate and empty genotypes are dropped as abstentions and counted.
   Scoring a refusal to answer against truth manufactures disagreements.

USAGE
    python real-genome-arm/scripts/13_calls_to_sample_table.py \
        --calls real-genome-arm/work/getrm-1000g/calls \
        --out   real-genome-arm/work/getrm-1000g/pypgx_calls.tsv
"""
from __future__ import annotations

import argparse
import csv
import io
import sys
import zipfile
from pathlib import Path

# PyPGx writes the called diplotype under "Genotype"; the sample id is the
# unnamed index column.
GENOTYPE_COL = "Genotype"
NO_CALL = {"", ".", "indeterminate", "none", "nan"}


def read_gene_archive(archive: Path) -> list[dict]:
    """Read one results.zip and return its rows as dicts.

    Raises on anything malformed. The caller decides what to do with that;
    this function never returns a partial result silently.
    """
    with zipfile.ZipFile(archive) as z:
        names = [n for n in z.namelist() if n.endswith("data.tsv")]
        if not names:
            raise ValueError("no data.tsv in archive")
        text = z.read(names[0]).decode("utf-8", "replace")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError("empty data.tsv")
    return list(reader)


def collect(calls_dir: Path) -> tuple[list[dict], list[str]]:
    """Walk a PyPGx calls directory, returning (rows, problems).

    `problems` is a human-readable list of every gene that could not be read.
    It is returned rather than logged so the caller cannot ignore it by accident.
    """
    rows: list[dict] = []
    problems: list[str] = []

    for gene_dir in sorted(p for p in Path(calls_dir).iterdir() if p.is_dir()):
        gene = gene_dir.name
        archive = gene_dir / "results.zip"
        if not archive.exists():
            problems.append(f"{gene}: no results.zip (gene was not called)")
            continue
        try:
            records = read_gene_archive(archive)
        except Exception as exc:  # noqa: BLE001 - the reason is reported verbatim
            problems.append(f"{gene}: unreadable results.zip ({exc})")
            continue

        if not records:
            problems.append(f"{gene}: results.zip contained no rows")
            continue

        sample_col = list(records[0].keys())[0]
        if GENOTYPE_COL not in records[0]:
            problems.append(f"{gene}: no '{GENOTYPE_COL}' column")
            continue

        for rec in records:
            sample = (rec.get(sample_col) or "").strip()
            diplotype = (rec.get(GENOTYPE_COL) or "").strip()
            if not sample:
                continue
            if diplotype.lower() in NO_CALL:
                continue  # abstention, not a call
            rows.append({"sample": sample, "gene": gene, "diplotype": diplotype})

    return rows, problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls", required=True, type=Path,
                    help="PyPGx calls directory, one subdirectory per gene")
    ap.add_argument("--out", required=True, type=Path,
                    help="output TSV (sample, gene, diplotype)")
    args = ap.parse_args(argv)

    rows, problems = collect(args.calls)

    for p in problems:
        print(f"  NOT CALLED  {p}", file=sys.stderr)

    if not rows:
        # Refuse rather than write an empty table. An evaluator handed an empty
        # call set reports a vacuous result, which is worse than an error.
        print("No calls extracted. Refusing to write an empty table.",
              file=sys.stderr)
        raise SystemExit(2)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "gene", "diplotype"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    genes = sorted({r["gene"] for r in rows})
    samples = {r["sample"] for r in rows}
    print(f"wrote {len(rows)} calls: {len(samples)} samples, {len(genes)} genes")
    print(f"genes: {', '.join(genes)}")
    if problems:
        print(f"{len(problems)} gene(s) produced no calls; see above.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
