#!/usr/bin/env python3
"""
Step 18 — rename sample identifiers inside a PyPGx archive.

Why this exists: the Corpasome BAM read groups carry sequencing run IDs
(e.g. CL100084707_L1_HUMbgjRAAAB-549) while the VCFs used for calling carry
container IDs (PT00007A). PyPGx matches the depth-of-coverage and control-
statistics archives to the variant VCF BY SAMPLE NAME, so the two must agree
or copy-number analysis silently applies to nobody.

The rename is deliberately strict. It refuses to write unless every name it was
asked to change is present and every name present is accounted for, because a
partially applied rename would attach one person's read depth to another
person's variants. That failure would not raise; it would produce a plausible
diplotype for the wrong individual.

Usage:
  18_rename_archive_samples.py <in.zip> <out.zip> OLD=NEW [OLD=NEW ...]
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path


def parse_mapping(pairs):
    mapping = {}
    for p in pairs:
        if "=" not in p:
            raise ValueError(f"mapping must be OLD=NEW, got {p!r}")
        old, new = p.split("=", 1)
        old, new = old.strip(), new.strip()
        if not old or not new:
            raise ValueError(f"mapping must be OLD=NEW, got {p!r}")
        # A shell that does not word-split (zsh, unquoted "$VAR") delivers every
        # pair as ONE argument. Splitting on the first '=' then yields a target
        # name containing the rest of the mapping, which is silently wrong.
        for part, label in ((old, "source"), (new, "target")):
            if any(ch.isspace() for ch in part):
                raise ValueError(
                    f"{label} name contains whitespace: {part!r}. "
                    "Pass each OLD=NEW as a separate argument."
                )
            if part.count("=") or "," in part:
                raise ValueError(f"{label} name contains a separator: {part!r}")
        if old in mapping:
            raise ValueError(f"duplicate source name {old!r}")
        mapping[old] = new
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("mapping is not one-to-one: duplicate target names")
    return mapping


def rename_header(header: str, mapping: dict, fixed: tuple) -> str:
    """Rewrite a TSV header, renaming sample columns only.

    `fixed` names leading non-sample columns (e.g. Chromosome, Position) that
    must be left alone. Raises if any mapped name is missing, or if a sample
    column has no mapping.
    """
    cols = header.rstrip("\n").split("\t")
    lead = [c for c in cols if c in fixed]
    samples = [c for c in cols if c not in fixed]
    missing = set(mapping) - set(samples)
    if missing:
        raise ValueError(f"mapping names absent from archive: {sorted(missing)}")
    unmapped = set(samples) - set(mapping)
    if unmapped:
        raise ValueError(f"sample columns with no mapping: {sorted(unmapped)}")
    out = [c if c in fixed else mapping[c] for c in cols]
    assert len(out) == len(cols) == len(lead) + len(samples)
    return "\t".join(out) + "\n"


FIXED = ("Chromosome", "Position", "Name", "")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("mapping", nargs="+")
    args = ap.parse_args()

    mapping = parse_mapping(args.mapping)
    src, dst = Path(args.src), Path(args.dst)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        with zipfile.ZipFile(src) as zf:
            names = zf.namelist()
            zf.extractall(tmp)
        data = [n for n in names if n.endswith("data.tsv")]
        if len(data) != 1:
            raise SystemExit(f"expected exactly one data.tsv, found {data}")
        target = tmp / data[0]
        lines = target.read_text().splitlines(keepends=True)
        if not lines:
            raise SystemExit("data.tsv is empty")
        first = lines[0]
        # A SampleTable puts samples in the FIRST COLUMN, one per row;
        # a CovFrame puts them in the header. Detect which.
        header_cols = first.rstrip("\n").split("\t")
        if set(mapping) & set(header_cols):
            lines[0] = rename_header(first, mapping, FIXED)
            renamed = len(mapping)
        else:
            renamed = 0
            for i in range(1, len(lines)):
                cells = lines[i].rstrip("\n").split("\t")
                if cells and cells[0] in mapping:
                    cells[0] = mapping[cells[0]]
                    lines[i] = "\t".join(cells) + "\n"
                    renamed += 1
            if renamed != len(mapping):
                raise SystemExit(
                    f"renamed {renamed} row labels but mapping has {len(mapping)}"
                )
        target.write_text("".join(lines))

        if dst.exists():
            dst.unlink()
        with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zf:
            for n in names:
                zf.write(tmp / n, n)

    print(f"renamed {renamed} sample identifiers -> {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
