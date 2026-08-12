#!/usr/bin/env python3
"""
Step 17 — aggregate the Corpas family WGS calls produced by step 16.

Reads   <calls>/<arm>/<sample>/<gene>/results.zip
Writes  a tidy (cohort, gene, diplotype, phenotype, n_carriers) table per arm,
        in the same schema as data/n5_four_cohorts.tsv, plus a QC report.

Three checks are run and printed, because this family has no external truth set
and these are the only independent evidence that the calls are coherent:

  1. Mendelian consistency of each child's diplotype against both parents.
     An inconsistency is a caller error, not a biological finding, and is
     reported as such rather than repaired or dropped.
  2. NGS versus chip pipeline concordance on identical input. This is the
     assay/calling-mode effect, measured rather than asserted.
  3. WGS versus the superseded 23andMe array calls, where those are supplied.

A gene that could not be called yields no verdict at all. It is never scored as
consistent, because "we could not check" and "we checked and it passed" are
different claims and only one of them is true.

Usage:
  17_aggregate_family_wgs.py <calls_dir> <out_dir> [--array-tsv PATH]
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

# Pedigree. Established by kinship on chr22, which corrected an earlier
# container list: PT00001A is unrelated to this family and PT00002A is the son.
# PT00010A (aunt) has only a 4.4x GRCh38 build and is not called.
FATHER = "PT00007A"
MOTHER = "PT00008A"
CHILDREN = ("PT00002A", "PT00009A")


def split_diplotype(dip: str):
    """('*1', '*80+*28') from '*1/*80+*28'. None if not a diplotype.

    Splits on the FIRST separator only: an allele may itself contain '+',
    parentheses or HGVS punctuation, e.g. 'c.85T>C (*9A)'.
    """
    if dip is None:
        return None
    dip = dip.strip()
    if not dip or "/" not in dip:
        return None
    a, b = dip.split("/", 1)
    a, b = a.strip(), b.strip()
    if not a or not b:
        return None
    return a, b


def mendelian_ok(child: str, father: str, mother: str):
    """True/False, or None when any of the three could not be parsed."""
    c, f, m = split_diplotype(child), split_diplotype(father), split_diplotype(mother)
    if c is None or f is None or m is None:
        return None
    # One child allele must come from each parent, in either orientation.
    return (c[0] in f and c[1] in m) or (c[1] in f and c[0] in m)


def read_calls(calls_dir: Path, arm: str):
    """{sample: {gene: (genotype, phenotype)}} for one arm."""
    out: dict[str, dict[str, tuple[str, str]]] = {}
    arm_dir = calls_dir / arm
    if not arm_dir.is_dir():
        return out
    for sample_dir in sorted(p for p in arm_dir.iterdir() if p.is_dir()):
        per_gene: dict[str, tuple[str, str]] = {}
        for gene_dir in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
            zpath = gene_dir / "results.zip"
            if not zpath.exists():
                continue
            with zipfile.ZipFile(zpath) as zf:
                names = [n for n in zf.namelist() if n.endswith("data.tsv")]
                if not names:
                    continue
                reader = csv.DictReader(
                    io.TextIOWrapper(io.BytesIO(zf.read(names[0])), encoding="utf-8"),
                    delimiter="\t",
                )
                for row in reader:
                    per_gene[gene_dir.name] = (
                        (row.get("Genotype") or "").strip(),
                        (row.get("Phenotype") or "").strip(),
                    )
        out[sample_dir.name] = per_gene
    return out


def tidy_rows(calls: dict, cohort: str):
    """(cohort, gene, diplotype, phenotype, n_carriers) per distinct state."""
    counts: dict[tuple[str, str], int] = defaultdict(int)
    phen: dict[tuple[str, str], str] = {}
    for _sample, per_gene in calls.items():
        for gene, (dip, ph) in per_gene.items():
            if not dip:
                continue
            counts[(gene, dip)] += 1
            phen[(gene, dip)] = ph
    return [
        (cohort, gene, dip, phen[(gene, dip)], n)
        for (gene, dip), n in sorted(counts.items())
    ]


def concordance(a: dict, b: dict):
    """(compared, agreeing, [(sample, gene, a_call, b_call), ...]) on shared pairs."""
    compared = agree = 0
    diffs = []
    for sample, per_gene in sorted(a.items()):
        for gene, (dip_a, _ph) in sorted(per_gene.items()):
            other = b.get(sample, {}).get(gene)
            if other is None:
                continue
            compared += 1
            if dip_a == other[0]:
                agree += 1
            else:
                diffs.append((sample, gene, dip_a, other[0]))
    return compared, agree, diffs


def mendelian_report(calls: dict):
    """(checked, inconsistent, [(child, gene, child_dip, dad_dip, mum_dip)], unchecked)."""
    checked = 0
    unchecked = 0
    bad = []
    for child in CHILDREN:
        for gene in sorted(calls.get(child, {})):
            cd = calls[child][gene][0]
            fd = calls.get(FATHER, {}).get(gene, ("", ""))[0]
            md = calls.get(MOTHER, {}).get(gene, ("", ""))[0]
            verdict = mendelian_ok(cd, fd, md)
            if verdict is None:
                unchecked += 1
                continue
            checked += 1
            if not verdict:
                bad.append((child, gene, cd, fd, md))
    return checked, len(bad), bad, unchecked


def read_array_tsv(path: Path, cohort_label: str = "CorpasFamily"):
    """Distinct (gene, diplotype) states in the superseded array table."""
    states = set()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("cohort") == cohort_label:
                states.add((row["gene"], row["diplotype"]))
    return states


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("calls_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--array-tsv", default=None)
    args = ap.parse_args()

    calls_dir = Path(args.calls_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    arms = {arm: read_calls(calls_dir, arm) for arm in ("ngs", "chip")}
    report: list[str] = []

    for arm, calls in arms.items():
        if not calls:
            report.append(f"{arm}: NO CALLS FOUND")
            continue
        label = "CorpasFamilyWGS" if arm == "ngs" else "CorpasFamilyWGSchip"
        rows = tidy_rows(calls, label)
        out_tsv = out_dir / f"family_wgs_{arm}.tsv"
        with open(out_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["cohort", "gene", "diplotype", "phenotype", "n_carriers"])
            w.writerows(rows)
        genes = {g for per in calls.values() for g in per}
        called = sum(1 for per in calls.values() for d, _ in per.values() if d)
        report.append(
            f"{arm}: {len(calls)} samples x {len(genes)} genes, "
            f"{called} diplotype calls, {len(rows)} distinct (gene, diplotype) states "
            f"-> {out_tsv.name}"
        )
        checked, n_bad, bad, unchecked = mendelian_report(calls)
        report.append(
            f"{arm}: Mendelian {n_bad} inconsistent of {checked} checked "
            f"({unchecked} not checkable)"
        )
        for child, gene, cd, fd, md in bad:
            report.append(f"    {child} {gene}: child {cd} | father {fd} | mother {md}")

    if arms["ngs"] and arms["chip"]:
        n, agree, diffs = concordance(arms["ngs"], arms["chip"])
        pct = f"{100.0 * agree / n:.1f}%" if n else "n/a"
        report.append(f"ngs vs chip: {agree}/{n} identical diplotypes ({pct})")
        for sample, gene, a, b in diffs:
            report.append(f"    {sample} {gene}: ngs {a} | chip {b}")

    if args.array_tsv:
        array_states = read_array_tsv(Path(args.array_tsv))
        wgs_states = {(g, d) for _c, g, d, _p, _n in tidy_rows(arms["ngs"], "x")}
        report.append(
            f"array vs WGS distinct states: array {len(array_states)}, "
            f"WGS {len(wgs_states)}, shared {len(array_states & wgs_states)}, "
            f"WGS-only {len(wgs_states - array_states)}, "
            f"array-only {len(array_states - wgs_states)}"
        )

    text = "\n".join(report)
    (out_dir / "family_wgs_qc.txt").write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
