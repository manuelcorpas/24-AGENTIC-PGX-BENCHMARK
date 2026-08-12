#!/usr/bin/env python3
"""Supplementary Table S1 (headline metrics, both scorers) and Table S3
(per-gene lethal-class errors), generated rather than transcribed.

WHY THIS EXISTS
Two supplementary tables had no generator. Table S1 still carried the earlier
three-arm protocol's denominators (8,738 and 1,096) and an answer-supplied
control block, none of which appears in the matched five-cell manuscript, so
the Results' claim that "every number is produced under both the locked
baseline scorer and a frozen clinical-equivalence scorer" had no supporting
table for the five cells. Table S3's per-gene lethal counts were correct but
unregenerable, which is the state a number is in just before it goes stale.

WHAT IT DOES NOT DO
It does not score anything. Scoring is done once, by the cell-blind scorer in
61-rescore-matched.py, and re-implementing it here would create a second
scorer that could silently disagree with the one the paper describes. Table S1
is read from that scorer's own report under both scorer sections; Table S3
aggregates its already-scored rows. The per-gene lethal totals are asserted
against the report's per-cell lethal counts, so the two tables cannot drift.

USAGE
    python code/89-supp-headline-and-per-gene.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
REPORT = DATA / "v3_five_cell_live_report.txt"
SCORED = DATA / "v3_matched_scored_rows_all5.json"

# Display order and the labels used in the manuscript, so the supplement and
# the main text name the same five cells in the same sequence.
CELLS = [
    ("free_generation", "free generation"),
    ("rag_generation", "RAG generation"),
    ("rag_execution", "RAG call + authored execution"),
    ("skill_generation", "authored rules, generated"),
    ("skill_execution", "authored rules, executed"),
]

LETHAL_GENE_ORDER = [
    "HLA-A*31:01", "HLA-B*15:02", "HLA-B*57:01", "HLA-B*58:01",
    "DPYD", "TPMT", "NUDT15", "G6PD", "CYP2D6", "MT-RNR1", "RYR1",
]


def parse_report(text: str) -> dict[str, dict[str, dict]]:
    """{scorer: {cell: metrics}} from the scorer's own report."""
    out: dict[str, dict[str, dict]] = {}
    scorer = None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("## "):
            scorer = s[3:].strip()
            out[scorer] = {}
            continue
        if not scorer or "n=" not in s:
            continue
        cell = s.split()[0]
        fields = {}
        for token in s.split():
            if "=" in token:
                key, value = token.split("=", 1)
                fields[key] = value
        errors = None
        if "(" in s and "errors)" in s:
            errors = int(s.split("(")[1].split()[0])
        out[scorer][cell] = {
            "n": int(fields["n"]),
            "coverage": float(fields["coverage"]),
            "A1": float(fields["A1"]),
            "A2": float(fields["A2"]),
            "lethal_accuracy": float(fields["lethal"]),
            "lethal_errors": errors,
            "parse_fail": int(fields.get("parse_fail", 0)),
        }
    return out


def pct(value: float) -> str:
    return f"{value * 100:.1f}%"


def pct2(value: float) -> str:
    """Coverage to two decimals.

    Authored-rule execution covers 0.999621, which is 99.96% in main-text
    Table 1 and rounds to 100.0% at one decimal. A supplementary table saying
    100.0% beside a main table saying 99.96% is the same number contradicting
    itself in print, which is the most damaging form of a rounding choice.
    """
    return f"{value * 100:.2f}%"


def build_headline(report: dict) -> list[str]:
    baseline = report["baseline"]
    equivalence = report["clinical_equivalence"]
    lines = [
        "Table S1. Headline metrics for the five matched cells under both scorers.",
        "",
        f"{'cell':<32}{'n':>7}{'coverage':>10}"
        f"{'A1 base':>9}{'A1 equiv':>10}{'A2 base':>9}{'A2 equiv':>10}"
        f"{'lethal errors':>15}{'parse fail':>12}",
    ]
    for key, label in CELLS:
        b, e = baseline[key], equivalence[key]
        if b["A2"] != e["A2"] or b["lethal_errors"] != e["lethal_errors"]:
            sys.stderr.write(
                f"{key}: the equivalence layer changed A2 or a lethal count; the "
                "supplement states it modifies phenotype identification only\n")
            raise SystemExit(1)
        lines.append(
            f"{label:<32}{b['n']:>7}{pct2(b['coverage']):>10}"
            f"{pct(b['A1']):>9}{pct(e['A1']):>10}"
            f"{pct(b['A2']):>9}{pct(e['A2']):>10}"
            f"{b['lethal_errors']:>15}{b['parse_fail']:>12}")
    return lines


def build_per_gene(rows: list[dict], report: dict) -> list[str]:
    lethal = [r for r in rows if r["lethal_action"] is not None]
    counts: dict[str, dict[str, list[int]]] = {}
    for row in lethal:
        bucket = counts.setdefault(row["gene"], {})
        errors, total = bucket.setdefault(row["cell"], [0, 0])
        bucket[row["cell"]] = [errors + (0 if row["lethal_action"] else 1),
                               total + 1]

    genes = [g for g in LETHAL_GENE_ORDER if g in counts]
    unexpected = sorted(set(counts) - set(genes))
    if unexpected:
        sys.stderr.write(f"lethal-class genes not in the display order: {unexpected}\n")
        raise SystemExit(1)

    lines = [
        "Table S3. Per-gene lethal-class errors across the five matched cells.",
        "",
        f"{'gene':<14}" + "".join(f"{label:>32}" for _, label in CELLS),
    ]
    totals = {key: [0, 0] for key, _ in CELLS}
    for gene in genes:
        cells = []
        for key, _ in CELLS:
            errors, total = counts[gene].get(key, [0, 0])
            totals[key][0] += errors
            totals[key][1] += total
            cells.append(f"{errors}/{total}")
        lines.append(f"{gene:<14}" + "".join(f"{c:>32}" for c in cells))
    lines.append(f"{'all loci':<14}"
                 + "".join(f"{f'{totals[k][0]}/{totals[k][1]}':>32}"
                           for k, _ in CELLS))

    # Bind to the scorer's own per-cell counts. If these disagree, one of the
    # two tables is describing a different run and neither should be published.
    for key, _ in CELLS:
        stated = report["baseline"][key]["lethal_errors"]
        if totals[key][0] != stated:
            sys.stderr.write(
                f"{key}: per-gene lethal errors sum to {totals[key][0]}, the "
                f"scorer's report says {stated}; refusing to write\n")
            raise SystemExit(1)
        if totals[key][1] != 336:
            sys.stderr.write(
                f"{key}: {totals[key][1]} lethal-class evaluations, expected 336\n")
            raise SystemExit(1)
    return lines


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DATA / "v3_supp_tables_s1_s3.txt")
    args = ap.parse_args(argv)

    for path in (REPORT, SCORED):
        if not path.exists():
            sys.stderr.write(
                f"missing {path.name}; raw rows are gitignored and come from the "
                "Zenodo version cited in the manuscript. INCOMPLETE, nothing written\n")
            return 2

    report = parse_report(REPORT.read_text())
    rows = json.loads(SCORED.read_text())

    lines = build_headline(report) + ["", ""] + build_per_gene(rows, report)
    text = "\n".join(lines) + "\n"
    args.out.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
