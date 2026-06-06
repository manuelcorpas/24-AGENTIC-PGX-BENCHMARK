#!/usr/bin/env python3
"""
A2 regression classifier, locked to the 3-of-3 replicate three-arm dataset.

For each (model, tc, pop) cell, the script:
  1. Aggregates A2 across the three replicates per condition (no_spec, cpic_rag).
  2. Identifies "A2 regression" cells where cpic_rag mean A2 < no_spec mean A2.
  3. Classifies each regression into one of four categories using a single
     representative cpic_rag parsed DRUG field (the first parsed run for the
     cell, matching the original 2-run-snapshot heuristic).

Categories (matching the original CSV taxonomy):
  a_drug_substitution: queried_drug name is not present (case-insensitive
                       substring) in the cpic_rag parsed DRUG field. The model
                       answered for a different drug than was queried (the
                       classic RAG chunk-substitution mode).
  b_phrasing_equiv:    queried_drug IS present AND cpic_rag direction equals
                       gt_direction (scorer regex missed an equivalent phrasing).
  c_wrong_direction:   queried_drug IS present AND cpic_rag direction is a
                       distinct, non-UNCLEAR direction that conflicts with gt.
  d_other:             queried_drug present but direction is UNCLEAR, or any
                       other case the above three do not match.

Direction extraction (gt_drug and cpic_rag DRUG) uses the same regex family as
the rigorous rescorer in 10-rescore-v3.py: AVOID, ALT, REDUCE, STANDARD, UNCLEAR.

Inputs:
  RESULTS/v3_raw_rescored_three_arm.json   (merged 3-run dataset)
  SPECS/test_cases_v3.json                 (queried_drug and gt_drug per tc)

Outputs:
  RESULTS/v3_three_arm_a2_regression_classified.csv  (overwrites existing)
  RESULTS/v3_three_arm_a2_regression_summary.txt     (human-readable counts)

The previous CSV (2-run-snapshot snapshot) should be archived before this
script overwrites it. The script will refuse to overwrite if PRESERVE_OLD=1
is set in the environment and an "_pre3run.csv" sibling does not exist; in
that case it writes the sibling first.

Usage:
  python3 33-classify-a2-regressions.py             # rerun + write CSV + summary
  python3 33-classify-a2-regressions.py --self-test # data sanity check only
  python3 33-classify-a2-regressions.py --compare   # diff vs existing CSV
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
MERGED = BASE / "RESULTS" / "v3_raw_rescored_three_arm.json"
CASES_FILE = BASE / "SPECS" / "test_cases_v3.json"
OUT_CSV = BASE / "RESULTS" / "v3_three_arm_a2_regression_classified.csv"
OUT_SUMMARY = BASE / "RESULTS" / "v3_three_arm_a2_regression_summary.txt"

CSV_COLUMNS = [
    "category", "model", "tc", "pop",
    "queried_drug", "gt_drug", "cpic_rag_DRUG",
    "gt_direction", "cpic_rag_direction",
    "A2_no_spec_mean", "A2_cpic_rag_mean",
]


_HAS_AVOID = re.compile(r"\bavoid\b|\bcontraindicat\w*\b|\bdo not use\b", re.I)
_HAS_ALT = re.compile(
    r"\balternative\b|\bswitch to\b|\binstead of\b|\bdifferent (drug|agent|therapy)\b"
    r"|\bnot indicated\b|\buse combination therapy\b",
    re.I,
)
_HAS_REDUCE = re.compile(
    r"\breduc\w*\b|\blower\b|\bdecreas\w*\b|\b\d+\s*%\b|\bsmaller\b", re.I
)
_HAS_STANDARD = re.compile(
    r"\bstandard\s+(dose|dosing|use|starting dose)\b|\busual\s+dose\b"
    r"|\bnormal\s+dose\b|\bregular\s+dose\b|^[^:]*:\s*standard\b"
    r"|\b(?<!not\s)indicated\b",
    re.I,
)


def extract_direction(text: str) -> str:
    if not text:
        return "UNCLEAR"
    has_avoid = bool(_HAS_AVOID.search(text))
    has_alt = bool(_HAS_ALT.search(text))
    has_reduce = bool(_HAS_REDUCE.search(text))
    has_standard = bool(_HAS_STANDARD.search(text))

    if has_alt:
        return "ALT"
    if has_avoid:
        return "AVOID"
    if has_reduce:
        return "REDUCE"
    if has_standard:
        return "STANDARD"
    return "UNCLEAR"


def aggregate_a2_per_cell(rows: list[dict]) -> dict:
    by_cell: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for r in rows:
        if r.get("scores", {}).get("format_fail"):
            continue
        if r["cond"] not in ("no_spec", "cpic_rag"):
            continue
        a2 = r.get("scores", {}).get("A2")
        if a2 is None:
            continue
        key = (r["model"], r["tc"], r["pop"], r["cond"])
        by_cell[key].append(r)

    means: dict[tuple[str, str, str], dict] = {}
    for (model, tc, pop, cond), cell_rows in by_cell.items():
        scores = [r["scores"]["A2"] for r in cell_rows]
        mean = sum(scores) / len(scores) if scores else None
        mk = (model, tc, pop)
        means.setdefault(mk, {})[cond] = {
            "mean": mean,
            "n": len(cell_rows),
            "rows": cell_rows,
        }
    return means


def first_parsed_drug(cell_rows: list[dict]) -> str:
    cell_rows = sorted(cell_rows, key=lambda r: r.get("run", 0))
    for r in cell_rows:
        drug = (r.get("parsed") or {}).get("DRUG")
        if drug:
            return drug.strip()
    return ""


def classify(queried_drug: str, gt_drug: str, cpic_rag_drug: str) -> tuple[str, str, str]:
    gt_dir = extract_direction(gt_drug)
    rag_dir = extract_direction(cpic_rag_drug)

    queried_lc = (queried_drug or "").lower()
    rag_lc = (cpic_rag_drug or "").lower()

    if queried_lc and queried_lc not in rag_lc:
        return "a_drug_substitution", gt_dir, rag_dir

    if gt_dir != "UNCLEAR" and rag_dir == gt_dir:
        return "b_phrasing_equiv", gt_dir, rag_dir

    if rag_dir != "UNCLEAR" and rag_dir != gt_dir:
        return "c_wrong_direction", gt_dir, rag_dir

    return "d_other", gt_dir, rag_dir


def build_classified_rows(merged_rows: list[dict], cases: list[dict]) -> list[dict]:
    cbi = {c["id"]: c for c in cases}
    cell_means = aggregate_a2_per_cell(merged_rows)

    out: list[dict] = []
    for (model, tc, pop), cond_data in cell_means.items():
        ns = cond_data.get("no_spec")
        rag = cond_data.get("cpic_rag")
        if not ns or not rag:
            continue
        if rag["mean"] is None or ns["mean"] is None:
            continue
        if rag["mean"] >= ns["mean"]:
            continue

        case = cbi.get(tc, {})
        queried_drug = case.get("drug", "")
        gt_drug = case.get("gt_drug", "")
        cpic_rag_drug = first_parsed_drug(rag["rows"])

        category, gt_dir, rag_dir = classify(queried_drug, gt_drug, cpic_rag_drug)
        out.append({
            "category": category,
            "model": model,
            "tc": tc,
            "pop": pop,
            "queried_drug": queried_drug,
            "gt_drug": gt_drug,
            "cpic_rag_DRUG": cpic_rag_drug,
            "gt_direction": gt_dir,
            "cpic_rag_direction": rag_dir,
            "A2_no_spec_mean": f"{ns['mean']:.4f}",
            "A2_cpic_rag_mean": f"{rag['mean']:.4f}",
        })
    return out


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def write_summary(rows: list[dict], path: Path) -> None:
    cats = Counter(r["category"] for r in rows)
    total = len(rows)
    drug_sub = cats.get("a_drug_substitution", 0)
    phr = cats.get("b_phrasing_equiv", 0)
    wd = cats.get("c_wrong_direction", 0)
    other = cats.get("d_other", 0)

    by_gene_total: Counter = Counter()
    by_gene_subs: Counter = Counter()
    tc2gene: dict[str, str] = {}
    main = json.loads(MERGED.read_text())
    for m in main:
        tc2gene[m["tc"]] = m["gene"]
    for r in rows:
        g = tc2gene.get(r["tc"], "UNKNOWN")
        by_gene_total[g] += 1
        if r["category"] == "a_drug_substitution":
            by_gene_subs[g] += 1

    lines = [
        "# A2 regression classification (locked 3-run dataset)",
        f"Source: {MERGED.name}",
        f"Total A2 regression combos (cpic_rag mean < no_spec mean): {total}",
        "",
        "== Aggregate by category ==",
        f"  a_drug_substitution: {drug_sub:4d}  ({drug_sub/total*100:.1f}%)",
        f"  b_phrasing_equiv:    {phr:4d}  ({phr/total*100:.1f}%)",
        f"  c_wrong_direction:   {wd:4d}  ({wd/total*100:.1f}%)",
        f"  d_other:             {other:4d}  ({other/total*100:.1f}%)",
        "",
        "== Combined chunk/multi-drug structural confusion ==",
        f"  drug_sub + wrong_dir + other = {drug_sub + wd + other:4d}  "
        f"({(drug_sub + wd + other)/total*100:.1f}%)",
        "",
        "== Drug-substitution rate per gene ==",
        "  gene                | regressions | drug-sub | rate",
    ]
    for gene in sorted(by_gene_total, key=lambda g: -by_gene_subs[g]):
        tot = by_gene_total[gene]
        sub = by_gene_subs[gene]
        rate = sub / tot * 100 if tot else 0.0
        lines.append(f"  {gene:20s}| {tot:11d} | {sub:8d} | {rate:5.1f}%")
    path.write_text("\n".join(lines) + "\n")


def self_test() -> int:
    failures: list[str] = []
    if not MERGED.exists():
        failures.append(f"missing input: {MERGED}")
    if not CASES_FILE.exists():
        failures.append(f"missing input: {CASES_FILE}")

    if not failures:
        rows = json.loads(MERGED.read_text())
        cases = json.loads(CASES_FILE.read_text())
        if len(rows) != 26730:
            failures.append(f"merged dataset row count: got {len(rows)}, want 26730")
        if len(cases) != 110:
            failures.append(f"test cases row count: got {len(cases)}, want 110")

        test_drug = "tamoxifen"
        test_response = "codeine: Use codeine label recommended age- or weight-specific dosing"
        cat, gt_dir, rag_dir = classify(
            test_drug,
            "tamoxifen: consider higher dose or alternative; reduced endoxifen formation",
            test_response,
        )
        if cat != "a_drug_substitution":
            failures.append(f"classify tamox/codeine: got {cat}, want a_drug_substitution")

        cat2, _, _ = classify(
            "codeine",
            "codeine: AVOID (no analgesic benefit)",
            "codeine: standard dosing",
        )
        if cat2 != "c_wrong_direction":
            failures.append(f"classify codeine wrong-dir: got {cat2}, want c_wrong_direction")

        cat3, _, _ = classify(
            "codeine",
            "codeine: AVOID (no analgesic benefit)",
            "codeine: contraindicated for this patient",
        )
        if cat3 != "b_phrasing_equiv":
            failures.append(f"classify codeine AVOID/contraindicated equiv: got {cat3}, want b_phrasing_equiv")

        cat4, _, _ = classify(
            "codeine",
            "codeine: AVOID",
            "codeine: response is unclear from available evidence",
        )
        if cat4 != "d_other":
            failures.append(f"classify codeine unclear-dir: got {cat4}, want d_other")

    if failures:
        print("SELF-TEST FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("SELF-TEST PASSED (inputs present; 4 classifier unit assertions OK).")
    return 0


def compare(existing_csv: Path) -> int:
    if not existing_csv.exists():
        print(f"no existing CSV at {existing_csv} to compare against")
        return 1
    existing = list(csv.DictReader(existing_csv.open()))
    rows = json.loads(MERGED.read_text())
    cases = json.loads(CASES_FILE.read_text())
    new_rows = build_classified_rows(rows, cases)

    ec = Counter(r["category"] for r in existing)
    nc = Counter(r["category"] for r in new_rows)
    et = len(existing)
    nt = len(new_rows)

    print(f"{'category':22s} | {'existing':>10s} | {'new (locked 3-run)':>20s}")
    print(f"{'total combos':22s} | {et:>10d} | {nt:>20d}")
    for cat in ("a_drug_substitution", "b_phrasing_equiv", "c_wrong_direction", "d_other"):
        eo = ec.get(cat, 0)
        no = nc.get(cat, 0)
        e_pct = eo / et * 100 if et else 0.0
        n_pct = no / nt * 100 if nt else 0.0
        print(f"{cat:22s} | {eo:>4d} ({e_pct:4.1f}%) | {no:>4d} ({n_pct:4.1f}%)")
    return 0


def archive_existing(target: Path) -> None:
    if not target.exists():
        return
    arch = target.with_name(target.stem + "_pre3run.csv")
    if arch.exists():
        return
    arch.write_bytes(target.read_bytes())
    print(f"archived existing CSV -> {arch.relative_to(BASE)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1] if __doc__ else None)
    parser.add_argument("--self-test", action="store_true",
                        help="Sanity check inputs and classifier unit assertions; do not write.")
    parser.add_argument("--compare", action="store_true",
                        help="Diff new classifier output against existing CSV; do not write.")
    args = parser.parse_args()

    if args.self_test:
        return self_test()

    rc = self_test()
    if rc != 0:
        print("Aborting: self-test failed.", file=sys.stderr)
        return rc

    if args.compare:
        return compare(OUT_CSV)

    rows = json.loads(MERGED.read_text())
    cases = json.loads(CASES_FILE.read_text())
    new_rows = build_classified_rows(rows, cases)

    archive_existing(OUT_CSV)
    write_csv(new_rows, OUT_CSV)
    write_summary(new_rows, OUT_SUMMARY)

    cats = Counter(r["category"] for r in new_rows)
    total = len(new_rows)
    drug_sub = cats.get("a_drug_substitution", 0)
    print(f"wrote {OUT_CSV.relative_to(BASE)} ({total} regression combos)")
    print(f"wrote {OUT_SUMMARY.relative_to(BASE)}")
    print(f"headline: drug-substitution {drug_sub}/{total} = {drug_sub/total*100:.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
