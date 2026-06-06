#!/usr/bin/env python3
"""
Three-arm analysis: produces the headline numbers and audit candidates from
the merged three-arm dataset.

Reads:
  ../RESULTS/v3_raw_rescored_three_arm.json  (output of 17-merge-rag-results.py)
  ../SPECS/test_cases_v3.json

Writes:
  ../RESULTS/v3_three_arm_analysis_report.txt
  ../RESULTS/v3_three_arm_analysis_summary.json
  ../RESULTS/v3_three_arm_per_case_a1.csv          (figure data)
  ../RESULTS/v3_three_arm_lethal_a3_errors.csv     (safety story)
  ../RESULTS/v3_three_arm_audit_candidates.csv     (per-gene cells where
                                                    cpic_rag < no_spec)

Sections produced (matching the manuscript questions):
  1. Per-condition aggregate A1 / A2 / A3 (parsed only)
  2. Lethal-class A1 / A3 per condition
  3. A1 < A3 misalignment ("correctness-by-coincidence") per condition
  4. Per-case A1 per condition (CSV for figure)
  5. Per-gene per-condition A1 means + cpic_rag vs no_spec gap (audit)
  6. Lethal A3 errors per gene per condition (safety table)
  7. Bucket placement verdict (cpic_rag aggregate strictly between)
  8. Per-model per-condition A1 (per-model cpic_rag breakdown)
"""
from __future__ import annotations
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
MERGED = BASE / "RESULTS" / "v3_raw_rescored_three_arm.json"
CASES_FILE = BASE / "SPECS" / "test_cases_v3.json"
OUT_TXT = BASE / "RESULTS" / "v3_three_arm_analysis_report.txt"
OUT_JSON = BASE / "RESULTS" / "v3_three_arm_analysis_summary.json"
OUT_PERCASE_CSV = BASE / "RESULTS" / "v3_three_arm_per_case_a1.csv"
OUT_LETHAL_CSV = BASE / "RESULTS" / "v3_three_arm_lethal_a3_errors.csv"
OUT_AUDIT_CSV = BASE / "RESULTS" / "v3_three_arm_audit_candidates.csv"

CONDS = ("no_spec", "cpic_rag", "with_spec")


def is_lethal(tc: dict) -> bool:
    return "lethal" in tc.get("gt_drug", "").lower()


def safe_mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def main() -> int:
    if not MERGED.exists():
        print(f"FATAL: {MERGED.name} not found. Run 17-merge-rag-results.py first.",
              file=sys.stderr)
        return 2

    rows = json.loads(MERGED.read_text())
    cases = json.loads(CASES_FILE.read_text())
    cbi = {c["id"]: c for c in cases}

    # Filter to parsed only (drop format_fail) once for reuse
    parsed = [r for r in rows if not r["scores"].get("format_fail")]
    conds_present = sorted({r["cond"] for r in rows},
                           key=lambda c: CONDS.index(c) if c in CONDS else 99)

    lines = [
        "# v3 THREE-ARM ANALYSIS REPORT",
        f"Source: {MERGED.name}",
        f"Total rows: {len(rows)}  (parsed-only, after format_fail filter: {len(parsed)})",
        f"Conditions present: {conds_present}",
        "",
    ]

    # ---- Section 1: Per-condition aggregate A1/A2/A3 -------------------------
    lines.append("=" * 70)
    lines.append("Section 1: Per-condition aggregate (parsed only)")
    lines.append("=" * 70)
    lines.append(f"{'condition':<12} {'n':>6}  {'A1 mean':>10} {'A2 mean':>10} {'A3 mean':>10}")
    aggregate = {}
    for cond in conds_present:
        rs = [r for r in parsed if r["cond"] == cond]
        n = len(rs)
        a1 = safe_mean([r["scores"].get("A1", 0) for r in rs])
        a2 = safe_mean([r["scores"].get("A2", 0) for r in rs])
        a3 = safe_mean([r["scores"].get("A3", 0) for r in rs])
        lines.append(f"{cond:<12} {n:>6}  {a1:>10.4f} {a2:>10.4f} {a3:>10.4f}")
        aggregate[cond] = {"n": n, "A1": a1, "A2": a2, "A3": a3}

    # ---- Section 2: Lethal-class A1/A3 per condition -------------------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 2: Lethal-class A1 / A3 per condition")
    lines.append("=" * 70)
    lines.append(f"{'condition':<12} {'n_lethal':>10}  {'A1 mean':>10} {'A3 mean':>10}  {'A3 errors':>10}")
    lethal = {}
    for cond in conds_present:
        rs = [r for r in parsed if r["cond"] == cond and is_lethal(cbi[r["tc"]])]
        n = len(rs)
        a1 = safe_mean([r["scores"].get("A1", 0) for r in rs])
        a3 = safe_mean([r["scores"].get("A3", 0) for r in rs])
        a3_errors = sum(1 for r in rs if r["scores"].get("A3", 1) < 1.0)
        lines.append(f"{cond:<12} {n:>10}  {a1:>10.4f} {a3:>10.4f}  {a3_errors:>10}")
        lethal[cond] = {"n": n, "A1": a1, "A3": a3, "A3_errors": a3_errors}

    # ---- Section 3: A1/A3 misalignment per condition ------------------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 3: A1 < 1 within A3 = 1 (correctness-by-coincidence)")
    lines.append("  Lethal-class subset only — the manuscript's safety story.")
    lines.append("=" * 70)
    lines.append(f"{'condition':<12} {'A3=1':>8} {'A1<1 in A3=1':>14} {'A1=0 in A3=1':>14} {'pct A1<1':>10}")
    misalign = {}
    for cond in conds_present:
        rs = [r for r in parsed if r["cond"] == cond and is_lethal(cbi[r["tc"]])]
        a3_eq_1 = [r for r in rs if r["scores"].get("A3", 1) >= 1.0]
        a1_lt_1 = [r for r in a3_eq_1 if r["scores"].get("A1", 0) < 1.0]
        a1_eq_0 = [r for r in a3_eq_1 if r["scores"].get("A1", 0) == 0.0]
        pct = (len(a1_lt_1) / len(a3_eq_1) * 100) if a3_eq_1 else 0.0
        lines.append(f"{cond:<12} {len(a3_eq_1):>8} {len(a1_lt_1):>14} {len(a1_eq_0):>14} {pct:>9.2f}%")
        misalign[cond] = {"A3_eq_1": len(a3_eq_1), "A1_lt_1": len(a1_lt_1),
                          "A1_eq_0": len(a1_eq_0), "pct_misaligned": pct}

    # ---- Section 4: Per-case A1 per condition (CSV for figure) --------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 4: Per-case A1 per condition")
    lines.append(f"  Full table written to: {OUT_PERCASE_CSV.name}")
    lines.append("=" * 70)
    case_ids = sorted({r["tc"] for r in parsed})
    with OUT_PERCASE_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case_id", "gene", "drug", "lethal"] + [f"A1_{c}" for c in conds_present])
        for tc_id in case_ids:
            tc = cbi[tc_id]
            row = [tc_id, tc["gene"], tc["drug"], int(is_lethal(tc))]
            for cond in conds_present:
                rs = [r for r in parsed if r["tc"] == tc_id and r["cond"] == cond]
                row.append(f"{safe_mean([r['scores'].get('A1', 0) for r in rs]):.4f}" if rs else "")
            w.writerow(row)

    # Print top 10 cases where cpic_rag has the largest gap below no_spec
    if "no_spec" in conds_present and "cpic_rag" in conds_present:
        gaps = []
        for tc_id in case_ids:
            tc = cbi[tc_id]
            ns = [r for r in parsed if r["tc"] == tc_id and r["cond"] == "no_spec"]
            cr = [r for r in parsed if r["tc"] == tc_id and r["cond"] == "cpic_rag"]
            if not ns or not cr:
                continue
            ns_a1 = safe_mean([r["scores"].get("A1", 0) for r in ns])
            cr_a1 = safe_mean([r["scores"].get("A1", 0) for r in cr])
            gaps.append((tc_id, tc["gene"], ns_a1, cr_a1, cr_a1 - ns_a1))
        gaps.sort(key=lambda x: x[4])  # ascending: most-negative first
        lines.append("")
        lines.append("  Cases where cpic_rag underperforms no_spec (top 10):")
        lines.append(f"  {'case':<28} {'gene':<14} {'no_spec':>9} {'cpic_rag':>9} {'delta':>9}")
        for tc_id, gene, ns_a1, cr_a1, d in gaps[:10]:
            if d >= 0:
                continue
            lines.append(f"  {tc_id:<28} {gene:<14} {ns_a1:>9.4f} {cr_a1:>9.4f} {d:>+9.4f}")

    # ---- Section 5: Per-gene per-condition A1 + cpic_rag vs no_spec gap -----
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 5: Per-gene per-condition A1 mean (audit candidates)")
    lines.append("  Genes where cpic_rag < no_spec are next-round equivalence audit")
    lines.append("  candidates or Limitations entries.")
    lines.append("=" * 70)
    genes = sorted({cbi[r["tc"]]["gene"] for r in parsed})
    audit_rows = []
    header = f"  {'gene':<14}"
    for cond in conds_present:
        header += f" {cond:>10}"
    header += f" {'rag-ns':>10}"
    lines.append(header)
    for gene in genes:
        gene_means = {}
        for cond in conds_present:
            rs = [r for r in parsed if cbi[r["tc"]]["gene"] == gene and r["cond"] == cond]
            gene_means[cond] = safe_mean([r["scores"].get("A1", 0) for r in rs])
        rag_ns = (gene_means.get("cpic_rag", 0) - gene_means.get("no_spec", 0)
                  if "cpic_rag" in gene_means and "no_spec" in gene_means else None)
        line = f"  {gene:<14}"
        for cond in conds_present:
            line += f" {gene_means[cond]:>10.4f}"
        line += f" {rag_ns:>+10.4f}" if rag_ns is not None else f" {'':>10}"
        if rag_ns is not None and rag_ns < 0:
            line += "  <-- audit"
            audit_rows.append((gene, gene_means, rag_ns))
        lines.append(line)

    # Audit candidate CSV
    with OUT_AUDIT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        cols = ["gene"] + [f"A1_{c}" for c in conds_present] + ["cpic_rag_minus_no_spec"]
        w.writerow(cols)
        for gene, gm, d in audit_rows:
            w.writerow([gene] + [f"{gm[c]:.4f}" for c in conds_present] + [f"{d:+.4f}"])

    # ---- Section 6: Lethal A3 errors per gene per condition -----------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 6: Lethal-class A3 errors per gene per condition")
    lines.append(f"  Full table also written to: {OUT_LETHAL_CSV.name}")
    lines.append("=" * 70)
    header = f"  {'gene':<14}"
    for cond in conds_present:
        header += f" {cond+' err/n':>16}"
    lines.append(header)
    lethal_csv_rows = []
    for gene in genes:
        line = f"  {gene:<14}"
        csv_row = [gene]
        any_lethal = False
        for cond in conds_present:
            rs = [r for r in parsed if cbi[r["tc"]]["gene"] == gene
                  and r["cond"] == cond and is_lethal(cbi[r["tc"]])]
            if rs:
                any_lethal = True
            errs = sum(1 for r in rs if r["scores"].get("A3", 1) < 1.0)
            line += f" {errs:>6}/{len(rs):<9}"
            csv_row.extend([errs, len(rs)])
        if any_lethal:
            lines.append(line)
            lethal_csv_rows.append(csv_row)
    with OUT_LETHAL_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        cols = ["gene"]
        for cond in conds_present:
            cols.extend([f"{cond}_a3_errors", f"{cond}_lethal_n"])
        w.writerow(cols)
        for r in lethal_csv_rows:
            w.writerow(r)

    # ---- Section 7: Bucket placement verdict --------------------------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 7: Bucket placement verdict")
    lines.append("=" * 70)
    bucket = {}
    if {"no_spec", "cpic_rag", "with_spec"}.issubset(conds_present):
        ns = aggregate["no_spec"]["A1"]
        cr = aggregate["cpic_rag"]["A1"]
        ws = aggregate["with_spec"]["A1"]
        gap_total = ws - ns
        gap_closed_by_rag = (cr - ns) / gap_total * 100 if gap_total > 0 else 0
        lines.append(f"  no_spec   A1: {ns:.4f}")
        lines.append(f"  cpic_rag  A1: {cr:.4f}")
        lines.append(f"  with_spec A1: {ws:.4f}")
        lines.append(f"  Total gap (with_spec - no_spec):  {gap_total:+.4f}")
        lines.append(f"  Gap closed by RAG (cpic_rag - no_spec): {cr - ns:+.4f}")
        lines.append(f"  Gap closure %: {gap_closed_by_rag:.1f}%")
        lines.append("")
        if cr <= ns + 0.001:
            lines.append("  VERDICT: cpic_rag <= no_spec — RAG provides no benefit. Investigate.")
            placement = "underperforms"
        elif cr >= ws - 0.001:
            lines.append("  VERDICT: cpic_rag matches with_spec — RAG sufficient modulo vocabulary.")
            placement = "matches_with_spec"
        else:
            lines.append("  VERDICT: cpic_rag strictly between no_spec and with_spec.")
            lines.append(f"           RAG closes {gap_closed_by_rag:.0f}% of the gap; SKILL.md "
                         f"closes the remaining {100 - gap_closed_by_rag:.0f}%.")
            placement = "between"
        bucket = {"no_spec": ns, "cpic_rag": cr, "with_spec": ws,
                  "gap_total": gap_total, "rag_closure_pct": gap_closed_by_rag,
                  "placement": placement}

    # ---- Section 8: Per-model per-condition A1 ------------------------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 8: Per-model per-condition A1 (parsed only)")
    lines.append("=" * 70)
    models = sorted({r["model"] for r in parsed})
    header = f"  {'model':<22}"
    for cond in conds_present:
        header += f" {cond:>10}"
    lines.append(header)
    per_model = {}
    for model in models:
        line = f"  {model:<22}"
        per_model[model] = {}
        for cond in conds_present:
            rs = [r for r in parsed if r["model"] == model and r["cond"] == cond]
            a1 = safe_mean([r["scores"].get("A1", 0) for r in rs])
            line += f" {a1:>10.4f}"
            per_model[model][cond] = a1
        lines.append(line)

    # ---- Write outputs ------------------------------------------------------
    text = "\n".join(lines)
    OUT_TXT.write_text(text)
    summary = {
        "n_total_rows": len(rows),
        "n_parsed_rows": len(parsed),
        "conditions_present": conds_present,
        "section_1_aggregate": aggregate,
        "section_2_lethal": lethal,
        "section_3_misalignment": misalign,
        "section_5_audit_candidate_genes": [g[0] for g in audit_rows],
        "section_7_bucket_placement": bucket,
        "section_8_per_model_a1": per_model,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    print(text)
    print()
    print(f"Wrote {OUT_TXT.name}, {OUT_JSON.name}, {OUT_PERCASE_CSV.name}, "
          f"{OUT_LETHAL_CSV.name}, {OUT_AUDIT_CSV.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
