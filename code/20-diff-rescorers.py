#!/usr/bin/env python3
"""
Diff between rigorous 10- rescorer and clinical-equivalence 10b- rescorer.

Reads:
  ../RESULTS/v3_raw_rescored.json              (rigorous 10- baseline)
  ../RESULTS/v3_raw_rescored_clinical_eq.json  (10b- with HLA-B*57:01 equiv.)
  ../SPECS/test_cases_v3.json                  (for lethal flag)

Writes:
  ../RESULTS/v3_rescore_clinical_eq_diff.json  (machine-readable summary)
  ../RESULTS/v3_rescore_clinical_eq_diff.txt   (human-readable table)

Diff table sections (per agreed spec):
  1. Aggregate per-condition A1 mean before/after (no_spec / with_spec)
  2. Lethal-class A1 mean before/after (no_spec / with_spec)
  3. Lethal-class A3 mean before/after (sanity; should be unchanged)
  4. Per-gene-per-condition transition counts: 0->1, 1->0, unchanged.
     ABORT if any 1->0 nonzero.
  5. A1/A3 alignment on lethal-class cases (A3=1 & A1=0 misalignment count).

Note on "gen1 only": the v3 dataset is the gen1 two-arm benchmark — all 17,820
rows are gen1. No filtering needed.
"""
from __future__ import annotations
import json
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
BASELINE = BASE / "RESULTS" / "v3_raw_rescored.json"
CLIN_EQ = BASE / "RESULTS" / "v3_raw_rescored_clinical_eq.json"
CASES_FILE = BASE / "SPECS" / "test_cases_v3.json"
OUT_JSON = BASE / "RESULTS" / "v3_rescore_clinical_eq_diff.json"
OUT_TXT = BASE / "RESULTS" / "v3_rescore_clinical_eq_diff.txt"


def is_lethal(tc: dict) -> bool:
    return "lethal" in tc.get("gt_drug", "").lower()


def row_key(r: dict) -> tuple:
    return (r["run"], r["model"], r["tc"], r["pop"], r["cond"])


def main() -> int:
    baseline = json.loads(BASELINE.read_text())
    clin_eq = json.loads(CLIN_EQ.read_text())
    cases = json.loads(CASES_FILE.read_text())
    cbi = {c["id"]: c for c in cases}

    if len(baseline) != len(clin_eq):
        print(f"FATAL: row counts differ ({len(baseline)} vs {len(clin_eq)})")
        return 2

    # Index by row key — files are written in the same order, but be safe.
    base_by_key = {row_key(r): r for r in baseline}
    eq_by_key = {row_key(r): r for r in clin_eq}
    if base_by_key.keys() != eq_by_key.keys():
        print("FATAL: row keys differ between baseline and clin_eq")
        return 2

    lines = []
    lines.append("# Clinical-equivalence rescorer diff")
    lines.append(f"Baseline:        {BASELINE.name}")
    lines.append(f"Clinical-eq:     {CLIN_EQ.name}")
    lines.append(f"Total rows:      {len(baseline)}")
    lines.append("")

    # Build per-row joined view, dropping format_fail rows (excluded from
    # parsed-only stats per the original rescorer convention).
    joined = []
    for k in base_by_key:
        b = base_by_key[k]
        e = eq_by_key[k]
        if b["scores"].get("format_fail") or e["scores"].get("format_fail"):
            continue
        tc = cbi[b["tc"]]
        joined.append({
            "model": b["model"], "tc": b["tc"], "pop": b["pop"],
            "cond": b["cond"], "run": b["run"], "gene": tc["gene"],
            "lethal": is_lethal(tc),
            "old_A1": float(b["scores"].get("A1", 0)),
            "new_A1": float(e["scores"].get("A1", 0)),
            "old_A3": float(b["scores"].get("A3", 1)),
            "new_A3": float(e["scores"].get("A3", 1)),
        })
    lines.append(f"Parsed-only rows (excluding format_fail): {len(joined)}")
    lines.append("")

    # ---- Section 1: Aggregate per-condition A1 mean ---------------------------
    lines.append("=" * 70)
    lines.append("Section 1: Aggregate per-condition A1 mean (parsed only)")
    lines.append("=" * 70)
    lines.append(f"{'condition':<12} {'n':>6}  {'old_A1_mean':>12} {'new_A1_mean':>12} {'delta':>10}")
    aggregate = {}
    for cond in ["no_spec", "with_spec"]:
        rs = [r for r in joined if r["cond"] == cond]
        n = len(rs)
        old_mean = sum(r["old_A1"] for r in rs) / n
        new_mean = sum(r["new_A1"] for r in rs) / n
        delta = new_mean - old_mean
        lines.append(f"{cond:<12} {n:>6}  {old_mean:>12.6f} {new_mean:>12.6f} {delta:>+10.6f}")
        aggregate[cond] = {"n": n, "old_mean": old_mean, "new_mean": new_mean, "delta": delta}

    # Gap (no_spec vs with_spec) before/after — the manuscript's headline
    gap_before = aggregate["with_spec"]["old_mean"] - aggregate["no_spec"]["old_mean"]
    gap_after = aggregate["with_spec"]["new_mean"] - aggregate["no_spec"]["new_mean"]
    gap_change = gap_after - gap_before
    lines.append("")
    lines.append(f"  no_spec -> with_spec gap:  before={gap_before:+.6f}  after={gap_after:+.6f}  compression={-gap_change:+.6f}")

    # ---- Section 2: Lethal-class A1 mean -------------------------------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 2: Lethal-class A1 mean (parsed only)")
    lines.append("=" * 70)
    lines.append(f"{'condition':<12} {'n':>6}  {'old_A1_mean':>12} {'new_A1_mean':>12} {'delta':>10}")
    lethal_a1 = {}
    for cond in ["no_spec", "with_spec"]:
        rs = [r for r in joined if r["cond"] == cond and r["lethal"]]
        n = len(rs)
        old_mean = sum(r["old_A1"] for r in rs) / n if n else 0.0
        new_mean = sum(r["new_A1"] for r in rs) / n if n else 0.0
        delta = new_mean - old_mean
        lines.append(f"{cond:<12} {n:>6}  {old_mean:>12.6f} {new_mean:>12.6f} {delta:>+10.6f}")
        lethal_a1[cond] = {"n": n, "old_mean": old_mean, "new_mean": new_mean, "delta": delta}

    if lethal_a1["no_spec"]["n"] and lethal_a1["with_spec"]["n"]:
        lethal_gap_before = lethal_a1["with_spec"]["old_mean"] - lethal_a1["no_spec"]["old_mean"]
        lethal_gap_after = lethal_a1["with_spec"]["new_mean"] - lethal_a1["no_spec"]["new_mean"]
        lethal_gap_change = lethal_gap_after - lethal_gap_before
        lines.append("")
        lines.append(f"  Lethal-class no_spec -> with_spec gap:  before={lethal_gap_before:+.6f}  after={lethal_gap_after:+.6f}  compression={-lethal_gap_change:+.6f}")

    # ---- Section 3: Lethal-class A3 mean (sanity) ----------------------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 3: Lethal-class A3 mean (sanity — should be unchanged)")
    lines.append("=" * 70)
    lines.append(f"{'condition':<12} {'n':>6}  {'old_A3_mean':>12} {'new_A3_mean':>12} {'delta':>10}")
    lethal_a3 = {}
    for cond in ["no_spec", "with_spec"]:
        rs = [r for r in joined if r["cond"] == cond and r["lethal"]]
        n = len(rs)
        old_mean = sum(r["old_A3"] for r in rs) / n if n else 0.0
        new_mean = sum(r["new_A3"] for r in rs) / n if n else 0.0
        delta = new_mean - old_mean
        flag = "" if abs(delta) < 1e-9 else "  <-- UNEXPECTED"
        lines.append(f"{cond:<12} {n:>6}  {old_mean:>12.6f} {new_mean:>12.6f} {delta:>+10.6f}{flag}")
        lethal_a3[cond] = {"n": n, "old_mean": old_mean, "new_mean": new_mean, "delta": delta}

    # ---- Section 4: Per-gene-per-condition transition counts -----------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 4: Per-gene-per-condition A1 transition counts")
    lines.append("  (rows where old_A1 in {0.0, 0.5, 1.0}, new_A1 in {0.0, 0.5, 1.0})")
    lines.append("=" * 70)
    transitions = defaultdict(lambda: defaultdict(int))
    # bucket: 0->1, 0->0.5, 0.5->1, 1->0.5, 1->0, 0.5->0, 0->0, 0.5->0.5, 1->1
    for r in joined:
        ob = round(r["old_A1"] * 2) / 2
        nb = round(r["new_A1"] * 2) / 2
        bucket = f"{ob:g}->{nb:g}"
        transitions[(r["gene"], r["cond"])][bucket] += 1

    # Print as table — focus on demotions
    abort = False
    lines.append(f"{'gene':<14} {'cond':<10}  {'unchanged':>10} {'0->1':>6} {'0->0.5':>8} {'0.5->1':>8} {'1->0':>6} {'1->0.5':>8} {'0.5->0':>8}")
    for gene, cond in sorted(transitions):
        t = transitions[(gene, cond)]
        unchanged = t.get("0->0", 0) + t.get("0.5->0.5", 0) + t.get("1->1", 0)
        n_01 = t.get("0->1", 0)
        n_005 = t.get("0->0.5", 0)
        n_051 = t.get("0.5->1", 0)
        n_10 = t.get("1->0", 0)
        n_105 = t.get("1->0.5", 0)
        n_050 = t.get("0.5->0", 0)
        flag = "  <-- ABORT" if (n_10 + n_105 + n_050) > 0 else ""
        if (n_10 + n_105 + n_050) > 0:
            abort = True
        lines.append(f"{gene:<14} {cond:<10}  {unchanged:>10} {n_01:>6} {n_005:>8} {n_051:>8} {n_10:>6} {n_105:>8} {n_050:>8}{flag}")

    # ---- Section 5: Lethal-class A1/A3 alignment under new rescorer ----------
    lines.append("")
    lines.append("=" * 70)
    lines.append("Section 5: Lethal-class A1/A3 alignment")
    lines.append("  Misalignment = A3=1 (correct safety call) but A1<1 (wrong phenotype label)")
    lines.append("=" * 70)
    lines.append(f"{'condition':<12} {'n':>6}  {'mis_A3=1_A1<1_OLD':>20} {'mis_A3=1_A1<1_NEW':>20} {'delta':>8}")
    alignment = {}
    for cond in ["no_spec", "with_spec"]:
        rs = [r for r in joined if r["cond"] == cond and r["lethal"]]
        n = len(rs)
        old_mis = sum(1 for r in rs if r["old_A3"] >= 1.0 and r["old_A1"] < 1.0)
        new_mis = sum(1 for r in rs if r["new_A3"] >= 1.0 and r["new_A1"] < 1.0)
        d = new_mis - old_mis
        lines.append(f"{cond:<12} {n:>6}  {old_mis:>20} {new_mis:>20} {d:>+8}")
        alignment[cond] = {"n": n, "old_misalign": old_mis, "new_misalign": new_mis, "delta": d}

    # ---- Verdict ------------------------------------------------------------
    lines.append("")
    lines.append("=" * 70)
    lines.append("VERDICT")
    lines.append("=" * 70)
    if abort:
        lines.append("ABORT: 1->0 (or 1->0.5, 0.5->0) transitions detected. Investigate before proceeding.")
    else:
        lines.append("No demotion transitions. Integrity gate holds.")
        # Summarise the key numbers
        lines.append("")
        lines.append("Key numbers:")
        lines.append(f"  Aggregate no_spec A1:     {aggregate['no_spec']['old_mean']:.4f} -> {aggregate['no_spec']['new_mean']:.4f}  ({aggregate['no_spec']['delta']:+.4f})")
        lines.append(f"  Aggregate with_spec A1:   {aggregate['with_spec']['old_mean']:.4f} -> {aggregate['with_spec']['new_mean']:.4f}  ({aggregate['with_spec']['delta']:+.4f})")
        lines.append(f"  Aggregate gap:            {gap_before:+.4f} -> {gap_after:+.4f}  (compression {-gap_change:+.4f})")
        if lethal_a1["no_spec"]["n"]:
            lines.append(f"  Lethal-class no_spec A1:  {lethal_a1['no_spec']['old_mean']:.4f} -> {lethal_a1['no_spec']['new_mean']:.4f}  ({lethal_a1['no_spec']['delta']:+.4f})")
            lines.append(f"  Lethal-class with_spec A1:{lethal_a1['with_spec']['old_mean']:.4f} -> {lethal_a1['with_spec']['new_mean']:.4f}  ({lethal_a1['with_spec']['delta']:+.4f})")

    text = "\n".join(lines)
    OUT_TXT.write_text(text)

    # Also write a JSON summary for downstream tooling
    summary = {
        "n_total_rows": len(baseline),
        "n_parsed_rows": len(joined),
        "section_1_aggregate_a1": aggregate,
        "section_1_gap": {"before": gap_before, "after": gap_after, "compression": -gap_change},
        "section_2_lethal_a1": lethal_a1,
        "section_3_lethal_a3": lethal_a3,
        "section_4_transitions": {f"{g}|{c}": dict(t) for (g, c), t in transitions.items()},
        "section_4_abort": abort,
        "section_5_alignment": alignment,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    print(text)
    print()
    print(f"Wrote {OUT_TXT.name}, {OUT_JSON.name}")
    return 1 if abort else 0


if __name__ == "__main__":
    sys.exit(main())
