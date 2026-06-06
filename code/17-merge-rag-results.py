#!/usr/bin/env python3
"""
Three-arm merge: combine the existing two-arm dataset (no_spec + with_spec,
already scored under 10b- in v3_raw_rescored_clinical_eq.json) with the
cpic_rag full-run output (v3_rag_raw.json), rescoring cpic_rag rows under
10b- for uniform treatment, and producing a single three-arm dataset.

Per Option B from HANDOFF-three-arm-rag.md: the clinical-equivalence rescorer
(10b-) is applied uniformly across all conditions. The original rigorous
rescorer (10-) is preserved as `scores_rigorous` on every row for transparency
and to allow downstream analyses (e.g. the substring false-positive screen in
19-validate-three-arm.py R3) to filter out 10b- promotions.

Reads:
  ../RESULTS/v3_raw_rescored_clinical_eq.json  (two-arm, already 10b- scored,
                                                 with scores_rigorous preserved)
  ../RESULTS/v3_rag_raw.json                   (cpic_rag full-run output,
                                                 preliminary scoring at runtime)
  ../SPECS/test_cases_v3.json

Writes:
  ../RESULTS/v3_raw_rescored_three_arm.json    (input to 18- and 19-)
  ../RESULTS/v3_three_arm_merge_report.txt
"""
from __future__ import annotations
import json
import sys
from collections import Counter, defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
PYDIR = BASE / "PYTHON"
TWO_ARM = BASE / "RESULTS" / "v3_raw_rescored_clinical_eq.json"
CPIC_RAG_RAW = BASE / "RESULTS" / "v3_rag_raw.json"
CASES = BASE / "SPECS" / "test_cases_v3.json"
OUT = BASE / "RESULTS" / "v3_raw_rescored_three_arm.json"
REPORT = BASE / "RESULTS" / "v3_three_arm_merge_report.txt"

# Import original rigorous rescorer (for scores_rigorous baseline) and
# clinical-equivalence rescorer (for headline scores).
_sp1 = spec_from_file_location("rescore", PYDIR / "10-rescore-v3.py")
rescore = module_from_spec(_sp1); _sp1.loader.exec_module(rescore)
_sp2 = spec_from_file_location("rescore_eq", PYDIR / "10b-rescore-v3-clinical-equivalence.py")
rescore_eq = module_from_spec(_sp2); _sp2.loader.exec_module(rescore_eq)


def rescore_pair(r: dict, tc: dict) -> tuple[dict, dict]:
    """Return (rigorous_scores, clinical_eq_scores) for a row's parsed output.
    Both share A2/A3/B1/B2/B3 — only A1 differs (10b- may promote).
    Preserves format_fail untouched."""
    if r["scores"].get("format_fail"):
        ff = dict(r["scores"])
        return ff, dict(ff)
    parsed = r.get("parsed", {}) or {}
    rigorous = rescore.rescore_row(r, tc)
    clin_eq = rescore_eq.rescore_row_clinical_eq(r, tc)
    return rigorous, clin_eq


def main() -> int:
    if not TWO_ARM.exists():
        print(f"FATAL: {TWO_ARM.name} not found. Run 10b- first.", file=sys.stderr)
        return 2
    if not CPIC_RAG_RAW.exists():
        print(f"FATAL: {CPIC_RAG_RAW.name} not found. Run the full cpic_rag arm first.",
              file=sys.stderr)
        return 2

    cases_by_id = {c["id"]: c for c in json.loads(CASES.read_text())}
    two_arm = json.loads(TWO_ARM.read_text())
    rag = json.loads(CPIC_RAG_RAW.read_text())

    out_rows = []

    # ---- Two-arm rows: already 10b- scored, scores_rigorous preserved ----
    # The clinical_eq output stores scores=10b-, scores_rigorous=rigorous-baseline.
    # We pass through as-is (do NOT rescore — these are locked headline numbers).
    cond_counts = Counter()
    for r in two_arm:
        cond_counts[r["cond"]] += 1
        out_rows.append(r)

    # ---- cpic_rag rows: rescore from parsed under both classifiers ----
    n_promoted = 0
    n_demoted = 0
    promotion_examples = []
    rag_format_fail = 0
    for r in rag:
        tc = cases_by_id[r["tc"]]
        rigorous, clin_eq = rescore_pair(r, tc)
        if r["scores"].get("format_fail"):
            rag_format_fail += 1
        # New row mirrors the two-arm shape:
        #   scores             = 10b- (headline)
        #   scores_rigorous    = rigorous 10- baseline
        #   scores_runtime     = preliminary score at API time (kept for audit)
        new_row = dict(r)
        new_row["scores_runtime"] = r["scores"]
        new_row["scores_rigorous"] = rigorous
        new_row["scores"] = clin_eq
        cond_counts[new_row["cond"]] += 1

        old_a1 = rigorous.get("A1", 0)
        new_a1 = clin_eq.get("A1", 0)
        if not new_row["scores"].get("format_fail"):
            if new_a1 > old_a1 + 0.01:
                n_promoted += 1
                if len(promotion_examples) < 30:
                    promotion_examples.append({
                        "model": r["model"], "tc": r["tc"], "pop": r["pop"],
                        "run": r["run"], "gene": tc["gene"],
                        "gt_phenotype": tc["gt_phenotype"],
                        "parsed_PHENOTYPE": (r["parsed"] or {}).get("PHENOTYPE", "")[:200],
                        "old_A1": old_a1, "new_A1": new_a1,
                    })
            if new_a1 < old_a1 - 0.01:
                n_demoted += 1
        out_rows.append(new_row)

    OUT.write_text(json.dumps(out_rows, indent=2))

    # ---- Report ----
    lines = [
        "# v3 three-arm merge report",
        f"Source two-arm:  {TWO_ARM.name}",
        f"Source cpic_rag: {CPIC_RAG_RAW.name}",
        f"Output:          {OUT.name}",
        "",
        f"Total rows in merged file:  {len(out_rows)}",
        f"  no_spec:    {cond_counts.get('no_spec', 0)}",
        f"  cpic_rag:   {cond_counts.get('cpic_rag', 0)}",
        f"  with_spec:  {cond_counts.get('with_spec', 0)}",
        "",
        f"cpic_rag format_fail count: {rag_format_fail} / {len(rag)}",
        "",
        f"cpic_rag rescoring (rigorous 10- vs 10b- clinical-equivalence):",
        f"  Promoted (A1 baseline<1 -> 1):  {n_promoted}",
        f"  Demoted  (A1 baseline=1 -> <1): {n_demoted}  (must be 0 by construction)",
        "",
    ]
    if promotion_examples:
        lines.append("## Sample cpic_rag promotions (up to 30)")
        for ex in promotion_examples:
            lines.append(
                f"  {ex['model']:<22} {ex['tc']:<28} {ex['pop']:<3} run={ex['run']} gene={ex['gene']}"
            )
            lines.append(f"    gt_phenotype: {ex['gt_phenotype']}")
            lines.append(f"    parsed:       {ex['parsed_PHENOTYPE']}")
            lines.append(f"    A1 {ex['old_A1']} -> {ex['new_A1']}")
            lines.append("")
    REPORT.write_text("\n".join(lines))
    print("\n".join(lines))
    print()
    print(f"Wrote {OUT.name}, {REPORT.name}")
    if n_demoted > 0:
        print(f"WARNING: {n_demoted} cpic_rag rows demoted; should be impossible.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
