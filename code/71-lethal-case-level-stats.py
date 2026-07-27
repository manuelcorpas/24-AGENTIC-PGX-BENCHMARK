#!/usr/bin/env python3
"""
Lethal-class comparison at the case level, with intervals (reviewer point 4).

WHY
The 336 lethal-class cells are 14 distinct cases replicated across eight models
and three runs, and the HLA stratum is 120 cells, that is five distinct cases.
A referee asked twice that the unit of analysis be the distinct case, not the
cell, and noted that the manuscript reported a doubling computed on
non-independent cells with no interval anywhere, while giving clustered
intervals for phenotype accuracy in Table 3.

This script does what was asked: resample the distinct cases, report the paired
difference in lethal-class error rate with a case-clustered interval, and do the
same for the HLA stratum separately. If the effect does not survive, it does not
belong in the Summary.

WHAT IT REPORTS AND WHY THAT SHAPE
The quantity is the difference in lethal-class error rate between the
retrieval-generated and free-generation cells. Cases are the resampling unit, so
the interval reflects the fact that fourteen clinical scenarios, not 336
independent observations, generated the evidence.

USAGE
    python code/71-lethal-case-level-stats.py
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
ROWS = BASE / "data" / "v3_five_cell_live_rows.json"
OUT = BASE / "data" / "v3_lethal_case_level.json"
REPORT = BASE / "data" / "v3_lethal_case_level.txt"

SEED = 20260727
N_BOOT = 10000


def error_rate(by_case, cases, cell):
    vals = [x for c in cases for x in by_case[c][cell]]
    return (sum(vals) / len(vals)) if vals else None


def bootstrap(by_case, cases, a, b, n_boot=N_BOOT, seed=SEED):
    rng = random.Random(seed)
    obs = error_rate(by_case, cases, b) - error_rate(by_case, cases, a)
    draws = []
    for _ in range(n_boot):
        s = [rng.choice(cases) for _ in cases]
        ra, rb = error_rate(by_case, s, a), error_rate(by_case, s, b)
        if ra is not None and rb is not None:
            draws.append(rb - ra)
    draws.sort()
    lo = draws[int(0.025 * len(draws))]
    hi = draws[int(0.975 * len(draws))]
    return {"difference": round(obs, 4),
            "ci_low": round(lo, 4), "ci_high": round(hi, 4),
            "excludes_zero": not (lo < 0 < hi),
            "n_cases": len(cases), "n_boot": len(draws)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--rows", type=Path, default=ROWS)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--report", type=Path, default=REPORT)
    a = ap.parse_args(argv)

    if not a.rows.exists():
        sys.stderr.write(f"missing {a.rows}; run 61-rescore-matched.py --rows\n")
        return 1
    rows = json.loads(a.rows.read_text())
    lethal = [r for r in rows if r.get("lethal_action") is not None]
    if not lethal:
        sys.stderr.write("no lethal-class rows found\n")
        return 1

    by_case = defaultdict(lambda: defaultdict(list))
    genes = {}
    for r in lethal:
        by_case[r["case_id"]][r["cell"]].append(0.0 if r["lethal_action"] == 1.0 else 1.0)
        genes[r["case_id"]] = r["gene"]

    cases = sorted(by_case)
    hla = sorted(c for c in cases if genes[c].startswith("HLA"))
    other = sorted(c for c in cases if not genes[c].startswith("HLA"))

    strata = {"all": cases, "HLA": hla, "non-HLA": other}
    result = {"n_distinct_cases": len(cases), "n_hla_cases": len(hla),
              "n_lethal_cells_per_cell": len(lethal) // len(by_case[cases[0]]),
              "strata": {}}
    for name, cs in strata.items():
        if not cs:
            continue
        result["strata"][name] = {
            "cases": cs,
            "free_generation": round(error_rate(by_case, cs, "free_generation"), 4),
            "rag_generation": round(error_rate(by_case, cs, "rag_generation"), 4),
            "retrieval_minus_free": bootstrap(by_case, cs, "free_generation",
                                              "rag_generation"),
        }
    a.out.write_text(json.dumps(result, indent=2))

    L = ["LETHAL-CLASS COMPARISON AT THE CASE LEVEL", ""]
    L.append(f"  distinct lethal-class cases: {result['n_distinct_cases']} "
             f"(HLA {result['n_hla_cases']}, non-HLA "
             f"{result['n_distinct_cases'] - result['n_hla_cases']})")
    L.append("  The unit of analysis is the distinct case, not the cell: the 336")
    L.append("  lethal-class cells are these cases replicated across models and runs.")
    L.append("")
    for name, s in result["strata"].items():
        d = s["retrieval_minus_free"]
        L.append(f"  {name} (n = {d['n_cases']} cases)")
        L.append(f"    free generation        {s['free_generation']}")
        L.append(f"    retrieval, generated   {s['rag_generation']}")
        L.append(f"    difference             {d['difference']:+.4f}   "
                 f"case-clustered 95% CI ({d['ci_low']:+.4f}, {d['ci_high']:+.4f})")
        L.append(f"    excludes zero          {d['excludes_zero']}")
        L.append("")
    text = "\n".join(L)
    a.report.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
