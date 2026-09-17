#!/usr/bin/env python3
"""Score every replicate with the locked v2.10 scorer and report the spread.

The spread is the point. Claude Science exposes no temperature or seed, so
replicate variance is a measurement of the non-determinism that a versioned
skill does not have. A cell with no answers reports "no answers found" rather
than zero, because an empty result is only news when nothing failed.

Usage: python PYTHON/03-aggregate_runs.py [--runs RUNS] [--repo PATH]
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ANSWERS = "claude_science_answers.csv"


def _load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _scorer(repo: Path):
    return _load(repo / "code" / "61-rescore-matched.py", "rescore_matched")


def score_replicate(answers: Path, cases: list[dict], id_map: dict, matched) -> dict:
    """One replicate. Missing cases count as failures, as in the paper."""
    cs = _load(HERE / "PYTHON" / "01-score_claude_science.py", "score_cs_inner")
    rows = cs.unblind(answers, id_map)
    by_id = {r["case_id"]: r for r in rows}
    scored = []
    for c in cases:
        r = by_id.get(c["id"])
        row = {"cell": "run", "case_id": c["id"], "raw": "", "parsed_phenotype": "",
               "parsed_drug": "", "abstained": False}
        if r is not None:
            row.update({
                "raw": json.dumps(r),
                "parsed_phenotype": (r.get("phenotype") or "").strip(),
                "parsed_drug": (r.get("recommendation") or "").strip(),
                "abstained": str(r.get("abstained", "")).strip().lower() in {"true", "1", "yes"},
            })
        scored.append(matched.score_row(row, c))
    n = len(scored)
    lethal = [s for s in scored if s["lethal_action"] is not None]
    return {
        "n_cases": n,
        "phenotype": sum(s["a1_phenotype"] for s in scored) / n,
        "recommendation": sum(s["a2_recommendation"] for s in scored) / n,
        "lethal_errors": sum(1 for s in lethal if s["lethal_action"] == 0.0),
        "lethal_n": len(lethal),
        "abstained": sum(s["abstained"] for s in scored),
        "parsed_ok": sum(s["parsed_ok"] for s in scored),
    }


def _spread(values: list[float]) -> dict:
    return {
        "mean": round(statistics.fmean(values), 4),
        "sd": round(statistics.pstdev(values), 4),
        "min": round(min(values), 4),
        "max": round(max(values), 4),
        "per_replicate": [round(v, 4) for v in values],
    }


def aggregate(runs_root: Path, cases: list[dict], id_map: dict,
              repo: Path | None = None) -> dict:
    runs_root = Path(runs_root)
    matched = _scorer(repo or (Path.home() / "dev" / "24-AGENTIC-PGX-BENCHMARK"))
    out = {}
    for cell in sorted(p.name for p in runs_root.iterdir()
                       if p.is_dir() and not p.name.startswith("_")):
        reps = sorted((runs_root / cell).glob("rep*/" + ANSWERS))
        if not reps:
            out[cell] = {"n_replicates": 0, "status": "no answers found"}
            continue
        per = [score_replicate(p, cases, id_map, matched) for p in reps]
        out[cell] = {
            "n_replicates": len(per),
            "status": "ok",
            "replicates": [p.parent.name for p in reps],
            "phenotype": _spread([x["phenotype"] for x in per]),
            "recommendation": _spread([x["recommendation"] for x in per]),
            "lethal_errors": [x["lethal_errors"] for x in per],
            "lethal_n": per[0]["lethal_n"],
            "abstained": [x["abstained"] for x in per],
            "parsed_ok": [x["parsed_ok"] for x in per],
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, default=HERE / "RUNS")
    ap.add_argument("--repo", type=Path, default=Path.home() / "dev" / "24-AGENTIC-PGX-BENCHMARK")
    args = ap.parse_args()

    cases = json.loads((HERE.parent / "SPECS" / "test_cases_v3.json").read_text())
    id_map = json.loads((args.runs / "_KEYS" / "id_map.json").read_text())
    res = aggregate(args.runs, cases, id_map, repo=args.repo)
    (args.runs / "AGGREGATE.json").write_text(json.dumps(res, indent=1))

    print(f"{'cell':14s} {'n':>2s}  phenotype              recommendation         lethal errors")
    for cell, v in res.items():
        if v["n_replicates"] == 0:
            print(f"{cell:14s}  -  {v['status']}")
            continue
        p, r = v["phenotype"], v["recommendation"]
        print(f"{cell:14s} {v['n_replicates']:>2d}  "
              f"{p['mean']:.3f} sd {p['sd']:.3f} [{p['min']:.3f}-{p['max']:.3f}]  "
              f"{r['mean']:.3f} sd {r['sd']:.3f} [{r['min']:.3f}-{r['max']:.3f}]  "
              f"{v['lethal_errors']} of {v['lethal_n']}")
    print("\nPaper, baseline scorer: skill generation 0.966/0.971, skill execution 0.967/0.973,"
          "\nfree generation 0.743/0.630, RAG generation 0.822/0.555.")


if __name__ == "__main__":
    main()
