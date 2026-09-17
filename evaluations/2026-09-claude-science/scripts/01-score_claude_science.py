#!/usr/bin/env python3
"""Score Claude Science's answers with the locked v2.10 matched scorer.

Reads RESULTS/claude_science_answers.csv (case_id, diplotype, phenotype,
recommendation, abstained), scores every one of the 110 cases with
code/61-rescore-matched.py score_row (missing cases count as failures), and
prints the result next to the manuscript's Table 1.

Usage: python PYTHON/01-score_claude_science.py [--repo PATH] [--answers CSV]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

# Manuscript v90 Results text and Table 1, baseline scorer, n = 2,640 per arm,
# lethal-class errors on 336 attempted evaluations. Numbers not stated in the
# text are left as None rather than guessed.
MANUSCRIPT_TABLE1 = {
    "free generation (no knowledge)": {"a2": 0.630, "lethal_errors": 46},
    "RAG generation": {"a2": 0.555, "lethal_errors": 137},
    "authored-rule generation (skill)": {"a2": 0.971, "lethal_errors": 24},
    "authored-rule execution (skill as code)": {"a2": 0.973, "lethal_errors": 23},
    "RAG-assisted execution": {"a2": 0.970, "lethal_errors": None},
}


def unblind(path: Path, surrogate_to_real: dict[str, str]) -> list[dict]:
    """Restore real case ids in an answers file produced from a hashed package.
    An unmapped id is left as it is, so a stray row is visible rather than lost."""
    with Path(path).open(newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        cid = r["case_id"].strip()
        r["case_id"] = surrogate_to_real.get(cid, cid)
    return rows


def read_answers(path: Path, cases: list[dict]) -> list[dict]:
    by_id = {}
    with Path(path).open(newline="") as f:
        for r in csv.DictReader(f):
            by_id[r["case_id"].strip()] = r
    rows = []
    for c in cases:
        r = by_id.get(c["id"])
        if r is None:
            rows.append({"cell": "claude_science", "case_id": c["id"], "raw": "",
                         "parsed_phenotype": "", "parsed_drug": "", "abstained": False})
            continue
        rows.append({
            "cell": "claude_science", "case_id": c["id"], "model": "claude-science",
            "raw": json.dumps(r), "parsed_diplotype": (r.get("diplotype") or "").strip(),
            "parsed_phenotype": (r.get("phenotype") or "").strip(),
            "parsed_drug": (r.get("recommendation") or "").strip(),
            "abstained": str(r.get("abstained", "")).strip().lower() in {"true", "1", "yes"},
        })
    return rows


def _load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=Path.home() / "dev" / "24-AGENTIC-PGX-BENCHMARK")
    ap.add_argument("--answers", type=Path, default=here / "RESULTS" / "claude_science_answers.csv")
    args = ap.parse_args()

    matched = _load(args.repo / "code" / "61-rescore-matched.py", "rescore_matched")
    cases = matched.CASES
    scored = [matched.score_row(r, matched.CASE_BY_ID[r["case_id"]])
              for r in read_answers(args.answers, cases)]

    n = len(scored)
    lethal = [s for s in scored if s["lethal_action"] is not None]
    summary = {
        "n_cases": n,
        "parsed_ok": sum(s["parsed_ok"] for s in scored),
        "abstained": sum(s["abstained"] for s in scored),
        "a1_phenotype": round(sum(s["a1_phenotype"] for s in scored) / n, 4),
        "a2_recommendation": round(sum(s["a2_recommendation"] for s in scored) / n, 4),
        "lethal_n": len(lethal),
        "lethal_errors": sum(1 for s in lethal if s["lethal_action"] == 0.0),
        "note": "one pass over 110 cases; manuscript arms are 8 models x 3 reps, so compare rates not counts",
    }
    out = here / "RESULTS"
    out.mkdir(exist_ok=True)
    (out / "claude_science_scored.json").write_text(json.dumps(scored, indent=1))
    (out / "claude_science_summary.json").write_text(
        json.dumps({"claude_science": summary, "manuscript_table1": MANUSCRIPT_TABLE1}, indent=1))
    print(json.dumps(summary, indent=1))
    for arm, v in MANUSCRIPT_TABLE1.items():
        print(f"{arm:42s} a2={v['a2']:.3f} lethal_errors/336={v['lethal_errors']}")


if __name__ == "__main__":
    main()
