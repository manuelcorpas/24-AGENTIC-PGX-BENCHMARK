#!/usr/bin/env python3
"""Recompute every number in the September 2026 Claude Science evaluation.

One command, repo-relative, no network. Rescores the six replication runs and
the pilot with the frozen matched scorer at this commit
(code/61-rescore-matched.py), reproduces the identifier-leak split, and exits
non-zero if anything differs from results/expected.json.

    python evaluations/2026-09-claude-science/verify.py
"""
from __future__ import annotations

import csv
import json
import re
import statistics
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ANSWERS = "claude_science_answers.csv"
_PHENOTYPE_TOKEN = re.compile(r"_(pm|im|um|nm|rm|norm|normal|pos|neg|def|fav|unf|sens|poor|decr)(_|$)")


def _load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def scorer():
    return _load(REPO / "code" / "61-rescore-matched.py", "rescore_matched")


def load_cases() -> list[dict]:
    return json.loads((REPO / "specs" / "test_cases_v3.json").read_text())


def read_answers(path: Path, id_map: dict | None) -> dict[str, dict]:
    with Path(path).open(newline="") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for r in rows:
        cid = r["case_id"].strip()
        out[id_map.get(cid, cid) if id_map else cid] = r
    return out


def score_file(path: Path, cases: list[dict], id_map: dict | None, matched=None) -> list[dict]:
    """Missing cases count as failures, as in the paper."""
    matched = matched or scorer()
    by_id = read_answers(path, id_map)
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
    return scored


def _summ(scored: list[dict]) -> dict:
    n = len(scored)
    lethal = [s for s in scored if s["lethal_action"] is not None]
    return {
        "phenotype": sum(s["a1_phenotype"] for s in scored) / n,
        "recommendation": sum(s["a2_recommendation"] for s in scored) / n,
        "lethal_errors": sum(1 for s in lethal if s["lethal_action"] == 0.0),
        "lethal_n": len(lethal),
        "abstained": sum(s["abstained"] for s in scored),
    }


def _spread(v: list[float]) -> dict:
    return {"mean": round(statistics.fmean(v), 4), "sd": round(statistics.pstdev(v), 4),
            "min": round(min(v), 4), "max": round(max(v), 4), "per_replicate": [round(x, 4) for x in v]}


def aggregate(runs_root: Path, cases: list[dict], id_map: dict) -> dict:
    matched = scorer()
    out = {}
    for cell in sorted(p.name for p in Path(runs_root).iterdir() if p.is_dir() and not p.name.startswith("_")):
        reps = sorted(Path(runs_root, cell).glob("rep*/" + ANSWERS))
        if not reps:
            out[cell] = {"n_replicates": 0, "status": "no answers found"}
            continue
        per = [_summ(score_file(p, cases, id_map, matched)) for p in reps]
        out[cell] = {
            "n_replicates": len(per), "status": "ok",
            "phenotype": _spread([x["phenotype"] for x in per]),
            "recommendation": _spread([x["recommendation"] for x in per]),
            "lethal_errors": [x["lethal_errors"] for x in per],
            "lethal_n": per[0]["lethal_n"],
            "abstained": [x["abstained"] for x in per],
        }
    return out


def leak_split(pilot_csv: Path, hashed_csv: Path, cases: list[dict], id_map: dict) -> dict:
    """Split cases by whether the ORIGINAL identifier encodes a phenotype token,
    then compare the pilot (semantic ids) against a hashed-id run case by case."""
    matched = scorer()
    pilot = {s["case_id"]: s["a1_phenotype"] for s in score_file(pilot_csv, cases, None, matched)}
    hashed = {s["case_id"]: s["a1_phenotype"] for s in score_file(hashed_csv, cases, id_map, matched)}
    groups = {"leaked": [], "control": []}
    for c in cases:
        groups["leaked" if _PHENOTYPE_TOKEN.search(c["id"] + "_") else "control"].append(c["id"])
    res = {}
    for g, ids in groups.items():
        res[g] = {"n": len(ids),
                  "pilot": round(sum(pilot[i] for i in ids) / len(ids), 4),
                  "hashed": round(sum(hashed[i] for i in ids) / len(ids), 4)}
        res[g]["delta"] = round(res[g]["hashed"] - res[g]["pilot"], 4)
    return res


def recompute() -> dict:
    cases = load_cases()
    id_map = json.loads((HERE / "keys" / "id_map.json").read_text())
    agg = aggregate(HERE / "runs", cases, id_map)
    pilot = _summ(score_file(HERE / "pilot" / ANSWERS, cases, None))
    leak = leak_split(HERE / "pilot" / ANSWERS, HERE / "runs" / "api_corpus" / "rep1" / ANSWERS, cases, id_map)
    return {
        "cells": {c: {"phenotype_mean": v["phenotype"]["mean"], "phenotype_sd": v["phenotype"]["sd"],
                      "abstained": v["abstained"], "lethal_errors": v["lethal_errors"]}
                  for c, v in agg.items() if v.get("status") == "ok"},
        "pilot": {"phenotype": round(pilot["phenotype"], 4), "abstained": pilot["abstained"],
                  "lethal_errors": pilot["lethal_errors"]},
        "leak_split": leak,
    }


def check(expected: dict) -> list[str]:
    """Compare a fresh recompute with results/expected.json. Empty list means agreement."""
    got = recompute()
    problems = []

    def walk(e, g, path):
        if isinstance(e, dict):
            for k, ev in e.items():
                if k not in g:
                    problems.append(f"{path}/{k}: missing in recompute")
                else:
                    walk(ev, g[k], f"{path}/{k}")
        elif isinstance(e, list):
            if list(e) != list(g):
                problems.append(f"{path}: expected {e}, got {g}")
        elif isinstance(e, float):
            if abs(e - float(g)) > 5e-4:
                problems.append(f"{path}: expected {e}, got {g}")
        elif e != g:
            problems.append(f"{path}: expected {e!r}, got {g!r}")

    walk(expected, got, "")
    return problems


def main() -> int:
    got = recompute()
    print(json.dumps(got, indent=1))
    exp_path = HERE / "results" / "expected.json"
    if not exp_path.exists():
        print("\nINCOMPLETE: results/expected.json missing; nothing to check against")
        return 2
    problems = check(json.loads(exp_path.read_text()))
    if problems:
        print("\nMISMATCH against results/expected.json:")
        for p in problems:
            print("  " + p)
        return 1
    print("\nOK: every number in results/expected.json reproduces from the shipped CSVs "
          "with code/61-rescore-matched.py at this commit.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
