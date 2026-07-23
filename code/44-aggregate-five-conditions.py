#!/usr/bin/env python3
"""
Unified five-condition aggregation on the SAME 3-pop x 3-rep x 9-model grid.

Conditions 1/2/5 (no_spec, cpic_rag, with_spec) come from the three-arm rescored
file (clinical-equivalence scorer = the `scores` dict, which reproduces the
manuscript headline: no_spec A1 80.6%, cpic_rag 89.5%, with_spec 100%).

Conditions 3/4 (skill-reasoning = A_reasoning, skill-execution = B_execution) come
from the new full-grid run (43-armAB-fullgrid.py). To match the headline scorer,
the skill phenotype (pphen) is re-scored with the SAME clinical-equivalence A1
layer (10b). A2/A3/lethal use the skill runner's deterministic scores.

Reports every metric two ways: ALL 9 MODELS and EX-MISTRAL (8 models, the headline),
the Mistral exclusion now applied uniformly to all five conditions.

Run: python3 code/44-aggregate-five-conditions.py
"""
from __future__ import annotations
import json
from pathlib import Path
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec

BASE = Path(__file__).resolve().parent.parent
THREE = BASE / "RESULTS" / "v3_raw_rescored_three_arm.json"
SKILL = BASE / "RESULTS" / "v3_armAB_fullgrid.json"
SKILL_JSONL = BASE / "RESULTS" / "v3_armAB_fullgrid.jsonl"
CASES = BASE / "specs" / "test_cases_v3.json"
OUT_REPORT = BASE / "RESULTS" / "v3_five_condition_matched_report.txt"
OUT_JSON = BASE / "RESULTS" / "v3_five_condition_matched.json"

MISTRAL = "Mistral Large 2"

# clinical-equivalence A1 scorer (same one used for the three-arm headline)
_sp = spec_from_file_location("ce", str(BASE / "code" / "10b-rescore-v3-clinical-equivalence.py"))
ce = module_from_spec(_sp); _sp.loader.exec_module(ce)

cases = {c["id"]: c for c in json.loads(CASES.read_text())}


def load_skill():
    if SKILL.exists():
        rows = json.loads(SKILL.read_text())
    else:
        rows = []
        for line in SKILL_JSONL.read_text().splitlines():
            if line.strip():
                rows.append(json.loads(line))
        # dedupe by cell, last wins
        dd = {}
        for r in rows:
            dd[(r["tc"], r["arm"], r["model"], r["rep"], r["pop"])] = r
        rows = list(dd.values())
    # re-score A1 with clinical equivalence for parity with the three-arm headline
    for r in rows:
        gt = cases[r["tc"]]["gt_phenotype"]
        r["A1ce"] = ce.score_a1_clinical_eq(r["pphen"], gt, r["gene"])
    return rows


def load_three():
    return json.loads(THREE.read_text())


def main():
    skill = load_skill()
    three = load_three()

    # normalise into a common record list: cond, model, pop, tc, rep/run, A1, A2, A3, lethal, lethal_err
    recs = []
    for r in three:
        if not r["parsed"]:
            continue
        s = r["scores"]  # clinical-equivalence
        gt_drug = cases[r["tc"]]["gt_drug"]
        recs.append({"cond": r["cond"], "model": r["model"], "pop": r["pop"], "tc": r["tc"],
                     "rep": r["run"], "A1": s["A1"], "A2": s["A2"], "A3": s["A3"],
                     "lethal": "lethal" in gt_drug.lower()})
    cond_name = {"A_reasoning": "skill_reasoning", "B_execution": "skill_execution"}
    for r in skill:
        recs.append({"cond": cond_name[r["arm"]], "model": r["model"], "pop": r["pop"], "tc": r["tc"],
                     "rep": r["rep"], "A1": r["A1ce"], "A2": r["A2"], "A3": r["A3"],
                     "lethal": r["lethal"]})

    CONDS = ["no_spec", "cpic_rag", "skill_reasoning", "skill_execution", "with_spec"]
    LABEL = {"no_spec": "1 free-prompted", "cpic_rag": "2 retrieval", "skill_reasoning": "3 skill-reasoning",
             "skill_execution": "4 skill-execution", "with_spec": "5 answer-supplied"}

    def sel(cond, pop=None, model=None, ex_mistral=False):
        return [x for x in recs if x["cond"] == cond
                and (pop is None or x["pop"] == pop)
                and (model is None or x["model"] == model)
                and (not ex_mistral or x["model"] != MISTRAL)]

    def mean(rows, dim):
        return 100 * sum(r[dim] for r in rows) / len(rows) if rows else float("nan")

    def lethal_err(rows):
        return sum(1 for r in rows if r["lethal"] and r["A3"] < 1)

    def repro(cond, ex_mistral=True):
        g = defaultdict(list)
        for x in sel(cond, ex_mistral=ex_mistral):
            g[(x["model"], x["tc"], x["pop"])].append(x["A1"])
        full = [v for v in g.values() if len(v) >= 3]
        return 100 * sum(1 for v in full if len(set(v[:3])) == 1) / len(full) if full else float("nan")

    out = ["FIVE-CONDITION MATCHED COMPARISON (3 populations x 3 replicates x 110 cases)",
           "Conditions 1/2/5 clinical-equivalence scorer; conditions 3/4 re-scored with the same A1 layer.",
           f"skill rows loaded: {len(skill)} | three-arm parsed rows: {sum(1 for r in three if r['parsed'])}", ""]

    for tag, exm in [("ALL 9 MODELS", False), ("EX-MISTRAL (8 models, HEADLINE)", True)]:
        out.append(f"================= {tag} =================")
        out.append(f"{'condition':20}{'A1%':>8}{'A2%':>8}{'A3%':>8}{'lethalErr':>11}{'repro3/3':>10}{'n':>8}")
        for c in CONDS:
            rows = sel(c, ex_mistral=exm)
            out.append(f"{LABEL[c]:20}{mean(rows,'A1'):8.1f}{mean(rows,'A2'):8.1f}{mean(rows,'A3'):8.1f}"
                       f"{lethal_err(rows):11d}{repro(c, ex_mistral=exm):9.1f}%{len(rows):8d}")
        out.append("")
        out.append(f"--- {tag}: phenotype A1 by population (EUR / AMR / AFR / spread) ---")
        for c in CONDS:
            v = {p: mean(sel(c, pop=p, ex_mistral=exm), 'A1') for p in ("EUR", "AMR", "AFR")}
            sp = max(v.values()) - min(v.values())
            out.append(f"{LABEL[c]:20}{v['EUR']:8.1f}{v['AMR']:8.1f}{v['AFR']:8.1f}{sp:8.1f}")
        out.append("")
        out.append(f"--- {tag}: lethal-class A3 errors by population (EUR / AMR / AFR) ---")
        for c in CONDS:
            out.append(f"{LABEL[c]:20}" + "".join(f"{lethal_err(sel(c, pop=p, ex_mistral=exm)):8d}" for p in ("EUR", "AMR", "AFR")))
        out.append("")

    txt = "\n".join(out)
    OUT_REPORT.write_text(txt)
    # machine-readable headline (ex-Mistral)
    head = {}
    for c in CONDS:
        rows = sel(c, ex_mistral=True)
        head[c] = {"A1": round(mean(rows, 'A1'), 1), "A2": round(mean(rows, 'A2'), 1),
                   "A3": round(mean(rows, 'A3'), 1), "lethal_err": lethal_err(rows),
                   "repro": round(repro(c), 1), "n": len(rows),
                   "by_pop": {p: round(mean(sel(c, pop=p, ex_mistral=True), 'A1'), 1) for p in ("EUR", "AMR", "AFR")}}
    OUT_JSON.write_text(json.dumps(head, indent=2))
    print(txt)
    print("\n-> wrote", OUT_REPORT.name, "and", OUT_JSON.name)


if __name__ == "__main__":
    main()
