#!/usr/bin/env python3
"""
Fabrication firewall for the Cell Genomics v16 manuscript.

Every quantitative claim in the manuscript is recomputed here from the raw result
files and checked against the value stated in the text. Each check prints the
manuscript value, the value recomputed from data, the source file, the model set,
and PASS/FAIL. Any claim that cannot be reproduced from data fails loudly.

This is the source of truth for every number. If the manuscript and this script
disagree, the manuscript is wrong, not this script.

  python3 code/56-validate-manuscript-numbers.py

Exit code 0 if all checks pass, 1 otherwise.
"""
from __future__ import annotations
import json
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
R = BASE / "RESULTS"
three = json.loads((R / "v3_raw_rescored_three_arm.json").read_text())
skill = json.loads((R / "v3_armA9_armBv2.json").read_text())
cases = {c["id"]: c for c in json.loads((BASE / "specs" / "test_cases_v3.json").read_text())}
pop_skill = json.loads((R / "v3_armA9_armBv2_POP.json").read_text())
adv_fwd = json.loads((R / "v3_adversarial_scrambled.json").read_text())
adv_rev = json.loads((R / "v3_adversarial_reverse.json").read_text())

MISTRAL = "Mistral Large 2"
CLEAN8 = lambda r: r["model"] != MISTRAL
LETHAL = {tc for tc, c in cases.items() if "lethal" in c["gt_drug"].lower()}

CHECKS = []
def check(name, computed, expected, tol, src, modelset):
    CHECKS.append((name, computed, expected, tol, src, modelset))

# ---- three-arm aggregates (ALL 9 models, scores field) ----
def three_mean(cond, dim, clean=False):
    rows = [r for r in three if r["cond"] == cond and (CLEAN8(r) if clean else True)
            and not r["scores"].get("format_fail")]
    return 100 * sum(r["scores"][dim] for r in rows) / len(rows)

def three_lethal_errs(cond):
    rows = [r for r in three if r["cond"] == cond and not r["scores"].get("format_fail")
            and r["tc"] in LETHAL]
    return sum(1 for r in rows if r["scores"]["A3"] < 1)

def three_consistency(cond):
    g = defaultdict(list)
    for r in three:
        if r["cond"] == cond and not r["scores"].get("format_fail"):
            g[(r["model"], r["tc"], r["pop"])].append(r["scores"]["A1"])
    f = [v for v in g.values() if len(v) == 3]
    return 100 * sum(1 for v in f if len(set(v)) == 1) / len(f)

def free_model_a1(m):
    rows = [r for r in three if r["cond"] == "no_spec" and r["model"] == m
            and not r["scores"].get("format_fail")]
    return 100 * sum(r["scores"]["A1"] for r in rows) / len(rows)

models9 = sorted({r["model"] for r in three})
spread = max(free_model_a1(m) for m in models9) - min(free_model_a1(m) for m in models9)

check("Free-prompt phenotype A1", three_mean("no_spec", "A1"), 80.6, 0.3, "three-arm", "all 9")
check("RAG phenotype A1", three_mean("cpic_rag", "A1"), 89.5, 0.3, "three-arm", "all 9")
check("Control phenotype A1", three_mean("with_spec", "A1"), 100.0, 0.1, "three-arm", "all 9")
check("Free-prompt drug-rec A2", three_mean("no_spec", "A2"), 61.6, 0.5, "three-arm", "all 9")
check("RAG drug-rec A2", three_mean("cpic_rag", "A2"), 53.0, 0.5, "three-arm", "all 9")
check("Free-prompt lethal-class errors", three_lethal_errs("no_spec"), 270, 0, "three-arm", "all 9")
check("RAG lethal-class errors", three_lethal_errs("cpic_rag"), 414, 0, "three-arm", "all 9")
check("Free-prompt 3/3 consistency", three_consistency("no_spec"), 82.7, 0.5, "three-arm", "all 9")
check("RAG 3/3 consistency", three_consistency("cpic_rag"), 93.8, 0.5, "three-arm", "all 9")
check("Free-prompt model spread", spread, 27.9, 0.3, "three-arm", "all 9")
check("Weakest free model (Gemini)", free_model_a1("Gemini 2.5 Flash"), 62.8, 0.3, "three-arm", "Gemini")
check("Strongest free model (o3)", free_model_a1("o3"), 90.7, 0.3, "three-arm", "o3")

# ---- correctness-by-coincidence (ALL 9) ----
def coincidence(cond):
    rows = [r for r in three if r["cond"] == cond and not r["scores"].get("format_fail")
            and r["tc"] in LETHAL]
    ca = [r for r in rows if r["scores"]["A3"] >= 1]
    return 100 * sum(1 for r in ca if r["scores"]["A1"] < 1) / len(ca)
check("Coincidence free-prompt", coincidence("no_spec"), 19.0, 0.5, "three-arm", "all 9")
check("Coincidence RAG", coincidence("cpic_rag"), 11.5, 0.5, "three-arm", "all 9")
check("Coincidence control", coincidence("with_spec"), 0.0, 0.1, "three-arm", "all 9")

# ---- information-without-action (per-gene lethal means, ALL models) ----
def gene_cond_mean(gene, cond, dim):
    rows = [r for r in three if r["gene"] == gene and r["cond"] == cond
            and not r["scores"].get("format_fail") and r["tc"] in LETHAL]
    return sum(r["scores"][dim] for r in rows) / len(rows)
check("HLA-B*57:01 RAG A1 (phenotype)", 100*gene_cond_mean("HLA-B*57:01", "cpic_rag", "A1"), 100.0, 0.5, "three-arm", "all 9")
check("HLA-B*57:01 RAG A3 (safety)", 100*gene_cond_mean("HLA-B*57:01", "cpic_rag", "A3"), 11.1, 0.5, "three-arm", "all 9")

# ---- skill arms (CLEAN 8) ----
def skill_mean(arm, dim):
    rows = [r for r in skill if r["arm"] == arm and CLEAN8(r)]
    return 100 * sum(r[dim] for r in rows) / len(rows)
def skill_repro(arm):
    g = defaultdict(list)
    for r in skill:
        if r["arm"] == arm and CLEAN8(r):
            g[(r["model"], r["tc"])].append(r["A1"])
    f = [v for v in g.values() if len(v) == 2]
    return 100 * sum(1 for v in f if len(set(v)) == 1) / len(f)
check("Skill-reasoning A1", skill_mean("A_reasoning", "A1"), 96.2, 0.3, "skill", "clean 8")
check("Skill-execution A1", skill_mean("B_execution", "A1"), 94.3, 0.3, "skill", "clean 8")
check("Skill-reasoning reproducibility", skill_repro("A_reasoning"), 97.8, 0.5, "skill", "clean 8")

import re
def canon(s, g):
    s = (s or "").lower().replace(g.lower(), ""); s = re.sub(r"\(.*?\)", "", s).replace("hla-", "")
    s = s.replace("positive", "pos").replace("negative", "neg").replace("carrier", "")
    s = re.sub(r"[^\w*/:.\- ]", "", s); p = [re.sub(r"\s+", "", x) for x in re.split(r"\s*/\s*", s) if x.strip()]
    if len(p) == 2: p = sorted(p)
    return "/".join(p) if p else re.sub(r"\s+", "", s)
B = [r for r in skill if r["arm"] == "B_execution" and CLEAN8(r)]
corr = [r for r in B if canon(r["called_diplotype"], r["gene"]) == canon(cases[r["tc"]]["gt_diplotype"], r["gene"])]
wrong = [r for r in B if r not in corr]
check("Correct-call rate", 100*len(corr)/len(B), 93.4, 0.3, "skill", "clean 8")
check("A1 given correct call", 100*sum(r["A1"] for r in corr)/len(corr), 100.0, 0.1, "skill", "clean 8")
check("A1 given wrong call (coincidence)", 100*sum(r["A1"] for r in wrong)/len(wrong), 14.5, 0.5, "skill", "clean 8")

# ---- population sweep (CLEAN 8) ----
def pop_acc(arm, pop):
    rows = [r for r in pop_skill if r["arm"] == arm and r["pop"] == pop and CLEAN8(r)]
    return 100 * sum(r["A1"] for r in rows) / len(rows)
for p, e in [("EUR", 95.6), ("AMR", 96.0), ("AFR", 95.8)]:
    check(f"Pop reasoning {p}", pop_acc("A_reasoning", p), e, 0.3, "pop sweep", "clean 8")
for p, e in [("EUR", 93.9), ("AMR", 93.9), ("AFR", 92.2)]:
    check(f"Pop execution {p}", pop_acc("B_execution", p), e, 0.3, "pop sweep", "clean 8")

# ---- adversarial ----
def adv_exec(rows):
    return sum(1 for r in rows if r["classification"] in ("ECHO_SCRAMBLED", "HEDGE"))
def adv_revert(rows):
    return sum(1 for r in rows if r["classification"] not in ("ECHO_SCRAMBLED", "HEDGE"))
check("Adversarial forward executed", adv_exec(adv_fwd), 45, 0, "adversarial", "3 models")
check("Adversarial reverse executed", adv_exec(adv_rev), 45, 0, "adversarial", "3 models")
check("Adversarial total reverted", adv_revert(adv_fwd) + adv_revert(adv_rev), 0, 0, "adversarial", "3 models")
hedged = sum(1 for r in adv_fwd + adv_rev if r["classification"] == "HEDGE")
check("Adversarial hedged", hedged, 2, 0, "adversarial", "3 models")

# ---- dataset scale ----
check("Three-arm evaluation count", len(three), 26730, 0, "three-arm", "all 9")
check("Test cases", len(cases), 110, 0, "cases", "n/a")

# ============================================================ report
def fmt(v):
    return f"{v:.1f}" if isinstance(v, float) else str(v)

print(f"{'CLAIM':38}{'MANUSCRIPT':>11}{'FROM DATA':>11}{'':>3}{'SOURCE':>12} {'MODELS'}")
print("-" * 96)
n_pass = 0
fails = []
for name, computed, expected, tol, src, ms in CHECKS:
    ok = abs(computed - expected) <= tol
    n_pass += ok
    if not ok:
        fails.append(name)
    flag = "OK " if ok else "XX "
    print(f"{name:38}{fmt(expected):>11}{fmt(computed):>11}  {flag}{src:>12} {ms}")
print("-" * 96)
print(f"{n_pass}/{len(CHECKS)} checks reproduce from raw data.")
if fails:
    print("FAILURES (manuscript value not reproduced from data):")
    for f in fails:
        print(f"  - {f}")
    sys.exit(1)
print("All manuscript numbers are reproducible from the raw result files. No fabrication.")
sys.exit(0)
