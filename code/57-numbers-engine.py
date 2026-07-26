#!/usr/bin/env python3
"""
v26 revision numbers engine.

STAGE 1  REPRODUCTION GATE   reproduce the v25 headline point estimates exactly
STAGE 2  CLUSTERED INFERENCE recompute CIs and the two key contrasts with a
                             case-level cluster bootstrap (+ model-level sensitivity),
                             report design effects and effective N
STAGE 3  MISSING COUNTS      nine-model skill aggregate (two honest treatments of
                             Mistral's rate-limit errors), skill-arm lethal-class
                             counts/denominators, Mistral usable-vs-error reconciliation

Headline phenotype accuracy = MEAN A1 over parsed cells (format_fail==False), which
reproduces v25 exactly from v3_raw_rescored_three_arm.json (free 80.6, RAG 89.5,
control 100) and from v3_armAB_fullgrid.json for the skill arms (95.5, 93.3).

Inputs (../RESULTS):  v3_raw_rescored_three_arm.json, v3_armAB_fullgrid.json,
                      v3_three_arm_per_case_a1.csv
Output: v26_numbers.json, v26_numbers.txt (this folder). No network. Seeded bootstrap.
"""
import os, json, csv, math, random
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "data")
SEED = 20260611
random.seed(SEED); np.random.seed(SEED)
B = 10000

def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)

def wilson(k, n, z=1.96):
    if n == 0: return (float("nan"), float("nan"))
    p = k/n; d = 1 + z*z/n; c = p + z*z/(2*n)
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return ((c-h)/d, (c+h)/d)

# ---- lethal flag per case ----
lethal_case = {}
with open(os.path.join(RES, "v3_three_arm_per_case_a1.csv")) as f:
    for row in csv.DictReader(f):
        lethal_case[row["case_id"]] = int(row["lethal"])
n_lethal_cases = sum(lethal_case.values())

# ---- core conditions from the canonical merged three-arm file ----
def core_cells(cond):
    out = []
    for r in load("v3_raw_rescored_three_arm.json"):
        if r["cond"] != cond: continue
        sc = r.get("scores") or {}
        out.append(dict(cond=cond, tc=r["tc"], model=r["model"], pop=r["pop"],
                        lethal=lethal_case.get(r["tc"], 0),
                        parsed=(not sc.get("format_fail", False)),
                        A1=sc.get("A1"), A3=sc.get("A3")))
    return out
free_ce    = core_cells("no_spec")
rag        = core_cells("cpic_rag")
control_ce = core_cells("with_spec")

# ---- skill arms from the full grid ----
fg = load("v3_armAB_fullgrid.json")
def skill_cells(arm, include_mistral):
    out = []
    for r in fg:
        if r["arm"] != arm: continue
        if (not include_mistral) and "istral" in r["model"]: continue
        raw = (r.get("raw_text") or "")
        out.append(dict(cond=arm, tc=r["tc"], model=r["model"], pop=r["pop"],
                        lethal=bool(r.get("lethal")),
                        parsed=(raw.strip() not in ("", "<error: 'choices'>")),
                        A1=r.get("A1"), A3=r.get("A3"), raw=raw))
    return out
reason8 = skill_cells("A_reasoning", False); exec8 = skill_cells("B_execution", False)
reason9 = skill_cells("A_reasoning", True);  exec9 = skill_cells("B_execution", True)

# ---- metrics (headline = mean A1 over parsed) ----
def mean_A1(cells):
    v = [c["A1"] for c in cells if c["parsed"] and c["A1"] is not None]
    return (sum(v)/len(v) if v else float("nan")), len(v)
def lethal_A3_errors(cells):
    leth = [c for c in cells if c["lethal"] and c["parsed"] and c["A3"] is not None]
    return sum(1 for c in leth if c["A3"] < 1.0), len(leth)

# ===================== STAGE 1: reproduction gate =====================
gate = {}
for name, cells, target in [("free_prompted_A1", free_ce, 80.6), ("retrieval_A1", rag, 89.5),
                            ("control_A1", control_ce, 100.0),
                            ("skill_reasoning_A1_8mdl", reason8, 95.5),
                            ("skill_execution_A1_8mdl", exec8, 93.3)]:
    m, n = mean_A1(cells)
    gate[name] = dict(pct=round(100*m, 1), n=n, target=target, ok=abs(100*m-target) <= 0.15)
for name, cells in [("free_prompted_lethal", free_ce), ("retrieval_lethal", rag)]:
    e, n = lethal_A3_errors(cells)
    gate[name] = dict(errors=e, n=n, rate=round(100*e/n, 1) if n else None)

# ===================== STAGE 2: clustered inference =====================
def by_cluster(cells, key):
    d = defaultdict(list)
    for c in cells: d[c[key]].append(c)
    return d

def boot_mean_A1(cells, key):
    groups = list(by_cluster(cells, key).values()); G = len(groups)
    ests = []
    for _ in range(B):
        idx = np.random.randint(0, G, G); s = 0.0; n = 0
        for i in idx:
            for c in groups[i]:
                if c["parsed"] and c["A1"] is not None:
                    s += c["A1"]; n += 1
        if n: ests.append(s/n)
    return np.percentile(ests, [2.5, 97.5])

def boot_lethal_diff(a, b, key):
    ga = by_cluster([c for c in a if c["lethal"]], key)
    gb = by_cluster([c for c in b if c["lethal"]], key)
    keys = sorted(set(ga) | set(gb)); K = len(keys)
    ea, na = lethal_A3_errors(a); eb, nb = lethal_A3_errors(b)
    diffs = []
    for _ in range(B):
        idx = np.random.randint(0, K, K); ea_=na_=eb_=nb_=0
        for i in idx:
            kk = keys[i]
            for c in ga.get(kk, []):
                if c["parsed"] and c["A3"] is not None: na_+=1; ea_+=(c["A3"]<1.0)
            for c in gb.get(kk, []):
                if c["parsed"] and c["A3"] is not None: nb_+=1; eb_+=(c["A3"]<1.0)
        if na_ and nb_: diffs.append(eb_/nb_ - ea_/na_)
    diffs = np.array(diffs); lo, hi = np.percentile(diffs, [2.5, 97.5])
    return dict(rate_a=round(100*ea/na,1), rate_b=round(100*eb/nb,1), err_a=ea, n_a=na,
                err_b=eb, n_b=nb, diff_pp=round(100*(eb/nb-ea/na),1),
                ci_pp=[round(100*lo,1), round(100*hi,1)],
                boot_p_two_sided=round(2*min(float(np.mean(diffs<=0)), float(np.mean(diffs>=0))),4),
                n_clusters=K)

clustered = {}
for name, cells in [("free_prompted", free_ce), ("retrieval", rag), ("control", control_ce),
                    ("skill_reasoning_8mdl", reason8), ("skill_execution_8mdl", exec8)]:
    m, n = mean_A1(cells)
    k = round(m*n)  # for a naive Wilson reference treating the mean as a binomial proportion
    nw = wilson(k, n)
    lo_c, hi_c = boot_mean_A1(cells, "tc")
    lo_m, hi_m = boot_mean_A1(cells, "model")
    nw_w = (nw[1]-nw[0]); cl_w = (hi_c-lo_c)
    clustered[name] = dict(point=round(100*m,1), n_cells=n,
                           n_cases=len(set(c["tc"] for c in cells)),
                           naive_wilson=[round(100*nw[0],1), round(100*nw[1],1)],
                           cluster_by_case=[round(100*lo_c,1), round(100*hi_c,1)],
                           cluster_by_model=[round(100*lo_m,1), round(100*hi_m,1)],
                           design_effect=round((cl_w/nw_w)**2,1) if nw_w>0 else None)
lethal_contrast = boot_lethal_diff(free_ce, rag, "tc")

# ===================== STAGE 3: missing counts =====================
def per_model(cells):
    d = defaultdict(lambda: dict(sum=0.0, n_parsed=0, n_total=0, err=0))
    for c in cells:
        if c["A1"] is None: continue
        d[c["model"]]["n_total"] += 1
        if c["parsed"]:
            d[c["model"]]["n_parsed"] += 1; d[c["model"]]["sum"] += c["A1"]
        else:
            d[c["model"]]["err"] += 1
    out = {}
    for m, v in sorted(d.items()):
        overall = v["sum"]/v["n_total"] if v["n_total"] else float("nan")   # error = wrong
        parsedonly = v["sum"]/v["n_parsed"] if v["n_parsed"] else float("nan")  # error = missing
        out[m] = dict(overall_pct=round(100*overall,1), parsed_only_pct=round(100*parsedonly,1),
                      n_parsed=v["n_parsed"], n_total=v["n_total"],
                      usable_pct=round(100*v["n_parsed"]/v["n_total"],1))
    return out

def agg(cells, parsed_only):
    v = [c["A1"] for c in cells if c["A1"] is not None and (c["parsed"] or not parsed_only)]
    return round(100*sum(v)/len(v),1), len(v)

def lethal_counts(cells):
    e, n = lethal_A3_errors(cells); return dict(errors=e, n=n, rate=round(100*e/n,1) if n else None)

stage3 = dict(
    nine_model_skill=dict(
        reasoning_error_as_wrong=dict(zip(("pct","n"), agg(reason9, False))),
        reasoning_parsed_only=dict(zip(("pct","n"), agg(reason9, True))),
        execution_error_as_wrong=dict(zip(("pct","n"), agg(exec9, False))),
        execution_parsed_only=dict(zip(("pct","n"), agg(exec9, True))),
        eight_model_reasoning=gate["skill_reasoning_A1_8mdl"]["pct"],
        eight_model_execution=gate["skill_execution_A1_8mdl"]["pct"],
        per_model_reasoning=per_model(reason9),
        per_model_execution=per_model(exec9),
    ),
    mistral=dict(
        reasoning_usable_pct=per_model(reason9)["Mistral Large 2"]["usable_pct"],
        execution_usable_pct=per_model(exec9)["Mistral Large 2"]["usable_pct"],
        reasoning_accuracy_when_usable=per_model(reason9)["Mistral Large 2"]["parsed_only_pct"],
        execution_accuracy_when_usable=per_model(exec9)["Mistral Large 2"]["parsed_only_pct"],
    ),
    skill_arm_lethal=dict(
        skill_reasoning_8mdl=lethal_counts(reason8),
        skill_execution_8mdl=lethal_counts(exec8),
    ),
)

out = dict(meta=dict(B=B, seed=SEED, n_total_cases=len(lethal_case), n_lethal_cases=n_lethal_cases),
           stage1_reproduction_gate=gate, stage2_clustered=clustered,
           stage2_lethal_contrast_free_to_rag=lethal_contrast, stage3_missing_counts=stage3)
json.dump(out, open(os.path.join(HERE, "v26_numbers.json"), "w"), indent=2)

# ---- human-readable ----
L = [f"V26 NUMBERS ENGINE  (mean-A1 metric; cluster bootstrap B={B}, seed={SEED})", "="*72]
L.append("\nSTAGE 1  REPRODUCTION GATE")
allok = True
for k, v in gate.items():
    if "ok" in v:
        allok &= v["ok"]
        L.append(f"  {k:30s} {v['pct']}%  target {v['target']}%  {'OK' if v['ok'] else 'MISMATCH'}  (n={v['n']})")
    else:
        L.append(f"  {k:30s} {v['errors']}/{v['n']} = {v['rate']}%")
L.append(f"  >>> GATE {'PASSED' if allok else 'FAILED'}")
L.append("\nSTAGE 2  ACCURACY CIs  (naive Wilson vs cluster bootstrap)")
for name, v in clustered.items():
    L.append(f"  {name}: {v['point']}%  (n_cells={v['n_cells']}, n_cases={v['n_cases']})")
    L.append(f"      naive Wilson     {v['naive_wilson']}  width {round(v['naive_wilson'][1]-v['naive_wilson'][0],1)}pp")
    L.append(f"      cluster-by-case  {v['cluster_by_case']}  -> design effect x{v['design_effect']}")
    L.append(f"      cluster-by-model {v['cluster_by_model']}")
lc = lethal_contrast
L.append("\nSTAGE 2  LETHAL-CLASS CONTRAST free -> retrieval (clustered by case)")
L.append(f"  free {lc['err_a']}/{lc['n_a']}={lc['rate_a']}%  retrieval {lc['err_b']}/{lc['n_b']}={lc['rate_b']}%")
L.append(f"  diff +{lc['diff_pp']}pp  95% CI {lc['ci_pp']}pp  two-sided bootstrap p={lc['boot_p_two_sided']}  ({lc['n_clusters']} lethal clusters; {n_lethal_cases} lethal cases)")
L.append("\nSTAGE 3  NINE-MODEL SKILL AGGREGATE")
nm = stage3["nine_model_skill"]
L.append(f"  reasoning: 8-model {nm['eight_model_reasoning']}% | 9-model error-as-wrong {nm['reasoning_error_as_wrong']['pct']}% | 9-model parsed-only {nm['reasoning_parsed_only']['pct']}%")
L.append(f"  execution: 8-model {nm['eight_model_execution']}% | 9-model error-as-wrong {nm['execution_error_as_wrong']['pct']}% | 9-model parsed-only {nm['execution_parsed_only']['pct']}%")
L.append("  per-model (reasoning) overall% / parsed-only% / usable%:")
for m, v in nm["per_model_reasoning"].items():
    L.append(f"      {m:20s} {v['overall_pct']:5}% / {v['parsed_only_pct']:5}% / {v['usable_pct']:5}%")
ms = stage3["mistral"]
L.append(f"  Mistral usable: reasoning {ms['reasoning_usable_pct']}% (acc when usable {ms['reasoning_accuracy_when_usable']}%); execution {ms['execution_usable_pct']}% (acc when usable {ms['execution_accuracy_when_usable']}%)")
L.append("  skill-arm lethal-class counts:")
for k, v in stage3["skill_arm_lethal"].items():
    L.append(f"      {k:24s} {v['errors']}/{v['n']} = {v['rate']}%")
open(os.path.join(HERE, "v26_numbers.txt"), "w").write("\n".join(L))
print("\n".join(L))
