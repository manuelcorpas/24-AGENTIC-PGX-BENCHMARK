#!/usr/bin/env python3
"""Paced Mistral-only run (serial, ~1s spacing to avoid 429 rate-limit artifact),
both arms, merged into v3_armA9_armBv2.json, full 9-model report recomputed."""
import json, time, threading, requests
from pathlib import Path
from collections import defaultdict
import importlib.util
BASE = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("m", str(BASE/"code"/"41-armA9-armBv2.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

_lock = threading.Lock(); _last = [0.0]
def mis_paced(p):
    with _lock:
        dt = time.time() - _last[0]
        if dt < 1.0: time.sleep(1.0 - dt)
        _last[0] = time.time()
    r = requests.post("https://api.mistral.ai/v1/chat/completions",
                      headers={"Authorization": f"Bearer {m.MKEY}", "Content-Type": "application/json"},
                      json={"model":"mistral-large-latest","messages":[{"role":"user","content":p}],"max_tokens":320}, timeout=120)
    return r.json()["choices"][0]["message"]["content"]
m._mis = mis_paced

M = "Mistral Large 2"
tasks = []
for c in m.cases:
    for arm, pf in [("A_reasoning", m.armA_prompt), ("B_execution", m.armB2_prompt)]:
        p = pf(c)
        for rep in range(m.N_REPS):
            tasks.append((c, arm, p, M, m.MODELS[M], rep))
print(f"paced Mistral run: {len(tasks)} calls (serial, ~1s spacing)", flush=True)
new = []
for i, t in enumerate(tasks, 1):
    new.append(m.run_one(t))
    if i % 50 == 0: print(f"  [{i}/{len(tasks)}]", flush=True)

old = json.loads((BASE/"RESULTS"/"v3_armA9_armBv2.json").read_text())
results = [r for r in old if r["model"] != M] + new
(BASE/"RESULTS"/"v3_armA9_armBv2.json").write_text(json.dumps(results, indent=2))

models = list(m.MODELS.keys())
def mean(arm, dim, model=None):
    rows=[r for r in results if r["arm"]==arm and (model is None or r["model"]==model)]
    return 100*sum(r[dim] for r in rows)/len(rows) if rows else float("nan")
def repro(arm, model=None):
    g=defaultdict(list)
    for r in results:
        if r["arm"]==arm and (model is None or r["model"]==model): g[(r["model"],r["tc"])].append(r["A1"])
    f=[v for v in g.values() if len(v)==m.N_REPS]; return 100*sum(1 for v in f if len(set(v))==1)/len(f) if f else float("nan")
def leth(arm, model=None):
    return sum(1 for r in results if r["arm"]==arm and r["lethal"] and r["A3"]<1 and (model is None or r["model"]==model))
def agg_ex_mistral(arm,dim):
    rows=[r for r in results if r["arm"]==arm and r["model"]!=M]; return 100*sum(r[dim] for r in rows)/len(rows)
out=["ARM A (skill-reasoning) + ARM B v2 (skill-execution) — FULL 9 MODELS (Mistral paced)",
     f"9 models x 110 cases x EUR x {m.N_REPS} reps per arm",""]
out.append("=== PER MODEL ===")
out.append(f"{'model':18}{'A_A1':>8}{'A_repro':>9}{'B_A1':>8}{'B_repro':>9}{'A_leth':>8}{'B_leth':>8}")
for mm in models:
    out.append(f"{mm:18}{mean('A_reasoning','A1',mm):7.1f}%{repro('A_reasoning',mm):8.1f}%{mean('B_execution','A1',mm):7.1f}%{repro('B_execution',mm):8.1f}%{leth('A_reasoning',mm):8d}{leth('B_execution',mm):8d}")
out.append("")
out.append(f"AGGREGATE all 9:        Arm A A1={mean('A_reasoning','A1'):.1f}%  Arm B A1={mean('B_execution','A1'):.1f}%")
out.append(f"AGGREGATE ex-Mistral:   Arm A A1={agg_ex_mistral('A_reasoning','A1'):.1f}%  Arm B A1={agg_ex_mistral('B_execution','A1'):.1f}%")
txt="\n".join(out); (BASE/"RESULTS"/"v3_armA9_armBv2_report.txt").write_text(txt); print("\n"+txt)
