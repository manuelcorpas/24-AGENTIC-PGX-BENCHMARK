#!/usr/bin/env python3
"""Re-run only Gemini (token-budget fix) and Mistral (hard throttle) and merge into
v3_armA9_armBv2.json, then recompute the full report. Fixes the two harness artifacts
(Gemini truncation, Mistral 429) without re-running the 7 clean models."""
import json, time, threading, requests
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import importlib.util
BASE = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("m", str(BASE/"code"/"41-armA9-armBv2.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

# fix 1: Gemini token budget (thinking model needs room for reasoning + output)
def gem_fixed(p):
    with m.sem["gemini"]:
        r = requests.post(f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={m.GKEY}",
                          json={"contents":[{"parts":[{"text":p}]}], "generationConfig":{"maxOutputTokens":2048}}, timeout=120)
        d = r.json()
        return d["candidates"][0]["content"]["parts"][0]["text"]
m._gem = gem_fixed
# fix 2: throttle Mistral to avoid 429s
m.sem["mistral"] = threading.Semaphore(2)

RERUN = {"Gemini 2.5 Flash", "Mistral Large 2"}
tasks = []
for c in m.cases:
    for arm, pf in [("A_reasoning", m.armA_prompt), ("B_execution", m.armB2_prompt)]:
        p = pf(c)
        for model in RERUN:
            for rep in range(m.N_REPS):
                tasks.append((c, arm, p, model, m.MODELS[model], rep))
print(f"re-running {len(tasks)} calls for {RERUN}", flush=True)
new = []; done = 0
with ThreadPoolExecutor(max_workers=12) as ex:
    for r in ex.map(m.run_one, tasks):
        new.append(r); done += 1
        if done % 100 == 0: print(f"  [{done}/{len(tasks)}]", flush=True)

# merge: drop old Gemini/Mistral rows, add fixed ones
old = json.loads((BASE/"RESULTS"/"v3_armA9_armBv2.json").read_text())
merged = [r for r in old if r["model"] not in RERUN] + new
(BASE/"RESULTS"/"v3_armA9_armBv2.json").write_text(json.dumps(merged, indent=2))
results = merged

# recompute report
models = list(m.MODELS.keys())
def mean(arm, dim, model=None):
    rows = [r for r in results if r["arm"]==arm and (model is None or r["model"]==model)]
    return 100*sum(r[dim] for r in rows)/len(rows) if rows else float("nan")
def repro(arm, model=None):
    g = defaultdict(list)
    for r in results:
        if r["arm"]==arm and (model is None or r["model"]==model): g[(r["model"],r["tc"])].append(r["A1"])
    full=[v for v in g.values() if len(v)==m.N_REPS]
    return 100*sum(1 for v in full if len(set(v))==1)/len(full) if full else float("nan")
def lethal(arm, model=None):
    return sum(1 for r in results if r["arm"]==arm and r["lethal"] and r["A3"]<1 and (model is None or r["model"]==model))
out=["ARM A (skill-reasoning, 9 models) + ARM B v2 (skill-execution) — Gemini/Mistral fixed",
     f"9 models x 110 cases x EUR x {m.N_REPS} reps per arm",""]
out.append("=== AGGREGATE (all 9) ===")
out.append(f"{'metric':24}{'Arm A':>12}{'Arm B v2':>12}")
for d in ["A1","A2","A3"]: out.append(f"{'aggregate '+d:24}{mean('A_reasoning',d):11.1f}%{mean('B_execution',d):11.1f}%")
out.append(f"{'reproducibility(A1)':24}{repro('A_reasoning'):11.1f}%{repro('B_execution'):11.1f}%")
out.append(f"{'lethal-class A3 errors':24}{lethal('A_reasoning'):12d}{lethal('B_execution'):12d}")
bmatch=[r for r in results if r["arm"]=="B_execution"]
out.append(f"Arm B v2 controlled-vocab match: {100*sum(1 for r in bmatch if r.get('matched'))/len(bmatch):.1f}%")
out.append("")
out.append("=== PER MODEL ===")
out.append(f"{'model':18}{'A_A1':>8}{'A_repro':>9}{'B_A1':>8}{'B_repro':>9}{'A_leth':>8}{'B_leth':>8}")
for mm in models:
    out.append(f"{mm:18}{mean('A_reasoning','A1',mm):7.1f}%{repro('A_reasoning',mm):8.1f}%{mean('B_execution','A1',mm):7.1f}%{repro('B_execution',mm):8.1f}%{lethal('A_reasoning',mm):8d}{lethal('B_execution',mm):8d}")
txt="\n".join(out); (BASE/"RESULTS"/"v3_armA9_armBv2_report.txt").write_text(txt); print("\n"+txt)
