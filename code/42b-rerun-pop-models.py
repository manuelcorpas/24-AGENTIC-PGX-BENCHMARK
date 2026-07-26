#!/usr/bin/env python3
"""
Targeted re-run of specific models for the population sweep, merging into
data/v3_armA9_armBv2_POP.json (replacing only the named models' rows).

Fixes carried over the main run that 42 missed:
  - Gemini 2.5 Flash is a thinking model: maxOutputTokens must be 2048, not 320
    (the 320 cap truncated all reasoning-arm output to empty). See 41b.
  - Empty/error responses are retried with longer backoff, and the run aborts
    cleanly if a provider returns a hard quota error (429 "exceeded your quota").
  - Task order is shuffled deterministically so any residual throttling is spread
    across populations rather than concentrated on the one that runs last.

Usage:
  python3 code/42b-rerun-pop-models.py "Gemini 2.5 Flash"
  python3 code/42b-rerun-pop-models.py "o3" "o4-mini" "GPT-5.2" "GPT-4.1"   # after OpenAI top-up
"""
from __future__ import annotations
import sys
import json
import time
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from importlib.util import spec_from_file_location, module_from_spec

BASE = Path(__file__).resolve().parent.parent
POP = BASE / "data" / "v3_armA9_armBv2_POP.json"

_sw = spec_from_file_location("sw", str(BASE / "code" / "42-armAB-population-sweep.py"))
sw = module_from_spec(_sw)
_sw.loader.exec_module(sw)          # loads sweep module
sw._sp.loader.exec_module(sw.h)     # loads underlying harness (clients, rules)
h = sw.h


def gem_fixed(p):
    # Gemini 2.5 Flash is a thinking model; 4096 leaves room for thinking plus the
    # short structured output so finishReason is STOP, not MAX_TOKENS (which truncates
    # the response before the PHENOTYPE line and yields an empty parse).
    with h.sem["gemini"]:
        r = requests.post(
            f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={h.GKEY}",
            json={"contents": [{"parts": [{"text": p}]}],
                  "generationConfig": {"maxOutputTokens": 4096}}, timeout=120)
        j = r.json()
        parts = j["candidates"][0]["content"]["parts"]
        return "".join(pt.get("text", "") for pt in parts)


# corrected per-model adapters (Gemini gets 4096 tokens, all parts concatenated)
ADAPTERS = dict(h.MODELS)
ADAPTERS["Gemini 2.5 Flash"] = gem_fixed


def run_cell(task):
    c, arm, prompt, model, fn, rep, pop = task
    text = ""
    for attempt in range(5):
        try:
            text = fn(prompt)
            if text and text.strip():
                break
        except Exception as e:
            msg = str(e)
            if "exceeded your current quota" in msg or "insufficient_quota" in msg:
                return {"_quota_dead": model, "_err": msg[:120]}
            text = f"<error: {e}>"
        time.sleep(3 * (attempt + 1))
    r = h.run_one((c, arm, prompt, model, fn, rep))
    r["pop"] = pop
    return r


def main(models):
    bad = [m for m in models if m not in ADAPTERS]
    if bad:
        sys.exit(f"unknown models: {bad}")
    tasks = []
    for pop in sw.POPULATIONS:
        for c in h.cases:
            for arm, pf in [("A_reasoning", sw.armA_prompt_pop), ("B_execution", sw.armB2_prompt_pop)]:
                p = pf(c, pop)
                for model in models:
                    tasks.append((c, arm, p, model, ADAPTERS[model], 0, pop["id"]))
    # deterministic shuffle (interleave populations); avoids confounding pop with time
    tasks.sort(key=lambda t: (hash((t[0]["id"], t[6], t[1])) & 0xffff))
    total = len(tasks)
    print(f"re-running {models}: {total} calls (shuffled, Gemini=2048tok)", flush=True)

    results, done, quota_dead = [], 0, set()
    with ThreadPoolExecutor(max_workers=12) as ex:
        for r in ex.map(run_cell, tasks):
            done += 1
            if r.get("_quota_dead"):
                quota_dead.add(r["_quota_dead"])
                if done % 50 == 0:
                    print(f"  [{done}/{total}] QUOTA DEAD: {r['_quota_dead']}", flush=True)
                continue
            results.append(r)
            if done % 200 == 0:
                print(f"  [{done}/{total}]", flush=True)

    if quota_dead:
        print(f"\nABORTED for quota-exhausted providers: {sorted(quota_dead)}. "
              f"No merge for those models.", flush=True)
        results = [r for r in results if r["model"] not in quota_dead]
    if not results:
        print("nothing to merge.", flush=True)
        return

    merged = json.loads(POP.read_text())
    rerun_models = {r["model"] for r in results}
    merged = [r for r in merged if r["model"] not in rerun_models]   # drop stale rows
    merged.extend(results)
    POP.write_text(json.dumps(merged, indent=2))
    print(f"merged {len(results)} fresh rows for {sorted(rerun_models)} -> {POP.name}", flush=True)


if __name__ == "__main__":
    models = sys.argv[1:] or ["Gemini 2.5 Flash"]
    main(models)
