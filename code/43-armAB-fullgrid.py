#!/usr/bin/env python3
"""
Full-grid skill arms (Arm A skill-reasoning, Arm B skill-execution v2) so that
conditions 3 and 4 are collected on the SAME grid as the three core conditions
(no_spec / cpic_rag / with_spec): 9 models x 3 populations x 3 replicates x 110
cases x 2 arms = 17,820 calls.

Why: the headline skill numbers were previously EUR-only / 2-rep / 8-model and the
population sweep was 3-pop / 1-rep / 8-model. Reviewers correctly note the skill
conditions were not on the same footing as the others. This run puts them on the
identical 3-pop x 3-rep x 9-model grid. Mistral Large 2 is INCLUDED in collection
(so collection matches the others exactly); the headline applies the Mistral
exclusion uniformly across all five conditions at the analysis stage.

Population framing is identical to 02-run-benchmark-v3.py / 42-armAB-population-sweep.py.
Model adapters, skill rules, controlled-vocabulary prompts and the baseline scorer
are reused from 41-armA9-armBv2.py.

Robustness: every completed call is appended to a JSONL checkpoint immediately, so a
crash/interruption loses nothing and a re-run resumes (skips completed cells).

Run:    python3 PYTHON/43-armAB-fullgrid.py
Resume: python3 PYTHON/43-armAB-fullgrid.py            (auto-skips checkpointed cells)
Smoke:  python3 PYTHON/43-armAB-fullgrid.py --smoke    (counts + sample prompts, no API)
Test:   python3 PYTHON/43-armAB-fullgrid.py --test1    (one live call per model)
"""
from __future__ import annotations
import json, sys, time, threading
from pathlib import Path
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from concurrent.futures import ThreadPoolExecutor

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "RESULTS" / "v3_armAB_fullgrid.json"
JSONL = BASE / "RESULTS" / "v3_armAB_fullgrid.jsonl"
REPORT = BASE / "RESULTS" / "v3_armAB_fullgrid_report.txt"
N_REPS = 3

# reuse the validated 41 harness (MODELS, cases, DIP2PHEN, REC, rules_dip,
# norm_dip, parse_field, skill_rules_text, rs scorer)
_sp = spec_from_file_location("h41", str(BASE / "PYTHON" / "41-armA9-armBv2.py"))
h = module_from_spec(_sp)
_sp.loader.exec_module(h)

# --- harness artifact fixes (identical to 41b-rerun-gemini-mistral.py) ---
# Gemini 2.5 Flash is a thinking model; 320 output tokens get consumed by reasoning
# and it returns empty text. Give it 2048. Mistral throttled to 2 to avoid 429s.
import requests as _rq
import threading as _th

def _gem_fixed(p):
    with h.sem["gemini"]:
        r = _rq.post(
            f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={h.GKEY}",
            json={"contents": [{"parts": [{"text": p}]}], "generationConfig": {"maxOutputTokens": 2048}},
            timeout=120)
        return r.json()["candidates"][0]["content"]["parts"][0]["text"]

h._gem = _gem_fixed
h.sem["mistral"] = _th.Semaphore(2)

# population framing, identical to 42-armAB-population-sweep.py / three-arm benchmark
POPULATIONS = [
    {"id": "EUR", "name": "European family cohort (Corpasome project)",
     "desc": "European ancestry, whole-genome sequencing"},
    {"id": "AMR", "name": "Peruvian Genome Project",
     "desc": "Admixed Latin American, 7 indigenous and mestizo subpopulations"},
    {"id": "AFR", "name": "Uganda Genome Resource",
     "desc": "East African, 6,407 whole-genome sequences"},
]
EXCLUDE_FROM_HEADLINE = {"Mistral Large 2"}


def armA_prompt_pop(c, pop):
    return f"""You are executing a ClawBio pharmacogenomics skill. Determine the patient's diplotype from the genotype, then APPLY the skill rules below EXACTLY to output the phenotype and drug recommendation.

{h.skill_rules_text(c['gene'])}

## Patient
Cohort: {pop['name']} ({pop['desc']})
Gene: {c['gene']}
Genotype: {c['genotype']}
Drug: {c['drug']}
Population context: {c['pop_note'][pop['id']]}

## Output (4 lines only)
DIPLOTYPE: [called diplotype]
PHENOTYPE: [from Rule 1]
DRUG: [from Rule 2]
HAZARD: [clinical hazard]"""


def armB2_prompt_pop(c, pop):
    valid = "\n".join(f"  - {d}" for d in sorted(h.rules_dip[c['gene']].keys()))
    return f"""You are the input-interpretation step of a ClawBio pharmacogenomics agent. Map the patient's genotype to exactly ONE diplotype from the controlled list below, copying its text VERBATIM. Do not invent notation. A downstream validated skill computes phenotype and recommendation.

Patient cohort: {pop['name']} ({pop['desc']})
Gene: {c['gene']}
Valid diplotypes (choose one, copy verbatim):
{valid}

Patient genotype: {c['genotype']}
Population context: {c['pop_note'][pop['id']]}

## Output (1 line only)
DIPLOTYPE: [exact text of one list item]"""


_wlock = threading.Lock()


def run_one_fg(task):
    c, arm, prompt, model, fn, rep, pop = task
    text = ""
    for attempt in range(4):
        try:
            text = fn(prompt)
            if text and text.strip():
                break  # got a real response
            text = "<empty>"  # empty (e.g. token-budget truncation) -> retry
        except Exception as e:
            text = f"<error: {e}>"
        time.sleep(2 * (attempt + 1))
    gt_phen, gt_drug = c["gt_phenotype"], c["gt_drug"]
    if arm == "A_reasoning":
        pphen = h.parse_field(text, "PHENOTYPE"); pdrug = h.parse_field(text, "DRUG")
        comp = h.parse_field(text, "DIPLOTYPE")
    else:
        comp = h.parse_field(text, "DIPLOTYPE"); nd = h.norm_dip(comp, c["gene"])
        pphen = h.DIP2PHEN.get((c["gene"], nd), "UNKNOWN")
        pdrug = h.REC.get((c["gene"], c["drug"], nd), "UNKNOWN")
    r = {"tc": c["id"], "gene": c["gene"], "drug": c["drug"], "arm": arm,
         "model": model, "rep": rep, "pop": pop, "called_diplotype": comp,
         "pphen": pphen, "pdrug": pdrug[:80], "raw_text": (text or "")[:1000],
         "A1": h.rs.score_a1(pphen, gt_phen), "A2": h.rs.score_a2(pdrug, gt_drug),
         "A3": h.rs.score_a3(pdrug, gt_drug), "lethal": "lethal" in gt_drug.lower(),
         "matched": (pphen != "UNKNOWN") if arm == "B_execution" else None}
    with _wlock:
        with open(JSONL, "a") as f:
            f.write(json.dumps(r) + "\n")
    return r


def build_tasks(done_keys):
    tasks = []
    for pop in POPULATIONS:
        for c in h.cases:
            for arm, pf in [("A_reasoning", armA_prompt_pop), ("B_execution", armB2_prompt_pop)]:
                p = pf(c, pop)
                for model, fn in h.MODELS.items():
                    for rep in range(N_REPS):
                        key = [c["id"], arm, model, rep, pop["id"]]
                        if tuple(key) in done_keys:
                            continue
                        tasks.append((c, arm, p, model, fn, rep, pop["id"]))
    return tasks


def load_done():
    done, rows = set(), []
    if JSONL.exists():
        for line in JSONL.read_text().splitlines():
            if not line.strip():
                continue
            try:
                r = json.loads(line)
            except Exception:
                continue
            done.add((r["tc"], r["arm"], r["model"], r["rep"], r["pop"])); rows.append(r)
    return done, rows


def write_report(results):
    models = list(h.MODELS.keys())
    keep = [m for m in models if m not in EXCLUDE_FROM_HEADLINE]

    def mean(arm, dim, pop=None, model=None, ex_mistral=False):
        rows = [r for r in results if r["arm"] == arm
                and (pop is None or r["pop"] == pop)
                and (model is None or r["model"] == model)
                and (not ex_mistral or r["model"] not in EXCLUDE_FROM_HEADLINE)]
        return 100 * sum(r[dim] for r in rows) / len(rows) if rows else float("nan")

    def lethal(arm, pop=None, ex_mistral=True):
        return sum(1 for r in results if r["arm"] == arm and r["lethal"] and r["A3"] < 1
                   and (pop is None or r["pop"] == pop)
                   and (not ex_mistral or r["model"] not in EXCLUDE_FROM_HEADLINE))

    def repro(arm, ex_mistral=True):
        g = defaultdict(list)
        for r in results:
            if r["arm"] == arm and (not ex_mistral or r["model"] not in EXCLUDE_FROM_HEADLINE):
                g[(r["model"], r["tc"], r["pop"])].append(r["A1"])
        full = [v for v in g.values() if len(v) == N_REPS]
        return 100 * sum(1 for v in full if len(set(v)) == 1) / len(full) if full else float("nan")

    o = ["SKILL ARMS FULL GRID (Arm A skill-reasoning, Arm B skill-execution v2)",
         f"9 models x 3 pops x {N_REPS} reps x 110 cases x 2 arms; rows={len(results)}",
         "Population framing identical to 02-run-benchmark-v3.py.", ""]
    for label, exm in [("ALL 9 MODELS", False), ("EX-MISTRAL (8 models, headline)", True)]:
        o.append(f"=== AGGREGATE {label} ===")
        o.append(f"{'metric':26}{'Arm A':>11}{'Arm B':>11}")
        for d in ["A1", "A2", "A3"]:
            o.append(f"{'aggregate '+d:26}{mean('A_reasoning', d, ex_mistral=exm):10.1f}%{mean('B_execution', d, ex_mistral=exm):10.1f}%")
        o.append("")
    o.append(f"reproducibility 3/3 (ex-Mistral)  Arm A {repro('A_reasoning'):.1f}%  Arm B {repro('B_execution'):.1f}%")
    o.append("")
    o.append("=== PHENOTYPE A1 BY POPULATION (ex-Mistral) ===")
    o.append(f"{'arm':16}{'EUR':>9}{'AMR':>9}{'AFR':>9}{'spread':>9}")
    for arm, lab in [("A_reasoning", "Skill-reason"), ("B_execution", "Skill-exec")]:
        v = {p: mean(arm, 'A1', pop=p, ex_mistral=True) for p in ("EUR", "AMR", "AFR")}
        o.append(f"{lab:16}{v['EUR']:>9.1f}{v['AMR']:>9.1f}{v['AFR']:>9.1f}{max(v.values())-min(v.values()):>9.1f}")
    o.append("")
    o.append("=== LETHAL-CLASS A3 ERRORS BY POPULATION (ex-Mistral) ===")
    o.append(f"{'arm':16}{'EUR':>9}{'AMR':>9}{'AFR':>9}")
    for arm, lab in [("A_reasoning", "Skill-reason"), ("B_execution", "Skill-exec")]:
        o.append(f"{lab:16}" + "".join(f"{lethal(arm, p):>9d}" for p in ("EUR", "AMR", "AFR")))
    o.append("")
    o.append("=== PER MODEL (Arm A A1 / Arm B A1, pooled over pops+reps) ===")
    o.append(f"{'model':18}{'A_A1':>9}{'B_A1':>9}")
    for m in models:
        tag = "  (excl)" if m in EXCLUDE_FROM_HEADLINE else ""
        o.append(f"{m:18}{mean('A_reasoning','A1',model=m):8.1f}%{mean('B_execution','A1',model=m):8.1f}%{tag}")
    txt = "\n".join(o)
    REPORT.write_text(txt)
    print("\n" + txt)


def main():
    done, rows = load_done()
    tasks = build_tasks(done)
    target = len(POPULATIONS) * len(h.cases) * 2 * len(h.MODELS) * N_REPS
    print(f"target {target} calls | cached {len(done)} | to run {len(tasks)}", flush=True)
    results = list(rows)
    d = 0
    with ThreadPoolExecutor(max_workers=64) as ex:
        for r in ex.map(run_one_fg, tasks):
            results.append(r); d += 1
            if d % 200 == 0:
                print(f"  [{d}/{len(tasks)}] (+{len(done)} cached)", flush=True)
    dedup = {}
    for r in results:
        dedup[(r["tc"], r["arm"], r["model"], r["rep"], r["pop"])] = r
    final = list(dedup.values())
    OUT.write_text(json.dumps(final, indent=2))
    write_report(final)
    print("\nDONE ->", OUT.name, len(final), "rows", flush=True)


if __name__ == "__main__":
    if "--smoke" in sys.argv:
        tasks = build_tasks(set())
        print(f"models ({len(h.MODELS)}):", ", ".join(h.MODELS))
        print(f"populations: {[p['id'] for p in POPULATIONS]}  reps={N_REPS}")
        print(f"task count: {len(tasks)} (expected {3*110*2*len(h.MODELS)*N_REPS})")
        c0 = h.cases[0]
        print("\n--- sample Arm A prompt (AFR) ---\n" + armA_prompt_pop(c0, POPULATIONS[2])[:700])
        print("\n--- sample Arm B prompt (AMR) ---\n" + armB2_prompt_pop(c0, POPULATIONS[1])[:700])
    elif "--test1" in sys.argv:
        c0 = h.cases[0]
        for model, fn in h.MODELS.items():
            t0 = time.time()
            try:
                txt = fn(armB2_prompt_pop(c0, POPULATIONS[0]))
                ok = h.parse_field(txt, "DIPLOTYPE")
                print(f"  {model:18} OK {time.time()-t0:5.1f}s -> DIPLOTYPE: {ok!r}")
            except Exception as e:
                print(f"  {model:18} FAIL {time.time()-t0:5.1f}s -> {e}")
    else:
        main()
