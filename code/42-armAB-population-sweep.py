#!/usr/bin/env python3
"""
Population sweep for the skill arms (Arm A reasoning, Arm B execution v2).

The main skill run (v3_armA9_armBv2.json) carried no population context, so the
population-invariance guarantee in F01 was only measured for the free/RAG/control
conditions. This script runs BOTH skill arms across all three populations with the
SAME population framing used by the three-arm benchmark (02-run-benchmark-v3.py),
so population-invariance of the skill arms is measured rather than asserted.

Design:
  - Populations: EUR (Corpasome), AMR (Peruvian Genome Project), AFR (Uganda Genome
    Resource, East African) -- identical framing to the three-arm experiment.
  - Clean 8 models (Mistral Large 2 excluded; handled by its own paced rerun).
  - 1 rep per cell (the EUR main run already established reproducibility; this sweep
    measures accuracy invariance across populations, a point estimate per population).
  - Writes data/v3_armA9_armBv2_POP.json + _POP_report.txt. Does NOT touch the
    existing EUR run file.

Call count: 110 cases x 2 arms x 8 models x 3 pops x 1 rep = 5,280 calls.

Run:   python3 code/42-armAB-population-sweep.py
Smoke: python3 code/42-armAB-population-sweep.py --smoke   (build tasks, no API calls)
"""
from __future__ import annotations
import sys
import json
import time
from pathlib import Path
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from concurrent.futures import ThreadPoolExecutor

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "data" / "v3_armA9_armBv2_POP.json"
REPORT = BASE / "data" / "v3_armA9_armBv2_POP_report.txt"

# reuse the validated harness (run_one, prompts, skill rules, model adapters)
_sp = spec_from_file_location("h41", str(BASE / "code" / "41-armA9-armBv2.py"))
h = module_from_spec(_sp)

POPULATIONS = [
    {"id": "EUR", "name": "European family cohort (Corpasome project)",
     "desc": "European ancestry, whole-genome sequencing"},
    {"id": "AMR", "name": "Peruvian Genome Project",
     "desc": "Admixed Latin American, 7 indigenous and mestizo subpopulations"},
    {"id": "AFR", "name": "Uganda Genome Resource",
     "desc": "East African, 6,407 whole-genome sequences"},
]

EXCLUDE = {"Mistral Large 2"}
N_REPS = 1


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


def build_tasks():
    models = {m: fn for m, fn in h.MODELS.items() if m not in EXCLUDE}
    tasks = []
    for pop in POPULATIONS:
        for c in h.cases:
            for arm, pf in [("A_reasoning", armA_prompt_pop), ("B_execution", armB2_prompt_pop)]:
                p = pf(c, pop)
                for model, fn in models.items():
                    for rep in range(N_REPS):
                        tasks.append((c, arm, p, model, fn, rep, pop["id"]))
    return tasks, list(models.keys())


def run_one_pop(task):
    c, arm, prompt, model, fn, rep, pop = task
    r = h.run_one((c, arm, prompt, model, fn, rep))
    r["pop"] = pop
    return r


def write_report(results, models):
    def acc(arm, pop=None, model=None):
        rows = [r for r in results if r["arm"] == arm
                and (pop is None or r["pop"] == pop)
                and (model is None or r["model"] == model)]
        return 100 * sum(r["A1"] for r in rows) / len(rows) if rows else float("nan")

    def leth(arm, pop):
        rows = [r for r in results if r["arm"] == arm and r["pop"] == pop and r["lethal"]]
        return sum(1 for r in rows if r["A3"] < 1)

    out = ["SKILL-ARM POPULATION SWEEP (Arm A reasoning, Arm B execution v2)",
           f"{len(models)} models x 110 cases x 3 pops x {N_REPS} rep per arm "
           f"(clean set; Mistral excluded)",
           "Population framing identical to 02-run-benchmark-v3.py.", ""]
    out.append("=== PHENOTYPE ACCURACY (A1) BY POPULATION ===")
    out.append(f"{'arm':16}{'EUR':>10}{'AMR':>10}{'AFR':>10}{'spread':>10}")
    for arm, label in [("A_reasoning", "Skill-reasoning"), ("B_execution", "Skill-execution")]:
        vals = {p: acc(arm, p) for p in ("EUR", "AMR", "AFR")}
        sp = max(vals.values()) - min(vals.values())
        out.append(f"{label:16}{vals['EUR']:>10.1f}{vals['AMR']:>10.1f}{vals['AFR']:>10.1f}{sp:>10.1f}")
    out.append("")
    out.append("=== LETHAL-CLASS SAFETY ERRORS BY POPULATION ===")
    out.append(f"{'arm':16}{'EUR':>10}{'AMR':>10}{'AFR':>10}")
    for arm, label in [("A_reasoning", "Skill-reasoning"), ("B_execution", "Skill-execution")]:
        out.append(f"{label:16}" + "".join(f"{leth(arm, p):>10d}" for p in ("EUR", "AMR", "AFR")))
    out.append("")
    out.append("=== PER MODEL x POPULATION (Skill-execution A1) ===")
    out.append(f"{'model':18}{'EUR':>10}{'AMR':>10}{'AFR':>10}{'spread':>10}")
    for m in models:
        vals = {p: acc("B_execution", p, m) for p in ("EUR", "AMR", "AFR")}
        sp = max(vals.values()) - min(vals.values())
        out.append(f"{m:18}{vals['EUR']:>10.1f}{vals['AMR']:>10.1f}{vals['AFR']:>10.1f}{sp:>10.1f}")
    text = "\n".join(out)
    REPORT.write_text(text)
    print("\n" + text)


def main():
    _sp.loader.exec_module(h)  # load harness (reads .env, builds clients + rules)
    tasks, models = build_tasks()
    total = len(tasks)
    print(f"population sweep: {total} calls "
          f"({len(models)} models x 3 pops x 2 arms x 110 cases x {N_REPS} rep)", flush=True)
    results = []
    done = 0
    with ThreadPoolExecutor(max_workers=64) as ex:
        for r in ex.map(run_one_pop, tasks):
            results.append(r)
            done += 1
            if done % 200 == 0:
                print(f"  [{done}/{total}]", flush=True)
    OUT.write_text(json.dumps(results, indent=2))
    write_report(results, models)
    print("\nDONE. ->", OUT.name)


if __name__ == "__main__":
    if "--smoke" in sys.argv:
        _sp.loader.exec_module(h)
        tasks, models = build_tasks()
        print(f"models ({len(models)}):", ", ".join(models))
        print(f"task count: {len(tasks)} (expected {110*2*len(models)*3*N_REPS})")
        c0 = h.cases[0]
        print("\n--- sample Arm A prompt (AFR) ---")
        print(armA_prompt_pop(c0, POPULATIONS[2])[:600])
        print("\n--- sample Arm B prompt (AMR) ---")
        print(armB2_prompt_pop(c0, POPULATIONS[1])[:600])
    else:
        main()
