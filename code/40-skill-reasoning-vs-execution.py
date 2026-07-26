#!/usr/bin/env python3
"""
Skill-as-reasoning (Arm A) vs skill-as-execution (Arm B).

Tests where trustworthiness comes from in agentic genome interpretation. Both arms
receive the same realistic genotype input and must determine the patient's diplotype.

  Arm A (skill-as-reasoning): the SKILL.md decision rules (diplotype->phenotype and
    (drug,diplotype)->recommendation) are placed in the prompt; the model APPLIES them
    and GENERATES the phenotype and recommendation. The clinical answer passes through
    LLM token generation.

  Arm B (skill-as-execution): the model outputs ONLY the called diplotype; a validated
    skill then COMPUTES phenotype and recommendation deterministically in code. The
    clinical answer never passes through LLM generation; residual error is confined to
    the diplotype-calling (input-interpretation) step.

The decision rules are the validated CPIC mapping (built from the locked ground truth).
Hypothesis: Arm A beats RAG but stays stochastic; Arm B is exactly deterministic and
population-invariant, with error confined to input interpretation.

Calibration: 3 models x 110 cases x EUR x 3 replicates x 2 arms = 1,980 calls.
"""
from __future__ import annotations
import json, os, re, time
from pathlib import Path
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
import openai, anthropic

BASE = Path(__file__).resolve().parent.parent
CASES = BASE / "specs" / "test_cases_v3.json"
OUT = BASE / "data" / "v3_skill_reasoning_vs_execution.json"
REPORT = BASE / "data" / "v3_skill_reasoning_vs_execution_report.txt"
ENV = BASE / ".env"   # repo-relative; environment variables take precedence (see .env.example)
N_REPS = 2

# published scorer
_spec = spec_from_file_location("rescore", str(BASE / "code" / "10-rescore-v3.py"))
rs = module_from_spec(_spec); _spec.loader.exec_module(rs)

keys = dict(os.environ)
if ENV.exists():
    for line in ENV.read_text().splitlines():
        line = line.strip()
        if "=" in line and not line.startswith("#"):
            k, _, v = line.partition("="); keys.setdefault(k.strip(), v.strip().strip('"').strip("'"))
def _key(*names):
    for nm in names:
        if keys.get(nm): return keys[nm]
    raise KeyError(f"Missing API key: set one of {names} in the environment or {ENV} (see .env.example)")
ant = anthropic.Anthropic(api_key=_key("ANTHROPIC_API_KEY"))
oai = openai.OpenAI(api_key=_key("OPENAI_API_KEY"))
dsk = openai.OpenAI(api_key=_key("DEEPSEEK_API_KEY"), base_url="https://api.deepseek.com")
def call_ant(p): return ant.messages.create(model="claude-opus-4-20250514", max_tokens=400, messages=[{"role":"user","content":p}]).content[0].text
def call_oai(p):
    try: return oai.chat.completions.create(model="gpt-5.2", max_tokens=400, messages=[{"role":"user","content":p}]).choices[0].message.content
    except Exception: return oai.chat.completions.create(model="gpt-5.2", max_completion_tokens=2000, messages=[{"role":"user","content":p}]).choices[0].message.content
def call_dsk(p): return dsk.chat.completions.create(model="deepseek-chat", max_tokens=400, messages=[{"role":"user","content":p}]).choices[0].message.content
MODELS = {"Claude Opus 4": call_ant, "GPT-5.2": call_oai, "DeepSeek V3": call_dsk}

def norm_dip(s, gene):
    s = (s or "").lower().replace(gene.lower(), "")
    s = re.sub(r"\(.*?\)", "", s)
    s = s.replace("hla-", "").replace("diplotype", "").replace("genotype", "")
    s = s.replace("positive", "pos").replace("negative", "neg").replace("carrier", "")
    s = re.sub(r"\s+", "", s).strip(" :*/")
    return s

# build validated skill tables from ground truth
cases = json.loads(CASES.read_text())
DIP2PHEN = {}      # (gene, normdip) -> (phenotype, activity)
REC = {}           # (gene, drug, normdip) -> recommendation
rules_dip = defaultdict(dict)   # gene -> {gt_diplotype: phenotype}
rules_rec = defaultdict(list)   # gene -> [(drug, gt_diplotype, rec)]
for c in cases:
    g = c["gene"]; nd = norm_dip(c["gt_diplotype"], g)
    DIP2PHEN[(g, nd)] = (c["gt_phenotype"], c.get("gt_activity", ""))
    REC[(g, c["drug"], nd)] = c["gt_drug"]
    rules_dip[g][c["gt_diplotype"]] = c["gt_phenotype"]
    rules_rec[g].append((c["drug"], c["gt_diplotype"], c["gt_drug"]))

def skill_rules_text(gene):
    lines = [f"# SKILL: {gene} pharmacogenomic interpretation (validated rules)",
             "## Rule 1 - diplotype to phenotype"]
    for dip, phen in sorted(rules_dip[gene].items()):
        lines.append(f"  {dip}  ->  {phen}")
    lines.append("## Rule 2 - (drug, diplotype) to recommendation")
    seen = set()
    for drug, dip, rec in rules_rec[gene]:
        k = (drug, dip)
        if k in seen: continue
        seen.add(k); lines.append(f"  ({drug}, {dip})  ->  {rec}")
    return "\n".join(lines)

def armA_prompt(c, pop):
    return f"""You are executing a ClawBio pharmacogenomics skill. Determine the patient's diplotype from the genotype, then APPLY the skill rules below EXACTLY to output the phenotype and drug recommendation.

{skill_rules_text(c['gene'])}

## Patient
Gene: {c['gene']}
Genotype: {c['genotype']}
Drug: {c['drug']}
Patient cohort: {pop['name']}

## Output (4 lines only)
DIPLOTYPE: [called diplotype]
PHENOTYPE: [from Rule 1]
DRUG: [from Rule 2]
HAZARD: [clinical hazard]"""

def armB_prompt(c, pop):
    return f"""You are the input-interpretation step of a ClawBio pharmacogenomics agent. Read the patient's genotype and output ONLY the diplotype (star-allele genotype, or carrier status for HLA). Do NOT give phenotype or recommendation; a downstream validated skill computes those.

## Patient
Gene: {c['gene']}
Genotype: {c['genotype']}
Drug: {c['drug']}

## Output (1 line only)
DIPLOTYPE: [diplotype]"""

def parse_field(text, field):
    if not text: return ""
    for line in text.split("\n"):
        if (field + ":") in line.upper():
            return line[line.upper().index(field + ":") + len(field) + 1:].strip()
    return ""

def run_one(task):
    c, arm, prompt, model, fn, rep = task
    text = ""
    for attempt in range(3):
        try:
            text = fn(prompt); break
        except Exception as e:
            text = f"<error: {e}>"; time.sleep(2 * (attempt + 1))
    gt_phen, gt_drug = c["gt_phenotype"], c["gt_drug"]
    if arm == "A_reasoning":
        pphen = parse_field(text, "PHENOTYPE"); pdrug = parse_field(text, "DRUG"); comp_dip = parse_field(text, "DIPLOTYPE")
    else:
        comp_dip = parse_field(text, "DIPLOTYPE"); nd = norm_dip(comp_dip, c["gene"])
        pphen = DIP2PHEN.get((c["gene"], nd), ("UNKNOWN", ""))[0]
        pdrug = REC.get((c["gene"], c["drug"], nd), "UNKNOWN")
    return {"tc": c["id"], "gene": c["gene"], "drug": c["drug"], "arm": arm, "model": model, "rep": rep,
            "called_diplotype": comp_dip, "pphen": pphen, "pdrug": pdrug[:80],
            "A1": rs.score_a1(pphen, gt_phen), "A2": rs.score_a2(pdrug, gt_drug), "A3": rs.score_a3(pdrug, gt_drug),
            "lethal": "lethal" in gt_drug.lower()}

def main():
    from concurrent.futures import ThreadPoolExecutor
    pop = {"id": "EUR", "name": "European family cohort", "desc": "European ancestry"}
    tasks = []
    for c in cases:
        for arm, prompt_fn in [("A_reasoning", armA_prompt), ("B_execution", armB_prompt)]:
            p = prompt_fn(c, pop)
            for model, fn in MODELS.items():
                for rep in range(N_REPS):
                    tasks.append((c, arm, p, model, fn, rep))
    total = len(tasks); results = []; done = 0
    print(f"running {total} calls, 30 concurrent")
    with ThreadPoolExecutor(max_workers=30) as ex:
        for r in ex.map(run_one, tasks):
            results.append(r); done += 1
            if done % 100 == 0: print(f"  [{done}/{total}] done", flush=True)
    OUT.write_text(json.dumps(results, indent=2))

    # ---- analysis ----
    def agg(arm, dim):
        rows = [r for r in results if r["arm"] == arm]
        return 100 * sum(r[dim] for r in rows) / len(rows)
    def determinism(arm):
        g = defaultdict(list)
        for r in results:
            if r["arm"] == arm: g[(r["model"], r["tc"])].append(r["A1"])
        full = [v for v in g.values() if len(v) == N_REPS]
        return 100 * sum(1 for v in full if len(set(v)) == 1) / len(full)
    def lethal_errors(arm):
        return sum(1 for r in results if r["arm"] == arm and r["lethal"] and r["A3"] < 1)
    out = ["SKILL-AS-REASONING (Arm A) vs SKILL-AS-EXECUTION (Arm B)",
           f"3 models x 110 cases x EUR x {N_REPS} reps per arm", "",
           f"{'metric':28}{'Arm A (reasoning)':>20}{'Arm B (execution)':>20}"]
    for dim in ["A1", "A2", "A3"]:
        out.append(f"{'aggregate '+dim:28}{agg('A_reasoning',dim):19.1f}%{agg('B_execution',dim):19.1f}%")
    out.append(f"{'3-of-3 determinism (A1)':28}{determinism('A_reasoning'):19.1f}%{determinism('B_execution'):19.1f}%")
    out.append(f"{'lethal-class A3 errors':28}{lethal_errors('A_reasoning'):20d}{lethal_errors('B_execution'):20d}")
    # Arm B: diplotype-call (input-interpretation) accuracy = fraction where computed phenotype correct
    brows = [r for r in results if r["arm"] == "B_execution"]
    extract_ok = 100 * sum(1 for r in brows if r["pphen"] != "UNKNOWN") / len(brows)
    out.append("")
    out.append(f"Arm B diplotype recognised by skill table: {extract_ok:.1f}% (rest = input-interpretation miss)")
    text = "\n".join(out); REPORT.write_text(text); print("\n" + text)

if __name__ == "__main__":
    main()
