#!/usr/bin/env python3
"""
Arm A (skill-as-reasoning) scaled to 9 models + Arm B v2 (skill-as-execution with
controlled-vocabulary input) across 9 models.

Arm A: model gets the validated decision rules and the genotype, reasons to phenotype+rec.
Arm B v2: model maps the genotype to ONE diplotype from the skill's controlled vocabulary
  (the list of valid diplotypes), copied verbatim; a validated skill then computes
  phenotype+recommendation in code. Controlled vocabulary removes the notation artifact
  seen in Arm B v1, isolating true input-interpretation accuracy.

9 models x 110 cases x EUR x 2 reps x 2 arms = 3,960 calls. Provider-level concurrency
caps + retry-backoff for maximum safe throughput.
"""
from __future__ import annotations
import json, os, re, time, threading, requests
from pathlib import Path
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from concurrent.futures import ThreadPoolExecutor
import openai, anthropic

BASE = Path(__file__).resolve().parent.parent
CASES = BASE / "specs" / "test_cases_v3.json"
OUT = BASE / "data" / "v3_armA9_armBv2.json"
REPORT = BASE / "data" / "v3_armA9_armBv2_report.txt"
ENV = BASE / ".env"   # repo-relative; environment variables take precedence (see .env.example)
N_REPS = 2

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
GKEY = _key("GEMINI_API_KEY", "GOOGLE_API_KEY"); MKEY = _key("MISTRAL_API_KEY")

# provider concurrency caps
sem = {"anthropic": threading.Semaphore(20), "openai": threading.Semaphore(20),
       "deepseek": threading.Semaphore(16), "gemini": threading.Semaphore(10),
       "mistral": threading.Semaphore(4)}

def _ant(model, p):
    with sem["anthropic"]:
        return ant.messages.create(model=model, max_tokens=320, messages=[{"role":"user","content":p}]).content[0].text
def _oai(model, p, reasoning=False):
    with sem["openai"]:
        if reasoning:
            return oai.chat.completions.create(model=model, max_completion_tokens=2000, messages=[{"role":"user","content":p}]).choices[0].message.content
        return oai.chat.completions.create(model=model, max_tokens=320, messages=[{"role":"user","content":p}]).choices[0].message.content
def _dsk(p):
    with sem["deepseek"]:
        return dsk.chat.completions.create(model="deepseek-chat", max_tokens=320, messages=[{"role":"user","content":p}]).choices[0].message.content
def _gem(p):
    with sem["gemini"]:
        r = requests.post(f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={GKEY}",
                          json={"contents":[{"parts":[{"text":p}]}], "generationConfig":{"maxOutputTokens":320}}, timeout=90)
        return r.json()["candidates"][0]["content"]["parts"][0]["text"]
# NOTE (peer review): Mistral is called through the moving 'mistral-large-latest'
# pointer, not a dated snapshot. See code/MODEL-VERSIONS.md for the reconciliation
# caveat; the headline is reported on the common eight-model set with Mistral
# excluded and reported separately, so no manuscript number depends on which
# snapshot this pointer resolved to.
def _mis(p):
    with sem["mistral"]:
        r = requests.post("https://api.mistral.ai/v1/chat/completions",
                          headers={"Authorization": f"Bearer {MKEY}", "Content-Type": "application/json"},
                          json={"model":"mistral-large-latest","messages":[{"role":"user","content":p}],"max_tokens":320}, timeout=90)
        return r.json()["choices"][0]["message"]["content"]

MODELS = {
    "Claude Opus 4":   lambda p: _ant("claude-opus-4-20250514", p),
    "Claude Sonnet 4": lambda p: _ant("claude-sonnet-4-20250514", p),
    "GPT-5.2":         lambda p: _oai("gpt-5.2", p, reasoning=True),
    "GPT-4.1":         lambda p: _oai("gpt-4.1", p),
    "o3":              lambda p: _oai("o3", p, reasoning=True),
    "o4-mini":         lambda p: _oai("o4-mini", p, reasoning=True),
    "Gemini 2.5 Flash":lambda p: _gem(p),
    "DeepSeek V3":     lambda p: _dsk(p),
    "Mistral Large 2": lambda p: _mis(p),
}

def norm_dip(s, gene):
    s = (s or "").lower().replace(gene.lower(), "")
    s = re.sub(r"\(.*?\)", "", s).replace("hla-", "")
    s = s.replace("positive", "pos").replace("negative", "neg").replace("carrier", "").replace("present", "pos").replace("absent", "neg")
    s = re.sub(r"[^\w*/:.\- ]", "", s)
    parts = [re.sub(r"\s+", "", p) for p in re.split(r"\s*/\s*", s) if p.strip()]
    if len(parts) == 2: parts = sorted(parts)
    return "/".join(parts) if parts else re.sub(r"\s+", "", s)

cases = json.loads(CASES.read_text()); cmap = {c["id"]: c for c in cases}
DIP2PHEN = {}; REC = {}; rules_dip = defaultdict(dict); rules_rec = defaultdict(list)
for c in cases:
    g = c["gene"]; nd = norm_dip(c["gt_diplotype"], g)
    DIP2PHEN[(g, nd)] = c["gt_phenotype"]; REC[(g, c["drug"], nd)] = c["gt_drug"]
    rules_dip[g][c["gt_diplotype"]] = c["gt_phenotype"]
    rules_rec[g].append((c["drug"], c["gt_diplotype"], c["gt_drug"]))

def skill_rules_text(gene):
    L = [f"# SKILL: {gene} (validated rules)", "## Rule 1 diplotype -> phenotype"]
    for dip, ph in sorted(rules_dip[gene].items()): L.append(f"  {dip} -> {ph}")
    L.append("## Rule 2 (drug, diplotype) -> recommendation"); seen=set()
    for drug, dip, rec in rules_rec[gene]:
        if (drug,dip) in seen: continue
        seen.add((drug,dip)); L.append(f"  ({drug}, {dip}) -> {rec}")
    return "\n".join(L)

def armA_prompt(c):
    return f"""You are executing a ClawBio pharmacogenomics skill. Determine the patient's diplotype from the genotype, then APPLY the skill rules below EXACTLY to output the phenotype and drug recommendation.

{skill_rules_text(c['gene'])}

## Patient
Gene: {c['gene']}
Genotype: {c['genotype']}
Drug: {c['drug']}

## Output (4 lines only)
DIPLOTYPE: [called diplotype]
PHENOTYPE: [from Rule 1]
DRUG: [from Rule 2]
HAZARD: [clinical hazard]"""

def armB2_prompt(c):
    valid = "\n".join(f"  - {d}" for d in sorted(rules_dip[c['gene']].keys()))
    return f"""You are the input-interpretation step of a ClawBio pharmacogenomics agent. Map the patient's genotype to exactly ONE diplotype from the controlled list below, copying its text VERBATIM. Do not invent notation. A downstream validated skill computes phenotype and recommendation.

Gene: {c['gene']}
Valid diplotypes (choose one, copy verbatim):
{valid}

Patient genotype: {c['genotype']}

## Output (1 line only)
DIPLOTYPE: [exact text of one list item]"""

def parse_field(text, field):
    if not text: return ""
    for line in text.split("\n"):
        if (field + ":") in line.upper():
            return line[line.upper().index(field + ":") + len(field) + 1:].strip()
    return ""

def run_one(task):
    c, arm, prompt, model, fn, rep = task
    text = ""
    for attempt in range(4):
        try: text = fn(prompt); break
        except Exception as e: text = f"<error: {e}>"; time.sleep(2 * (attempt + 1))
    gt_phen, gt_drug = c["gt_phenotype"], c["gt_drug"]
    if arm == "A_reasoning":
        pphen = parse_field(text, "PHENOTYPE"); pdrug = parse_field(text, "DRUG"); comp = parse_field(text, "DIPLOTYPE")
    else:
        comp = parse_field(text, "DIPLOTYPE"); nd = norm_dip(comp, c["gene"])
        pphen = DIP2PHEN.get((c["gene"], nd), "UNKNOWN"); pdrug = REC.get((c["gene"], c["drug"], nd), "UNKNOWN")
    return {"tc": c["id"], "gene": c["gene"], "drug": c["drug"], "arm": arm, "model": model, "rep": rep,
            "called_diplotype": comp, "pphen": pphen, "pdrug": pdrug[:80],
            "A1": rs.score_a1(pphen, gt_phen), "A2": rs.score_a2(pdrug, gt_drug), "A3": rs.score_a3(pdrug, gt_drug),
            "lethal": "lethal" in gt_drug.lower(), "matched": pphen != "UNKNOWN" if arm == "B_execution" else None}

def main():
    tasks = []
    for c in cases:
        for arm, pf in [("A_reasoning", armA_prompt), ("B_execution", armB2_prompt)]:
            p = pf(c)
            for model, fn in MODELS.items():
                for rep in range(N_REPS):
                    tasks.append((c, arm, p, model, fn, rep))
    total = len(tasks); results = []; done = 0
    print(f"running {total} calls", flush=True)
    with ThreadPoolExecutor(max_workers=72) as ex:
        for r in ex.map(run_one, tasks):
            results.append(r); done += 1
            if done % 200 == 0: print(f"  [{done}/{total}]", flush=True)
    OUT.write_text(json.dumps(results, indent=2))

    models = list(MODELS.keys())
    def mean(arm, dim, model=None):
        rows = [r for r in results if r["arm"] == arm and (model is None or r["model"] == model)]
        return 100 * sum(r[dim] for r in rows) / len(rows) if rows else float("nan")
    def repro(arm, model=None):
        g = defaultdict(list)
        for r in results:
            if r["arm"] == arm and (model is None or r["model"] == model): g[(r["model"], r["tc"])].append(r["A1"])
        full = [v for v in g.values() if len(v) == N_REPS]
        return 100 * sum(1 for v in full if len(set(v)) == 1) / len(full) if full else float("nan")
    def lethal(arm, model=None):
        return sum(1 for r in results if r["arm"] == arm and r["lethal"] and r["A3"] < 1 and (model is None or r["model"] == model))
    out = ["ARM A (skill-reasoning, 9 models) + ARM B v2 (skill-execution, controlled vocabulary)",
           f"9 models x 110 cases x EUR x {N_REPS} reps per arm", ""]
    out.append("=== AGGREGATE ===")
    out.append(f"{'metric':24}{'Arm A':>12}{'Arm B v2':>12}")
    for d in ["A1","A2","A3"]: out.append(f"{'aggregate '+d:24}{mean('A_reasoning',d):11.1f}%{mean('B_execution',d):11.1f}%")
    out.append(f"{'reproducibility(A1)':24}{repro('A_reasoning'):11.1f}%{repro('B_execution'):11.1f}%")
    out.append(f"{'lethal-class A3 errors':24}{lethal('A_reasoning'):12d}{lethal('B_execution'):12d}")
    bmatch = [r for r in results if r["arm"]=="B_execution"]
    out.append(f"Arm B v2 diplotype matched controlled vocab: {100*sum(1 for r in bmatch if r['matched'])/len(bmatch):.1f}%")
    out.append("")
    out.append("=== PER MODEL ===")
    out.append(f"{'model':18}{'A_A1':>8}{'A_repro':>9}{'B_A1':>8}{'B_repro':>9}{'A_leth':>8}{'B_leth':>8}")
    for m in models:
        out.append(f"{m:18}{mean('A_reasoning','A1',m):7.1f}%{repro('A_reasoning',m):8.1f}%{mean('B_execution','A1',m):7.1f}%{repro('B_execution',m):8.1f}%{lethal('A_reasoning',m):8d}{lethal('B_execution',m):8d}")
    text = "\n".join(out); REPORT.write_text(text); print("\n" + text)

if __name__ == "__main__":
    main()
