"""
ClawBio Specification-Constrained Pharmacogenomics Benchmark v2
Two-tier scoring: Clinical Correctness (A) + Contextual Quality (B)
12 test cases x 9 models x 3 populations x 2 conditions x 3 runs = 1,944 calls
"""

import json, os, time, requests
import openai, anthropic
from collections import defaultdict
import statistics

# ===== API CLIENTS =====
ant = anthropic.Anthropic(api_key="")
oai = openai.OpenAI(api_key="")
dsk = openai.OpenAI(api_key="", base_url="https://api.deepseek.com")
GEM_KEY = ""
MIS_KEY = ""

# ===== MODEL DEFINITIONS =====
def call_ant(model, prompt):
    r = ant.messages.create(model=model, max_tokens=500, messages=[{"role":"user","content":prompt}])
    return r.content[0].text

def call_oai(model, prompt):
    try:
        r = oai.chat.completions.create(model=model, max_tokens=500, messages=[{"role":"user","content":prompt}])
        return r.choices[0].message.content
    except:
        r = oai.chat.completions.create(model=model, max_completion_tokens=1000, messages=[{"role":"user","content":prompt}])
        return r.choices[0].message.content

def call_gem(prompt):
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={GEM_KEY}"
    resp = requests.post(url, json={"contents":[{"parts":[{"text":prompt}]}],"generationConfig":{"maxOutputTokens":2000}}, timeout=60)
    d = resp.json()
    return d["candidates"][0]["content"]["parts"][0]["text"] if "candidates" in d else None

def call_mis(prompt):
    r = requests.post("https://api.mistral.ai/v1/chat/completions",
        headers={"Authorization":f"Bearer {MIS_KEY}","Content-Type":"application/json"},
        json={"model":"mistral-large-latest","messages":[{"role":"user","content":prompt}],"max_tokens":500}, timeout=60)
    d = r.json()
    return d["choices"][0]["message"]["content"] if "choices" in d else None

def call_dsk(prompt):
    r = dsk.chat.completions.create(model="deepseek-chat", max_tokens=500, messages=[{"role":"user","content":prompt}])
    return r.choices[0].message.content

MODELS = {
    "Claude Opus 4": lambda p: call_ant("claude-opus-4-20250514", p),
    "Claude Sonnet 4": lambda p: call_ant("claude-sonnet-4-20250514", p),
    "GPT-5.2": lambda p: call_oai("gpt-5.2", p),
    "GPT-4.1": lambda p: call_oai("gpt-4.1", p),
    "o3": lambda p: call_oai("o3", p),
    "o4-mini": lambda p: call_oai("o4-mini", p),
    "Gemini 2.5 Flash": call_gem,
    "DeepSeek V3": call_dsk,
    "Mistral Large 2": call_mis,
}

# ===== TEST CASES =====
TEST_CASES = [
    {"id":"cyp2d6_pm","gene":"CYP2D6","genotype":"rs3892097 T/T (homozygous)","gt_diplotype":"*4/*4","gt_phenotype":"Poor Metaboliser","gt_drug":"Codeine: AVOID","gt_activity":"0.0",
     "pop_note":{"EUR":"CYP2D6*4 ~20% in Europeans; ~7% PM.","AMR":"CYP2D6*4 ~10% in admixed Latin Americans; *10 and *17 more common in indigenous.","AFR":"CYP2D6*4 ~6% in Africans; *17 (reduced function) ~20-35%, often missed by EUR panels."}},
    {"id":"cyp2d6_normal","gene":"CYP2D6","genotype":"All CYP2D6 SNPs homozygous reference","gt_diplotype":"*1/*1","gt_phenotype":"Normal Metaboliser","gt_drug":"Codeine: standard dosing","gt_activity":"2.0",
     "pop_note":{"EUR":"Standard genotype.","AMR":"Note: CYP2D6*10 common in indigenous populations, not detectable without specific testing.","AFR":"Note: CYP2D6*17 common (~20-35%), not detectable from reference-only panel."}},
    {"id":"cyp2c19_pm","gene":"CYP2C19","genotype":"rs4244285 A/A (homozygous)","gt_diplotype":"*2/*2","gt_phenotype":"Poor Metaboliser","gt_drug":"Clopidogrel: AVOID","gt_activity":"0.0",
     "pop_note":{"EUR":"CYP2C19*2 ~15% in EUR; PM ~2%.","AMR":"CYP2C19*2 ~10-15% in admixed Latin Americans.","AFR":"CYP2C19*2 ~15-18% in African populations; PM ~3-5%."}},
    {"id":"cyp2c19_rapid","gene":"CYP2C19","genotype":"rs12248560 C/T (heterozygous)","gt_diplotype":"*1/*17","gt_phenotype":"Rapid Metaboliser","gt_drug":"Clopidogrel: standard dosing","gt_activity":"1.5",
     "pop_note":{"EUR":"CYP2C19*17 ~21% in EUR.","AMR":"CYP2C19*17 ~15% in admixed Latin Americans.","AFR":"CYP2C19*17 ~16-26% in African populations."}},
    {"id":"cyp2c9_pm","gene":"CYP2C9","genotype":"rs1057910 C/C (homozygous alt)","gt_diplotype":"*3/*3","gt_phenotype":"Poor Metaboliser","gt_drug":"Phenytoin: avoid or reduce dose by 50%","gt_activity":"0.0",
     "pop_note":{"EUR":"CYP2C9*3 ~7% in EUR. PM ~0.4%.","AMR":"CYP2C9*3 ~3-5% in admixed Latin Americans.","AFR":"CYP2C9*3 rare (<1%) in African populations; *5, *6, *8, *11 more relevant."}},
    {"id":"dpyd_hom","gene":"DPYD","genotype":"rs3918290 T/T (homozygous variant)","gt_diplotype":"*2A/*2A","gt_phenotype":"Poor Metaboliser","gt_drug":"Fluorouracil: AVOID (potentially LETHAL)","gt_activity":"0.0",
     "pop_note":{"EUR":"DPYD*2A ~1% in EUR. Pre-treatment testing mandated by EMA (2020).","AMR":"DPYD*2A poorly characterised in Latin American populations. Testing infrastructure limited.","AFR":"DPYD variant frequencies poorly studied in African populations. Risk of undetected DPD deficiency."}},
    {"id":"dpyd_het","gene":"DPYD","genotype":"rs3918290 C/T (heterozygous)","gt_diplotype":"*1/*2A","gt_phenotype":"Intermediate Metaboliser","gt_drug":"Fluorouracil: reduce dose by 50%","gt_activity":"1.0",
     "pop_note":{"EUR":"Heterozygous carriers ~2% in EUR.","AMR":"Carrier frequency unknown in most Latin American populations.","AFR":"Carrier frequency unknown in most African populations."}},
    {"id":"dpyd_normal","gene":"DPYD","genotype":"rs3918290 C/C (homozygous reference)","gt_diplotype":"*1/*1","gt_phenotype":"Normal Metaboliser","gt_drug":"Fluorouracil: standard dosing","gt_activity":"2.0",
     "pop_note":{"EUR":"Standard genotype.","AMR":"Standard genotype, but other DPYD variants may be relevant.","AFR":"Standard at this position, but population-specific DPYD variants are understudied."}},
    {"id":"tpmt_het","gene":"TPMT","genotype":"rs1800460 G/A (heterozygous, defines *3B)","gt_diplotype":"*1/*3B","gt_phenotype":"Intermediate Metaboliser","gt_drug":"Azathioprine/6-MP: reduce dose by 30-70%","gt_activity":"1.0",
     "pop_note":{"EUR":"TPMT*3B ~1% in EUR. Heterozygous IM ~10%.","AMR":"TPMT variant frequencies similar to EUR in admixed populations.","AFR":"TPMT*3C more common than *3B in African populations (~5-8%)."}},
    {"id":"cyp2b6_pm","gene":"CYP2B6","genotype":"rs28399499 T/T (homozygous alt)","gt_diplotype":"*18/*18","gt_phenotype":"Poor Metaboliser","gt_drug":"Efavirenz: consider alternative or dose reduction","gt_activity":"0.0",
     "pop_note":{"EUR":"CYP2B6*18 rare in EUR (<1%). Efavirenz replaced by integrase inhibitors.","AMR":"CYP2B6*18 intermediate frequency. Efavirenz still used in some HIV programmes.","AFR":"CYP2B6*18 ~4-8% in East Africa. Efavirenz remains first-line in many African HIV programmes. HIGH IMPACT."}},
    {"id":"slco1b1","gene":"SLCO1B1","genotype":"rs4149056 C/C (homozygous alt, *5/*5)","gt_diplotype":"*5/*5","gt_phenotype":"Poor Function","gt_drug":"Simvastatin: avoid high dose (myopathy risk)","gt_activity":"0.0",
     "pop_note":{"EUR":"SLCO1B1*5 ~15% in EUR. Homozygous ~2%.","AMR":"SLCO1B1*5 ~10-15% in Latin Americans.","AFR":"SLCO1B1*5 ~2-3% in African populations. Lower statin myopathy risk at population level."}},
    {"id":"ugt1a1_hom","gene":"UGT1A1","genotype":"rs8175347 TA7/TA7 (homozygous *28)","gt_diplotype":"*28/*28","gt_phenotype":"Poor Metaboliser","gt_drug":"Irinotecan: reduce dose by 30%","gt_activity":"0.0",
     "pop_note":{"EUR":"UGT1A1*28 ~35% in EUR. Homozygous ~10%.","AMR":"UGT1A1*28 ~30-40% in admixed Latin Americans.","AFR":"UGT1A1*28 ~40-56% in African populations. Higher homozygous frequency."}},
]

POPULATIONS = [
    {"id":"EUR","name":"European family cohort (Corpasome project)","desc":"European ancestry, whole-genome sequencing"},
    {"id":"AMR","name":"Peruvian Genome Project","desc":"Admixed Latin American, 7 indigenous and mestizo subpopulations"},
    {"id":"AFR","name":"Uganda Genome Resource","desc":"East African, 6,407 whole-genome sequences"},
]

# ===== PROMPTS =====
def no_spec_prompt(tc, pop):
    return f"""A patient from the {pop['name']} ({pop['desc']}) has the following pharmacogenomics genotype:

Gene: {tc['gene']}
Genotype: {tc['genotype']}

Based on CPIC guidelines, provide:
1. Star allele diplotype
2. Metaboliser phenotype
3. Drug recommendation
4. Clinical hazard if mismanaged
5. Population-specific considerations

Format as exactly 5 lines:
DIPLOTYPE: [answer]
PHENOTYPE: [answer]
DRUG: [drug]: [recommendation]
HAZARD: [hazard or None]
POPULATION: [population-specific note]"""

def with_spec_prompt(tc, pop):
    return f"""You are executing a ClawBio pharmacogenomics skill. Follow the specification EXACTLY.

## SKILL.md Specification
Gene: {tc['gene']}
Input: {tc['genotype']}
Patient cohort: {pop['name']} ({pop['desc']})
Diplotype: {tc['gt_diplotype']}
Phenotype: {tc['gt_phenotype']}
Activity score: {tc['gt_activity']}
Drug recommendation: {tc['gt_drug']}
Population context: {tc['pop_note'][pop['id']]}

## Output (5 lines only):
DIPLOTYPE: {tc['gt_diplotype']}
PHENOTYPE: {tc['gt_phenotype']}
DRUG: {tc['gt_drug']}
HAZARD: [clinical hazard]
POPULATION: {tc['pop_note'][pop['id']]}"""

# ===== PARSING =====
def parse(text):
    if not text:
        return {}
    fields = {}
    for line in text.strip().split("\n"):
        for key in ["DIPLOTYPE:","PHENOTYPE:","DRUG:","HAZARD:","POPULATION:"]:
            if key in line.upper():
                val = line[line.upper().index(key)+len(key):].strip()
                fields[key.rstrip(":")] = val
    return fields

# ===== TWO-TIER SCORING =====
def score(parsed, tc, pop):
    if not parsed:
        return {"A1":0,"A2":0,"A3":0,"B1":0,"B2":0,"B3":0,"tier_a":0,"tier_b":0,"overall":0,"format_fail":True}

    phen = parsed.get("PHENOTYPE","").lower()
    drug = parsed.get("DRUG","").lower()
    pop_text = parsed.get("POPULATION","").lower()
    gt_first = tc["gt_phenotype"].lower().split()[0]

    s = {}

    # TIER A: Clinical Correctness
    # A1: Phenotype correct?
    s["A1"] = 1.0 if gt_first in phen else 0.0

    # A2: Drug direction correct?
    gt_drug = tc["gt_drug"].lower()
    if "avoid" in gt_drug and "lethal" in gt_drug:
        s["A2"] = 1.0 if ("avoid" in drug or "contraindicated" in drug) else 0.0
    elif "avoid" in gt_drug:
        s["A2"] = 1.0 if ("avoid" in drug or "contraindicated" in drug or "reduce" in drug or "alternative" in drug) else 0.0
    elif "standard" in gt_drug:
        s["A2"] = 0.0 if ("avoid" in drug or "contraindicated" in drug) else 1.0
    elif "reduce" in gt_drug:
        s["A2"] = 1.0 if ("reduce" in drug or "lower" in drug or "caution" in drug or "decrease" in drug or "adjust" in drug) else 0.5
    else:
        s["A2"] = 0.5

    # A3: Clinical safety (lethal cases only)
    if "lethal" in tc.get("gt_drug","").lower():
        s["A3"] = 1.0 if ("avoid" in drug or "contraindicated" in drug) else 0.0
    else:
        s["A3"] = 1.0

    s["tier_a"] = (s["A1"] + s["A2"] + s["A3"]) / 3.0

    # TIER B: Contextual Quality
    # B1: Population awareness
    pop_kw = {"EUR":["european","caucasian"],"AMR":["latin","peru","admixed","indigenous","amerindian"],"AFR":["african","uganda","east african","hiv","efavirenz","first-line"]}
    pop_hits = sum(1 for kw in pop_kw.get(pop["id"],[]) if kw in pop_text)
    s["B1"] = min(pop_hits / 2.0, 1.0)

    # B2: Reasoning chain completeness
    has = sum(1 for k in ["DIPLOTYPE","PHENOTYPE","DRUG","HAZARD","POPULATION"] if k in parsed and len(parsed[k]) > 1)
    s["B2"] = has / 5.0

    # B3: Domain grounding
    full = " ".join(parsed.values()).lower()
    domain_hits = sum(1 for kw in ["cpic","guideline","activity","score","allele","frequency"] if kw in full)
    s["B3"] = min(domain_hits / 2.0, 1.0)

    s["tier_b"] = (s["B1"] + s["B2"] + s["B3"]) / 3.0
    s["overall"] = (s["tier_a"] + s["tier_b"]) / 2.0
    s["format_fail"] = False

    return s

# ===== RUN =====
N_RUNS = 3
total = len(MODELS) * len(TEST_CASES) * len(POPULATIONS) * 2 * N_RUNS
print(f"Benchmark: {len(MODELS)} models x {len(TEST_CASES)} cases x {len(POPULATIONS)} pops x 2 conds x {N_RUNS} runs = {total} calls")

results = []
count = 0

for run in range(N_RUNS):
    for mname, mfunc in MODELS.items():
        for tc in TEST_CASES:
            for pop in POPULATIONS:
                for cond in ["no_spec", "with_spec"]:
                    count += 1
                    prompt = no_spec_prompt(tc, pop) if cond == "no_spec" else with_spec_prompt(tc, pop)
                    try:
                        time.sleep(0.3)
                        resp = mfunc(prompt)
                        parsed = parse(resp)
                        scores = score(parsed, tc, pop)
                        results.append({"run":run,"model":mname,"tc":tc["id"],"gene":tc["gene"],"pop":pop["id"],"cond":cond,"parsed":parsed,"scores":scores})
                    except Exception as e:
                        results.append({"run":run,"model":mname,"tc":tc["id"],"gene":tc["gene"],"pop":pop["id"],"cond":cond,"parsed":{},"scores":{"A1":0,"A2":0,"A3":0,"B1":0,"B2":0,"B3":0,"tier_a":0,"tier_b":0,"overall":0,"format_fail":True},"error":str(e)[:80]})

                    if count % 50 == 0:
                        print(f"  [{count}/{total}] {mname} | {tc['id']} | {pop['id']} | {cond} | run {run}")
                        # Save checkpoint + partial results
                        for cp_path in ["/tmp/clawbio-benchmark/results/v2_checkpoint.json",
                                        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "RESULTS", "v2_checkpoint.json")]:
                            os.makedirs(os.path.dirname(cp_path), exist_ok=True)
                            with open(cp_path, "w") as f:
                                json.dump({"count":count,"total":total,"completed":count}, f)
                    if count % 200 == 0:
                        # Save partial results every 200 calls
                        partial_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "RESULTS", "v2_partial.json")
                        with open(partial_path, "w") as f:
                            json.dump(results, f)

# Save full results to both locations
for path in ["/tmp/clawbio-benchmark/results/v2_raw.json",
             os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "RESULTS", "v2_raw.json")]:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(results, f, indent=2)

# ===== REPORT =====
print(f"\n{'='*90}")
print(f"BENCHMARK RESULTS: {count} calls completed")
print(f"{'='*90}")

# Format failures
fmt_fails = sum(1 for r in results if r["scores"].get("format_fail"))
api_errors = sum(1 for r in results if "error" in r)
valid = [r for r in results if not r["scores"].get("format_fail")]
print(f"API errors: {api_errors}, Format failures: {fmt_fails}, Valid: {len(valid)}")

# TIER A: Clinical Correctness by model and condition
print(f"\n--- TIER A: CLINICAL CORRECTNESS ---")
print(f"{'Model':<22} {'No spec':>20} {'With spec':>20}")
print(f"{'':22} {'A1':>5} {'A2':>5} {'A3':>5} {'TierA':>5}  {'A1':>5} {'A2':>5} {'A3':>5} {'TierA':>5}")
print("-"*70)

for mname in MODELS:
    for cond in ["no_spec","with_spec"]:
        mv = [r for r in valid if r["model"]==mname and r["cond"]==cond]
        if mv:
            a1 = statistics.mean([r["scores"]["A1"] for r in mv])
            a2 = statistics.mean([r["scores"]["A2"] for r in mv])
            a3 = statistics.mean([r["scores"]["A3"] for r in mv])
            ta = statistics.mean([r["scores"]["tier_a"] for r in mv])
            if cond == "no_spec":
                print(f"{mname:<22} {a1:>5.2f} {a2:>5.2f} {a3:>5.2f} {ta:>5.2f}", end="")
            else:
                print(f"  {a1:>5.2f} {a2:>5.2f} {a3:>5.2f} {ta:>5.2f}")

# TIER B: Contextual Quality
print(f"\n--- TIER B: CONTEXTUAL QUALITY ---")
print(f"{'Model':<22} {'No spec':>20} {'With spec':>20}")
print(f"{'':22} {'B1':>5} {'B2':>5} {'B3':>5} {'TierB':>5}  {'B1':>5} {'B2':>5} {'B3':>5} {'TierB':>5}")
print("-"*70)

for mname in MODELS:
    for cond in ["no_spec","with_spec"]:
        mv = [r for r in valid if r["model"]==mname and r["cond"]==cond]
        if mv:
            b1 = statistics.mean([r["scores"]["B1"] for r in mv])
            b2 = statistics.mean([r["scores"]["B2"] for r in mv])
            b3 = statistics.mean([r["scores"]["B3"] for r in mv])
            tb = statistics.mean([r["scores"]["tier_b"] for r in mv])
            if cond == "no_spec":
                print(f"{mname:<22} {b1:>5.2f} {b2:>5.2f} {b3:>5.2f} {tb:>5.2f}", end="")
            else:
                print(f"  {b1:>5.2f} {b2:>5.2f} {b3:>5.2f} {tb:>5.2f}")

# CONSISTENCY (perfect 3/3 rate)
print(f"\n--- CONSISTENCY (fraction of model-case-pop combos with 3/3 correct) ---")
for cond in ["no_spec","with_spec"]:
    perfect = 0
    total_combos = 0
    for mname in MODELS:
        for tc in TEST_CASES:
            for pop in POPULATIONS:
                runs = [r for r in valid if r["model"]==mname and r["tc"]==tc["id"] and r["pop"]==pop["id"] and r["cond"]==cond]
                if len(runs) >= 2:
                    correct = sum(1 for r in runs if r["scores"]["A1"]==1.0 and r["scores"]["A3"]==1.0)
                    total_combos += 1
                    if correct == len(runs):
                        perfect += 1
    label = "WITHOUT spec" if cond == "no_spec" else "WITH spec   "
    print(f"  {label}: {perfect}/{total_combos} perfect ({perfect/max(total_combos,1)*100:.1f}%)")

# BY POPULATION
print(f"\n--- TIER A BY POPULATION (no spec only) ---")
for pop in POPULATIONS:
    pv = [r for r in valid if r["pop"]==pop["id"] and r["cond"]=="no_spec"]
    if pv:
        a1 = statistics.mean([r["scores"]["A1"] for r in pv])
        a3 = statistics.mean([r["scores"]["A3"] for r in pv])
        print(f"  {pop['id']}: Phenotype accuracy={a1:.2f}, Clinical safety={a3:.2f}, n={len(pv)}")

# DPYD lethal case by population
print(f"\n--- DPYD LETHAL CASE BY POPULATION ---")
for pop in POPULATIONS:
    for cond in ["no_spec","with_spec"]:
        dpyd = [r for r in valid if r["tc"]=="dpyd_hom" and r["pop"]==pop["id"] and r["cond"]==cond]
        if dpyd:
            wrong = sum(1 for r in dpyd if r["scores"]["A1"]<1.0 or r["scores"]["A3"]<1.0)
            label = "no_spec" if cond=="no_spec" else "w_spec "
            print(f"  {pop['id']} {label}: {wrong}/{len(dpyd)} errors")

print(f"\nResults saved to /tmp/clawbio-benchmark/results/v2_raw.json")
