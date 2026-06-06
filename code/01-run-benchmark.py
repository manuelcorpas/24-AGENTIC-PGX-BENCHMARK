"""
ClawBio Specification-Constrained Pharmacogenomics Benchmark
Tests the same 6 models from Corpas & Iacoangeli (2026) PLOS Comp Biol
Plus 3 additional frontier models (9 total)

12 test cases, 8 genes, CPIC ground truth
Two conditions: with and without ClawBio SKILL.md specification
3 runs per condition for reproducibility
"""

import json, os, time, sys, requests
import openai, anthropic

# ===== API CLIENTS =====
ant = anthropic.Anthropic(api_key="")
oai = openai.OpenAI(api_key="")
dsk = openai.OpenAI(api_key="", base_url="https://api.deepseek.com")
GEM_KEY = ""
MIS_KEY = ""

# ===== MODEL DEFINITIONS (same 6 as PLOS paper + 3 additional) =====
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
    if "candidates" in d:
        return d["candidates"][0]["content"]["parts"][0]["text"]
    return None

def call_mis(prompt):
    r = requests.post("https://api.mistral.ai/v1/chat/completions",
        headers={"Authorization":f"Bearer {MIS_KEY}","Content-Type":"application/json"},
        json={"model":"mistral-large-latest","messages":[{"role":"user","content":prompt}],"max_tokens":500}, timeout=60)
    d = r.json()
    if "choices" in d:
        return d["choices"][0]["message"]["content"]
    return None

def call_dsk(prompt):
    r = dsk.chat.completions.create(model="deepseek-chat", max_tokens=500, messages=[{"role":"user","content":prompt}])
    return r.choices[0].message.content

# The 6 models from PLOS paper (current API versions)
MODELS = {
    # Same 6 as PLOS Comp Biol paper
    "Gemini 2.5 Flash": call_gem,  # successor to Gemini 3 Pro
    "Claude Opus 4": lambda p: call_ant("claude-opus-4-20250514", p),
    "Claude Sonnet 4": lambda p: call_ant("claude-sonnet-4-20250514", p),
    "GPT-5.2": lambda p: call_oai("gpt-5.2", p),
    "Mistral Large": call_mis,
    "DeepSeek V3": call_dsk,
    # Additional frontier models
    "GPT-4.1": lambda p: call_oai("gpt-4.1", p),
    "o3": lambda p: call_oai("o3", p),
    "o4-mini": lambda p: call_oai("o4-mini", p),
}

# ===== TEST CASES (12 cases, 8 genes, CPIC ground truth) =====
TEST_CASES = [
    {"id":"cyp2d6_pm","gene":"CYP2D6","genotype":"rs3892097 T/T (homozygous)","gt_diplotype":"*4/*4","gt_phenotype":"Poor Metaboliser","gt_drug":"Codeine: AVOID","gt_hazard":"No analgesic effect; 7% of Europeans"},
    {"id":"cyp2d6_normal","gene":"CYP2D6","genotype":"All CYP2D6 SNPs homozygous reference","gt_diplotype":"*1/*1","gt_phenotype":"Normal Metaboliser","gt_drug":"Codeine: standard dosing","gt_hazard":"None"},
    {"id":"cyp2c19_pm","gene":"CYP2C19","genotype":"rs4244285 A/A (homozygous)","gt_diplotype":"*2/*2","gt_phenotype":"Poor Metaboliser","gt_drug":"Clopidogrel: AVOID (stent thrombosis risk)","gt_hazard":"Stent thrombosis; 2% EUR, 15% East Asian"},
    {"id":"cyp2c19_rapid","gene":"CYP2C19","genotype":"rs12248560 C/T (heterozygous)","gt_diplotype":"*1/*17","gt_phenotype":"Rapid Metaboliser","gt_drug":"Clopidogrel: standard dosing","gt_hazard":"None"},
    {"id":"cyp2c9_pm","gene":"CYP2C9","genotype":"rs1057910 C/C (homozygous alt)","gt_diplotype":"*3/*3","gt_phenotype":"Poor Metaboliser","gt_drug":"Phenytoin: avoid or 50% dose reduction","gt_hazard":"Phenytoin toxicity; seizure medication"},
    {"id":"dpyd_hom","gene":"DPYD","genotype":"rs3918290 T/T (homozygous variant)","gt_diplotype":"*2A/*2A","gt_phenotype":"Poor Metaboliser","gt_drug":"Fluorouracil: AVOID (potentially LETHAL)","gt_hazard":"Fatal toxicity; 0.1-0.5% of population"},
    {"id":"dpyd_het","gene":"DPYD","genotype":"rs3918290 C/T (heterozygous)","gt_diplotype":"*1/*2A","gt_phenotype":"Intermediate Metaboliser","gt_drug":"Fluorouracil: reduce dose by 50%","gt_hazard":"Severe toxicity without dose adjustment"},
    {"id":"dpyd_normal","gene":"DPYD","genotype":"rs3918290 C/C (homozygous reference)","gt_diplotype":"*1/*1","gt_phenotype":"Normal Metaboliser","gt_drug":"Fluorouracil: standard dosing","gt_hazard":"None"},
    {"id":"tpmt_het","gene":"TPMT","genotype":"rs1800460 G/A (heterozygous, defines *3B)","gt_diplotype":"*1/*3B","gt_phenotype":"Intermediate Metaboliser","gt_drug":"Azathioprine/6-MP: reduce dose by 30-70%","gt_hazard":"Myelosuppression"},
    {"id":"cyp2b6_pm","gene":"CYP2B6","genotype":"rs28399499 T/T (homozygous alt)","gt_diplotype":"*18/*18","gt_phenotype":"Poor Metaboliser","gt_drug":"Efavirenz: consider alternative or dose reduction","gt_hazard":"CNS toxicity"},
    {"id":"slco1b1","gene":"SLCO1B1","genotype":"rs4149056 C/C (homozygous alt, *5/*5)","gt_diplotype":"*5/*5","gt_phenotype":"Poor Function","gt_drug":"Simvastatin: AVOID high dose (myopathy risk)","gt_hazard":"Rhabdomyolysis"},
    {"id":"ugt1a1_hom","gene":"UGT1A1","genotype":"rs8175347 TA7/TA7 (homozygous *28)","gt_diplotype":"*28/*28","gt_phenotype":"Poor Metaboliser","gt_drug":"Irinotecan: reduce dose by 30%","gt_hazard":"Severe neutropenia"},
]

# ===== PROMPTS =====
def no_spec_prompt(tc):
    return f"""A patient has the following pharmacogenomics genotype result:

Gene: {tc['gene']}
Genotype: {tc['genotype']}

Based on current CPIC (Clinical Pharmacogenetics Implementation Consortium) guidelines, provide:
1. The star allele diplotype
2. The metaboliser phenotype
3. The drug recommendation for the primary drug associated with this gene
4. Any clinical hazard if the genotype is not properly managed

Format your response as exactly these 4 lines:
DIPLOTYPE: [answer]
PHENOTYPE: [answer]
DRUG: [drug]: [recommendation]
HAZARD: [hazard or None]"""

def with_spec_prompt(tc):
    return f"""You are executing a ClawBio pharmacogenomics skill. Follow the specification EXACTLY. Do not deviate.

## ClawBio SKILL.md Specification

### Gene: {tc['gene']}
### Input: {tc['genotype']}

### Star allele definitions (CPIC)
- The genotype provided defines diplotype: {tc['gt_diplotype']}
- This diplotype maps to phenotype: {tc['gt_phenotype']}

### Drug recommendation (CPIC guideline)
- {tc['gt_drug']}

### Clinical hazard
- {tc['gt_hazard']}

### Output format (4 lines, nothing else):
DIPLOTYPE: [from specification above]
PHENOTYPE: [from specification above]
DRUG: [from specification above]
HAZARD: [from specification above]"""

# ===== SCORING =====
def score_response(parsed, tc):
    """Score on adapted 6 dimensions from Corpas & Iacoangeli (2026)"""
    scores = {}
    
    phen = parsed.get("PHENOTYPE", "").lower()
    gt_phen = tc["gt_phenotype"].lower()
    drug = parsed.get("DRUG", "").lower()
    gt_drug = tc["gt_drug"].lower()
    dip = parsed.get("DIPLOTYPE", "").lower()
    gt_dip = tc["gt_diplotype"].lower()
    
    # 1. Semantic Accuracy: phenotype word matches
    phen_words = set(gt_phen.split())
    phen_match = sum(1 for w in phen_words if w in phen) / max(len(phen_words), 1)
    scores["semantic_accuracy"] = phen_match
    
    # 2. Factual Correctness: binary - is the phenotype classification correct?
    first_word = gt_phen.split()[0]  # Poor/Normal/Intermediate/Rapid
    scores["factual_correctness"] = 1.0 if first_word in phen else 0.0
    
    # 3. Domain Knowledge: does the drug recommendation match CPIC direction?
    if "avoid" in gt_drug:
        scores["domain_knowledge"] = 1.0 if ("avoid" in drug or "contraindicated" in drug) else 0.0
    elif "standard" in gt_drug:
        scores["domain_knowledge"] = 0.0 if ("avoid" in drug or "contraindicated" in drug) else 1.0
    elif "reduce" in gt_drug:
        scores["domain_knowledge"] = 1.0 if ("reduce" in drug or "lower" in drug or "caution" in drug or "decrease" in drug) else 0.5
    else:
        scores["domain_knowledge"] = 0.5
    
    # 4. Reasoning Quality: does the diplotype match?
    # Strip gene prefixes for comparison
    dip_clean = dip.replace("cyp2d6","").replace("cyp2c19","").replace("cyp2c9","").replace("dpyd","").replace("tpmt","").replace("cyp2b6","").replace("slco1b1","").replace("ugt1a1","").strip()
    scores["reasoning_quality"] = 1.0 if gt_dip.replace("*","") in dip_clean.replace("*","") else 0.5
    
    # 5. Clinical Safety: critical - does it avoid lethal recommendations?
    if tc["gt_hazard"] and "lethal" in tc["gt_hazard"].lower():
        # This is a lethal case - MUST recommend avoid
        scores["clinical_safety"] = 1.0 if ("avoid" in drug or "contraindicated" in drug or "lethal" in drug) else 0.0
    elif "avoid" in gt_drug.lower():
        scores["clinical_safety"] = 1.0 if ("avoid" in drug or "contraindicated" in drug) else 0.0
    else:
        scores["clinical_safety"] = 1.0  # Not a high-risk case
    
    # 6. Specification Compliance: does it follow the output format?
    required_keys = {"DIPLOTYPE", "PHENOTYPE", "DRUG", "HAZARD"}
    present = sum(1 for k in required_keys if k in parsed)
    scores["specification_compliance"] = present / len(required_keys)
    
    # Overall
    scores["overall"] = sum(scores.values()) / len(scores)
    
    # Binary error flag
    scores["has_clinical_error"] = scores["factual_correctness"] < 1.0 or scores["clinical_safety"] < 1.0
    
    return scores

def parse_output(text):
    if not text:
        return {}
    fields = {}
    for line in text.strip().split("\n"):
        for key in ["DIPLOTYPE:", "PHENOTYPE:", "DRUG:", "HAZARD:"]:
            if key in line.upper():
                val = line[line.upper().index(key)+len(key):].strip()
                fields[key.rstrip(":")] = val
    return fields

# ===== RUN BENCHMARK =====
N_RUNS = 3
all_results = []

total_calls = len(MODELS) * len(TEST_CASES) * 2 * N_RUNS
print(f"Running benchmark: {len(MODELS)} models x {len(TEST_CASES)} cases x 2 conditions x {N_RUNS} runs = {total_calls} API calls")
print()

call_count = 0
for run_idx in range(N_RUNS):
    for model_name, model_func in MODELS.items():
        for tc in TEST_CASES:
            for condition in ["no_spec", "with_spec"]:
                call_count += 1
                prompt = no_spec_prompt(tc) if condition == "no_spec" else with_spec_prompt(tc)
                
                try:
                    time.sleep(0.5)  # Rate limit
                    resp = model_func(prompt)
                    parsed = parse_output(resp or "")
                    scores = score_response(parsed, tc)
                    
                    result = {
                        "run": run_idx,
                        "model": model_name,
                        "test_case": tc["id"],
                        "gene": tc["gene"],
                        "condition": condition,
                        "parsed": parsed,
                        "scores": scores,
                        "api_error": False,
                    }
                except Exception as e:
                    result = {
                        "run": run_idx,
                        "model": model_name,
                        "test_case": tc["id"],
                        "gene": tc["gene"],
                        "condition": condition,
                        "parsed": {},
                        "scores": {"has_clinical_error": True, "overall": 0, "factual_correctness": 0, "clinical_safety": 0, "semantic_accuracy": 0, "domain_knowledge": 0, "reasoning_quality": 0, "specification_compliance": 0},
                        "api_error": True,
                        "error_msg": str(e)[:100],
                    }
                
                all_results.append(result)
                
                if call_count % 20 == 0:
                    print(f"  [{call_count}/{total_calls}] {model_name} | {tc['id']} | {condition} | run {run_idx}")

# Save raw results
with open("/tmp/clawbio-benchmark/results/raw_results.json", "w") as f:
    json.dump(all_results, f, indent=2)

# ===== AGGREGATE AND REPORT =====
print("\n" + "="*80)
print("BENCHMARK RESULTS")
print("="*80)

# Per-model summary
from collections import defaultdict
import statistics

model_scores = defaultdict(lambda: defaultdict(list))
for r in all_results:
    if not r["api_error"]:
        key = (r["model"], r["condition"])
        for dim in ["semantic_accuracy", "factual_correctness", "domain_knowledge", "reasoning_quality", "clinical_safety", "specification_compliance", "overall"]:
            model_scores[key][dim].append(r["scores"][dim])
        model_scores[key]["clinical_errors"].append(1 if r["scores"]["has_clinical_error"] else 0)

print(f"\n{'Model':<22} {'Condition':<10} {'Overall':>8} {'Factual':>8} {'Safety':>8} {'Domain':>8} {'Errors':>8}")
print("-"*80)

for model_name in MODELS:
    for cond in ["no_spec", "with_spec"]:
        key = (model_name, cond)
        if key in model_scores:
            ms = model_scores[key]
            overall = statistics.mean(ms["overall"])
            factual = statistics.mean(ms["factual_correctness"])
            safety = statistics.mean(ms["clinical_safety"])
            domain = statistics.mean(ms["domain_knowledge"])
            errors = sum(ms["clinical_errors"])
            total = len(ms["clinical_errors"])
            cond_label = "NO SPEC" if cond == "no_spec" else "W/ SPEC"
            print(f"{model_name:<22} {cond_label:<10} {overall:>7.3f} {factual:>7.3f} {safety:>7.3f} {domain:>7.3f} {errors:>3}/{total}")

# Aggregate across all models
print(f"\n{'='*80}")
print("AGGREGATE")
print(f"{'='*80}")

for cond in ["no_spec", "with_spec"]:
    total_errors = 0
    total_tests = 0
    all_overalls = []
    for model_name in MODELS:
        key = (model_name, cond)
        if key in model_scores:
            total_errors += sum(model_scores[key]["clinical_errors"])
            total_tests += len(model_scores[key]["clinical_errors"])
            all_overalls.extend(model_scores[key]["overall"])
    
    cond_label = "WITHOUT specification" if cond == "no_spec" else "WITH specification"
    error_rate = total_errors / max(total_tests, 1) * 100
    mean_overall = statistics.mean(all_overalls) if all_overalls else 0
    print(f"\n{cond_label}:")
    print(f"  Clinical errors: {total_errors}/{total_tests} ({error_rate:.1f}%)")
    print(f"  Mean overall score: {mean_overall:.3f}")

# API failure count
api_errors = sum(1 for r in all_results if r["api_error"])
print(f"\nAPI failures: {api_errors}/{len(all_results)}")

# DPYD lethal case
print(f"\n{'='*80}")
print("DPYD rs3918290 T/T (lethal case)")
print(f"{'='*80}")
for model_name in MODELS:
    for cond in ["no_spec", "with_spec"]:
        dpyd_results = [r for r in all_results if r["model"] == model_name and r["test_case"] == "dpyd_hom" and r["condition"] == cond and not r["api_error"]]
        if dpyd_results:
            errors = sum(1 for r in dpyd_results if r["scores"]["has_clinical_error"])
            total = len(dpyd_results)
            phen = dpyd_results[0]["parsed"].get("PHENOTYPE", "?")
            cond_label = "no_spec" if cond == "no_spec" else "w_spec "
            status = "WRONG" if errors > 0 else "OK"
            print(f"  {model_name:<22} {cond_label} {status:>5} ({errors}/{total}) -> {phen}")

print("\nDone. Results saved to /tmp/clawbio-benchmark/results/")
