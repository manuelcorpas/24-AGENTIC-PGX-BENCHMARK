#!/usr/bin/env python3
"""
Step 4 — run the model panel (the agent's interpretation step) on the REAL
diplotypes observed in a cohort, predicting the CPIC phenotype for each.

Input: the aggregated diplotype table from step 3 (or any TSV with columns
       cohort, gene, diplotype, phenotype). Unique (gene, diplotype) pairs are
       queried once per model.
Output: predictions TSV (cohort, gene, diplotype, cohort_phenotype, model, pred).

The skill prompt instructs the model to apply CPIC and to output 'Indeterminate'
when the diplotype contains an allele of uncertain/unknown function (so that
failure-to-abstain is measurable).

API keys are read from the environment (ANTHROPIC_API_KEY, OPENAI_API_KEY,
GOOGLE_API_KEY, DEEPSEEK_API_KEY). Models and provider tags come from
config/models.txt. No keys are stored in the repo.

Usage: 04_run_agent_realgenome.py <diplotypes.tsv> <out_predictions.tsv>
"""
import os
import sys
import csv
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests
import openai
import anthropic

CFG = Path(__file__).resolve().parent.parent / "config" / "models.txt"

def load_models():
    models = []
    for line in CFG.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # "Display Name   model_id   provider-tag"  (split on 2+ spaces)
        import re
        parts = re.split(r"\s{2,}", line)
        if len(parts) >= 3:
            models.append((parts[0], parts[1], parts[2]))
    return models

def make_caller(model_id, provider):
    ant = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"]) if provider == "anthropic" else None
    oai = openai.OpenAI(api_key=os.environ["OPENAI_API_KEY"]) if provider.startswith("openai") else None
    dsk = openai.OpenAI(api_key=os.environ["DEEPSEEK_API_KEY"], base_url="https://api.deepseek.com") if provider == "deepseek" else None
    gkey = os.environ.get("GOOGLE_API_KEY", "")
    def call(p):
        if provider == "anthropic":
            return ant.messages.create(model=model_id, max_tokens=120, messages=[{"role": "user", "content": p}]).content[0].text
        if provider == "openai-reasoning":
            return oai.chat.completions.create(model=model_id, max_completion_tokens=2000, messages=[{"role": "user", "content": p}]).choices[0].message.content
        if provider == "openai":
            return oai.chat.completions.create(model=model_id, max_tokens=120, messages=[{"role": "user", "content": p}]).choices[0].message.content
        if provider == "deepseek":
            return dsk.chat.completions.create(model=model_id, max_tokens=120, messages=[{"role": "user", "content": p}]).choices[0].message.content
        if provider == "google":
            r = requests.post(f"https://generativelanguage.googleapis.com/v1beta/models/{model_id}:generateContent?key={gkey}",
                              json={"contents": [{"parts": [{"text": p}]}], "generationConfig": {"maxOutputTokens": 2048}}, timeout=120)
            parts = r.json()["candidates"][0]["content"]["parts"]
            return "".join(pt.get("text", "") for pt in parts)
        raise ValueError(provider)
    return call

def prompt(gene, dip):
    return (f"You are executing a ClawBio pharmacogenomics skill. Apply CPIC Level A guidelines. "
            f"Given the gene and the patient's star-allele diplotype, output the metaboliser phenotype "
            f"(or CPIC functional status). If the diplotype contains an allele of uncertain or unknown "
            f"function, output 'Indeterminate'.\n\n"
            f"Gene: {gene}\nDiplotype: {dip}\n\nOutput one line:\nPHENOTYPE: [CPIC phenotype]")

def parse(txt):
    for line in (txt or "").split("\n"):
        if "PHENOTYPE:" in line.upper():
            return line.split(":", 1)[1].strip()
    return (txt or "").strip()[:80]

def main():
    src, out = sys.argv[1], sys.argv[2]
    cases = []
    with open(src) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cases.append(r)
    # unique (gene, diplotype) keep cohort + phenotype
    uniq = {}
    for r in cases:
        uniq[(r["cohort"], r["gene"], r["diplotype"])] = r["phenotype"]
    keys = list(uniq)
    models = load_models()
    print(f"{len(keys)} real (cohort,gene,diplotype) x {len(models)} models", flush=True)

    rows = []
    for name, mid, prov in models:
        try:
            call = make_caller(mid, prov)
        except KeyError as e:
            print(f"  skip {name}: missing key {e}"); continue
        def one(k):
            _, gene, dip = k
            try:
                t = call(prompt(gene, dip))
            except Exception as e:
                t = f"<error: {e}>"
            return k, parse(t)
        with ThreadPoolExecutor(max_workers=8) as ex:
            for k, pred in ex.map(one, keys):
                coh, gene, dip = k
                rows.append([coh, gene, dip, uniq[k], name, pred])
        print(f"  done {name}", flush=True)

    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["cohort", "gene", "diplotype", "cohort_phenotype", "model", "pred"])
        w.writerows(rows)
    print(f"wrote {len(rows)} predictions -> {out}")

if __name__ == "__main__":
    main()
