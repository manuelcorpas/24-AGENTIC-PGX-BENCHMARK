#!/usr/bin/env python3
"""
(gene, drug)-keyed chunking comparator for the cpic_rag arm.

Reviewer concern: the drug-substitution failure (Section 4) is conditional on
gene-keyed chunking. This script re-runs the retrieval-augmented condition on the
six worst-affected genes (CYP2D6, CYP2C19, CYP2C9, UGT1A1, SLCO1B1, IFNL3) using
finer (gene, drug)-keyed chunks: the model sees ONLY the CPIC annotation section
for the queried drug, not the full multi-drug gene excerpt. We then measure the
drug-substitution rate under the finer chunking and compare it to the gene-keyed
baseline on the same cells.

Same drug-substitution definition as 33-classify-a2-regressions.py:
queried drug name absent (case-insensitive substring) from the parsed DRUG field.

Run: 3 models x 63 cases x EUR x 2 reps = 378 calls.
"""
from __future__ import annotations
import json, os, re, time
from pathlib import Path
import openai, anthropic

BASE = Path(__file__).resolve().parent.parent
CORPUS = BASE / "SPECS" / "cpic_rag_corpus_v3.json"
CASES = BASE / "SPECS" / "test_cases_v3.json"
BASELINE = BASE / "RESULTS" / "v3_raw_rescored_three_arm.json"
OUT = BASE / "RESULTS" / "v3_rag_genedrug_chunking.json"
REPORT = BASE / "RESULTS" / "v3_rag_genedrug_chunking_report.txt"
ENV = Path("/Users/manuelcorpas1/dev/AGENTIC-AI/.env")

WORST = ["CYP2D6", "CYP2C19", "CYP2C9", "UGT1A1", "SLCO1B1", "IFNL3"]
N_REPS = 2

keys = {}
for line in ENV.read_text().splitlines():
    line = line.strip()
    if "=" in line and not line.startswith("#"):
        k, _, v = line.partition("="); keys[k.strip()] = v.strip().strip('"').strip("'")

ant = anthropic.Anthropic(api_key=keys["ANTHROPIC_API_KEY"])
oai = openai.OpenAI(api_key=keys["OPENAI_API_KEY"])
dsk = openai.OpenAI(api_key=keys["DEEPSEEK_API_KEY"], base_url="https://api.deepseek.com")

def call_ant(p):
    return ant.messages.create(model="claude-opus-4-20250514", max_tokens=600,
                               messages=[{"role": "user", "content": p}]).content[0].text
def call_oai(p):
    try:
        return oai.chat.completions.create(model="gpt-5.2", max_tokens=600,
                                           messages=[{"role": "user", "content": p}]).choices[0].message.content
    except Exception:
        return oai.chat.completions.create(model="gpt-5.2", max_completion_tokens=2000,
                                           messages=[{"role": "user", "content": p}]).choices[0].message.content
def call_dsk(p):
    return dsk.chat.completions.create(model="deepseek-chat", max_tokens=600,
                                       messages=[{"role": "user", "content": p}]).choices[0].message.content
MODELS = {"Claude Opus 4": call_ant, "GPT-5.2": call_oai, "DeepSeek V3": call_dsk}

DRUG_ALIASES = {  # case drug -> tokens that identify its corpus chunk
    "peg-ifn-α": ["peginterferon"],
}

def split_chunks(excerpt):
    parts = re.split(r"(?=^### Annotation)", excerpt, flags=re.M)
    chunks = {}
    for p in parts:
        if not p.strip():
            continue
        m = re.search(r"Drug scope:\s*([^\n]+)", p)
        if m:
            for d in re.split(r"[;,/]| and ", m.group(1)):
                d = d.strip().lower()
                if d:
                    chunks.setdefault(d, []).append(p.strip())
    return {d: "\n\n".join(v) for d, v in chunks.items()}

def drug_chunk(gene_chunks, drug):
    dl = drug.lower()
    toks = DRUG_ALIASES.get(dl, [dl])
    sel = []
    for scope, text in gene_chunks.items():
        if any(t in scope or scope in t for t in toks):
            sel.append(text)
    return "\n\n".join(sel) if sel else None

def prompt(tc, pop, excerpt):
    return f"""You are a pharmacogenomics interpretation system. Use the CPIC guideline excerpt below to decide the diplotype, phenotype, drug recommendation, hazard, and population-specific note for the patient.

## CPIC guideline excerpt (retrieved for gene: {tc['gene']}, drug: {tc['drug']})

{excerpt}

## Patient

Gene: {tc['gene']}
Genotype: {tc['genotype']}
Patient cohort: {pop['name']} ({pop['desc']})

## Output (exactly 5 lines, no preamble)

DIPLOTYPE: [your answer]
PHENOTYPE: [your answer]
DRUG: [drug]: [recommendation]
HAZARD: [clinical hazard]
POPULATION: [population-specific note]"""

def parse_drug(text):
    if not text:
        return ""
    for line in text.split("\n"):
        if "DRUG:" in line.upper():
            return line[line.upper().index("DRUG:") + 5:].strip()
    return ""

def is_substitution(queried, drug_field):
    q = (queried or "").lower()
    return bool(q) and q not in (drug_field or "").lower()

def main():
    corpus = json.loads(CORPUS.read_text())["genes"]
    cases = [c for c in json.loads(CASES.read_text()) if c["gene"] in WORST]
    pop = {"id": "EUR", "name": "European family cohort (Corpasome project)", "desc": "European ancestry, whole-genome sequencing"}
    gene_chunks = {g: split_chunks(corpus[g]["guideline_excerpt"]) for g in WORST}

    results = []
    total = len(cases) * len(MODELS) * N_REPS
    n = 0
    for tc in cases:
        ex = drug_chunk(gene_chunks[tc["gene"]], tc["drug"])
        if ex is None:
            print(f"  WARN no chunk for {tc['id']} ({tc['drug']}); skipping")
            continue
        p = prompt(tc, pop, ex)
        for model, fn in MODELS.items():
            for rep in range(N_REPS):
                n += 1
                try:
                    time.sleep(0.4)
                    text = fn(p)
                except Exception as e:
                    text = f"<error: {e}>"
                df = parse_drug(text)
                sub = is_substitution(tc["drug"], df)
                results.append({"tc": tc["id"], "gene": tc["gene"], "drug": tc["drug"],
                                "model": model, "rep": rep, "parsed_drug": df,
                                "substitution": sub, "chunk_chars": len(ex)})
                print(f"  [{n}/{total}] {model} {tc['id']} r{rep}: sub={sub}")
    OUT.write_text(json.dumps(results, indent=2))

    # gene-keyed baseline on the same cells (EUR, these genes) from merged dataset
    base = json.loads(BASELINE.read_text())
    casemap = {c["id"]: c for c in json.loads(CASES.read_text())}
    base_rows = [r for r in base if r["cond"] == "cpic_rag" and r["gene"] in WORST and r["pop"] == "EUR"
                 and not r["scores"].get("format_fail", False)]
    def base_sub(r):
        return is_substitution(casemap[r["tc"]]["drug"], r["parsed"].get("DRUG", ""))

    out = ["(gene, drug)-KEYED CHUNKING COMPARATOR", "",
           f"New (gene,drug) chunking run: {len(results)} calls, {len(MODELS)} models, EUR, {N_REPS} reps", ""]
    out.append(f"{'gene':<10}{'genekey_sub%':>14}{'genedrug_sub%':>15}{'n_new':>7}")
    from collections import defaultdict
    def rate(rows, fn):
        rows = list(rows)
        return (100 * sum(1 for r in rows if fn(r)) / len(rows)) if rows else float("nan"), len(rows)
    for g in WORST:
        gk, gkn = rate([r for r in base_rows if r["gene"] == g], base_sub)
        gd, gdn = rate([r for r in results if r["gene"] == g], lambda r: r["substitution"])
        out.append(f"{g:<10}{gk:13.1f}%{gd:14.1f}%{gdn:7d}")
    gk_all, _ = rate(base_rows, base_sub)
    gd_all, gdn_all = rate(results, lambda r: r["substitution"])
    out.append(f"{'ALL 6':<10}{gk_all:13.1f}%{gd_all:14.1f}%{gdn_all:7d}")
    out.append("")
    out.append("genekey_sub%  = drug-substitution rate under gene-keyed chunks (published baseline, EUR)")
    out.append("genedrug_sub% = drug-substitution rate under (gene,drug)-keyed chunks (this run)")
    text = "\n".join(out)
    REPORT.write_text(text)
    print("\n" + text)

if __name__ == "__main__":
    main()
