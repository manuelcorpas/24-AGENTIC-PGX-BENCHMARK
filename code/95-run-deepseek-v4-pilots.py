#!/usr/bin/env python3
"""Run one pinned DeepSeek V4 Flash request per paper prompt family."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
CODE = BASE / "code"
AUDIT = BASE / "docs" / "deepseek-v4-flash-rerun-audit-2026-08-30"
DEFAULT_OUT = AUDIT / "pilots" / "deepseek-v4-flash-pilots.jsonl"


def load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def family(family_id: str, dependency_class: str, runner: str, prompt: str,
           max_tokens: int) -> dict:
    return {
        "family_id": family_id,
        "dependency_class": dependency_class,
        "runner": runner,
        "prompt": prompt,
        "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
        "max_tokens": max_tokens,
    }


def build_pilot_prompts() -> list[dict]:
    mf = load(CODE / "60-matched-factorial.py", "pilot_matched_factorial")
    case = mf.CASES[0]
    rows = [
        family(f"matched_{cell}", "matched_five_cell", "code/60-matched-factorial.py",
               mf.build_prompt(cell, case), 320)
        for cell in mf.CELLS
    ]

    rules = load(CODE / "67-extract-rules-from-prose.py", "pilot_extract_rules")
    corpus = json.loads(rules.CORPUS_FILE.read_text())["genes"]
    gene = sorted(corpus)[0]
    rows.append(family(
        "extracted_rule_table", "extracted_rules", "code/67-extract-rules-from-prose.py",
        rules.PROMPT.format(gene=gene, excerpt=corpus[gene]["guideline_excerpt"]), 8000,
    ))

    forward = load(CODE / "14-adversarial-scrambled-spec.py", "pilot_adversarial_forward")
    reverse = load(CODE / "14b-adversarial-reverse-spec.py", "pilot_adversarial_reverse")
    cases = {c["id"]: c for c in json.loads(forward.SPECS.read_text())}
    pop = {"id": "EUR", "name": "European family cohort (Corpasome project)",
           "desc": "European ancestry, whole-genome sequencing"}
    s = forward.SCRAMBLES[0]
    c = cases[s["case_id"]]
    rows.append(family(
        "adversarial_forward", "adversarial_contract", "code/14-adversarial-scrambled-spec.py",
        forward.with_spec_prompt(s["case_id"], c["gene"], c["genotype"], pop["name"], pop["desc"],
                                 s["scrambled_diplotype"], s["scrambled_phen"],
                                 s["scrambled_activity"], s["scrambled_drug"],
                                 c["pop_note"]["EUR"]), 600,
    ))
    s = reverse.SCRAMBLES[0]
    c = cases[s["case_id"]]
    rows.append(family(
        "adversarial_reverse", "adversarial_contract", "code/14b-adversarial-reverse-spec.py",
        reverse.with_spec_prompt(c["gene"], c["genotype"], pop["name"], pop["desc"],
                                 s["scrambled_diplotype"], s["scrambled_phen"],
                                 s["scrambled_activity"], s["scrambled_drug"],
                                 c["pop_note"]["EUR"]), 600,
    ))

    previous = mf.RETRIEVAL_MODE
    mf.RETRIEVAL_MODE = "drug"
    try:
        for cell in ("rag_generation", "rag_execution"):
            rows.append(family(
                f"drug_keyed_{cell}", "drug_keyed_retrieval", "code/62-rag-drugkeyed-fullgrid.py",
                mf.build_prompt(cell, case), 320,
            ))
    finally:
        mf.RETRIEVAL_MODE = previous

    population = load(CODE / "02-run-benchmark-v3.py", "pilot_population")
    pcase = json.loads(population.SPECS_FILE.read_text())[0]
    ppop = population.POPULATIONS[0]
    rows.extend([
        family("population_no_spec", "population_framing", "code/02-run-benchmark-v3.py",
               population.no_spec_prompt(pcase, ppop), 500),
        family("population_with_spec", "population_framing", "code/02-run-benchmark-v3.py",
               population.with_spec_prompt(pcase, ppop), 500),
    ])
    rag = load(CODE / "16-run-rag-condition.py", "pilot_population_rag")
    rag_corpus = json.loads(rag.CORPUS_FILE.read_text())
    excerpt = rag_corpus["genes"][pcase["gene"]]["guideline_excerpt"]
    rows.append(family(
        "population_cpic_rag", "population_framing", "code/16-run-rag-condition.py",
        rag.cpic_rag_prompt(pcase, ppop, excerpt), 600,
    ))

    placeholder_keys = ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "DEEPSEEK_API_KEY",
                        "GEMINI_API_KEY", "MISTRAL_API_KEY")
    inserted = [name for name in placeholder_keys if name not in os.environ]
    for name in inserted:
        os.environ[name] = "pilot-prompt-construction-only"
    try:
        skills = load(CODE / "43-armAB-fullgrid.py", "pilot_population_skills")
    finally:
        for name in inserted:
            os.environ.pop(name, None)
    scase = skills.h.cases[0]
    spop = skills.POPULATIONS[0]
    rows.extend([
        family("population_skill_reasoning", "population_framing", "code/43-armAB-fullgrid.py",
               skills.armA_prompt_pop(scase, spop), 320),
        family("population_skill_execution", "population_framing", "code/43-armAB-fullgrid.py",
               skills.armB2_prompt_pop(scase, spop), 320),
    ])

    normal = load(CODE / "74-input-normalisation.py", "pilot_normalisation")
    inputs = json.loads(normal.INPUTS_OUT.read_text())["inputs"]
    nkey = sorted(inputs)[0]
    rec = inputs[nkey]
    rows.extend([
        family("normalisation_without_definitions", "input_normalisation",
               "code/74-input-normalisation.py", normal.build_prompt(
                   rec["gene"], rec["variants"], "vcf", rec["truncated"], None), 6000),
        family("normalisation_with_definitions", "input_normalisation",
               "code/74-input-normalisation.py", normal.build_prompt(
                   rec["gene"], rec["variants"], "vcf", rec["truncated"],
                   normal.load_definitions(rec["gene"])), 6000),
    ])

    real = load(BASE / "real-genome-arm" / "scripts" / "21_run_agent_four_cohorts.py",
                "pilot_real_genome")
    state = real.read_states(BASE / "data" / "n5_four_cohorts.tsv")[0]
    rows.append(family(
        "real_genome_interpretation", "real_genome",
        "real-genome-arm/scripts/21_run_agent_four_cohorts.py",
        real.build_prompt(state["gene"], state["diplotype"]), 1024,
    ))
    return rows


def read_key() -> str:
    if os.environ.get("DEEPSEEK_API_KEY"):
        return os.environ["DEEPSEEK_API_KEY"]
    env_path = BASE / ".env"
    for line in env_path.read_text().splitlines() if env_path.exists() else []:
        if line.strip().startswith("DEEPSEEK_API_KEY="):
            return line.partition("=")[2].strip().strip('"').strip("'")
    raise SystemExit("DEEPSEEK_API_KEY is not set")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="list families without API calls")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    prompts = build_pilot_prompts()
    if args.list:
        for row in prompts:
            print(row["family_id"], row["runner"], row["prompt_sha256"])
        return 0

    import openai
    client = openai.OpenAI(api_key=read_key(), base_url="https://api.deepseek.com")
    ds_v4 = load(CODE / "94-deepseek-v4-client.py", "deepseek_v4_client_pilots")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise SystemExit(f"refusing to overwrite existing pilot evidence: {args.out}")
    with args.out.open("x") as sink:
        for index, row in enumerate(prompts, 1):
            result = ds_v4.call_deepseek_v4(
                client, row["prompt"], max_tokens=row["max_tokens"],
                max_attempts=3, retry_base_seconds=2.0,
            )
            record = {**row, "deepseek_provenance": result}
            sink.write(json.dumps(record, ensure_ascii=False) + "\n")
            sink.flush()
            os.fsync(sink.fileno())
            print(f"[{index}/{len(prompts)}] {row['family_id']}: {result['stop_reason']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
