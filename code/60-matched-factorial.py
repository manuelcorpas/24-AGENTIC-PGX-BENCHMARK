#!/usr/bin/env python3
"""
Matched factorial rerun (revision item N1; Reviewer 1 point 3, Reviewer 2 point 1).

THE PROBLEM THIS FIXES
The original five conditions differed in more than the variable of interest.
The skill arms were told the target drug and were scored drug-specifically; the
free-prompted and retrieval arms were not, and had to infer which drug the
question was about. So the skill arms' accuracy lead confounds "executed rules
help" with "an easier task". Reviewers 1.3 and 2.1 are correct about this and
the fix is to hold the task identical and vary only two things:

    knowledge representation  x  decision mechanism

                    | model generation      | deterministic execution
    ----------------+-----------------------+--------------------------
    no knowledge    | free_generation       | UNDEFINED, omitted
    retrieved prose | rag_generation        | rag_execution
    structured rules| skill_generation      | skill_execution

The (no knowledge, execution) cell is undefined, because execution needs a
structured input call and there is nothing to produce one from. We omit it and
say so rather than inventing a cell.

rag_execution is the informative new cell: the model extracts a structured
diplotype call from retrieved prose, and the validated skill executes the
mapping. It separates extraction quality from rule quality, which nothing in
the original design did.

WHAT IS RERUN AND WHAT IS REUSED
Only the three cells marked rerun=True issue new model calls. The two skill
arms are reused unchanged from the existing 17,820-evaluation dataset, because
the plan of revision declares their output schema canonical and the matched
design adopts that schema rather than re-running them.

Single framing: the three population framings are the same genotypes with
different wording, so they are not crossed into the matched grid. Framing
invariance is a separate, smaller check (R2.1).

USAGE
    python code/60-matched-factorial.py --dry-run            # free; writes prompts only
    python code/60-matched-factorial.py --estimate           # print the call budget
    python code/60-matched-factorial.py --cells rag_execution --limit 5
    python code/60-matched-factorial.py                      # full run (COSTS MONEY)

Scored by code/61-rescore-matched.py under both the baseline scorer and the
frozen clinical-equivalence scorer (see SCORING-PREREG.md).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
CORPUS_FILE = BASE / "specs" / "cpic_rag_corpus_v3.json"
OUT = BASE / "data" / "v3_matched_factorial.json"
JSONL = BASE / "data" / "v3_matched_factorial.jsonl"
PROMPTS_OUT = BASE / "data" / "v3_matched_factorial_prompts.json"

_spec = spec_from_file_location("_pgx_rules", Path(__file__).resolve().parent / "_pgx_rules.py")
rules = module_from_spec(_spec)
_spec.loader.exec_module(rules)

CASES = rules.load_cases()
DIP2PHEN, REC, RULES_DIP, RULES_REC = rules.build_rules(CASES)
execute_skill = rules.make_executor(DIP2PHEN, REC)

N_REPS = 3
N_MODELS = 9

# Populated by load_models() only when a live run is requested, so that import,
# --dry-run and the tests need no API keys at all.
MODELS: dict = {}


# =============================================================================
# The factorial
# =============================================================================

CELLS: dict[str, dict] = {
    "free_generation": {
        "knowledge": "none", "mechanism": "model_generation", "rerun": True,
        "label": "free-prompted",
    },
    "rag_generation": {
        "knowledge": "retrieved_prose", "mechanism": "model_generation", "rerun": True,
        "label": "retrieval-augmented",
    },
    "rag_execution": {
        "knowledge": "retrieved_prose", "mechanism": "deterministic_execution", "rerun": True,
        "label": "extract-then-execute",
    },
    "skill_generation": {
        "knowledge": "structured_rules", "mechanism": "model_generation", "rerun": False,
        "label": "skill-reasoning (reused)",
    },
    "skill_execution": {
        "knowledge": "structured_rules", "mechanism": "deterministic_execution", "rerun": False,
        "label": "skill-execution (reused)",
    },
}

# One instruction and one schema for every generation cell. Identical text, so
# the only thing that varies is the knowledge block. This is the matched design.
INSTRUCTION_GENERATION = (
    "You are a pharmacogenomics interpretation system. Determine the patient's "
    "diplotype from the genotype below, then give the phenotype and the "
    "recommendation for the target drug named in the Patient section."
)

INSTRUCTION_EXECUTION = (
    "You are the input-interpretation step of a pharmacogenomics agent. Map the "
    "patient's genotype to exactly ONE diplotype from the controlled list below, "
    "copying its text VERBATIM. Do not invent notation. A downstream validated "
    "skill computes the phenotype and recommendation."
)

SCHEMA_GENERATION = """## Output (4 lines only)
DIPLOTYPE: [called diplotype]
PHENOTYPE: [phenotype]
DRUG: [recommendation for the target drug]
HAZARD: [clinical hazard]"""

SCHEMA_EXECUTION = """## Output (1 line only)
DIPLOTYPE: [exact text of one list item]"""


def _corpus_excerpt(gene: str) -> str:
    corpus = json.loads(CORPUS_FILE.read_text())["genes"]
    entry = corpus.get(gene)
    if not entry:
        return f"(no CPIC guideline excerpt retrieved for {gene})"
    return entry["guideline_excerpt"]


def _valid_diplotype_list(gene: str) -> str:
    return "\n".join(f"  - {d}" for d in sorted(RULES_DIP[gene].keys()))


def knowledge_block(cell: str, case: dict) -> str:
    """The ONLY part of a matched prompt that varies between cells."""
    gene = case["gene"]
    kind = CELLS[cell]["knowledge"]
    if kind == "none":
        return "## Knowledge provided\n\n(none: answer from your own knowledge)"
    if kind == "retrieved_prose":
        return (f"## Knowledge provided\n\nCPIC guideline excerpt "
                f"(retrieved for gene: {gene})\n\n{_corpus_excerpt(gene)}")
    return (f"## Knowledge provided\n\n"
            f"{rules.skill_rules_text(gene, RULES_DIP, RULES_REC)}")


def build_prompt(cell: str, case: dict) -> str:
    """Compose a prompt for one cell.

    Generation cells are byte-identical apart from knowledge_block(), which is
    what makes the comparison matched. Execution cells share one instruction,
    one interface (the controlled diplotype vocabulary) and one schema; they
    differ in the knowledge they are given, which is the point of the cell.
    """
    if CELLS[cell]["mechanism"] == "model_generation":
        return "\n\n".join([
            INSTRUCTION_GENERATION,
            knowledge_block(cell, case),
            "## Patient\n"
            f"Gene: {case['gene']}\n"
            f"Genotype: {case['genotype']}\n"
            f"Drug: {case['drug']}",
            SCHEMA_GENERATION,
        ])
    return "\n\n".join([
        INSTRUCTION_EXECUTION,
        knowledge_block(cell, case),
        "## Patient\n"
        f"Gene: {case['gene']}\n"
        f"Genotype: {case['genotype']}\n"
        f"Drug: {case['drug']}\n\n"
        f"Valid diplotypes (choose one, copy verbatim):\n{_valid_diplotype_list(case['gene'])}",
        SCHEMA_EXECUTION,
    ])


def planned_calls(n_models: int = N_MODELS, n_reps: int = N_REPS,
                  n_cases: int | None = None) -> int:
    """New model calls this runner issues. The figure quoted to the editor is
    computed here rather than typed into a document by hand."""
    rerun_cells = sum(1 for c in CELLS.values() if c["rerun"])
    return rerun_cells * (n_cases or len(CASES)) * n_models * n_reps


# =============================================================================
# Execution
# =============================================================================

def load_models() -> dict:
    """Construct provider clients. Called only for a live run, never at import."""
    import anthropic
    import openai

    env_path = BASE / ".env"
    keys = dict(os.environ)
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, _, v = line.partition("=")
                keys.setdefault(k.strip(), v.strip().strip('"').strip("'"))

    def key(*names):
        for nm in names:
            if keys.get(nm):
                return keys[nm]
        raise KeyError(f"Missing API key: set one of {names} in the environment or {env_path}")

    ant = anthropic.Anthropic(api_key=key("ANTHROPIC_API_KEY"))
    oai = openai.OpenAI(api_key=key("OPENAI_API_KEY"))
    dsk = openai.OpenAI(api_key=key("DEEPSEEK_API_KEY"), base_url="https://api.deepseek.com")
    gkey = key("GEMINI_API_KEY", "GOOGLE_API_KEY")
    mkey = key("MISTRAL_API_KEY")

    sem = {"anthropic": threading.Semaphore(20), "openai": threading.Semaphore(20),
           "deepseek": threading.Semaphore(16), "gemini": threading.Semaphore(10),
           "mistral": threading.Semaphore(4)}

    def _ant(model, p):
        with sem["anthropic"]:
            return ant.messages.create(model=model, max_tokens=320,
                                       messages=[{"role": "user", "content": p}]).content[0].text

    def _oai(model, p, reasoning=False):
        with sem["openai"]:
            if reasoning:
                return oai.chat.completions.create(
                    model=model, max_completion_tokens=2000,
                    messages=[{"role": "user", "content": p}]).choices[0].message.content
            return oai.chat.completions.create(
                model=model, max_tokens=320,
                messages=[{"role": "user", "content": p}]).choices[0].message.content

    def _dsk(p):
        with sem["deepseek"]:
            return dsk.chat.completions.create(
                model="deepseek-chat", max_tokens=320,
                messages=[{"role": "user", "content": p}]).choices[0].message.content

    def _gem(p):
        import urllib.request
        with sem["gemini"]:
            req = urllib.request.Request(
                "https://generativelanguage.googleapis.com/v1beta/models/"
                f"gemini-2.5-flash:generateContent?key={gkey}",
                data=json.dumps({"contents": [{"parts": [{"text": p}]}]}).encode(),
                headers={"Content-Type": "application/json"})
            body = json.loads(urllib.request.urlopen(req, timeout=120).read())
            return body["candidates"][0]["content"]["parts"][0]["text"]

    def _mis(p):
        import urllib.request
        with sem["mistral"]:
            req = urllib.request.Request(
                "https://api.mistral.ai/v1/chat/completions",
                data=json.dumps({"model": "mistral-large-2411", "max_tokens": 320,
                                 "messages": [{"role": "user", "content": p}]}).encode(),
                headers={"Content-Type": "application/json",
                         "Authorization": f"Bearer {mkey}"})
            body = json.loads(urllib.request.urlopen(req, timeout=120).read())
            return body["choices"][0]["message"]["content"]

    return {
        "Claude Opus 4":    lambda p: _ant("claude-opus-4-20250514", p),
        "Claude Sonnet 4":  lambda p: _ant("claude-sonnet-4-20250514", p),
        "GPT-5.2":          lambda p: _oai("gpt-5.2", p, reasoning=True),
        "GPT-4.1":          lambda p: _oai("gpt-4.1", p),
        "o3":               lambda p: _oai("o3", p, reasoning=True),
        "o4-mini":          lambda p: _oai("o4-mini", p, reasoning=True),
        "Gemini 2.5 Flash": lambda p: _gem(p),
        "DeepSeek V3":      lambda p: _dsk(p),
        "Mistral Large 2":  lambda p: _mis(p),
    }


def run_one(cell: str, case: dict, model_name: str, fn, rep: int) -> dict:
    """One evaluation. Execution cells run the model for the CALL only, then the
    validated skill computes phenotype and recommendation in code."""
    prompt = build_prompt(cell, case)
    text, error = "", None
    for attempt in range(3):
        try:
            text = fn(prompt)
            break
        except Exception as exc:                        # noqa: BLE001 - recorded, not swallowed
            error = f"{type(exc).__name__}: {exc}"
    row = {
        "cell": cell, "knowledge": CELLS[cell]["knowledge"],
        "mechanism": CELLS[cell]["mechanism"],
        "case_id": case["id"], "gene": case["gene"], "drug": case["drug"],
        "model": model_name, "rep": rep,
        "raw": text, "error": error,
        "called_diplotype": rules.parse_field(text, "DIPLOTYPE"),
    }
    if CELLS[cell]["mechanism"] == "deterministic_execution":
        phen, rec = execute_skill(case["gene"], case["drug"], row["called_diplotype"])
        row["parsed_phenotype"] = phen or ""
        row["parsed_drug"] = rec or ""
        row["abstained"] = phen is None
    else:
        row["parsed_phenotype"] = rules.parse_field(text, "PHENOTYPE")
        row["parsed_drug"] = rules.parse_field(text, "DRUG")
        row["abstained"] = False
    return row


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--dry-run", action="store_true",
                    help="Build every prompt and write them out; issue no API calls")
    ap.add_argument("--estimate", action="store_true", help="Print the call budget and exit")
    ap.add_argument("--limit", type=int, default=0, help="Use only the first N cases")
    ap.add_argument("--reps", type=int, default=N_REPS)
    ap.add_argument("--cells", nargs="*", default=None,
                    help="Restrict to these cells (default: all rerun cells)")
    ap.add_argument("--models", nargs="*", default=None, help="Restrict to these models")
    args = ap.parse_args(argv)

    cases = CASES[:args.limit] if args.limit else CASES
    cells = args.cells or [c for c, v in CELLS.items() if v["rerun"]]
    for c in cells:
        if c not in CELLS:
            sys.stderr.write(f"unknown cell: {c}\n")
            return 2

    if args.estimate:
        print(f"cases={len(cases)} models={N_MODELS} reps={args.reps}")
        print(f"rerun cells: {[c for c, v in CELLS.items() if v['rerun']]}")
        print(f"planned new model calls: {planned_calls(N_MODELS, args.reps, len(cases))}")
        return 0

    if args.dry_run:
        prompts = [
            {"cell": cell, "case_id": c["id"], "prompt": build_prompt(cell, c)}
            for cell in cells for c in cases
        ]
        PROMPTS_OUT.parent.mkdir(exist_ok=True)
        PROMPTS_OUT.write_text(json.dumps(prompts, indent=2))
        print(f"DRY RUN: built {len(prompts)} prompts across {len(cells)} cells, "
              f"issued 0 API calls, wrote {PROMPTS_OUT.name}")
        return 0

    models = MODELS or load_models()
    if args.models:
        models = {k: v for k, v in models.items() if k in args.models}

    tasks = [(cell, c, mn, fn, rep)
             for cell in cells for c in cases
             for mn, fn in models.items() for rep in range(args.reps)]
    print(f"issuing {len(tasks)} model calls across {len(cells)} cells")

    results = []
    JSONL.parent.mkdir(exist_ok=True)
    with JSONL.open("a") as sink, ThreadPoolExecutor(max_workers=24) as pool:
        for row in pool.map(lambda t: run_one(*t), tasks):
            results.append(row)
            sink.write(json.dumps(row) + "\n")
            sink.flush()
    OUT.write_text(json.dumps(results, indent=2))
    print(f"wrote {len(results)} evaluations to {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
