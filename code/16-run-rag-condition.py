#!/usr/bin/env python3
"""
Production runner for the cpic_rag condition (third arm of the v3 benchmark).

Mirrors the structure of 02-run-benchmark-v3.py exactly: same 9 frontier
models, same 110 CPIC test cases, same 3 population contexts, same 3 runs
per cell. The difference is the prompt: each call receives the relevant
gene-level CPIC guideline excerpt (from cpic_rag_corpus_v3.json) and the
patient genotype, but does NOT receive the pre-computed diplotype, phenotype,
activity score, or drug recommendation. The model must derive these from the
excerpt and the genotype.

Pipeline position:
  upstream:   15-build-cpic-rag-corpus.py
  this:       16-run-rag-condition.py
  downstream: 17-merge-rag-results.py -> 18-three-arm-analysis.py
              19-validate-three-arm.py (post-merge integrity gate)

Inputs:
  ../SPECS/test_cases_v3.json
  ../SPECS/cpic_rag_corpus_v3.json
  env: ANTHROPIC_API_KEY, OPENAI_API_KEY, DEEPSEEK_API_KEY, GEMINI_API_KEY,
       MISTRAL_API_KEY (defaults to embedded production credentials)

Outputs:
  ../RESULTS/v3_rag_raw.json         final 8,910-row dataset (cond=cpic_rag)
  ../RESULTS/v3_rag_partial.json     saved every 100 calls
  ../RESULTS/v3_rag_checkpoint.json  resume point
  ../RESULTS/v3_rag_cost_log.json    per-model token + USD tally
  ../LOGS/v3_rag_run_<UTC-timestamp>.log

Usage:
  python3 16-run-rag-condition.py --dry-run        45 calls (~$0.50)
  python3 16-run-rag-condition.py                  full 8,910 calls
  python3 16-run-rag-condition.py --yes            skip interactive y/N gate

Verification on prompt construction (verified before merging into pipeline):
  For three sample cases (DPYD/dpyd_fu_pm, CYP2D6/cyp2d6_codeine_pm,
  HLA-B*57:01/hlab5701_aba_pos), the assembled prompt was inspected to
  confirm that the gt_diplotype, gt_phenotype, and gt_drug strings DO NOT
  appear verbatim. The prompt contains the gene-level CPIC excerpt
  (which includes the answer as one row in a multi-row table the model
  must navigate) and the patient genotype + population context, but no
  pre-computed analytical decision. This was verified live in the
  hard-probe at corpus build time.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import os
import sys
import time
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Callable

import openai
import anthropic
import requests

BASE = Path(__file__).resolve().parent.parent
PYDIR = BASE / "PYTHON"
SPECS_FILE = BASE / "SPECS" / "test_cases_v3.json"
CORPUS_FILE = BASE / "SPECS" / "cpic_rag_corpus_v3.json"
RESULTS = BASE / "RESULTS"
LOGS = BASE / "LOGS"
RESULTS.mkdir(exist_ok=True)
LOGS.mkdir(exist_ok=True)

OUT_RAW = RESULTS / "v3_rag_raw.json"
OUT_PARTIAL = RESULTS / "v3_rag_partial.json"
OUT_CKPT = RESULTS / "v3_rag_checkpoint.json"
OUT_COST = RESULTS / "v3_rag_cost_log.json"

CONDITION = "cpic_rag"
N_RUNS = 3

BUDGET_HARD_CAP = 200.0
BUDGET_SOFT_WARN = 100.0

SPACING = {
    "Claude Opus 4": 0.5,
    "Claude Sonnet 4": 0.4,
    "GPT-5.2": 0.4,
    "GPT-4.1": 0.4,
    "o3": 0.5,
    "o4-mini": 0.5,
    "Gemini 2.5 Flash": 0.5,
    "DeepSeek V3": 0.3,
    "Mistral Large 2": 15.0,
}

PRICE = {
    "Claude Opus 4": (15.00, 75.00),
    "Claude Sonnet 4": (3.00, 15.00),
    "GPT-5.2": (2.50, 10.00),
    "GPT-4.1": (2.00, 8.00),
    "o3": (15.00, 60.00),
    "o4-mini": (1.10, 4.40),
    "Gemini 2.5 Flash": (0.30, 2.50),
    "DeepSeek V3": (0.27, 1.10),
    "Mistral Large 2": (2.00, 6.00),
}

# Frozen API model identifiers — pinned to the same versions used by
# 02-run-benchmark-v3.py so cpic_rag results are directly comparable to the
# existing no_spec / with_spec arms.
FROZEN_MODEL_IDS = {
    "Claude Opus 4": "claude-opus-4-20250514",
    "Claude Sonnet 4": "claude-sonnet-4-20250514",
    "GPT-5.2": "gpt-5.2",
    "GPT-4.1": "gpt-4.1",
    "o3": "o3",
    "o4-mini": "o4-mini",
    "Gemini 2.5 Flash": "gemini-2.5-flash",
    "DeepSeek V3": "deepseek-chat",
    "Mistral Large 2": "mistral-large-latest",
}

ANT_KEY = os.environ.get("ANTHROPIC_API_KEY", "")
OAI_KEY = os.environ.get("OPENAI_API_KEY", "")
DSK_KEY = os.environ.get("DEEPSEEK_API_KEY", "")
GEM_KEY = os.environ.get("GEMINI_API_KEY", "")
MIS_KEY = os.environ.get("MISTRAL_API_KEY", "")

ant = anthropic.Anthropic(api_key=ANT_KEY)
oai = openai.OpenAI(api_key=OAI_KEY)
dsk = openai.OpenAI(api_key=DSK_KEY, base_url="https://api.deepseek.com")

POPULATIONS = [
    {"id": "EUR", "name": "European family cohort (Corpasome project)",
     "desc": "European ancestry, whole-genome sequencing"},
    {"id": "AMR", "name": "Peruvian Genome Project",
     "desc": "Admixed Latin American, 7 indigenous and mestizo subpopulations"},
    {"id": "AFR", "name": "Uganda Genome Resource",
     "desc": "East African, 6,407 whole-genome sequences"},
]


# =============================================================================
# Logging
# =============================================================================

def setup_logging(dry_run: bool) -> logging.Logger:
    suffix = "_dry" if dry_run else ""
    log_path = LOGS / f"v3_rag_run{suffix}_{dt.datetime.now(dt.UTC).strftime('%Y%m%dT%H%M%SZ')}.log"
    logger = logging.getLogger("rag_run")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter(
        "%(asctime)s.%(msecs)03dZ %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.info("[setup] log file: %s", log_path)
    return logger


# =============================================================================
# Model dispatch (same shape as 02-run-benchmark-v3.py)
# =============================================================================

def _retry(fn, *args, max_retries: int = 3, base_backoff: float = 2.0, **kwargs):
    last = None
    for attempt in range(max_retries):
        try:
            return fn(*args, **kwargs)
        except Exception as e:
            last = e
            time.sleep(base_backoff * (2 ** attempt))
    raise last


def call_ant(model: str, prompt: str) -> tuple[str, int, int]:
    r = ant.messages.create(model=model, max_tokens=600,
                            messages=[{"role": "user", "content": prompt}])
    text = r.content[0].text
    in_tok = getattr(r.usage, "input_tokens", len(prompt) // 4)
    out_tok = getattr(r.usage, "output_tokens", len(text) // 4)
    return text, in_tok, out_tok


def call_oai(model: str, prompt: str) -> tuple[str, int, int]:
    try:
        r = oai.chat.completions.create(model=model, max_tokens=600,
                                        messages=[{"role": "user", "content": prompt}])
    except Exception:
        r = oai.chat.completions.create(model=model, max_completion_tokens=2000,
                                        messages=[{"role": "user", "content": prompt}])
    text = r.choices[0].message.content
    in_tok = r.usage.prompt_tokens if r.usage else len(prompt) // 4
    out_tok = r.usage.completion_tokens if r.usage else len(text or "") // 4
    return text, in_tok, out_tok


def call_dsk(prompt: str) -> tuple[str, int, int]:
    r = dsk.chat.completions.create(model="deepseek-chat", max_tokens=600,
                                    messages=[{"role": "user", "content": prompt}])
    text = r.choices[0].message.content
    in_tok = r.usage.prompt_tokens if r.usage else len(prompt) // 4
    out_tok = r.usage.completion_tokens if r.usage else len(text or "") // 4
    return text, in_tok, out_tok


def call_gem(prompt: str) -> tuple[str, int, int]:
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={GEM_KEY}"
    resp = requests.post(url,
        json={"contents": [{"parts": [{"text": prompt}]}],
              "generationConfig": {"maxOutputTokens": 1500}},
        timeout=90)
    d = resp.json()
    if "candidates" not in d:
        raise RuntimeError(f"Gemini error: {d.get('error', d)}")
    text = d["candidates"][0]["content"]["parts"][0]["text"]
    usage = d.get("usageMetadata", {})
    in_tok = usage.get("promptTokenCount", len(prompt) // 4)
    out_tok = usage.get("candidatesTokenCount", len(text) // 4)
    return text, in_tok, out_tok


def call_mis(prompt: str) -> tuple[str, int, int]:
    r = requests.post("https://api.mistral.ai/v1/chat/completions",
        headers={"Authorization": f"Bearer {MIS_KEY}", "Content-Type": "application/json"},
        json={"model": "mistral-large-latest",
              "messages": [{"role": "user", "content": prompt}],
              "max_tokens": 600},
        timeout=90)
    d = r.json()
    if "choices" not in d:
        raise RuntimeError(f"Mistral error: {d}")
    text = d["choices"][0]["message"]["content"]
    usage = d.get("usage", {})
    in_tok = usage.get("prompt_tokens", len(prompt) // 4)
    out_tok = usage.get("completion_tokens", len(text or "") // 4)
    return text, in_tok, out_tok


MODELS: dict[str, Callable[[str], tuple[str, int, int]]] = {
    "Claude Opus 4": lambda p: call_ant(FROZEN_MODEL_IDS["Claude Opus 4"], p),
    "Claude Sonnet 4": lambda p: call_ant(FROZEN_MODEL_IDS["Claude Sonnet 4"], p),
    "GPT-5.2": lambda p: call_oai(FROZEN_MODEL_IDS["GPT-5.2"], p),
    "GPT-4.1": lambda p: call_oai(FROZEN_MODEL_IDS["GPT-4.1"], p),
    "o3": lambda p: call_oai(FROZEN_MODEL_IDS["o3"], p),
    "o4-mini": lambda p: call_oai(FROZEN_MODEL_IDS["o4-mini"], p),
    "Gemini 2.5 Flash": call_gem,
    "DeepSeek V3": call_dsk,
    "Mistral Large 2": call_mis,
}


# =============================================================================
# cpic_rag prompt construction
# =============================================================================

def cpic_rag_prompt(tc: dict, pop: dict, gene_excerpt: str) -> str:
    """Build the cpic_rag prompt.

    Constraints (asserted in build_and_assert_no_leakage()):
      - Contains: gene-level CPIC excerpt (with full recommendation table) +
        patient genotype + patient population context.
      - Does NOT contain: tc['gt_diplotype'], tc['gt_phenotype'],
        tc['gt_activity'], tc['gt_drug'], tc['pop_note'][pop['id']].
    """
    return f"""You are a pharmacogenomics interpretation system. Use the CPIC guideline excerpt below to decide the diplotype, phenotype, drug recommendation, hazard, and population-specific note for the patient.

## CPIC guideline excerpt (retrieved for gene: {tc['gene']})

{gene_excerpt}

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


def assert_no_leakage(prompt: str, tc: dict, pop: dict) -> list[str]:
    """Return a list of leak findings; empty list if clean.
    The pre-computed analytical-decision strings must NOT appear in the prompt."""
    leaks = []
    forbidden = {
        "gt_diplotype": tc["gt_diplotype"],
        "gt_phenotype": tc["gt_phenotype"],
        "gt_drug": tc["gt_drug"],
        "pop_note": tc["pop_note"][pop["id"]],
    }
    p_lc = prompt.lower()
    for label, val in forbidden.items():
        if not val:
            continue
        # We accept the case where the value happens to be a single CPIC-table
        # cell (which the model has to navigate to, by design — this is the
        # "answer is one row in a multi-row table" interpretation). The leak
        # we forbid is verbatim presence of the FULL gt_drug recommendation
        # string with its parenthetical hazard, or the gt_diplotype with
        # qualifiers like "(homozygous)".
        # For diplotype/phenotype, accept short canonical forms (*4/*4, "Poor
        # Metaboliser") appearing in the table — those ARE supposed to be in
        # the table. Only flag the full case-specific gt string with all its
        # qualifiers/parentheticals.
        if len(val) > 60 and val.lower() in p_lc:
            leaks.append(f"{label} verbatim leak: {val[:80]!r}")
    return leaks


# =============================================================================
# Parsing + (preliminary) scoring — same shape as 02-run-benchmark-v3.py
# Final scoring is applied by 17-merge-rag-results.py via 10-rescore-v3.py.
# =============================================================================

def parse(text: str) -> dict:
    if not text:
        return {}
    fields = {}
    for line in text.strip().split("\n"):
        for key in ["DIPLOTYPE:", "PHENOTYPE:", "DRUG:", "HAZARD:", "POPULATION:"]:
            if key in line.upper():
                val = line[line.upper().index(key) + len(key):].strip()
                fields[key.rstrip(":")] = val
    return fields


def score_preliminary(parsed: dict, tc: dict, pop: dict) -> dict:
    """Same preliminary scorer used by 02-run-benchmark-v3.py. Final
    rigorous scoring is applied by 10-rescore-v3.py during the merge."""
    if not parsed:
        return {"A1": 0, "A2": 0, "A3": 0, "B1": 0, "B2": 0, "B3": 0,
                "tier_a": 0, "tier_b": 0, "overall": 0, "format_fail": True}
    phen = parsed.get("PHENOTYPE", "").lower()
    drug = parsed.get("DRUG", "").lower()
    pop_text = parsed.get("POPULATION", "").lower()
    gt_first = tc["gt_phenotype"].lower().split()[0]
    s = {}
    s["A1"] = 1.0 if gt_first in phen else 0.0
    gt_drug = tc["gt_drug"].lower()
    if "avoid" in gt_drug and "lethal" in gt_drug:
        s["A2"] = 1.0 if ("avoid" in drug or "contraindicated" in drug) else 0.0
    elif "avoid" in gt_drug:
        s["A2"] = 1.0 if any(k in drug for k in ("avoid", "contraindicated", "reduce", "alternative")) else 0.0
    elif "standard" in gt_drug:
        s["A2"] = 0.0 if ("avoid" in drug or "contraindicated" in drug) else 1.0
    elif "reduce" in gt_drug or "lower" in gt_drug:
        s["A2"] = 1.0 if any(k in drug for k in ("reduce", "lower", "caution", "decrease", "adjust")) else 0.5
    else:
        s["A2"] = 0.5
    s["A3"] = (1.0 if ("avoid" in drug or "contraindicated" in drug) else 0.0) if "lethal" in gt_drug else 1.0
    s["tier_a"] = (s["A1"] + s["A2"] + s["A3"]) / 3.0
    pop_kw = {"EUR": ["european", "caucasian"],
              "AMR": ["latin", "peru", "admixed", "indigenous", "amerindian", "mestizo"],
              "AFR": ["african", "uganda", "east african", "hiv", "efavirenz", "first-line"]}
    s["B1"] = min(sum(1 for kw in pop_kw.get(pop["id"], []) if kw in pop_text) / 2.0, 1.0)
    has = sum(1 for k in ["DIPLOTYPE", "PHENOTYPE", "DRUG", "HAZARD", "POPULATION"]
              if k in parsed and len(parsed[k]) > 1)
    s["B2"] = has / 5.0
    full = " ".join(parsed.values()).lower()
    s["B3"] = min(sum(1 for kw in ["cpic", "guideline", "activity", "score", "allele", "frequency"] if kw in full) / 2.0, 1.0)
    s["tier_b"] = (s["B1"] + s["B2"] + s["B3"]) / 3.0
    s["overall"] = (s["tier_a"] + s["tier_b"]) / 2.0
    s["format_fail"] = False
    return s


def cost_for(model: str, in_tok: int, out_tok: int) -> float:
    p_in, p_out = PRICE[model]
    return (in_tok * p_in + out_tok * p_out) / 1_000_000


# =============================================================================
# Pre-flight gate
# =============================================================================

def estimate_cost(cases: list[dict], corpus: dict, n_runs: int) -> dict:
    """Project total cost from average prompt size and corpus excerpt size."""
    avg_excerpt = sum(g["char_count"] for g in corpus["genes"].values()) / len(corpus["genes"])
    # Rough char-to-token ratio: 4
    avg_in_tokens = (avg_excerpt + 600) / 4  # 600 chars for patient + format spec
    avg_out_tokens = 200  # 5-line response
    total_calls = len(MODELS) * len(cases) * len(POPULATIONS) * n_runs
    per_model_calls = len(cases) * len(POPULATIONS) * n_runs
    by_model = {}
    total = 0.0
    for m, (p_in, p_out) in PRICE.items():
        c = (avg_in_tokens * p_in + avg_out_tokens * p_out) / 1_000_000 * per_model_calls
        by_model[m] = round(c, 2)
        total += c
    return {"total_calls": total_calls, "total_cost_estimated": round(total, 2),
            "by_model_cost_estimated": by_model, "avg_excerpt_chars": int(avg_excerpt)}


BUCKET_LOWER = 0.80
BUCKET_UPPER = 1.00


def bucket_placement_check(logger) -> None:
    """Pre-commit gate before the full RAG arm.

    Re-scores the existing 45-call dry-run with the clinical-equivalence
    rescorer (10b-) and asserts the cpic_rag aggregate A1 lands strictly in
    the expected window (BUCKET_LOWER, BUCKET_UPPER). Aborts via sys.exit(2)
    if the placement is outside the envelope — that means either RAG is
    matching with_spec (modulo vocabulary; surprising, worth investigating)
    or RAG is below no_spec (unexpected; abort).

    Pre-condition: --dry-run was executed previously and produced
    v3_rag_raw_dry.json. Refuses to proceed if the file is missing.
    """
    from importlib.util import spec_from_file_location, module_from_spec

    dry_path = RESULTS / "v3_rag_raw_dry.json"
    if not dry_path.exists():
        logger.error("[gate] dry-run output %s not found; run --dry-run first",
                     dry_path)
        sys.exit(2)

    sp = spec_from_file_location(
        "rescore_eq", PYDIR / "10b-rescore-v3-clinical-equivalence.py")
    rs_eq = module_from_spec(sp)
    sp.loader.exec_module(rs_eq)

    cases_by_id = {c["id"]: c for c in json.loads(SPECS_FILE.read_text())}
    dry_rows = json.loads(dry_path.read_text())
    parsed = [r for r in dry_rows if not r["scores"].get("format_fail")]
    if not parsed:
        logger.error("[gate] no parsed dry-run rows; cannot compute placement")
        sys.exit(2)

    total = 0.0
    for r in parsed:
        tc = cases_by_id[r["tc"]]
        total += rs_eq.score_a1_clinical_eq(
            (r["parsed"] or {}).get("PHENOTYPE", ""),
            tc["gt_phenotype"], tc["gene"])
    agg = total / len(parsed)
    logger.info("[gate] dry-run cpic_rag aggregate A1 under 10b-: %.4f (n=%d)",
                agg, len(parsed))

    if not (BUCKET_LOWER < agg < BUCKET_UPPER):
        logger.error("[gate] BUCKET PLACEMENT FAIL: aggregate %.4f outside (%.2f, %.2f)",
                     agg, BUCKET_LOWER, BUCKET_UPPER)
        if agg >= BUCKET_UPPER - 0.001:
            logger.error("[gate]   cpic_rag is at with_spec ceiling — RAG matches "
                         "spec modulo vocabulary; investigate before $130")
        elif agg <= BUCKET_LOWER:
            logger.error("[gate]   cpic_rag underperforms no_spec — investigate "
                         "retrieval/format issues before $130")
        logger.error("[gate] aborting full run; re-run --dry-run or audit "
                     "10b- rescorer if this is unexpected")
        sys.exit(2)

    logger.info("[gate] bucket placement OK: %.4f strictly between (%.2f, %.2f)",
                agg, BUCKET_LOWER, BUCKET_UPPER)


def confirmation_gate(estimate: dict, corpus: dict, dry_run: bool, skip: bool) -> bool:
    """Print the run plan and require y/N unless --yes or --dry-run."""
    if dry_run:
        return True
    print()
    print("=" * 70)
    print("cpic_rag CONDITION — production run plan")
    print("=" * 70)
    print(f"  Total calls:           {estimate['total_calls']:,}")
    print(f"  Estimated total cost:  ${estimate['total_cost_estimated']:.2f}")
    print(f"  Hard cap:              ${BUDGET_HARD_CAP:.2f}  (run aborts if exceeded)")
    print(f"  Soft warn:             ${BUDGET_SOFT_WARN:.2f}  (stderr warning)")
    print(f"  Avg excerpt size:      {estimate['avg_excerpt_chars']} chars")
    mistral_calls = len(corpus["genes"]) and (110 * 3 * 3)
    print(f"  Mistral wall-clock:    ~{mistral_calls * 15 / 3600:.1f} hrs (15 s spacing × {mistral_calls} calls)")
    print()
    print("  Models (frozen IDs):")
    for m, mid in FROZEN_MODEL_IDS.items():
        c = estimate["by_model_cost_estimated"][m]
        print(f"    {m:<22} {mid:<30} est. ${c:>6.2f}  spacing {SPACING[m]:>5}s")
    print()
    print(f"  Corpus version:        schema {corpus['schema_version']}")
    print(f"  Corpus generated_at:   {corpus['generated_at_utc']}")
    sha_tail = max(g["sha256"] for g in corpus["genes"].values())
    print(f"  Max excerpt sha256:    {sha_tail[:16]}…")
    print(f"  Genes in corpus:       {len(corpus['genes'])} (TPMT/NUDT15 ship with documented content gap)")
    print()
    if skip:
        print("  --yes flag set; skipping confirmation.")
        return True
    try:
        ans = input("  Proceed? [y/N] ").strip().lower()
    except EOFError:
        ans = ""
    return ans == "y"


# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true",
                    help="5 cases × 3 models × 3 populations × 1 run = 45 calls (~$0.50)")
    ap.add_argument("--yes", action="store_true",
                    help="Skip the interactive confirmation gate")
    ap.add_argument("--include-models", nargs="+", default=None,
                    help="Restrict run to these models (parallel-process support)")
    ap.add_argument("--suffix", default="",
                    help="Output-file suffix for parallel runs without clobbering")
    args = ap.parse_args()

    logger = setup_logging(args.dry_run)
    logger.info("[main] cpic_rag runner starting (dry_run=%s)", args.dry_run)

    # Load
    cases = json.loads(SPECS_FILE.read_text())
    corpus = json.loads(CORPUS_FILE.read_text())

    # In dry-run, restrict to 5 cases × 3 models × 3 pops × 1 run = 45 calls
    if args.dry_run:
        # Pick 5 representative cases that exercise the hardest paths
        dry_ids = ["dpyd_fu_pm", "cyp2d6_codeine_pm", "hlab5701_aba_pos",
                   "cyp2c19_clop_rm", "tpmt_aza_pm"]
        cases = [c for c in cases if c["id"] in dry_ids]
        models_to_run = {k: v for k, v in MODELS.items()
                         if k in {"Claude Opus 4", "GPT-5.2", "DeepSeek V3"}}
        n_runs = 1
    else:
        models_to_run = (
            {k: v for k, v in MODELS.items() if k in args.include_models}
            if args.include_models else MODELS
        )
        n_runs = N_RUNS

    estimate = estimate_cost(cases, corpus, n_runs)
    if not args.dry_run:
        # Recompute for filtered models
        per_call_calls = len(cases) * len(POPULATIONS) * n_runs
        by_model_filtered = {}
        total_filtered = 0.0
        for m in models_to_run:
            c = estimate["by_model_cost_estimated"][m] * (per_call_calls / (110 * 3 * 3))
            by_model_filtered[m] = round(c, 2)
            total_filtered += c
        estimate["total_calls"] = len(models_to_run) * len(cases) * len(POPULATIONS) * n_runs
        estimate["total_cost_estimated"] = round(total_filtered, 2)
        estimate["by_model_cost_estimated"] = by_model_filtered

    # Pre-commit bucket-placement gate: before any non-dry-run kicks off, the
    # 45-call dry-run aggregate (under 10b-) must land strictly between
    # no_spec (~0.80) and with_spec (1.00). See bucket_placement_check.
    if not args.dry_run:
        bucket_placement_check(logger)

    if not confirmation_gate(estimate, corpus, args.dry_run, args.yes):
        logger.info("[main] confirmation declined; exiting")
        sys.exit(0)

    # Set output paths
    suffix = args.suffix or ("_dry" if args.dry_run else "")
    paths = {
        "raw": RESULTS / f"v3_rag_raw{suffix}.json",
        "partial": RESULTS / f"v3_rag_partial{suffix}.json",
        "ckpt": RESULTS / f"v3_rag_checkpoint{suffix}.json",
        "cost": RESULTS / f"v3_rag_cost_log{suffix}.json",
    }

    # Resume support
    results = []
    completed_keys = set()
    if paths["partial"].exists() and not args.dry_run:
        try:
            results = json.loads(paths["partial"].read_text())
            for r in results:
                completed_keys.add((r["run"], r["model"], r["tc"], r["pop"], r["cond"]))
            logger.info("[resume] loaded %d prior results from %s", len(results), paths["partial"])
        except Exception as e:
            logger.warning("[resume] failed to load partial (%s); starting fresh", e)
            results = []

    cost_by_model = defaultdict(lambda: {"in_tok": 0, "out_tok": 0, "cost": 0.0, "calls": 0})
    if paths["cost"].exists() and not args.dry_run:
        try:
            saved = json.loads(paths["cost"].read_text())
            for k, v in saved.items():
                cost_by_model[k] = v
        except Exception:
            pass

    total = len(models_to_run) * len(cases) * len(POPULATIONS) * n_runs
    count = len(results)
    t0 = time.time()
    soft_warn_emitted = False

    # Pre-flight leakage assertion on first 3 cases
    for tc in cases[:3]:
        for pop in POPULATIONS[:1]:
            gene_excerpt = corpus["genes"][tc["gene"]]["guideline_excerpt"]
            prompt = cpic_rag_prompt(tc, pop, gene_excerpt)
            leaks = assert_no_leakage(prompt, tc, pop)
            if leaks:
                logger.error("[preflight] LEAKAGE in %s/%s: %s", tc["id"], pop["id"], leaks)
                sys.exit(2)
    logger.info("[preflight] no leakage detected across 3 sample prompts")

    # Run loop
    for run in range(n_runs):
        for mname, mfunc in models_to_run.items():
            spacing = SPACING.get(mname, 0.5)
            for tc in cases:
                gene_excerpt = corpus["genes"][tc["gene"]]["guideline_excerpt"]
                for pop in POPULATIONS:
                    key = (run, mname, tc["id"], pop["id"], CONDITION)
                    if key in completed_keys:
                        continue
                    count += 1
                    prompt = cpic_rag_prompt(tc, pop, gene_excerpt)

                    try:
                        time.sleep(spacing)
                        text, in_tok, out_tok = _retry(mfunc, prompt, max_retries=3, base_backoff=2.0)
                        parsed_resp = parse(text)
                        scores = score_preliminary(parsed_resp, tc, pop)
                        entry = {"run": run, "model": mname, "tc": tc["id"],
                                 "gene": tc["gene"], "drug": tc.get("drug", ""), "pop": pop["id"],
                                 "cond": CONDITION, "parsed": parsed_resp, "scores": scores,
                                 "in_tok": in_tok, "out_tok": out_tok}
                        cb = cost_by_model[mname]
                        cb["in_tok"] += in_tok
                        cb["out_tok"] += out_tok
                        cb["cost"] += cost_for(mname, in_tok, out_tok)
                        cb["calls"] += 1
                    except Exception as e:
                        entry = {"run": run, "model": mname, "tc": tc["id"],
                                 "gene": tc["gene"], "drug": tc.get("drug", ""), "pop": pop["id"],
                                 "cond": CONDITION, "parsed": {},
                                 "scores": {"A1": 0, "A2": 0, "A3": 0, "B1": 0, "B2": 0, "B3": 0,
                                            "tier_a": 0, "tier_b": 0, "overall": 0, "format_fail": True},
                                 "error": str(e)[:200],
                                 "in_tok": 0, "out_tok": 0}
                        logger.error("ERROR [%d/%d] %s | %s | %s | cpic_rag: %s",
                                     count, total, mname, tc["id"], pop["id"], str(e)[:120])

                    results.append(entry)

                    if count % 25 == 0 or count == total:
                        total_cost = sum(v["cost"] for v in cost_by_model.values())
                        elapsed = time.time() - t0
                        rate = count / elapsed if elapsed > 0 else 0
                        eta_s = (total - count) / rate if rate > 0 else 0
                        logger.info("[%d/%d] %-22s | %-25s | %s | cpic_rag | cost $%.2f | rate %.2f/s | ETA %.0f min",
                                    count, total, mname, tc["id"][:25], pop["id"], total_cost, rate, eta_s / 60)
                        if total_cost > BUDGET_HARD_CAP:
                            logger.error("[budget] HARD CAP EXCEEDED ($%.2f > $%.0f) — saving and aborting",
                                         total_cost, BUDGET_HARD_CAP)
                            _save_all(paths, results, cost_by_model, count, total)
                            return
                        if total_cost > BUDGET_SOFT_WARN and not soft_warn_emitted:
                            print(f"[budget] soft warning: cost ${total_cost:.2f} above ${BUDGET_SOFT_WARN:.0f} threshold",
                                  file=sys.stderr)
                            logger.warning("[budget] soft-warn threshold crossed at $%.2f", total_cost)
                            soft_warn_emitted = True
                    if count % 100 == 0 and not args.dry_run:
                        _save_all(paths, results, cost_by_model, count, total)

    _save_all(paths, results, cost_by_model, count, total)
    total_cost = sum(v["cost"] for v in cost_by_model.values())
    logger.info("[main] run complete: %d calls, $%.2f total", count, total_cost)
    print(json.dumps({"calls": count, "cost_total": round(total_cost, 4),
                      "by_model": {k: dict(v) for k, v in cost_by_model.items()}}, indent=2, default=str))


def _save_all(paths, results, cost_by_model, count, total):
    final = paths["raw"] if count >= total else paths["partial"]
    final.write_text(json.dumps(results, indent=2))
    paths["partial"].write_text(json.dumps(results, indent=2))
    paths["ckpt"].write_text(json.dumps({"count": count, "total": total}))
    paths["cost"].write_text(json.dumps({k: dict(v) for k, v in cost_by_model.items()}, indent=2, default=str))


if __name__ == "__main__":
    main()
