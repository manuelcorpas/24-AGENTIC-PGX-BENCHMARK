#!/usr/bin/env python3
"""
Extract a structured rule table from CPIC guideline prose (reviewer point 2).

WHY THIS EXISTS
The manuscript claimed that extracting decision logic from a guideline reaches
the same place as a hand-authored rule table. It did not test that. In the
matched factorial BOTH execution cells are handed the authored controlled
vocabulary and both call the same make_executor(dip2phen, rec), so the retrieved
prose only changes what sits in the context window while the model picks a
diplotype. Their agreement was guaranteed by construction and carried no
information about extraction.

This script performs the experiment the reviewer asked for. A model reads only
the guideline prose and emits a rule table: diplotype -> phenotype, and
(drug, diplotype) -> recommendation. Nothing from the benchmark's ground truth
enters the prompt. The output is frozen with a SHA-256 so the table that gets
diffed and executed is provably the table that was extracted.

WHAT IT DELIBERATELY DOES NOT DO
It does not see gt_diplotype, gt_phenotype or gt_drug, and it is not shown the
controlled vocabulary. Supplying either would reintroduce exactly the circularity
that made the original claim vacuous.

USAGE
    python code/67-extract-rules-from-prose.py --estimate
    python code/67-extract-rules-from-prose.py --max-spend 10 --models "GPT-5.2" "Claude Opus 4.5" "o3"
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
CORPUS_FILE = BASE / "specs" / "cpic_rag_corpus_v3.json"
OUT_DEFAULT = BASE / "data" / "v3_extracted_rules.json"

DEFAULT_MODELS = ("GPT-5.2", "Claude Opus 4.5", "o3")


def _load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


runner = _load(CODE / "60-matched-factorial.py", "matched_factorial")

PROMPT = """You are converting a clinical pharmacogenomics guideline into a machine-executable rule table.

Below is the CPIC guideline text for {gene}. Read ONLY this text.

{excerpt}

Produce a JSON object with exactly two keys:

"diplotype_to_phenotype": an object mapping each diplotype named in the guideline
  to its metaboliser or function phenotype, using the guideline's own wording.
  Example: {{"*1/*1": "Normal Metabolizer", "*4/*4": "Poor Metabolizer"}}

"recommendations": an array of objects, one per (drug, diplotype) pair for which
  the guideline states a therapeutic recommendation, each with keys
  "drug", "diplotype" and "recommendation". Give the recommendation as the
  guideline states it, including any explicit avoid or contraindication wording.

Rules:
- Use only what the text supports. Do not add diplotypes or drugs it does not mention.
- Use the diplotype notation the guideline uses.
- Output the JSON object and nothing else. No commentary, no code fences.
"""


def extract_json(text: str) -> dict | None:
    """Pull the JSON object out of a model reply, tolerating code fences."""
    if not text:
        return None
    t = text.strip()
    t = re.sub(r"^```(?:json)?\s*", "", t)
    t = re.sub(r"\s*```$", "", t)
    start = t.find("{")
    if start < 0:
        return None
    depth = 0
    for i, ch in enumerate(t[start:], start):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(t[start:i + 1])
                except json.JSONDecodeError:
                    return None
    return None


def normalise(obj: dict | None) -> dict:
    """Coerce a parsed reply into the frozen schema, dropping anything malformed.

    A malformed entry is dropped rather than guessed at: a rule table assembled
    partly from repair heuristics is not the table the model produced.
    """
    if not isinstance(obj, dict):
        return {"diplotype_to_phenotype": {}, "recommendations": []}
    d2p = obj.get("diplotype_to_phenotype") or {}
    if not isinstance(d2p, dict):
        d2p = {}
    d2p = {str(k).strip(): str(v).strip() for k, v in d2p.items()
           if isinstance(k, str) and isinstance(v, (str, int, float))}
    recs = obj.get("recommendations") or []
    out = []
    if isinstance(recs, list):
        for r in recs:
            if not isinstance(r, dict):
                continue
            drug, dip, rec = r.get("drug"), r.get("diplotype"), r.get("recommendation")
            if all(isinstance(x, str) and x.strip() for x in (drug, dip, rec)):
                out.append({"drug": drug.strip(), "diplotype": dip.strip(),
                            "recommendation": rec.strip()})
    return {"diplotype_to_phenotype": d2p, "recommendations": out}


def freeze(table: dict) -> str:
    return hashlib.sha256(
        json.dumps(table, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build_extractors() -> dict:
    """Model callers with an output budget large enough for a rule table.

    The matched-factorial runner caps output at 320 tokens, which is right for a
    phenotype-and-recommendation answer and far too small here: a truncated table
    would look like a model that omitted rules, which is precisely the kind of
    harness artefact this project keeps catching.
    """
    import os
    import anthropic
    import openai

    def key(name):
        v = os.environ.get(name)
        if not v:
            raise SystemExit(f"{name} not set; export keys before running")
        return v

    ant = anthropic.Anthropic(api_key=key("ANTHROPIC_API_KEY"))
    oai = openai.OpenAI(api_key=key("OPENAI_API_KEY"))
    sem = {"anthropic": threading.Semaphore(4), "openai": threading.Semaphore(4)}

    def _ant(model, p):
        with sem["anthropic"]:
            r = ant.messages.create(model=model, max_tokens=8000,
                                    messages=[{"role": "user", "content": p}])
            return r.content[0].text, r.usage.input_tokens, r.usage.output_tokens

    def _oai(model, p, reasoning=False):
        with sem["openai"]:
            kw = ({"max_completion_tokens": 24000} if reasoning
                  else {"max_tokens": 8000})
            r = oai.chat.completions.create(
                model=model, messages=[{"role": "user", "content": p}], **kw)
            return (r.choices[0].message.content,
                    r.usage.prompt_tokens, r.usage.completion_tokens)

    return {
        "Claude Opus 4.5":   lambda p: _ant("claude-opus-4-5-20251101", p),
        "Claude Sonnet 4.5": lambda p: _ant("claude-sonnet-4-5-20250929", p),
        "GPT-5.2":           lambda p: _oai("gpt-5.2", p, reasoning=True),
        "GPT-4.1":           lambda p: _oai("gpt-4.1", p),
        "o3":                lambda p: _oai("o3", p, reasoning=True),
        "o4-mini":           lambda p: _oai("o4-mini", p, reasoning=True),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--models", nargs="*", default=list(DEFAULT_MODELS))
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--max-spend", type=float, default=None, metavar="USD")
    ap.add_argument("--estimate", action="store_true")
    ap.add_argument("--pace", type=float, default=0.0)
    a = ap.parse_args(argv)

    corpus = json.loads(CORPUS_FILE.read_text())["genes"]
    genes = a.genes or sorted(corpus)
    jobs = [(g, m) for g in genes for m in a.models]

    if a.estimate:
        chars = sum(len(corpus[g]["guideline_excerpt"]) for g in genes)
        print(f"{len(genes)} genes x {len(a.models)} models = {len(jobs)} calls")
        print(f"~{chars // 4:,} input tokens per model pass, "
              f"~{(chars // 4) * len(a.models):,} total in")
        return 0

    fns = build_extractors()
    spend = runner.Spend(a.max_spend) if a.max_spend else None
    lock = threading.Lock()
    results: list[dict] = []

    def one(job):
        gene, model = job
        fn = fns.get(model)
        if fn is None:
            return {"gene": gene, "model": model, "error": f"unknown model {model}"}
        prompt = PROMPT.format(gene=gene, excerpt=corpus[gene]["guideline_excerpt"])
        if spend is not None:
            spend.check()
        try:
            text, in_tok, out_tok = fn(prompt)
            if spend is not None:
                spend.add(model, in_tok, out_tok)
        except Exception as exc:                                   # noqa: BLE001
            return {"gene": gene, "model": model, "error": str(exc)}
        parsed_obj = extract_json(text)
        table = normalise(parsed_obj)
        # A truncated reply yields unbalanced JSON. Recording it as an empty table
        # would read as "the model extracted no rules", which is a different and
        # much more interesting claim than "we cut it off".
        truncated = bool(text) and parsed_obj is None
        row = {
            "gene": gene, "model": model,
            "table": table, "sha256": freeze(table),
            "n_diplotypes": len(table["diplotype_to_phenotype"]),
            "n_recommendations": len(table["recommendations"]),
            "parsed": bool(table["diplotype_to_phenotype"] or table["recommendations"]),
            "truncated_or_unparseable": truncated,
            "in_tokens": in_tok, "out_tokens": out_tok,
            "cost_usd": round(runner.Spend.cost(model, in_tok, out_tok), 6),
            "raw": text,
        }
        with lock:
            results.append(row)
            print(f"  {gene:12s} {model:18s} "
                  f"{row['n_diplotypes']:3d} diplotypes, "
                  f"{row['n_recommendations']:3d} recommendations"
                  + ("  TRUNCATED/UNPARSEABLE" if truncated else ""))
        return row

    print(f"extracting rule tables: {len(jobs)} calls"
          + (f", hard cap ${a.max_spend:.2f}" if a.max_spend else ""))
    with ThreadPoolExecutor(max_workers=6) as ex:
        list(ex.map(one, jobs))

    errs = [r for r in results if r.get("error")]
    total = sum(r.get("cost_usd") or 0 for r in results)
    a.out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {a.out}  ({len(results)} tables, {len(errs)} errors, ${total:.2f})")
    if errs:
        for r in errs[:5]:
            print(f"  ERROR {r['gene']} {r['model']}: {r['error'][:100]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
