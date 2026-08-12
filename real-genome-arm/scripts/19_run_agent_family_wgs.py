#!/usr/bin/env python3
"""
Step 19 — run the model panel over a table of observed (gene, diplotype) states.

Supersedes 04_run_agent_realgenome.py for the family arm. Three differences,
each of them a defect this project has already paid for:

  1. PROVENANCE. Every row records in/out tokens, computed cost, a UTC timestamp
     and the provider's own response id. Script 04 recorded none of these, which
     is why the reproducibility answer to Reviewer 2 could only be half given:
     the raw rows carried no timestamps or response ids and they were not
     reconstructable after the fact.
  2. FAIL LOUDLY. A model named in the panel file but not implemented raises
     before any call is made, and the script exits non-zero if any call errored
     or if the row count does not equal states x models. A short run can no
     longer terminate as a success (CORRECTIONS.md C16).
  3. NO SILENT ABSTENTION. A provider error is recorded as an error, never as a
     model declining to answer. An HTTP 429 read as an abstention is C12, the
     sharpest entry in the corrections log.

Usage:
  19_run_agent_family_wgs.py <states.tsv> <out.json> [--models PATH] [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import urllib.error
import urllib.request

HERE = Path(__file__).resolve().parent
ARM = HERE.parent

# $ per million tokens (input, output).
# These are not list prices recalled from memory: they were recovered by least
# squares from the 13,199 rows of data/v3_five_cell_live.json, which record
# in_tokens, out_tokens and cost_usd per call. The fit is exact to 1e-6, so
# these are the rates this project was actually charged for that run.
RATES = {
    "Claude Opus 4.5": (5.00, 25.00),
    "Claude Sonnet 4.5": (3.00, 15.00),
    "GPT-5.2": (1.25, 10.00),
    "GPT-4.1": (2.00, 8.00),
    "o3": (2.00, 8.00),
    "o4-mini": (1.10, 4.40),
    "Gemini 2.5 Flash": (0.30, 2.50),
    "DeepSeek V3": (0.27, 1.10),
}

PROVIDERS = {"anthropic", "openai", "openai-reasoning", "google", "deepseek"}


def build_prompt(gene: str, diplotype: str) -> str:
    """Identical wording to 04_run_agent_realgenome.py, so predictions from the
    two scripts are comparable and existing rows remain valid."""
    return (
        "You are executing a ClawBio pharmacogenomics skill. Apply CPIC Level A guidelines. "
        "Given the gene and the patient's star-allele diplotype, output the metaboliser phenotype "
        "(or CPIC functional status). If the diplotype contains an allele of uncertain or unknown "
        "function, output 'Indeterminate'.\n\n"
        f"Gene: {gene}\nDiplotype: {diplotype}\n\n"
        "Output one line:\nPHENOTYPE: [CPIC phenotype]"
    )


def parse_phenotype(text: str | None) -> str:
    for line in (text or "").split("\n"):
        if "PHENOTYPE:" in line.upper():
            return line.split(":", 1)[1].strip()
    return (text or "").strip()[:80]


def cost_usd(model: str, in_tok: int, out_tok: int) -> float:
    ri, ro = RATES[model]
    return round(in_tok * ri / 1e6 + out_tok * ro / 1e6, 8)


def load_panel(path: Path):
    """[(display_name, model_id, provider)]. Raises on an unknown provider or an
    unpriced model BEFORE any call is made."""
    panel = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p for p in line.split("  ") if p.strip()]
        if len(parts) < 3:
            continue
        name, mid, prov = (p.strip() for p in parts[:3])
        if prov not in PROVIDERS:
            raise SystemExit(f"unknown provider {prov!r} for model {name!r}")
        if name not in RATES:
            raise SystemExit(f"no price recorded for model {name!r}; refusing to run")
        panel.append((name, mid, prov))
    if not panel:
        raise SystemExit(f"no models parsed from {path}")
    return panel


def _post(url: str, payload: dict, headers: dict, timeout: int = 180):
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", **headers},
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def call_model(model_id: str, provider: str, prompt: str):
    """(text, in_tokens, out_tokens, response_id). Raises on any provider error."""
    if provider == "anthropic":
        b = _post(
            "https://api.anthropic.com/v1/messages",
            {"model": model_id, "max_tokens": 120,
             "messages": [{"role": "user", "content": prompt}]},
            {"x-api-key": os.environ["ANTHROPIC_API_KEY"],
             "anthropic-version": "2023-06-01"},
        )
        u = b.get("usage", {})
        return (b["content"][0]["text"], u.get("input_tokens", 0),
                u.get("output_tokens", 0), b.get("id"))

    if provider in ("openai", "openai-reasoning"):
        payload = {"model": model_id,
                   "messages": [{"role": "user", "content": prompt}]}
        # Reasoning models bill their hidden reasoning as output tokens and
        # truncate at a low cap, so they get the same 2000 budget as before.
        if provider == "openai-reasoning":
            payload["max_completion_tokens"] = 2000
        else:
            payload["max_tokens"] = 120
        b = _post("https://api.openai.com/v1/chat/completions", payload,
                  {"Authorization": "Bearer " + os.environ["OPENAI_API_KEY"]})
        u = b.get("usage", {})
        return (b["choices"][0]["message"]["content"], u.get("prompt_tokens", 0),
                u.get("completion_tokens", 0), b.get("id"))

    if provider == "deepseek":
        b = _post(
            "https://api.deepseek.com/chat/completions",
            {"model": model_id, "max_tokens": 120,
             "messages": [{"role": "user", "content": prompt}]},
            {"Authorization": "Bearer " + os.environ["DEEPSEEK_API_KEY"]},
        )
        u = b.get("usage", {})
        return (b["choices"][0]["message"]["content"], u.get("prompt_tokens", 0),
                u.get("completion_tokens", 0), b.get("id"))

    if provider == "google":
        key = os.environ["GOOGLE_API_KEY"]
        b = _post(
            f"https://generativelanguage.googleapis.com/v1beta/models/{model_id}:generateContent?key={key}",
            {"contents": [{"parts": [{"text": prompt}]}],
             "generationConfig": {"maxOutputTokens": 2048}},
            {},
        )
        parts = b["candidates"][0]["content"]["parts"]
        u = b.get("usageMetadata", {})
        return ("".join(p.get("text", "") for p in parts),
                u.get("promptTokenCount", 0), u.get("candidatesTokenCount", 0),
                b.get("responseId"))

    raise ValueError(provider)


def read_states(path: Path):
    seen, states = set(), []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["cohort"], r["gene"], r["diplotype"])
            if key in seen:
                continue
            seen.add(key)
            states.append(
                {"cohort": r["cohort"], "gene": r["gene"],
                 "diplotype": r["diplotype"], "caller_phenotype": r["phenotype"]}
            )
    return states


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("states_tsv")
    ap.add_argument("out_json")
    ap.add_argument("--models", default=str(ARM / "config" / "models.txt"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    states = read_states(Path(args.states_tsv))
    if args.limit:
        states = states[: args.limit]
    panel = load_panel(Path(args.models))
    expected = len(states) * len(panel)
    print(f"{len(states)} states x {len(panel)} models = {expected} calls", flush=True)

    rows, errors = [], []

    def one(job):
        name, mid, prov, st = job
        prompt = build_prompt(st["gene"], st["diplotype"])
        started = datetime.now(timezone.utc).isoformat()
        for attempt in range(3):
            try:
                text, tin, tout, rid = call_model(mid, prov, prompt)
                return {
                    **st, "model": name, "model_id": mid, "provider": prov,
                    "raw": text, "pred": parse_phenotype(text),
                    "in_tokens": tin, "out_tokens": tout,
                    "cost_usd": cost_usd(name, tin, tout),
                    "requested_at_utc": started, "response_id": rid,
                    "attempt": attempt + 1, "error": None,
                }
            except Exception as exc:  # noqa: BLE001
                if attempt == 2:
                    return {
                        **st, "model": name, "model_id": mid, "provider": prov,
                        "raw": None, "pred": None, "in_tokens": 0, "out_tokens": 0,
                        "cost_usd": 0.0, "requested_at_utc": started,
                        "response_id": None, "attempt": attempt + 1,
                        "error": f"{type(exc).__name__}: {exc}"[:300],
                    }
                time.sleep(2 * (attempt + 1))

    jobs = [(n, m, p, s) for (n, m, p) in panel for s in states]
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for row in ex.map(one, jobs):
            rows.append(row)
            if row["error"]:
                errors.append(row)
                print(f"  ERROR {row['model']} {row['gene']} {row['diplotype']}: "
                      f"{row['error']}", flush=True)

    total = sum(r["cost_usd"] for r in rows)
    Path(args.out_json).write_text(json.dumps(rows, indent=1))
    print(f"wrote {len(rows)} rows, total cost ${total:.4f} -> {args.out_json}")

    # A short or partly failed run must not exit 0.
    if len(rows) != expected:
        print(f"FAIL: {len(rows)} rows but expected {expected}", file=sys.stderr)
        return 1
    if errors:
        print(f"FAIL: {len(errors)} calls errored", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
