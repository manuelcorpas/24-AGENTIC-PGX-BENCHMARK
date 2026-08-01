#!/usr/bin/env python3
"""Freeze the seven-model input-normalisation analysis and its uncertainty.

The raw responses remain in their original source files. This script defines
the exact confirmatory cohort, hashes every source, checks completeness and
produces one small canonical analysis file used by Figure 8, Table S9, the
manuscript validator and the Zenodo manifest.

No API calls are made. Provider errors and output-budget truncations are
reported but excluded from behavioural denominators. With/without comparisons
use the intersection of (sample, gene) units available for that model; o3 is
therefore paired on its prespecified 150-pair no-definition subsample.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from datetime import date
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
OUT = DATA / "v3_input_normalisation_seven_model_freeze.json"

MODELS = [
    "Claude Opus 4.5",
    "Claude Sonnet 4.5",
    "GPT-5.2",
    "GPT-4.1",
    "o3",
    "o4-mini",
    "DeepSeek V3",
]

DEFINITION_FILES = {
    "Claude Opus 4.5": ["v3_input_normalisation_defs.json",
                         "v3_input_normalisation_defs_tail.json"],
    "Claude Sonnet 4.5": ["v3_input_normalisation_defs_sonnet_gpt41_o4mini.json",
                           "v3_input_normalisation_defs_sonnet_gpt41_o4mini_rest.json"],
    "GPT-5.2": ["v3_input_normalisation_defs_gpt52.json"],
    "GPT-4.1": ["v3_input_normalisation_defs_sonnet_gpt41_o4mini.json",
                "v3_input_normalisation_defs_sonnet_gpt41_o4mini_rest.json"],
    "o3": ["v3_input_normalisation_defs_o3.json"],
    "o4-mini": ["v3_input_normalisation_defs_sonnet_gpt41_o4mini.json",
                 "v3_input_normalisation_defs_sonnet_gpt41_o4mini_rest.json"],
    "DeepSeek V3": ["v3_input_normalisation_defs_deepseek.json"],
}

NO_DEFINITION_FILES = {
    "Claude Opus 4.5": ["v3_input_normalisation_main.json"],
    "Claude Sonnet 4.5": ["v3_input_normalisation_nodefs_four.json"],
    "GPT-5.2": ["v3_input_normalisation_main.json"],
    "GPT-4.1": ["v3_input_normalisation_nodefs_four.json"],
    "o3": ["v3_input_normalisation_o3.json"],
    "o4-mini": ["v3_input_normalisation_nodefs_four.json"],
    "DeepSeek V3": ["v3_input_normalisation_nodefs_four.json"],
}

EXPECTED_NO_DEFINITION = {m: 527 for m in MODELS}
EXPECTED_NO_DEFINITION["o3"] = 150
BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_SEED = 20260801


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_blob(name: str) -> dict:
    path = DATA / name
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def model_rows(files: list[str], model: str) -> list[dict]:
    rows = []
    for name in files:
        rows.extend(r for r in load_blob(name)["rows"] if r["model"] == model)
    return rows


def nd(value: str) -> tuple[str, ...]:
    return tuple(sorted(part.strip() for part in value.split("/")))


def operational(rows: list[dict]) -> dict:
    counts = {k: 0 for k in ("call", "abstain", "error", "truncated_output")}
    for row in rows:
        status = row.get("status")
        if status not in counts:
            raise ValueError(f"unexpected status {status!r}")
        counts[status] += 1
    return {
        "attempted": len(rows),
        **counts,
        "scored": counts["call"] + counts["abstain"],
    }


def percentile(values: list[float], p: float) -> float:
    values = sorted(values)
    pos = (len(values) - 1) * p
    lo = int(pos)
    hi = min(lo + 1, len(values) - 1)
    frac = pos - lo
    return values[lo] * (1 - frac) + values[hi] * frac


def cluster_interval(rows: list[dict], numerator) -> dict:
    """Sample-cluster bootstrap interval for a binary row-level proportion."""
    if not rows:
        return {"estimate": None, "ci95": [None, None], "clusters": 0}
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["sample"]].append(row)
    clusters = sorted(grouped)
    point_num = sum(numerator(row) for row in rows)
    point = point_num / len(rows)
    cluster_den = np.array([len(grouped[c]) for c in clusters], dtype=float)
    cluster_num = np.array(
        [sum(numerator(row) for row in grouped[c]) for c in clusters], dtype=float
    )
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    picks = rng.integers(0, len(clusters),
                         size=(BOOTSTRAP_REPLICATES, len(clusters)))
    draws = (cluster_num[picks].sum(axis=1) /
             cluster_den[picks].sum(axis=1)).tolist()
    return {
        "estimate": point,
        "ci95": [percentile(draws, 0.025), percentile(draws, 0.975)],
        "clusters": len(clusters),
    }


def cluster_difference(rows: list[dict]) -> dict:
    """Paired model-minus-caller accuracy against GeT-RM on emitted calls."""
    if not rows:
        return {"estimate": None, "ci95": [None, None], "clusters": 0}
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["sample"]].append(row)
    clusters = sorted(grouped)

    cluster_den = np.array([len(grouped[c]) for c in clusters], dtype=float)
    cluster_diff = np.array([
        sum((nd(row["call"]) == nd(row["getrm"])) -
            (nd(row["reference"]) == nd(row["getrm"]))
            for row in grouped[c])
        for c in clusters
    ], dtype=float)
    point = cluster_diff.sum() / cluster_den.sum()
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    picks = rng.integers(0, len(clusters),
                         size=(BOOTSTRAP_REPLICATES, len(clusters)))
    draws = (cluster_diff[picks].sum(axis=1) /
             cluster_den[picks].sum(axis=1)).tolist()
    return {
        "estimate": point,
        "ci95": [percentile(draws, 0.025), percentile(draws, 0.975)],
        "clusters": len(clusters),
    }


def metrics(rows: list[dict]) -> dict:
    op = operational(rows)
    scored = [r for r in rows if r["status"] in ("call", "abstain")]
    calls = [r for r in scored if r["status"] == "call"]
    coverage = cluster_interval(scored, lambda r: r["status"] == "call")
    pypgx = cluster_interval(calls, lambda r: nd(r["call"]) == nd(r["reference"]))
    getrm = cluster_interval(calls, lambda r: nd(r["call"]) == nd(r["getrm"]))
    caller = cluster_interval(calls, lambda r: nd(r["reference"]) == nd(r["getrm"]))
    return {
        "operational": op,
        "coverage": coverage,
        "accuracy_vs_pypgx_among_calls": pypgx,
        "accuracy_vs_getrm_among_calls": getrm,
        "caller_vs_getrm_on_model_calls": caller,
        "model_minus_caller_vs_getrm": cluster_difference(calls),
    }


def main() -> int:
    source_names = sorted({f for fs in DEFINITION_FILES.values() for f in fs} |
                          {f for fs in NO_DEFINITION_FILES.values() for f in fs})
    blobs = {name: load_blob(name) for name in source_names}
    input_hashes = {blob.get("inputs_sha256") for blob in blobs.values()}
    if len(input_hashes) != 1:
        raise AssertionError(f"input hashes disagree: {input_hashes}")

    result = {
        "freeze_version": "2026-08-01-c15",
        "frozen_on": str(date.today()),
        "models": MODELS,
        "gemini_policy": (
            "Gemini 2.5 Flash is not assigned a performance estimate: a complete "
            "confirmatory run was not collected. Its boxed-LaTeX pilot response is "
            "accepted by the same presentation-neutral parser used for every model; "
            "no model is excluded because of response wrapping."
        ),
        "parser_policy": (
            "Use the last explicit DIPLOTYPE marker; ignore Markdown/LaTeX presentation "
            "wrappers; strip an optional gene prefix; validate the first value token as "
            "one complete schema-supported diplotype; never search surrounding reasoning."
        ),
        "inputs_sha256": next(iter(input_hashes)),
        "bootstrap": {
            "unit": "sample cluster",
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": BOOTSTRAP_SEED,
            "interval": "percentile 95%",
        },
        "source_files": {
            name: {"sha256": sha256(DATA / name), "rows": len(blobs[name]["rows"])}
            for name in source_names
        },
        "models_analysis": {},
    }

    for model in MODELS:
        with_rows = model_rows(DEFINITION_FILES[model], model)
        without_rows = [r for r in model_rows(NO_DEFINITION_FILES[model], model)
                        if r["form"] == "vcf"]
        if len(with_rows) != 527:
            raise AssertionError(f"{model}: expected 527 definition rows, got {len(with_rows)}")
        expected = EXPECTED_NO_DEFINITION[model]
        if len(without_rows) != expected:
            raise AssertionError(f"{model}: expected {expected} no-definition rows, got {len(without_rows)}")
        for label, rows in (("with", with_rows), ("without", without_rows)):
            keys = [(r["sample"], r["gene"]) for r in rows]
            if len(keys) != len(set(keys)):
                raise AssertionError(f"{model} {label}: duplicate sample-gene rows")
        common = ({(r["sample"], r["gene"]) for r in with_rows} &
                  {(r["sample"], r["gene"]) for r in without_rows})
        if len(common) != expected:
            raise AssertionError(f"{model}: paired key count {len(common)} != {expected}")
        paired_with = [r for r in with_rows if (r["sample"], r["gene"]) in common]
        paired_without = [r for r in without_rows if (r["sample"], r["gene"]) in common]
        result["models_analysis"][model] = {
            "paired_units": len(common),
            "without_definitions": metrics(paired_without),
            "with_definitions": metrics(paired_with),
            "with_definitions_full_527": metrics(with_rows),
        }

    all_with = []
    for model in MODELS:
        all_with.extend(model_rows(DEFINITION_FILES[model], model))
    result["definition_arm_operational_total"] = operational(all_with)
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {OUT}")
    print(f"input SHA-256 {result['inputs_sha256']}")
    print(f"definition arm {result['definition_arm_operational_total']}")
    for model in MODELS:
        m = result["models_analysis"][model]
        wo = m["without_definitions"]
        wi = m["with_definitions"]
        diff = wi["model_minus_caller_vs_getrm"]
        print(
            f"{model:20s} n={m['paired_units']:3d} "
            f"coverage {wo['coverage']['estimate']:.3f}->{wi['coverage']['estimate']:.3f} "
            f"GeT-RM {wi['accuracy_vs_getrm_among_calls']['estimate']:.3f} "
            f"caller {wi['caller_vs_getrm_on_model_calls']['estimate']:.3f} "
            f"diff {diff['estimate']:+.3f} "
            f"CI [{diff['ci95'][0]:+.3f}, {diff['ci95'][1]:+.3f}]"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
