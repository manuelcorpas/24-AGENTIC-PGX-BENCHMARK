#!/usr/bin/env python3
"""Build the replication runs: 2 cells x 3 replicates of a hashed blind package.

Cells (the two that are worth the wall-clock, see NEXT-THREAD-BRIEF.md):
  api_corpus   CPIC API allowed plus the corpus. Tests whether 97.3% holds.
  corpus_only  network denied, corpus only. The like-for-like contrast with the
               paper's RAG arm, which never had API access.

Every replicate gets its own folder with hashed case ids and a corpus stripped
of the metadata that revealed the package was built around this case set. The
surrogate-to-real id map is written to _KEYS/, outside every handed folder.

Usage: python PYTHON/02-build_runs.py [--reps 3] [--out RUNS]
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]

# Per-gene keys that describe the case set or the evaluation setup rather than
# the guideline. n_test_cases and verification were the leak Claude Science
# found on 2026-09-15; verbatim_fraction_estimate reveals that the corpus was
# scored for an evaluation, which is setup the agent should not see either.
CASE_SET_METADATA = ("n_test_cases", "verification", "verbatim_fraction_estimate")

CELLS = {
    "api_corpus": {
        "brief_rule": (
            "You may use the CPIC data API (api.cpicpgx.org), PharmGKB or PharmVar in addition "
            "to the supplied corpus. Record which source and version you used for each gene."
        ),
    },
    "corpus_only": {
        "brief_rule": (
            "Do not use any network source. Work only from the supplied corpus and your own "
            "reasoning. If a gene's guideline text is insufficient, say so and abstain rather "
            "than fetching it. Deny any network request this session offers."
        ),
    },
}

BRIEF = """# Blinded pharmacogenomics study

You are running an independent pharmacogenomics study. Work end to end and keep every step
reproducible.

## Question

Given a patient genotype description and a named drug, what is the patient's diplotype, CPIC
metaboliser (or function) phenotype, and CPIC Level A drug recommendation? How accurately can
this be done, and where do errors occur, especially in cases where standard dosing can cause
severe or life-threatening harm?

## Inputs

- `cases_blinded.json`: 110 cases, each with `id`, `gene`, `drug`, `genotype`. No answers are
  included, and the identifiers are opaque hashes that carry no information.
- `cpic_guideline_corpus.json`: CPIC guideline text for the 21 genes (snapshot of 2026-05-09).

## Rules

1. {rule}
2. Do not look for or use any published benchmark, paper or repository that reports answers to
   these cases. If you come across one, stop and say so.
3. Build whatever you think is best: lookup tables, code, or reasoning. Record which method you
   used.
4. Abstain when the input cannot be mapped with confidence, and explain why. An abstention is
   better than a guess.

## Required outputs

1. `claude_science_answers.csv` with exactly these columns, one row per case id:
   `case_id,diplotype,phenotype,recommendation,abstained`
   `recommendation` should name the drug and the action, e.g. "codeine: avoid; use a
   non-CYP2D6 opioid". `abstained` is true or false. Use the case ids exactly as given.
2. Any diplotype-to-phenotype and phenotype-to-recommendation tables you built, as CSV.
3. A short methods note: approach, sources and versions, which cases you consider lethal-class
   and why, and your own estimate of where you are most likely wrong.
4. The code and notebook that produced 1 and 2.
"""


def _builder():
    spec = spec_from_file_location("build_pkg_inner", HERE / "PYTHON" / "00-build_package.py")
    mod = module_from_spec(spec)
    sys.modules["build_pkg_inner"] = mod
    spec.loader.exec_module(mod)
    return mod


def strip_corpus(corpus: dict) -> dict:
    """Remove the per-gene metadata that reveals the corpus was assembled around
    this case set. Returns a copy; the input is not mutated."""
    out = copy.deepcopy(corpus)
    for gene in out.get("genes", {}).values():
        for k in CASE_SET_METADATA:
            gene.pop(k, None)
    return out


def build(cases: list[dict], corpus: dict, root: Path, reps: int = 3) -> list[tuple[Path, dict]]:
    """Write every cell x replicate folder. Returns [(folder, id_map), ...]."""
    b = _builder()
    root = Path(root)
    blind = b.blind_cases(cases, hash_ids=True)
    id_map = b.id_map(cases)
    clean = strip_corpus(corpus)

    keys = root / "_KEYS"
    keys.mkdir(parents=True, exist_ok=True)
    (keys / "id_map.json").write_text(json.dumps(id_map, indent=1))

    made = []
    for cell, cfg in CELLS.items():
        for rep in range(1, reps + 1):
            d = root / cell / f"rep{rep}"
            d.mkdir(parents=True, exist_ok=True)
            (d / "cases_blinded.json").write_text(json.dumps(blind, indent=1))
            (d / "cpic_guideline_corpus.json").write_text(json.dumps(clean, indent=1))
            (d / "BRIEF.md").write_text(BRIEF.format(rule=cfg["brief_rule"]))
            made.append((d, id_map))
    return made


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--out", type=Path, default=HERE / "RUNS")
    args = ap.parse_args()

    specs = HERE.parent / "SPECS"
    cases = json.loads((specs / "test_cases_v3.json").read_text())
    # The source is cpic_rag_corpus_v3.json; cpic_guideline_corpus.json is only
    # the name given to the stripped copy inside each handed folder.
    corpus = json.loads((specs / "cpic_rag_corpus_v3.json").read_text())
    made = build(cases, corpus, args.out, reps=args.reps)
    for d, _ in made:
        print(d)
    print(f"\n{len(made)} run folders, {len(cases)} cases each. Key: {args.out / '_KEYS'}")


if __name__ == "__main__":
    main()
