#!/usr/bin/env python3
"""Per-gene input-normalisation table for the supplementary (Table S7).

WHY THIS EXISTS
The supplementary caption for this table asserted that its values were
"computed from the deposited raw rows by build-supp-v28.py and are not
transcribed". No such script existed at any tag, and the three surviving
copies of the table disagreed with each other: the caption quoted CYP2D6
accuracy as 0.000 to 0.825, the table body printed 0.030 and 0.789, and
data/v3_normalisation_by_gene.txt printed a third set on pre-C15 coverage.
A caption that claims a generator is the strongest possible statement that
the numbers were not typed by hand, so it has to be true. This is that
generator.

WHAT IT COMPUTES
One row per gene for Claude Opus 4.5 on the variant-call rendering, over the
527 (sample, gene) pairs, in both arms:

  coverage  call / (call + abstention), the behavioural denominator. Provider
            errors and output-budget truncations are excluded from both
            denominators, because neither is the model declining to answer.
  accuracy  among emitted calls only, against the deterministic caller
            (PyPGx). Diplotypes are unordered, so *1/*4 and *4/*1 match.

Source files and the arm-to-file mapping are taken from the eight-model
freeze (85-freeze-eight-model-normalisation.py) so that this table and
Figure 9 cannot drift apart. The script asserts its own column totals against
the freeze's Claude Opus 4.5 operational counts and fails loudly if they
disagree, rather than emitting a plausible table nobody checks.

USAGE
    python code/88-supp-normalisation-by-gene.py
    python code/88-supp-normalisation-by-gene.py --out data/v3_normalisation_by_gene.txt
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
FREEZE = DATA / "v3_input_normalisation_eight_model_freeze.json"
OUT_DEFAULT = DATA / "v3_normalisation_by_gene.txt"
JSON_DEFAULT = DATA / "v3_normalisation_by_gene.json"

MODEL = "Claude Opus 4.5"
FORM = "vcf"

# Identical to the mapping in 85-freeze-eight-model-normalisation.py. Kept
# explicit rather than imported, because the freeze module has no importable
# name (it starts with a digit) and a silent divergence would be invisible;
# the totals assertion at the end is what actually binds the two together.
NO_DEFINITION_FILES = ["v3_input_normalisation_main.json"]
DEFINITION_FILES = ["v3_input_normalisation_defs.json",
                    "v3_input_normalisation_defs_tail.json"]


def unordered(diplotype: str) -> tuple[str, ...]:
    return tuple(sorted(part.strip() for part in diplotype.split("/")))


def load_rows(names: list[str]) -> list[dict]:
    rows: list[dict] = []
    for name in names:
        blob = json.loads((DATA / name).read_text())
        rows.extend(r for r in blob["rows"]
                    if r["model"] == MODEL and r.get("form", FORM) == FORM)
    return rows


def summarise(rows: list[dict]) -> dict[str, dict]:
    """Per-gene call/abstain/error/truncated counts and correct-call counts."""
    per: dict[str, dict] = {}
    for row in rows:
        gene = row["gene"]
        bucket = per.setdefault(gene, {"call": 0, "abstain": 0, "error": 0,
                                       "truncated_output": 0, "correct": 0,
                                       "scorable": 0})
        status = row["status"]
        if status not in bucket:
            raise ValueError(f"unexpected status {status!r} for {gene}")
        bucket[status] += 1
        if status != "call":
            continue
        reference = row.get("reference")
        if not reference:
            # No caller answer means the call cannot be scored either way; it
            # still counts as behaviour, so it stays in the coverage numerator.
            continue
        bucket["scorable"] += 1
        if unordered(row["call"]) == unordered(reference):
            bucket["correct"] += 1
    return per


def rate(numerator: int, denominator: int) -> float | None:
    return round(numerator / denominator, 3) if denominator else None


def fmt(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.3f}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--json-out", type=Path, default=JSON_DEFAULT)
    args = ap.parse_args(argv)

    nodefs = summarise(load_rows(NO_DEFINITION_FILES))
    withdefs = summarise(load_rows(DEFINITION_FILES))

    genes = sorted(set(nodefs) | set(withdefs))
    empty = {"call": 0, "abstain": 0, "error": 0, "truncated_output": 0,
             "correct": 0, "scorable": 0}

    table = []
    for gene in genes:
        a = nodefs.get(gene, empty)
        b = withdefs.get(gene, empty)
        pairs = a["call"] + a["abstain"] + a["error"] + a["truncated_output"]
        table.append({
            "gene": gene,
            "pairs": pairs,
            "coverage_no_definitions": rate(a["call"], a["call"] + a["abstain"]),
            "accuracy_no_definitions": rate(a["correct"], a["scorable"]),
            "coverage_with_definitions": rate(b["call"], b["call"] + b["abstain"]),
            "accuracy_with_definitions": rate(b["correct"], b["scorable"]),
        })

    # Bind this table to the frozen analysis. If the arm-to-file mapping above
    # ever drifts from the freeze, these totals stop matching and the script
    # refuses to write a table that would quietly contradict Figure 9.
    freeze = json.loads(FREEZE.read_text())["models_analysis"][MODEL]
    checks = [
        ("no definitions", nodefs, freeze["without_definitions"]["operational"]),
        ("with definitions", withdefs, freeze["with_definitions"]["operational"]),
    ]
    for label, per, expected in checks:
        for key in ("call", "abstain", "error", "truncated_output"):
            got = sum(v[key] for v in per.values())
            if got != expected[key]:
                sys.stderr.write(
                    f"{label}: {key} total {got} does not match the eight-model "
                    f"freeze ({expected[key]}); refusing to write\n")
                return 1

    lines = [
        "Table S7. Input normalisation by gene, Claude Opus 4.5, variant-call rendering,",
        "527 (sample, gene) pairs. Coverage is call/(call + abstention); accuracy is among",
        "emitted calls, against the deterministic caller, with diplotypes unordered.",
        "Provider failures and output-budget truncations are excluded from both denominators.",
        f"Generated by code/{Path(__file__).name} from the deposited raw rows.",
        "",
        f"{'gene':<12}{'pairs':>7}{'cov -defs':>11}{'acc -defs':>11}"
        f"{'cov +defs':>11}{'acc +defs':>11}",
    ]
    for row in table:
        lines.append(
            f"{row['gene']:<12}{row['pairs']:>7}"
            f"{fmt(row['coverage_no_definitions']):>11}"
            f"{fmt(row['accuracy_no_definitions']):>11}"
            f"{fmt(row['coverage_with_definitions']):>11}"
            f"{fmt(row['accuracy_with_definitions']):>11}")

    text = "\n".join(lines) + "\n"
    args.out.write_text(text)
    args.json_out.write_text(json.dumps(
        {"model": MODEL, "rendering": FORM, "rows": table}, indent=2) + "\n")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
