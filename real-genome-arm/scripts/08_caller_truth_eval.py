#!/usr/bin/env python3
"""
Validate the CALLING step against external truth (revision item N5).

N0 measured the mapping step: PyPGx called the diplotype and also supplied the
phenotype it was scored against, so input-call error was unobservable. This
script closes that loop. Given a truth table of consensus diplotypes (GeT-RM)
and a table of called diplotypes, it reports per-gene call concordance and,
crucially, the quantity the paper's error-localisation claim rests on: how much
end-to-end error is attributable to the input call rather than to the mapping.

It is deliberately source-agnostic. GeT-RM is the truth set the plan commits
to, but any (sample, gene, diplotype) table works, so the harness does not rot
if the truth source changes.

USAGE
    python real-genome-arm/scripts/08_caller_truth_eval.py \
        --truth real-genome-arm/getrm/getrm_consensus.tsv \
        --called path/to/pypgx_calls.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
REQUIRED = ("sample", "gene", "diplotype")


def _load_n0():
    spec = spec_from_file_location("n0_pipeline", HERE / "07_executed_pipeline_n0.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_calls(path: Path) -> list[dict]:
    """Read a (sample, gene, diplotype) TSV, refusing anything malformed.

    An empty or wrong-shaped truth set must stop the run. Silently evaluating
    against nothing would report a vacuous 100% concordance, which is the exact
    failure mode this whole revision is guarding against.
    """
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit(f"{path} contains no rows; refusing to report a vacuous result")
    missing = [c for c in REQUIRED if c not in rows[0]]
    if missing:
        raise SystemExit(f"{path} is missing required columns: {missing}")
    return rows


def evaluate(truth: list[dict], called: list[dict]) -> dict:
    n0 = _load_n0()

    def key(r):
        return (r["sample"].strip(), r["gene"].strip())

    truth_map = {key(r): r["diplotype"].strip() for r in truth}
    called_map = {key(r): r["diplotype"].strip() for r in called}
    shared = sorted(set(truth_map) & set(called_map))

    per_gene: dict[str, dict] = defaultdict(lambda: {"n": 0, "correct": 0, "errors": []})
    for k in shared:
        sample, gene = k
        t = n0.norm_dip(truth_map[k], gene)
        c = n0.norm_dip(called_map[k], gene)
        g = per_gene[gene]
        g["n"] += 1
        if t == c:
            g["correct"] += 1
        else:
            g["errors"].append({"sample": sample, "truth": truth_map[k],
                                "called": called_map[k]})

    n = sum(g["n"] for g in per_gene.values())
    correct = sum(g["correct"] for g in per_gene.values())
    in_vocab = sum(1 for k in shared
                   if (k[1], n0.norm_dip(truth_map[k], k[1])) in n0.SKILL_MAP)
    return {
        "truth_rows": len(truth_map),
        "called_rows": len(called_map),
        "evaluable": n,
        "truth_without_call": len(set(truth_map) - set(called_map)),
        "call_concordance": round(correct / n, 4) if n else None,
        "call_errors": n - correct,
        "truth_states_in_skill_vocabulary": in_vocab,
        "per_gene": {g: {"n": v["n"], "correct": v["correct"],
                         "concordance": round(v["correct"] / v["n"], 4) if v["n"] else None,
                         "errors": v["errors"][:20]}
                     for g, v in sorted(per_gene.items())},
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", type=Path, required=True)
    ap.add_argument("--called", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=BASE / "data" / "v3_caller_truth_eval.json")
    args = ap.parse_args(argv)

    for p in (args.truth, args.called):
        if not p.exists():
            sys.stderr.write(
                f"missing {p}\nSee real-genome-arm/getrm/README.md for how to obtain "
                "the GeT-RM consensus table.\n")
            return 1

    result = evaluate(load_calls(args.truth), load_calls(args.called))
    args.out.write_text(json.dumps(result, indent=2))
    print(f"evaluable (sample, gene) pairs   {result['evaluable']}")
    print(f"call concordance                 {result['call_concordance']}")
    print(f"call errors                      {result['call_errors']}")
    print(f"truth states in skill vocabulary {result['truth_states_in_skill_vocabulary']}")
    for gene, g in result["per_gene"].items():
        print(f"   {gene:<12} {g['correct']}/{g['n']}  {g['concordance']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
