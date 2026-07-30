#!/usr/bin/env python3
"""
Evaluate the input-normalisation experiment (reviewer point R1.1).

Reports coverage, abstention and accuracy among emitted calls, against two
references kept strictly apart:

  vs_caller   PyPGx, the deterministic caller the paper compares against
  vs_getrm    GeT-RM external consensus

They are never pooled. The paper already reports that the two agree on only
0.761 of these pairs, so a blended accuracy would belong to neither.

Responses that exhausted the output budget without reaching the DIPLOTYPE line
are counted as `truncated_output`, reported separately, and excluded from the
scored denominator. They are a harness limit, not a model declining to answer,
and the first run of this experiment would have published 100 per cent o3
abstention if they had been absorbed into the abstention rate.

USAGE
    python code/75-evaluate-input-normalisation.py \
        data/v3_input_normalisation_main.json data/v3_input_normalisation_o3.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
REPORT_OUT = BASE / "data" / "v3_input_normalisation.txt"
JSON_OUT = BASE / "data" / "v3_input_normalisation_eval.json"


def _norm_dip(d: str) -> tuple:
    return tuple(sorted(p.strip() for p in d.split("/")))


def _arm(rows: list[dict], ref_key: str) -> dict:
    rows = [r for r in rows if r.get(ref_key)]
    n = len(rows)
    emitted = [r for r in rows if r.get("call")]
    correct = [r for r in emitted
               if _norm_dip(r["call"]) == _norm_dip(r[ref_key])]
    return {
        "n": n,
        "emitted": len(emitted),
        "coverage": (len(emitted) / n) if n else 0.0,
        "abstention": (1 - len(emitted) / n) if n else 0.0,
        "correct": len(correct),
        "accuracy_among_emitted": (len(correct) / len(emitted)) if emitted else None,
    }


def _cell(rows: list[dict]) -> dict:
    return {"vs_caller": _arm(rows, "reference"), "vs_getrm": _arm(rows, "getrm")}


def evaluate(rows: list[dict]) -> dict:
    """Overall and per model / rendering / gene.

    Rows whose response was truncated by the output budget are excluded from
    the scored set and counted separately, so they neither inflate abstention
    nor shrink the denominator silently.
    """
    truncated = [r for r in rows if r.get("status") == "truncated_output"]
    scored = [r for r in rows if r.get("status") != "truncated_output"]

    def group(key):
        out = {}
        for r in scored:
            out.setdefault(r[key], []).append(r)
        return {k: _cell(v) for k, v in sorted(out.items())}

    return {
        "overall": _cell(scored),
        "by_model": group("model"),
        "by_form": group("form"),
        "by_gene": group("gene"),
        "truncated_output": len(truncated),
        "errors": sum(1 for r in rows if r.get("error")),
        "n_rows": len(rows),
        "n_scored": len(scored),
    }


def _fmt(c: dict) -> str:
    a = c["accuracy_among_emitted"]
    return (f"n={c['n']:5d}  cov={c['coverage']:.3f}  abst={c['abstention']:.3f}  "
            f"acc={'  n/a' if a is None else f'{a:.3f}'}  ({c['correct']}/{c['emitted']})")


def report(ev: dict) -> str:
    L = ["INPUT NORMALISATION (R1.1): can the model do the job the paper assigns it?", ""]
    L.append(f"rows {ev['n_rows']}   scored {ev['n_scored']}   "
             f"truncated_output {ev['truncated_output']}   errors {ev['errors']}")
    L.append("")
    for ref, label in (("vs_caller", "against PyPGx (deterministic caller)"),
                       ("vs_getrm", "against GeT-RM (external consensus)")):
        L.append(f"  {label}")
        L.append(f"    OVERALL              {_fmt(ev['overall'][ref])}")
        for section, title in (("by_model", "model"), ("by_form", "rendering"),
                               ("by_gene", "gene")):
            L.append(f"    by {title}:")
            for k, v in ev[section].items():
                L.append(f"      {k:20s} {_fmt(v[ref])}")
        L.append("")
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputs", nargs="+", type=Path)
    args = ap.parse_args(argv)

    rows = []
    for p in args.inputs:
        if not p.exists():
            print(f"missing {p}", file=sys.stderr)
            return 2
        rows += json.loads(p.read_text())["rows"]

    ev = evaluate(rows)
    text = report(ev)
    print(text)
    REPORT_OUT.write_text(text + "\n")
    JSON_OUT.write_text(json.dumps(ev, indent=2))
    print(f"\nwrote {REPORT_OUT}\nwrote {JSON_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
