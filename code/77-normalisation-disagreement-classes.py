#!/usr/bin/env python3
"""
Decompose the model's normalisation disagreements by mechanism (R1.1).

WHY THIS EXISTS
`14_getrm_disagreement_classes.py` decomposes the deterministic caller's
disagreements with GeT-RM rather than quoting its raw 0.761 unqualified, on the
grounds that most of them are nomenclature-version differences or variation the
input physically cannot carry. That was the right standard to apply to PyPGx.
Applying a weaker standard to the model would be choosing the denominator that
flatters our own argument, which is the error C6 in CORRECTIONS.md ran the other
way and which this project has now made four times.

So the model's wrong answers get the same decomposition, using the same classes
and the same code path where possible. If the model's errors are largely legacy
nomenclature, that is a materially different finding from the model calling the
wrong haplotype, and the paper must say which.

WHAT THIS IS NOT
The classification is by pattern, not per-sample curation. A classified
disagreement is attributable to a mechanism, not proven benign. The unexplained
remainder is the honest residual and is reported separately.

Note that coverage is unaffected by any of this. A model that abstains on 76 per
cent of inputs has not been rescued by a nomenclature argument about the
remaining 24.

USAGE
    python code/77-normalisation-disagreement-classes.py \
        data/v3_input_normalisation_main.json data/v3_input_normalisation_o3.json
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "data" / "v3_normalisation_disagreement_classes.json"
TXT = BASE / "data" / "v3_normalisation_disagreement_classes.txt"

# Legacy SLCO1B1 haplotype naming: GeT-RM and older references use *1A/*1B where
# current PharmVar definitions renumber the same haplotypes (*1B against *37 is
# the dominant pair). UGT1A1 differs in which alleles are defined at all.
LEGACY_SLCO1B1 = re.compile(r"\*1[AB]\b")
STRUCTURAL = re.compile(r"\*5\b|x\d")


def classify(gene: str, call: str, reference: str) -> str:
    """Attribute one disagreement to a mechanism, or leave it unexplained."""
    if gene == "SLCO1B1" and (LEGACY_SLCO1B1.search(call) or LEGACY_SLCO1B1.search(reference)):
        return "nomenclature_version"
    if gene == "UGT1A1":
        return "nomenclature_version"
    if gene == "CYP2D6" and (STRUCTURAL.search(call) or STRUCTURAL.search(reference)):
        return "structural_variation"
    if gene == "CYP2B6" and {"*6", "*9"} & (set(call.split("/")) | set(reference.split("/"))):
        return "missing_site_6_vs_9"
    return "unexplained"


def _norm(d: str) -> tuple:
    return tuple(sorted(p.strip() for p in d.split("/")))


def decompose(rows: list[dict], ref_key: str = "reference") -> dict:
    emitted = [r for r in rows if r.get("status") == "call" and r.get(ref_key)]
    wrong = [r for r in emitted if _norm(r["call"]) != _norm(r[ref_key])]
    classes = Counter(classify(r["gene"], r["call"], r[ref_key]) for r in wrong)
    attributable = sum(v for k, v in classes.items() if k != "unexplained")
    return {
        "emitted": len(emitted),
        "correct": len(emitted) - len(wrong),
        "wrong": len(wrong),
        "classes": dict(classes),
        "attributable": attributable,
        "unexplained": classes.get("unexplained", 0),
        "accuracy_raw": (len(emitted) - len(wrong)) / len(emitted) if emitted else None,
        "accuracy_if_nomenclature_forgiven": (
            (len(emitted) - classes.get("unexplained", 0)
             - classes.get("structural_variation", 0)
             - classes.get("missing_site_6_vs_9", 0)) / len(emitted)
            if emitted else None),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputs", nargs="+", type=Path)
    args = ap.parse_args(argv)

    rows = []
    for p in args.inputs:
        rows += json.loads(p.read_text())["rows"]

    out = {"vs_caller": decompose(rows, "reference"),
           "vs_getrm": decompose(rows, "getrm")}

    L = ["NORMALISATION DISAGREEMENTS, DECOMPOSED BY MECHANISM", ""]
    L.append("The same decomposition applied to the deterministic caller in")
    L.append("14_getrm_disagreement_classes.py, applied to the model.")
    L.append("")
    for k, d in out.items():
        L.append(f"  {k}")
        L.append(f"    emitted {d['emitted']}  correct {d['correct']}  wrong {d['wrong']}")
        for cls, n in sorted(d["classes"].items(), key=lambda x: -x[1]):
            L.append(f"      {cls:24s} {n:5d}")
        L.append(f"    attributable to a mechanism : {d['attributable']}")
        L.append(f"    unexplained (honest residual): {d['unexplained']}")
        L.append(f"    accuracy, raw                : {d['accuracy_raw']:.3f}")
        L.append(f"    accuracy, nomenclature forgiven: "
                 f"{d['accuracy_if_nomenclature_forgiven']:.3f}")
        L.append("")
    L.append("Coverage is untouched by this decomposition. It concerns only the")
    L.append("answers the model chose to emit.")
    text = "\n".join(L)
    print(text)
    TXT.write_text(text + "\n")
    OUT.write_text(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
