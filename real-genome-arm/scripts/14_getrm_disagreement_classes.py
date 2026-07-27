#!/usr/bin/env python3
"""
Decompose GeT-RM caller disagreements by mechanism (revision item N5b).

WHY THIS EXISTS
08_caller_truth_eval.py reports a single external concordance figure. On the
1000 Genomes phase 3 cohort that figure is 0.7609, and quoting it as "PyPGx call
accuracy" would be wrong in a way that matters: most disagreements are not the
caller getting a haplotype wrong. They are the truth set and the caller speaking
different versions of star-allele nomenclature, or the caller being asked to
resolve variation that the input data physically cannot carry.

Reporting 0.7609 unqualified would understate a deterministic caller and hand a
reviewer an easy objection. Reporting only the clean subset would be the
opposite error: choosing the denominator that flatters. This script reports
both, plus the mechanism split, so the reader can see what the number is made
of.

THE CLASSES, and what each claim rests on
  SLCO1B1, UGT1A1 -- allele-definition / nomenclature version.
      GeT-RM's SLCO1B1 calls use the legacy *1A/*1B haplotype naming; PyPGx uses
      current PharmVar definitions, in which the same haplotypes carry different
      numbers (*1B against *37 is the dominant pair here). UGT1A1 differs in
      which alleles are defined at all: GeT-RM reports *60, PyPGx does not
      define it and reports the linked *80+*28 haplotype.
  CYP2D6 -- structural variation.
      Whole-gene deletions (*5) and duplications (xN) require read depth.
      1000G phase 3 is low-coverage with imputation and carries no per-sample
      depth, so these are unresolvable from this input by construction.
  CYP2B6 -- *6 against *9.
      *6 is *9 plus 785A>G. Where that site is absent from the phase 3 release,
      *9 is the correct call from the data available.

WHAT THIS IS NOT
The classification is by pattern, not by per-sample curation. A classified
disagreement is *attributable to* a mechanism, not proven to be benign. The
unexplained remainder is reported separately and is the honest residual: it is
the part that may be genuine calling error.

USAGE
    python real-genome-arm/scripts/14_getrm_disagreement_classes.py
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARM = HERE.parent
REPO = ARM.parent

TRUTH_DEFAULT = ARM / "getrm" / "getrm_consensus.tsv"
CALLED_DEFAULT = ARM / "work" / "getrm-1000g" / "pypgx_calls.tsv"
OUT_DEFAULT = REPO / "data" / "v3_getrm_disagreement_classes.json"
REPORT_DEFAULT = REPO / "data" / "v3_getrm_disagreement_classes.txt"

# Genes whose truth and caller definitions align and which carry no CNV calls
# here. Fixed in this file rather than chosen after seeing the per-gene results.
SHARED_DEFINITION = ("CYP2C19", "CYP2C9", "CYP3A5", "TPMT")

NOMENCLATURE_GENES = ("SLCO1B1", "UGT1A1")
CNV_RE = re.compile(r"x\d|(?:^|/)\*5(?:$|/)")
CYP2B6_RE = re.compile(r"\*6|\*9")


def _norm():
    spec = spec_from_file_location("n0", HERE / "07_executed_pipeline_n0.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.norm_dip


def load(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}")
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    missing = [c for c in ("sample", "gene", "diplotype") if rows and c not in rows[0]]
    if missing:
        raise SystemExit(f"{path} is missing columns: {missing}")
    return {(r["sample"], r["gene"]): r["diplotype"] for r in rows}


def classify(gene: str, truth: str, called: str) -> str:
    if gene == "CYP2D6" and (CNV_RE.search(truth) or CNV_RE.search(called)):
        return "CYP2D6 structural variant (needs read depth)"
    if gene in NOMENCLATURE_GENES:
        return f"{gene} allele-definition / nomenclature version"
    if gene == "CYP2B6" and CYP2B6_RE.search(truth + called):
        return "CYP2B6 *6 against *9 (785A>G site coverage)"
    return f"{gene} unexplained"


def analyse(truth: dict, called: dict) -> dict:
    norm_dip = _norm()
    keys = sorted(set(truth) & set(called))
    per = defaultdict(lambda: {"n": 0, "correct": 0})
    classes = Counter()
    examples = defaultdict(Counter)

    for k in keys:
        sample, gene = k
        t, c = norm_dip(truth[k], gene), norm_dip(called[k], gene)
        per[gene]["n"] += 1
        if t == c:
            per[gene]["correct"] += 1
            continue
        label = classify(gene, truth[k], called[k])
        classes[label] += 1
        examples[label][(truth[k], called[k])] += 1

    n = sum(v["n"] for v in per.values())
    correct = sum(v["correct"] for v in per.values())
    unexplained = sum(v for k, v in classes.items() if k.endswith("unexplained"))
    sn = sum(per[g]["n"] for g in SHARED_DEFINITION if g in per)
    sc = sum(per[g]["correct"] for g in SHARED_DEFINITION if g in per)

    return {
        "evaluable": n,
        "concordant": correct,
        "concordance": round(correct / n, 4) if n else None,
        "disagreements": n - correct,
        "attributable": (n - correct) - unexplained,
        "unexplained": unexplained,
        "classes": dict(classes.most_common()),
        "examples": {k: [{"truth": t, "called": c, "n": v}
                         for (t, c), v in ex.most_common(3)]
                     for k, ex in examples.items()},
        "shared_definition_genes": list(SHARED_DEFINITION),
        "shared_definition_n": sn,
        "shared_definition_correct": sc,
        "shared_definition_concordance": round(sc / sn, 4) if sn else None,
        "per_gene": {g: {**v, "concordance": round(v["correct"] / v["n"], 4)}
                     for g, v in sorted(per.items())},
    }


def render(r: dict) -> str:
    L = ["GeT-RM CALLER CONCORDANCE, DECOMPOSED", ""]
    L.append(f"  evaluable (sample, gene) pairs   {r['evaluable']}")
    L.append(f"  concordant                       {r['concordant']}")
    L.append(f"  concordance                      {r['concordance']}")
    L.append(f"  disagreements                    {r['disagreements']}")
    L.append(f"    attributable to a mechanism    {r['attributable']}")
    L.append(f"    unexplained                    {r['unexplained']}")
    L.append("")
    L.append("  disagreements by mechanism:")
    for k, v in r["classes"].items():
        L.append(f"    {v:4d}  {k}")
    L.append("")
    L.append(f"  genes with shared definitions and no CNV "
             f"({', '.join(r['shared_definition_genes'])}):")
    L.append(f"    {r['shared_definition_correct']}/{r['shared_definition_n']} "
             f"= {r['shared_definition_concordance']}")
    L.append("")
    L.append("  per gene:")
    for g, v in r["per_gene"].items():
        L.append(f"    {g:10s} {v['correct']:4d}/{v['n']:<4d} {v['concordance']}")
    L.append("")
    L.append("  Input is 1000 Genomes phase 3: low-coverage, imputed, no per-sample")
    L.append("  read depth. This evaluates the caller on genotype data of that kind,")
    L.append("  the same mode used for the Peru, Uganda and IBS arms. It does not")
    L.append("  evaluate it on high-depth sequencing, and CYP2D6 structural variation")
    L.append("  cannot be resolved from this input at all.")
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--truth", type=Path, default=TRUTH_DEFAULT)
    ap.add_argument("--called", type=Path, default=CALLED_DEFAULT)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--report", type=Path, default=REPORT_DEFAULT)
    a = ap.parse_args(argv)

    result = analyse(load(a.truth), load(a.called))
    if not result["evaluable"]:
        sys.stderr.write("no overlapping (sample, gene) pairs; refusing to report\n")
        return 1
    a.out.write_text(json.dumps(result, indent=2))
    text = render(result)
    a.report.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
