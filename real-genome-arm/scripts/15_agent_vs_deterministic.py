#!/usr/bin/env python3
"""
Agent against a deterministic pipeline on real cohorts (reviewer point 8).

THE OBJECTION
The referee asked twice for a direct comparison against a conventional
pipeline, to establish what the agent contributes beyond input normalisation
and orchestration. The revision compared against PharmCAT only on the 110
curated synthetic cases, using a vocabulary we defined ourselves. On the real
cohorts, where the question is live, there was no non-LLM baseline at all.

WHAT THIS COMPARES
Both arms are scored against the same external reference: the CPIC phenotype
returned by the deterministic caller (PyPGx) for each observed diplotype.

  deterministic arm   the validated skill's mapping applied in code to the
                      called diplotype; abstains when the diplotype is outside
                      the skill's vocabulary
  agent arm           each panel model asked to interpret the same diplotype

The reference is independent of the skill's authored table: PyPGx derives its
phenotype from its own allele-function definitions, so agreement between the
executed skill and the reference is a measurement rather than a tautology. This
is the distinction the referee drew for the 110-case comparison, and it is
respected here.

WHAT THE COMPARISON CANNOT SETTLE
The deterministic arm is handed a diplotype that a deterministic caller already
produced. It is not doing input interpretation, which is exactly the step the
agent exists to perform on heterogeneous input. So this measures the
interpretation step alone, and a deterministic pipeline winning it does not show
the agent is useless; it shows the agent should not be the one doing the
mapping. That is the paper's claim, and this is the test of it on real data.

USAGE
    python real-genome-arm/scripts/15_agent_vs_deterministic.py \
        --preds data/v3_realgenome_preds_4cohorts.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARM = HERE.parent
REPO = ARM.parent

STATES = REPO / "data" / "n5_four_cohorts.tsv"
PREDS = REPO / "data" / "v3_realgenome_preds_4cohorts.tsv"
OUT = REPO / "data" / "v3_agent_vs_deterministic.json"
REPORT = REPO / "data" / "v3_agent_vs_deterministic.txt"


def _mod(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    m = module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


rules = _mod(REPO / "code" / "_pgx_rules.py", "_pgx_rules")


def phen_key(text: str) -> str:
    """Compare phenotype tiers, not wording, spelling or capitalisation.

    Hyphens and spaces are stripped before matching. Without that,
    "Ultra-rapid Metaboliser" misses the "ultrarapid" needle, falls through to
    "rapid", and is scored as disagreeing with "Ultrarapid Metabolizer" -- a
    spelling difference reported as a clinical one.
    """
    t = (text or "").strip().lower().replace("-", "").replace(" ", "")
    if not t:
        return ""
    likely = "likely" in t
    for needle, tier in (("ultrarapid", "ultrarapid"), ("rapid", "rapid"),
                         ("normal", "normal"), ("extensive", "normal"),
                         ("intermediate", "intermediate"), ("poor", "poor"),
                         ("indeterminate", "indeterminate"),
                         ("uncertain", "indeterminate"),
                         ("increased", "increased"), ("decreased", "decreased"),
                         ("nofunction", "poor"), ("positive", "positive"),
                         ("negative", "negative")):
        if needle in t:
            return ("likely " if likely else "") + tier
    return t


def is_reference_abstention(phenotype: str) -> bool:
    """True when the CALLER declined to assign a phenotype.

    PyPGx returns Indeterminate for genes whose alleles it does not assign a
    function to (CYP4F2 here). Those states carry no reference answer, so
    scoring either arm against them measures the reference's silence. N0 treated
    them as abstention targets for exactly this reason; doing otherwise counted
    eight reference abstentions as skill errors and dropped the deterministic
    arm from ~1.00 to 0.82.
    """
    return phen_key(phenotype) in ("", "indeterminate")


def load_states(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def deterministic_arm(states: list[dict], dip2phen: dict) -> dict:
    per = defaultdict(lambda: {"n": 0, "emitted": 0, "correct": 0})
    for s in states:
        gene, cohort = s["gene"], s["cohort"]
        key = (gene, rules.norm_dip(s["diplotype"], gene))
        if is_reference_abstention(s["phenotype"]):
            continue
        got = dip2phen.get(key)
        c = per[cohort]
        c["n"] += 1
        if got is None:
            continue
        c["emitted"] += 1
        if phen_key(got) == phen_key(s["phenotype"]):
            c["correct"] += 1
    return {k: summarise(v) for k, v in per.items()}


def agent_arm(preds: list[dict]) -> dict:
    per = defaultdict(lambda: {"n": 0, "emitted": 0, "correct": 0})
    for r in preds:
        cohort = r.get("cohort", "")
        pred = r.get("pred", "")
        ref = r.get("cohort_phenotype", "")
        if is_reference_abstention(ref):
            continue
        c = per[cohort]
        c["n"] += 1
        pk = phen_key(pred)
        if not pk or pk == "indeterminate":
            continue          # abstention, counted but not scored as an error
        c["emitted"] += 1
        if pk == phen_key(ref):
            c["correct"] += 1
    return {k: summarise(v) for k, v in per.items()}


def summarise(v: dict) -> dict:
    return {
        "n": v["n"], "emitted": v["emitted"], "correct": v["correct"],
        "abstained": v["n"] - v["emitted"],
        "coverage": round(v["emitted"] / v["n"], 4) if v["n"] else None,
        "accuracy_among_emitted": (round(v["correct"] / v["emitted"], 4)
                                   if v["emitted"] else None),
        "accuracy_overall": round(v["correct"] / v["n"], 4) if v["n"] else None,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--states", type=Path, default=STATES)
    ap.add_argument("--preds", type=Path, default=PREDS)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--report", type=Path, default=REPORT)
    a = ap.parse_args(argv)

    if not a.preds.exists():
        sys.stderr.write(f"missing {a.preds}; run 04_run_agent_realgenome.py first\n")
        return 1
    states = load_states(a.states)
    with a.preds.open() as fh:
        preds = list(csv.DictReader(fh, delimiter="\t"))

    dip2phen, _, _, _ = rules.build_rules(rules.load_cases())
    det = deterministic_arm(states, dip2phen)

    by_model = defaultdict(list)
    for r in preds:
        by_model[r.get("model", "?")].append(r)
    agent_all = agent_arm(preds)
    agent_by_model = {m: agent_arm(rows) for m, rows in sorted(by_model.items())}

    result = {"deterministic": det, "agent_pooled": agent_all,
              "agent_by_model": agent_by_model}
    a.out.write_text(json.dumps(result, indent=2))

    cohorts = sorted(set(det) | set(agent_all))
    L = ["AGENT AGAINST A DETERMINISTIC PIPELINE, ON REAL COHORTS", ""]
    L.append("  Both arms scored against the same external reference: the CPIC phenotype")
    L.append("  returned by the deterministic caller for each observed diplotype.")
    L.append("")
    L.append(f"  {'cohort':14s} {'arm':14s} {'n':>5} {'cov':>7} {'acc|emit':>9} {'acc':>7}")
    for c in cohorts:
        for name, src in (("deterministic", det), ("agent (pooled)", agent_all)):
            v = src.get(c)
            if not v:
                continue
            L.append(f"  {c:14s} {name:14s} {v['n']:5d} {str(v['coverage']):>7} "
                     f"{str(v['accuracy_among_emitted']):>9} {str(v['accuracy_overall']):>7}")
        L.append("")
    L.append("  The deterministic arm receives a diplotype the caller already produced, so")
    L.append("  it performs no input interpretation. This compares the mapping step only.")
    L.append("  A deterministic pipeline winning it does not show the agent is redundant;")
    L.append("  it shows the agent should not be the component doing the mapping, which is")
    L.append("  the architecture this paper argues for.")
    text = "\n".join(L)
    a.report.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
