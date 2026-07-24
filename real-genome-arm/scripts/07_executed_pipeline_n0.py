#!/usr/bin/env python3
"""
N0: end-to-end executed pipeline on real genomes (deterministic caller -> executed
skill -> abstention), reported as accuracy among emitted answers alongside coverage.

This is the experiment that turns the manuscript's central prediction into a result.
The deployment architecture the paper proposes is: replace the model's diplotype call
with a validated deterministic caller (PyPGx 0.26.0, GRCh37), execute the validated
skill's CPIC mapping in code, and abstain on no-calls, ambiguous (uncertain-function)
and out-of-scope diplotypes. The claim (Cell Genomics reviewer 2, point 2) is that this
yields 100% correctness among emitted in-scope answers at a measured coverage; here we
measure it rather than predict it.

Ground truth for each emitted answer is the deterministic caller's own CPIC phenotype
(PyPGx), so accuracy among emitted is the agreement of our executed skill with an
INDEPENDENT implementation of CPIC. Any disagreement is recorded, not hidden: it is the
falsifier for the claim.

No new model calls: PyPGx and the executed skill are deterministic.

Input: the per-cohort aggregated diplotype table from 03_aggregate_diplotypes.py,
       columns: cohort, gene, diplotype, phenotype (PyPGx CPIC phenotype), n_carriers.

Usage:
  07_executed_pipeline_n0.py <diplotypes.tsv> <out.json>   # one or many cohorts in the tsv
  07_executed_pipeline_n0.py --demo                         # runs on a built-in fixture, no data
"""
from __future__ import annotations
import csv
import json
import re
import sys
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
CASES = BASE / "specs" / "test_cases_v3.json"

# Reuse the manuscript's own phenotype-tier scorer (pure: json/re/pathlib, no API).
_spec = spec_from_file_location("rescore", str(BASE / "code" / "10-rescore-v3.py"))
_rs = module_from_spec(_spec)
_spec.loader.exec_module(_rs)
score_a1 = _rs.score_a1  # score_a1(parsed_phen, gt_phen) -> 1.0 / 0.5 / 0.0

NO_CALL_MARKERS = {"", "indeterminate", "no call", "no result", "n/a", "na",
                   "not available", "no data", "none", "unknown"}
AMBIGUOUS_MARKERS = ["uncertain", "unknown function", "possible ", "ambiguous",
                     "indeterminate function"]


def norm_dip(s: str, gene: str) -> str:
    """Normalise a diplotype string the same way the skill runners do."""
    s = (s or "").lower().replace(gene.lower(), "")
    s = re.sub(r"\(.*?\)", "", s).replace("hla-", "")
    s = s.replace("positive", "pos").replace("negative", "neg").replace("carrier", "").replace("present", "pos").replace("absent", "neg")
    s = re.sub(r"[^\w*/:.\- ]", "", s)
    parts = [re.sub(r"\s+", "", p) for p in re.split(r"\s*/\s*", s) if p.strip()]
    if len(parts) == 2:
        parts = sorted(parts)
    return "/".join(parts) if parts else re.sub(r"\s+", "", s)


def build_skill_map(cases) -> dict:
    """The skill's validated (gene, normalised diplotype) -> canonical phenotype map."""
    m = {}
    for c in cases:
        m[(c["gene"], norm_dip(c["gt_diplotype"], c["gene"]))] = c["gt_phenotype"]
    return m


def _load_cases():
    return json.loads(CASES.read_text())


# Built at import so the skill vocabulary is available to callers and tests.
SKILL_MAP = build_skill_map(_load_cases()) if CASES.exists() else {}


def classify(gene: str, diplotype: str, pypgx_phenotype: str, skill_map: dict):
    """Route a real diplotype to 'emitted' or 'abstain' with a reason.

    Abstention reasons: no_call (caller returned nothing usable), ambiguous
    (uncertain-function allele), out_of_scope (diplotype not in the validated skill).
    """
    dl = (diplotype or "").strip()
    pl = (pypgx_phenotype or "").strip().lower()
    if not dl or pl in NO_CALL_MARKERS or "no call" in pl or "indeterminate" in pl:
        return ("abstain", "no_call")
    if any(m in pl for m in AMBIGUOUS_MARKERS):
        return ("abstain", "ambiguous")
    if (gene, norm_dip(dl, gene)) not in skill_map:
        return ("abstain", "out_of_scope")
    return ("emitted", "")


def execute(gene: str, diplotype: str, skill_map: dict) -> str:
    """The executed skill: deterministic (gene, diplotype) -> phenotype in code."""
    return skill_map[(gene, norm_dip(diplotype, gene))]


def run(rows, skill_map: dict, score_fn) -> dict:
    """Aggregate per cohort. rows: dicts with cohort/gene/diplotype/phenotype/n_carriers."""
    agg = defaultdict(lambda: {
        "n_states": 0, "emitted": 0, "correct_emitted": 0,
        "abstain": {"no_call": 0, "ambiguous": 0, "out_of_scope": 0},
        "carriers_total": 0, "carriers_emitted": 0, "disagreements": [],
    })
    for r in rows:
        coh = r["cohort"]; gene = r["gene"]; dip = r["diplotype"]
        ref_phen = r.get("phenotype", "")
        try:
            nc = int(r.get("n_carriers", 0) or 0)
        except (TypeError, ValueError):
            nc = 0
        a = agg[coh]
        a["n_states"] += 1
        a["carriers_total"] += nc
        status, reason = classify(gene, dip, ref_phen, skill_map)
        if status == "abstain":
            a["abstain"][reason] += 1
            continue
        a["emitted"] += 1
        a["carriers_emitted"] += nc
        executed_phen = execute(gene, dip, skill_map)
        if score_fn(executed_phen, ref_phen) == 1.0:
            a["correct_emitted"] += 1
        else:
            a["disagreements"].append({
                "gene": gene, "diplotype": dip,
                "executed_phenotype": executed_phen, "pypgx_phenotype": ref_phen,
            })
    out = {}
    for coh, a in agg.items():
        emitted = a["emitted"]
        out[coh] = {
            **a,
            "coverage_states": (emitted / a["n_states"]) if a["n_states"] else 0.0,
            "coverage_carriers": (a["carriers_emitted"] / a["carriers_total"]) if a["carriers_total"] else 0.0,
            "accuracy_among_emitted": (a["correct_emitted"] / emitted) if emitted else None,
        }
    return out


def _read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _demo_rows():
    # Small synthetic fixture spanning emitted, out_of_scope and no_call, three cohorts.
    return [
        {"cohort": "EUR", "gene": "CYP2D6", "diplotype": "*4/*4", "phenotype": "Poor Metabolizer", "n_carriers": 2},
        {"cohort": "EUR", "gene": "CYP2C19", "diplotype": "*2/*2", "phenotype": "Poor Metabolizer", "n_carriers": 1},
        {"cohort": "AMR", "gene": "CYP2D6", "diplotype": "*1/*99", "phenotype": "Normal Metabolizer", "n_carriers": 4},
        {"cohort": "AFR", "gene": "CYP2C19", "diplotype": "", "phenotype": "Indeterminate", "n_carriers": 3},
    ]


def _print_summary(summary):
    for coh, s in summary.items():
        acc = s["accuracy_among_emitted"]
        acc_s = "n/a" if acc is None else f"{100 * acc:.1f}%"
        print(f"[{coh}] states={s['n_states']} emitted={s['emitted']} "
              f"({100 * s['coverage_states']:.0f}% coverage; carriers {100 * s['coverage_carriers']:.0f}%) "
              f"accuracy-among-emitted={acc_s} "
              f"abstain(no_call={s['abstain']['no_call']}, ambiguous={s['abstain']['ambiguous']}, "
              f"out_of_scope={s['abstain']['out_of_scope']}) disagreements={len(s['disagreements'])}")


def main(argv):
    if not argv or argv[0] == "--demo":
        summary = run(_demo_rows(), SKILL_MAP, score_a1)
        _print_summary(summary)
        return 0
    if len(argv) < 2:
        print(__doc__)
        return 2
    tsv, out_json = argv[0], argv[1]
    rows = _read_tsv(tsv)
    summary = run(rows, SKILL_MAP, score_a1)
    Path(out_json).write_text(json.dumps(summary, indent=2))
    _print_summary(summary)
    print(f"-> {out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
