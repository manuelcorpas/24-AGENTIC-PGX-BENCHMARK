#!/usr/bin/env python3
"""
One scorer for all five matched cells (revision items N1, N3, N6).

WHY ONE SCORER
The original design scored the skill arms drug-specifically and the other arms
not at all drug-specifically, so the arms were not comparable (R1.3, R2.1).
Here a single function scores every cell and does not know which cell it is
looking at; a test asserts no cell name appears in its source.

WHAT IS SCORED
  a1_phenotype       phenotype tier match, via the locked baseline scorer
  a2_recommendation  recommendation direction, conditional on the right drug
  drug_match         did the answer address the TARGET drug at all
  lethal_action      on lethal-class cases only: was the avoid/contraindicate
                     action taken (None elsewhere, never 1.0-by-default)
  parsed_ok          missing, empty and unparsed outputs are failures (R2.5)
  abstained          abstention is tracked separately from a wrong answer

WHAT IS NOT SCORED
Aggregate A3 is retired. It returns 1.0 by definition for the 96 non-lethal
cases, so its aggregate measures case composition rather than model behaviour
(R2.5). Lethal-class performance is reported directly instead.

DUAL SCORING
Every number is produced twice, under the locked baseline scorer and under the
frozen clinical-equivalence scorer, so that no conclusion can rest on the
equivalence layer. The equivalence path calls verify_frozen() first, so a
silently widened pattern table cannot produce a reported number
(see SCORING-PREREG.md).

USAGE
    python code/61-rescore-matched.py                 # score data/v3_matched_factorial.json
    python code/61-rescore-matched.py --input FILE --out FILE
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
IN_DEFAULT = BASE / "data" / "v3_matched_factorial.json"
OUT_DEFAULT = BASE / "data" / "v3_matched_factorial_scored.json"
REPORT_DEFAULT = BASE / "data" / "v3_matched_factorial_report.txt"


def _load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rules = _load(CODE / "_pgx_rules.py", "_pgx_rules")
baseline = _load(CODE / "10-rescore-v3.py", "rescore_v3")

CASES = rules.load_cases()
CASE_BY_ID = {c["id"]: c for c in CASES}

# Every drug named anywhere in the benchmark. Used to detect drug substitution:
# an answer that names a different benchmark drug is answering a different
# question, however fluent it is.
ALL_DRUGS = sorted({c["drug"].lower() for c in CASES})


def load_equivalence_scorer():
    """Import the clinical-equivalence scorer and refuse to proceed if its
    pattern tables have moved since registration (SCORING-PREREG.md)."""
    eq = _load(CODE / "10b-rescore-v3-clinical-equivalence.py", "rescore_eq")
    eq.verify_frozen()
    return eq


def is_lethal(case: dict) -> bool:
    """Lethal-class membership comes from the CPIC guideline text, not from a
    post-hoc judgment about which cases models found hard (SCORING-PREREG.md)."""
    return "LETHAL" in str(case.get("gt_drug", "")).upper()


def _substitutes_another_drug(text: str, target: str) -> bool:
    """True only when the answer is about a DIFFERENT drug than the one asked.

    Deliberately not "does the answer repeat the target drug name". Models
    frequently give a correct recommendation without restating the drug ("Use
    label recommended age- or weight-specific dosing"), while the execution
    cells emit canonical rule text that always contains it. Requiring the name
    therefore measured output style, not drug substitution, and penalised
    generation cells for prose that was correct. Substitution means naming a
    different benchmark drug while not naming the target.
    """
    low = (text or "").lower()
    if not low:
        return False
    names_target = bool(re.search(rf"\b{re.escape(target.lower())}\b", low))
    others = [d for d in ALL_DRUGS
              if d != target.lower() and re.search(rf"\b{re.escape(d)}\b", low)]
    return bool(others) and not names_target


def score_row(row: dict, case: dict, a1_scorer=None) -> dict:
    """Score one evaluation. Deliberately blind to which cell produced it."""
    a1_fn = a1_scorer or (lambda phen, gt, gene: baseline.score_a1(phen, gt))

    parsed_phen = (row.get("parsed_phenotype") or "").strip()
    parsed_drug = (row.get("parsed_drug") or "").strip()
    # An abstention is NOT a parse failure. The model produced a usable output
    # and the validated skill declined to map it, which is the designed
    # behaviour. Conflating the two would report deliberate caution as broken
    # infrastructure and would make the abstention argument unreadable.
    abstained = bool(row.get("abstained"))
    parsed_ok = (
        not row.get("error")
        and bool((row.get("raw") or "").strip())
        and (abstained or bool(parsed_phen or parsed_drug))
    )

    substituted = _substitutes_another_drug(parsed_drug, case["drug"])
    drug_match = 0.0 if substituted else 1.0

    a1 = a1_fn(parsed_phen, case["gt_phenotype"], case["gene"]) if parsed_phen else 0.0
    # Recommendation credit is conditional on having answered about the right
    # drug: a correct-sounding recommendation for another drug is not partial
    # credit for this patient.
    # Credit is withheld only for genuine substitution, not for terse phrasing.
    a2 = 0.0 if substituted else baseline.score_a2(parsed_drug, case["gt_drug"])

    lethal_action = None
    if is_lethal(case):
        avoided = bool(re.search(r"\bavoid\b|\bcontraindicat", parsed_drug, re.IGNORECASE))
        lethal_action = 1.0 if (avoided and not substituted) else 0.0

    return {
        "cell": row["cell"], "case_id": case["id"], "gene": case["gene"],
        "drug": case["drug"], "model": row.get("model"), "rep": row.get("rep"),
        "a1_phenotype": float(a1),
        "a2_recommendation": float(a2),
        "drug_match": drug_match, "substituted": substituted,
        "lethal_action": lethal_action,
        "parsed_ok": parsed_ok,
        "abstained": abstained,
    }


def summarise(rows: list[dict], case_by_id: dict | None = None, a1_scorer=None) -> dict:
    """Per-cell summary. Coverage, abstention and parse failure are reported
    beside accuracy, never conditioned away (R2.2, R2.5)."""
    cases = case_by_id or CASE_BY_ID
    by_cell: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_cell[r["cell"]].append(score_row(r, cases[r["case_id"]], a1_scorer))

    out = {}
    for cell, scored in by_cell.items():
        n = len(scored)
        emitted = [s for s in scored if s["parsed_ok"] and not s["abstained"]]
        lethal = [s for s in scored if s["lethal_action"] is not None]
        out[cell] = {
            "n": n,
            "parse_failures": sum(1 for s in scored if not s["parsed_ok"]),
            "abstentions": sum(1 for s in scored if s["abstained"]),
            "coverage": round(len(emitted) / n, 4) if n else 0.0,
            "a1_phenotype_mean": round(sum(s["a1_phenotype"] for s in scored) / n, 4) if n else 0.0,
            "a2_recommendation_mean": round(sum(s["a2_recommendation"] for s in scored) / n, 4) if n else 0.0,
            "no_substitution_rate": round(sum(s["drug_match"] for s in scored) / n, 4) if n else 0.0,
            "a1_among_emitted": round(sum(s["a1_phenotype"] for s in emitted) / len(emitted), 4) if emitted else None,
            "lethal_n": len(lethal),
            "lethal_action_accuracy": round(sum(s["lethal_action"] for s in lethal) / len(lethal), 4) if lethal else None,
            "lethal_errors": sum(1 for s in lethal if s["lethal_action"] == 0.0),
        }
    return out


def score_both(rows: list[dict], case_by_id: dict | None = None) -> dict:
    """Report under both scorers. If a conclusion holds under only one of them,
    it is not a conclusion (SCORING-PREREG.md section 1)."""
    eq = load_equivalence_scorer()
    return {
        "baseline": summarise(rows, case_by_id),
        "clinical_equivalence": summarise(
            rows, case_by_id,
            a1_scorer=lambda phen, gt, gene: eq.score_a1_clinical_eq(phen, gt, gene)),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=IN_DEFAULT)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--report", type=Path, default=REPORT_DEFAULT)
    args = ap.parse_args(argv)

    if not args.input.exists():
        sys.stderr.write(
            f"no input at {args.input}\n"
            "Run code/60-matched-factorial.py first (or pass --input).\n")
        return 1

    rows = json.loads(args.input.read_text())
    result = score_both(rows)
    args.out.write_text(json.dumps(result, indent=2))

    lines = ["MATCHED FACTORIAL: reported under both scorers", ""]
    for scorer, cells in result.items():
        lines.append(f"## {scorer}")
        for cell, s in sorted(cells.items()):
            lines.append(
                f"  {cell:<18} n={s['n']:<6} coverage={s['coverage']:<7} "
                f"A1={s['a1_phenotype_mean']:<7} A2={s['a2_recommendation_mean']:<7} "
                f"no_sub={s['no_substitution_rate']:<7} "
                f"lethal={s['lethal_action_accuracy']} ({s['lethal_errors']} errors) "
                f"parse_fail={s['parse_failures']}")
        lines.append("")
    report = "\n".join(lines)
    args.report.write_text(report)
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
