#!/usr/bin/env python3
"""
Ancestry reanalysis with matched diplotype states (revision item N5;
Reviewer 2 point 4).

WHAT THE REVIEWER SAID
R2.4: the real-genome experiment does not identify an ancestry effect. Three
cohorts of wildly different size, carrying different diplotype distributions,
were compared directly, so any cohort difference confounds ancestry with
composition and with cohort size. The paper's "degrades with ancestry" language
was not supportable, and has already been removed from the working manuscript.

WHAT THIS DOES
1. Restricts to diplotype states OBSERVED IN MORE THAN ONE COHORT, so the
   comparison is like-for-like rather than a comparison of different states.
2. Reports coverage on BOTH denominators, distinct states and carrier-weighted,
   because the steep state gradient is largely the rare-state tail exposed by a
   larger cohort. Applying the reviewer's own denominator correction to the
   coverage report is the point of the exercise.
3. Stratifies by known versus uncertain allele function, since "Indeterminate"
   and "Uncertain Susceptibility" states are a different phenomenon from
   genuine coverage gaps.
4. Decomposes the cohort difference into a composition part and a residual
   cohort part, by reweighting each cohort to a common state distribution
   (direct standardisation). If the difference vanishes under standardisation,
   it was composition, not ancestry.

WHAT IT DELIBERATELY DOES NOT DO
It does not fit a mixed model with a cohort term and call the coefficient an
ancestry effect. Three cohorts is not a sample of ancestries; it is three
cohorts. The standardisation above answers the reviewer's question without
dressing a descriptive comparison as an inferential one.

USAGE
    python code/64-ancestry-matched-analysis.py \
        --input real-genome-arm/n0/n0_input_3cohorts.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
IN_DEFAULT = BASE / "real-genome-arm" / "n0" / "n0_input_3cohorts.tsv"
OUT_DEFAULT = BASE / "data" / "v3_ancestry_matched.json"
REPORT_DEFAULT = BASE / "data" / "v3_ancestry_matched.txt"

# States whose phenotype is itself an expression of uncertainty rather than a
# call. Counting these as coverage failures would conflate "the allele tables
# do not reach this state" with "this state has no agreed clinical meaning".
UNCERTAIN_PHENOTYPES = {
    "indeterminate", "uncertain susceptibility", "unknown", "n/a", "no result",
    "uncertain function", "not available",
}


def _load_n0():
    """The N0 pipeline module, for its diplotype normalisation and skill vocabulary."""
    from importlib.util import spec_from_file_location, module_from_spec
    path = BASE / "real-genome-arm" / "scripts" / "07_executed_pipeline_n0.py"
    if not path.exists():
        return None
    spec = spec_from_file_location("n0_pipeline", path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_states(path: Path) -> list[dict]:
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["n_carriers"] = int(r.get("n_carriers") or 0)
    return rows


def is_uncertain(phenotype: str) -> bool:
    return (phenotype or "").strip().lower() in UNCERTAIN_PHENOTYPES


def matched_states(rows: list[dict], min_cohorts: int = 2) -> set[tuple[str, str]]:
    """(gene, diplotype) states observed in at least `min_cohorts` cohorts.

    This is the like-for-like set: comparing cohorts on states only one of them
    carries measures which states they carry, not how well anything performs.
    """
    seen: dict[tuple[str, str], set[str]] = defaultdict(set)
    for r in rows:
        seen[(r["gene"], r["diplotype"])].add(r["cohort"])
    return {state for state, cohorts in seen.items() if len(cohorts) >= min_cohorts}


def coverage_by_cohort(rows: list[dict], covered: set[tuple[str, str]],
                       restrict: set[tuple[str, str]] | None = None) -> dict:
    """Coverage on both denominators, per cohort."""
    out: dict[str, dict] = {}
    by_cohort: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        state = (r["gene"], r["diplotype"])
        if restrict is not None and state not in restrict:
            continue
        by_cohort[r["cohort"]].append(r)
    for cohort, rs in by_cohort.items():
        states = {(r["gene"], r["diplotype"]) for r in rs}
        n_states = len(states)
        n_covered = len(states & covered)
        carriers = sum(r["n_carriers"] for r in rs)
        carriers_covered = sum(r["n_carriers"] for r in rs
                               if (r["gene"], r["diplotype"]) in covered)
        uncertain = {(r["gene"], r["diplotype"]) for r in rs if is_uncertain(r["phenotype"])}
        out[cohort] = {
            "distinct_states": n_states,
            "states_covered": n_covered,
            "coverage_states": round(n_covered / n_states, 4) if n_states else None,
            "carriers": carriers,
            "carriers_covered": carriers_covered,
            "coverage_carriers": round(carriers_covered / carriers, 4) if carriers else None,
            "uncertain_function_states": len(uncertain),
            "coverage_states_excluding_uncertain": (
                round(len(states & covered) / len(states - uncertain), 4)
                if (states - uncertain) else None),
        }
    return out


def standardised_coverage(rows: list[dict], covered: set[tuple[str, str]]) -> dict:
    """Direct standardisation to a common state distribution.

    Each cohort is reweighted so that every cohort is evaluated on the SAME set
    of states with the SAME weights (the pooled carrier distribution over
    matched states). If the cohort differences disappear here, the raw gradient
    was composition rather than an ancestry effect, which is precisely the
    confound R2.4 identifies.
    """
    matched = matched_states(rows)
    pooled: dict[tuple[str, str], int] = defaultdict(int)
    for r in rows:
        state = (r["gene"], r["diplotype"])
        if state in matched:
            pooled[state] += r["n_carriers"]
    total_weight = sum(pooled.values())

    out: dict[str, dict] = {}
    cohorts = sorted({r["cohort"] for r in rows})
    for cohort in cohorts:
        present = {(r["gene"], r["diplotype"]) for r in rows if r["cohort"] == cohort}
        num = sum(w for state, w in pooled.items()
                  if state in present and state in covered)
        den = sum(w for state, w in pooled.items() if state in present)
        out[cohort] = {
            "standardised_coverage": round(num / den, 4) if den else None,
            "matched_states_present": len(present & matched),
            "weight_share": round(den / total_weight, 4) if total_weight else None,
        }
    return out


def cohort_specific_coverage(rows: list[dict], covered: set[tuple[str, str]]) -> dict:
    """Coverage on states UNIQUE to each cohort.

    This exists because matched-state standardisation is blind by construction
    to states only one cohort carries, and that is exactly where the coverage
    gap lives: a cohort's untypeable rare alleles are cohort-specific. Reporting
    only the standardised figure would announce that the disparity disappears
    when in fact the analysis had stopped looking at it.
    """
    matched = matched_states(rows)
    out: dict[str, dict] = {}
    by_cohort: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_cohort[r["cohort"]].append(r)
    for cohort, rs in by_cohort.items():
        specific = {(r["gene"], r["diplotype"]) for r in rs} - matched
        carriers = sum(r["n_carriers"] for r in rs
                       if (r["gene"], r["diplotype"]) in specific)
        carriers_covered = sum(r["n_carriers"] for r in rs
                               if (r["gene"], r["diplotype"]) in specific
                               and (r["gene"], r["diplotype"]) in covered)
        out[cohort] = {
            "cohort_specific_states": len(specific),
            "covered": len(specific & covered),
            "coverage_states": (round(len(specific & covered) / len(specific), 4)
                                if specific else None),
            "coverage_carriers": (round(carriers_covered / carriers, 4)
                                  if carriers else None),
        }
    return out


def analyse(rows: list[dict], covered: set[tuple[str, str]]) -> dict:
    matched = matched_states(rows)
    raw = coverage_by_cohort(rows, covered)
    restricted = coverage_by_cohort(rows, covered, restrict=matched)
    standardised = standardised_coverage(rows, covered)
    specific = cohort_specific_coverage(rows, covered)

    raw_vals = [v["coverage_states"] for v in raw.values() if v["coverage_states"] is not None]
    std_vals = [v["standardised_coverage"] for v in standardised.values()
                if v["standardised_coverage"] is not None]
    return {
        "n_states": len({(r["gene"], r["diplotype"]) for r in rows}),
        "n_matched_states": len(matched),
        "cohorts": sorted({r["cohort"] for r in rows}),
        "raw_by_cohort": raw,
        "matched_states_only": restricted,
        "standardised": standardised,
        "cohort_specific": specific,
        "raw_state_coverage_spread": (round(max(raw_vals) - min(raw_vals), 4)
                                      if raw_vals else None),
        "standardised_coverage_spread": (round(max(std_vals) - min(std_vals), 4)
                                         if std_vals else None),
        "cohort_specific_coverage_spread": _spread(
            [v["coverage_states"] for v in specific.values()]),
    }


def _spread(values) -> float | None:
    vals = [v for v in values if v is not None]
    return round(max(vals) - min(vals), 4) if len(vals) > 1 else None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=IN_DEFAULT)
    ap.add_argument("--covered", type=Path, default=None,
                    help="TSV of gene<TAB>diplotype states the validated skill covers; "
                         "defaults to the states with a non-uncertain phenotype")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--report", type=Path, default=REPORT_DEFAULT)
    args = ap.parse_args(argv)

    if not args.input.exists():
        sys.stderr.write(f"no input at {args.input}\n")
        return 1
    rows = load_states(args.input)

    # Coverage means "in the validated skill's vocabulary", exactly as N0 defines
    # it. Diplotypes are normalised with the same function N0 uses, because a
    # cohort writes "*1/*1xN" where the skill stores "*1/*1xn"; matching raw
    # strings would report a vocabulary miss that is really a spelling
    # difference, and would inflate the coverage gap this analysis exists to
    # measure honestly.
    n0 = _load_n0()
    if args.covered and args.covered.exists():
        with args.covered.open() as fh:
            covered = {(a.strip(), n0.norm_dip(b.strip(), a.strip()))
                       for a, b, *_ in csv.reader(fh, delimiter="\t")}
    elif n0 is not None:
        covered = set(n0.SKILL_MAP)
    else:
        covered = {(r["gene"], r["diplotype"]) for r in rows if not is_uncertain(r["phenotype"])}

    if n0 is not None:
        for r in rows:
            r["diplotype"] = n0.norm_dip(r["diplotype"], r["gene"])

    result = analyse(rows, covered)
    args.out.write_text(json.dumps(result, indent=2))

    lines = ["ANCESTRY REANALYSIS ON MATCHED DIPLOTYPE STATES", "",
             f"  distinct states {result['n_states']}, "
             f"observed in >1 cohort: {result['n_matched_states']}", ""]
    lines.append(f"  {'cohort':<16}{'states':>8}{'cov(state)':>12}{'cov(carrier)':>14}"
                 f"{'standardised':>14}{'own-states':>12}")
    for cohort in result["cohorts"]:
        r = result["raw_by_cohort"][cohort]
        st = result["standardised"][cohort]
        sp = result["cohort_specific"][cohort]
        lines.append(f"  {cohort:<16}{r['distinct_states']:>8}"
                     f"{str(r['coverage_states']):>12}{str(r['coverage_carriers']):>14}"
                     f"{str(st['standardised_coverage']):>14}"
                     f"{str(sp['coverage_states']):>12}")
    lines += [
        "",
        f"  spread in raw state coverage       {result['raw_state_coverage_spread']}",
        f"  spread after standardisation       {result['standardised_coverage_spread']}",
        f"  spread on cohort-specific states   {result['cohort_specific_coverage_spread']}",
        "",
        "  Standardisation compares cohorts on states they SHARE, so a small spread",
        "  there means composition explains the raw gradient on common ground. It is",
        "  blind to states only one cohort carries, which is where the untypeable",
        "  rare-allele tail sits; the last column reports those directly. Read the",
        "  two together: the disparity is relocated and measured, not dissolved.",
    ]
    text = "\n".join(lines)
    args.report.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
