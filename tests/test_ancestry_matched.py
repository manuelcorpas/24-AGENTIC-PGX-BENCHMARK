#!/usr/bin/env python3
"""
Tests for the ancestry reanalysis (revision item N5; Reviewer 2 point 4).

R2.4's objection is that three cohorts of different sizes carrying different
diplotype distributions cannot identify an ancestry effect. The analysis must
therefore be able to SHOW that a raw difference is composition, and must not
manufacture a difference where none exists.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def am():
    spec = spec_from_file_location("ancestry", REPO / "code" / "64-ancestry-matched-analysis.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def row(cohort, gene, dip, phen="Normal Metaboliser", n=1):
    return {"cohort": cohort, "gene": gene, "diplotype": dip,
            "phenotype": phen, "n_carriers": n}


def test_matched_states_require_more_than_one_cohort(am):
    rows = [row("A", "CYP2D6", "*1/*1"), row("B", "CYP2D6", "*1/*1"),
            row("A", "CYP2D6", "*4/*4")]
    assert am.matched_states(rows) == {("CYP2D6", "*1/*1")}


def test_coverage_reported_on_both_denominators(am):
    rows = [row("A", "CYP2D6", "*1/*1", n=100), row("A", "CYP2D6", "*9/*9", n=1)]
    cov = am.coverage_by_cohort(rows, covered={("CYP2D6", "*1/*1")})["A"]
    assert cov["coverage_states"] == pytest.approx(0.5)
    assert cov["coverage_carriers"] == pytest.approx(100 / 101, abs=1e-4)


def test_rare_state_tail_depresses_state_coverage_but_not_carrier_coverage(am):
    """The mechanism behind the apparent ancestry gradient, isolated.

    A larger cohort exposes a long tail of rare states. Counted per distinct
    state, that looks like collapsing coverage; counted per carrier, it barely
    moves. Reporting only the first denominator is what produced the original
    'degrades with ancestry' reading.
    """
    small = [row("small", "CYP2D6", "*1/*1", n=100)]
    large = ([row("large", "CYP2D6", "*1/*1", n=10000)]
             + [row("large", "CYP2D6", f"*{i}/*{i}", n=1) for i in range(20, 40)])
    cov = am.coverage_by_cohort(small + large, covered={("CYP2D6", "*1/*1")})
    assert cov["small"]["coverage_states"] == 1.0
    assert cov["large"]["coverage_states"] < 0.1
    assert cov["large"]["coverage_carriers"] > 0.99


def test_standardisation_removes_a_purely_compositional_difference(am):
    """If cohorts differ only in HOW MANY carriers sit on each shared state,
    standardising to a common distribution must equalise them."""
    rows = [row("A", "G", "x", n=90), row("A", "G", "y", n=10),
            row("B", "G", "x", n=10), row("B", "G", "y", n=90)]
    out = am.analyse(rows, covered={("G", "x")})
    a, b = (out["standardised"]["A"]["standardised_coverage"],
            out["standardised"]["B"]["standardised_coverage"])
    assert a == pytest.approx(b)
    assert out["raw_by_cohort"]["A"]["coverage_carriers"] != \
           out["raw_by_cohort"]["B"]["coverage_carriers"]


def test_standardisation_is_blind_to_cohort_specific_states(am):
    """The method's own limitation, pinned so it cannot be forgotten.

    Standardisation compares cohorts on states they SHARE. Cohort B below
    carries an uncovered state that A does not, but that state is unmatched, so
    standardisation cannot see it and reports the cohorts as equal. Reporting
    only the standardised figure would announce that the disparity vanished
    when the analysis had merely stopped looking at it.
    """
    rows = [row("A", "G", "x", n=50), row("A", "G", "y", n=50),
            row("B", "G", "x", n=50), row("B", "G", "z", n=50)]
    out = am.analyse(rows, covered={("G", "x"), ("G", "y")})
    assert out["standardised"]["A"]["standardised_coverage"] == \
           out["standardised"]["B"]["standardised_coverage"]


def test_cohort_specific_coverage_catches_what_standardisation_misses(am):
    """The companion metric that makes the pair honest.

    Same data as above: the difference invisible to standardisation must show
    up here, because that is where the untypeable rare-allele tail lives.
    """
    rows = [row("A", "G", "x", n=50), row("A", "G", "y", n=50),
            row("B", "G", "x", n=50), row("B", "G", "z", n=50)]
    out = am.analyse(rows, covered={("G", "x"), ("G", "y")})
    assert out["cohort_specific"]["A"]["coverage_states"] == 1.0
    assert out["cohort_specific"]["B"]["coverage_states"] == 0.0


def test_uncertain_phenotypes_are_identified(am):
    assert am.is_uncertain("Indeterminate")
    assert am.is_uncertain("Uncertain Susceptibility")
    assert not am.is_uncertain("Poor Metaboliser")


def test_real_cohort_file_runs_and_reports_all_three_cohorts(am):
    rows = am.load_states(REPO / "real-genome-arm" / "n0" / "n0_input_3cohorts.tsv")
    n0 = am._load_n0()
    for r in rows:
        r["diplotype"] = n0.norm_dip(r["diplotype"], r["gene"])
    out = am.analyse(rows, covered=set(n0.SKILL_MAP))
    assert set(out["cohorts"]) == {"CorpasFamily", "Peru", "UGR"}
    assert out["n_matched_states"] > 0
    # The finding: standardisation shrinks the raw gradient substantially.
    assert out["standardised_coverage_spread"] < out["raw_state_coverage_spread"]
