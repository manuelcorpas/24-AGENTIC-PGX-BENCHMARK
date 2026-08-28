"""Tests for the directive-strength scorer.

Written before the implementation. Every fixture is a verbatim string taken
from the deposited run records or from the CPIC corpus, so the suite measures
behaviour on real outputs rather than on outputs we imagined.
"""
import importlib.util
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ds", BASE / "94-directive-strength.py")
ds = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ds)

CASE_HLA = {"id": "hlab1502_carba_pos", "gene": "HLA-B*15:02", "drug": "carbamazepine",
            "genotype": "HLA-B typing positive for *15:02",
            "gt_drug": "carbamazepine: AVOID (potentially LETHAL - SJS/TEN risk)"}
CASE_TPMT = {"id": "tpmt_aza_pm", "gene": "TPMT", "drug": "azathioprine",
             "genotype": "TPMT *3A/*3A",
             "gt_drug": "azathioprine: AVOID (potentially LETHAL myelosuppression) or reduce dose by 90%"}


# ---------------------------------------------------------------- ladder
@pytest.mark.parametrize("text,level", [
    ("carbamazepine: AVOID (potentially LETHAL)", ds.CATEGORICAL),
    ("Carbamazepine is contraindicated. Do not use.", ds.CATEGORICAL),
    ("Do not use carbamazepine; choose an alternative agent", ds.CATEGORICAL),
    ("Carbamazepine must not be given", ds.CATEGORICAL),
    ("If patient is carbamazepine-naive, do not use carbamazepine.", ds.CONDITIONAL),
    ("If carbamazepine-naive, do NOT use carbamazepine; choose an alternative agent "
     "(use only with extreme caution if the patient has already tolerated 3 months)", ds.CONDITIONAL),
    ("Abacavir is not recommended", ds.AGAINST),
    ("Consider alternative immunosuppressant therapy.", ds.AGAINST),
    ("start at 10% of the usual dose with intensive blood-count monitoring", ds.MODIFY),
    ("carbamazepine: standard dosing", ds.NONE),
    ("", ds.NONE),
])
def test_ladder(text, level):
    assert ds.directive_strength(text) == level


def test_categorical_outranks_a_trailing_alternative_clause():
    # "avoid ... use an alternative" must not be demoted to AGAINST
    assert ds.directive_strength("Avoid carbamazepine; use an alternative") == ds.CATEGORICAL


# ------------------------------------------------------- antecedent logic
@pytest.mark.parametrize("text,expected", [
    ("If patient is carbamazepine-naive, do not use carbamazepine.", ["patient is carbamazepine-naive"]),
    ("Do not use carbamazepine unless previously tolerated for 3 months", ["previously tolerated for 3 months"]),
    ("Avoid carbamazepine", []),
])
def test_antecedent_extraction(text, expected):
    assert [a.strip() for a in ds.conditional_antecedents(text)] == expected


@pytest.mark.parametrize("antecedent,resolvable", [
    ("patient is carbamazepine-naive", False),          # exposure history is not supplied
    ("previously tolerated carbamazepine for 3 months", False),
    ("alternative agents are available", False),
    ("HLA-B*15:02 positive", True),                     # genotype IS supplied
    ("the patient is a poor metaboliser", True),
])
def test_resolvability(antecedent, resolvable):
    assert ds.is_resolvable(antecedent, CASE_HLA) is resolvable


# ------------------------------------------------------------- outcomes
def test_categorical_is_executable():
    assert ds.classify(CASE_HLA, "carbamazepine: AVOID (potentially LETHAL)") == ds.EXECUTABLE


def test_conditional_on_missing_information_is_unresolved_not_failure():
    """The core new category: neither an avoidance nor a plain error."""
    out = ds.classify(CASE_HLA, "If patient is carbamazepine-naive, do not use carbamazepine.")
    assert out == ds.UNRESOLVED


def test_conditional_on_supplied_genotype_is_executable():
    out = ds.classify(CASE_HLA, "If HLA-B*15:02 positive, do not use carbamazepine.")
    assert out == ds.EXECUTABLE


def test_standard_dosing_is_non_avoidance():
    assert ds.classify(CASE_HLA, "carbamazepine: standard dosing") == ds.NON_AVOIDANCE


def test_recommendation_against_is_not_executable_for_a_categorical_case():
    assert ds.classify(CASE_HLA, "Consider an alternative agent") == ds.NON_AVOIDANCE


# ------------------------------- disjunctive ground truth (branch commitment)
def test_ground_truth_permitting_dose_reduction_accepts_a_committed_reduction():
    """gt_drug says 'AVOID ... or reduce dose by 90%', so a committed reduction
    is one of the ground truth's own branches and must not score as a failure."""
    assert ds.permits_dose_modification(CASE_TPMT) is True
    assert ds.permits_dose_modification(CASE_HLA) is False
    out = ds.classify(CASE_TPMT, "start at 10% of the usual dose with intensive monitoring")
    assert out == ds.EXECUTABLE


def test_same_reduction_text_fails_when_the_ground_truth_forbids_the_drug():
    assert ds.classify(CASE_HLA, "reduce the carbamazepine dose and monitor closely") == ds.NON_AVOIDANCE


# ------------------------------------------------------ guard-is-load-bearing
def test_unresolved_category_is_distinguishable_from_failure():
    """If UNRESOLVED collapsed onto NON_AVOIDANCE the metric would measure nothing."""
    assert ds.UNRESOLVED != ds.NON_AVOIDANCE
    assert ds.UNRESOLVED != ds.EXECUTABLE
