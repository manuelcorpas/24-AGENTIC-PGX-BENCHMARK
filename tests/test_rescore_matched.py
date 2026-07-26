#!/usr/bin/env python3
"""
Tests for the matched scorer (revision items N1, N3, N6; R2.1, R2.2, R2.5).

One scorer judges all five cells. Any per-cell branch in the scoring path would
reintroduce the confound the matched design exists to remove, so these tests
pin: drug-specificity (an answer about the wrong drug is wrong), dual scoring
under baseline and frozen equivalence, parse failures counted as failures
rather than dropped, lethal-class reported directly, and no aggregate A3.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


def _load(name):
    spec = spec_from_file_location(name.replace("-", "_"), REPO / "code" / name)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def rs():
    return _load("61-rescore-matched.py")


@pytest.fixture(scope="module")
def case(rs):
    return next(c for c in rs.CASES if c["id"] == "cyp2d6_codeine_pm")


@pytest.fixture(scope="module")
def lethal_case(rs):
    return next(c for c in rs.CASES if c["id"] == "cyp2d6_codeine_um")


def row(case, **kw):
    base = {
        "cell": "free_generation", "case_id": case["id"], "gene": case["gene"],
        "drug": case["drug"], "model": "m", "rep": 0, "raw": "x",
        "parsed_phenotype": case["gt_phenotype"], "parsed_drug": case["gt_drug"],
        "abstained": False, "error": None,
    }
    base.update(kw)
    return base


# ---------------------------------------------------------------- drug specificity

def test_correct_drug_answer_scores_full_marks(rs, case):
    s = rs.score_row(row(case), case)
    assert s["a1_phenotype"] == 1.0
    assert s["drug_match"] == 1.0


def test_answer_about_a_different_drug_is_scored_wrong(rs, case):
    """Drug substitution is the failure the original scoring could not see.

    A fluent, clinically sensible recommendation for the WRONG drug is not a
    partially correct answer; for this patient it is an incorrect one.
    """
    s = rs.score_row(row(case, parsed_drug="warfarin: reduce dose by 50%"), case)
    assert s["drug_match"] == 0.0
    assert s["a2_recommendation"] == 0.0


def test_terse_correct_answers_are_not_scored_as_substitution(rs, case):
    """A correct recommendation that does not restate the drug name is not a
    substitution. Requiring the name measured output style and penalised prose
    cells against execution cells, which emit canonical rule text."""
    s = rs.score_row(row(case, parsed_drug="Use label recommended age- or weight-specific dosing."), case)
    assert s["substituted"] is False
    assert s["drug_match"] == 1.0


def test_empty_answer_is_not_a_substitution(rs, case):
    assert rs.score_row(row(case, parsed_drug=""), case)["substituted"] is False


def test_scoring_has_no_per_cell_branch(rs):
    """The scorer must not know which cell it is scoring: that is what 'matched' means."""
    import inspect
    src = inspect.getsource(rs.score_row)
    for cell in ("free_generation", "rag_generation", "rag_execution",
                 "skill_generation", "skill_execution"):
        assert cell not in src


# ---------------------------------------------------------------- failures are failures

def test_empty_output_is_a_failure_not_a_dropped_row(rs, case):
    s = rs.score_row(row(case, parsed_phenotype="", parsed_drug="", raw=""), case)
    assert s["parsed_ok"] is False
    assert s["a1_phenotype"] == 0.0
    assert s["a2_recommendation"] == 0.0


def test_api_error_is_a_failure(rs, case):
    s = rs.score_row(row(case, raw="", parsed_phenotype="", parsed_drug="",
                         error="RateLimitError: 429"), case)
    assert s["parsed_ok"] is False


def test_abstention_is_distinct_from_a_wrong_answer(rs, case):
    """Abstaining is not correct, but it is not a dangerous error either; the
    safety analysis must be able to tell them apart."""
    s = rs.score_row(row(case, abstained=True, parsed_phenotype="", parsed_drug=""), case)
    assert s["abstained"] is True
    assert s["a1_phenotype"] == 0.0


# ---------------------------------------------------------------- lethal class

def test_lethal_case_is_flagged_from_the_guideline_text(rs, lethal_case, case):
    assert rs.is_lethal(lethal_case) is True
    assert rs.is_lethal(case) is False


def test_lethal_action_scored_only_on_lethal_cases(rs, case, lethal_case):
    assert rs.score_row(row(case), case)["lethal_action"] is None
    assert rs.score_row(row(lethal_case), lethal_case)["lethal_action"] == 1.0


def test_failing_to_avoid_on_a_lethal_case_is_a_lethal_error(rs, lethal_case):
    s = rs.score_row(row(lethal_case, parsed_drug="codeine: standard dosing"), lethal_case)
    assert s["lethal_action"] == 0.0


def test_no_aggregate_a3_anywhere(rs, case):
    """R2.5: A3 is 1.0 by definition for non-lethal cases, so its aggregate
    measures case composition rather than performance. Retired."""
    s = rs.score_row(row(case), case)
    assert "a3" not in s
    summary = rs.summarise([row(case)], {case["id"]: case})
    assert not any("a3" in k.lower() for k in summary["free_generation"])


# ---------------------------------------------------------------- dual scoring + freeze

def test_summary_reports_coverage_and_parse_failure(rs, case):
    rows = [row(case), row(case, parsed_phenotype="", parsed_drug="", raw=""),
            row(case, abstained=True, parsed_phenotype="", parsed_drug="")]
    s = rs.summarise(rows, {case["id"]: case})["free_generation"]
    assert s["n"] == 3
    assert s["parse_failures"] == 1
    assert s["abstentions"] == 1
    assert 0.0 <= s["coverage"] <= 1.0


def test_equivalence_scoring_requires_the_frozen_patterns(rs):
    """Scoring under the equivalence layer must run the freeze check, so a
    silently widened table cannot produce a reported number."""
    import inspect
    src = inspect.getsource(rs.load_equivalence_scorer)
    assert "verify_frozen" in src


def test_dual_scoring_reports_both_scorers(rs, case):
    out = rs.score_both([row(case)], {case["id"]: case})
    assert set(out) == {"baseline", "clinical_equivalence"}
