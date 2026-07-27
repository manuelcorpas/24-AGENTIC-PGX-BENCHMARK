#!/usr/bin/env python3
"""
Tests for the extracted-rule evaluator (reviewer point 2).

The reviewer's objection was that the manuscript claimed extraction matches
authorship without ever extracting anything: both execution cells were handed the
authored vocabulary and ran the same executor, so their agreement was guaranteed.
The evaluator under test closes that gap by diffing a model-extracted table
against the authored one and executing the extracted table.

What these tests pin is that the comparison cannot flatter itself:
  - an extracted table that omits a diplotype must abstain, never fall back to
    the authored mapping, because a silent fallback would measure the authored
    table twice and reproduce exactly the circularity being corrected;
  - abstention is counted separately from a wrong answer;
  - the diff reports both directions, missing and extra.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "code" / "68-evaluate-extracted-rules.py"


@pytest.fixture(scope="module")
def ev():
    spec = spec_from_file_location("eval_extracted", SCRIPT)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AUTHORED = {
    "CYP2C19": {"*1/*1": "Normal Metabolizer",
                "*2/*2": "Poor Metabolizer",
                "*17/*17": "Ultrarapid Metabolizer"},
}


def table(d2p, recs=None):
    return {"diplotype_to_phenotype": d2p, "recommendations": recs or []}


# ---------------------------------------------------------------- diff

def test_diff_reports_missing_and_extra_in_both_directions(ev):
    extracted = table({"*1/*1": "Normal Metabolizer", "*3/*3": "Poor Metabolizer"})
    d = ev.diff_gene("CYP2C19", AUTHORED["CYP2C19"], extracted)
    assert d["missing"] == ["*2/*2", "*17/*17"] or set(d["missing"]) == {"*2/*2", "*17/*17"}
    assert set(d["extra"]) == {"*3/*3"}
    assert d["shared"] == 1


def test_diff_counts_phenotype_disagreement_on_shared_diplotypes(ev):
    extracted = table({"*1/*1": "Poor Metabolizer"})
    d = ev.diff_gene("CYP2C19", AUTHORED["CYP2C19"], extracted)
    assert d["shared"] == 1
    assert d["agree"] == 0
    assert d["disagree_examples"]


def test_phenotype_comparison_tolerates_wording_not_meaning(ev):
    """Guidelines say 'CYP2C19 normal metabolizer'; the authored table says
    'Normal Metabolizer'. Those are the same call and must not be scored as a
    disagreement, or the diff measures capitalisation."""
    extracted = table({"*1/*1": "CYP2C19 normal metabolizer"})
    d = ev.diff_gene("CYP2C19", AUTHORED["CYP2C19"], extracted)
    assert d["agree"] == 1


# ---------------------------------------------------------------- execution

def test_missing_diplotype_abstains_and_never_falls_back(ev):
    """A silent fallback to the authored mapping would measure the authored
    table twice and rebuild the circularity this experiment exists to remove."""
    extracted = table({"*1/*1": "Normal Metabolizer"})
    out = ev.execute_with("CYP2C19", "*2/*2", extracted)
    assert out is None


def test_present_diplotype_executes_from_the_extracted_table(ev):
    extracted = table({"*2/*2": "Poor Metabolizer"})
    assert ev.execute_with("CYP2C19", "*2/*2", extracted) == "Poor Metabolizer"


def test_scoring_separates_abstention_from_error(ev):
    rows = [
        {"gene": "CYP2C19", "called_diplotype": "*1/*1", "gt_phenotype": "Normal Metabolizer"},
        {"gene": "CYP2C19", "called_diplotype": "*2/*2", "gt_phenotype": "Poor Metabolizer"},
        {"gene": "CYP2C19", "called_diplotype": "*17/*17", "gt_phenotype": "Ultrarapid Metabolizer"},
    ]
    extracted = table({"*1/*1": "Normal Metabolizer", "*2/*2": "Normal Metabolizer"})
    r = ev.score_rows(rows, {"CYP2C19": extracted})
    assert r["n"] == 3
    assert r["correct"] == 1        # *1/*1
    assert r["wrong"] == 1          # *2/*2 mapped to the wrong phenotype
    assert r["abstained"] == 1      # *17/*17 absent from the extracted table
    assert r["accuracy_among_emitted"] == pytest.approx(0.5)


def test_empty_extracted_table_abstains_on_everything_rather_than_scoring_zero(ev):
    rows = [{"gene": "CYP2C19", "called_diplotype": "*1/*1",
             "gt_phenotype": "Normal Metabolizer"}]
    r = ev.score_rows(rows, {"CYP2C19": table({})})
    assert r["abstained"] == 1
    assert r["accuracy_among_emitted"] is None
