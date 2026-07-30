"""Tests for 75-evaluate-input-normalisation.py.

The evaluation answers R1.1: can a model do the input-normalisation job the
paper argues belongs to it? Three properties matter enough to test.

1. The two references are not interchangeable. PyPGx is the deterministic
   caller the paper compares against; GeT-RM is external consensus. The paper
   already reports they agree only 0.761 of the time, so an evaluation that
   silently pooled them would report a number belonging to neither.
2. Accuracy must stay conditional on emission, per gene and per form as well as
   overall, or a model that abstains on the hard genes looks better than one
   that answers them.
3. A cell with no emitted calls has no accuracy. It must report None, not 0.0
   and not 1.0. Both defaults have appeared in this project's CORRECTIONS.md.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent


def _load(name, path):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ev = _load("evaluate_normalisation", BASE / "code" / "75-evaluate-input-normalisation.py")


ROWS = [
    # agrees with both references
    {"sample": "A", "gene": "CYP2D6", "form": "vcf", "model": "M1",
     "call": "*1/*1", "reference": "*1/*1", "getrm": "*1/*1"},
    # agrees with the caller, disagrees with GeT-RM
    {"sample": "B", "gene": "CYP2D6", "form": "vcf", "model": "M1",
     "call": "*4/*4", "reference": "*4/*4", "getrm": "*1/*4"},
    # abstained
    {"sample": "C", "gene": "CYP2D6", "form": "vcf", "model": "M1",
     "call": None, "reference": "*1/*1", "getrm": "*1/*1"},
    # wrong against both
    {"sample": "D", "gene": "TPMT", "form": "prose", "model": "M2",
     "call": "*1/*1", "reference": "*1/*3A", "getrm": "*1/*3A"},
]


def test_reports_both_references_separately():
    out = ev.evaluate(ROWS)
    assert out["overall"]["vs_caller"]["accuracy_among_emitted"] == pytest.approx(2 / 3)
    assert out["overall"]["vs_getrm"]["accuracy_among_emitted"] == pytest.approx(1 / 3)


def test_references_are_never_pooled():
    """A single blended accuracy would belong to neither reference."""
    out = ev.evaluate(ROWS)
    assert "accuracy" not in out["overall"]
    assert set(out["overall"]) >= {"vs_caller", "vs_getrm"}


def test_coverage_is_reference_independent():
    out = ev.evaluate(ROWS)
    assert out["overall"]["vs_caller"]["coverage"] == pytest.approx(0.75)
    assert out["overall"]["vs_getrm"]["coverage"] == pytest.approx(0.75)


def test_breaks_down_by_model_form_and_gene():
    out = ev.evaluate(ROWS)
    assert set(out["by_model"]) == {"M1", "M2"}
    assert set(out["by_form"]) == {"vcf", "prose"}
    assert set(out["by_gene"]) == {"CYP2D6", "TPMT"}


def test_empty_cell_has_no_accuracy():
    rows = [{"sample": "Z", "gene": "CYP2D6", "form": "vcf", "model": "M9",
             "call": None, "reference": "*1/*1", "getrm": "*1/*1"}]
    out = ev.evaluate(rows)
    assert out["overall"]["vs_caller"]["accuracy_among_emitted"] is None
    assert out["overall"]["vs_caller"]["coverage"] == 0.0


def test_rows_missing_getrm_are_excluded_from_the_getrm_arm_only():
    rows = ROWS + [{"sample": "E", "gene": "CYP2D6", "form": "vcf", "model": "M1",
                    "call": "*1/*1", "reference": "*1/*1", "getrm": None}]
    out = ev.evaluate(rows)
    assert out["overall"]["vs_caller"]["n"] == 5
    assert out["overall"]["vs_getrm"]["n"] == 4


def test_errored_rows_are_counted_as_abstentions_not_dropped():
    """Dropping an API failure shrinks the denominator and inflates coverage."""
    rows = ROWS + [{"sample": "F", "gene": "CYP2D6", "form": "vcf", "model": "M1",
                    "call": None, "reference": "*1/*1", "getrm": "*1/*1",
                    "error": "APIError: 500"}]
    out = ev.evaluate(rows)
    assert out["overall"]["vs_caller"]["n"] == 5
    assert out["overall"]["vs_caller"]["emitted"] == 3
    assert out["errors"] == 1
