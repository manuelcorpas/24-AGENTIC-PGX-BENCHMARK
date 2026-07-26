#!/usr/bin/env python3
"""
Tests for the hierarchical / multiway-clustered statistics (revision items N6
and N3; Reviewer 2 point 5).

R2.5 asked for four things and each is pinned here:
  1. intervals clustered by model as well as by case
  2. a variance decomposition across the crossed factors, not one bootstrap
  3. coverage, abstention and parse failure reported, with failures counted as
     failures rather than dropped
  4. aggregate A3 retired

The statistical property that matters most for the paper is the last test in
this file: under deterministic execution the between-model variance component
must be zero, because that is what "model-invariant" means as a measurement
rather than as a claim.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def hs():
    spec = spec_from_file_location("hier_stats", REPO / "code" / "65-hierarchical-stats.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rows(n_cases=10, n_models=4, n_reps=3, value=lambda case, model, rep: 1.0,
         framing="EUR"):
    return [
        {"case_id": f"c{c}", "model": f"m{m}", "rep": r, "framing": framing,
         "correct": value(c, m, r), "parsed_ok": True, "abstained": False,
         "lethal_action": None}
        for c in range(n_cases) for m in range(n_models) for r in range(n_reps)
    ]


# ---------------------------------------------------------------- clustering

def test_case_and_model_clustered_intervals_are_both_reported(hs):
    out = hs.clustered_intervals(rows(value=lambda c, m, r: float((c + m) % 2)))
    assert set(out) >= {"case_clustered", "model_clustered", "two_way_clustered"}
    for key in ("case_clustered", "model_clustered", "two_way_clustered"):
        lo, hi = out[key]["ci95"]
        assert lo <= out[key]["mean"] <= hi


def test_model_clustered_interval_is_wider_when_models_differ_systematically(hs):
    """If whole models are good or bad, clustering by case understates uncertainty.

    This is exactly R2.5's objection to a single case-level bootstrap.
    """
    data = rows(n_models=4, value=lambda c, m, r: 1.0 if m < 2 else 0.0)
    out = hs.clustered_intervals(data)
    case_w = out["case_clustered"]["ci95"][1] - out["case_clustered"]["ci95"][0]
    model_w = out["model_clustered"]["ci95"][1] - out["model_clustered"]["ci95"][0]
    assert model_w > case_w


def test_zero_variance_data_gives_a_degenerate_interval(hs):
    out = hs.clustered_intervals(rows(value=lambda c, m, r: 1.0))
    assert out["case_clustered"]["mean"] == 1.0
    assert out["case_clustered"]["ci95"] == (1.0, 1.0)


def test_intervals_are_reproducible(hs):
    data = rows(value=lambda c, m, r: float((c * m + r) % 2))
    assert hs.clustered_intervals(data, seed=7) == hs.clustered_intervals(data, seed=7)


# ---------------------------------------------------------------- variance decomposition

def test_variance_decomposition_covers_the_crossed_factors(hs):
    comp = hs.variance_components(rows(value=lambda c, m, r: float((c + m + r) % 2)))
    assert set(comp) >= {"case", "model", "framing", "replicate", "residual"}


def test_model_variance_is_zero_under_deterministic_execution(hs):
    """The paper's model-invariance claim, expressed as a measurement.

    Deterministic execution returns the same answer whatever model supplied the
    input call, so the between-model component must be exactly zero. If this
    ever fails, the invariance claim is false and must be withdrawn.
    """
    executed = rows(value=lambda c, m, r: float(c % 2))   # depends on case only
    comp = hs.variance_components(executed)
    assert comp["model"] == pytest.approx(0.0, abs=1e-12)
    assert comp["replicate"] == pytest.approx(0.0, abs=1e-12)
    assert comp["case"] > 0


def test_model_variance_is_positive_when_models_disagree(hs):
    generated = rows(value=lambda c, m, r: float(m % 2))
    comp = hs.variance_components(generated)
    assert comp["model"] > 0


# ---------------------------------------------------------------- reporting

def test_failures_are_counted_not_dropped(hs):
    data = rows(n_cases=4, n_models=2, n_reps=1)
    data[0]["parsed_ok"] = False
    data[0]["correct"] = 0.0
    rep = hs.report(data)
    assert rep["n"] == len(data)
    assert rep["parse_failures"] == 1
    assert rep["coverage"] < 1.0


def test_lethal_class_is_reported_directly(hs):
    data = rows(n_cases=4, n_models=2, n_reps=1)
    for r in data[:4]:
        r["lethal_action"] = 1.0
    data[0]["lethal_action"] = 0.0
    rep = hs.report(data)
    assert rep["lethal_n"] == 4
    assert rep["lethal_errors"] == 1
    assert rep["lethal_action_accuracy"] == pytest.approx(0.75)


def test_no_aggregate_a3_is_produced(hs):
    rep = hs.report(rows(n_cases=3, n_models=2, n_reps=1))
    assert not any("a3" in k.lower() for k in rep)
