#!/usr/bin/env python3
"""
Tests for the frozen clinical-equivalence scorer (revision item N7, R2.5).

Reviewer 2 point 5: the clinical-equivalence patterns may have been widened
after model outputs were inspected, so the scorer could have been tuned to the
data. The remedy is mechanical rather than promissory: the pattern tables carry
a checksum, --frozen refuses to score if that checksum has moved, and every
headline number is reported under both this scorer and the untouched baseline.

These tests pin the three properties that make the freeze meaningful:
  1. the fingerprint covers the semantics of the tables (adding a pattern moves it)
  2. --frozen refuses to run once the fingerprint moves
  3. the equivalence layer can only promote a score, never demote one
"""
import re
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCORER = REPO / "code" / "10b-rescore-v3-clinical-equivalence.py"


def load_scorer():
    spec = spec_from_file_location("scorer_10b", SCORER)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def scorer():
    return load_scorer()


def test_fingerprint_is_deterministic(scorer):
    """Same tables, same hash: the freeze must not depend on dict ordering or run."""
    assert scorer.pattern_fingerprint() == scorer.pattern_fingerprint()
    assert re.fullmatch(r"[0-9a-f]{64}", scorer.pattern_fingerprint())


def test_declared_frozen_hash_matches_the_tables(scorer):
    """The hash recorded in the file (and in SCORING-PREREG.md) is the live one.

    This is the test that fails the moment anyone edits a pattern without
    re-registering the change, which is the whole point of the freeze.
    """
    assert scorer.FROZEN_PATTERN_SHA256 == scorer.pattern_fingerprint(), (
        "equivalence patterns changed without updating FROZEN_PATTERN_SHA256; "
        "if the change is intended, re-register it in SCORING-PREREG.md with a "
        "dated entry saying what changed and why"
    )


def test_fingerprint_moves_when_a_pattern_is_added(scorer, monkeypatch):
    """A widened table must not be able to masquerade as the frozen one."""
    before = scorer.pattern_fingerprint()
    widened = dict(scorer.CLINICAL_EQUIVALENCE)
    widened["CYP2D6"] = list(widened.get("CYP2D6", [])) + [
        (re.compile(r"anything at all"), "NORMAL", "test-only")
    ]
    monkeypatch.setattr(scorer, "CLINICAL_EQUIVALENCE", widened)
    assert scorer.pattern_fingerprint() != before


def test_fingerprint_moves_when_a_tier_label_changes(scorer, monkeypatch):
    """Not just the regex: the tier a pattern maps to is part of the contract."""
    before = scorer.pattern_fingerprint()
    generic = [(p, "SOMETHING_ELSE", tag) for (p, tier, tag) in scorer.GENERIC_HLA_RISK_ALLELE_PATTERNS]
    monkeypatch.setattr(scorer, "GENERIC_HLA_RISK_ALLELE_PATTERNS", generic)
    assert scorer.pattern_fingerprint() != before


def test_verify_frozen_passes_on_unmodified_tables(scorer):
    scorer.verify_frozen()  # must not raise


def test_verify_frozen_refuses_on_modified_tables(scorer, monkeypatch):
    monkeypatch.setattr(scorer, "FROZEN_PATTERN_SHA256", "0" * 64)
    with pytest.raises(SystemExit):
        scorer.verify_frozen()


def test_equivalence_layer_only_promotes(scorer):
    """The stated safety property, encoded as a test rather than a docstring claim.

    score_a1_clinical_eq must never return less than the baseline scorer, so the
    equivalence layer cannot manufacture a 1 -> 0 transition in either direction.
    """
    baseline = scorer.rescore.score_a1
    cases = [
        ("Poor Metaboliser", "Poor Metaboliser", "CYP2D6"),
        ("Ultra-rapid Metaboliser", "Poor Metaboliser", "CYP2D6"),
        ("positive for HLA-B*57:01", "HLA-B*57:01 positive", "HLA-B*57:01"),
        ("normal risk of hypersensitivity", "HLA-B*57:01 negative", "HLA-B*57:01"),
        ("", "Poor Metaboliser", "CYP2D6"),
        ("complete gibberish tokens", "Normal Metaboliser", "TPMT"),
    ]
    for parsed, gt, gene in cases:
        eq = scorer.score_a1_clinical_eq(parsed, gt, gene)
        base = baseline(parsed, gt)
        assert eq >= base, f"equivalence layer demoted {parsed!r} vs {gt!r}: {eq} < {base}"
