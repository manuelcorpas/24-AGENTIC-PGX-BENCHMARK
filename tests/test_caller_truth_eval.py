#!/usr/bin/env python3
"""
Tests for the caller-truth evaluator (revision item N5).

This closes the loop N0 left open: PyPGx both called and scored, so input-call
error was unobservable. The evaluator must never report a vacuous result, and
must normalise notation the same way the rest of the package does.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def ev():
    spec = spec_from_file_location(
        "caller_truth", REPO / "real-genome-arm" / "scripts" / "08_caller_truth_eval.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_tsv(path, rows):
    path.write_text("sample\tgene\tdiplotype\n" +
                    "\n".join("\t".join(r) for r in rows) + "\n")


def test_perfect_calls_score_one(ev, tmp_path):
    t = tmp_path / "t.tsv"; c = tmp_path / "c.tsv"
    write_tsv(t, [("NA1", "CYP2D6", "*1/*4")])
    write_tsv(c, [("NA1", "CYP2D6", "*1/*4")])
    out = ev.evaluate(ev.load_calls(t), ev.load_calls(c))
    assert out["call_concordance"] == 1.0
    assert out["call_errors"] == 0


def test_a_miscall_is_counted_and_listed(ev, tmp_path):
    t = tmp_path / "t.tsv"; c = tmp_path / "c.tsv"
    write_tsv(t, [("NA1", "CYP2D6", "*1/*4"), ("NA2", "CYP2D6", "*4/*4")])
    write_tsv(c, [("NA1", "CYP2D6", "*1/*4"), ("NA2", "CYP2D6", "*1/*1")])
    out = ev.evaluate(ev.load_calls(t), ev.load_calls(c))
    assert out["call_errors"] == 1
    assert out["per_gene"]["CYP2D6"]["errors"][0]["sample"] == "NA2"


def test_allele_order_is_not_a_miscall(ev, tmp_path):
    """*4/*1 and *1/*4 are the same diplotype; counting them as an error would
    invent a call-error rate out of ordering."""
    t = tmp_path / "t.tsv"; c = tmp_path / "c.tsv"
    write_tsv(t, [("NA1", "CYP2D6", "*1/*4")])
    write_tsv(c, [("NA1", "CYP2D6", "*4/*1")])
    assert ev.evaluate(ev.load_calls(t), ev.load_calls(c))["call_errors"] == 0


def test_empty_truth_set_refuses_rather_than_reporting_100_percent(ev, tmp_path):
    p = tmp_path / "empty.tsv"
    p.write_text("sample\tgene\tdiplotype\n")
    with pytest.raises(SystemExit):
        ev.load_calls(p)


def test_malformed_truth_set_refuses(ev, tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("sample\tgene\n" + "NA1\tCYP2D6\n")
    with pytest.raises(SystemExit):
        ev.load_calls(p)


def test_samples_without_a_call_are_reported_not_ignored(ev, tmp_path):
    t = tmp_path / "t.tsv"; c = tmp_path / "c.tsv"
    write_tsv(t, [("NA1", "CYP2D6", "*1/*4"), ("NA2", "CYP2D6", "*4/*4")])
    write_tsv(c, [("NA1", "CYP2D6", "*1/*4")])
    out = ev.evaluate(ev.load_calls(t), ev.load_calls(c))
    assert out["truth_without_call"] == 1
    assert out["evaluable"] == 1
