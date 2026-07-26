#!/usr/bin/env python3
"""
Tests for drug-keyed retrieval (revision item N4; R1.2, R2.3).

The strengthened retrieval arm must retrieve the queried drug's section rather
than the whole gene's annotation set, must never fall back silently, and must
declare its class-to-member mapping explicitly.
"""
import json
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


def _load(name, mod_name):
    spec = spec_from_file_location(mod_name, REPO / "code" / name)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def r():
    return _load("_pgx_rules.py", "_pgx_rules")


@pytest.fixture(scope="module")
def mf():
    return _load("60-matched-factorial.py", "matched_factorial")


@pytest.fixture(scope="module")
def corpus():
    return json.loads((REPO / "specs" / "cpic_rag_corpus_v3.json").read_text())["genes"]


def test_drug_chunk_returns_only_the_queried_drug(r, corpus):
    chunk, route = r.drug_chunk(corpus["CYP2D6"]["guideline_excerpt"], "codeine")
    assert route == "exact"
    assert "codeine" in chunk.lower()
    assert "tamoxifen" not in chunk.lower()
    assert "ondansetron" not in chunk.lower()


def test_drug_chunk_is_much_smaller_than_the_gene_excerpt(r, corpus):
    full = corpus["CYP2D6"]["guideline_excerpt"]
    chunk, _ = r.drug_chunk(full, "codeine")
    assert len(chunk) < len(full) / 2


def test_drug_class_queries_resolve_to_member_sections(r, corpus):
    chunk, route = r.drug_chunk(corpus["NUDT15"]["guideline_excerpt"], "thiopurines")
    assert route == "class"
    assert "azathioprine" in chunk.lower()


def test_unknown_drug_reports_none_rather_than_guessing(r, corpus):
    chunk, route = r.drug_chunk(corpus["CYP2D6"]["guideline_excerpt"], "not-a-real-drug")
    assert chunk is None
    assert route == "none"


def test_every_benchmark_case_resolves(r, corpus):
    """All 110 cases must retrieve something, by an explicitly named route."""
    unresolved = []
    for c in r.load_cases():
        excerpt = corpus.get(c["gene"], {}).get("guideline_excerpt", "")
        chunk, route = r.drug_chunk(excerpt, c["drug"])
        if chunk is None:
            unresolved.append((c["gene"], c["drug"]))
    assert not unresolved, f"cases with no drug-keyed chunk: {unresolved}"


def test_class_aliases_are_declared_not_inferred(r):
    """A fuzzy matcher here would silently change what the model was shown."""
    assert set(r.DRUG_CLASS_ALIASES) >= {"thiopurines", "aminoglycosides",
                                         "volatile-anaesthetics"}


def test_retrieval_route_is_recorded_for_every_retrieval_row(mf):
    case = next(c for c in mf.CASES if c["id"] == "cyp2d6_codeine_pm")
    previous = mf.RETRIEVAL_MODE
    mf.RETRIEVAL_MODE = "drug"
    try:
        row = mf.run_one("rag_generation", case, "fake",
                         lambda p: ("DIPLOTYPE: *4/*4\nPHENOTYPE: Poor Metaboliser\n"
                                    "DRUG: codeine: AVOID\nHAZARD: none", 10, 10), 0)
    finally:
        mf.RETRIEVAL_MODE = previous
    assert row["retrieval_mode"] == "drug"
    assert row["retrieval_route"].startswith("drug_keyed")


def test_non_retrieval_rows_carry_no_retrieval_route(mf):
    case = next(c for c in mf.CASES if c["id"] == "cyp2d6_codeine_pm")
    row = mf.run_one("free_generation", case, "fake", lambda p: ("x", 1, 1), 0)
    assert row["retrieval_route"] is None


def test_gene_keyed_mode_still_available_as_the_ablation(mf):
    case = next(c for c in mf.CASES if c["id"] == "cyp2d6_codeine_pm")
    previous = mf.RETRIEVAL_MODE
    mf.RETRIEVAL_MODE = "gene"
    try:
        excerpt, route = mf.retrieve(case["gene"], case["drug"])
    finally:
        mf.RETRIEVAL_MODE = previous
    assert route == "gene_keyed"
    assert "tamoxifen" in excerpt.lower()


def test_route_census_covers_all_cases(mf):
    driver = _load("62-rag-drugkeyed-fullgrid.py", "drugkeyed")
    census = driver.route_census()
    assert sum(census.values()) == len(mf.CASES)
    assert census.get("fallback_gene_keyed", 0) == 0
