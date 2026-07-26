#!/usr/bin/env python3
"""
Tests for the PharmCAT comparator (revision item N2; R2.2, R2.7).

The comparator's danger is not that it crashes; it is that it reports notation
mismatches as scientific disagreement. Every translation below was established
by probing the pinned jar, and these tests stop them silently changing.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def pc():
    spec = spec_from_file_location("pharmcat_cmp", REPO / "code" / "63-pharmcat-comparator.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- normalisation

def test_us_and_uk_spelling_are_the_same_phenotype(pc):
    assert pc.normalise_phenotype("Intermediate Metabolizer") == \
           pc.normalise_phenotype("Intermediate Metaboliser")


def test_parenthetical_qualifiers_do_not_create_disagreement(pc):
    """'Normal Metaboliser (Expressor)' and 'Normal Metabolizer' are one call."""
    assert pc.normalise_phenotype("Normal Metaboliser (Expressor)") == \
           pc.normalise_phenotype("Normal Metabolizer")


def test_no_result_is_not_a_phenotype(pc):
    """Reading PharmCAT's 'No Result' as a call would fabricate disagreements."""
    assert pc.normalise_phenotype("No Result") == ""
    assert pc.normalise_phenotype("n/a") == ""


def test_hla_status_normalises_to_positive_negative(pc):
    assert pc.normalise_phenotype("*57:01 positive") == "positive"
    assert pc.normalise_phenotype("Negative") == "negative"


# ---------------------------------------------------------------- translation

def test_dpyd_reference_allele_is_translated(pc):
    """PharmCAT names DPYD's reference allele 'Reference'; '*1' returns Indeterminate."""
    call, route = pc.to_outside_call("DPYD", "*1/*1")
    assert route == "diplotype"
    assert call == "Reference/Reference"


def test_cyp2d6_duplication_gets_an_explicit_copy_number(pc):
    """'*1xN' does not parse in PharmCAT at all."""
    call, _ = pc.to_outside_call("CYP2D6", "*1/*1xN (gene duplication)")
    assert "xN" not in call
    assert call == "*1/*1x2"


def test_hla_positive_uses_pharmcat_status_notation(pc):
    call, route = pc.to_outside_call("HLA-B*57:01", "HLA-B*57:01 positive")
    assert route == "hla_status"
    assert call == "*57:01 positive"


def test_non_diplotype_states_are_declared_not_expressible(pc):
    for gene, dip in [("G6PD", "G6PD Mediterranean (hemizygous male or homozygous female)"),
                      ("MT-RNR1", "m.1555A>G"),
                      ("RYR1", "RYR1 pathogenic variant present"),
                      ("CFTR", "G551D mutation present")]:
        call, route = pc.to_outside_call(gene, dip)
        assert call is None and route == "not_expressible", (gene, dip)


def test_translations_are_declared_not_inferred(pc):
    assert pc.ALLELE_TRANSLATIONS["DPYD"]["*1"] == "Reference"
    assert pc.ALLELE_TRANSLATIONS["CYP2D6"]["*1xN"] == "*1x2"


# ---------------------------------------------------------------- batching

def test_batches_never_repeat_a_gene(pc):
    """One phenotyper run reports each gene once, so a repeated gene loses data."""
    pairs = [("CYP2D6", "*1/*1"), ("CYP2D6", "*4/*4"), ("CYP2C19", "*2/*2")]
    for batch in pc.batch_states(pairs):
        genes = [g for g, _ in batch]
        assert len(genes) == len(set(genes))


def test_batching_preserves_every_state(pc):
    pairs = [("CYP2D6", f"*{i}/*1") for i in range(5)] + [("TPMT", "*1/*1")]
    flat = [p for batch in pc.batch_states(pairs) for p in batch]
    assert sorted(flat) == sorted(pairs)


# ---------------------------------------------------------------- classification

def test_indeterminate_is_an_abstention_not_a_disagreement(pc):
    """PharmCAT refusing to call is not PharmCAT making a rival call."""
    case = next(c for c in pc.CASES if c["gene"] == "DPYD")
    call, _ = pc.to_outside_call(case["gene"], case["gt_diplotype"])
    rows = pc.compare([case], {(case["gene"], call): "Indeterminate"})
    assert rows[0]["pharmcat_abstained"] is True
    assert rows[0]["ours_vs_pharmcat"] is None
    summary = pc.summarise(rows)
    assert summary["pharmcat_abstained"] == 1
    assert summary["ours_vs_pharmcat_disagree"] == 0


def test_a_genuine_disagreement_is_still_counted(pc):
    """Guard against the abstention rule swallowing real discordance."""
    case = next(c for c in pc.CASES if c["id"] == "cyp2d6_codeine_pm")
    call, _ = pc.to_outside_call(case["gene"], case["gt_diplotype"])
    rows = pc.compare([case], {(case["gene"], call): "Ultrarapid Metabolizer"})
    assert rows[0]["ours_vs_pharmcat"] is False
    assert pc.summarise(rows)["ours_vs_pharmcat_disagree"] == 1
