"""Tests for the Corpas family WGS aggregation and Mendelian check.

The Mendelian check is the load-bearing part: it is the only independent
evidence that the WGS calls are internally coherent, since there is no external
truth set for this family. A wrong implementation would either hide real
inconsistencies or invent them, and both would reach the manuscript as a
statement about caller reliability.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

_p = Path(__file__).resolve().parents[1] / "real-genome-arm" / "scripts" / "17_aggregate_family_wgs.py"
_spec = spec_from_file_location("agg", _p)
agg = module_from_spec(_spec)
_spec.loader.exec_module(agg)


class TestSplitDiplotype:
    def test_plain_star_alleles(self):
        assert agg.split_diplotype("*1/*2") == ("*1", "*2")

    def test_reference_homozygote(self):
        assert agg.split_diplotype("Reference/Reference") == ("Reference", "Reference")

    def test_allele_containing_a_plus(self):
        # UGT1A1 *80+*28 is ONE allele, not two.
        assert agg.split_diplotype("*1/*80+*28") == ("*1", "*80+*28")

    def test_allele_containing_hgvs_and_parentheses(self):
        assert agg.split_diplotype("Reference/c.85T>C (*9A)") == (
            "Reference",
            "c.85T>C (*9A)",
        )

    def test_rsid_named_alleles(self):
        assert agg.split_diplotype("rs12979860/rs12980275") == (
            "rs12979860",
            "rs12980275",
        )

    def test_empty_is_none(self):
        assert agg.split_diplotype("") is None
        assert agg.split_diplotype("   ") is None

    def test_missing_separator_is_none(self):
        assert agg.split_diplotype("*1") is None


class TestMendelian:
    def test_simple_consistent_transmission(self):
        # father *1/*2, mother *1/*3, child *2/*3
        assert agg.mendelian_ok("*2/*3", "*1/*2", "*1/*3") is True

    def test_consistent_in_either_orientation(self):
        # child *3/*2 is the same transmission written the other way round
        assert agg.mendelian_ok("*3/*2", "*1/*2", "*1/*3") is True

    def test_homozygous_parent_forces_an_allele(self):
        # father *6/*6 must transmit *6; a child without *6 is inconsistent
        assert agg.mendelian_ok("*2/*9", "*6/*6", "*2/*5") is False

    def test_child_allele_absent_from_mother(self):
        assert agg.mendelian_ok("*1/*4", "*1/*2", "*10/*41") is False

    def test_homozygous_child_needs_the_allele_from_both(self):
        assert agg.mendelian_ok("*2/*2", "*1/*2", "*1/*2") is True
        assert agg.mendelian_ok("*2/*2", "*1/*2", "*1/*1") is False

    def test_reference_alleles_behave_like_any_other(self):
        assert agg.mendelian_ok(
            "Reference/Reference", "Reference/Reference", "Reference/*2"
        ) is True
        assert agg.mendelian_ok(
            "Reference/Reference", "rs1/rs2", "Reference/Reference"
        ) is False

    def test_unparseable_input_is_not_a_verdict(self):
        # An uncallable gene must not be silently scored as consistent.
        assert agg.mendelian_ok("", "*1/*2", "*1/*3") is None
        assert agg.mendelian_ok("*1/*2", "", "*1/*3") is None


class TestTidyRows:
    def test_counts_carriers_per_distinct_state(self):
        calls = {
            "PT1": {"CYP2C9": ("*1/*2", "Intermediate Metabolizer")},
            "PT2": {"CYP2C9": ("*1/*2", "Intermediate Metabolizer")},
            "PT3": {"CYP2C9": ("*2/*2", "Poor Metabolizer")},
        }
        rows = agg.tidy_rows(calls, "CorpasFamily")
        by_dip = {r[2]: r for r in rows}
        assert by_dip["*1/*2"][4] == 2
        assert by_dip["*2/*2"][4] == 1
        assert by_dip["*1/*2"][3] == "Intermediate Metabolizer"

    def test_blank_genotypes_are_dropped_not_counted(self):
        calls = {"PT1": {"CYP2C9": ("", "")}}
        assert agg.tidy_rows(calls, "CorpasFamily") == []

    def test_cohort_label_is_carried(self):
        calls = {"PT1": {"TPMT": ("*1/*1", "Normal Metabolizer")}}
        assert agg.tidy_rows(calls, "CorpasFamilyWGS")[0][0] == "CorpasFamilyWGS"


class TestConcordance:
    def test_identical_call_sets_are_fully_concordant(self):
        a = {"PT1": {"TPMT": ("*1/*1", "Normal Metabolizer")}}
        n, agree, diffs = agg.concordance(a, a)
        assert (n, agree, diffs) == (1, 1, [])

    def test_disagreement_is_reported_with_both_values(self):
        a = {"PT1": {"TPMT": ("*1/*1", "Normal Metabolizer")}}
        b = {"PT1": {"TPMT": ("*1/*9", "Indeterminate")}}
        n, agree, diffs = agg.concordance(a, b)
        assert (n, agree) == (1, 0)
        assert diffs == [("PT1", "TPMT", "*1/*1", "*1/*9")]

    def test_only_shared_sample_gene_pairs_are_compared(self):
        a = {"PT1": {"TPMT": ("*1/*1", "x")}, "PT2": {"TPMT": ("*1/*1", "x")}}
        b = {"PT1": {"TPMT": ("*1/*1", "x")}}
        n, agree, diffs = agg.concordance(a, b)
        assert n == 1
