"""Tests for N0: end-to-end executed pipeline (PyPGx call -> executed skill -> abstain).

The claim under test (Cell Genomics revision, R2.2): on real diplotypes, replacing the
model's input call with a deterministic caller and abstaining on no-calls / ambiguous /
out-of-scope diplotypes yields 100% correctness among emitted in-scope answers, at a
measured coverage. These tests pin the deterministic logic (skill map, abstention
classifier, aggregation) so the live run on real cohorts is trustworthy.
"""
import importlib.util
from pathlib import Path

MOD_PATH = Path(__file__).resolve().parents[1] / "real-genome-arm" / "scripts" / "07_executed_pipeline_n0.py"
spec = importlib.util.spec_from_file_location("n0", MOD_PATH)
n0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(n0)


def test_skill_map_loads_from_cases():
    # 110 canonical CPIC cases -> a substantial (gene, diplotype) vocabulary
    assert len(n0.SKILL_MAP) > 50


def test_norm_dip_strips_gene_and_orders():
    assert n0.norm_dip("CYP2D6 *4/*4", "CYP2D6") == "*4/*4"
    assert n0.norm_dip("*17/*1", "CYP2C19") == "*1/*17"  # sorted


def test_classify_emitted_known_diplotype():
    sm = {("CYP2D6", "*4/*4"): "Poor Metaboliser"}
    assert n0.classify("CYP2D6", "*4/*4", "Poor Metabolizer", sm) == ("emitted", "")


def test_classify_out_of_scope():
    sm = {("CYP2D6", "*4/*4"): "Poor Metaboliser"}
    status, reason = n0.classify("CYP2D6", "*1/*99", "Normal Metabolizer", sm)
    assert (status, reason) == ("abstain", "out_of_scope")


def test_classify_no_call_on_blank_or_indeterminate():
    sm = {("CYP2D6", "*4/*4"): "Poor Metaboliser"}
    assert n0.classify("CYP2D6", "", "Indeterminate", sm)[1] == "no_call"
    assert n0.classify("CYP2D6", "*4/*4", "Indeterminate", sm)[1] == "no_call"


def test_classify_ambiguous_uncertain_function():
    sm = {("CYP2D6", "*4/*4"): "Poor Metaboliser"}
    # in-vocab diplotype but PyPGx flags uncertain function -> abstain ambiguous
    assert n0.classify("CYP2D6", "*4/*4", "Poor Metabolizer (uncertain function allele)", sm)[1] == "ambiguous"


def test_run_aggregation_coverage_and_accuracy():
    sm = {("CYP2D6", "*4/*4"): "Poor Metaboliser", ("CYP2C19", "*2/*2"): "Poor Metaboliser"}
    rows = [
        {"cohort": "EUR", "gene": "CYP2D6", "diplotype": "*4/*4", "phenotype": "Poor Metabolizer", "n_carriers": 3},
        {"cohort": "EUR", "gene": "CYP2C19", "diplotype": "*2/*2", "phenotype": "Poor Metabolizer", "n_carriers": 2},
        {"cohort": "EUR", "gene": "CYP2D6", "diplotype": "*1/*99", "phenotype": "Normal Metabolizer", "n_carriers": 5},
        {"cohort": "EUR", "gene": "CYP2C19", "diplotype": "", "phenotype": "Indeterminate", "n_carriers": 1},
    ]
    summ = n0.run(rows, sm, n0.score_a1)["EUR"]
    assert summ["n_states"] == 4
    assert summ["emitted"] == 2
    assert summ["abstain"]["out_of_scope"] == 1
    assert summ["abstain"]["no_call"] == 1
    assert summ["correct_emitted"] == 2
    assert summ["accuracy_among_emitted"] == 1.0
    assert round(summ["coverage_states"], 3) == 0.5


def test_run_flags_a_disagreement_as_falsifier():
    # If the executed skill disagrees with PyPGx on an in-scope diplotype, it must be
    # counted as incorrect and recorded, not hidden. This is the falsifier for the claim.
    sm = {("CYP2D6", "*4/*4"): "Poor Metaboliser"}
    rows = [{"cohort": "AFR", "gene": "CYP2D6", "diplotype": "*4/*4",
             "phenotype": "Normal Metabolizer", "n_carriers": 1}]  # PyPGx says NM, skill says PM
    summ = n0.run(rows, sm, n0.score_a1)["AFR"]
    assert summ["emitted"] == 1
    assert summ["correct_emitted"] == 0
    assert summ["accuracy_among_emitted"] == 0.0
    assert len(summ["disagreements"]) == 1
