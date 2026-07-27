#!/usr/bin/env python3
"""
Tests for the GeT-RM consolidated-table ingester (revision item N5).

The consolidated table is a WIDE sheet: one row per Coriell sample, one column
per gene, diplotypes as star alleles. The caller-truth evaluator wants a LONG
(sample, gene, diplotype) table. This ingester converts one to the other.

What these tests pin is mostly refusal. A truth set that silently loses rows,
or that admits "Not tested" as though it were a diplotype, produces a concordance
number that looks fine and means nothing, which is the failure mode the whole
revision is guarding against.

The fixtures here are synthetic and are NOT a truth set. They exercise the
parser's shape handling only; the real table is a documented manual download
(real-genome-arm/getrm/README.md).
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "real-genome-arm" / "scripts" / "11_ingest_getrm_consolidated.py"


@pytest.fixture(scope="module")
def ing():
    spec = spec_from_file_location("getrm_ingest", SCRIPT)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_csv(tmp_path, rows):
    p = tmp_path / "consolidated.csv"
    p.write_text("\n".join(",".join(r) for r in rows) + "\n")
    return p


# ---------------------------------------------------------------- shape

def test_wide_table_becomes_long_rows(ing, tmp_path):
    src = write_csv(tmp_path, [
        ["Coriell #", "CYP2C19", "CYP2D6", "TPMT"],
        ["NA12878", "*1/*2", "*3/*4", "*1/*1"],
        ["NA18563", "*1/*1", "*1/*2", "*1/*3"],
    ])
    rows = ing.to_long(ing.read_table(src))
    assert len(rows) == 6
    assert {"sample", "gene", "diplotype"} <= set(rows[0])
    got = {(r["sample"], r["gene"]): r["diplotype"] for r in rows}
    assert got[("NA12878", "CYP2D6")] == "*3/*4"
    assert got[("NA18563", "TPMT")] == "*1/*3"


def test_sample_column_is_detected_under_several_headings(ing, tmp_path):
    for heading in ("Coriell #", "Coriell ID", "Sample", "Coriell Sample ID"):
        src = write_csv(tmp_path, [[heading, "TPMT"], ["NA12878", "*1/*1"]])
        rows = ing.to_long(ing.read_table(src))
        assert rows and rows[0]["sample"] == "NA12878", heading


def test_unparseable_sample_column_is_an_error_not_a_guess(ing, tmp_path):
    src = write_csv(tmp_path, [["mystery", "TPMT"], ["NA12878", "*1/*1"]])
    with pytest.raises(SystemExit):
        ing.to_long(ing.read_table(src))


# ---------------------------------------------------------------- refusal

@pytest.mark.parametrize("blank", ["", " ", "Not tested", "not tested", "ND", "N/A", "-"])
def test_non_genotype_placeholders_are_dropped_not_admitted(ing, tmp_path, blank):
    src = write_csv(tmp_path, [
        ["Coriell #", "TPMT"],
        ["NA12878", blank],
    ])
    assert ing.to_long(ing.read_table(src)) == []


def test_a_table_with_no_usable_genotypes_raises(ing, tmp_path):
    src = write_csv(tmp_path, [["Coriell #", "TPMT"], ["NA12878", "Not tested"]])
    with pytest.raises(SystemExit):
        ing.build(src)


def test_non_gene_metadata_columns_are_ignored(ing, tmp_path):
    src = write_csv(tmp_path, [
        ["Coriell #", "TPMT", "Reference", "BAM available", "1000 Genomes"],
        ["NA12878", "*1/*1", "Pratt 2016", "yes", "yes"],
    ])
    rows = ing.to_long(ing.read_table(src))
    assert [r["gene"] for r in rows] == ["TPMT"]


def test_diplotype_must_look_like_a_diplotype(ing, tmp_path):
    """A free-text cell is not a call. Admitting it would put prose into the
    truth set and score it as a mismatch, understating concordance."""
    src = write_csv(tmp_path, [
        ["Coriell #", "TPMT", "CYP2D6"],
        ["NA12878", "*1/*1", "see reference"],
    ])
    rows = ing.to_long(ing.read_table(src))
    assert [r["gene"] for r in rows] == ["TPMT"]


def test_gene_restriction_keeps_only_requested_genes(ing, tmp_path):
    src = write_csv(tmp_path, [
        ["Coriell #", "TPMT", "CYP2D6", "CYP2C19"],
        ["NA12878", "*1/*1", "*1/*2", "*1/*3"],
    ])
    rows = ing.to_long(ing.read_table(src), genes={"TPMT", "CYP2C19"})
    assert sorted(r["gene"] for r in rows) == ["CYP2C19", "TPMT"]


def test_output_is_the_schema_the_evaluator_requires(ing, tmp_path):
    src = write_csv(tmp_path, [["Coriell #", "TPMT"], ["NA12878", "*1/*1"]])
    out = tmp_path / "getrm_consensus.tsv"
    ing.build(src, out)
    header = out.read_text().splitlines()[0].split("\t")
    assert header == ["sample", "gene", "diplotype"]
