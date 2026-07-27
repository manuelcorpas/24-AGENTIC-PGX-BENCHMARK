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
import csv
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
    """Write a fixture with real CSV quoting.

    Naive comma-joining silently splits a header like
    'VKORC1 NM_024006.6:c.196G>A, rs72547529', which the published sheet
    genuinely contains, and shifts every later column by one.
    """
    p = tmp_path / "consolidated.csv"
    with p.open("w", newline="") as fh:
        csv.writer(fh).writerows(rows)
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


# ---------------------------------------------------------------- real sheet shape

def test_header_row_is_found_below_a_blank_preamble(ing, tmp_path):
    """The published sheet starts with a blank row and puts the header on the
    second. Treating row 0 as the header yields a table whose columns are all
    empty strings, which silently produces zero genotypes."""
    src = write_csv(tmp_path, [
        ["", "", ""],
        ["GeT-RM Characterization", "Coriell ID #", "CYP2C19"],
        ["3", "HG00111", "*1/*2"],
    ])
    rows = ing.to_long(ing.read_table(src))
    assert [(r["sample"], r["gene"], r["diplotype"]) for r in rows] == \
        [("HG00111", "CYP2C19", "*1/*2")]


def test_per_gene_reference_columns_are_not_treated_as_genes(ing, tmp_path):
    """The sheet interleaves '<GENE> References' after each gene column."""
    src = write_csv(tmp_path, [
        ["Coriell ID #", "CYP1A2", "CYP1A2 References", "CYP2B6", "CYP2B6 References"],
        ["HG00111", "*1/*1", "Pratt 2016", "*1/*6", "Gaedigk 2019"],
    ])
    rows = ing.to_long(ing.read_table(src))
    assert sorted(r["gene"] for r in rows) == ["CYP1A2", "CYP2B6"]


def test_accession_and_ftp_columns_do_not_become_genotypes(ing, tmp_path):
    src = write_csv(tmp_path, [
        ["Coriell ID #", "run_accession", "sra_ftp", "TPMT"],
        ["HG00111", "ERR000123", "ftp.sra.ebi.ac.uk/vol1/x.fastq.gz", "*1/*1"],
    ])
    rows = ing.to_long(ing.read_table(src))
    assert [r["gene"] for r in rows] == ["TPMT"]


def test_header_names_are_stripped_of_whitespace_and_newlines(ing, tmp_path):
    src = tmp_path / "h.csv"
    src.write_text('"Coriell ID #  \nhttps://example","CYP2C19  "\n"HG00111","*1/*2"\n')
    rows = ing.to_long(ing.read_table(src))
    assert rows[0]["sample"] == "HG00111"
    assert rows[0]["gene"] == "CYP2C19"


def test_nucleotide_genotypes_are_not_star_allele_diplotypes(ing, tmp_path):
    """The consolidated sheet carries rsID-specific SNP columns alongside the
    star-allele columns, with values like 'G/A'. Those match the shape of a
    diplotype but are not one. Admitted into the truth set they would be scored
    against PyPGx star-allele calls, mismatch every time, and understate
    concordance: a confident wrong result."""
    src = write_csv(tmp_path, [
        ["Coriell ID #", "VKORC1 NM_024006.6:c.196G>A, rs72547529", "CYP2C19"],
        ["HG00111", "G/A", "*1/*2"],
    ])
    rows = ing.to_long(ing.read_table(src))
    assert [r["gene"] for r in rows] == ["CYP2C19"]


@pytest.mark.parametrize("cell", ["G/G", "A/A", "C/G", "T/T", "A", "G>A"])
def test_single_base_calls_are_rejected(ing, cell, tmp_path):
    src = write_csv(tmp_path, [["Coriell ID #", "GGCX"], ["HG00111", cell]])
    assert ing.to_long(ing.read_table(src)) == []


def test_wild_type_shorthand_is_rejected_as_out_of_vocabulary(ing, tmp_path):
    """SLCO2B1 is reported as 'WT/WT'. That is a real call, but it is not in the
    star-allele vocabulary PyPGx emits, so comparing it would produce spurious
    mismatches rather than information. Excluded deliberately, and reported."""
    src = write_csv(tmp_path, [["Coriell ID #", "SLCO2B1"], ["HG00111", "WT/WT"]])
    assert ing.to_long(ing.read_table(src)) == []


def test_dropped_columns_are_reported_not_silently_discarded(ing, tmp_path):
    src = write_csv(tmp_path, [
        ["Coriell ID #", "CYP2C19", "SLCO2B1"],
        ["HG00111", "*1/*2", "WT/WT"],
    ])
    rows, dropped = ing.to_long(ing.read_table(src), report=True)
    assert [r["gene"] for r in rows] == ["CYP2C19"]
    assert "SLCO2B1" in dropped


def test_output_is_the_schema_the_evaluator_requires(ing, tmp_path):
    src = write_csv(tmp_path, [["Coriell #", "TPMT"], ["NA12878", "*1/*1"]])
    out = tmp_path / "getrm_consensus.tsv"
    ing.build(src, out)
    header = out.read_text().splitlines()[0].split("\t")
    assert header == ["sample", "gene", "diplotype"]

def test_hla_alleles_with_colons_are_kept(ing, tmp_path):
    """GeT-RM reports HLA as colon-separated star alleles (*57:01:01). Four of
    the benchmark's lethal-class loci are HLA, so silently dropping them would
    remove exactly the cases where a wrong call is most dangerous."""
    src = write_csv(tmp_path, [
        ["Coriell ID #", "HLA-A", "HLA-B"],
        ["HG00111", "*30:02:01", "*57:01:01"],
    ])
    rows = ing.to_long(ing.read_table(src))
    got = {r["gene"]: r["diplotype"] for r in rows}
    assert got == {"HLA-A": "*30:02:01", "HLA-B": "*57:01:01"}


def test_hla_diplotype_pairs_are_kept(ing, tmp_path):
    src = write_csv(tmp_path, [
        ["Coriell ID #", "HLA-B"], ["HG00111", "*57:01:01/*44:02:01"]])
    rows = ing.to_long(ing.read_table(src))
    assert rows[0]["diplotype"] == "*57:01:01/*44:02:01"
