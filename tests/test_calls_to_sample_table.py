"""Tests for 13_calls_to_sample_table.py.

The extractor turns PyPGx per-gene results archives into the (sample, gene,
diplotype) table that 08_caller_truth_eval.py consumes. The failure mode worth
testing is silence: a gene whose archive is missing or malformed must be
reported, not skipped into a smaller and apparently cleaner table.
"""
import io
import zipfile
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "real-genome-arm"
    / "scripts"
    / "13_calls_to_sample_table.py"
)


def load():
    spec = spec_from_file_location("calls_to_table", SCRIPT)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_results(gene_dir: Path, rows, header="\tGenotype\tPhenotype"):
    """Write a minimal PyPGx results.zip: metadata.txt plus data.tsv."""
    gene_dir.mkdir(parents=True, exist_ok=True)
    body = io.StringIO()
    body.write(header + "\n")
    for sample, genotype, phenotype in rows:
        body.write(f"{sample}\t{genotype}\t{phenotype}\n")
    with zipfile.ZipFile(gene_dir / "results.zip", "w") as z:
        z.writestr("tmp/metadata.txt", "Gene=X\nAssembly=GRCh37\n")
        z.writestr("tmp/data.tsv", body.getvalue())


def test_extracts_sample_gene_diplotype(tmp_path):
    mod = load()
    write_results(
        tmp_path / "CYP2C19",
        [("HG00111", "*1/*2", "Intermediate Metabolizer"),
         ("NA12878", "*1/*1", "Normal Metabolizer")],
    )
    rows, problems = mod.collect(tmp_path)
    assert problems == []
    assert sorted(rows, key=lambda r: r["sample"]) == [
        {"sample": "HG00111", "gene": "CYP2C19", "diplotype": "*1/*2"},
        {"sample": "NA12878", "gene": "CYP2C19", "diplotype": "*1/*1"},
    ]


def test_gene_name_comes_from_the_directory(tmp_path):
    mod = load()
    write_results(tmp_path / "TPMT", [("HG00111", "*1/*3A", "IM")])
    rows, _ = mod.collect(tmp_path)
    assert rows[0]["gene"] == "TPMT"


def test_missing_archive_is_reported_not_skipped(tmp_path):
    mod = load()
    write_results(tmp_path / "TPMT", [("HG00111", "*1/*1", "NM")])
    (tmp_path / "CYP2D6").mkdir()  # called, but produced no results archive
    rows, problems = mod.collect(tmp_path)
    assert len(rows) == 1
    assert any("CYP2D6" in p for p in problems)


def test_malformed_archive_is_reported_not_skipped(tmp_path):
    mod = load()
    (tmp_path / "DPYD").mkdir(parents=True)
    (tmp_path / "DPYD" / "results.zip").write_bytes(b"not a zip file")
    rows, problems = mod.collect(tmp_path)
    assert rows == []
    assert any("DPYD" in p for p in problems)


def test_indeterminate_calls_are_dropped_and_counted(tmp_path):
    """PyPGx emits an empty or Indeterminate genotype where it cannot call.

    Those are abstentions, not calls, and scoring them against truth would
    manufacture disagreements out of a refusal to answer.
    """
    mod = load()
    write_results(
        tmp_path / "CYP2C9",
        [("HG00111", "*1/*2", "IM"),
         ("NA12878", "", "Indeterminate"),
         ("NA18980", "Indeterminate", "Indeterminate")],
    )
    rows, _ = mod.collect(tmp_path)
    assert [r["sample"] for r in rows] == ["HG00111"]


def test_refuses_an_empty_call_set(tmp_path):
    """An empty table must raise, not write a file that scores vacuously."""
    mod = load()
    (tmp_path / "empty").mkdir()
    with pytest.raises(SystemExit):
        mod.main(["--calls", str(tmp_path), "--out", str(tmp_path / "o.tsv")])
