"""Tests for PyPGx archive sample renaming.

A partially applied or misparsed rename attaches one person's read depth to
another person's variants. That does not raise anywhere downstream; it produces
a plausible diplotype for the wrong individual. Every check here exists to make
that failure loud.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

_p = (
    Path(__file__).resolve().parents[1]
    / "real-genome-arm"
    / "scripts"
    / "18_rename_archive_samples.py"
)
_spec = spec_from_file_location("ren", _p)
ren = module_from_spec(_spec)
_spec.loader.exec_module(ren)


class TestParseMapping:
    def test_single_pair(self):
        assert ren.parse_mapping(["A=B"]) == {"A": "B"}

    def test_several_pairs(self):
        assert ren.parse_mapping(["A=B", "C=D"]) == {"A": "B", "C": "D"}

    def test_rejects_missing_separator(self):
        with pytest.raises(ValueError):
            ren.parse_mapping(["AB"])

    def test_rejects_empty_side(self):
        with pytest.raises(ValueError):
            ren.parse_mapping(["=B"])
        with pytest.raises(ValueError):
            ren.parse_mapping(["A="])

    def test_rejects_duplicate_source(self):
        with pytest.raises(ValueError):
            ren.parse_mapping(["A=B", "A=C"])

    def test_rejects_duplicate_target(self):
        # two people renamed to the same identifier collapses two genomes into one
        with pytest.raises(ValueError):
            ren.parse_mapping(["A=X", "B=X"])

    def test_rejects_whole_mapping_passed_as_one_argument(self):
        # zsh does not word-split an unquoted "$VAR"; the whole mapping arrives
        # as a single argv entry and the naive split produces a target name
        # carrying the remaining pairs.
        with pytest.raises(ValueError, match="whitespace"):
            ren.parse_mapping(["A=B C=D"])

    def test_rejects_embedded_separator(self):
        with pytest.raises(ValueError):
            ren.parse_mapping(["A=B=C"])


class TestRenameHeader:
    FIXED = ("Chromosome", "Position")

    def test_renames_sample_columns_only(self):
        header = "Chromosome\tPosition\tRUN1\tRUN2\n"
        out = ren.rename_header(header, {"RUN1": "PT1", "RUN2": "PT2"}, self.FIXED)
        assert out == "Chromosome\tPosition\tPT1\tPT2\n"

    def test_preserves_column_order(self):
        header = "Chromosome\tPosition\tRUNB\tRUNA\n"
        out = ren.rename_header(header, {"RUNA": "PT_A", "RUNB": "PT_B"}, self.FIXED)
        assert out.rstrip().split("\t") == ["Chromosome", "Position", "PT_B", "PT_A"]

    def test_raises_when_a_sample_column_has_no_mapping(self):
        header = "Chromosome\tPosition\tRUN1\tRUN2\n"
        with pytest.raises(ValueError, match="no mapping"):
            ren.rename_header(header, {"RUN1": "PT1"}, self.FIXED)

    def test_raises_when_a_mapped_name_is_absent(self):
        header = "Chromosome\tPosition\tRUN1\n"
        with pytest.raises(ValueError, match="absent"):
            ren.rename_header(header, {"RUN1": "PT1", "GHOST": "PT9"}, self.FIXED)
