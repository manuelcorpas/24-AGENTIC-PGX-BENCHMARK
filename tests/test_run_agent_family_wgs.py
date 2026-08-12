"""Tests for the family-arm agent runner.

The guards here are the ones the corrections log says this project needs:
an unpriced or unimplemented model must stop the run before any call is made,
and a provider error must never be recorded as a model abstention.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

_p = (
    Path(__file__).resolve().parents[1]
    / "real-genome-arm"
    / "scripts"
    / "19_run_agent_family_wgs.py"
)
_spec = spec_from_file_location("runner", _p)
runner = module_from_spec(_spec)
_spec.loader.exec_module(runner)


class TestPrompt:
    def test_contains_gene_and_diplotype(self):
        p = runner.build_prompt("CYP2D6", "*1/*68+*4")
        assert "Gene: CYP2D6" in p
        assert "Diplotype: *1/*68+*4" in p

    def test_wording_matches_the_original_runner(self):
        # Predictions from script 04 are reused alongside these, so a changed
        # prompt would silently mix two different questions in one table.
        p = runner.build_prompt("TPMT", "*1/*9")
        assert p.startswith("You are executing a ClawBio pharmacogenomics skill.")
        assert p.endswith("Output one line:\nPHENOTYPE: [CPIC phenotype]")


class TestParse:
    def test_extracts_labelled_phenotype(self):
        assert runner.parse_phenotype("PHENOTYPE: Poor Metabolizer") == "Poor Metabolizer"

    def test_is_case_insensitive_on_the_label(self):
        assert runner.parse_phenotype("phenotype: Normal Metabolizer") == "Normal Metabolizer"

    def test_finds_the_line_among_others(self):
        txt = "Some preamble\nPHENOTYPE: Intermediate Metabolizer\ntrailing"
        assert runner.parse_phenotype(txt) == "Intermediate Metabolizer"

    def test_none_does_not_raise(self):
        assert runner.parse_phenotype(None) == ""

    def test_unlabelled_output_is_truncated_not_invented(self):
        assert runner.parse_phenotype("x" * 200) == "x" * 80


class TestCost:
    def test_uses_the_fitted_rate(self):
        # Claude Opus 4.5: $5/M in, $25/M out
        assert runner.cost_usd("Claude Opus 4.5", 1_000_000, 0) == pytest.approx(5.00)
        assert runner.cost_usd("Claude Opus 4.5", 0, 1_000_000) == pytest.approx(25.00)

    def test_small_call_is_priced(self):
        assert runner.cost_usd("DeepSeek V3", 95, 47) == pytest.approx(
            95 * 0.27 / 1e6 + 47 * 1.10 / 1e6
        )

    def test_unpriced_model_raises(self):
        with pytest.raises(KeyError):
            runner.cost_usd("Mistral Large 2512", 10, 10)


class TestPanelLoading:
    def _write(self, tmp_path, text):
        p = tmp_path / "models.txt"
        p.write_text(text)
        return p

    def test_parses_a_well_formed_panel(self, tmp_path):
        p = self._write(tmp_path, "# comment\nClaude Opus 4.5      claude-opus-4-5      anthropic\n")
        assert runner.load_panel(p) == [("Claude Opus 4.5", "claude-opus-4-5", "anthropic")]

    def test_rejects_an_unknown_provider(self, tmp_path):
        p = self._write(tmp_path, "Claude Opus 4.5      claude-opus-4-5      mistralai\n")
        with pytest.raises(SystemExit, match="unknown provider"):
            runner.load_panel(p)

    def test_rejects_a_model_with_no_price(self, tmp_path):
        # Running a model whose cost cannot be computed produces a spend with no
        # record of what it cost, which is the provenance gap this script exists
        # to close.
        p = self._write(tmp_path, "Mistral Large 2512      mistral-large      openai\n")
        with pytest.raises(SystemExit, match="no price recorded"):
            runner.load_panel(p)

    def test_empty_panel_raises(self, tmp_path):
        p = self._write(tmp_path, "# only a comment\n")
        with pytest.raises(SystemExit, match="no models parsed"):
            runner.load_panel(p)


class TestReadStates:
    def test_deduplicates_states(self, tmp_path):
        p = tmp_path / "s.tsv"
        p.write_text(
            "cohort\tgene\tdiplotype\tphenotype\tn_carriers\n"
            "F\tTPMT\t*1/*1\tNormal Metabolizer\t2\n"
            "F\tTPMT\t*1/*1\tNormal Metabolizer\t2\n"
            "F\tTPMT\t*1/*9\tIndeterminate\t2\n"
        )
        states = runner.read_states(p)
        assert len(states) == 2
        assert states[0]["caller_phenotype"] == "Normal Metabolizer"
