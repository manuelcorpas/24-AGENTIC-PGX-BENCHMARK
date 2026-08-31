"""Tests for provenance-preserving real-genome scoring.

These tests define the correction contract before the replacement runner and
scorer are implemented. A response can be clinically wrong only after it has
produced one valid phenotype under the supplied one-line output contract.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


contract = load(
    REPO / "real-genome-arm" / "scripts" / "20_realgenome_contract.py",
    "realgenome_contract",
)
runner = load(
    REPO / "real-genome-arm" / "scripts" / "21_run_agent_four_cohorts.py",
    "realgenome_runner",
)
scorer = load(
    REPO / "real-genome-arm" / "scripts" / "22_score_agent_realgenome.py",
    "realgenome_scorer",
)


class TestContractParser:
    def test_accepts_one_exact_labelled_line(self):
        parsed = contract.parse_contract_output(
            "PHENOTYPE: Poor Metabolizer", stop_reason="end_turn"
        )
        assert parsed["parse_status"] == "valid"
        assert parsed["phenotype_key"] == "poor metabolizer"
        assert parsed["invalid_reason"] is None

    def test_accepts_british_spelling_without_substring_matching(self):
        parsed = contract.parse_contract_output(
            "phenotype: Intermediate Metaboliser", stop_reason="stop"
        )
        assert parsed["parse_status"] == "valid"
        assert parsed["phenotype_key"] == "intermediate metabolizer"

    def test_missing_marker_is_invalid_not_a_prediction(self):
        parsed = contract.parse_contract_output(
            "I need to inspect the allele functions first", stop_reason="max_tokens"
        )
        assert parsed["parse_status"] == "invalid"
        assert parsed["phenotype_key"] is None
        assert parsed["invalid_reason"] == "truncated_without_phenotype"

    def test_extra_reasoning_is_scorable_but_flagged_as_nonconformant(self):
        parsed = contract.parse_contract_output(
            "Reasoning first\nPHENOTYPE: Normal Metabolizer", stop_reason="end_turn"
        )
        assert parsed["parse_status"] == "valid"
        assert parsed["invalid_reason"] is None
        assert parsed["contract_status"] == "extra_text"
        assert parsed["extracted_phenotype"] == "Normal Metabolizer"

    def test_unrecognised_value_is_parseable_but_not_normalised(self):
        parsed = contract.parse_contract_output(
            "PHENOTYPE: probably fine", stop_reason="end_turn"
        )
        assert parsed["parse_status"] == "valid"
        assert parsed["invalid_reason"] is None
        assert parsed["value_status"] == "unrecognised"
        assert parsed["phenotype_key"] is None

    def test_explicit_indeterminate_is_a_valid_abstention(self):
        parsed = contract.parse_contract_output(
            "PHENOTYPE: Indeterminate", stop_reason="end_turn"
        )
        assert parsed["parse_status"] == "valid"
        assert parsed["phenotype_key"] == "indeterminate"

    def test_uncertain_susceptibility_is_not_an_abstention(self):
        assert contract.normalise_phenotype("Uncertain Susceptibility") == (
            "uncertain susceptibility"
        )
        assert contract.reference_is_scorable("Uncertain Susceptibility") is True
        assert contract.reference_is_scorable("Indeterminate") is False

    def test_gene_prefix_and_standard_abbreviation_do_not_change_the_label(self):
        assert contract.normalise_phenotype(
            "CYP2D6 Intermediate Metabolizer", gene="CYP2D6"
        ) == "intermediate metabolizer"
        assert contract.normalise_phenotype(
            "Normal metabolizer (NM)", gene="CYP2D6"
        ) == "normal metabolizer"

    def test_cyp3a5_expressor_qualifier_does_not_change_the_metabolizer_tier(self):
        assert contract.normalise_phenotype(
            "Normal Metaboliser (Expressor)", gene="CYP3A5"
        ) == "normal metabolizer"
        assert contract.normalise_phenotype(
            "Poor Metaboliser (Non-expressor)", gene="CYP3A5"
        ) == "poor metabolizer"

    def test_clear_nonanswer_is_an_abstention(self):
        assert contract.normalise_phenotype("Not applicable") == "indeterminate"


class TestProviderExtraction:
    def test_reasoning_and_nonreasoning_openai_share_one_rate_limit(self):
        assert runner.provider_group("openai") == "openai"
        assert runner.provider_group("openai-reasoning") == "openai"

    def test_google_receives_a_separate_thinking_aware_budget(self):
        payload = runner.request_payload(
            "gemini-2.5-flash",
            "google",
            "prompt",
            standard_tokens=1024,
            reasoning_tokens=4096,
            google_tokens=4096,
        )
        assert payload["generationConfig"]["maxOutputTokens"] == 4096

    def test_anthropic_raw_text_and_stop_reason_are_retained(self):
        body = {
            "id": "msg_123",
            "model": "claude-opus-4-5-20251101",
            "content": [{"type": "text", "text": "PHENOTYPE: Normal Metabolizer"}],
            "stop_reason": "end_turn",
            "usage": {"input_tokens": 91, "output_tokens": 8},
        }
        result = runner.extract_provider_result("anthropic", body)
        assert result["raw_text"] == "PHENOTYPE: Normal Metabolizer"
        assert result["provider_stop_reason"] == "end_turn"
        assert result["response_id"] == "msg_123"
        assert result["input_tokens"] == 91
        assert result["output_tokens"] == 8

    def test_openai_finish_reason_is_retained(self):
        body = {
            "id": "chatcmpl_123",
            "choices": [
                {
                    "finish_reason": "length",
                    "message": {"content": "I need to determine the phenotype"},
                }
            ],
            "usage": {"prompt_tokens": 90, "completion_tokens": 1024},
        }
        result = runner.extract_provider_result("openai", body)
        assert result["provider_stop_reason"] == "length"
        assert result["raw_text"] == "I need to determine the phenotype"

    def test_google_finish_reason_is_retained(self):
        body = {
            "responseId": "google_123",
            "candidates": [
                {
                    "finishReason": "STOP",
                    "content": {"parts": [{"text": "PHENOTYPE: Poor Function"}]},
                }
            ],
            "usageMetadata": {"promptTokenCount": 90, "candidatesTokenCount": 7},
        }
        result = runner.extract_provider_result("google", body)
        assert result["provider_stop_reason"] == "STOP"
        assert result["raw_text"] == "PHENOTYPE: Poor Function"


class TestFourWayScoring:
    def record(self, raw_text, reference="Normal Metabolizer", stop="end_turn", error=None):
        return {
            "cohort": "F",
            "gene": "CYP2D6",
            "diplotype": "*1/*1",
            "caller_phenotype": reference,
            "raw_text": raw_text,
            "provider_stop_reason": stop,
            "error": error,
        }

    def test_correct(self):
        assert scorer.score_record(self.record("PHENOTYPE: Normal Metabolizer"))["outcome"] == "correct"

    def test_clinically_wrong_requires_a_valid_output(self):
        assert scorer.score_record(self.record("PHENOTYPE: Poor Metabolizer"))["outcome"] == "clinically_wrong"

    def test_unlabelled_text_is_invalid_not_clinically_wrong(self):
        scored = scorer.score_record(
            self.record("I need to determine the phenotype", stop="max_tokens")
        )
        assert scored["outcome"] == "invalid_output"
        assert scored["invalid_reason"] == "truncated_without_phenotype"

    def test_labelled_but_unrecognised_value_is_clinically_wrong(self):
        scored = scorer.score_record(self.record("PHENOTYPE: probably fine"))
        assert scored["outcome"] == "clinically_wrong"

    def test_provider_error_is_invalid_not_abstention(self):
        scored = scorer.score_record(
            self.record(None, stop=None, error={"type": "HTTPError", "status": 429})
        )
        assert scored["outcome"] == "invalid_output"
        assert scored["invalid_reason"] == "provider_error"

    def test_explicit_indeterminate_is_abstention(self):
        assert scorer.score_record(self.record("PHENOTYPE: Indeterminate"))["outcome"] == "abstained"

    def test_reference_indeterminate_is_excluded(self):
        scored = scorer.score_record(
            self.record("PHENOTYPE: Normal Metabolizer", reference="Indeterminate")
        )
        assert scored["outcome"] == "excluded_reference"

    def test_summary_has_four_mutually_exclusive_categories(self):
        records = [
            scorer.score_record(self.record("PHENOTYPE: Normal Metabolizer")),
            scorer.score_record(self.record("PHENOTYPE: Poor Metabolizer")),
            scorer.score_record(self.record("PHENOTYPE: Indeterminate")),
            scorer.score_record(self.record("unlabelled")),
        ]
        summary = scorer.summarise(records)
        assert summary["n"] == 4
        assert summary["correct"] == 1
        assert summary["abstained"] == 1
        assert summary["clinically_wrong"] == 1
        assert summary["invalid_output"] == 1
        assert sum(summary[key] for key in scorer.OUTCOMES) == summary["n"]

    def test_deterministic_control_retains_valid_uncertain_susceptibility(self):
        states = [
            {
                "cohort": "Family",
                "gene": "RYR1",
                "diplotype": "Negative/Negative",
                "phenotype": "Uncertain Susceptibility",
            },
            {
                "cohort": "Family",
                "gene": "CYP4F2",
                "diplotype": "*1/*1",
                "phenotype": "Indeterminate",
            },
        ]
        mapping = {("RYR1", "neg/neg"): "Uncertain Susceptibility"}
        result = scorer.deterministic_arm(states, mapping)
        assert result["Family"] == {
            "n": 1,
            "emitted": 1,
            "correct": 1,
            "abstained": 0,
            "coverage": 1.0,
            "accuracy_among_emitted": 1.0,
            "accuracy_overall": 1.0,
        }

    def test_deterministic_coverage_is_not_rounded_before_reporting(self):
        summary = scorer.summarise_deterministic(
            {"n": 95, "emitted": 20, "correct": 20}
        )
        assert summary["coverage"] == 20 / 95
        assert f"{summary['coverage']:.3f}" == "0.211"


class TestModelVersionAudit:
    def test_dated_openai_snapshot_matches_its_requested_alias(self):
        assert scorer.returned_model_matches("gpt-4.1", "gpt-4.1-2025-04-14") is True

    def test_deepseek_family_change_is_not_silently_accepted(self):
        assert scorer.returned_model_matches("deepseek-chat", "deepseek-v4-flash") is False
