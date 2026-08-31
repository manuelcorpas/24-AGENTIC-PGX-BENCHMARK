from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLIENT_PATH = ROOT / "code" / "94-deepseek-v4-client.py"
PILOT_PATH = ROOT / "code" / "95-run-deepseek-v4-pilots.py"


def load_client_module():
    spec = importlib.util.spec_from_file_location("deepseek_v4_client", CLIENT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FakeUsage:
    prompt_tokens = 12
    completion_tokens = 7
    total_tokens = 19


class FakeMessage:
    content = "PHENOTYPE: Normal Metaboliser"
    reasoning_content = None


class FakeChoice:
    message = FakeMessage()
    finish_reason = "stop"


class FakeResponse:
    id = "resp-test"
    model = "deepseek-v4-flash"
    created = 123456
    choices = [FakeChoice()]
    usage = FakeUsage()

    def model_dump(self):
        return {
            "id": self.id,
            "model": self.model,
            "created": self.created,
            "choices": [{"finish_reason": "stop"}],
            "usage": {"prompt_tokens": 12, "completion_tokens": 7, "total_tokens": 19},
        }


class FakeCompletions:
    def __init__(self):
        self.kwargs = None

    def create(self, **kwargs):
        self.kwargs = kwargs
        return FakeResponse()


class FakeClient:
    def __init__(self):
        self.chat = type("Chat", (), {})()
        self.chat.completions = FakeCompletions()


def test_pinned_non_thinking_request_and_complete_provenance():
    mod = load_client_module()
    fake = FakeClient()
    result = mod.call_deepseek_v4(fake, "test prompt", max_tokens=50)

    assert fake.chat.completions.kwargs["model"] == "deepseek-v4-flash"
    assert fake.chat.completions.kwargs["extra_body"] == {
        "thinking": {"type": "disabled"}
    }
    assert result["requested_model"] == "deepseek-v4-flash"
    assert result["returned_model"] == "deepseek-v4-flash"
    assert result["thinking_mode"] == "disabled"
    assert result["response_id"] == "resp-test"
    assert result["stop_reason"] == "stop"
    assert result["usage"]["input_tokens"] == 12
    assert result["usage"]["output_tokens"] == 7
    assert result["provider_response"]["model"] == "deepseek-v4-flash"
    assert result["request_payload"]["messages"][0]["content"] == "test prompt"
    assert result["attempt_log"][0]["status"] == "success"
    assert result["requested_at_utc"]
    assert result["completed_at_utc"]
    assert result["latency_ms"] >= 0


def test_returned_model_mismatch_is_a_hard_failure():
    mod = load_client_module()
    fake = FakeClient()
    fake.chat.completions.create = lambda **kwargs: type(
        "WrongResponse",
        (),
        {
            "id": "wrong",
            "model": "deepseek-chat",
            "created": 123,
            "choices": [FakeChoice()],
            "usage": FakeUsage(),
            "model_dump": lambda self: {"id": "wrong", "model": "deepseek-chat"},
        },
    )()

    try:
        mod.call_deepseek_v4(fake, "test prompt", max_tokens=50)
    except mod.ReturnedModelMismatch as exc:
        assert "deepseek-chat" in str(exc)
        assert exc.provenance["returned_model"] == "deepseek-chat"
        assert exc.provenance["provider_response"]["model"] == "deepseek-chat"
    else:
        raise AssertionError("returned-model mismatch did not fail")


def test_all_manifest_runners_are_pinned_and_do_not_request_alias():
    manifest_runners = [
        "code/60-matched-factorial.py",
        "code/67-extract-rules-from-prose.py",
        "code/14-adversarial-scrambled-spec.py",
        "code/14b-adversarial-reverse-spec.py",
        "code/02-run-benchmark-v3.py",
        "code/16-run-rag-condition.py",
        "code/43-armAB-fullgrid.py",
        "code/74-input-normalisation.py",
        "real-genome-arm/scripts/21_run_agent_four_cohorts.py",
    ]
    for relative in manifest_runners:
        source = (ROOT / relative).read_text()
        assert 'model="deepseek-chat"' not in source, relative
        assert '"DeepSeek V3": call_dsk' not in source, relative
        assert '"DeepSeek V3":       _dsk' not in source, relative


def test_v4_output_paths_are_separate_from_historical_files():
    expected = {
        "code/60-matched-factorial.py": "v4_flash_five_cell",
        "code/67-extract-rules-from-prose.py": "v4_flash_extracted_rules",
        "code/14-adversarial-scrambled-spec.py": "v4_flash_adversarial_forward",
        "code/14b-adversarial-reverse-spec.py": "v4_flash_adversarial_reverse",
        "code/02-run-benchmark-v3.py": "v4_flash_population_three_arm",
        "code/16-run-rag-condition.py": "v4_flash_population_three_arm",
        "code/43-armAB-fullgrid.py": "v4_flash_population_skill_arms",
        "code/74-input-normalisation.py": "v4_flash_input_normalisation",
        "real-genome-arm/scripts/21_run_agent_four_cohorts.py": "deepseek-v4-flash",
    }
    for relative, marker in expected.items():
        assert marker in (ROOT / relative).read_text(), relative


def test_pilot_manifest_has_one_request_per_distinct_prompt_family():
    spec = importlib.util.spec_from_file_location("deepseek_v4_pilots", PILOT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    families = module.build_pilot_prompts()

    assert len(families) == 18
    assert len({row["family_id"] for row in families}) == 18
    assert all(row["prompt"].strip() for row in families)
    assert {row["dependency_class"] for row in families} == {
        "matched_five_cell",
        "extracted_rules",
        "adversarial_contract",
        "drug_keyed_retrieval",
        "population_framing",
        "input_normalisation",
        "real_genome",
    }
