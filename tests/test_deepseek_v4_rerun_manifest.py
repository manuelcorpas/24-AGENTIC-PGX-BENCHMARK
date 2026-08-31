import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs" / "deepseek-v4-flash-rerun-audit-2026-08-30"
MANIFEST = AUDIT / "RERUN-MANIFEST.json"


def test_manifest_has_complete_pinned_rerun_scope():
    data = json.loads(MANIFEST.read_text())

    assert data["target"]["manuscript_label"] == "DeepSeek V4 Flash"
    assert data["target"]["request_model_id"] == "deepseek-v4-flash"
    assert data["target"]["thinking_mode"] == "disabled"
    assert data["target"]["returned_model_gate"] == "deepseek-v4-flash"

    jobs = data["jobs"]
    assert len(jobs) == 10
    assert sum(job["expected_calls"] for job in jobs) == 8901
    assert len({job["id"] for job in jobs}) == len(jobs)
    assert all(job["preserve_existing"] for job in jobs)
    assert all(job["output_path"].startswith("data/v4_flash_") for job in jobs)


def test_manifest_covers_every_reported_dependency_class():
    data = json.loads(MANIFEST.read_text())
    classes = {job["dependency_class"] for job in data["jobs"]}
    assert classes == {
        "matched_five_cell",
        "extracted_rules",
        "adversarial_contract",
        "drug_keyed_retrieval",
        "population_framing",
        "input_normalisation",
        "real_genome",
    }

    required_items = {
        "Table 1", "Figure 2", "Table 2", "Figure 3", "Figure 4",
        "Figure S2", "Figure S3", "Figure S4", "Figure S5",
        "Figure S7", "Figure S8", "Figure S9",
        "Table S1", "Table S2", "Table S3", "Table S4", "Table S5",
        "Table S6", "Table S7", "Table S8", "Table S9",
    }
    assert required_items <= set(data["downstream"]["display_items"])


def test_manifest_has_non_averageable_provenance_gates():
    data = json.loads(MANIFEST.read_text())
    gates = set(data["acceptance_gates"])
    assert "all_rows_return_requested_model_id" in gates
    assert "no_deepseek_chat_request_alias" in gates
    assert "all_expected_units_exactly_once" in gates
    assert "no_overwrite_of_historical_raw_data" in gates
    assert "all_manuscript_deepseek_v3_labels_removed" in gates
