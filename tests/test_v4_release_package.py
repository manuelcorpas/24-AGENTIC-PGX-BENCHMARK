import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_stage_module():
    path = ROOT / "code" / "76-stage-zenodo-revision.py"
    spec = importlib.util.spec_from_file_location("stage_zenodo_revision", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_release_identity_is_new_and_v4_bound():
    stage = load_stage_module()
    assert stage.RELEASE_VERSION == "1.8.0"
    assert stage.RELEASE_TAG == "agentic-pgx-benchmark-v2.9"
    assert stage.MODEL_LABEL == "DeepSeek V4 Flash"
    assert stage.REQUEST_MODEL == "deepseek-v4-flash"
    assert stage.RETURNED_MODEL == "deepseek-v4-flash"
    assert stage.THINKING_MODE == "disabled"


def test_ten_canonical_v4_sources_are_explicit_and_hash_bound():
    stage = load_stage_module()
    validation = json.loads(
        (ROOT / "docs" / "deepseek-v4-flash-rerun-audit-2026-08-30"
         / "FULL-RERUN-VALIDATION.json").read_text()
    )
    expected = {job["path"]: job["sha256"] for job in validation["jobs"]}
    assert stage.CANONICAL_V4_FILES == expected
    assert len(expected) == 10
    for relative, expected_sha in expected.items():
        path = ROOT / relative
        assert path.exists(), relative
        assert stage.sha256(path) == expected_sha


def test_release_stages_all_canonical_v4_sources(tmp_path):
    stage = load_stage_module()
    result = stage.stage_release(tmp_path, tag=None, version_doi="ZENODO_VERSION_DOI_PENDING")
    assert not result.missing
    staged = {item.path for item in result.files}
    assert {Path(path).name for path in stage.CANONICAL_V4_FILES} <= staged
    assert "RELEASE-MANIFEST.json" in staged
    assert "CHECKSUMS.sha256" in staged
    manifest = json.loads((tmp_path / "RELEASE-MANIFEST.json").read_text())
    assert manifest["release"]["tag"] == "agentic-pgx-benchmark-v2.9"
    assert manifest["release"]["model_label"] == "DeepSeek V4 Flash"
    assert manifest["canonical_v4_validation"]["calls"] == 8901


def test_current_metadata_names_the_reserved_v4_release_identity():
    combined = "\n".join(
        (ROOT / path).read_text()
        for path in ("README.md", "data/README.md", "data/MANIFEST.md")
    )
    assert "zenodo.22141182" in combined
    assert "agentic-pgx-benchmark-v2.9" in combined
    assert "agentic-pgx-benchmark-v2.8" not in combined
    assert "DeepSeek V3" not in (ROOT / "MODEL-VERSIONS.md").read_text()
    assert "deepseek-chat" not in (ROOT / "MODEL-VERSIONS.md").read_text()
