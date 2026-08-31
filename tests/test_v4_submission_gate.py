import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_gate():
    path = ROOT / "code" / "81-validate-submission-package.py"
    spec = importlib.util.spec_from_file_location("submission_gate", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_gate_defaults_to_v4_release_identity():
    gate = load_gate()
    assert gate.DEFAULT_TAG == "agentic-pgx-benchmark-v2.10"
    assert gate.RELEASE_VERSION == "1.8.0"
    assert gate.REQUIRED_MODEL == "DeepSeek V4 Flash"


def test_gate_resolves_current_submission_names(tmp_path):
    gate = load_gate()
    for name in (
        "00-Manuscript-v90.docx", "01-Supplementary-v48.docx",
        "02-Response-to-Reviewers-v68.docx", "03-Cover-Letter-v12.docx",
        "04-Declaration-of-Interests-v32.docx",
    ):
        (tmp_path / name).touch()
    paths = gate.resolve_submission_documents(tmp_path)
    assert paths["manuscript"].name == "00-Manuscript-v90.docx"
    assert paths["supplement"].name == "01-Supplementary-v48.docx"
    assert paths["response"].name == "02-Response-to-Reviewers-v68.docx"
    assert paths["cover"].name == "03-Cover-Letter-v12.docx"
    assert paths["declaration"].name == "04-Declaration-of-Interests-v32.docx"


def test_release_manifest_gate_rejects_wrong_hash(tmp_path):
    gate = load_gate()
    (tmp_path / "payload.txt").write_text("changed")
    manifest = {
        "release": {"tag": gate.DEFAULT_TAG, "version": gate.RELEASE_VERSION,
                    "version_doi": "10.5281/zenodo.123"},
        "files": [{"path": "payload.txt", "sha256": "0" * 64, "size": 7}],
        "missing": [],
    }
    failures = gate.validate_release_manifest(tmp_path, manifest, "123")
    assert any("checksum" in failure for failure in failures)
