#!/usr/bin/env python3
"""Stage the hash-bound DeepSeek V4 Flash Zenodo release package."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
AUDIT = BASE / "docs" / "deepseek-v4-flash-rerun-audit-2026-08-30"
RELEASE_VERSION = "1.8.0"
RELEASE_TAG = "agentic-pgx-benchmark-v2.10"
MODEL_LABEL = "DeepSeek V4 Flash"
REQUEST_MODEL = "deepseek-v4-flash"
RETURNED_MODEL = "deepseek-v4-flash"
THINKING_MODE = "disabled"

_VALIDATION = json.loads((AUDIT / "FULL-RERUN-VALIDATION.json").read_text())
CANONICAL_V4_FILES = {job["path"]: job["sha256"] for job in _VALIDATION["jobs"]}

# Historical non-DeepSeek arms reconstruct the common panel. V4 merged and
# derived files replace only the DeepSeek arm in the reported analyses.
SUPPORTING_DATA_FILES = [
    "data/v3_five_cell_live.json", "data/v3_five_cell_live_rows.json",
    "data/v3_five_cell_live_scored.json", "data/v3_five_cell_live_stats.json",
    "data/v3_cell_provenance.json", "data/v3_model_caller_eval.json",
    "data/v3_matched_scored_rows.json", "data/v3_matched_scored_rows_all5.json",
    "data/v3_extracted_rules.json", "data/v3_getrm_caller_truth.json",
    "data/v3_getrm_disagreement_classes.json", "data/v3_caller_truth_eval.json",
    "data/v3_ancestry_four_cohorts.json", "data/v3_ancestry_matched.json",
    "data/v3_pharmcat_comparison.json", "data/v3_normalisation_inputs.json",
    "data/v3_input_normalisation_main.json", "data/v3_input_normalisation_o3.json",
    "data/v3_input_normalisation_defs.json", "data/v3_input_normalisation_o3_retry.json",
    "data/v3_input_normalisation_eval.json", "data/v3_normalisation_disagreement_classes.json",
    "data/v3_input_normalisation_defs_tail.json", "data/v3_input_normalisation_defs_gpt52.json",
    "data/v3_input_normalisation_defs_o3.json",
    "data/v3_input_normalisation_defs_sonnet_gpt41_o4mini.json",
    "data/v3_input_normalisation_defs_sonnet_gpt41_o4mini_rest.json",
    "data/v3_input_normalisation_nodefs_four.json", "data/v3_input_normalisation_nodefs_gemini.json",
    "data/v3_input_normalisation_nodefs_o3_527.json", "data/v3_input_normalisation_defs_gemini.json",
    "data/v3_rag_genedrug_chunking.json", "data/v3_armA9_armBv2_POP.json",
    "data/v4_flash_eight_model_five_cell.json", "data/v4_flash_five_cell_scored.json",
    "data/v4_flash_matched_scored_rows_all5.json", "data/v4_flash_hierarchical_stats.json",
    "data/v4_flash_crossed_mixed_model.json", "data/v4_flash_eight_model_extracted_rules.json",
    "data/v4_flash_lethal_case_stats.json",
    "data/v4_flash_extracted_rules_eval.json", "data/v4_flash_eight_model_rag_drugkeyed.json",
    "data/v4_flash_rag_drugkeyed_scored_rows.json", "data/v4_flash_rag_drugkeyed_summary.json",
    "data/v4_flash_adversarial_forward_panel.json", "data/v4_flash_adversarial_reverse_panel.json",
    "data/v4_flash_population_three_arm_panel.json", "data/v4_flash_population_skill_arms_panel.json",
    "data/v4_flash_input_normalisation_eight_model_freeze.json",
    "data/v4_flash_eight_model_realgenome_provider_responses.jsonl",
    "data/v4_flash_realgenome_scored.jsonl", "data/v4_flash_realgenome_model_audit.json",
    "data/v4_flash_agent_vs_deterministic.json",
]
TEXT_FILES = [
    "CORRECTIONS.md", "MODEL-VERSIONS.md", "data/MANIFEST.md", "data/README.md",
    "data/v4_flash_five_cell_report.txt", "data/v4_flash_extracted_rules_eval.txt",
    "data/v4_flash_hierarchical_stats.txt", "data/v4_flash_crossed_mixed_model.txt",
    "data/v4_flash_lethal_case_stats.txt", "data/v4_flash_rag_drugkeyed_report.txt",
    "data/v4_flash_realgenome_model_audit.txt", "data/v4_flash_supp_tables_s1_s3.txt",
    "docs/deepseek-v4-flash-rerun-audit-2026-08-30/RERUN-MANIFEST.json",
    "docs/deepseek-v4-flash-rerun-audit-2026-08-30/FULL-RERUN-VALIDATION.json",
    "docs/deepseek-v4-flash-rerun-audit-2026-08-30/SOURCE-HASHES.tsv",
]

@dataclass(frozen=True)
class StagedFile:
    path: str
    size: int
    sha256: str
    role: str

@dataclass(frozen=True)
class StageResult:
    files: tuple[StagedFile, ...]
    missing: tuple[str, ...]

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def _copy(relative: str, out: Path, role: str) -> StagedFile | None:
    source = BASE / relative
    if not source.exists():
        return None
    destination = out / Path(relative).name
    shutil.copy2(source, destination)
    return StagedFile(destination.name, destination.stat().st_size, sha256(destination), role)

def _verify_canonical_sources() -> None:
    failures = []
    for relative, expected in CANONICAL_V4_FILES.items():
        path = BASE / relative
        actual = sha256(path) if path.exists() else None
        if actual != expected:
            failures.append(f"{relative}: expected {expected}, found {actual}")
    if failures:
        raise ValueError("Canonical V4 source gate failed:\n" + "\n".join(failures))

def stage_release(out: Path, tag: str | None, version_doi: str) -> StageResult:
    _verify_canonical_sources()
    out.mkdir(parents=True, exist_ok=True)
    staged: list[StagedFile] = []
    missing: list[str] = []
    for relative in CANONICAL_V4_FILES:
        item = _copy(relative, out, "canonical V4 raw provider output")
        (staged if item else missing).append(item or relative)
    for relative in SUPPORTING_DATA_FILES:
        if relative in CANONICAL_V4_FILES:
            continue
        item = _copy(relative, out, "supporting or derived analysis input")
        (staged if item else missing).append(item or relative)
    for relative in TEXT_FILES:
        item = _copy(relative, out, "release metadata or readable result")
        (staged if item else missing).append(item or relative)
    if tag:
        archive_name = f"24-AGENTIC-PGX-BENCHMARK_{tag}.zip"
        archive = out / archive_name
        subprocess.run(["git", "archive", "--format=zip", f"--output={archive}", tag],
                       cwd=BASE, check=True)
        staged.append(StagedFile(archive_name, archive.stat().st_size,
                                 sha256(archive), "immutable tagged source archive"))
    provenance = out / "PROVENANCE.md"
    provenance.write_text(
        "# Provenance\n\n"
        "The ten canonical DeepSeek V4 Flash files are listed with SHA-256 hashes "
        "in FULL-RERUN-VALIDATION.json. All 8,901 requested units are present "
        "exactly once, with zero provider errors. Every canonical request used "
        "`deepseek-v4-flash`, required the returned identifier "
        "`deepseek-v4-flash`, and disabled thinking.\n\n"
        "The raw records retain the exact request payload, complete provider "
        "response, returned model identifier, raw text, stop reason, token usage, "
        "response identifier, timestamps, latency, retry history and error body "
        "where supplied by the provider and runner. SOURCE-HASHES.tsv freezes the "
        "benchmark inputs. RELEASE-MANIFEST.json maps every deposited file to its "
        "size and SHA-256 checksum.\n"
    )
    staged.append(StagedFile(provenance.name, provenance.stat().st_size,
                             sha256(provenance), "V4 run provenance"))
    notes = out / f"RELEASE-NOTES-v{RELEASE_VERSION}.md"
    notes.write_text(
        f"# Version {RELEASE_VERSION}\n\n"
        f"This release replaces every reported DeepSeek V3 arm with {MODEL_LABEL}, "
        f"requested and returned as `{REQUEST_MODEL}` with thinking disabled. The "
        f"ten canonical raw files contain {_VALIDATION['canonical_units_unique']:,} "
        f"unique calls with zero provider errors. It is paired with immutable tag "
        f"`{RELEASE_TAG}` and version DOI `10.5281/zenodo.{version_doi}`.\n\n"
        "Figure 3A contains all eight requested models. Figure 3B pools the seven "
        "models with version-matched curated and real-genome results because no "
        "matched curated V4 Flash result exists.\n"
    )
    staged.append(StagedFile(notes.name, notes.stat().st_size, sha256(notes), "release notes"))
    manifest_path = out / "RELEASE-MANIFEST.json"
    manifest = {
        "schema_version": "agentic-pgx-release-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "release": {"version": RELEASE_VERSION, "tag": RELEASE_TAG,
                    "version_doi": f"10.5281/zenodo.{version_doi}",
                    "model_label": MODEL_LABEL, "request_model": REQUEST_MODEL,
                    "returned_model": RETURNED_MODEL, "thinking_mode": THINKING_MODE},
        "canonical_v4_validation": {
            "calls": _VALIDATION["canonical_units_unique"],
            "provider_errors": _VALIDATION["provider_errors"],
            "source": "docs/deepseek-v4-flash-rerun-audit-2026-08-30/FULL-RERUN-VALIDATION.json"},
        "files": [item.__dict__ for item in staged], "missing": missing,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    staged.append(StagedFile(manifest_path.name, manifest_path.stat().st_size,
                             sha256(manifest_path), "machine-readable release manifest"))
    checksums = out / "CHECKSUMS.sha256"
    checksums.write_text("".join(f"{item.sha256}  {item.path}\n" for item in staged))
    staged.append(StagedFile(checksums.name, checksums.stat().st_size,
                             sha256(checksums), "SHA-256 checksum ledger"))
    return StageResult(tuple(staged), tuple(missing))

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--tag")
    parser.add_argument("--version-doi", required=True)
    args = parser.parse_args(argv)
    result = stage_release(args.out, args.tag, args.version_doi)
    total = sum(item.size for item in result.files)
    print(f"staged {len(result.files)} files, {total / 1e6:.1f} MB -> {args.out}")
    if result.missing:
        for relative in result.missing:
            print(f"MISSING: {relative}")
        return 1
    print(f"PASS: {len(CANONICAL_V4_FILES)} canonical V4 files are hash-bound")
    print(f"PASS: {RELEASE_TAG} / Zenodo {RELEASE_VERSION} metadata staged")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
