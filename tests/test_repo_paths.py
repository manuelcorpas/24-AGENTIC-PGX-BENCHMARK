#!/usr/bin/env python3
"""
Path-integrity tests for the reproducibility package.

Why these exist: every script derives its paths from the repository root, but
several referenced the directories in a different case than they exist on disk
(BASE / "RESULTS" against a real data/, BASE / "FIGURES" against a real
figures/). On macOS the filesystem is case-insensitive, so this resolved
silently on the authors' machines and failed with FileNotFoundError on a
reviewer's case-sensitive Linux checkout.

These tests are deliberately case-sensitive even when run on a case-insensitive
filesystem: they compare against os.listdir() of the real directory rather than
asking the filesystem whether a path exists.
"""
import os
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

# Directories the scripts create themselves at run time (mkdir), so they are not
# tracked in git and will not be present in a fresh clone.
RUNTIME_CREATED = {"logs"}

# Tokens that are module/attribute names rather than directories.
NOT_DIRECTORIES = {"__file__"}

# BASE / "X", ROOT / "X", PROJECT_ROOT / "X", ...
DIR_TOKEN = re.compile(r'(?:BASE|ROOT|REPO|PROJECT_ROOT)\s*/\s*"([A-Za-z0-9_.\-]+)"')

# The legacy uppercase names this normalisation removed. Guards against regression.
LEGACY_UPPERCASE = {"RESULTS", "FIGURES", "LOGS", "DOCS"}


def python_files():
    """Every script in the package, excluding this test (it names the legacy strings)."""
    for sub in ("code", "real-genome-arm", "tests"):
        for path in sorted((REPO / sub).rglob("*.py")):
            if path.resolve() != Path(__file__).resolve():
                yield path


def referenced_tokens():
    """(token, relative path of the file that referenced it) for every script."""
    for path in python_files():
        for token in DIR_TOKEN.findall(path.read_text(encoding="utf-8", errors="replace")):
            if token in NOT_DIRECTORIES or "." in token:
                continue
            yield token, path.relative_to(REPO)


def test_referenced_directories_exist_with_exact_case():
    """Every root-level directory a script names must exist with that exact case.

    Compares against os.listdir() so the assertion is case-sensitive even on
    macOS, which is the whole point: this must fail here if it would fail on a
    reviewer's Linux checkout.
    """
    actual = set(os.listdir(REPO))
    missing = sorted(
        {
            (token, str(src))
            for token, src in referenced_tokens()
            if token not in actual and token not in RUNTIME_CREATED
        }
    )
    assert not missing, "directories referenced by scripts but absent (exact case): " + repr(missing)


def test_no_legacy_uppercase_directory_names():
    """The uppercase names are gone, so the case mismatch cannot silently return.

    Scans every quoted occurrence, not just the BASE / "X" form, so os.path.join
    call sites and docstring path listings are covered too.
    """
    offenders = []
    for path in python_files():
        for i, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            for name in LEGACY_UPPERCASE:
                if re.search(rf'["\']{name}["\']|(?<![A-Za-z0-9_-]){name}/', line):
                    offenders.append((name, f"{path.relative_to(REPO)}:{i}"))
    assert not offenders, "legacy uppercase directory names still referenced: " + repr(offenders)


def test_runtime_created_directories_are_actually_created():
    """Anything exempted as runtime-created must really be mkdir'd by a script."""
    for name in RUNTIME_CREATED:
        creators = [
            p.relative_to(REPO)
            for p in python_files()
            if re.search(rf'"{name}"', p.read_text(encoding="utf-8", errors="replace"))
            and "mkdir" in p.read_text(encoding="utf-8", errors="replace")
        ]
        assert creators, f"{name}/ is exempted as runtime-created but no script mkdirs it"


@pytest.mark.parametrize("path", list(python_files()), ids=lambda p: p.name)
def test_no_hardcoded_absolute_home_paths(path):
    """No script may point at a personal machine; a reviewer's checkout has no /Users/<name>."""
    text = path.read_text(encoding="utf-8", errors="replace")
    hits = re.findall(r'["\'](/Users/[^"\']+|/home/[^"\']+)["\']', text)
    assert not hits, f"{path.name} hardcodes an absolute home path: {hits}"
