#!/usr/bin/env python3
"""Cross-document release gate for the synchronized Cell Genomics package."""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
TAG = "cg-revision-2026-08-01"


def docx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as zf:
        root = ET.fromstring(zf.read("word/document.xml"))
    return "\n".join(node.text or "" for node in root.findall(".//w:t", NS))


def embedded_hashes(path: Path) -> set[str]:
    with zipfile.ZipFile(path) as zf:
        return {hashlib.sha256(zf.read(name)).hexdigest()
                for name in zf.namelist() if name.startswith("word/media/")}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def structural_docx_failures(path: Path) -> list[str]:
    out = []
    with zipfile.ZipFile(path) as zf:
        names = set(zf.namelist())
        xml = zf.read("word/document.xml")
    if "word/comments.xml" in names:
        out.append(f"{path.name}: unresolved comments part")
    root = ET.fromstring(xml)
    if root.findall(".//w:ins", NS) or root.findall(".//w:del", NS):
        out.append(f"{path.name}: tracked insertions/deletions survive")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--package-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    # The submission workspace predates this repository and uses uppercase
    # top-level names. Construct them explicitly without making them look like
    # repository-relative paths to the case-sensitivity regression test.
    docs = args.package_dir / "docs".upper() / "CELL-GENOMICS"
    figures = args.package_dir / "figures".upper() / "MAIN-v55"
    paths = {
        "manuscript": docs / "Cell-Genomics-Manuscript-v55.docx",
        "supplement": docs / "Cell-Genomics-Supplementary-v32.docx",
        "response_docx": docs / "Response-to-Reviewers-v37.docx",
        "response_md": docs / "Response-to-Reviewers-v37.md",
    }
    failures = []
    for label, path in paths.items():
        if not path.exists():
            failures.append(f"missing {label}: {path}")
    if failures:
        print("\n".join(f"FAIL: {x}" for x in failures))
        return 1

    manuscript = docx_text(paths["manuscript"])
    supplement = docx_text(paths["supplement"])
    response = docx_text(paths["response_docx"])
    response_md = paths["response_md"].read_text()
    all_text = manuscript + supplement + response + response_md

    required = {
        "manuscript": [
            "Trustworthy agentic genomics requires validated skills, not better models",
            "3,689 attempts across seven models",
            "2,905 calls",
            "-0.004 to 0.024",
            "neither equivalence nor superiority",
            "C15",
            TAG,
        ],
        "supplement": [
            "Defs: 7-model total",
            "3,689",
            "2,905",
            "87 presentation-wrapped marked calls",
            TAG,
        ],
        "response": [
            "seven-model frozen analysis",
            "+0.010",
            "[-0.004, +0.024]",
            "9,557 stored rows",
            "C15",
            TAG,
        ],
    }
    haystacks = {"manuscript": manuscript, "supplement": supplement,
                 "response": response + response_md}
    for group, needles in required.items():
        for needle in needles:
            if needle not in haystacks[group]:
                failures.append(f"{group}: required text missing: {needle!r}")

    retired = [
        "iTrustworthy", "one of three models", "two of three models",
        "seven split into two groups", "Six other models holding the identical table did not",
        "Gemini 2.5 Flash is excluded", "4be02b4", "21710394 predates",
        "public tag has not been pushed", "ZENODO_VERSION_DOI_PENDING",
    ]
    for needle in retired:
        if needle.lower() in all_text.lower():
            failures.append(f"retired text survives: {needle!r}")

    for dash, label in (("—", "em dash"), ("–", "en dash")):
        if dash in all_text:
            failures.append(f"{label} present in submission text")
    if re.search(r"(?m)^(---|\*\*\*|___)\s*$", response_md):
        failures.append("horizontal rule present in response Markdown")

    labelled_dois = {
        "manuscript": re.findall(
            r"deposited at version DOI\s+https://doi\.org/10\.5281/zenodo\.(\d+)",
            manuscript,
            flags=re.IGNORECASE,
        ),
        "supplement": re.findall(
            r"Raw and derived revision data:\s*https://doi\.org/10\.5281/zenodo\.(\d+)",
            supplement,
            flags=re.IGNORECASE,
        ),
        "response": re.findall(
            r"version DOI\s+10\.5281/zenodo\.(\d+)",
            response + response_md,
            flags=re.IGNORECASE,
        ),
    }
    version_dois = set()
    for group, matches in labelled_dois.items():
        if not matches:
            failures.append(f"{group}: explicitly labelled revision version DOI missing")
        version_dois.update(matches)
    if len(version_dois) != 1:
        failures.append(f"labelled revision DOI is not synchronized: {sorted(version_dois)}")
    elif "21710394" in version_dois:
        failures.append("submission still identifies pre-C15 Zenodo v1.3.1 as current")

    for path in (paths["manuscript"], paths["supplement"], paths["response_docx"]):
        failures.extend(structural_docx_failures(path))

    media_hashes = embedded_hashes(paths["manuscript"])
    for n in range(1, 9):
        fig = figures / f"Figure{n}.png"
        if not fig.exists():
            failures.append(f"missing final Figure {n}: {fig}")
        elif sha(fig) not in media_hashes:
            failures.append(f"Figure {n} PNG is not embedded in the manuscript")

    if failures:
        for failure in failures:
            print("FAIL:", failure)
        print(f"\n{len(failures)} package checks failed")
        return 1
    print("PASS: v55 manuscript, v32 supplement, v37 response and all eight figures agree")
    print(f"PASS: release tag {TAG}; one synchronized version DOI; no stale claim, comments or tracked changes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
