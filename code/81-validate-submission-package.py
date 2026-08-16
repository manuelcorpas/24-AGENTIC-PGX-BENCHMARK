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
# The tag and version DOI are arguments, not constants. They were constants,
# and the file therefore validated a package that no longer existed (v55 / v32 /
# v37) while the submitted package was v70 / v41 / v63. A validator pinned to
# superseded filenames does not fail loudly, it fails to run, which is worse.
DEFAULT_TAG = "agentic-pgx-benchmark-v2.4"


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


def latest(directory: Path, stem: str) -> Path:
    """Resolve the highest version of a document, rather than a pinned name."""
    hits = sorted(directory.glob(f"{stem}-v*.docx"),
                  key=lambda q: int(re.search(r"-v(\d+)\.docx$", q.name).group(1)))
    return hits[-1] if hits else directory / f"{stem}-vNONE.docx"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--package-dir", type=Path, required=True,
                    help="the submission pack itself, holding the uploaded files")
    ap.add_argument("--tag", default=DEFAULT_TAG)
    ap.add_argument("--version-doi", required=True,
                    help="Zenodo version DOI number, e.g. 21905387. Required "
                         "because the whole point is that the archive and the "
                         "package name the same thing.")
    args = ap.parse_args(argv)
    pack = args.package_dir
    figures = pack / "figures".upper()
    paths = {
        "manuscript": latest(pack, "Cell-Genomics-Manuscript"),
        "supplement": latest(pack, "Cell-Genomics-Supplementary"),
        "response_docx": latest(pack, "Response-to-Reviewers"),
    }
    failures = []
    for label, path in paths.items():
        if not path.exists():
            failures.append(f"missing {label}: {path}")
    if failures:
        print("\n".join(f"FAIL: {x}" for x in failures))
        return 1
    print("validating: " + ", ".join(p.name for p in paths.values()))

    manuscript = docx_text(paths["manuscript"])
    supplement = docx_text(paths["supplement"])
    response = docx_text(paths["response_docx"])
    response_md = ""
    all_text = manuscript + supplement + response

    required = {
        "manuscript": [
            "Validated skills are necessary but not sufficient for trustworthy agentic genomics",
            "13,200 attempted evaluations",
            "13,199",
            args.tag,
        ],
        "supplement": ["configuration", args.tag],
        "response": ["matched five-configuration comparison", args.tag],
    }
    haystacks = {"manuscript": manuscript, "supplement": supplement,
                 "response": response}
    for group, needles in required.items():
        for needle in needles:
            if needle not in haystacks[group]:
                failures.append(f"{group}: required text missing: {needle!r}")

    retired = [
        "iTrustworthy", "one of three models", "two of three models",
        "Gemini 2.5 Flash is excluded", "ZENODO_VERSION_DOI_PENDING",
        "public tag has not been pushed",
        # The vocabulary retired on 16 August. "cell" survives legitimately in
        # the journal name and in one deposited filename, so this looks only for
        # the arm sense, which is the one that was renamed.
        "five-cell comparison", "per cell", "execution cells",
        # Superseded release identifiers: every benchmark tag except the one
        # this run is validating against. Hardcoding v2.3 here made the check
        # contradict its own --tag argument.
        *[t for t in ("agentic-pgx-benchmark-v2.2", "agentic-pgx-benchmark-v2.3",
                      "agentic-pgx-benchmark-v2.4") if t != args.tag],
    ]
    for needle in retired:
        if needle and needle.lower() in all_text.lower():
            failures.append(f"retired text survives: {needle!r}")

    for dash, label in (("\u2014", "em dash"), ("\u2013", "en dash")):
        if dash in all_text:
            failures.append(f"{label} present in submission text")

    version_dois = set(re.findall(r"zenodo\.(\d+)", all_text))
    # Deposits that are legitimately cited and are not this paper's archive.
    # Sweeping every zenodo.NNN in the text flagged the ClawBio skill library as
    # a synchronisation failure, which is a false accusation of a correct
    # citation, so each exclusion is named rather than pattern-matched.
    version_dois.discard("20567742")          # this paper's concept DOI
    version_dois.discard("19420648")          # ClawBio skill library v0.5.0
    if version_dois != {args.version_doi}:
        failures.append(f"version DOI not synchronized: found {sorted(version_dois)}, "
                        f"expected {args.version_doi!r}")

    for path in paths.values():
        failures.extend(structural_docx_failures(path))

    ms_media = embedded_hashes(paths["manuscript"])
    supp_media = embedded_hashes(paths["supplement"])
    for name, media in (("Figure{}", ms_media), ("FigureS{}", supp_media)):
        count = 4 if name == "Figure{}" else 8
        for n in range(1, count + 1):
            candidates = [figures / (name.format(n) + ext)
                          for ext in (".png", ".tiff")]
            found = [c for c in candidates if c.exists()]
            if not found:
                failures.append(f"missing separate copy of {name.format(n)}")
            elif not any(sha(c) in media for c in found):
                failures.append(f"{name.format(n)} separate copy is not the embedded one")

    if failures:
        for failure in failures:
            print("FAIL:", failure)
        print(f"\n{len(failures)} package checks failed")
        return 1
    print(f"PASS: {paths['manuscript'].name}, {paths['supplement'].name}, "
          f"{paths['response_docx'].name} and all twelve figures agree")
    print(f"PASS: release tag {args.tag}; version DOI zenodo.{args.version_doi}; "
          "no stale claim, comment or tracked change")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
