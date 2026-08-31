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
DEFAULT_TAG = "agentic-pgx-benchmark-v2.10"
RELEASE_VERSION = "1.8.0"
REQUIRED_MODEL = "DeepSeek V4 Flash"


def docx_text(path: Path, prose_only: bool = False) -> str:
    """Visible text. With prose_only, skip monospaced runs.

    Why the option exists: the dash and style gates below police authored prose.
    A verbatim listing quoted from a deposited artefact (the authored rule table,
    set in a monospaced face) must be reproduced character for character, so
    normalising its punctuation to satisfy a prose rule would falsify the quote.
    Scoping the gate is correct; loosening it for prose would not be.
    """
    with zipfile.ZipFile(path) as zf:
        root = ET.fromstring(zf.read("word/document.xml"))
    if not prose_only:
        return "\n".join(node.text or "" for node in root.findall(".//w:t", NS))
    MONO = ("Consolas", "Courier New", "Courier", "Menlo", "Monaco")
    out = []
    for run in root.findall(".//w:r", NS):
        fonts = run.findall("./w:rPr/w:rFonts", NS)
        ascii_font = fonts[0].get(f"{{{NS['w']}}}ascii") if fonts else None
        if ascii_font in MONO:
            continue
        out.extend(node.text or "" for node in run.findall("./w:t", NS))
    return "\n".join(out)


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


def resolve_submission_documents(directory: Path) -> dict[str, Path]:
    """Resolve the numbered filenames used by the live submission package."""
    stems = {
        "manuscript": "00-Manuscript",
        "supplement": "01-Supplementary",
        "response": "02-Response-to-Reviewers",
        "cover": "03-Cover-Letter",
        "declaration": "04-Declaration-of-Interests",
    }
    return {label: latest(directory, stem) for label, stem in stems.items()}


def validate_release_manifest(directory: Path, manifest: dict,
                              version_doi: str) -> list[str]:
    failures = []
    release = manifest.get("release", {})
    expected = {
        "tag": DEFAULT_TAG,
        "version": RELEASE_VERSION,
        "version_doi": f"10.5281/zenodo.{version_doi}",
    }
    for key, value in expected.items():
        if release.get(key) != value:
            failures.append(f"release manifest {key}: {release.get(key)!r} != {value!r}")
    if manifest.get("missing"):
        failures.append(f"release manifest records missing files: {manifest['missing']}")
    for item in manifest.get("files", []):
        path = directory / item["path"]
        if not path.exists():
            failures.append(f"release file missing: {item['path']}")
        elif sha(path) != item["sha256"]:
            failures.append(f"release checksum mismatch: {item['path']}")
        elif path.stat().st_size != item["size"]:
            failures.append(f"release size mismatch: {item['path']}")
    return failures


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--package-dir", type=Path, required=True,
                    help="the submission pack itself, holding the uploaded files")
    ap.add_argument("--tag", default=DEFAULT_TAG)
    ap.add_argument("--version-doi", required=True,
                    help="Zenodo version DOI number, e.g. 21905387. Required "
                         "because the whole point is that the archive and the "
                         "package name the same thing.")
    ap.add_argument("--release-dir", type=Path, required=True,
                    help="staged Zenodo directory containing RELEASE-MANIFEST.json")
    args = ap.parse_args(argv)
    pack = args.package_dir
    figures = pack / "figures".upper()
    paths = resolve_submission_documents(pack)
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
    response = docx_text(paths["response"])
    cover = docx_text(paths["cover"])
    declaration = docx_text(paths["declaration"])
    all_text = manuscript + supplement + response + cover + declaration

    required = {
        "manuscript": [
            "Validated skills are necessary but not sufficient for trustworthy agentic genomics",
            "13,200 attempted evaluations",
            "13,199",
            args.tag, REQUIRED_MODEL,
        ],
        "supplement": ["configuration", args.tag],
        "response": ["matched five-configuration comparison", args.tag, REQUIRED_MODEL],
        "cover": [args.tag, REQUIRED_MODEL],
        "declaration": [args.tag],
    }
    haystacks = {"manuscript": manuscript, "supplement": supplement,
                 "response": response, "cover": cover, "declaration": declaration}
    for group, needles in required.items():
        for needle in needles:
            if needle not in haystacks[group]:
                failures.append(f"{group}: required text missing: {needle!r}")

    retired = [
        "iTrustworthy", "one of three models", "two of three models",
        "Gemini 2.5 Flash is excluded", "ZENODO_VERSION_DOI_PENDING",
        "public tag has not been pushed",
        "DeepSeek V3", "deepseek-chat",
        # The vocabulary retired on 16 August. "cell" survives legitimately in
        # the journal name and in one deposited filename, so this looks only for
        # the arm sense, which is the one that was renamed.
        "five-cell comparison", "per cell", "execution cells",
        # Superseded release identifiers: every benchmark tag except the one
        # this run is validating against. Hardcoding v2.3 here made the check
        # contradict its own --tag argument.
        *[t for t in ("agentic-pgx-benchmark-v2.2", "agentic-pgx-benchmark-v2.3",
                      "agentic-pgx-benchmark-v2.4", "agentic-pgx-benchmark-v2.8",
                      "agentic-pgx-benchmark-v2.9") if t != args.tag],
    ]
    for needle in retired:
        if needle and needle.lower() in all_text.lower():
            failures.append(f"retired text survives: {needle!r}")

    prose_text = (docx_text(paths["manuscript"], prose_only=True)
                  + docx_text(paths["supplement"], prose_only=True)
                  + docx_text(paths["response"], prose_only=True)
                  + docx_text(paths["cover"], prose_only=True)
                  + docx_text(paths["declaration"], prose_only=True))
    for dash, label in (("\u2014", "em dash"), ("\u2013", "en dash")):
        if dash in prose_text:
            failures.append(f"{label} present in submission prose")

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

    release_manifest = args.release_dir / "RELEASE-MANIFEST.json"
    if not release_manifest.exists():
        failures.append(f"missing release manifest: {release_manifest}")
    else:
        import json
        failures.extend(validate_release_manifest(
            args.release_dir, json.loads(release_manifest.read_text()),
            args.version_doi))

    # The separate figure copies are publication TIFFs built from the
    # full-resolution originals, so they are deliberately NOT byte-identical to
    # the downscaled images embedded in the documents. Byte-identity was the
    # right check while both were the same PNG and is the wrong one now, so
    # this compares the picture instead: same image, and at least 300 dpi
    # across a 180 mm page.
    # docxlib lives beside the submission pack, not in this repository, so it
    # is imported from there at run time rather than assumed to be importable.
    sys.path.insert(0, str(pack.parent))
    import docx
    from docxlib import pair_images_to_legends as dl_pair

    from io import BytesIO
    import zipfile as _zip

    from PIL import Image, ImageFilter
    import numpy as _np

    def _thumb(img):
        small = img.convert("L").resize((96, 96), Image.LANCZOS)
        return _np.asarray(small.filter(ImageFilter.GaussianBlur(1.2)), dtype=float)

    for doc_key, prefix, count in (("manuscript", "Figure", 4),
                                   ("supplement", "FigureS", 9)):
        doc = docx.Document(str(paths[doc_key]))
        pairs = dl_pair(doc)
        with _zip.ZipFile(paths[doc_key]) as z:
            for n in range(1, count + 1):
                legend = f"Figure {n}." if prefix == "Figure" else f"Figure S{n}."
                tiff = figures / f"{prefix}{n}.tiff"
                if not tiff.exists():
                    failures.append(f"missing publication TIFF: {tiff.name}")
                    continue
                if (figures / f"{prefix}{n}.png").exists():
                    failures.append(f"a PNG copy of {prefix}{n} survives beside the TIFF")
                embedded = Image.open(BytesIO(z.read("word/" + pairs[legend][1])))
                separate = Image.open(tiff)
                diff = float(_np.abs(_thumb(embedded) - _thumb(separate)).mean())
                if diff > 3.0:
                    failures.append(f"{prefix}{n}.tiff is not the figure embedded "
                                    f"under {legend} (difference {diff:.2f})")
                dpi = separate.size[0] / (180.0 / 25.4)
                if dpi < 300.0:
                    failures.append(f"{prefix}{n}.tiff is {dpi:.0f} dpi across 180 mm, under 300")

    if failures:
        for failure in failures:
            print("FAIL:", failure)
        print(f"\n{len(failures)} package checks failed")
        return 1
    print(f"PASS: {paths['manuscript'].name}, {paths['supplement'].name}, "
          f"{paths['response'].name}, {paths['cover'].name}, "
          f"{paths['declaration'].name} and all thirteen publication TIFFs agree")
    print(f"PASS: release tag {args.tag}; version DOI zenodo.{args.version_doi}; "
          "no stale claim, comment or tracked change")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
