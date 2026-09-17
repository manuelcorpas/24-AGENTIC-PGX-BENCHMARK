#!/usr/bin/env python3
"""Build the blinded input package handed to Claude Science.

Removes every ground-truth field from the 110 CPIC cases and copies the same
CPIC guideline corpus the RAG arm used. Answers stay behind in SPECS/.

Usage: python PYTHON/00-build_package.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

ALLOWED = ("id", "gene", "drug", "genotype")
KNOWN_HIDDEN = {"gt_diplotype", "gt_phenotype", "gt_drug", "gt_activity", "pop_note"}


def hashed_id(case_id: str) -> str:
    """Stable 12-hex-character surrogate. Hex only, so no identifier can spell
    a phenotype hint such as 'pm' or 'um'."""
    return hashlib.sha256(f"cpic-study::{case_id}".encode()).hexdigest()[:12]


def id_map(cases: list[dict]) -> dict[str, str]:
    """Surrogate id -> real case id, kept out of the handed-over package."""
    return {hashed_id(c["id"]): c["id"] for c in cases}


def blind_cases(cases: list[dict], hash_ids: bool = False) -> list[dict]:
    """Keep only the inputs. Any field not known to be safe or hidden raises,
    so a new answer field added upstream cannot leak silently.

    The real case identifiers encode the intended phenotype (cyp2d6_codeine_pm),
    which the 2026-09-15 Claude Science run found and disclosed. hash_ids
    replaces them with surrogates so a rerun is blind on that channel too.
    """
    out = []
    for c in cases:
        unknown = set(c) - set(ALLOWED) - KNOWN_HIDDEN
        if unknown:
            raise ValueError(f"unclassified case fields {sorted(unknown)} in {c.get('id')}")
        row = {k: c[k] for k in ALLOWED}
        if hash_ids:
            row["id"] = hashed_id(c["id"])
        out.append(row)
    return out


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    here = Path(__file__).resolve().parents[1]
    specs = here.parent / "SPECS"
    pkg = here / "PACKAGE"
    pkg.mkdir(exist_ok=True)

    cases = json.loads((specs / "test_cases_v3.json").read_text())
    (pkg / "cases_blinded.json").write_text(json.dumps(blind_cases(cases), indent=1))
    shutil.copyfile(specs / "cpic_rag_corpus_v3.json", pkg / "cpic_guideline_corpus.json")

    manifest = {p.name: sha256(p) for p in sorted(pkg.glob("*.json")) if p.name != "MANIFEST.json"}
    manifest["n_cases"] = len(cases)
    (pkg / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    print(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
