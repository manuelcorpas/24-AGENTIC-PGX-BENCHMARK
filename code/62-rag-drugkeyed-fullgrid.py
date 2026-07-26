#!/usr/bin/env python3
"""
Drug-keyed retrieval on the full grid (revision item N4; R1.2, R2.3).

WHAT THE REVIEWERS SAID
R1.2 and R2.3: the paper's retrieval conclusions are broader than its baseline
supports. The retrieval arm was indexed gene-keyed, so a query about codeine
returned the entire CYP2D6 annotation set (24,781 characters covering codeine,
tamoxifen, ondansetron, paroxetine and more) and the model had to select the
right drug's table itself. The drug-substitution failure is therefore partly a
property of that indexing choice, not of retrieval augmentation as such.

WHAT THIS DOES
Runs the retrieval cells with drug-keyed chunking on the full grid, so the
strengthened configuration becomes the primary retrieval arm and the gene-keyed
configuration becomes a declared ablation. 16b-rag-genedrug-chunking.py
established the effect on six genes and 378 cells; this extends it to all 110
cases, all models and both retrieval cells (generation and extract-then-execute).

RETRIEVAL ROUTES, REPORTED NOT HIDDEN
Every row records how its excerpt was retrieved:
    drug_keyed_exact  101 cases: the corpus has a section for that drug
    drug_keyed_class    9 cases: the case queries a drug CLASS (thiopurines,
                        aminoglycosides, volatile-anaesthetics, PEG-IFN-alpha)
                        and the member sections are concatenated
    fallback_gene_keyed no section found; gene-level excerpt used
A run that reported drug-keyed coverage without these counts would be
misdescribing the experiment, so the report prints the breakdown.

COST
Drug-keyed chunks average 5,148 characters against 12,915 gene-level, so the
stronger arm is also the cheaper one.

USAGE
    python code/62-rag-drugkeyed-fullgrid.py --estimate
    python code/62-rag-drugkeyed-fullgrid.py --pilot 3
    python code/62-rag-drugkeyed-fullgrid.py --max-spend 40
"""
from __future__ import annotations

import collections
import json
import sys
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

CODE = Path(__file__).resolve().parent
BASE = CODE.parent

_spec = spec_from_file_location("matched_factorial", CODE / "60-matched-factorial.py")
mf = module_from_spec(_spec)
_spec.loader.exec_module(mf)

RETRIEVAL_CELLS = ["rag_generation", "rag_execution"]
OUT = BASE / "data" / "v3_rag_drugkeyed_fullgrid.json"


def route_census() -> dict[str, int]:
    """How the 110 cases resolve under drug-keyed retrieval. No API calls."""
    previous = mf.RETRIEVAL_MODE
    mf.RETRIEVAL_MODE = "drug"
    try:
        census = collections.Counter(
            mf.retrieve(c["gene"], c["drug"])[1] for c in mf.CASES)
    finally:
        mf.RETRIEVAL_MODE = previous
    return dict(census)


def main(argv: list[str] | None = None) -> int:
    argv = list(argv if argv is not None else sys.argv[1:])

    census = route_census()
    print("drug-keyed retrieval routes across the 110 cases:")
    for route, n in sorted(census.items(), key=lambda kv: -kv[1]):
        print(f"    {route:<22} {n}")
    if census.get("fallback_gene_keyed"):
        print(f"    NOTE: {census['fallback_gene_keyed']} cases fell back to the "
              "gene-level excerpt and are flagged in every row.")
    print()

    # Force drug-keyed retrieval and the two retrieval cells; everything else
    # (matched prompts, spend cap, resume, pilot) is inherited from 60-.
    if "--retrieval" not in argv:
        argv += ["--retrieval", "drug"]
    if "--cells" not in argv:
        argv += ["--cells", *RETRIEVAL_CELLS]

    mf.OUT = OUT
    mf.JSONL = OUT.with_suffix(".jsonl")
    return mf.main(argv)


if __name__ == "__main__":
    sys.exit(main())
