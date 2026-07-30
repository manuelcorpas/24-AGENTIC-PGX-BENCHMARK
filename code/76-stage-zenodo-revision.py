#!/usr/bin/env python3
"""
Stage the Zenodo v1.3.0 deposit for the revision, and record its provenance.

WHY THIS EXISTS
Two outstanding items close together here.

R2.6 asked that the reproducibility package match the manuscript. The raw
evaluation rows are gitignored by design (they are large, and data/README.md
sends a reader to Zenodo for them), so "the repository matches the manuscript"
is only true once the revision's raw data is actually deposited. Until then the
tag carries the code and the DOI carries the previous version's data.

The second item is provenance. The raw rows carry per-call input tokens, output
tokens and cost, but no timestamps and no provider response IDs: they were never
captured and cannot be reconstructed. A Zenodo deposit does not invent them, but
it does supply something the project currently lacks entirely, which is an
EXTERNALLY ATTESTED DATE. After deposit, the raw rows are bounded above by a
publication date Zenodo records and we do not control. Combined with the
provider billing period, which bounds the spend, that is the strongest
provenance claim the surviving evidence supports, and it is weaker than
per-call timestamps would have been. This script writes that statement rather
than leaving the letter to overclaim.

USAGE
    python code/76-stage-zenodo-revision.py --out ~/Desktop/ZENODO-v1.3.0
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"

# The raw and derived files the revision's reported numbers are computed from.
# Anything a manuscript number depends on belongs here; anything else does not.
FILES = [
    "v3_five_cell_live.json",            # matched factorial, raw rows
    "v3_five_cell_live_rows.json",       # scored rows
    "v3_five_cell_live_scored.json",
    "v3_five_cell_live_stats.json",
    "v3_matched_scored_rows.json",
    "v3_extracted_rules.json",           # R2.7 extraction experiment
    "v3_extracted_rules_eval.json",
    "v3_getrm_caller_truth.json",        # GeT-RM external truth
    "v3_getrm_disagreement_classes.json",
    "v3_caller_truth_eval.json",
    "v3_ancestry_four_cohorts.json",     # four-cohort rerun
    "v3_ancestry_matched.json",
    "v3_agent_vs_deterministic.json",
    "v3_pharmcat_comparison.json",
    "v3_crossed_mixed_model.json",       # R2.5 mixed model
    "v3_hierarchical_stats.json",
    "v3_lethal_case_level.json",
    "v3_cell_provenance.json",
    "v3_normalisation_inputs.json",      # R1.1 input normalisation, frozen inputs
    "v3_input_normalisation_main.json",  # arm 1, Claude + GPT-5.2
    "v3_input_normalisation_o3.json",    # arm 1, o3 subsample
    "v3_input_normalisation_defs.json",  # arm 2, the result that changed the conclusion
    "v3_input_normalisation_o3_retry.json",  # the 429 batch: all error, deposited anyway
    "v3_input_normalisation_eval.json",
    "v3_normalisation_disagreement_classes.json",
]

PROVENANCE = """# Provenance of the raw evaluation rows

## What every row carries

Each evaluation row records the model it was sent to, the exact prompt inputs
(or the hash that freezes them), the raw response text, input tokens, output
tokens and cost in USD. Per-cell totals regenerate from these rows via
`code/69-cell-provenance.py`.

## What no row carries, and why it cannot be added

No row carries a wall-clock timestamp or a provider response ID. These were not
captured at run time. They are not reconstructable after the fact: the provider
APIs used here do not expose a per-response lookup keyed by content, and the
run logs that might have carried them were not retained. We state this rather
than reconstructing an approximate time from file modification dates, which
would be a local artefact presented as evidence.

## What bounds the runs instead

Two external bounds, neither of which we control:

1. **This Zenodo deposit.** Its publication date is recorded by Zenodo and
   places an upper bound on when the raw rows were produced.
2. **The provider billing period.** The per-row costs sum to the totals in
   `data/v3_cell_provenance.txt`; the corresponding provider invoices place the
   spend, and therefore the calls, inside a billing window that the provider
   dates independently of us.

Together these bound the runs from outside the project. They are weaker than
per-call timestamps, and we do not present them as equivalent.
"""


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    staged, missing = [], []
    for name in FILES:
        src = DATA / name
        if not src.exists():
            missing.append(name)
            continue
        shutil.copy2(src, args.out / name)
        staged.append((name, src.stat().st_size, sha256(src)))

    lines = ["# SHA-256 checksums, Zenodo v1.3.0 (revision dataset)", ""]
    lines += [f"{h}  {n}" for n, _, h in staged]
    (args.out / "CHECKSUMS.sha256").write_text("\n".join(lines) + "\n")
    (args.out / "PROVENANCE.md").write_text(PROVENANCE)

    total = sum(s for _, s, _ in staged)
    print(f"staged {len(staged)} files, {total/1e6:.1f} MB -> {args.out}")
    for n, s, h in staged:
        print(f"  {n:44s} {s/1e6:7.2f} MB  {h[:16]}")
    if missing:
        print(f"\nMISSING ({len(missing)}), not staged:")
        for n in missing:
            print(f"  {n}")
        print("\nDeposit is incomplete until these exist. Not writing a manifest "
              "that claims otherwise.")
        return 1
    print("\nAll files present. Upload this directory to Zenodo as a new version.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
