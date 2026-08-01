#!/usr/bin/env python3
"""
Stage the synchronized Zenodo revision deposit and record its provenance.

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
    python code/76-stage-zenodo-revision.py --out ~/Desktop/ZENODO-v1.4.0 \
        --tag cg-revision-2026-08-01
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
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
    "v3_realgenome_preds_4cohorts.tsv", # Figure 7, four-cohort raw predictions
    "v3_pharmcat_comparison.json",
    "v3_crossed_mixed_model.json",       # R2.5 mixed model
    "v3_hierarchical_stats.json",
    "v3_lethal_case_level.json",
    "v3_cell_provenance.json",
    "v3_adversarial_scrambled.json",    # Figure 6, forward corruption
    "v3_adversarial_reverse.json",      # Figure 6, reverse corruption
    "v3_normalisation_inputs.json",      # R1.1 input normalisation, frozen inputs
    "v3_input_normalisation_main.json",  # arm 1, Claude + GPT-5.2
    "v3_input_normalisation_o3.json",    # arm 1, o3 subsample
    "v3_input_normalisation_defs.json",  # arm 2, the result that changed the conclusion
    "v3_input_normalisation_o3_retry.json",  # the 429 batch: all error, deposited anyway
    "v3_input_normalisation_eval.json",
    "v3_normalisation_disagreement_classes.json",
    "v3_input_normalisation_defs_tail.json",   # arm 2, the final 45 pairs
    "v3_input_normalisation_defs_gpt52.json",  # arm 2, GPT-5.2
    "v3_input_normalisation_defs_o3.json",     # arm 2, o3
    "v3_input_normalisation_defs_deepseek.json",
    "v3_input_normalisation_defs_sonnet_gpt41_o4mini.json",
    "v3_input_normalisation_defs_sonnet_gpt41_o4mini_rest.json",
    "v3_input_normalisation_nodefs_four.json",
    "v3_input_normalisation_seven_model_freeze.json",
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

RELEASE_NOTES = """# Version 1.4.0 release notes

Prepared 1 August 2026 for the Cell Genomics revision titled "Trustworthy
agentic genomics requires validated skills, not better models".

This version supersedes v1.3.1 for the input-normalisation analysis. Correction
C14 makes optional gene prefixes presentation-neutral. Correction C15 accepts an
explicit DIPLOTYPE marker inside Markdown, LaTeX or prose wrapping, applies the
same rule to every model, and uses the last explicit marker as the final answer.
Re-parsing 9,557 stored rows without new API calls recovered 87 marked calls and
removed two provisional calls followed by an explicit final ABSTAIN. No raw
response text was changed.

The definition-supplied arm is frozen across seven models: 3,689 attempts, 2,905
calls, 780 abstentions and four output-budget truncations. The frozen analysis
records source-file hashes, exact paired units and 10,000-replicate
sample-cluster bootstrap intervals. Gemini is assigned no performance estimate
because a complete comparable run was not collected; its boxed pilot response
is accepted by the same model-neutral parser.

Code is pinned by the source archive for tag cg-revision-2026-08-01. See
CORRECTIONS.md, MODEL-VERSIONS.md, DATA-MANIFEST.md and PROVENANCE.md for the
auditable scope and limitations.
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
    ap.add_argument(
        "--tag",
        help="annotated release tag to archive; omit only for a pre-tag dry run",
    )
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

    supplementary = {
        "CORRECTIONS.md": BASE / "CORRECTIONS.md",
        "MODEL-VERSIONS.md": BASE / "MODEL-VERSIONS.md",
        "DATA-MANIFEST.md": DATA / "MANIFEST.md",
        "DATA-README.md": DATA / "README.md",
    }
    for name, src in supplementary.items():
        if not src.exists():
            missing.append(str(src.relative_to(BASE)))
            continue
        dst = args.out / name
        shutil.copy2(src, dst)
        staged.append((name, dst.stat().st_size, sha256(dst)))

    provenance = args.out / "PROVENANCE.md"
    provenance.write_text(PROVENANCE)
    staged.append((provenance.name, provenance.stat().st_size, sha256(provenance)))

    release_notes = args.out / "RELEASE-NOTES-v1.4.0.md"
    release_notes.write_text(RELEASE_NOTES)
    staged.append((release_notes.name, release_notes.stat().st_size,
                   sha256(release_notes)))

    if args.tag:
        archive_name = f"24-AGENTIC-PGX-BENCHMARK_{args.tag}.zip"
        archive = args.out / archive_name
        try:
            subprocess.run(
                ["git", "archive", "--format=zip", f"--output={archive}", args.tag],
                cwd=BASE,
                check=True,
            )
        except subprocess.CalledProcessError:
            print(f"Unable to archive tag {args.tag!r}; has it been created?")
            return 1
        staged.append((archive_name, archive.stat().st_size, sha256(archive)))

    lines = ["# SHA-256 checksums, synchronized Cell Genomics revision dataset", ""]
    lines += [f"{h}  {n}" for n, _, h in staged]
    (args.out / "CHECKSUMS.sha256").write_text("\n".join(lines) + "\n")

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
    suffix = " including the tagged source archive" if args.tag else " (pre-tag dry run)"
    print(f"\nAll files present{suffix}. Upload this directory to Zenodo as a new version.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
