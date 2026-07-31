#!/usr/bin/env python3
"""
Re-parse the stored input-normalisation responses after the C14 parser fix.

WHY
The gene-prefix strip in 74-input-normalisation.py required a space or hyphen
after the gene symbol, so the standard PharmVar rendering "CYP2D6*4/*4" did not
parse and was recorded as an abstention. Coverage is this experiment's headline
metric, so the artefact moved the headline: 744 well-formed diplotypes across
the deposited rows were scored as refusals to answer.

This re-parses the stored `raw` text with the corrected parser and rewrites
`call` and `status`. It issues no API calls and changes no response text, which
is the same remedy applied to C1. Every row's raw text is left byte-identical,
so the correction is auditable by diffing call/status alone.

WHAT IT DOES NOT DO
It does not touch rows that errored at the provider or were truncated by the
output budget. Those are not the model declining to answer and they stay out of
the scored denominator.

USAGE
    python code/78-reparse-normalisation.py            # report only
    python code/78-reparse-normalisation.py --write    # rewrite the data files
"""
from __future__ import annotations

import argparse
import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent

FILES = [
    "v3_input_normalisation.json",
    "v3_input_normalisation_main.json",
    "v3_input_normalisation_o3.json",
    "v3_input_normalisation_o3_retry.json",
    "v3_input_normalisation_defs.json",
    "v3_input_normalisation_defs_o3.json",
    "v3_input_normalisation_defs_gpt52.json",
    "v3_input_normalisation_defs_tail.json",
]


def _load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def reparse(rows: list[dict], norm) -> dict:
    """Rewrite call/status in place. Returns a per-file tally."""
    tally = {"rows": len(rows), "recovered": 0, "lost": 0, "changed_status": 0}
    for r in rows:
        if r.get("error") or r.get("truncated"):
            continue
        before = r.get("call")
        after = norm.parse_call(r.get("raw"), r.get("gene"))
        if before == after:
            continue
        if after and not before:
            tally["recovered"] += 1
        elif before and not after:
            tally["lost"] += 1
        r["call"] = after
        status = norm.classify(r.get("raw"), r.get("out_tokens") or 0,
                               error=r.get("error"), gene=r.get("gene"))
        if status != r.get("status"):
            r["status"] = status
            tally["changed_status"] += 1
    return tally


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true",
                    help="rewrite the data files; otherwise report only")
    args = ap.parse_args(argv)

    norm = _load("norm", CODE / "74-input-normalisation.py")

    total = {"rows": 0, "recovered": 0, "lost": 0, "changed_status": 0}
    for name in FILES:
        path = BASE / "data" / name
        if not path.exists():
            print(f"  {name:44s} absent, skipped")
            continue
        blob = json.loads(path.read_text())
        rows = blob["rows"] if isinstance(blob, dict) else blob
        tally = reparse(rows, norm)
        for k in total:
            total[k] += tally[k]
        print(f"  {name:44s} rows={tally['rows']:5d} "
              f"recovered={tally['recovered']:4d} lost={tally['lost']:3d}")
        if args.write:
            path.write_text(json.dumps(blob, indent=2))

    print(f"\n  {'TOTAL':44s} rows={total['rows']:5d} "
          f"recovered={total['recovered']:4d} lost={total['lost']:3d}")
    if total["lost"]:
        print("  WARNING: the fix removed calls it should only have added. "
              "Inspect before accepting.")
    if not args.write:
        print("\n  report only; re-run with --write to rewrite the data files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
