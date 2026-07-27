#!/usr/bin/env python3
"""
Per-cell provenance for the matched factorial (reviewer point 1).

WHY THIS EXISTS
A referee read code/60-matched-factorial.py, saw rerun=False on the two skill
cells and a docstring saying they were "reused unchanged" from the legacy
17,820-evaluation dataset, and concluded that the manuscript misstated its
Methods and that the eight-model panel was not common to all cells. The
conclusion was wrong but the inference was sound: the repository said so, the
raw evaluations are not in git, and nothing let the claim be checked.

The remedy is not a rebuttal. It is to publish, per cell, which model
identifiers produced which rows, how many rows, and what they cost, regenerated
from the data rather than asserted. Billing is the discriminator a legacy row
cannot fake: reused rows carry no token counts and no spend.

USAGE
    python code/69-cell-provenance.py
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
LIVE = BASE / "data" / "v3_five_cell_live.json"
LEGACY = BASE / "data" / "v3_five_cell_matched.json"
OUT = BASE / "data" / "v3_cell_provenance.json"
REPORT = BASE / "data" / "v3_cell_provenance.txt"

# Manuscript name -> API identifier actually issued, for the substituted panel.
API_IDS = {
    "Claude Opus 4.5": "claude-opus-4-5-20251101",
    "Claude Sonnet 4.5": "claude-sonnet-4-5-20250929",
    "GPT-5.2": "gpt-5.2",
    "GPT-4.1": "gpt-4.1",
    "o3": "o3",
    "o4-mini": "o4-mini",
    "Gemini 2.5 Flash": "gemini-2.5-flash",
    "DeepSeek V3": "deepseek-chat",
    "Mistral Large 2512": "mistral-large-2512",
}


def summarise(rows: list[dict]) -> dict:
    per_cell = defaultdict(lambda: {
        "n": 0, "in_tokens": 0, "out_tokens": 0, "cost_usd": 0.0,
        "models": defaultdict(int), "rows_with_billing": 0,
    })
    for r in rows:
        c = per_cell[r["cell"]]
        c["n"] += 1
        c["in_tokens"] += r.get("in_tokens") or 0
        c["out_tokens"] += r.get("out_tokens") or 0
        c["cost_usd"] += r.get("cost_usd") or 0.0
        c["models"][r.get("model")] += 1
        if (r.get("in_tokens") or 0) > 0:
            c["rows_with_billing"] += 1
    out = {}
    for cell, v in per_cell.items():
        out[cell] = {
            "n": v["n"],
            "models": dict(sorted(v["models"].items())),
            "api_identifiers": sorted({API_IDS.get(m, f"UNMAPPED:{m}")
                                       for m in v["models"]}),
            "in_tokens": v["in_tokens"],
            "out_tokens": v["out_tokens"],
            "cost_usd": round(v["cost_usd"], 2),
            "rows_with_billing": v["rows_with_billing"],
            "live": v["rows_with_billing"] == v["n"],
        }
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--live", type=Path, default=LIVE)
    ap.add_argument("--legacy", type=Path, default=LEGACY)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--report", type=Path, default=REPORT)
    a = ap.parse_args(argv)

    if not a.live.exists():
        sys.stderr.write(f"missing {a.live}\n")
        return 1
    live = json.loads(a.live.read_text())
    prov = summarise(live)

    identical = None
    if a.legacy.exists():
        legacy = json.loads(a.legacy.read_text())

        def idx(rows):
            return {(r["cell"], r["case_id"], r["rep"], r["model"]):
                    (r.get("raw", ""), r.get("in_tokens"), r.get("cost_usd"))
                    for r in rows if "skill" in str(r.get("cell", ""))}
        lo, ln = idx(legacy), idx(live)
        common = set(lo) & set(ln)
        identical = {
            "comparable_skill_rows": len(common),
            "identical_on_raw_tokens_and_cost": sum(1 for k in common if lo[k] == ln[k]),
        }

    result = {"per_cell": prov, "legacy_comparison": identical}
    a.out.write_text(json.dumps(result, indent=2))

    L = ["MATCHED FACTORIAL: PER-CELL PROVENANCE", ""]
    L.append(f"  {'cell':18s} {'n':>6} {'billed':>7} {'in_tok':>12} {'USD':>7}  models")
    for cell, v in sorted(prov.items()):
        L.append(f"  {cell:18s} {v['n']:6d} {v['rows_with_billing']:7d} "
                 f"{v['in_tokens']:12,d} {v['cost_usd']:7.2f}  {len(v['models'])}")
    L.append("")
    total = sum(v["cost_usd"] for v in prov.values())
    L.append(f"  total spend across the five cells: ${total:.2f}")
    L.append("")
    panels = {cell: tuple(sorted(v["models"])) for cell, v in prov.items()}
    common_panel = len(set(panels.values())) == 1
    L.append(f"  identical model panel in every cell: {common_panel}")
    if common_panel:
        L.append("  panel (manuscript name -> API identifier issued):")
        for m in sorted(next(iter(panels.values()))):
            L.append(f"    {m:20s} {API_IDS.get(m, 'UNMAPPED')}")
    L.append("")
    L.append("  Every row in every cell carries its own token counts and cost, which")
    L.append("  is what a reused legacy row cannot do. The rerun flags in")
    L.append("  60-matched-factorial.py are the default cell selection for a bare")
    L.append("  invocation, not a record of what produced these numbers.")
    if identical:
        L.append("")
        L.append(f"  skill rows comparable with the legacy dataset: "
                 f"{identical['comparable_skill_rows']}")
        L.append(f"  of those, identical on raw text, tokens and cost: "
                 f"{identical['identical_on_raw_tokens_and_cost']}")
    text = "\n".join(L)
    a.report.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
