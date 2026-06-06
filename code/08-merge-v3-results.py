#!/usr/bin/env python3
"""
Merge v3_raw_mistral.json + v3_raw_others.json -> v3_raw.json.

Used when the v3 benchmark was run as two parallel processes (Mistral on a
longer rate-limit interval and the other 8 models on tight spacing).

Reads:
  ../RESULTS/v3_raw_mistral.json
  ../RESULTS/v3_raw_others.json
Writes:
  ../RESULTS/v3_raw.json          — combined list, sorted (run, model, tc, pop, cond)
  ../RESULTS/v3_cost_log.json     — combined cost-by-model
  ../RESULTS/v3_summary.txt       — quick QC report
"""
from __future__ import annotations
import json
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parent.parent
R = BASE / "RESULTS"

MISTRAL = R / "v3_raw_mistral.json"
OTHERS = R / "v3_raw_others.json"
MISTRAL_COST = R / "v3_cost_log_mistral.json"
OTHERS_COST = R / "v3_cost_log_others.json"

OUT_RAW = R / "v3_raw.json"
OUT_COST = R / "v3_cost_log.json"
OUT_SUMMARY = R / "v3_summary.txt"


def main():
    if not MISTRAL.exists():
        raise SystemExit(f"Missing {MISTRAL}; run Mistral process first")
    if not OTHERS.exists():
        raise SystemExit(f"Missing {OTHERS}; run Others process first")

    m = json.loads(MISTRAL.read_text())
    o = json.loads(OTHERS.read_text())
    combined = m + o
    combined.sort(key=lambda r: (r["run"], r["model"], r["tc"], r["pop"], r["cond"]))
    OUT_RAW.write_text(json.dumps(combined, indent=2))

    cost = {}
    for path in (MISTRAL_COST, OTHERS_COST):
        if path.exists():
            data = json.loads(path.read_text())
            for k, v in data.items():
                cost[k] = v
    OUT_COST.write_text(json.dumps(cost, indent=2, default=str))

    # QC report
    lines = [
        f"# v3 benchmark merged",
        f"Total rows: {len(combined)} (mistral={len(m)} + others={len(o)})",
        "",
        "## Per-model summary",
    ]
    by_model = defaultdict(lambda: {"calls": 0, "errors": 0, "parse_fail": 0, "a1_correct": 0})
    for r in combined:
        m_ = by_model[r["model"]]
        m_["calls"] += 1
        if "error" in r:
            m_["errors"] += 1
        if r["scores"].get("format_fail"):
            m_["parse_fail"] += 1
        if r["scores"].get("A1") == 1.0:
            m_["a1_correct"] += 1

    for model, s in sorted(by_model.items()):
        parsed = s["calls"] - s["parse_fail"]
        a1_pct = 100 * s["a1_correct"] / parsed if parsed else 0
        lines.append(
            f"  {model:<22} calls={s['calls']:>5} errors={s['errors']:>3} "
            f"parse_fail={s['parse_fail']:>3} A1_correct={s['a1_correct']}/{parsed} ({a1_pct:.1f}%)"
        )

    lines.append("")
    lines.append("## Cost summary")
    total_cost = 0
    for k, v in sorted(cost.items()):
        c = v.get("cost", 0)
        total_cost += c
        lines.append(f"  {k:<22} ${c:.2f}  (in={v.get('in_tok',0):>9}, out={v.get('out_tok',0):>9}, calls={v.get('calls',0):>5})")
    lines.append(f"  {'TOTAL':<22} ${total_cost:.2f}")

    OUT_SUMMARY.write_text("\n".join(lines))
    print("\n".join(lines))
    print()
    print(f"Wrote {OUT_RAW.name}")
    print(f"Wrote {OUT_COST.name}")
    print(f"Wrote {OUT_SUMMARY.name}")


if __name__ == "__main__":
    main()
