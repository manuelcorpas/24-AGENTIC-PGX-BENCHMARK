#!/usr/bin/env python3
"""
Recompute every quantitative claim in the BiB manuscript with Mistral Large 2 excluded.

Mirrors 05-compute-paper-stats.py exactly but applies a model filter on load.
Writes paper_stats_no_mistral.{txt,json} for direct comparison.

Run: python3 06-stats-no-mistral.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Callable

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "RESULTS" / "v2_raw.json"
OUT_TXT = BASE / "RESULTS" / "paper_stats_no_mistral.txt"
OUT_JSON = BASE / "RESULTS" / "paper_stats_no_mistral.json"

EXCLUDE_MODELS = {"Mistral Large 2"}


def is_parsed(r: dict) -> bool:
    return bool(r.get("parsed"))


def by_combo(rows, key):
    out = defaultdict(list)
    for r in rows:
        out[key(r)].append(r)
    return out


def per_combo_mean(rows, key, score):
    combos = by_combo(rows, key)
    means = [sum(r["scores"][score] for r in v) / len(v) for v in combos.values()]
    return sum(means) / len(means) if means else 0.0


@dataclass
class Stat:
    section: str
    name: str
    filter_def: str
    computed: Any

    def render(self) -> str:
        comp = self.computed
        if isinstance(comp, float):
            comp = f"{comp:.4f}"
        return f"  {self.section:25s} {self.name:55s} = {comp}"


def compute_all(data: list) -> list[Stat]:
    parsed = [r for r in data if is_parsed(r)]
    no_spec = [r for r in parsed if r["cond"] == "no_spec"]
    with_spec = [r for r in parsed if r["cond"] == "with_spec"]
    out: list[Stat] = []

    # Overview
    out.append(Stat("Overview", "total evaluations", "all rows", len(data)))
    out.append(Stat("Overview", "parseable evaluations", "parsed=true", len(parsed)))
    out.append(Stat("Overview", "parse rate %", "parsed/total*100", round(len(parsed) / len(data) * 100, 1)))
    api_errors = sum(1 for r in data if r.get("api_error"))
    out.append(Stat("Overview", "API errors", "api_error truthy", api_errors))
    out.append(Stat("Overview", "n models", "distinct model values", len({r["model"] for r in data})))

    # Per-model parse rates
    parse_rates = {}
    for m in sorted({r["model"] for r in data}):
        tot = sum(1 for r in data if r["model"] == m)
        prs = sum(1 for r in parsed if r["model"] == m)
        parse_rates[m] = (prs, tot, round(prs / tot * 100, 1))
    for m, (p, t, pct) in parse_rates.items():
        out.append(Stat("Parse-rate", m, "parsed/total*100", f"{p}/{t} = {pct}%"))
    rates = [v[2] for v in parse_rates.values()]
    out.append(Stat("Parse-rate", "min..max across all models", "", f"{min(rates)}..{max(rates)}"))

    # Tier A overview (V3 method)
    a1_ns = per_combo_mean(no_spec, lambda r: (r["model"], r["tc"], r["pop"]), "A1") * 100
    a1_ws = per_combo_mean(with_spec, lambda r: (r["model"], r["tc"], r["pop"]), "A1") * 100
    out.append(Stat("Tier A", "mean A1 no_spec %", "V3 per-(m,tc,pop) mean", round(a1_ns, 1)))
    out.append(Stat("Tier A", "mean A1 with_spec %", "V3 per-(m,tc,pop) mean", round(a1_ws, 1)))

    # Perfect consistency
    combos_ns = by_combo(no_spec, lambda r: (r["model"], r["tc"], r["pop"]))
    eligible_ns = {k: v for k, v in combos_ns.items() if len(v) >= 2}
    perfect_ns = sum(1 for v in eligible_ns.values() if all(r["scores"]["A1"] >= 1.0 for r in v))
    out.append(Stat("Tier A", "perfect-consistency no_spec", ">=2 parsed, all correct",
                    f"{perfect_ns}/{len(eligible_ns)} ({100*perfect_ns/len(eligible_ns):.1f}%)"))

    combos_ws = by_combo(with_spec, lambda r: (r["model"], r["tc"], r["pop"]))
    eligible_ws = {k: v for k, v in combos_ws.items() if len(v) >= 2}
    perfect_ws = sum(1 for v in eligible_ws.values() if all(r["scores"]["A1"] >= 1.0 for r in v))
    out.append(Stat("Tier A", "perfect-consistency with_spec", ">=2 parsed, all correct",
                    f"{perfect_ws}/{len(eligible_ws)} ({100*perfect_ws/len(eligible_ws):.1f}%)"))

    # Stochastic decomposition
    imperfect_ns = [v for v in eligible_ns.values() if not all(r["scores"]["A1"] >= 1.0 for r in v)]
    out.append(Stat("Tier A", "imperfect-combos no_spec", "denom - perfect", len(imperfect_ns)))

    # Tier A by dimension (cross-model means)
    for dim in ["A1", "A2", "A3"]:
        ns = per_combo_mean(no_spec, lambda r: (r["model"], r["tc"], r["pop"]), dim)
        ws = per_combo_mean(with_spec, lambda r: (r["model"], r["tc"], r["pop"]), dim)
        out.append(Stat("Tier A dim", f"{dim} no_spec", "V3", round(ns, 2)))
        out.append(Stat("Tier A dim", f"{dim} with_spec", "V3", round(ws, 2)))

    # Per-model Tier A (for Table 1)
    for m in sorted({r["model"] for r in parsed}):
        for dim in ["A1", "A2", "A3"]:
            for cond, rows in [("no_spec", no_spec), ("with_spec", with_spec)]:
                rs = [r for r in rows if r["model"] == m]
                if not rs:
                    continue
                v = per_combo_mean(rs, lambda r: (r["tc"], r["pop"]), dim)
                out.append(Stat("Table 1", f"{m} {dim} {cond}", "V3 per-(tc,pop)", round(v, 2)))

    # Per-model Tier B (for Table 2)
    for m in sorted({r["model"] for r in parsed}):
        for dim in ["B1", "B2", "B3"]:
            for cond, rows in [("no_spec", no_spec), ("with_spec", with_spec)]:
                rs = [r for r in rows if r["model"] == m]
                if not rs:
                    continue
                v = per_combo_mean(rs, lambda r: (r["tc"], r["pop"]), dim)
                out.append(Stat("Table 2", f"{m} {dim} {cond}", "V3 per-(tc,pop)", round(v, 2)))

    # Tier B aggregate
    for dim in ["B1", "B2", "B3"]:
        ns = per_combo_mean(no_spec, lambda r: (r["model"], r["tc"], r["pop"]), dim)
        ws = per_combo_mean(with_spec, lambda r: (r["model"], r["tc"], r["pop"]), dim)
        out.append(Stat("Tier B agg", f"{dim} no_spec", "V3", round(ns, 2)))
        out.append(Stat("Tier B agg", f"{dim} with_spec", "V3", round(ws, 2)))

    # Gemini specifics
    gem_ns = [r for r in no_spec if r["model"] == "Gemini 2.5 Flash"]
    if gem_ns:
        gem_a1 = sum(r["scores"]["A1"] for r in gem_ns) / len(gem_ns) * 100
        out.append(Stat("Gemini", "A1 no_spec %", "simple mean", round(gem_a1, 0)))
        for pop in ["EUR", "AMR", "AFR"]:
            rs = [r for r in gem_ns if r["pop"] == pop]
            if rs:
                v = sum(r["scores"]["A1"] for r in rs) / len(rs) * 100
                out.append(Stat("Gemini", f"A1 no_spec {pop} %", "simple mean", round(v, 0)))

    # A3 cross-model mean
    a3_means = {}
    for m in {r["model"] for r in no_spec}:
        rs = [r for r in no_spec if r["model"] == m]
        a3_means[m] = sum(r["scores"]["A3"] for r in rs) / len(rs)
    if a3_means:
        out.append(Stat("Tier A", "A3 cross-model mean no_spec", "mean over per-model A3", round(sum(a3_means.values()) / len(a3_means), 2)))
    if "Gemini 2.5 Flash" in a3_means:
        out.append(Stat("Tier A", "Gemini A3 no_spec", "", round(a3_means["Gemini 2.5 Flash"], 2)))

    # DPYD lethal case
    dpyd_ns = [r for r in no_spec if r["tc"] == "dpyd_hom"]
    dpyd_ws = [r for r in with_spec if r["tc"] == "dpyd_hom"]
    out.append(Stat("DPYD", "total parsed no_spec", "tc=dpyd_hom", len(dpyd_ns)))
    for pop in ["EUR", "AMR", "AFR"]:
        rs = [r for r in dpyd_ns if r["pop"] == pop]
        errs = [r for r in rs if r["scores"]["A1"] < 1.0]
        out.append(Stat("DPYD", f"errors {pop} no_spec", "A1<1.0", f"{len(errs)}/{len(rs)}"))
    total_err = sum(1 for r in dpyd_ns if r["scores"]["A1"] < 1.0)
    non_eur = sum(1 for r in dpyd_ns if r["scores"]["A1"] < 1.0 and r["pop"] != "EUR")
    out.append(Stat("DPYD", "total errors no_spec", "A1<1.0", total_err))
    out.append(Stat("DPYD", "non-European errors no_spec", "pop!=EUR, A1<1.0", non_eur))
    out.append(Stat("DPYD", "with_spec errors", "tc=dpyd_hom, A1<1.0",
                    f"{sum(1 for r in dpyd_ws if r['scores']['A1'] < 1.0)}/{len(dpyd_ws)}"))

    # Gemini perfect-consistency rate
    gem_combos = by_combo(gem_ns, lambda r: (r["tc"], r["pop"]))
    gem_eligible = {k: v for k, v in gem_combos.items() if len(v) >= 2}
    if gem_eligible:
        gem_perfect = sum(1 for v in gem_eligible.values() if all(r["scores"]["A1"] >= 1.0 for r in v))
        out.append(Stat("Consistency", "Gemini perfect-consistency rate %",
                        f"{gem_perfect}/{len(gem_eligible)}",
                        round(100 * gem_perfect / len(gem_eligible), 0)))

    # Population effects
    for pop in ["EUR", "AMR", "AFR"]:
        rs = [r for r in no_spec if r["pop"] == pop]
        v = per_combo_mean(rs, lambda r: (r["model"], r["tc"]), "A1") * 100
        out.append(Stat("Population", f"A1 {pop} no_spec %", "V3 per-(m,tc) within pop", round(v, 0)))

    return out


def main():
    data_all = json.loads(RAW.read_text())
    data = [r for r in data_all if r["model"] not in EXCLUDE_MODELS]

    stats = compute_all(data)

    excluded_n = len(data_all) - len(data)
    header = [
        "# Paper statistics report (Mistral excluded)",
        f"Source: {RAW.name}",
        f"Excluded models: {sorted(EXCLUDE_MODELS)}",
        f"Excluded rows: {excluded_n}",
        f"Remaining rows: {len(data)}",
        f"Models remaining: {sorted({r['model'] for r in data})}",
        "",
    ]
    body = [s.render() for s in stats]
    text = "\n".join(header + body)

    OUT_TXT.write_text(text + "\n")
    OUT_JSON.write_text(json.dumps([asdict(s) for s in stats], indent=2, default=str))
    print(text)
    print(f"\nWrote {OUT_TXT.name} and {OUT_JSON.name}")


if __name__ == "__main__":
    main()
