#!/usr/bin/env python3
"""
Recompute every quantitative claim in the BiB manuscript body from v2_raw.json.

Each metric is computed under a documented filter, compared to the claimed
manuscript value, and reported as PASS / MISMATCH.

Reads:  ../RESULTS/v2_raw.json
Writes: ../RESULTS/paper_stats.txt   (human-readable report)
        ../RESULTS/paper_stats.json  (machine-readable values)

Run from anywhere:  python3 05-compute-paper-stats.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Callable

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "RESULTS" / "v2_raw.json"
OUT_TXT = BASE / "RESULTS" / "paper_stats.txt"
OUT_JSON = BASE / "RESULTS" / "paper_stats.json"

# ----------------------------------------------------------------------------
# Filter primitives
# ----------------------------------------------------------------------------

def is_parsed(r: dict) -> bool:
    return bool(r.get("parsed"))


def by_combo(rows: list, key: Callable[[dict], tuple]) -> dict[tuple, list]:
    out: dict[tuple, list] = defaultdict(list)
    for r in rows:
        out[key(r)].append(r)
    return out


def per_combo_mean(rows: list, key: Callable[[dict], tuple], score: str) -> float:
    """V3 method: mean of per-combo means. This is the manuscript's A1 aggregation."""
    combos = by_combo(rows, key)
    means = [sum(r["scores"][score] for r in v) / len(v) for v in combos.values()]
    return sum(means) / len(means) if means else 0.0


# ----------------------------------------------------------------------------
# Result record
# ----------------------------------------------------------------------------

@dataclass
class Stat:
    paragraph: int
    name: str
    filter_def: str
    computed: Any
    claim: Any = None
    match: str = ""

    def status(self) -> str:
        if self.claim is None:
            return "INFO"
        # Tolerance for floats
        if isinstance(self.computed, float) and isinstance(self.claim, (int, float)):
            return "PASS" if abs(self.computed - self.claim) < 0.005 else "MISMATCH"
        if self.computed == self.claim:
            return "PASS"
        # Allow string equality on rendered percentage
        return "MISMATCH"

    def render_line(self) -> str:
        s = self.status()
        comp = self.computed
        if isinstance(comp, float):
            comp = f"{comp:.4f}"
        claim = self.claim
        if claim is None:
            claim = "(no claim)"
        return f"[{s:8s}] §{self.paragraph:3d} {self.name}\n           filter: {self.filter_def}\n           computed: {comp}    claim: {claim}\n"


# ----------------------------------------------------------------------------
# Metric computations
# ----------------------------------------------------------------------------

def compute_all(data: list) -> list[Stat]:
    parsed = [r for r in data if is_parsed(r)]
    no_spec = [r for r in parsed if r["cond"] == "no_spec"]
    with_spec = [r for r in parsed if r["cond"] == "with_spec"]
    out: list[Stat] = []

    # ---- §29 Overview ----
    out.append(Stat(29, "total evaluations", "len(data)", len(data), 1944))
    out.append(Stat(29, "parseable evaluations", "parsed=true", len(parsed), 1681))
    out.append(Stat(29, "parse rate %", "len(parsed)/len(data)*100", round(len(parsed) / len(data) * 100, 1), 86.5))

    api_errors = sum(1 for r in data if r.get("api_error"))
    out.append(Stat(29, "API errors", "row.api_error truthy", api_errors, 1))

    mistral_total = sum(1 for r in data if r["model"] == "Mistral Large 2")
    mistral_parsed = sum(1 for r in parsed if r["model"] == "Mistral Large 2")
    out.append(Stat(29, "Mistral Large 2 total runs", "model='Mistral Large 2'", mistral_total, 216))
    out.append(Stat(29, "Mistral Large 2 parsed runs", "parsed and model=Mistral", mistral_parsed, 12))
    mistral_pct = round(mistral_parsed / mistral_total * 100, 1) if mistral_total else 0
    out.append(Stat(29, "Mistral parse rate %", "parsed/total *100", mistral_pct, 5.6))

    # Other models 95-100% parse rates
    nonm_rates = sorted({
        m: round(sum(1 for r in parsed if r["model"] == m) /
                 sum(1 for r in data if r["model"] == m) * 100, 1)
        for m in {r["model"] for r in data} if m != "Mistral Large 2"
    }.items(), key=lambda kv: kv[1])
    nonm_min = min(v for _, v in nonm_rates)
    nonm_max = max(v for _, v in nonm_rates)
    out.append(Stat(29, "non-Mistral parse-rate range",
                    "min..max across non-Mistral models",
                    f"{nonm_min}..{nonm_max}",
                    "85..100"))

    # ---- §30 Tier A overview ----
    a1_ns_v3 = per_combo_mean(no_spec, lambda r: (r["model"], r["tc"], r["pop"]), "A1") * 100
    a1_ws_v3 = per_combo_mean(with_spec, lambda r: (r["model"], r["tc"], r["pop"]), "A1") * 100
    out.append(Stat(30, "mean A1 no_spec (%)", "per-(model,tc,pop) mean, then mean of combo means", round(a1_ns_v3, 1), 92.4))
    out.append(Stat(30, "mean A1 with_spec (%)", "per-(model,tc,pop) mean, then mean of combo means", round(a1_ws_v3, 1), 100.0))

    # Perfect consistency under manuscript filter (>=2 runs parsed, all parsed correct)
    combos_ns = by_combo(no_spec, lambda r: (r["model"], r["tc"], r["pop"]))
    eligible_ns = {k: v for k, v in combos_ns.items() if len(v) >= 2}
    perfect_ns = sum(1 for v in eligible_ns.values() if all(r["scores"]["A1"] >= 1.0 for r in v))
    out.append(Stat(30, "perfect-consistency denominator no_spec",
                    "(model,tc,pop) with >=2 of 3 runs parsed",
                    len(eligible_ns), 275))
    out.append(Stat(30, "perfect-consistency numerator no_spec",
                    "all parsed runs A1>=1.0, denominator above",
                    perfect_ns, 241))
    out.append(Stat(30, "perfect-consistency % no_spec",
                    "numerator/denominator * 100",
                    round(perfect_ns / len(eligible_ns) * 100, 1), 87.6))

    combos_ws = by_combo(with_spec, lambda r: (r["model"], r["tc"], r["pop"]))
    eligible_ws = {k: v for k, v in combos_ws.items() if len(v) >= 2}
    perfect_ws = sum(1 for v in eligible_ws.values() if all(r["scores"]["A1"] >= 1.0 for r in v))
    out.append(Stat(30, "perfect-consistency with_spec",
                    "(model,tc,pop) with >=2 runs parsed, all correct",
                    f"{perfect_ws}/{len(eligible_ws)}",
                    "290/290"))
    out.append(Stat(30, "perfect-consistency % with_spec", "x/y *100",
                    round(perfect_ws / len(eligible_ws) * 100, 1), 100.0))

    # Stochastic-failure decomposition (manuscript treats as one number; we split)
    imperfect_ns = [v for v in eligible_ns.values() if not all(r["scores"]["A1"] >= 1.0 for r in v)]
    truly_stochastic = sum(1 for v in imperfect_ns
                           if any(r["scores"]["A1"] >= 1.0 for r in v))
    consistent_fail = sum(1 for v in imperfect_ns
                          if not any(r["scores"]["A1"] >= 1.0 for r in v))
    total_imperfect = len(imperfect_ns)
    out.append(Stat(30, "imperfect total no_spec", "denom - perfect",
                    total_imperfect,
                    "12.4% (manuscript labels all as 'stochastic')"))
    out.append(Stat(30, "  truly stochastic (runs disagree)", "imperfect with at least one correct run",
                    f"{truly_stochastic}/{len(eligible_ns)} = {100*truly_stochastic/len(eligible_ns):.1f}%",
                    "(not separately reported)"))
    out.append(Stat(30, "  consistent failure (all parsed wrong)", "imperfect with no correct run",
                    f"{consistent_fail}/{len(eligible_ns)} = {100*consistent_fail/len(eligible_ns):.1f}%",
                    "(not separately reported)"))

    # ---- §32 Tier A by dimension ----
    for dim in ["A1", "A2", "A3"]:
        ns = per_combo_mean(no_spec, lambda r: (r["model"], r["tc"], r["pop"]), dim)
        ws = per_combo_mean(with_spec, lambda r: (r["model"], r["tc"], r["pop"]), dim)
        out.append(Stat(32, f"{dim} cross-model mean no_spec", "V3 method", round(ns, 2), None))
        out.append(Stat(32, f"{dim} cross-model mean with_spec", "V3 method", round(ws, 2), None))
    # Manuscript-specific A2 claim: 0.87 -> 0.97
    a2_ns = per_combo_mean(no_spec, lambda r: (r["model"], r["tc"], r["pop"]), "A2")
    a2_ws = per_combo_mean(with_spec, lambda r: (r["model"], r["tc"], r["pop"]), "A2")
    out.append(Stat(32, "A2 no_spec (claim: 0.87)", "V3 method", round(a2_ns, 2), 0.87))
    out.append(Stat(32, "A2 with_spec (claim: 0.96)", "V3 method", round(a2_ws, 2), 0.96))

    # ---- §35 Gemini & A3 ----
    gemini_ns = [r for r in no_spec if r["model"] == "Gemini 2.5 Flash"]
    gemini_a1_mean = sum(r["scores"]["A1"] for r in gemini_ns) / len(gemini_ns) * 100
    out.append(Stat(35, "Gemini 2.5 Flash A1 % no_spec",
                    "simple mean A1 over parsed",
                    round(gemini_a1_mean, 0), 61))

    # A3 cross-model mean no_spec
    a3_means = {}
    for m in {r["model"] for r in no_spec}:
        rows_m = [r for r in no_spec if r["model"] == m]
        a3_means[m] = sum(r["scores"]["A3"] for r in rows_m) / len(rows_m)
    overall_a3 = sum(a3_means.values()) / len(a3_means)
    out.append(Stat(35, "A3 mean across models no_spec", "mean over per-model A3 means", round(overall_a3, 2), 0.99))
    out.append(Stat(35, "Gemini A3 no_spec", "mean A3 over Gemini parsed rows",
                    round(a3_means.get("Gemini 2.5 Flash", -1), 2), 0.95))

    # ---- §37 DPYD lethal case ----
    dpyd_ns = [r for r in no_spec if r["tc"] == "dpyd_hom"]
    dpyd_ws = [r for r in with_spec if r["tc"] == "dpyd_hom"]
    out.append(Stat(37, "DPYD lethal: total parsed no_spec", "tc='dpyd_hom' parsed", len(dpyd_ns), None))
    for pop, claim_n, claim_d in [("EUR", 1, 23), ("AMR", 3, 20), ("AFR", 3, 24)]:
        rows = [r for r in dpyd_ns if r["pop"] == pop]
        errs = [r for r in rows if r["scores"]["A1"] < 1.0]
        out.append(Stat(37, f"DPYD lethal {pop} errors", f"pop={pop}, A1<1.0",
                        f"{len(errs)}/{len(rows)}", f"{claim_n}/{claim_d}"))
    total_err_ns = sum(1 for r in dpyd_ns if r["scores"]["A1"] < 1.0)
    out.append(Stat(37, "DPYD lethal total errors no_spec", "A1<1.0", total_err_ns, 7))
    non_eur_err = sum(1 for r in dpyd_ns if r["scores"]["A1"] < 1.0 and r["pop"] != "EUR")
    out.append(Stat(37, "DPYD lethal errors in non-European",
                    "pop in {AMR,AFR}, A1<1.0", non_eur_err, 6))
    total_ws = len([r for r in with_spec if r["tc"] == "dpyd_hom"])
    err_ws = sum(1 for r in dpyd_ws if r["scores"]["A1"] < 1.0)
    out.append(Stat(37, "DPYD lethal errors with_spec",
                    f"tc='dpyd_hom' with_spec parsed (n={total_ws})",
                    f"{err_ws}/{total_ws}", "0/72"))

    # ---- §39 Consistency analysis ----
    # Already covered above for 241/275; plus Gemini perfect-consistency rate
    gem_combos = by_combo([r for r in no_spec if r["model"] == "Gemini 2.5 Flash"],
                          lambda r: (r["tc"], r["pop"]))
    gem_eligible = {k: v for k, v in gem_combos.items() if len(v) >= 2}
    gem_perfect = sum(1 for v in gem_eligible.values() if all(r["scores"]["A1"] >= 1.0 for r in v))
    gem_pct = round(100 * gem_perfect / len(gem_eligible), 0)
    out.append(Stat(39, "Gemini perfect-consistency rate %",
                    "(tc,pop) with >=2 parsed, all correct, /eligible",
                    int(gem_pct), 36))

    # ---- §42 Tier B ----
    for dim, claim_ns, claim_ws in [("B1", 0.69, 0.38), ("B2", 0.99, 1.00), ("B3", 0.63, 0.09)]:
        b_ns = per_combo_mean(no_spec, lambda r: (r["model"], r["tc"], r["pop"]), dim)
        b_ws = per_combo_mean(with_spec, lambda r: (r["model"], r["tc"], r["pop"]), dim)
        out.append(Stat(42, f"{dim} mean no_spec", "V3 method", round(b_ns, 2), claim_ns))
        out.append(Stat(42, f"{dim} mean with_spec", "V3 method", round(b_ws, 2), claim_ws))

    # ---- §47 Population effect ----
    for pop, claim in [("EUR", 93), ("AMR", 92), ("AFR", 92)]:
        rows = [r for r in no_spec if r["pop"] == pop]
        # V3 method per population
        combos = by_combo(rows, lambda r: (r["model"], r["tc"]))
        means = [sum(r["scores"]["A1"] for r in v) / len(v) for v in combos.values()]
        pct = round(sum(means) / len(means) * 100, 0)
        out.append(Stat(47, f"A1 % no_spec {pop}", "V3 method per population", pct, claim))

    # Gemini European vs East African A1 (claim 67% -> 53%)
    for pop, claim in [("EUR", 67), ("AFR", 53)]:
        rows = [r for r in no_spec if r["model"] == "Gemini 2.5 Flash" and r["pop"] == pop]
        if not rows:
            out.append(Stat(47, f"Gemini A1 {pop}", "Gemini, pop, parsed", "no rows", claim))
            continue
        v = sum(r["scores"]["A1"] for r in rows) / len(rows) * 100
        out.append(Stat(47, f"Gemini A1 % no_spec {pop}", "simple mean over parsed", round(v, 0), claim))

    # DeepSeek V3 by population (manuscript example removed; report as INFO only)
    for pop in ["EUR", "AMR", "AFR"]:
        rows = [r for r in no_spec if r["model"] == "DeepSeek V3" and r["pop"] == pop]
        if not rows:
            continue
        v = sum(r["scores"]["A1"] for r in rows) / len(rows) * 100
        out.append(Stat(47, f"DeepSeek V3 A1 % no_spec {pop} (no claim, example removed)",
                        "simple mean over parsed",
                        round(v, 0),
                        None))

    return out


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main() -> None:
    assert RAW.exists(), f"Raw data not found: {RAW}"
    data = json.loads(RAW.read_text())
    stats = compute_all(data)

    lines = ["# Paper statistics report", f"Source: {RAW.name}", f"Total rows: {len(data)}", ""]
    pass_n = mismatch_n = info_n = 0
    for s in stats:
        if s.status() == "PASS":
            pass_n += 1
        elif s.status() == "MISMATCH":
            mismatch_n += 1
        else:
            info_n += 1
        lines.append(s.render_line())
    lines.append("")
    lines.append(f"Summary: PASS={pass_n}, MISMATCH={mismatch_n}, INFO={info_n}, total={len(stats)}")

    OUT_TXT.write_text("\n".join(lines))
    OUT_JSON.write_text(json.dumps([{**asdict(s), "status": s.status()} for s in stats], indent=2, default=str))

    print("\n".join(lines))
    print(f"\nWrote {OUT_TXT.name} and {OUT_JSON.name}")


if __name__ == "__main__":
    main()
