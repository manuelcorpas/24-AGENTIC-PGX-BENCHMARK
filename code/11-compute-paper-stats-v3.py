#!/usr/bin/env python3
"""
Compute every quantitative claim that goes into the BiB manuscript v3, from
the rescored v3 raw data. Single source of truth for numbers in abstract,
key points, results, and tables.

Reads:
  ../RESULTS/v3_raw_rescored.json  (output of 10-rescore-v3.py)
  ../SPECS/test_cases_v3.json

Writes:
  ../RESULTS/v3_paper_stats.txt    (human-readable)
  ../RESULTS/v3_paper_stats.json   (machine-readable, keyed by claim ID)

Claim IDs are referenced from the manuscript so reviewers can map any number
back to its computation.
"""
from __future__ import annotations
import json
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "RESULTS" / "v3_raw_rescored.json"
CASES = BASE / "SPECS" / "test_cases_v3.json"
OUT_TXT = BASE / "RESULTS" / "v3_paper_stats.txt"
OUT_JSON = BASE / "RESULTS" / "v3_paper_stats.json"


def parsed_only(rows):
    return [r for r in rows if not r["scores"].get("format_fail")]


def per_combo_mean(rows, key, score):
    combos = defaultdict(list)
    for r in rows:
        combos[key(r)].append(r)
    if not combos:
        return 0.0
    means = [sum(r["scores"][score] for r in v) / len(v) for v in combos.values()]
    return sum(means) / len(means)


def main():
    rows = json.loads(RAW.read_text())
    cases = json.loads(CASES.read_text())
    cbi = {c["id"]: c for c in cases}

    parsed = parsed_only(rows)
    no_spec = [r for r in parsed if r["cond"] == "no_spec"]
    with_spec = [r for r in parsed if r["cond"] == "with_spec"]
    models = sorted({r["model"] for r in rows})
    pops = sorted({r["pop"] for r in rows})

    claims = {}

    # -------- OVERVIEW --------
    claims["overview.total_evaluations"] = len(rows)
    claims["overview.total_parsed"] = len(parsed)
    claims["overview.parse_rate_pct"] = round(100 * len(parsed) / len(rows), 2)
    claims["overview.api_errors"] = sum(1 for r in rows if "error" in r)
    claims["overview.n_models"] = len(models)
    claims["overview.n_cases"] = len({r["tc"] for r in rows})
    claims["overview.n_populations"] = len(pops)
    claims["overview.n_runs"] = max(r["run"] for r in rows) + 1
    claims["overview.n_genes"] = len({c["gene"] for c in cases})
    claims["overview.n_gene_drug_pairs"] = len({(c["gene"], c["drug"]) for c in cases})

    # Per-model parse rate
    parse_rates = {}
    for m in models:
        n = sum(1 for r in rows if r["model"] == m)
        p = sum(1 for r in parsed if r["model"] == m)
        parse_rates[m] = {"parsed": p, "total": n, "pct": round(100 * p / n, 2)}
    claims["overview.parse_rate_per_model"] = parse_rates

    rates = [v["pct"] for v in parse_rates.values()]
    claims["overview.parse_rate_min"] = min(rates)
    claims["overview.parse_rate_max"] = max(rates)

    # -------- TIER A AGGREGATE --------
    for cond_label, rs in [("no_spec", no_spec), ("with_spec", with_spec)]:
        for dim in ["A1", "A2", "A3"]:
            correct = sum(1 for r in rs if r["scores"].get(dim) == 1.0)
            partial = sum(1 for r in rs if r["scores"].get(dim) == 0.5)
            claims[f"tier_a.{cond_label}.{dim}.correct"] = correct
            claims[f"tier_a.{cond_label}.{dim}.partial"] = partial
            claims[f"tier_a.{cond_label}.{dim}.n"] = len(rs)
            claims[f"tier_a.{cond_label}.{dim}.pct_correct"] = round(100 * correct / len(rs), 2)
            mean = sum(r["scores"].get(dim, 0) for r in rs) / len(rs)
            claims[f"tier_a.{cond_label}.{dim}.mean"] = round(mean, 4)

    # -------- TIER A: per-model x condition (Table 1) --------
    table1 = {}
    for m in models:
        table1[m] = {}
        for dim in ["A1", "A2", "A3"]:
            for cond in ["no_spec", "with_spec"]:
                rs = [r for r in parsed if r["model"] == m and r["cond"] == cond]
                v = per_combo_mean(rs, lambda r: (r["tc"], r["pop"]), dim)
                table1[m][f"{dim}_{cond}"] = round(v, 2)
    claims["table1_tier_a_by_model"] = table1

    # -------- TIER B: per-model x condition (Table 2) --------
    table2 = {}
    for m in models:
        table2[m] = {}
        for dim in ["B1", "B2", "B3"]:
            for cond in ["no_spec", "with_spec"]:
                rs = [r for r in parsed if r["model"] == m and r["cond"] == cond]
                v = per_combo_mean(rs, lambda r: (r["tc"], r["pop"]), dim)
                table2[m][f"{dim}_{cond}"] = round(v, 2)
    claims["table2_tier_b_by_model"] = table2

    # -------- CONSISTENCY (perfect 3/3) --------
    for cond_label, rs in [("no_spec", no_spec), ("with_spec", with_spec)]:
        combos = defaultdict(list)
        for r in rs:
            combos[(r["model"], r["tc"], r["pop"])].append(r)
        eligible = {k: v for k, v in combos.items() if len(v) >= 2}
        perfect = sum(1 for v in eligible.values() if all(x["scores"]["A1"] == 1.0 for x in v))
        claims[f"consistency.{cond_label}.numerator"] = perfect
        claims[f"consistency.{cond_label}.denominator"] = len(eligible)
        claims[f"consistency.{cond_label}.pct"] = round(100 * perfect / len(eligible), 2) if eligible else 0

    # Stochastic failure rate (no_spec only)
    combos_ns = defaultdict(list)
    for r in no_spec:
        combos_ns[(r["model"], r["tc"], r["pop"])].append(r)
    eligible_ns = {k: v for k, v in combos_ns.items() if len(v) >= 2}
    imperfect = [v for v in eligible_ns.values() if not all(x["scores"]["A1"] == 1.0 for x in v)]
    truly_stochastic = sum(1 for v in imperfect if any(x["scores"]["A1"] == 1.0 for x in v))
    consistent_fail = sum(1 for v in imperfect if not any(x["scores"]["A1"] == 1.0 for x in v))
    claims["consistency.no_spec.imperfect_total"] = len(imperfect)
    claims["consistency.no_spec.truly_stochastic"] = truly_stochastic
    claims["consistency.no_spec.consistent_fail"] = consistent_fail
    claims["consistency.no_spec.imperfect_pct"] = round(100 * len(imperfect) / len(eligible_ns), 2)

    # -------- POPULATION EFFECT (Tier A by population, no_spec) --------
    for pop in pops:
        rs = [r for r in no_spec if r["pop"] == pop]
        v = per_combo_mean(rs, lambda r: (r["model"], r["tc"]), "A1") * 100
        claims[f"population.no_spec.{pop}.A1_pct"] = round(v, 2)

    # -------- DPYD LETHAL CASE --------
    dpyd_ns = [r for r in no_spec if r["tc"] == "dpyd_fu_pm"]
    claims["dpyd.no_spec.parsed_total"] = len(dpyd_ns)
    claims["dpyd.no_spec.A3_errors_total"] = sum(1 for r in dpyd_ns if r["scores"].get("A3", 1.0) < 1.0)
    for pop in pops:
        ps = [r for r in dpyd_ns if r["pop"] == pop]
        pe = [r for r in ps if r["scores"].get("A3", 1.0) < 1.0]
        claims[f"dpyd.no_spec.{pop}.A3_errors"] = len(pe)
        claims[f"dpyd.no_spec.{pop}.parsed"] = len(ps)
    dpyd_ws = [r for r in with_spec if r["tc"] == "dpyd_fu_pm"]
    claims["dpyd.with_spec.parsed_total"] = len(dpyd_ws)
    claims["dpyd.with_spec.A3_errors_total"] = sum(1 for r in dpyd_ws if r["scores"].get("A3", 1.0) < 1.0)

    # -------- ALL LETHAL-CASE A3 ERRORS BY GENE --------
    lethal_by_gene = defaultdict(int)
    for r in no_spec:
        c = cbi[r["tc"]]
        if "lethal" not in c["gt_drug"].lower(): continue
        if r["scores"].get("A3", 1.0) < 1.0:
            lethal_by_gene[c["gene"]] += 1
    claims["lethal_errors.no_spec.by_gene"] = dict(lethal_by_gene)
    claims["lethal_errors.no_spec.total"] = sum(lethal_by_gene.values())
    claims["lethal_errors.with_spec.total"] = sum(
        1 for r in with_spec if "lethal" in cbi[r["tc"]]["gt_drug"].lower() and r["scores"].get("A3", 1.0) < 1.0
    )

    # Non-European share of lethal errors
    lethal_eur = sum(1 for r in no_spec if "lethal" in cbi[r["tc"]]["gt_drug"].lower()
                     and r["scores"].get("A3", 1.0) < 1.0 and r["pop"] == "EUR")
    lethal_non_eur = sum(1 for r in no_spec if "lethal" in cbi[r["tc"]]["gt_drug"].lower()
                          and r["scores"].get("A3", 1.0) < 1.0 and r["pop"] != "EUR")
    claims["lethal_errors.no_spec.EUR"] = lethal_eur
    claims["lethal_errors.no_spec.non_EUR"] = lethal_non_eur

    # -------- TIER B SUMMARY --------
    for cond_label, rs in [("no_spec", no_spec), ("with_spec", with_spec)]:
        for dim in ["B1", "B2", "B3"]:
            mean = sum(r["scores"].get(dim, 0) for r in rs) / len(rs)
            claims[f"tier_b.{cond_label}.{dim}.mean"] = round(mean, 3)

    # -------- COST + RUNTIME (from cost log) --------
    cost_path = BASE / "RESULTS" / "v3_cost_log.json"
    if cost_path.exists():
        cost_log = json.loads(cost_path.read_text())
        total_cost = sum(v["cost"] for v in cost_log.values())
        claims["benchmark.total_cost_usd"] = round(total_cost, 2)
        claims["benchmark.cost_by_model"] = {k: round(v["cost"], 2) for k, v in cost_log.items()}

    # -------- HUMAN-READABLE REPORT --------
    lines = ["# BiB v3 paper statistics", f"Source: {RAW.name}", ""]

    lines.append("## Overview")
    lines.append(f"  Total evaluations: {claims['overview.total_evaluations']}")
    lines.append(f"  Parsed:            {claims['overview.total_parsed']} ({claims['overview.parse_rate_pct']}%)")
    lines.append(f"  API errors:        {claims['overview.api_errors']}")
    lines.append(f"  Models:            {claims['overview.n_models']}")
    lines.append(f"  Cases:             {claims['overview.n_cases']} across {claims['overview.n_genes']} genes / {claims['overview.n_gene_drug_pairs']} gene-drug pairs")
    lines.append(f"  Populations:       {claims['overview.n_populations']} (EUR/AMR/AFR)")
    lines.append(f"  Runs:              {claims['overview.n_runs']}")
    lines.append(f"  Per-model parse-rate range: {claims['overview.parse_rate_min']}–{claims['overview.parse_rate_max']}%")

    lines.append("")
    lines.append("## Tier A aggregate (parsed only)")
    for cond in ["no_spec", "with_spec"]:
        for dim in ["A1", "A2", "A3"]:
            c = claims[f"tier_a.{cond}.{dim}.correct"]
            n = claims[f"tier_a.{cond}.{dim}.n"]
            pct = claims[f"tier_a.{cond}.{dim}.pct_correct"]
            lines.append(f"  {cond} {dim}: {c}/{n} = {pct}%")

    lines.append("")
    lines.append("## Table 1 (Tier A by model)")
    lines.append(f"  {'Model':<22} {'A1ns':>6} {'A1ws':>6} {'A2ns':>6} {'A2ws':>6} {'A3ns':>6} {'A3ws':>6}")
    for m in models:
        t = table1[m]
        lines.append(f"  {m:<22} {t['A1_no_spec']:>6.2f} {t['A1_with_spec']:>6.2f} {t['A2_no_spec']:>6.2f} {t['A2_with_spec']:>6.2f} {t['A3_no_spec']:>6.2f} {t['A3_with_spec']:>6.2f}")

    lines.append("")
    lines.append("## Table 2 (Tier B by model)")
    lines.append(f"  {'Model':<22} {'B1ns':>6} {'B1ws':>6} {'B2ns':>6} {'B2ws':>6} {'B3ns':>6} {'B3ws':>6}")
    for m in models:
        t = table2[m]
        lines.append(f"  {m:<22} {t['B1_no_spec']:>6.2f} {t['B1_with_spec']:>6.2f} {t['B2_no_spec']:>6.2f} {t['B2_with_spec']:>6.2f} {t['B3_no_spec']:>6.2f} {t['B3_with_spec']:>6.2f}")

    lines.append("")
    lines.append("## Consistency (perfect 3/3)")
    for cond in ["no_spec", "with_spec"]:
        n = claims[f"consistency.{cond}.numerator"]
        d = claims[f"consistency.{cond}.denominator"]
        p = claims[f"consistency.{cond}.pct"]
        lines.append(f"  {cond}: {n}/{d} ({p}%)")
    lines.append(f"  no_spec stochastic-failure decomposition: imperfect={claims['consistency.no_spec.imperfect_total']}, "
                 f"truly stochastic={claims['consistency.no_spec.truly_stochastic']}, "
                 f"consistent failure={claims['consistency.no_spec.consistent_fail']}")

    lines.append("")
    lines.append("## Population effect (no_spec A1 by pop)")
    for pop in pops:
        lines.append(f"  {pop}: {claims[f'population.no_spec.{pop}.A1_pct']}%")

    lines.append("")
    lines.append("## DPYD lethal case (rs3918290 T/T)")
    lines.append(f"  no_spec: {claims['dpyd.no_spec.A3_errors_total']}/{claims['dpyd.no_spec.parsed_total']} A3 errors")
    for pop in pops:
        lines.append(f"    {pop}: {claims[f'dpyd.no_spec.{pop}.A3_errors']}/{claims[f'dpyd.no_spec.{pop}.parsed']}")
    lines.append(f"  with_spec: {claims['dpyd.with_spec.A3_errors_total']}/{claims['dpyd.with_spec.parsed_total']}")

    lines.append("")
    lines.append("## All lethal-case A3 errors by gene (no_spec)")
    for gene, n in sorted(claims["lethal_errors.no_spec.by_gene"].items(), key=lambda x: -x[1]):
        lines.append(f"  {gene:<20} {n}")
    lines.append(f"  TOTAL no_spec lethal errors: {claims['lethal_errors.no_spec.total']} (EUR={claims['lethal_errors.no_spec.EUR']}, non-EUR={claims['lethal_errors.no_spec.non_EUR']})")
    lines.append(f"  TOTAL with_spec lethal errors: {claims['lethal_errors.with_spec.total']}")

    lines.append("")
    lines.append("## Tier B drift")
    for cond in ["no_spec", "with_spec"]:
        for dim in ["B1", "B2", "B3"]:
            v = claims[f"tier_b.{cond}.{dim}.mean"]
            lines.append(f"  {cond} {dim}: {v}")

    if "benchmark.total_cost_usd" in claims:
        lines.append("")
        lines.append("## Benchmark cost")
        lines.append(f"  Total: ${claims['benchmark.total_cost_usd']}")
        for k, v in sorted(claims["benchmark.cost_by_model"].items()):
            lines.append(f"    {k:<22} ${v}")

    text = "\n".join(lines)
    OUT_TXT.write_text(text)
    OUT_JSON.write_text(json.dumps(claims, indent=2, default=str))
    print(text)
    print(f"\nWrote {OUT_TXT.name}, {OUT_JSON.name}")


if __name__ == "__main__":
    main()
