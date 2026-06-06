#!/usr/bin/env python3
"""
Figure 6 for Cell Genomics Manuscript v8.

The population-equity finding: aggregate population gaps under free-prompted and
retrieval-augmented conditions are small but non-zero, and per-locus population
gaps can be substantial (up to ~15 percentage points on individual genes).
Specification-constrained execution eliminates all population variation by
construction: every (case, replicate) cell receives identical output regardless
of the population annotation.

Two-panel design:
  Panel A - aggregate A1 by (condition, population). Grouped bars showing the
            EUR / AMR / AFR means under each condition. With_spec collapses to
            a uniform 100% across all three populations.
  Panel B - per-locus EUR-vs-AFR A1 deltas by condition. Strip plot showing the
            distribution of per-gene EUR-AFR differences. Free-prompted and
            retrieval-augmented conditions show heterogeneous per-locus gaps
            (different genes favour different populations, no monotonic
            EUR > AMR > AFR pattern); with_spec is a delta of 0 on every gene.

Inputs:
  RESULTS/v3_raw_rescored_three_arm.json
  RESULTS/v3_three_arm_analysis_summary.json (cross-check section_1 aggregate)
  SPECS/test_cases_v3.json

Outputs (PNG 300 DPI + TIFF 600 DPI LZW):
  FIGURES/Figure6_population_equity.png
  FIGURES/Figure6_population_equity.tiff

Usage:
  python3 36-figure6-population-equity.py             # render
  python3 36-figure6-population-equity.py --self-test # verify inputs
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


BASE = Path(__file__).resolve().parent.parent
MAIN = BASE / "RESULTS" / "v3_raw_rescored_three_arm.json"
SUMMARY = BASE / "RESULTS" / "v3_three_arm_analysis_summary.json"
FIGDIR = BASE / "FIGURES"

CONDITIONS = ["no_spec", "cpic_rag", "with_spec"]
CONDITION_LABELS = {
    "no_spec": "Free-prompted",
    "cpic_rag": "Retrieval-augmented",
    "with_spec": "Specification-constrained",
}
CONDITION_LABELS_SHORT = {
    "no_spec": "Free-\nprompted",
    "cpic_rag": "Retrieval-\naugmented",
    "with_spec": "Specification-\nconstrained",
}
POPULATIONS = ["EUR", "AMR", "AFR"]
POPULATION_COLOURS = {
    "EUR": "#4a78a8",
    "AMR": "#c47a3a",
    "AFR": "#5a8f5a",
}
CONDITION_COLOURS = {
    "no_spec": "#d4756f",
    "cpic_rag": "#e0a96d",
    "with_spec": "#5a8f5a",
}

EXPECTED = {
    "n_per_cond_per_pop_approx": 2900,
    "with_spec_pop_a1_exact": 1.0,
    "no_spec_eur_a1_approx": 0.810,
    "no_spec_afr_a1_approx": 0.803,
    "cpic_rag_eur_a1_approx": 0.900,
}


def per_pop_means(main: list[dict]) -> dict[str, dict[str, dict]]:
    by: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in main:
        if r["scores"].get("format_fail"):
            continue
        a1 = r["scores"].get("A1")
        if a1 is None:
            continue
        by[r["cond"]][r["pop"]].append(a1)
    out: dict[str, dict[str, dict]] = {}
    for cond, pops in by.items():
        out[cond] = {}
        for pop, vals in pops.items():
            out[cond][pop] = {"mean": sum(vals) / len(vals) if vals else 0.0,
                              "n": len(vals)}
    return out


def per_gene_pop_means(main: list[dict]) -> dict[str, dict[str, dict[str, float]]]:
    by: dict[str, dict[str, dict[str, list[float]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    for r in main:
        if r["scores"].get("format_fail"):
            continue
        a1 = r["scores"].get("A1")
        if a1 is None:
            continue
        by[r["gene"]][r["cond"]][r["pop"]].append(a1)
    out: dict[str, dict[str, dict[str, float]]] = {}
    for gene, conds in by.items():
        out[gene] = {}
        for cond, pops in conds.items():
            out[gene][cond] = {p: sum(v) / len(v) for p, v in pops.items() if v}
    return out


def eur_afr_deltas(per_gene: dict) -> dict[str, list[tuple[str, float]]]:
    out: dict[str, list[tuple[str, float]]] = {c: [] for c in CONDITIONS}
    for gene, conds in per_gene.items():
        for cond in CONDITIONS:
            pops = conds.get(cond, {})
            if "EUR" in pops and "AFR" in pops:
                out[cond].append((gene, pops["EUR"] - pops["AFR"]))
    return out


def self_test() -> int:
    failures: list[str] = []

    main = json.loads(MAIN.read_text())
    pp = per_pop_means(main)

    for cond in CONDITIONS:
        for pop in POPULATIONS:
            n = pp.get(cond, {}).get(pop, {}).get("n", 0)
            if abs(n - EXPECTED["n_per_cond_per_pop_approx"]) > 100:
                failures.append(f"{cond}/{pop} n: got {n}, want ~{EXPECTED['n_per_cond_per_pop_approx']}")

    for pop in POPULATIONS:
        m = pp["with_spec"][pop]["mean"]
        if m != EXPECTED["with_spec_pop_a1_exact"]:
            failures.append(f"with_spec/{pop} A1: got {m}, want exactly 1.0")

    if abs(pp["no_spec"]["EUR"]["mean"] - EXPECTED["no_spec_eur_a1_approx"]) > 0.005:
        failures.append(f"no_spec/EUR A1: got {pp['no_spec']['EUR']['mean']:.4f}, "
                        f"want ~{EXPECTED['no_spec_eur_a1_approx']}")
    if abs(pp["no_spec"]["AFR"]["mean"] - EXPECTED["no_spec_afr_a1_approx"]) > 0.005:
        failures.append(f"no_spec/AFR A1: got {pp['no_spec']['AFR']['mean']:.4f}, "
                        f"want ~{EXPECTED['no_spec_afr_a1_approx']}")
    if abs(pp["cpic_rag"]["EUR"]["mean"] - EXPECTED["cpic_rag_eur_a1_approx"]) > 0.005:
        failures.append(f"cpic_rag/EUR A1: got {pp['cpic_rag']['EUR']['mean']:.4f}, "
                        f"want ~{EXPECTED['cpic_rag_eur_a1_approx']}")

    summ = json.loads(SUMMARY.read_text())
    for cond in CONDITIONS:
        s_a1 = summ["section_1_aggregate"][cond]["A1"]
        all_pop_vals = [pp[cond][p]["mean"] for p in POPULATIONS]
        all_pop_ns = [pp[cond][p]["n"] for p in POPULATIONS]
        weighted = sum(v * n for v, n in zip(all_pop_vals, all_pop_ns)) / sum(all_pop_ns)
        if abs(weighted - s_a1) > 0.001:
            failures.append(f"{cond} weighted-mean cross-check: figure={weighted:.4f}, summary={s_a1:.4f}")

    if failures:
        print("SELF-TEST FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"SELF-TEST PASSED (per-pop n ~{EXPECTED['n_per_cond_per_pop_approx']}; "
          f"with_spec uniform 1.0 across EUR/AMR/AFR; aggregate cross-check vs summary green).")
    return 0


def save_fig(fig: plt.Figure, name: str) -> None:
    png = FIGDIR / f"{name}.png"
    tiff = FIGDIR / f"{name}.tiff"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(tiff, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"wrote {png.relative_to(BASE)}")
    print(f"wrote {tiff.relative_to(BASE)}")


def render() -> None:
    FIGDIR.mkdir(exist_ok=True)
    main = json.loads(MAIN.read_text())
    pp = per_pop_means(main)
    per_gene = per_gene_pop_means(main)
    deltas = eur_afr_deltas(per_gene)

    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 300,
    })

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 5.2),
                                   gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.28})

    x = np.arange(len(CONDITIONS))
    width = 0.26

    for i, pop in enumerate(POPULATIONS):
        offsets = (i - 1) * width
        vals = [100 * pp[c][pop]["mean"] for c in CONDITIONS]
        bars = axA.bar(x + offsets, vals, width,
                       label=pop, color=POPULATION_COLOURS[pop],
                       edgecolor="black", linewidth=0.4)
        for rect, v in zip(bars, vals):
            axA.text(rect.get_x() + rect.get_width() / 2, v + 1.2,
                     f"{v:.1f}", ha="center", va="bottom", fontsize=7)

    axA.set_xticks(x)
    axA.set_xticklabels([CONDITION_LABELS_SHORT[c] for c in CONDITIONS], fontweight="bold")
    axA.set_ylabel("Aggregate phenotype accuracy A1 (%)")
    axA.set_ylim(0, 115)
    axA.set_yticks([0, 25, 50, 75, 100])
    axA.axhline(100, color="grey", linewidth=0.5, linestyle=":")
    axA.set_title("A  Aggregate A1 by condition and population",
                  loc="left", fontweight="bold")
    axA.legend(loc="lower center", bbox_to_anchor=(0.5, -0.36), ncol=3,
               frameon=False, title="Population")
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)

    for i, cond in enumerate(CONDITIONS):
        eur = pp[cond]["EUR"]["mean"]
        afr = pp[cond]["AFR"]["mean"]
        gap = abs(eur - afr) * 100
        axA.text(i, 110, f"EUR-AFR\n{gap:.2f} pp",
                 ha="center", va="bottom", fontsize=7, color="#555555", style="italic")

    rng = np.random.default_rng(7)
    for i, cond in enumerate(CONDITIONS):
        d_pairs = deltas[cond]
        vals = np.array([d * 100 for _, d in d_pairs])
        jitter = rng.uniform(-0.18, 0.18, size=len(vals))
        axB.scatter(np.full_like(vals, i, dtype=float) + jitter, vals,
                    s=44, color=CONDITION_COLOURS[cond],
                    edgecolors="black", linewidths=0.4, alpha=0.85, zorder=3)
        if cond != "with_spec":
            ranked = sorted(d_pairs, key=lambda t: -abs(t[1]))
            for gene, d in ranked[:3]:
                jval = rng.uniform(-0.18, 0.18)
                axB.annotate(gene, (i + jval, d * 100), xytext=(7, 0),
                             textcoords="offset points", fontsize=7, color="#333333",
                             va="center")

    axB.axhline(0, color="#888888", linewidth=0.8, linestyle="--", zorder=1)
    axB.set_xticks(np.arange(len(CONDITIONS)))
    axB.set_xticklabels([CONDITION_LABELS_SHORT[c] for c in CONDITIONS], fontweight="bold")
    axB.set_ylabel("Per-gene EUR - AFR A1 delta (percentage points)")
    axB.set_ylim(-18, 18)
    axB.set_title("B  Per-locus EUR vs AFR A1 difference",
                  loc="left", fontweight="bold")
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    axB.text(0.02, 0.96, "EUR higher", fontsize=7, color="#4a78a8",
             style="italic", transform=axB.transAxes, va="top", ha="left")
    axB.text(0.02, 0.04, "AFR higher", fontsize=7, color="#5a8f5a",
             style="italic", transform=axB.transAxes, va="bottom", ha="left")

    fig.suptitle("Aggregate population gaps are small; per-locus gaps are real and heterogeneous; specification eliminates both",
                 fontsize=10.5, fontweight="bold", color="#333333", y=1.005)

    fig.text(0.5, -0.04,
             "Per-condition n approximately 2,900 parsed evaluations per population (EUR / AMR / AFR). "
             "Per-locus deltas computed across 21 genes; the three largest absolute deltas per condition are annotated.",
             ha="center", fontsize=7, color="#555555")

    save_fig(fig, "Figure6_population_equity")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1] if __doc__ else None)
    parser.add_argument("--self-test", action="store_true",
                        help="Verify inputs; do not render.")
    args = parser.parse_args()

    if args.self_test:
        return self_test()

    rc = self_test()
    if rc != 0:
        print("Aborting render: self-test failed.", file=sys.stderr)
        return rc
    render()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
