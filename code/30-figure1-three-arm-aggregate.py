#!/usr/bin/env python3
"""
Figure 1 for Cell Genomics Manuscript v8.

Two-panel figure illustrating the three-arm headline:
  Panel A: aggregate phenotype/recommendation/safety accuracy (A1, A2, A3) by condition.
  Panel B: total lethal-class A3 error count by condition.

The headline story: free-prompted -> retrieval-augmented -> specification-constrained
moves accuracy from 80.6% to 100% on A1, and (paradoxically) increases lethal-class
A3 errors from 270 to 414 under RAG before collapsing to 0 under spec.

Inputs:
  RESULTS/v3_three_arm_analysis_summary.json  (machine-readable summary; primary)
  RESULTS/v3_three_arm_per_case_a1.csv        (used for self-test cross-check)
  RESULTS/v3_three_arm_lethal_a3_errors.csv   (per-gene, used by Figure 1 supplement)

Outputs (PNG 300 DPI + TIFF 600 DPI LZW, Cell Press requirements):
  FIGURES/Figure1_three_arm_aggregate.png
  FIGURES/Figure1_three_arm_aggregate.tiff

Usage:
  python3 30-figure1-three-arm-aggregate.py             # render figure
  python3 30-figure1-three-arm-aggregate.py --self-test # verify headline numbers
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASE = Path(__file__).resolve().parent.parent
SUMMARY = BASE / "RESULTS" / "v3_three_arm_analysis_summary.json"
PER_CASE_A1 = BASE / "RESULTS" / "v3_three_arm_per_case_a1.csv"
LETHAL_BY_GENE = BASE / "RESULTS" / "v3_three_arm_lethal_a3_errors.csv"
FIGDIR = BASE / "FIGURES"

CONDITIONS = ["no_spec", "cpic_rag", "with_spec"]
CONDITION_LABELS = {
    "no_spec": "Free-prompted",
    "cpic_rag": "Retrieval-augmented",
    "with_spec": "Specification-constrained",
}
CONDITION_COLOURS = {
    "no_spec": "#d4756f",
    "cpic_rag": "#e0a96d",
    "with_spec": "#5a8f5a",
}

EXPECTED_HEADLINES = {
    "A1": {"no_spec": 0.806, "cpic_rag": 0.895, "with_spec": 1.000},
    "A2": {"no_spec": 0.616, "cpic_rag": 0.530, "with_spec": 1.000},
    "A3": {"no_spec": 0.969, "cpic_rag": 0.953, "with_spec": 1.000},
    "lethal_a3_errors": {"no_spec": 270, "cpic_rag": 414, "with_spec": 0},
}


def load_summary() -> dict:
    return json.loads(SUMMARY.read_text())


def self_test() -> int:
    summary = load_summary()
    failures: list[str] = []

    agg = summary["section_1_aggregate"]
    for metric in ("A1", "A2", "A3"):
        for cond in CONDITIONS:
            got = round(agg[cond][metric], 3)
            want = EXPECTED_HEADLINES[metric][cond]
            if abs(got - want) > 0.002:
                failures.append(f"aggregate {metric}/{cond}: got {got}, want {want}")

    lethal = summary["section_2_lethal"]
    for cond in CONDITIONS:
        got = lethal[cond]["A3_errors"]
        want = EXPECTED_HEADLINES["lethal_a3_errors"][cond]
        if got != want:
            failures.append(f"lethal A3 errors/{cond}: got {got}, want {want}")

    with PER_CASE_A1.open() as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != 110:
        failures.append(f"per-case A1 row count: got {len(rows)}, want 110")

    with LETHAL_BY_GENE.open() as fh:
        gene_rows = list(csv.DictReader(fh))
    total_no_spec = sum(int(r["no_spec_a3_errors"]) for r in gene_rows)
    total_cpic_rag = sum(int(r["cpic_rag_a3_errors"]) for r in gene_rows)
    total_with_spec = sum(int(r["with_spec_a3_errors"]) for r in gene_rows)
    if total_no_spec != 270:
        failures.append(f"per-gene no_spec lethal total: got {total_no_spec}, want 270")
    if total_cpic_rag != 414:
        failures.append(f"per-gene cpic_rag lethal total: got {total_cpic_rag}, want 414")
    if total_with_spec != 0:
        failures.append(f"per-gene with_spec lethal total: got {total_with_spec}, want 0")

    if failures:
        print("SELF-TEST FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("SELF-TEST PASSED (12 assertions: 9 aggregate means, 3 lethal counts, row totals).")
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
    summary = load_summary()
    agg = summary["section_1_aggregate"]
    lethal = summary["section_2_lethal"]

    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 300,
    })

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.5, 4.0), gridspec_kw={"width_ratios": [1.6, 1.0]})

    metrics = ["A1", "A2", "A3"]
    metric_labels = [
        "A1\nPhenotype",
        "A2\nDrug recommendation",
        "A3\nSafety action",
    ]
    x = np.arange(len(metrics))
    width = 0.26

    for i, cond in enumerate(CONDITIONS):
        offsets = (i - 1) * width
        vals = [100 * agg[cond][m] for m in metrics]
        bars = axA.bar(x + offsets, vals, width,
                       label=CONDITION_LABELS[cond],
                       color=CONDITION_COLOURS[cond],
                       edgecolor="black", linewidth=0.4)
        for rect, v in zip(bars, vals):
            axA.text(rect.get_x() + rect.get_width() / 2, v + 1.2,
                     f"{v:.1f}", ha="center", va="bottom", fontsize=7)

    axA.set_xticks(x)
    axA.set_xticklabels(metric_labels)
    axA.set_ylabel("Mean accuracy (%)")
    axA.set_ylim(0, 115)
    axA.set_yticks([0, 25, 50, 75, 100])
    axA.axhline(100, color="grey", linewidth=0.5, linestyle=":")
    axA.set_title("A  Aggregate accuracy by condition", loc="left", fontweight="bold")
    axA.legend(loc="lower center", bbox_to_anchor=(0.5, -0.32), ncol=3, frameon=False)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)

    cond_x = np.arange(len(CONDITIONS))
    err_vals = [lethal[c]["A3_errors"] for c in CONDITIONS]
    err_colours = [CONDITION_COLOURS[c] for c in CONDITIONS]
    bars = axB.bar(cond_x, err_vals, color=err_colours,
                   edgecolor="black", linewidth=0.4, width=0.62)
    for rect, v in zip(bars, err_vals):
        axB.text(rect.get_x() + rect.get_width() / 2, v + 12,
                 str(v), ha="center", va="bottom", fontsize=9, fontweight="bold")

    axB.set_xticks(cond_x)
    axB.set_xticklabels([CONDITION_LABELS[c].replace("-", "-\n", 1) for c in CONDITIONS],
                        fontsize=8)
    axB.set_ylabel("Lethal-class A3 errors (count)")
    axB.set_ylim(0, max(err_vals) * 1.25 + 10)
    axB.set_title("B  Lethal-class clinical-safety errors", loc="left", fontweight="bold")
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    # Sample-size footnote removed per co-author review (the n is stated in the figure legend).

    fig.tight_layout()
    save_fig(fig, "Figure1_three_arm_aggregate")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1] if __doc__ else None)
    parser.add_argument("--self-test", action="store_true",
                        help="Verify headline numbers in inputs match locked values; do not render.")
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
