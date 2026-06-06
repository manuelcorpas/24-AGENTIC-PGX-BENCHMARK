#!/usr/bin/env python3
"""
Figure 4 for Cell Genomics Manuscript v8.

The "information-without-action" finding: under retrieval-augmented generation,
several HLA loci show high phenotype-identification accuracy (A1 ~ 1.0) but low
clinical-safety action accuracy (A3 << 1.0). The retrieved CPIC text surfaces the
correct status, but the model fails to translate it into the canonical AVOID
recommendation. Specification-constrained execution is the only condition that
aligns A1 and A3 across every lethal-class gene.

Three-panel scatter, one per condition. Each point is a lethal-class gene
plotted at (lethal-class mean A1, lethal-class mean A3). The diagonal x = y is
the alignment line:
  - top-right (1, 1)            = clinical-grade alignment
  - below the diagonal           = information without action (high A1, low A3)
  - above the diagonal           = action without information (low A1, high A3)

Inputs:
  RESULTS/v3_raw_rescored_three_arm.json    (per-gene lethal-class means)
  RESULTS/v3_three_arm_lethal_a3_errors.csv (sanity-check counts)
  SPECS/test_cases_v3.json                  (lethal-flag per case)

Outputs (PNG 300 DPI + TIFF 600 DPI LZW):
  FIGURES/Figure4_information_without_action.png
  FIGURES/Figure4_information_without_action.tiff

Usage:
  python3 34-figure4-information-without-action.py             # render
  python3 34-figure4-information-without-action.py --self-test # verify inputs
"""
from __future__ import annotations

import argparse
import csv
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
LETHAL_CSV = BASE / "RESULTS" / "v3_three_arm_lethal_a3_errors.csv"
CASES_FILE = BASE / "SPECS" / "test_cases_v3.json"
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

INFO_WITHOUT_ACTION_HIGHLIGHT = {"HLA-B*57:01", "HLA-B*15:02", "HLA-A*31:01"}

EXPECTED = {
    "n_lethal_genes": 11,
    "n_lethal_cases": 14,
    "hlab5701_cpic_rag_a1_approx": 1.000,
    "hlab5701_cpic_rag_a3_approx": 0.111,
    "hlab1502_cpic_rag_a3_approx": 0.130,
    "hlaa3101_cpic_rag_a3_approx": 0.321,
    "mtrnr1_a1_zero": True,
}


def load_lethal_cases() -> set[str]:
    cases = json.loads(CASES_FILE.read_text())
    return {c["id"] for c in cases if "lethal" in (c.get("gt_drug") or "").lower()}


def lethal_means_per_gene(main: list[dict], lethal_tcs: set[str]) -> dict:
    by: dict[str, dict[str, dict[str, list[float]]]] = defaultdict(
        lambda: defaultdict(lambda: {"A1": [], "A3": []})
    )
    for r in main:
        if r["scores"].get("format_fail"):
            continue
        if r["tc"] not in lethal_tcs:
            continue
        cond = r["cond"]
        gene = r["gene"]
        a1 = r["scores"].get("A1")
        a3 = r["scores"].get("A3")
        if a1 is not None:
            by[gene][cond]["A1"].append(a1)
        if a3 is not None:
            by[gene][cond]["A3"].append(a3)

    out: dict[str, dict[str, dict[str, float]]] = {}
    for gene, cmap in by.items():
        out[gene] = {}
        for cond, vals in cmap.items():
            a1_list = vals["A1"]
            a3_list = vals["A3"]
            out[gene][cond] = {
                "A1": sum(a1_list) / len(a1_list) if a1_list else float("nan"),
                "A3": sum(a3_list) / len(a3_list) if a3_list else float("nan"),
                "n": len(a3_list),
            }
    return out


def self_test() -> int:
    failures: list[str] = []

    lethal_tcs = load_lethal_cases()
    if len(lethal_tcs) != EXPECTED["n_lethal_cases"]:
        failures.append(f"lethal case count: got {len(lethal_tcs)}, want {EXPECTED['n_lethal_cases']}")

    main = json.loads(MAIN.read_text())
    means = lethal_means_per_gene(main, lethal_tcs)
    if len(means) != EXPECTED["n_lethal_genes"]:
        failures.append(f"lethal gene count: got {len(means)}, want {EXPECTED['n_lethal_genes']}")

    rag = means.get("HLA-B*57:01", {}).get("cpic_rag", {})
    if abs(rag.get("A1", -1) - EXPECTED["hlab5701_cpic_rag_a1_approx"]) > 0.005:
        failures.append(f"HLA-B*57:01 cpic_rag A1: got {rag.get('A1')}, want ~{EXPECTED['hlab5701_cpic_rag_a1_approx']}")
    if abs(rag.get("A3", -1) - EXPECTED["hlab5701_cpic_rag_a3_approx"]) > 0.005:
        failures.append(f"HLA-B*57:01 cpic_rag A3: got {rag.get('A3')}, want ~{EXPECTED['hlab5701_cpic_rag_a3_approx']}")

    rag1502 = means.get("HLA-B*15:02", {}).get("cpic_rag", {})
    if abs(rag1502.get("A3", -1) - EXPECTED["hlab1502_cpic_rag_a3_approx"]) > 0.005:
        failures.append(f"HLA-B*15:02 cpic_rag A3: got {rag1502.get('A3')}, want ~{EXPECTED['hlab1502_cpic_rag_a3_approx']}")

    rag3101 = means.get("HLA-A*31:01", {}).get("cpic_rag", {})
    if abs(rag3101.get("A3", -1) - EXPECTED["hlaa3101_cpic_rag_a3_approx"]) > 0.005:
        failures.append(f"HLA-A*31:01 cpic_rag A3: got {rag3101.get('A3')}, want ~{EXPECTED['hlaa3101_cpic_rag_a3_approx']}")

    mt = means.get("MT-RNR1", {})
    for cond in ("no_spec", "cpic_rag"):
        if abs(mt.get(cond, {}).get("A1", 1.0)) > 0.01:
            failures.append(f"MT-RNR1 {cond} A1 (expected ~0.000): got {mt.get(cond, {}).get('A1')}")

    for gene, cmap in means.items():
        ws = cmap.get("with_spec", {})
        if abs(ws.get("A1", 0) - 1.0) > 0.001 or abs(ws.get("A3", 0) - 1.0) > 0.001:
            failures.append(f"with_spec ceiling A1=A3=1.0 violated at {gene}: A1={ws.get('A1')}, A3={ws.get('A3')}")

    sum_a3_err_no_spec = 0
    sum_a3_err_cpic_rag = 0
    with LETHAL_CSV.open() as fh:
        for row in csv.DictReader(fh):
            sum_a3_err_no_spec += int(row["no_spec_a3_errors"])
            sum_a3_err_cpic_rag += int(row["cpic_rag_a3_errors"])
    if sum_a3_err_no_spec != 270:
        failures.append(f"lethal A3 error sum no_spec: got {sum_a3_err_no_spec}, want 270")
    if sum_a3_err_cpic_rag != 414:
        failures.append(f"lethal A3 error sum cpic_rag: got {sum_a3_err_cpic_rag}, want 414")

    if failures:
        print("SELF-TEST FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"SELF-TEST PASSED ({EXPECTED['n_lethal_genes']} lethal-class genes; "
          f"HLA-B*57:01 cpic_rag A1/A3 = {rag['A1']:.3f}/{rag['A3']:.3f}; "
          f"MT-RNR1 A1 ~ 0; with_spec A1=A3=1.0 on every lethal gene; "
          f"lethal A3 error sums 270 / 414 / 0).")
    return 0


def save_fig(fig: plt.Figure, name: str) -> None:
    png = FIGDIR / f"{name}.png"
    tiff = FIGDIR / f"{name}.tiff"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(tiff, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"wrote {png.relative_to(BASE)}")
    print(f"wrote {tiff.relative_to(BASE)}")


def label_offsets(cond: str) -> dict:
    """Per-condition, per-gene manual offsets to avoid overlap. (dx, dy, ha, va)."""
    if cond == "no_spec":
        return {
            "CYP2D6":     (0.028, 0.020, "left",  "bottom"),
            "HLA-B*58:01":(-0.025, -0.030, "right", "top"),
            "RYR1":       (0.028, -0.005, "left",  "center"),
            "G6PD":       (-0.025, 0.035, "right", "bottom"),
            "HLA-B*57:01":(0.028, 0.020, "left",  "bottom"),
            "HLA-A*31:01":(-0.025, -0.020, "right", "top"),
            "HLA-B*15:02":(0.028, -0.020, "left",  "top"),
            "MT-RNR1":    (0.030, 0.000, "left",  "center"),
            "DPYD":       (0.028, 0.000, "left",  "center"),
            "NUDT15":     (0.028, 0.000, "left",  "center"),
            "TPMT":       (0.028, 0.000, "left",  "center"),
        }
    if cond == "cpic_rag":
        return {
            "CYP2D6":     (0.028, 0.025, "left",  "bottom"),
            "HLA-B*58:01":(-0.025, 0.025, "right", "bottom"),
            "RYR1":       (0.028, 0.000, "left",  "center"),
            "G6PD":       (-0.025, -0.025, "right", "top"),
            "HLA-B*57:01":(0.030, -0.020, "left",  "top"),
            "HLA-A*31:01":(0.030, 0.020, "left",  "bottom"),
            "HLA-B*15:02":(0.030, 0.020, "left",  "bottom"),
            "MT-RNR1":    (0.030, 0.000, "left",  "center"),
            "DPYD":       (0.028, 0.000, "left",  "center"),
            "NUDT15":     (-0.028, 0.000, "right", "center"),
            "TPMT":       (0.028, 0.000, "left",  "center"),
        }
    return {}


def render() -> None:
    FIGDIR.mkdir(exist_ok=True)
    lethal_tcs = load_lethal_cases()
    main = json.loads(MAIN.read_text())
    means = lethal_means_per_gene(main, lethal_tcs)

    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 300,
    })

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.0), sharex=True, sharey=True)

    for ax, cond in zip(axes, CONDITIONS):
        ax.plot([-0.05, 1.05], [-0.05, 1.05], color="#bbbbbb",
                linewidth=0.8, linestyle="--", zorder=1)
        ax.fill_between([-0.05, 1.05], [-0.05, 1.05], [-0.10, -0.10],
                        color="#f5e6e2", alpha=0.4, zorder=0)
        ax.fill_between([-0.05, 1.05], [-0.05, 1.05], [1.15, 1.15],
                        color="#e6efe6", alpha=0.4, zorder=0)

        offsets = label_offsets(cond)

        for gene, cmap in means.items():
            d = cmap.get(cond)
            if not d:
                continue
            x, y = d["A1"], d["A3"]
            base_colour = CONDITION_COLOURS[cond]
            highlight = (cond == "cpic_rag" and gene in INFO_WITHOUT_ACTION_HIGHLIGHT)
            marker_colour = "#7a2820" if highlight else base_colour
            edge = "black" if highlight else "white"
            size = 130 if highlight else 80
            ax.scatter([x], [y], s=size, c=marker_colour, edgecolors=edge,
                       linewidths=1.0 if highlight else 0.6, zorder=3)

            if cond != "with_spec":
                dx, dy, ha, va = offsets.get(gene, (0.020, 0.020, "left", "bottom"))
                ax.annotate(gene, (x + dx, y + dy),
                            fontsize=7.5, color="#333333", ha=ha, va=va,
                            fontweight="bold" if highlight else "normal")

        ax.set_xlim(-0.08, 1.18)
        ax.set_ylim(-0.10, 1.15)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlabel("Lethal-class phenotype accuracy (A1)")
        ax.set_title(CONDITION_LABELS[cond], color=CONDITION_COLOURS[cond], fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_aspect("equal", adjustable="box")

    axes[0].set_ylabel("Lethal-class safety action (A3)")

    axes[2].annotate("All 11 lethal-class genes\nA1 = A3 = 1.000\n(every cell, every replicate,\nevery population)",
                     xy=(1.0, 1.0), xytext=(0.45, 0.50),
                     fontsize=8.5, color="#2a5a2a", fontweight="bold",
                     ha="center", va="center",
                     arrowprops=dict(arrowstyle="->", color="#2a5a2a", lw=1.0))

    fig.suptitle("Same phenotype, different action: retrieval surfaces correct information; specification translates it into correct action",
                 fontsize=10.5, fontweight="bold", color="#333333", y=1.015)

    fig.text(0.5, 0.96,
             "Above diagonal: action without information   |   Below diagonal: information without action",
             ha="center", fontsize=8, style="italic", color="#666666")

    fig.text(0.5, -0.05,
             "11 lethal-class genes (14 lethal cases). Per-gene means computed across all parsed lethal-class cells "
             "(populations and replicates pooled). Highlighted points in the Retrieval-augmented panel "
             "(HLA-B*57:01, HLA-B*15:02, HLA-A*31:01) show the canonical information-without-action pattern: "
             "phenotype correctly identified, AVOID action not produced.",
             ha="center", fontsize=7, color="#555555")

    save_fig(fig, "Figure4_information_without_action")


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
