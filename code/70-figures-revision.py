#!/usr/bin/env python3
"""Revision figures computed from preserved benchmark and cohort data.

Figure R1 reads the R3 attempted-evaluation summary: original-reference
recommendation concordance above both wording-rule counts. R2 localises
phenotype disagreement against input-call concordance. R3 reports between-model
phenotype concordance. R4 retains the four ancestry-coverage denominators.
Use docs/R3-2026-09-16/SOURCE/01-regenerate.py for the complete current pipeline.
Historical output paths remain available; no public figure is overwritten by R3.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
FIVE_CELL = BASE / "data" / "v4_flash_eight_model_five_cell.json"
CALLER = BASE / "data" / "v3_model_caller_eval.json"
ANCESTRY = BASE / "data" / "v3_ancestry_four_cohorts.json"
FIGDIR = BASE / "figures"

CELLS = ["free_generation", "rag_generation", "rag_execution",
         "skill_generation", "skill_execution"]
LABELS = {
    "free_generation": "No knowledge",
    "rag_generation": "RAG generation",
    "rag_execution": "RAG call + authored execution",
    "skill_generation": "Authored rules, generated",
    "skill_execution": "Authored rules, executed",
}
FIGURE_R3_LABELS = [
    "No\nknowledge",
    "RAG\ngeneration",
    "RAG call +\nauthored execution",
    "Authored rules,\ngenerated",
    "Authored rules,\nexecuted",
]
FIGURE_R2_LABEL_OFFSETS = {
    "HLA-B*15:02": (-8, 10),
    "NUDT15": (-38, -14),
    "CYP2C19": (-48, -3),
    "CYP4F2": (6, -12),
    "CYP2D6": (5, 8),
    "CYP3A5": (5, -8),
    "G6PD": (5, 8),
    "DPYD": (5, 8),
    "SLCO1B1": (5, -8),
}
FIGURE_R2_CLUSTER_CUTOFF = (0.97, 0.98)
GENERATION = "#d4756f"
EXECUTION = "#5a8f5a"
NEUTRAL = "#8a8a8a"
COLOURS = {c: (EXECUTION if c.endswith("execution") else GENERATION) for c in CELLS}


def _load(name, mod):
    spec = spec_from_file_location(mod, CODE / name)
    m = module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def scored_cells() -> dict:
    """Recompute per-cell headline numbers from the raw scored rows."""
    rs = _load("61-rescore-matched.py", "rs")
    rows = json.loads(FIVE_CELL.read_text())
    return rs.summarise(rows), rows, rs


def per_model_accuracy(rows, rs) -> dict:
    out = defaultdict(dict)
    models = sorted({r["model"] for r in rows})
    for cell in CELLS:
        for m in models:
            sub = [r for r in rows if r["model"] == m and r["cell"] == cell]
            if not sub:
                continue
            sc = [rs.score_row(r, rs.CASE_BY_ID[r["case_id"]])["a1_phenotype"] for r in sub]
            out[cell][m] = sum(sc) / len(sc)
    return out


def style() -> None:
    plt.rcParams.update({
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "figure.dpi": 300,
    })


def save(fig, stem: str) -> None:
    FIGDIR.mkdir(exist_ok=True)
    fig.savefig(FIGDIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / f"{stem}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  wrote figures/{stem}.png and .tiff")


# ---------------------------------------------------------------- FigureR1

def figure_r1(summary, save_output: bool = True):
    """Original-reference concordance and both wording rules on attempted units."""
    blue, orange = "#3775BA", "#D95F02"
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(7.0, 6.6))
    labels = ["Free generation", "RAG generation", "RAG-assisted execution",
              "Authored-rule generation", "Authored-rule execution"]
    ys = list(range(len(CELLS)))
    values = [100 * summary[c]["a2_recommendation_mean"] for c in CELLS]
    bars = axA.barh(ys, values, height=0.56, color=blue, edgecolor="#222222", linewidth=0.6)
    for y, value in zip(ys, values):
        axA.text(value + 1.5, y, f"{value:.1f}%", va="center", fontsize=8.5)
    axA.set_xlim(0, 113)
    axA.set_xticks([0, 25, 50, 75, 100])
    axA.set_xlabel("Recommendation concordance with authored references (%)", fontsize=8.5)
    axA.set_title("A  Original references, 110 cases", loc="left", fontsize=10, weight="bold")
    for offset, key, colour, hatch, label in [
        (-0.17, "strict_failures", blue, "", "Original rule"),
        (0.17, "widened_failures", orange, "///", 'Also accepts “do not use” (post hoc)'),
    ]:
        vals = [summary[c][key] for c in CELLS]
        axB.barh([y + offset for y in ys], vals, height=0.29, color=colour,
                 edgecolor="#222222", linewidth=0.6, hatch=hatch, label=label)
        for y, value in zip(ys, vals):
            axB.text(value + 2, y + offset, str(value), va="center", fontsize=8.5)
    axB.set_xlim(0, max(summary[c]["strict_failures"] for c in CELLS) * 1.16)
    axB.set_xlabel("Wording-based failures in high-severity cases (of 336)", fontsize=8.5)
    axB.set_title("B  Sensitivity to accepted wording", loc="left", fontsize=10, weight="bold")
    for ax in (axA, axB):
        ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=8.5)
        ax.invert_yaxis(); ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_axisbelow(True); ax.grid(axis="x", color="#e6e6e6", linewidth=0.5)
    fig.legend(*axB.get_legend_handles_labels(), loc="lower center", ncol=1,
               frameon=False, fontsize=8.5, bbox_to_anchor=(0.62, 0.002))
    fig.subplots_adjust(left=0.29, right=0.97, top=0.94, bottom=0.16, hspace=0.64)
    if save_output:
        save(fig, "FigureR1_matched_factorial_r3")
    return fig


# ---------------------------------------------------------------- FigureR2

def figure_r2(summary, rows=None, already_scored=False, save_output=True, caller=None):
    caller = json.loads(CALLER.read_text()) if caller is None else caller
    genes = {g: v for g, v in caller["per_gene"].items() if v["n"] >= 60}
    xs = [v["concordance"] for v in genes.values()]
    ys = []
    # end-to-end accuracy per gene, execution cells only, from the raw rows
    rows = json.loads(FIVE_CELL.read_text()) if rows is None else rows
    rs = _load("61-rescore-matched.py", "rs")
    for g in genes:
        sub = [r for r in rows if r["gene"] == g and r["cell"].endswith("execution")]
        sc = [r['a1_phenotype'] if already_scored else rs.score_row(r, rs.CASE_BY_ID[r["case_id"]])["a1_phenotype"] for r in sub]
        ys.append(sum(sc) / len(sc) if sc else 0.0)

    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    lo = min(min(xs), min(ys)) - 0.05
    ax.plot([lo, 1.0], [lo, 1.0], color=NEUTRAL, lw=1, ls="--", zorder=1)
    ax.scatter(xs, ys, s=44, color=EXECUTION, zorder=3, edgecolor="white", linewidth=0.6)
    high_cluster = []
    for g, x, y in zip(genes, xs, ys):
        if x >= FIGURE_R2_CLUSTER_CUTOFF[0] and y >= FIGURE_R2_CLUSTER_CUTOFF[1]:
            high_cluster.append(g)
            continue
        offset = FIGURE_R2_LABEL_OFFSETS.get(g, (5, -3))
        ax.annotate(g, (x, y), textcoords="offset points", xytext=offset,
                    fontsize=6.8)
    if high_cluster:
        ax.text(0.98, 0.08, f"{len(high_cluster)} genes in high-concordance cluster",
                transform=ax.transAxes, ha="right", fontsize=6.8)
    ax.set_xlabel("Input-call concordance (model as caller)")
    ax.set_ylabel("Phenotype concordance with authored references")
    ax.set_xlim(lo, 1.02)
    ax.set_ylim(lo, 1.02)
    ax.set_title("Within authored execution, phenotype concordance\ntracks input-call concordance", loc="left")
    ax.text(0.02, 0.96, f"overall call concordance {caller['call_concordance']:.3f}",
            transform=ax.transAxes, fontsize=8, va="top")
    fig.tight_layout()
    if save_output:save(fig, "FigureR2_error_localisation")
    return fig


# ---------------------------------------------------------------- FigureR3

def figure_r3(by_model, save_output=True):
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for i, cell in enumerate(CELLS):
        vals = list(by_model[cell].values())
        ax.scatter([i] * len(vals), vals, s=26, color=COLOURS[cell],
                   alpha=0.75, edgecolor="white", linewidth=0.5, zorder=3)
        ax.plot([i - 0.24, i + 0.24], [statistics.fmean(vals)] * 2,
                color="black", lw=1.4, zorder=4)
        ax.text(i, 0.34, f"spread\n{max(vals) - min(vals):.3f}",
                ha="center", fontsize=7.5)
    ax.set_xticks(range(len(CELLS)))
    ax.set_xticklabels(FIGURE_R3_LABELS, fontsize=7.5, linespacing=1.0)
    ax.set_ylabel("Phenotype concordance, one point per model")
    ax.set_ylim(0.3, 1.04)
    ax.set_title("Between-model spread across the five matched configurations",
                 loc="left")
    fig.subplots_adjust(left=0.12, right=0.98, top=0.90, bottom=0.20)
    if save_output:save(fig, "FigureR3_model_spread")
    return fig


# ---------------------------------------------------------------- FigureR4

def figure_r4(save_output=True):
    a = json.loads(ANCESTRY.read_text())
    cohorts = a["cohorts"]
    series = [
        ("distinct states", [a["raw_by_cohort"][c]["coverage_states"] for c in cohorts], "#d4756f"),
        ("carriers", [a["raw_by_cohort"][c]["coverage_carriers"] for c in cohorts], "#e0a96d"),
        ("standardised", [a["standardised"][c]["standardised_coverage"] for c in cohorts], "#5a8f5a"),
        ("cohort-specific states", [a["cohort_specific"][c]["coverage_states"] for c in cohorts], "#4a6f9a"),
    ]
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    width = 0.2
    for k, (label, vals, colour) in enumerate(series):
        xs = [i + (k - 1.5) * width for i in range(len(cohorts))]
        ax.bar(xs, vals, width=width, label=label, color=colour)
        for x, v in zip(xs, vals):
            ax.text(x, v + 0.012, f"{v:.2f}", ha="center", fontsize=6.6)
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels(cohorts)
    ax.set_ylabel("Coverage of the authored skill vocabulary")
    ax.set_ylim(0, max(max(v for _, v, _ in series for v in v) * 1.25, 0.6))
    ax.legend(frameon=False, ncol=2, fontsize=7.5)
    ax.set_title("Coverage on four denominators: the disparity is relocated, not dissolved",
                 loc="left", fontsize=9.5)
    fig.tight_layout()
    if save_output:save(fig, "FigureR4_ancestry_coverage")
    return fig


# ---------------------------------------------------------------- entry

def self_test() -> int:
    missing = [p.name for p in (FIVE_CELL, CALLER, ANCESTRY) if not p.exists()]
    if missing:
        sys.stderr.write(f"missing input data: {missing}\n")
        return 1
    summary, rows, _ = scored_cells()
    failures = []
    for cell in CELLS:
        if cell not in summary:
            failures.append(f"{cell} absent from the scored data")
    exec_cells = [c for c in CELLS if c.endswith("execution")]
    if not all(summary[c]["lethal_errors"] < summary["rag_generation"]["lethal_errors"]
               for c in exec_cells):
        failures.append("execution cells no longer have fewer lethal errors than "
                        "rag_generation; the figure caption would be wrong")
    caller = json.loads(CALLER.read_text())
    if not 0.0 <= caller["call_concordance"] <= 1.0:
        failures.append("call concordance out of range")
    if failures:
        for f in failures:
            sys.stderr.write(f"SELF-TEST FAIL: {f}\n")
        return 1
    print("SELF-TEST PASSED: inputs present and headline relations hold")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        return self_test()
    if self_test() != 0:
        return 1
    style()
    summary, rows, rs = scored_cells()
    print("rendering revision figures from data files:")
    figure_r1(json.loads((BASE / "docs/R3-2026-09-16/DATA/analysis.json").read_text())["baseline"])
    figure_r2(summary)
    figure_r3(per_model_accuracy(rows, rs))
    figure_r4()
    return 0


if __name__ == "__main__":
    sys.exit(main())
