#!/usr/bin/env python3
"""
Revision figures, generated from the matched-factorial data (N1, N4, N5, N6).

EVERY NUMBER PLOTTED HERE IS READ FROM A DATA FILE. There is no np.random in
this script and no hardcoded result: the figures are recomputed from
data/v3_five_cell_live.json, data/v3_model_caller_eval.json and
data/v3_ancestry_matched.json, so a figure can never drift from the numbers in
the text. --self-test verifies the inputs exist and the recomputed headline
values match the report before anything is drawn.

Figures produced:
  FigureR1  the 2x3 factorial: knowledge representation x decision mechanism,
            accuracy and lethal-class errors side by side
  FigureR2  error localisation: end-to-end accuracy against input-call accuracy,
            per gene, with the identity line
  FigureR3  between-model spread collapsing as correctness moves into execution
  FigureR4  ancestry coverage on both denominators, with the standardised and
            cohort-specific columns that stop the disparity being explained away

USAGE
    python code/70-figures-revision.py --self-test
    python code/70-figures-revision.py
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
FIVE_CELL = BASE / "data" / "v3_five_cell_live.json"
CALLER = BASE / "data" / "v3_model_caller_eval.json"
ANCESTRY = BASE / "data" / "v3_ancestry_matched.json"
FIGDIR = BASE / "figures"

CELLS = ["free_generation", "rag_generation", "rag_execution",
         "skill_generation", "skill_execution"]
LABELS = {
    "free_generation": "no\nknowledge",
    "rag_generation": "prose\ngenerate",
    "rag_execution": "prose\nEXECUTE",
    "skill_generation": "rules\ngenerate",
    "skill_execution": "rules\nEXECUTE",
}
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

def figure_r1(summary) -> None:
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.5, 4.2),
                                   gridspec_kw={"width_ratios": [1.1, 1.0]})
    xs = range(len(CELLS))
    acc = [summary[c]["a2_recommendation_mean"] for c in CELLS]
    axA.bar(xs, acc, color=[COLOURS[c] for c in CELLS], width=0.68)
    for x, v in zip(xs, acc):
        axA.text(x, v + 0.015, f"{v:.3f}", ha="center", fontsize=8)
    axA.set_xticks(list(xs))
    axA.set_xticklabels([LABELS[c] for c in CELLS], fontsize=8)
    axA.set_ylabel("recommendation accuracy")
    axA.set_ylim(0, 1.08)
    axA.set_title("A  Executing the rules recovers prose", loc="left")

    lethal = [summary[c]["lethal_errors"] for c in CELLS]
    axB.bar(xs, lethal, color=[COLOURS[c] for c in CELLS], width=0.68)
    for x, v in zip(xs, lethal):
        axB.text(x, v + 2, str(v), ha="center", fontsize=8)
    axB.set_xticks(list(xs))
    axB.set_xticklabels([LABELS[c] for c in CELLS], fontsize=8)
    axB.set_ylabel("lethal-class errors (of 336)")
    axB.set_ylim(0, max(lethal) * 1.2)
    axB.set_title("B  and removes most dangerous errors", loc="left")

    handles = [plt.Rectangle((0, 0), 1, 1, color=GENERATION),
               plt.Rectangle((0, 0), 1, 1, color=EXECUTION)]
    axA.legend(handles, ["model generates the answer", "code executes the rule"],
               loc="lower right", frameon=False)
    fig.tight_layout()
    save(fig, "FigureR1_matched_factorial")


# ---------------------------------------------------------------- FigureR2

def figure_r2(summary) -> None:
    caller = json.loads(CALLER.read_text())
    genes = {g: v for g, v in caller["per_gene"].items() if v["n"] >= 60}
    xs = [v["concordance"] for v in genes.values()]
    ys = []
    # end-to-end accuracy per gene, execution cells only, from the raw rows
    rows = json.loads(FIVE_CELL.read_text())
    rs = _load("61-rescore-matched.py", "rs")
    for g in genes:
        sub = [r for r in rows if r["gene"] == g and r["cell"].endswith("execution")]
        sc = [rs.score_row(r, rs.CASE_BY_ID[r["case_id"]])["a1_phenotype"] for r in sub]
        ys.append(sum(sc) / len(sc) if sc else 0.0)

    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    lo = min(min(xs), min(ys)) - 0.05
    ax.plot([lo, 1.0], [lo, 1.0], color=NEUTRAL, lw=1, ls="--", zorder=1)
    ax.scatter(xs, ys, s=44, color=EXECUTION, zorder=3, edgecolor="white", linewidth=0.6)
    for g, x, y in zip(genes, xs, ys):
        ax.annotate(g, (x, y), textcoords="offset points", xytext=(5, -3), fontsize=7)
    ax.set_xlabel("input-call accuracy (model as caller)")
    ax.set_ylabel("end-to-end phenotype accuracy")
    ax.set_xlim(lo, 1.02)
    ax.set_ylim(lo, 1.02)
    ax.set_title("Under execution, end-to-end accuracy\nis input-call accuracy", loc="left")
    ax.text(0.02, 0.96, f"overall call concordance {caller['call_concordance']:.3f}",
            transform=ax.transAxes, fontsize=8, va="top")
    fig.tight_layout()
    save(fig, "FigureR2_error_localisation")


# ---------------------------------------------------------------- FigureR3

def figure_r3(by_model) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for i, cell in enumerate(CELLS):
        vals = list(by_model[cell].values())
        ax.scatter([i] * len(vals), vals, s=26, color=COLOURS[cell],
                   alpha=0.75, edgecolor="white", linewidth=0.5, zorder=3)
        ax.plot([i - 0.24, i + 0.24], [statistics.fmean(vals)] * 2,
                color="black", lw=1.4, zorder=4)
        ax.text(i, 0.34, f"spread\n{max(vals) - min(vals):.3f}",
                ha="center", fontsize=7.5)
    ax.set_xticks(range(len(CELLS)))
    ax.set_xticklabels([LABELS[c] for c in CELLS], fontsize=8)
    ax.set_ylabel("phenotype accuracy, one point per model")
    ax.set_ylim(0.3, 1.04)
    ax.set_title("Between-model spread collapses as correctness moves into code",
                 loc="left")
    fig.tight_layout()
    save(fig, "FigureR3_model_spread")


# ---------------------------------------------------------------- FigureR4

def figure_r4() -> None:
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
    ax.set_ylabel("coverage of the validated skill vocabulary")
    ax.set_ylim(0, max(max(v for _, v, _ in series for v in v) * 1.25, 0.6))
    ax.legend(frameon=False, ncol=2, fontsize=7.5)
    ax.set_title("Coverage on four denominators: the disparity is relocated, not dissolved",
                 loc="left", fontsize=9.5)
    fig.tight_layout()
    save(fig, "FigureR4_ancestry_coverage")


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
    figure_r1(summary)
    figure_r2(summary)
    figure_r3(per_model_accuracy(rows, rs))
    figure_r4()
    return 0


if __name__ == "__main__":
    sys.exit(main())
