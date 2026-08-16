#!/usr/bin/env python3
"""Render Figure 7 from the same baseline scoring used in the manuscript text."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data" / "v3_agent_vs_deterministic.json"
FIGURES = BASE / "figures"

ORDER = [
    ("CorpasFamily", "European\n(family WGS, n=4)"),
    ("Peru", "Latin American\n(Peru)"),
    ("1000G_IBS", "European\n(Iberian, n=93)"),
    ("UGR", "East African\n(Uganda)"),
]


SCORED = BASE / "data" / "v3_matched_scored_rows_all5.json"
# The benchmark bar is authored-rule generation: the cell where the model holds
# the rule table and produces the phenotype itself, which is the mechanism the
# cohort bars test. The execution cells are the opposite handoff (model supplies
# the diplotype, code maps it), so they are not the analogue.
BENCHMARK_CELL = "skill_generation"
BENCHMARK_ATTEMPTED = 2640


def benchmark_bar():
    """Score the curated benchmark under the same three-way scheme as the cohorts.

    These three values were literals: 0.96, 0.0 and 0.04. Nothing in the data is
    0.96, the two round numbers beside it were placeholders, and the panel was
    therefore unreproducible from its own generator for the one bar every other
    bar is compared against. Computing it here is what makes the figure's
    headline claim checkable.
    """
    rows = [r for r in json.loads(SCORED.read_text())
            if r["cell"] == BENCHMARK_CELL]
    if not rows:
        raise AssertionError(f"no scored rows for {BENCHMARK_CELL}")
    n = BENCHMARK_ATTEMPTED
    correct = sum(r["a1_phenotype"] for r in rows) / n
    abstain = sum(1 for r in rows if r["abstained"]) / n
    wrong = 1.0 - correct - abstain
    if wrong < 0:
        raise AssertionError("negative wrong category on the benchmark bar")
    return correct, abstain, wrong


def values(data):
    b_correct, b_abstain, b_wrong = benchmark_bar()
    correct = [b_correct]
    abstain = [b_abstain]
    wrong = [b_wrong]
    for key, _ in ORDER:
        row = data["agent_pooled"][key]
        c = row["correct"] / row["n"]
        a = row["abstained"] / row["n"]
        w = (row["emitted"] - row["correct"]) / row["n"]
        if abs(c + a + w - 1.0) > 1e-12:
            raise AssertionError(f"categories do not sum to one for {key}")
        correct.append(c)
        abstain.append(a)
        wrong.append(w)
    return correct, abstain, wrong


def self_test(data):
    correct, abstain, wrong = values(data)
    expected = [0.9644, 0.6823, 0.5344, 0.6062, 0.3717]
    for got, want in zip(correct, expected):
        if round(got, 4) != want:
            raise AssertionError(f"correctness {got:.4f} != {want:.4f}")
    expected_abstain = [0.0, 0.1354, 0.2656, 0.1371, 0.3723]
    for got, want in zip(abstain, expected_abstain):
        if round(got, 4) != want:
            raise AssertionError(f"abstention {got:.4f} != {want:.4f}")
    if any(v < 0 for v in wrong):
        raise AssertionError("negative wrong category")
    print("SELF-TEST PASSED: baseline Figure 7 values match manuscript text")


def draw(data):
    correct, abstain, wrong = values(data)
    labels = ["Curated\nbenchmark"] + [label for _, label in ORDER]
    x = np.arange(len(labels))
    colors = {"correct": "#2A7B7B", "abstain": "#E0A96D", "wrong": "#D1495B"}

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 13, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.bar(x, np.array(correct) * 100, 0.62, color=colors["correct"],
           edgecolor="white", linewidth=0.7, label="Correct (matches CPIC)")
    ax.bar(x, np.array(abstain) * 100, 0.62, bottom=np.array(correct) * 100,
           color=colors["abstain"], edgecolor="white", linewidth=0.7,
           label="Abstained")
    ax.bar(x, np.array(wrong) * 100, 0.62,
           bottom=(np.array(correct) + np.array(abstain)) * 100,
           color=colors["wrong"], edgecolor="white", linewidth=0.7,
           label="Clinically wrong")
    for i, value in enumerate(correct):
        ax.text(i, value * 50, f"{100 * value:.1f}%", ha="center", va="center",
                color="white", fontweight="bold", fontsize=9.5)
    ax.set_title("Curated accuracy does not transfer to real genomes", loc="left", pad=10)
    ax.set_ylabel("Agent responses (%)")
    ax.set_ylim(0, 100)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.yaxis.grid(True, color="#ECECEC")
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.24))
    fig.tight_layout()
    fig.savefig(FIGURES / "Figure7_real_genome.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "Figure7_real_genome.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    print("wrote Figure7_real_genome.png and Figure7_real_genome.tiff")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    data = json.loads(DATA.read_text())
    self_test(data)
    if not args.self_test:
        draw(data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
