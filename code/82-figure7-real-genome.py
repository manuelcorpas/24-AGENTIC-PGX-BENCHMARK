#!/usr/bin/env python3
"""Render the model-level real-genome figure from manuscript source data."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data" / "v4_flash_agent_vs_deterministic.json"
FIGURES = BASE / "figures"

ORDER = [
    ("CorpasFamily", "Family"),
    ("1000G_IBS", "IBS"),
    ("Peru", "Peru"),
    ("UGR", "Uganda"),
]

MODELS = [
    "Claude Opus 4.5",
    "Claude Sonnet 4.5",
    "GPT-5.2",
    "GPT-4.1",
    "o3",
    "o4-mini",
    "Gemini 2.5 Flash",
    "DeepSeek V4 Flash",
]

REAL_MODEL_KEYS = {model: model for model in MODELS}
BENCHMARK_MODEL_KEYS = {model: model for model in MODELS}
POOLED_MODELS = MODELS

COLUMN_KEYS = ["benchmark"] + [key for key, _ in ORDER]
COLUMN_LABELS = ["Curated\ncases"] + [label for _, label in ORDER]


SCORED = BASE / "data" / "v4_flash_matched_scored_rows_all5.json"
# The benchmark bar is authored-rule generation: the cell where the model holds
# the rule table and produces the phenotype itself, which is the mechanism the
# cohort bars test. The execution cells are the opposite handoff (model supplies
# the diplotype, code maps it), so they are not the analogue.
BENCHMARK_CELL = "skill_generation"
BENCHMARK_ATTEMPTED_PER_MODEL = 330


def benchmark_by_model():
    """Return four-way curated outcomes for the eight pinned models."""
    rows = [
        row
        for row in json.loads(SCORED.read_text())
        if row["cell"] == BENCHMARK_CELL
    ]
    source_models = {
        source for source in BENCHMARK_MODEL_KEYS.values() if source is not None
    }
    grouped = {model: [] for model in source_models}
    for row in rows:
        if row["model"] in grouped:
            grouped[row["model"]].append(row)

    result = {}
    for model in MODELS:
        source_model = BENCHMARK_MODEL_KEYS[model]
        if source_model is None:
            result[model] = (np.nan, np.nan, np.nan, np.nan)
            continue
        model_rows = grouped[source_model]
        if len(model_rows) != BENCHMARK_ATTEMPTED_PER_MODEL:
            raise AssertionError(
                f"{model} has {len(model_rows)} benchmark rows, "
                f"expected {BENCHMARK_ATTEMPTED_PER_MODEL}"
            )
        n = BENCHMARK_ATTEMPTED_PER_MODEL
        correct = sum(float(row.get("a1_phenotype") or 0) for row in model_rows)
        invalid = sum(not bool(row.get("parsed_ok")) for row in model_rows)
        abstained = sum(
            bool(row.get("parsed_ok")) and bool(row.get("abstained"))
            for row in model_rows
        )
        wrong = n - correct - abstained - invalid
        if wrong < 0:
            raise AssertionError(f"negative wrong category for {model}")
        result[model] = tuple(
            100 * value / n for value in (correct, abstained, wrong, invalid)
        )
    return result


def benchmark_bar():
    """Score the seven version-matched benchmark models in four categories."""
    benchmark_sources = {BENCHMARK_MODEL_KEYS[model] for model in POOLED_MODELS}
    rows = [
        row
        for row in json.loads(SCORED.read_text())
        if row["cell"] == BENCHMARK_CELL and row["model"] in benchmark_sources
    ]
    if not rows:
        raise AssertionError(f"no scored rows for {BENCHMARK_CELL}")
    n = len(POOLED_MODELS) * BENCHMARK_ATTEMPTED_PER_MODEL
    if len(rows) != n:
        raise AssertionError(f"matched benchmark has {len(rows)} rows, expected {n}")
    correct = sum(row["a1_phenotype"] for row in rows) / n
    invalid = sum(not bool(row.get("parsed_ok")) for row in rows) / n
    abstain = sum(
        bool(row.get("parsed_ok")) and bool(row.get("abstained"))
        for row in rows
    ) / n
    wrong = 1.0 - correct - abstain - invalid
    if wrong < 0:
        raise AssertionError("negative wrong category on the benchmark bar")
    return correct, abstain, wrong, invalid


def model_level_matrices(data):
    """Build four model by setting matrices for panel A."""
    benchmark = benchmark_by_model()
    correct_rows = []
    abstained_rows = []
    wrong_rows = []
    invalid_rows = []

    deepseek = data["model_provenance"]["DeepSeek V4 Flash"]
    if deepseek["returned_model_id"] != ["deepseek-v4-flash"]:
        raise AssertionError("DeepSeek returned-model provenance has changed")

    for model in MODELS:
        model_correct = [benchmark[model][0]]
        model_abstained = [benchmark[model][1]]
        model_wrong = [benchmark[model][2]]
        model_invalid = [benchmark[model][3]]
        for key, _ in ORDER:
            row = data["agent_by_model"][REAL_MODEL_KEYS[model]][key]
            n = row["n"]
            correct = 100 * row["correct"] / n
            abstained = 100 * row["abstained"] / n
            wrong = 100 * row["clinically_wrong"] / n
            invalid = 100 * row["invalid_output"] / n
            if abs(correct + abstained + wrong + invalid - 100.0) > 1e-10:
                raise AssertionError(f"categories do not sum to 100 for {model}, {key}")
            model_correct.append(correct)
            model_abstained.append(abstained)
            model_wrong.append(wrong)
            model_invalid.append(invalid)
        correct_rows.append(model_correct)
        abstained_rows.append(model_abstained)
        wrong_rows.append(model_wrong)
        invalid_rows.append(model_invalid)

    labels = MODELS
    return (
        labels,
        np.asarray(correct_rows, dtype=float),
        np.asarray(abstained_rows, dtype=float),
        np.asarray(wrong_rows, dtype=float),
        np.asarray(invalid_rows, dtype=float),
    )


def pooled_percentages(data):
    """Return the equal-weight pooled percentages shown in panel B."""
    correct, abstained, wrong, invalid = values(data)
    return tuple(
        np.asarray([100 * value for value in outcome], dtype=float)
        for outcome in (correct, abstained, wrong, invalid)
    )


def values(data):
    b_correct, b_abstain, b_wrong, b_invalid = benchmark_bar()
    correct = [b_correct]
    abstain = [b_abstain]
    wrong = [b_wrong]
    invalid = [b_invalid]
    for key, _ in ORDER:
        rows = [
            data["agent_by_model"][REAL_MODEL_KEYS[model]][key]
            for model in POOLED_MODELS
        ]
        n = sum(row["n"] for row in rows)
        c = sum(row["correct"] for row in rows) / n
        a = sum(row["abstained"] for row in rows) / n
        w = sum(row["clinically_wrong"] for row in rows) / n
        i = sum(row["invalid_output"] for row in rows) / n
        if abs(c + a + w + i - 1.0) > 1e-12:
            raise AssertionError(f"categories do not sum to one for {key}")
        correct.append(c)
        abstain.append(a)
        wrong.append(w)
        invalid.append(i)
    return correct, abstain, wrong, invalid


def self_test(data):
    labels, correct_matrix, abstained_matrix, wrong_matrix, invalid_matrix = model_level_matrices(data)
    if labels != MODELS:
        raise AssertionError("model labels are not in the manuscript order")
    if correct_matrix.shape != (8, 5):
        raise AssertionError(f"unexpected matrix shape: {correct_matrix.shape}")
    total = correct_matrix + abstained_matrix + wrong_matrix + invalid_matrix
    if not np.allclose(total[~np.isnan(total)], 100.0):
        raise AssertionError("model-level categories do not sum to 100")

    pooled_correct, pooled_abstained, pooled_wrong, pooled_invalid = pooled_percentages(data)
    expected = [
        96.5530303030,
        63.9423076923,
        61.1842105263,
        55.3571428571,
        38.3552631579,
    ]
    if not np.allclose(pooled_correct, expected):
        raise AssertionError("pooled correctness values do not reproduce")
    expected_abstain = [
        0.0,
        14.9038461538,
        15.6578947368,
        28.2738095238,
        39.6710526316,
    ]
    if not np.allclose(pooled_abstained, expected_abstain):
        raise AssertionError("pooled abstention values do not reproduce")
    if np.any(wrong_matrix[:, 1:] < 0) or np.any(pooled_wrong < 0):
        raise AssertionError("negative wrong category")
    if not np.allclose(
        pooled_correct + pooled_abstained + pooled_wrong + pooled_invalid,
        100.0,
    ):
        raise AssertionError("pooled categories do not sum to 100")
    print("SELF-TEST PASSED: four-outcome Figure 7 values match source data")


def _heatmap(ax, matrix, title, color):
    cmap = LinearSegmentedColormap.from_list(f"white_to_{title}", ["#FFFFFF", color])
    cmap.set_bad("#E6E6E6")
    ax.imshow(
        np.ma.masked_invalid(matrix),
        cmap=cmap,
        vmin=0,
        vmax=100,
        aspect="auto",
        interpolation="none",
    )
    ax.set_title(title, fontsize=13, fontweight="bold", pad=6)
    ax.set_xticks(np.arange(len(COLUMN_LABELS)))
    ax.set_xticklabels(COLUMN_LABELS, fontsize=9.5)
    ax.tick_params(axis="both", which="both", length=0)
    ax.set_xticks(np.arange(-0.5, len(COLUMN_LABELS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            if np.isnan(value):
                text_color = "#4D4D4D"
                label = "NA"
            else:
                text_color = "white" if value >= 52 else "#202020"
                label = "0" if abs(value) < 0.05 else f"{value:.1f}"
            ax.text(
                col,
                row,
                label,
                ha="center",
                va="center",
                fontsize=10.5,
                fontweight="normal",
                color=text_color,
            )


def _pooled_bars(ax, data, colors):
    correct, abstained, wrong, invalid = pooled_percentages(data)
    positions = np.arange(len(COLUMN_LABELS))
    outcomes = [
        ("Correct", correct, colors["correct"]),
        ("Abstained", abstained, colors["abstain"]),
        ("Clinically wrong", wrong, colors["wrong"]),
        ("Invalid output", invalid, colors["invalid"]),
    ]
    bottoms = np.zeros(len(COLUMN_LABELS), dtype=float)

    for label, values_, color in outcomes:
        ax.bar(
            positions,
            values_,
            width=0.68,
            bottom=bottoms,
            label=label,
            color=color,
            edgecolor="white",
            linewidth=1.2,
        )
        for position, value, bottom in zip(positions, values_, bottoms):
            if value < 1.0:
                continue
            text_color = (
                "white"
                if label in {"Correct", "Clinically wrong", "Invalid output"}
                else "#202020"
            )
            text_y = bottom + value / 2
            if bottom + value > 99.5 and value < 7:
                text_y = 102.2
                text_color = "#202020"
            ax.text(
                position,
                text_y,
                f"{value:.1f}",
                ha="center",
                va="center",
                fontsize=10.5,
                fontweight="bold",
                color=text_color,
            )
        bottoms += values_

    ax.set_title("B  Equal-weight pooled summary", loc="left", fontsize=14, fontweight="bold", pad=18)
    ax.set_ylabel("Pooled evaluations (%)", fontsize=10.5)
    ax.set_xticks(positions)
    ax.set_xticklabels(COLUMN_LABELS, fontsize=10.5)
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.tick_params(axis="both", labelsize=10, length=0)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.67, 1.20),
        ncol=4,
        frameon=False,
        fontsize=9.5,
        handlelength=1.25,
        columnspacing=1.0,
    )
    ax.spines["left"].set_color("#767676")
    ax.spines["bottom"].set_color("#767676")


def make_figure(data):
    labels, correct, abstained, wrong, invalid = model_level_matrices(data)
    colors = {
        "correct": "#2A7B7B",
        "abstain": "#E0A96D",
        "wrong": "#D1495B",
        "invalid": "#5B6573",
    }

    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 12,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(10.5, 9.0), facecolor="white")
    grid = fig.add_gridspec(3, 2, height_ratios=[2.2, 2.2, 1.45])
    axes = [fig.add_subplot(grid[0, 0])]
    axes.append(fig.add_subplot(grid[0, 1], sharey=axes[0]))
    axes.append(fig.add_subplot(grid[1, 0], sharey=axes[0]))
    axes.append(fig.add_subplot(grid[1, 1], sharey=axes[0]))
    _heatmap(axes[0], correct, "Correct", colors["correct"])
    _heatmap(axes[1], abstained, "Abstained", colors["abstain"])
    _heatmap(axes[2], wrong, "Clinically wrong", colors["wrong"])
    _heatmap(axes[3], invalid, "Invalid output", colors["invalid"])

    for index in (0, 2):
        axes[index].set_yticks(np.arange(len(labels)))
        axes[index].set_yticklabels(labels, fontsize=10.5)
    for ax in (axes[1], axes[3]):
        ax.tick_params(labelleft=False)

    pooled_axis = fig.add_subplot(grid[2, :])
    _pooled_bars(pooled_axis, data, colors)

    fig.suptitle(
        "Curated accuracy does not transfer to real genomes",
        x=0.04,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
    )
    fig.text(0.04, 0.925, "A", ha="left", va="center", fontsize=14, fontweight="bold")
    fig.text(
        0.075,
        0.925,
        "Model-level outcomes",
        ha="left",
        va="center",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.99,
        0.015,
        "Panel B pools the seven version-matched models",
        ha="right",
        va="bottom",
        fontsize=9.5,
        color="#4D4D4D",
    )
    fig.subplots_adjust(
        left=0.22,
        right=0.99,
        top=0.875,
        bottom=0.075,
        wspace=0.12,
        hspace=0.48,
    )
    return fig


def draw(data):
    fig = make_figure(data)
    fig.savefig(FIGURES / "Figure7_real_genome.png", dpi=300)
    fig.savefig(FIGURES / "Figure7_real_genome.tiff", dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
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
