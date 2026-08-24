#!/usr/bin/env python3
"""Draw main-manuscript Figure 4 from the frozen eight-model analysis.

Supersedes 78-figure-input-normalisation.py, which reads the seven-model
freeze and asserts o3 has 150 paired units. That figure was last regenerated on
1 August and was still in the manuscript on 9 August under a caption reading
"Paired n = 527 for every model", i.e. a seven-model figure beneath an
eight-model caption. Same defect class as the Figure S3 quadrant error: the
artefact and its legend asserting different things.

The redesigned figure uses two vertically stacked horizontal dumbbell plots.
Panel A shows the coverage gain produced by supplying validated allele
definitions. Panel B separates the remaining model-dependent accuracy from the
PyPGx caller's concordance on the same model-emitted pairs. Both panels use the
same paired variant-call units and the same model order. Provider errors and
output-budget truncations remain visible in the operational counts but leave
behavioural denominators. Error bars are percentile 95% intervals from the
frozen sample-cluster bootstrap. No model calls are issued here.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
FIGDIR = BASE / "figures"
FREEZE = DATA / "v3_input_normalisation_eight_model_freeze.json"

WITHOUT = "#3775BA"
WITH = "#D95F02"
CALLER = "#4D4D4D"
INK = "#111111"
CONNECTOR = "#C7C7C7"
GRID = "#E6E6E6"
MODEL_LABELS = {
    "Claude Opus 4.5": "Claude\nOpus 4.5",
    "Claude Sonnet 4.5": "Claude\nSonnet 4.5",
    "GPT-5.2": "GPT-5.2",
    "GPT-4.1": "GPT-4.1",
    "o3": "o3",
    "o4-mini": "o4-mini",
    "DeepSeek V3": "DeepSeek\nV3",
    "Gemini 2.5 Flash": "Gemini\n2.5 Flash",
}
MIN_CALLS_FOR_ACCURACY = 10
FIGURE_SIZE = (6.3, 7.0)
LEGEND_FONT_SIZE = 9
ANNOTATION_FONT_SIZE = 9
MODEL_FONT_SIZE = 9
PANEL_FONT_SIZE = 11
TITLE_FONT_SIZE = 14
TITLE = (
    "Definitions unlock calls; accuracy remains model-dependent"
)
PANEL_A_TITLE = "A  Coverage increases for every model"
PANEL_B_TITLE = "B  External-consensus accuracy still spans 46% to 81%"
OUTPUT_STEM = "Figure4-input-normalisation-stacked-v1"


def load_freeze() -> dict:
    if not FREEZE.exists():
        raise FileNotFoundError(
            f"missing {FREEZE}; run code/85-freeze-eight-model-normalisation.py"
        )
    return json.loads(FREEZE.read_text())


def build_plot_rows(data: dict) -> list[dict]:
    """Return one plotting record per model, using only frozen values."""
    rows: list[dict] = []
    for model in data["models"]:
        cell = data["models_analysis"][model]
        without = cell["without_definitions"]
        with_definitions = cell["with_definitions"]
        calls_without = without["operational"]["call"]
        accuracy_without = (
            without["accuracy_vs_getrm_among_calls"]
            if calls_without >= MIN_CALLS_FOR_ACCURACY
            else None
        )
        rows.append(
            {
                "model": model,
                "label": MODEL_LABELS[model].replace("\n", " "),
                "paired_units": cell["paired_units"],
                "calls_without": calls_without,
                "coverage_without": without["coverage"]["estimate"],
                "coverage_without_ci": without["coverage"]["ci95"],
                "coverage_with": with_definitions["coverage"]["estimate"],
                "coverage_with_ci": with_definitions["coverage"]["ci95"],
                "coverage_delta": (
                    with_definitions["coverage"]["estimate"]
                    - without["coverage"]["estimate"]
                ),
                "accuracy_without": (
                    accuracy_without["estimate"] if accuracy_without else None
                ),
                "accuracy_without_ci": (
                    accuracy_without["ci95"] if accuracy_without else None
                ),
                "accuracy_with": with_definitions[
                    "accuracy_vs_getrm_among_calls"
                ]["estimate"],
                "accuracy_with_ci": with_definitions[
                    "accuracy_vs_getrm_among_calls"
                ]["ci95"],
                "caller_accuracy": with_definitions[
                    "caller_vs_getrm_on_model_calls"
                ]["estimate"],
                "caller_accuracy_ci": with_definitions[
                    "caller_vs_getrm_on_model_calls"
                ]["ci95"],
            }
        )
    return sorted(rows, key=lambda row: row["accuracy_with"], reverse=True)


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": MODEL_FONT_SIZE,
            "axes.linewidth": 1.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def horizontal_interval(
    ax,
    *,
    estimate: float,
    interval: list[float],
    y: float,
    colour: str,
    marker: str,
    filled: bool,
    zorder: int = 3,
) -> None:
    estimate_pct = 100 * estimate
    errors = [[100 * (estimate - interval[0])], [100 * (interval[1] - estimate)]]
    ax.errorbar(
        estimate_pct,
        y,
        xerr=errors,
        fmt=marker,
        markersize=6.2,
        markerfacecolor=colour if filled else "white",
        markeredgecolor=colour,
        markeredgewidth=1.2,
        ecolor=colour,
        elinewidth=1.2,
        capsize=2.5,
        capthick=1.2,
        zorder=zorder,
    )


def style_axis(ax, rows: list[dict], title: str) -> None:
    positions = list(range(len(rows)))
    ax.set_yticks(positions, [row["label"] for row in rows])
    ax.tick_params(axis="y", labelsize=MODEL_FONT_SIZE, length=0, pad=7)
    ax.tick_params(axis="x", labelsize=MODEL_FONT_SIZE)
    ax.set_xlim(-3, 113)
    ax.set_xticks([0, 25, 50, 75, 100], ["0", "25", "50", "75", "100%"])
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    ax.set_title(title, loc="left", fontsize=PANEL_FONT_SIZE, fontweight="bold", pad=9)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#666666")


def build_figure(data: dict):
    """Build the publication figure and return ``(figure, axes)`` for QA."""
    apply_style()
    rows = build_plot_rows(data)
    figure, axes = plt.subplots(
        2,
        1,
        figsize=FIGURE_SIZE,
        gridspec_kw={"height_ratios": [1, 1.18], "hspace": 0.39},
    )
    figure.suptitle(TITLE, fontsize=TITLE_FONT_SIZE, fontweight="bold", y=0.985)
    y_positions = list(range(len(rows)))

    coverage_ax = axes[0]
    style_axis(coverage_ax, rows, PANEL_A_TITLE)
    for y, row in zip(y_positions, rows):
        coverage_ax.hlines(
            y,
            100 * row["coverage_without"],
            100 * row["coverage_with"],
            color=CONNECTOR,
            linewidth=2.0,
            zorder=1,
        )
        horizontal_interval(
            coverage_ax,
            estimate=row["coverage_without"],
            interval=row["coverage_without_ci"],
            y=y,
            colour=WITHOUT,
            marker="o",
            filled=False,
        )
        horizontal_interval(
            coverage_ax,
            estimate=row["coverage_with"],
            interval=row["coverage_with_ci"],
            y=y,
            colour=WITH,
            marker="o",
            filled=True,
        )
        coverage_ax.text(
            102,
            y,
            f"+{100 * row['coverage_delta']:.0f} pp",
            va="center",
            ha="left",
            fontsize=ANNOTATION_FONT_SIZE,
            color=WITH,
            fontweight="bold",
        )
    coverage_ax.text(
        102,
        -0.78,
        "gain",
        va="center",
        ha="left",
        fontsize=ANNOTATION_FONT_SIZE,
        color="#555555",
    )
    coverage_legend = [
        Line2D(
            [0], [0], marker="o", color="none", markeredgecolor=WITHOUT,
            markerfacecolor="white", markeredgewidth=1.2, markersize=6.2,
            label="No definitions",
        ),
        Line2D(
            [0], [0], marker="o", color="none", markeredgecolor=WITH,
            markerfacecolor=WITH, markeredgewidth=1.2, markersize=6.2,
            label="Validated definitions supplied",
        ),
    ]
    accuracy_ax = axes[1]
    style_axis(accuracy_ax, rows, PANEL_B_TITLE)
    without_offset = 0.21
    caller_offset = -0.21
    for y, row in zip(y_positions, rows):
        accuracy_ax.plot(
            [100 * row["accuracy_with"], 100 * row["caller_accuracy"]],
            [y, y + caller_offset],
            color=CONNECTOR,
            linewidth=1.2,
            zorder=1,
        )
        if row["accuracy_without"] is None:
            accuracy_ax.text(
                1.5,
                y + without_offset,
                f"n={row['calls_without']}",
                va="center",
                ha="left",
                fontsize=ANNOTATION_FONT_SIZE,
                color="#666666",
            )
        else:
            horizontal_interval(
                accuracy_ax,
                estimate=row["accuracy_without"],
                interval=row["accuracy_without_ci"],
                y=y + without_offset,
                colour=WITHOUT,
                marker="o",
                filled=False,
            )
        horizontal_interval(
            accuracy_ax,
            estimate=row["accuracy_with"],
            interval=row["accuracy_with_ci"],
            y=y,
            colour=WITH,
            marker="o",
            filled=True,
        )
        horizontal_interval(
            accuracy_ax,
            estimate=row["caller_accuracy"],
            interval=row["caller_accuracy_ci"],
            y=y + caller_offset,
            colour=CALLER,
            marker="D",
            filled=False,
            zorder=4,
        )

    accuracy_legend = coverage_legend + [
        Line2D(
            [0], [0], marker="D", color="none", markeredgecolor=CALLER,
            markerfacecolor="white", markeredgewidth=1.2, markersize=5.8,
            label="PyPGx on the same emitted pairs",
        )
    ]
    figure.legend(
        handles=accuracy_legend,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.942),
        ncol=3,
        fontsize=LEGEND_FONT_SIZE,
        handletextpad=0.45,
        columnspacing=1.0,
        borderaxespad=0.2,
    )

    figure.text(
        0.5,
        0.012,
        (
            "Paired variant-call inputs: n=527 per model; identical units for all eight models.\n"
            "Error bars: sample-cluster bootstrap 95% intervals; accuracy omitted below 10 calls."
        ),
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#444444",
    )
    figure.subplots_adjust(left=0.265, right=0.97, top=0.82, bottom=0.105)
    return figure, axes


def draw(data: dict, stem: str = OUTPUT_STEM) -> list[Path]:
    figure, _ = build_figure(data)
    FIGDIR.mkdir(exist_ok=True)
    outputs = [
        FIGDIR / f"{stem}.png",
        FIGDIR / f"{stem}.tiff",
        FIGDIR / f"{stem}.pdf",
        FIGDIR / f"{stem}.svg",
    ]
    figure.savefig(outputs[0], dpi=300, facecolor="white")
    figure.savefig(
        outputs[1],
        dpi=600,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    figure.savefig(outputs[2], facecolor="white")
    figure.savefig(outputs[3], facecolor="white")
    plt.close(figure)
    print("wrote " + ", ".join(str(path) for path in outputs))
    return outputs


def self_test(data: dict) -> None:
    assert len(data["models"]) == 8
    assert data["definition_arm_operational_total"] == {
        "attempted": 4216, "call": 3002, "abstain": 1210,
        "error": 0, "truncated_output": 4, "scored": 4212,
    }
    opus = data["models_analysis"]["Claude Opus 4.5"]["with_definitions_full_527"]
    assert round(opus["accuracy_vs_pypgx_among_calls"]["estimate"], 3) == 0.967
    diff = opus["model_minus_caller_vs_getrm"]
    assert diff["ci95"][0] < 0 < diff["ci95"][1]
    assert data["models_analysis"]["o3"]["paired_units"] == 527
    # Every model is paired on the same units; this is what the caption asserts.
    assert all(m["paired_units"] == 527 for m in data["models_analysis"].values())
    print("self-test OK")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    try:
        data = load_freeze()
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 2
    self_test(data)
    if not args.self_test:
        draw(data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
