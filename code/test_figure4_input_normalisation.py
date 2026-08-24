#!/usr/bin/env python3
"""Regression gates for the redesigned main-manuscript Figure 4.

The submitted landscape figure becomes illegible at manuscript width and asks
the reader to decode three comparisons at once.  These gates require a
portrait, vertically stacked redesign whose two panels share one model order
and one typography system.

Run:  python3 code/test_figure4_input_normalisation.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
GEN = Path(__file__).resolve().parent / "87-figure4-input-normalisation-stacked.py"


def load_generator():
    spec = importlib.util.spec_from_file_location("figure4", GEN)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    failures: list[str] = []

    def check(name: str, condition: bool, detail: str = "") -> None:
        if condition:
            print(f"  pass  {name}")
        else:
            print(f"  FAIL  {name}{('  ' + detail) if detail else ''}")
            failures.append(name)

    fig4 = load_generator()
    data = fig4.load_freeze()

    print("Figure 4 redesign gates")

    check(
        "figure is portrait for two vertically stacked panels",
        fig4.FIGURE_SIZE[1] > fig4.FIGURE_SIZE[0],
        str(fig4.FIGURE_SIZE),
    )
    check(
        "figure uses the manuscript's native 6.3-inch image width",
        fig4.FIGURE_SIZE == (6.3, 7.0),
        str(fig4.FIGURE_SIZE),
    )
    check(
        "legend and annotation text use one uniform size",
        fig4.LEGEND_FONT_SIZE == fig4.ANNOTATION_FONT_SIZE,
        f"legend={fig4.LEGEND_FONT_SIZE}, annotation={fig4.ANNOTATION_FONT_SIZE}",
    )

    rows = fig4.build_plot_rows(data)
    check("all eight models are represented", len(rows) == 8, str(len(rows)))
    check(
        "every model uses the same 527 paired inputs",
        all(row["paired_units"] == 527 for row in rows),
    )
    check(
        "every coverage dumbbell points toward higher coverage",
        all(row["coverage_with"] > row["coverage_without"] for row in rows),
        str([(row["model"], row["coverage_delta"]) for row in rows]),
    )
    check(
        "models are ordered by accuracy with definitions",
        [row["accuracy_with"] for row in rows]
        == sorted((row["accuracy_with"] for row in rows), reverse=True),
    )

    gpt52 = next(row for row in rows if row["model"] == "GPT-5.2")
    check(
        "accuracy without definitions is suppressed below ten calls",
        gpt52["accuracy_without"] is None and gpt52["calls_without"] == 2,
        str((gpt52["accuracy_without"], gpt52["calls_without"])),
    )

    check(
        "accuracy range with definitions is data-derived",
        round(min(row["accuracy_with"] for row in rows), 3) == 0.464
        and round(max(row["accuracy_with"] for row in rows), 3) == 0.814,
    )
    check(
        "title states the necessary-but-not-sufficient message",
        "unlock calls" in fig4.TITLE.lower()
        and "accuracy" in fig4.TITLE.lower()
        and "model-dependent" in fig4.TITLE.lower(),
        fig4.TITLE,
    )

    figure, axes = fig4.build_figure(data)
    check("exactly two data panels are drawn", len(axes) == 2, str(len(axes)))
    if len(axes) == 2:
        check(
            "coverage panel is above accuracy panel",
            axes[0].get_position().y0 > axes[1].get_position().y0,
        )
        check(
            "zero-value markers have visible left padding",
            all(axis.get_xlim()[0] <= -2 for axis in axes),
            str([axis.get_xlim() for axis in axes]),
        )
        check(
            "one shared legend replaces colliding panel legends",
            len(figure.legends) == 1
            and all(axis.get_legend() is None for axis in axes),
            f"figure legends={len(figure.legends)}",
        )
        legend_sizes = {
            text.get_fontsize()
            for legend in figure.legends
            for text in legend.get_texts()
        }
        check(
            "all rendered legend labels use the uniform size",
            legend_sizes == {float(fig4.LEGEND_FONT_SIZE)},
            str(legend_sizes),
        )

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    width, _ = figure.canvas.get_width_height()
    bounded = []
    for text in figure.texts:
        box = text.get_window_extent(renderer=renderer)
        bounded.append(box.x0 >= 0 and box.x1 <= width)
    check("figure-level text stays inside the canvas", all(bounded), str(bounded))
    methods_notes = [text.get_text() for text in figure.texts if "Paired variant-call" in text.get_text()]
    check(
        "methods note is split across lines for manuscript width",
        len(methods_notes) == 1 and "\n" in methods_notes[0],
        str(methods_notes),
    )

    strings = [fig4.TITLE, fig4.PANEL_A_TITLE, fig4.PANEL_B_TITLE]
    check(
        "artwork strings contain no long dashes",
        all(not any(ch in value for ch in "—–") for value in strings),
    )

    import matplotlib.pyplot as plt

    plt.close(figure)
    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("all Figure 4 redesign gates pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
