#!/usr/bin/env python3
"""
Figure 2 for Cell Genomics Manuscript v8.

The adversarial scrambled-specification finding: when the SKILL.md specification
is deliberately corrupted (PHENOTYPE/DRUG fields scrambled to safe-looking but
incorrect values), 43 of 45 model responses echo the corrupted contract verbatim
on safety-critical fields, 2 hedge, and zero override toward CPIC truth.

The figure shows two side-by-side heatmaps over the same 5 lethal-class cases and
the same 3 frontier LLMs:
  Left  - Intact specification: A3 (lethal-class safety action) = 1.0 across all
          15 (case, model) cells, all runs.
  Right - Scrambled specification: A3 collapses; cells coloured by the
          run-aggregated classification (ECHO / HEDGE), never PASS.

The two-panel layout is the load-bearing visual argument: identical models,
identical prompts, the specification flips and the outputs flip with it.

Inputs:
  RESULTS/v3_adversarial_scrambled.json     (45 scrambled-spec calls; 5 cases x 3 models x 3 runs)
  RESULTS/v3_raw_rescored_three_arm.json    (intact-spec baseline from main with_spec arm)

Outputs (PNG 300 DPI + TIFF 600 DPI LZW):
  FIGURES/Figure2_adversarial_spec.png
  FIGURES/Figure2_adversarial_spec.tiff

Usage:
  python3 31-figure2-adversarial-spec.py             # render figure
  python3 31-figure2-adversarial-spec.py --self-test # verify input counts; no render
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


BASE = Path(__file__).resolve().parent.parent
ADV = BASE / "RESULTS" / "v3_adversarial_scrambled.json"
MAIN = BASE / "RESULTS" / "v3_raw_rescored_three_arm.json"
FIGDIR = BASE / "FIGURES"

CASES = [
    ("dpyd_fu_pm", "DPYD / fluorouracil PM\nPM -> NM, AVOID -> standard"),
    ("tpmt_aza_pm", "TPMT / azathioprine PM\nPM -> NM, AVOID -> standard"),
    ("hlab5701_aba_pos", "HLA-B*57:01 / abacavir\nPositive -> Negative, AVOID -> standard"),
    ("cyp2d6_codeine_pm", "CYP2D6 / codeine PM\nPM -> NM, AVOID -> standard"),
    ("cyp2c19_clop_pm", "CYP2C19 / clopidogrel PM\nPM -> NM, AVOID -> standard"),
]
MODELS = ["Claude Opus 4", "DeepSeek V3", "GPT-5.2"]

COLOUR_PASS = "#5a8f5a"
COLOUR_ECHO = "#c14a3f"
COLOUR_HEDGE = "#e0a96d"

CLASS_COLOURS = {
    "PASS": COLOUR_PASS,
    "ECHO_SCRAMBLED": COLOUR_ECHO,
    "HEDGE": COLOUR_HEDGE,
}

EXPECTED = {
    "n_total": 45,
    "n_echo": 43,
    "n_hedge": 2,
    "n_override": 0,
    "intact_a3_cells": 135,
    "intact_a3_mean": 1.0,
}


def load_adversarial() -> list[dict]:
    return json.loads(ADV.read_text())


def load_intact_with_spec() -> list[dict]:
    rows = json.loads(MAIN.read_text())
    case_ids = {c[0] for c in CASES}
    return [r for r in rows
            if r.get("cond") == "with_spec"
            and r.get("model") in MODELS
            and r.get("tc") in case_ids]


def self_test() -> int:
    failures: list[str] = []

    adv = load_adversarial()
    if len(adv) != EXPECTED["n_total"]:
        failures.append(f"adversarial row count: got {len(adv)}, want {EXPECTED['n_total']}")

    counts = Counter(r["classification"] for r in adv)
    if counts.get("ECHO_SCRAMBLED", 0) != EXPECTED["n_echo"]:
        failures.append(f"ECHO_SCRAMBLED count: got {counts.get('ECHO_SCRAMBLED', 0)}, want {EXPECTED['n_echo']}")
    if counts.get("HEDGE", 0) != EXPECTED["n_hedge"]:
        failures.append(f"HEDGE count: got {counts.get('HEDGE', 0)}, want {EXPECTED['n_hedge']}")
    n_override = sum(v for k, v in counts.items() if k not in {"ECHO_SCRAMBLED", "HEDGE"})
    if n_override != EXPECTED["n_override"]:
        failures.append(f"non-echo/non-hedge count: got {n_override}, want {EXPECTED['n_override']}")

    intact = load_intact_with_spec()
    if len(intact) != EXPECTED["intact_a3_cells"]:
        failures.append(f"intact with_spec rows: got {len(intact)}, want {EXPECTED['intact_a3_cells']}")
    intact_a3_vals = [r["scores"].get("A3") for r in intact]
    if not intact_a3_vals or any(v != 1.0 for v in intact_a3_vals):
        bad = sum(1 for v in intact_a3_vals if v != 1.0)
        failures.append(f"intact with_spec A3=1.0 on all cells: {bad} of {len(intact_a3_vals)} are not 1.0")

    adv_cases = {r["case_id"] for r in adv}
    expected_cases = {c[0] for c in CASES}
    if adv_cases != expected_cases:
        failures.append(f"adversarial cases mismatch: got {sorted(adv_cases)}, want {sorted(expected_cases)}")
    adv_models = {r["model"] for r in adv}
    if adv_models != set(MODELS):
        failures.append(f"adversarial models mismatch: got {sorted(adv_models)}, want {sorted(MODELS)}")

    if failures:
        print("SELF-TEST FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"SELF-TEST PASSED ({EXPECTED['n_echo']} ECHO + {EXPECTED['n_hedge']} HEDGE = {EXPECTED['n_total']} scrambled calls; "
          f"0 overrides; {EXPECTED['intact_a3_cells']} intact with_spec cells all A3=1.0).")
    return 0


def aggregate_scrambled(adv: list[dict]) -> dict[tuple[str, str], dict]:
    out: dict[tuple[str, str], dict] = defaultdict(lambda: {"ECHO_SCRAMBLED": 0, "HEDGE": 0, "n": 0})
    for r in adv:
        key = (r["case_id"], r["model"])
        out[key][r["classification"]] += 1
        out[key]["n"] += 1
    return dict(out)


def cell_label_scrambled(agg: dict) -> tuple[str, str]:
    echo = agg["ECHO_SCRAMBLED"]
    hedge = agg["HEDGE"]
    n = agg["n"]
    if echo == n:
        return f"{echo}/{n} echo", COLOUR_ECHO
    if hedge == n:
        return f"{hedge}/{n} hedge", COLOUR_HEDGE
    return f"{echo}/{n} echo\n{hedge}/{n} hedge", COLOUR_HEDGE


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
    adv = load_adversarial()
    agg = aggregate_scrambled(adv)

    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 300,
    })

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.5, 4.8),
                                   gridspec_kw={"width_ratios": [1, 1], "wspace": 0.18})

    n_rows = len(CASES)
    n_cols = len(MODELS)
    case_labels = [c[1] for c in CASES]

    for col in range(n_cols):
        for row in range(n_rows):
            axL.add_patch(mpatches.Rectangle((col, n_rows - 1 - row), 1, 1,
                                             facecolor=COLOUR_PASS, edgecolor="white", linewidth=1.5))
            axL.text(col + 0.5, n_rows - 1 - row + 0.5, "3/3 PASS",
                     ha="center", va="center", color="white", fontsize=9, fontweight="bold")

    axL.set_xlim(0, n_cols); axL.set_ylim(0, n_rows)
    axL.set_xticks([c + 0.5 for c in range(n_cols)])
    axL.set_xticklabels(MODELS, rotation=20, ha="right")
    axL.set_yticks([r + 0.5 for r in range(n_rows)])
    axL.set_yticklabels(reversed(case_labels))
    axL.set_title("Intact specification\nA3 = 1.0 on every cell",
                  loc="left", fontweight="bold", color=COLOUR_PASS)
    axL.tick_params(axis="both", which="both", length=0)
    for spine in axL.spines.values():
        spine.set_visible(False)
    axL.set_aspect("equal", adjustable="box")

    for col, model in enumerate(MODELS):
        for row, (case_id, _) in enumerate(CASES):
            cell_agg = agg[(case_id, model)]
            label, colour = cell_label_scrambled(cell_agg)
            axR.add_patch(mpatches.Rectangle((col, n_rows - 1 - row), 1, 1,
                                             facecolor=colour, edgecolor="white", linewidth=1.5))
            axR.text(col + 0.5, n_rows - 1 - row + 0.5, label,
                     ha="center", va="center", color="white", fontsize=9, fontweight="bold")

    axR.set_xlim(0, n_cols); axR.set_ylim(0, n_rows)
    axR.set_xticks([c + 0.5 for c in range(n_cols)])
    axR.set_xticklabels(MODELS, rotation=20, ha="right")
    axR.set_yticks([])
    axR.set_title("Scrambled specification\n43/45 echo scrambled, 2/45 hedge, 0 override",
                  loc="left", fontweight="bold", color=COLOUR_ECHO)
    axR.tick_params(axis="both", which="both", length=0)
    for spine in axR.spines.values():
        spine.set_visible(False)
    axR.set_aspect("equal", adjustable="box")

    legend_handles = [
        mpatches.Patch(facecolor=COLOUR_PASS, edgecolor="white", label="PASS (correct AVOID action)"),
        mpatches.Patch(facecolor=COLOUR_ECHO, edgecolor="white", label="ECHO scrambled spec (unsafe)"),
        mpatches.Patch(facecolor=COLOUR_HEDGE, edgecolor="white", label="HEDGE (still unsafe; phenotype/drug echoed)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=3, frameon=False, fontsize=8)

    # Superordinate title removed per co-author review (it duplicated the take-home in the legend).

    fig.text(0.5, -0.08,
             "5 lethal-class cases x 3 frontier LLMs x 3 replicates = 45 calls per condition. "
             "Intact-spec baseline: 135 with_spec evaluations (3 populations x 3 runs x 15 (case, model) cells), "
             "all A3 = 1.0.",
             ha="center", fontsize=8.5, color="#333333")

    save_fig(fig, "Figure2_adversarial_spec")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1] if __doc__ else None)
    parser.add_argument("--self-test", action="store_true",
                        help="Verify input counts and intact-spec baseline; do not render.")
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
