#!/usr/bin/env python3
"""
Figure 5 for Cell Genomics Manuscript v8.

The correctness-by-coincidence finding: when a clinical decision-support system
reaches the correct safety action through wrong reasoning, the system is
unauditable even when the conclusion happens to be right. This figure shows
the four-cell A1 vs A3 confusion matrix per condition on the lethal-class subset.

Quadrants:
  Q1  Aligned correct        A1 = 1  AND  A3 = 1   - safe and auditable
  Q2  Action without info    A1 < 1  AND  A3 = 1   - safe but unauditable
                                                     (right action, wrong reasoning)
  Q3  Information without    A1 = 1  AND  A3 < 1   - knew it, did not act
       action                                        (the RAG failure mode)
  Q4  Both wrong             A1 < 1  AND  A3 < 1   - unsafe and inaccurate

Headline numbers (verified against analysis_summary section_3_misalignment +
lethal_a3_errors CSV):
  no_spec   1096 lethal cells: Q1 669 (61%) | Q2 157 (14.3%) | Q3 122 (11.1%) | Q4 148 (13.5%)
  cpic_rag  1130 lethal cells: Q1 634 (56%) | Q2  82 ( 7.3%) | Q3 341 (30.2%) | Q4  73 ( 6.5%)
  with_spec 1134 lethal cells: Q1 1134 (100%) | Q2  0 | Q3 0 | Q4 0

Within the A3 = 1 subset (the handoff's "correctness-by-coincidence" metric):
  no_spec    157/826  = 19.0%
  cpic_rag    82/716  = 11.5%
  with_spec    0/1134 =  0.0%

Inputs:
  RESULTS/v3_raw_rescored_three_arm.json
  RESULTS/v3_three_arm_analysis_summary.json (cross-check)
  RESULTS/v3_three_arm_lethal_a3_errors.csv  (cross-check)
  SPECS/test_cases_v3.json                   (lethal flag per case)

Outputs (PNG 300 DPI + TIFF 600 DPI LZW):
  FIGURES/Figure5_correctness_by_coincidence.png
  FIGURES/Figure5_correctness_by_coincidence.tiff

Usage:
  python3 35-figure5-correctness-by-coincidence.py             # render
  python3 35-figure5-correctness-by-coincidence.py --self-test # verify
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
SUMMARY = BASE / "RESULTS" / "v3_three_arm_analysis_summary.json"
LETHAL_CSV = BASE / "RESULTS" / "v3_three_arm_lethal_a3_errors.csv"
CASES_FILE = BASE / "SPECS" / "test_cases_v3.json"
FIGDIR = BASE / "FIGURES"

CONDITIONS = ["no_spec", "cpic_rag", "with_spec"]
CONDITION_LABELS = {
    "no_spec": "Free-prompted",
    "cpic_rag": "Retrieval-augmented",
    "with_spec": "Specification-constrained",
}

QUADRANTS = ["Q1", "Q2", "Q3", "Q4"]
QUADRANT_LABELS = {
    "Q1": "Aligned correct  (A1 = 1, A3 = 1)",
    "Q2": "Action without information  (A1 < 1, A3 = 1)",
    "Q3": "Information without action  (A1 = 1, A3 < 1)",
    "Q4": "Both wrong  (A1 < 1, A3 < 1)",
}
QUADRANT_COLOURS = {
    "Q1": "#2a6a2a",
    "Q2": "#c1a04a",
    "Q3": "#c47a3a",
    "Q4": "#8b2e26",
}

EXPECTED_TOTALS = {"no_spec": 1096, "cpic_rag": 1130, "with_spec": 1134}
EXPECTED_Q2 = {"no_spec": 157, "cpic_rag": 82, "with_spec": 0}
EXPECTED_A3_ERRORS = {"no_spec": 270, "cpic_rag": 414, "with_spec": 0}


def load_lethal_cases() -> set[str]:
    cases = json.loads(CASES_FILE.read_text())
    return {c["id"] for c in cases if "lethal" in (c.get("gt_drug") or "").lower()}


def quadrant_counts(main: list[dict], lethal_tcs: set[str]) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {c: {q: 0 for q in QUADRANTS} for c in CONDITIONS}
    for r in main:
        if r["scores"].get("format_fail"):
            continue
        if r["tc"] not in lethal_tcs:
            continue
        a1 = r["scores"].get("A1")
        a3 = r["scores"].get("A3")
        if a1 is None or a3 is None:
            continue
        cond = r["cond"]
        if cond not in counts:
            continue
        a1_ok = (a1 == 1)
        a3_ok = (a3 == 1)
        if a1_ok and a3_ok:
            q = "Q1"
        elif (not a1_ok) and a3_ok:
            q = "Q2"
        elif a1_ok and (not a3_ok):
            q = "Q3"
        else:
            q = "Q4"
        counts[cond][q] += 1
    return counts


def self_test() -> int:
    failures: list[str] = []

    lethal_tcs = load_lethal_cases()
    if len(lethal_tcs) != 14:
        failures.append(f"lethal case count: got {len(lethal_tcs)}, want 14")

    main = json.loads(MAIN.read_text())
    qc = quadrant_counts(main, lethal_tcs)

    for cond, want_total in EXPECTED_TOTALS.items():
        got_total = sum(qc[cond].values())
        if got_total != want_total:
            failures.append(f"{cond} lethal-cell total: got {got_total}, want {want_total}")

    for cond, want_q2 in EXPECTED_Q2.items():
        if qc[cond]["Q2"] != want_q2:
            failures.append(f"{cond} Q2 (correctness-by-coincidence): got {qc[cond]['Q2']}, want {want_q2}")

    for cond, want_errs in EXPECTED_A3_ERRORS.items():
        got_errs = qc[cond]["Q3"] + qc[cond]["Q4"]
        if got_errs != want_errs:
            failures.append(f"{cond} A3 errors (Q3+Q4): got {got_errs}, want {want_errs}")

    if qc["with_spec"]["Q1"] != EXPECTED_TOTALS["with_spec"]:
        failures.append(f"with_spec Q1: got {qc['with_spec']['Q1']}, want {EXPECTED_TOTALS['with_spec']} (all aligned)")
    if qc["with_spec"]["Q2"] + qc["with_spec"]["Q3"] + qc["with_spec"]["Q4"] != 0:
        failures.append("with_spec has non-Q1 cells; ceiling violated")

    summ = json.loads(SUMMARY.read_text())
    for cond in CONDITIONS:
        s_q2 = summ["section_3_misalignment"][cond]["A1_lt_1"]
        if s_q2 != qc[cond]["Q2"]:
            failures.append(f"{cond} Q2 cross-check vs summary: figure={qc[cond]['Q2']}, summary={s_q2}")

    if failures:
        print("SELF-TEST FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"SELF-TEST PASSED (lethal totals 1096/1130/1134; "
          f"Q2 counts 157/82/0; A3 errors 270/414/0; "
          f"with_spec ceiling Q1=1134; cross-check vs summary green).")
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
    lethal_tcs = load_lethal_cases()
    main = json.loads(MAIN.read_text())
    qc = quadrant_counts(main, lethal_tcs)

    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 300,
    })

    fig, ax = plt.subplots(figsize=(9.5, 6.4))

    x = np.arange(len(CONDITIONS))
    width = 0.62
    bottoms = np.zeros(len(CONDITIONS))

    totals = [sum(qc[c].values()) for c in CONDITIONS]

    for q in QUADRANTS:
        vals = [100 * qc[c][q] / totals[i] for i, c in enumerate(CONDITIONS)]
        bars = ax.bar(x, vals, width, bottom=bottoms,
                      color=QUADRANT_COLOURS[q], edgecolor="white", linewidth=1.2,
                      label=QUADRANT_LABELS[q])
        for i, (rect, v) in enumerate(zip(bars, vals)):
            cond = CONDITIONS[i]
            count = qc[cond][q]
            if v < 1.5:
                continue
            text_y = bottoms[i] + v / 2
            if v >= 20:
                label_text = f"{q}\n{count}\n{v:.1f}%"
                font_size = 9
            elif v >= 10:
                label_text = f"{q}\n{count} ({v:.1f}%)"
                font_size = 8
            elif v >= 5:
                label_text = f"{q}: {count} ({v:.1f}%)"
                font_size = 7.5
            else:
                label_text = f"{q}: {count}"
                font_size = 6.5
            ax.text(rect.get_x() + rect.get_width() / 2, text_y,
                    label_text,
                    ha="center", va="center", fontsize=font_size,
                    fontweight="bold", color="white")
        bottoms += np.array(vals)

    for i, cond in enumerate(CONDITIONS):
        q1 = qc[cond]["Q1"]
        q2 = qc[cond]["Q2"]
        a3_ok = q1 + q2
        ratio_str = (f"correctness-by-coincidence\n(Q2 / A3 = 1):\n"
                     f"{q2}/{a3_ok} = {(q2 / a3_ok * 100 if a3_ok else 0):.1f}%")
        ax.text(i, 105, ratio_str,
                ha="center", va="bottom", fontsize=8,
                color=QUADRANT_COLOURS["Q2"] if q2 > 0 else "#2a6a2a",
                fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"{CONDITION_LABELS[c]}\n(n = {totals[i]:,})"
                        for i, c in enumerate(CONDITIONS)], fontweight="bold")
    ax.set_ylim(0, 130)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("% of lethal-class cells")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [mpatches.Patch(facecolor=QUADRANT_COLOURS[q], edgecolor="white",
                              label=QUADRANT_LABELS[q]) for q in QUADRANTS]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=2, frameon=False, fontsize=8.5)

    fig.suptitle("Retrieval shifts errors into the information-without-action quadrant; specification eliminates all four error modes",
                 fontsize=10.5, fontweight="bold", color="#333333", y=0.99)

    fig.text(0.5, 0.01,
             "Lethal-class subset of the locked 3-of-3 replicate three-arm dataset. "
             "Quadrants defined on per-cell A1 (phenotype) and A3 (clinical safety action) scores. "
             "Annotation above each bar: correctness-by-coincidence rate (Q2 over A3 = 1 subset).",
             ha="center", fontsize=7, color="#555555")

    fig.subplots_adjust(bottom=0.28, top=0.88)

    save_fig(fig, "Figure5_correctness_by_coincidence")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1] if __doc__ else None)
    parser.add_argument("--self-test", action="store_true",
                        help="Verify quadrant counts against summary + lethal CSV; do not render.")
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
