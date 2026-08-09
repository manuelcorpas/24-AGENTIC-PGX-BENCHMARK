#!/usr/bin/env python3
"""Regenerate supplementary Figures S2 and S3 from the five-cell scored rows.

Both figures previously derived from the retired three-arm protocol and were
plotted against aggregate A3, which the manuscript retires because it returns
1.0 by definition on the 96 non-lethal cases and therefore measured case
composition rather than model behaviour. Citing them from the main text meant
the paper repudiated a metric and then leaned on two figures built from it, one
of which carries the central safety mechanism.

Nothing about the underlying phenomenon needed A3. Both figures are restated against the lethal-class safety action, restricted to
the 336 lethal-class evaluations per cell where a required action is defined.

The action axis is `lethal_action`, not `a2_recommendation`. They are different
measures and only the first reproduces the manuscript: 137 errors under RAG
generation and 54 under free generation. Plotting a2 instead put one cell in the
"phenotype right, action wrong" quadrant for RAG generation, which would have
contradicted the paper's own central mechanism.

  S2  per lethal-class gene, phenotype accuracy against recommendation
      accuracy, one panel per cell. Below the diagonal is information without
      action: the model identifies the phenotype and does not act on it.
  S3  the same lethal-class cells decomposed into the four
      phenotype-by-recommendation quadrants, so correctness by coincidence
      (right action, wrong phenotype) is separated from correctness proper.

Usage:
    python3 code/61-rescore-matched.py --input data/v3_five_cell_live.json \
        --rows data/v3_matched_scored_rows_all5.json --out /tmp/o.json --report /tmp/r.txt
    python3 code/86-figures-s2-s3-five-cell.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
ROWS = BASE / "data" / "v3_matched_scored_rows_all5.json"
OUTDIR = BASE / "figures"

CELLS = [
    ("free_generation", "Free generation"),
    ("rag_generation", "RAG generation"),
    ("rag_execution", "RAG-assisted call\n+ authored execution"),
    ("skill_generation", "Authored rules,\ngenerated"),
    ("skill_execution", "Authored rules,\nexecuted"),
]
# Blue/orange/grey: the red/green pair used elsewhere fails protan separation.
BLUE, ORANGE, GREY, DARK = "#2E6FB7", "#E8873A", "#9E9E9E", "#333333"


def load():
    rows = json.loads(ROWS.read_text())
    return [r for r in rows if r.get("lethal_action") is not None]


def figure_s2(lethal):
    per = defaultdict(lambda: defaultdict(list))
    for r in lethal:
        per[r["cell"]][r["gene"]].append((r["a1_phenotype"], r["lethal_action"]))

    fig, axes = plt.subplots(1, len(CELLS), figsize=(17.5, 4.0), sharex=True, sharey=True)
    for ax, (key, label) in zip(axes, CELLS):
        ax.plot([0, 1], [0, 1], color=GREY, lw=1, ls="--", zorder=1)
        ax.fill_between([0, 1], [0, 1], [0, 0], color=ORANGE, alpha=0.07, zorder=0)
        for gene, pts in sorted(per[key].items()):
            x = sum(p[0] for p in pts) / len(pts)
            y = sum(p[1] for p in pts) / len(pts)
            below = y < x - 0.02
            ax.scatter(x, y, s=52, zorder=3, color=ORANGE if below else BLUE,
                       edgecolor="white", linewidth=0.8)
            if below:
                ax.annotate(gene, (x, y), fontsize=6.4, color=DARK,
                            xytext=(3, -8), textcoords="offset points")
        ax.set_title(label, fontsize=9.5)
        ax.set_xlim(-0.04, 1.04); ax.set_ylim(-0.04, 1.04)
        ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "50", "100"]); ax.set_yticklabels(["0", "50", "100"])
        ax.tick_params(labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("Correct lethal-class action (%)", fontsize=9)
    fig.text(0.5, 0.015, "Phenotype accuracy (%)", ha="center", fontsize=9)
    fig.suptitle("Figure S2. Information without action, lethal-class loci",
                 fontsize=12.5, fontweight="bold", y=0.99)
    fig.text(0.5, 0.925,
             "Each point is one lethal-class gene. The shaded region below the diagonal is "
             "where the phenotype is identified and the required action is not taken.",
             ha="center", fontsize=8.6, color="#555555")
    fig.tight_layout(rect=[0, 0.04, 1, 0.90])
    for ext in ("png", "tiff"):
        fig.savefig(OUTDIR / f"FigureS2-information-without-action.{ext}", dpi=300,
                    facecolor="white")
    plt.close(fig)


def figure_s3(lethal):
    quad = defaultdict(lambda: [0, 0, 0, 0])
    for r in lethal:
        p_ok = r["a1_phenotype"] >= 0.5
        a_ok = r["lethal_action"] >= 0.5
        idx = 0 if (p_ok and a_ok) else 1 if (not p_ok and a_ok) else 2 if (p_ok and not a_ok) else 3
        quad[r["cell"]][idx] += 1

    names = ["phenotype and action both correct",
             "action correct, phenotype wrong (correctness by coincidence)",
             "phenotype correct, action wrong (information without action)",
             "both wrong"]
    colours = [BLUE, ORANGE, "#7A2306", GREY]

    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    ys = range(len(CELLS))
    for i, (key, label) in enumerate(CELLS):
        counts = quad[key]
        total = sum(counts) or 1
        left = 0.0
        for j, c in enumerate(counts):
            frac = c / total * 100
            ax.barh(i, frac, left=left, color=colours[j], edgecolor="white", height=0.62)
            if frac > 4:
                ax.text(left + frac / 2, i, f"{c}", va="center", ha="center",
                        fontsize=8, color="white", fontweight="bold")
            left += frac
    ax.set_yticks(list(ys))
    ax.set_yticklabels([l.replace("\n", " ") for _, l in CELLS], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Percentage of the 336 lethal-class evaluations per cell", fontsize=9)
    ax.set_xlim(0, 100)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colours]
    ax.legend(handles, names, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.20), ncol=2, frameon=False)
    ax.set_title("Figure S3. Correctness by coincidence, lethal-class evaluations",
                 fontsize=12.5, fontweight="bold", loc="left")
    fig.tight_layout()
    for ext in ("png", "tiff"):
        fig.savefig(OUTDIR / f"FigureS3-correctness-by-coincidence.{ext}", dpi=300,
                    facecolor="white")
    plt.close(fig)
    return quad


def main() -> int:
    if not ROWS.exists():
        raise SystemExit(f"missing {ROWS}; run 61-rescore-matched.py --rows first")
    lethal = load()
    expected = 336 * len(CELLS)
    print(f"lethal-class rows: {len(lethal)} (expected about {expected})")
    figure_s2(lethal)
    quad = figure_s3(lethal)
    print("\nquadrant counts per cell (both-correct, coincidence, no-action, both-wrong):")
    for key, label in CELLS:
        c = quad[key]
        print(f"  {label.replace(chr(10), ' '):<34}{c}  n={sum(c)}")
    print(f"\nwrote {OUTDIR}/FigureS2-information-without-action.png")
    print(f"wrote {OUTDIR}/FigureS3-correctness-by-coincidence.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
