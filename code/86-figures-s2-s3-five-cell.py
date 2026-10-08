#!/usr/bin/env python3
"""Supplementary phenotype and wording-rule plots for the five configurations.

S2 plots per-gene phenotype concordance against original wording-rule acceptance
among the 336 high-severity evaluations per configuration. S3 partitions those
scores at 0.5 and shows the joint scoring categories. These lexical endpoints do
not independently establish clinical action, safety or harm. Both wording rules
are reported in main Figure 2; conclusions depending on the strict rule alone
are withdrawn. R3 uses the complete attempted grid from the regeneration driver.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
ROWS = BASE / "data" / "v4_flash_matched_scored_rows_all5.json"
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


def figure_s2(lethal, save_output=True):
    per = defaultdict(lambda: defaultdict(list))
    for r in lethal:
        per[r["cell"]][r["gene"]].append((r["a1_phenotype"], r["lethal_action"]))

    fig, grid = plt.subplots(3, 2, figsize=(7.2, 7.4), sharex=True, sharey=True)
    axes=list(grid.flat);axes[-1].set_visible(False)
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
                ax.annotate(gene, (x, y), fontsize=7.2, color=DARK,
                            xytext=(3, 8) if gene == 'HLA-A*31:01' else (3, -8),
                            textcoords="offset points")
        ax.set_title(label.replace('\n',' '), fontsize=8.5)
        ax.set_xlim(-0.04, 1.04); ax.set_ylim(-0.04, 1.04)
        ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "50", "100"]); ax.set_yticklabels(["0", "50", "100"])
        ax.tick_params(labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.supylabel("Original wording rule accepted (%)", fontsize=9,x=0.02)
    fig.supxlabel("Phenotype concordance (%)", fontsize=9,y=0.015)
    fig.suptitle("Figure S2. Phenotype concordance and original wording rule",
                 fontsize=10, fontweight="bold", y=0.99)
    fig.text(0.5, 0.945,
             "Each point is one high-severity gene. Shading marks lower wording acceptance.\nThese lexical classifications do not establish clinical harm.",
             ha="center", va="top", fontsize=8, color="#555555")
    fig.tight_layout(rect=[0.04, 0.04, 1, 0.925])
    if save_output:
        for ext in ("png", "tiff"):
            fig.savefig(OUTDIR / f"FigureS2-information-without-action.{ext}", dpi=300, facecolor="white")
        plt.close(fig)
    return fig


def figure_s3(lethal, save_output=True):
    quad = defaultdict(lambda: [0, 0, 0, 0])
    for r in lethal:
        p_ok = r["a1_phenotype"] >= 0.5
        a_ok = r["lethal_action"] >= 0.5
        idx = 0 if (p_ok and a_ok) else 1 if (not p_ok and a_ok) else 2 if (p_ok and not a_ok) else 3
        quad[r["cell"]][idx] += 1

    names = ["phenotype score >=0.5; wording accepted",
             "phenotype score <0.5; wording accepted",
             "phenotype score >=0.5; wording rejected",
             "phenotype score <0.5; wording rejected"]
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
    ax.set_xlabel("Percentage of 336 high-severity evaluations per configuration", fontsize=9)
    ax.set_xlim(0, 100)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colours]
    ax.legend(handles, names, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.20), ncol=2, frameon=False)
    ax.set_title("Figure S3. Joint phenotype scores and original wording classifications",
                 fontsize=12.5, fontweight="bold", loc="left")
    fig.tight_layout()
    if not save_output:return fig
    for ext in ("png", "tiff"):
        fig.savefig(OUTDIR / f"FigureS3-correctness-by-coincidence.{ext}", dpi=300, facecolor="white")
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
