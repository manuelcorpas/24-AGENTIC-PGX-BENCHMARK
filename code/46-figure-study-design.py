#!/usr/bin/env python3
"""
Study-design schematic (Figure 1). Three panels:
  A. The matched five-cell comparison: knowledge representation and decision mechanism.
     Two of the six combinations have no knowledge to execute, so five cells
     are defined and one is drawn explicitly as undefined.
  B. The evaluation grid, identical for every cell
     (110 cases x 8 models x 3 replicates = 2,640 per cell, 13,200 attempted;
     13,199 records returned and one absent record counted as failure).
  C. The two arms: curated analytical benchmark + real-genome validation.

The mechanism colours match Figure 2 (red = the model generates the answer,
green = code executes the rule) so the design figure and the result figure
read as one system.

Output: figures/Figure1_study_design.png (300 dpi) + .tiff (600 dpi).
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = Path(__file__).resolve().parent.parent
PNG = BASE / "figures" / "Figure1_study_design.png"
TIFF = BASE / "figures" / "Figure1_study_design.tiff"

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.titlesize": 10,
})

GEN = "#d4756f"      # the model generates the answer   (matches Figure 2)
EXE = "#5a8f5a"      # code executes the rule           (matches Figure 2)
UNDEF = "#e4e4e4"
INK = "#222222"


def box(ax, x, y, w, h, fc, text, ec="black", lw=0.7, fs=8.0, tc="white",
        weight="normal", round=0.04, hatch=None, ls="-"):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0.0,rounding_size={round}",
        linewidth=lw, edgecolor=ec, facecolor=fc, mutation_aspect=1,
        hatch=hatch, linestyle=ls))
    if text:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                color=tc, weight=weight, linespacing=1.25, zorder=5)


def arrow(ax, x1, y1, x2, y2, color=INK, lw=1.4, style="-|>", ms=9):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=ms, lw=lw, color=color,
                                 shrinkA=0, shrinkB=0))


fig, ax = plt.subplots(figsize=(7.4, 8.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")


def ptitle(y, t):
    ax.text(0, y, t, ha="left", va="bottom", fontsize=10, weight="bold", color=INK)


# ---------------- Panel A: the matched factorial ----------------
ptitle(97.3, "A   Matched five-cell comparison")

KNOW = [
    ("none", "no knowledge\nsupplied"),
    ("retrieved prose", "CPIC guideline text\nretrieved into context"),
    ("validated rules", "curated CPIC decision\ntable, versioned"),
]
COL_X0, COL_W, COL_GAP = 20.0, 25.0, 2.5
ROWS = [
    ("the model\ngenerates\nthe answer", 75.5, GEN),
    ("code\nexecutes\nthe rule", 63.0, EXE),
]
ROW_H = 10.0

# column headers. The subtitle sits clear of the box tops: at two lines it
# occupies roughly 87.6 to 90.4, and the first row of boxes tops out at 85.5.
for j, (name, sub) in enumerate(KNOW):
    cx = COL_X0 + j * (COL_W + COL_GAP)
    ax.text(cx + COL_W / 2, 92.2, name, ha="center", va="center",
            fontsize=8.4, weight="bold", color=INK)
    ax.text(cx + COL_W / 2, 89.0, sub, ha="center", va="center",
            fontsize=6.0, color="#6b6b6b", linespacing=1.2)
ax.text(COL_X0 + (3 * COL_W + 2 * COL_GAP) / 2, 94.8,
        "knowledge representation", ha="center", va="center",
        fontsize=7.0, style="italic", color="#8a8a8a")
# Sits above the row labels rather than rotated beside them, where it collided.
ax.text(9.5, 92.2, "decision\nmechanism", ha="center", va="center",
        fontsize=7.0, style="italic", color="#8a8a8a", linespacing=1.2)

CELLS = {
    (0, 0): ("Free\ngeneration", GEN),
    (0, 1): ("RAG\ngeneration", GEN),
    (0, 2): ("Authored rules,\ngenerated", GEN),
    (1, 0): (None, None),                      # nothing to execute
    (1, 1): ("RAG-assisted call\n+ authored execution", EXE),
    (1, 2): ("Authored rules,\nexecuted", EXE),
}

for i, (rlab, ry, _) in enumerate(ROWS):
    ax.text(17.5, ry + ROW_H / 2, rlab, ha="right", va="center",
            fontsize=7.4, weight="bold", color=INK, linespacing=1.2)
    for j in range(3):
        cx = COL_X0 + j * (COL_W + COL_GAP)
        name, fc = CELLS[(i, j)]
        if name is None:
            box(ax, cx, ry, COL_W, ROW_H, UNDEF, "", ec="#b4b4b4", lw=0.8,
                round=0.6, ls=(0, (3, 2)))
            ax.text(cx + COL_W / 2, ry + ROW_H / 2, "undefined\nno knowledge to execute",
                    ha="center", va="center", fontsize=6.2, color="#8a8a8a",
                    style="italic", linespacing=1.25)
        else:
            box(ax, cx, ry, COL_W, ROW_H, fc, name, fs=7.2, weight="bold",
                round=0.6)

ax.text(50.0, 60.4,
        "Every cell receives identical patient text, names the same target drug, emits one output schema,\n"
        "and is scored by one function that is not told which cell produced the row.",
        ha="center", va="top", fontsize=6.6, color=INK, linespacing=1.35)

# ---------------- Panel B: evaluation grid ----------------
ptitle(52.8, "B   Evaluation grid, identical for every cell")
gx, gy0, gw, gh = 1.5, 30.6, 97, 20.5
ax.add_patch(FancyBboxPatch((gx, gy0), gw, gh, boxstyle="round,pad=0,rounding_size=1.2",
                            facecolor="#fafafa", edgecolor="#cfcfcf", lw=0.8))

facs = [("110", "CPIC Level A\ncases"), ("8", "frontier LLMs\ncommon to all cells"),
        ("3", "replicates")]
cw, ch = 21.0, 7.6
cy = gy0 + gh - 10.0
cx0 = gx + 9
for i, (big, lab) in enumerate(facs):
    cx = cx0 + i * (cw + 6.0)
    box(ax, cx, cy, cw, ch, "white", "", ec="#9bb0c6", lw=1.0, round=0.8)
    ax.text(cx + cw / 2, cy + ch * 0.64, big, ha="center", va="center",
            fontsize=13, weight="bold", color="#3f6695")
    ax.text(cx + cw / 2, cy + ch * 0.20, lab, ha="center", va="center",
            fontsize=5.9, color=INK, linespacing=1.05)
    if i < len(facs) - 1:
        ax.text(cx + cw + 3.0, cy + ch / 2, "x", ha="center", va="center",
                fontsize=12, color="#888888")

ax.text(gx + gw / 2, cy - 2.6,
        "= 2,640 attempted evaluations per cell; 13,200 attempted, 13,199 records returned",
        ha="center", va="center", fontsize=7.4, weight="bold", color=INK)
ax.text(gx + gw / 2, gy0 + 3.2,
        "Scored on phenotype and on the drug-specific recommendation, with coverage, abstention and parse failure\n"
        "reported alongside; lethal-class errors counted directly on 336 lethal-class cells.",
        ha="center", va="center", fontsize=5.9, color="#6b6b6b", linespacing=1.3)

# ---------------- Panel C: two arms ----------------
ptitle(27.6, "C   Two arms: controlled benchmark and real-genome validation")
lx, ly, lw_, lh = 1.5, 1.2, 46, 24.0
ax.add_patch(FancyBboxPatch((lx, ly), lw_, lh, boxstyle="round,pad=0,rounding_size=1.2",
                            facecolor="#f6f8fa", edgecolor="#c4d0db", lw=0.9))
ax.text(lx + lw_ / 2, ly + lh - 3.0, "Curated analytical benchmark", ha="center",
        va="center", fontsize=8.6, weight="bold", color=INK)
ax.text(lx + lw_ / 2, ly + lh - 7.0,
        "known CPIC ground truth; balanced\nlethal-class coverage; the five cells above",
        ha="center", va="center", fontsize=6.6, color=INK, linespacing=1.25)
for i, fc in enumerate([GEN, GEN, GEN, EXE, EXE]):
    ax.add_patch(plt.Circle((lx + 10 + i * 6.5, ly + 10.0), 1.9, color=fc,
                            ec="black", lw=0.4))
ax.text(lx + lw_ / 2, ly + 4.6,
        "compares model generation with authored\nexecution under explicitly matched inputs",
        ha="center", va="center", fontsize=6.6, color=INK, linespacing=1.25)

rx = lx + lw_ + 5
ax.add_patch(FancyBboxPatch((rx, ly), lw_, lh, boxstyle="round,pad=0,rounding_size=1.2",
                            facecolor="#f8f6f4", edgecolor="#dbcdc4", lw=0.9))
ax.text(rx + lw_ / 2, ly + lh - 3.0, "Real-genome validation", ha="center",
        va="center", fontsize=8.6, weight="bold", color=INK)
# The three-line description that used to sit here was removed: at this width it
# ran the full box and butted straight into the cohort row, so the panel read as
# crammed. Nothing checkable is lost. The GeT-RM arm (113 reference samples, 527
# pairs) is stated five times in the manuscript body, the per-cohort counts are
# in the boxes directly below, and "7,240 individuals" appeared nowhere in the
# text and was not a registered claim, so the figure was the only thing
# asserting it.
ax.text(rx + lw_ / 2, ly + lh - 7.0,
        "PyPGx-called diplotypes, GRCh37, four cohorts",
        ha="center", va="center", fontsize=6.6, color=INK, linespacing=1.25)
# Every cohort carries its sample count. Two of the four previously showed an
# ancestry label where the other two showed an n, so the panel could not be read
# as a statement of scale.
cohorts = [("1000G IBS", "EUR, n = 93"), ("Corpas\nfamily (WGS)", "EUR, n = 4"),
           ("Peruvian\nGenome Project", "AMR, n = 736"),
           ("Uganda\nGenome Resource", "AFR, n = 6,407")]
for i, (nm, anc) in enumerate(cohorts):
    cxx = rx + 2.2 + i * 10.8
    box(ax, cxx, ly + 8.4, 9.6, 6.6, "white", nm, ec="#cbb6a8", lw=0.8, fs=5.0,
        tc=INK, round=0.6)
    ax.text(cxx + 4.8, ly + 7.6, anc, ha="center", va="top", fontsize=5.2, color="#8a7a6a")
ax.text(rx + lw_ / 2, ly + 4.0,
        "tests whether curated accuracy transfers,\nand measures vocabulary coverage by cohort",
        ha="center", va="center", fontsize=6.6, color=INK, linespacing=1.25)

arrow(ax, lx + lw_ + 0.5, ly + lh / 2, rx - 0.5, ly + lh / 2, color="#888888", lw=1.3, ms=10)

fig.savefig(PNG, dpi=300, bbox_inches="tight")
fig.savefig(TIFF, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
print("wrote", PNG.name, "and", TIFF.name)
