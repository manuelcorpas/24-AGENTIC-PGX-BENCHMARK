#!/usr/bin/env python3
"""
Study-design schematic (new Figure 1). Three panels:
  A. The five-condition constraint gradient (correctness moves model -> executed skill).
  B. The evaluation grid, identical for all five conditions
     (110 cases x 9 models x 3 ancestry framings x 3 replicates; headline ex-Mistral).
  C. The two arms: curated analytical benchmark + real-genome validation across 3 cohorts.

Output: FIGURES/Figure1_study_design.png (300 dpi) + .tiff (600 dpi).
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = Path(__file__).resolve().parent.parent
PNG = BASE / "FIGURES" / "Figure1_study_design.png"
TIFF = BASE / "FIGURES" / "Figure1_study_design.tiff"

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.titlesize": 10,
})

# condition palette (first three match the existing figures)
C = {
    "free":   "#d4756f",
    "rag":    "#e0a96d",
    "reason": "#7aa0c4",
    "exec":   "#3f6695",
    "ctrl":   "#5a8f5a",
}
INK = "#222222"


def box(ax, x, y, w, h, fc, text, ec="black", lw=0.7, fs=8.0, tc="white", weight="normal", round=0.04):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.0,rounding_size={round}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, mutation_aspect=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc,
            weight=weight, linespacing=1.25, zorder=5)


def arrow(ax, x1, y1, x2, y2, color=INK, lw=1.4, style="-|>", ms=9):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms,
                                 lw=lw, color=color, shrinkA=0, shrinkB=0))


fig, ax = plt.subplots(figsize=(7.4, 8.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")


def ptitle(y, t):
    ax.text(0, y, t, ha="left", va="bottom", fontsize=10, weight="bold", color=INK)


# ---------------- Panel A: constraint gradient ----------------
ptitle(95.5, "A   The five-condition constraint gradient")
conds = [
    (C["free"],   "Free-prompted",     "model reasons\nfrom its prior"),
    (C["rag"],    "Retrieval-\naugmented", "model reasons over\nretrieved CPIC text"),
    (C["reason"], "Skill-\nreasoning",  "model applies the\nskill's rules"),
    (C["exec"],   "Skill-\nexecution",  "model supplies input;\nskill computes in code"),
    (C["ctrl"],   "Answer-\nsupplied",  "deterministic\npositive control"),
]
n = len(conds); x0, gap, bw = 1.5, 1.5, (97 - (len(conds) - 1) * 1.5) / 5
by, bh = 85.5, 7.0
for i, (fc, name, desc) in enumerate(conds):
    x = x0 + i * (bw + gap)
    box(ax, x, by, bw, bh, fc, name, fs=8.2, weight="bold", round=0.6)
    ax.text(x + bw / 2, by - 1.2, desc, ha="center", va="top", fontsize=6.2, color=INK, linespacing=1.15)
# gradient bar beneath, clear of the descriptors above
gy = 79.3
ax.add_patch(FancyBboxPatch((x0, gy - 1.0), 97, 2.0, boxstyle="round,pad=0,rounding_size=1.0",
                            facecolor="#ededed", edgecolor="#bcbcbc", lw=0.6))
arrow(ax, x0 + 2, gy, x0 + 95, gy, color="#555555", lw=1.6, ms=11)
# both label pairs sit BELOW the bar so they never collide with the box descriptors
ax.text(x0 + 1, gy - 2.0, "correctness in the stochastic model\nstochastic, unauditable",
        ha="left", va="top", fontsize=6.3, style="italic", color="#666666", linespacing=1.3)
ax.text(x0 + 96, gy - 2.0, "correctness in the executed skill\ndeterministic, auditable",
        ha="right", va="top", fontsize=6.3, style="italic", color="#666666", linespacing=1.3)

# ---------------- Panel B: evaluation grid ----------------
ptitle(71.6, "B   Evaluation grid: 8,910 evaluations per condition, 44,550 in total")
gx, gy0, gw, gh = 1.5, 45.0, 97, 24.5
ax.add_patch(FancyBboxPatch((gx, gy0), gw, gh, boxstyle="round,pad=0,rounding_size=1.2",
                            facecolor="#fafafa", edgecolor="#cfcfcf", lw=0.8))
# factorial line of chips (how each condition's count is built)
facs = [("110", "CPIC Level A\ncases"), ("9", "frontier LLMs"),
        ("3", "ancestry framings\nEUR / AMR / AFR"), ("3", "replicates")]
cw, ch = 18.5, 7.2; cy = gy0 + gh - 9.0; cx0 = gx + 6
for i, (big, lab) in enumerate(facs):
    cx = cx0 + i * (cw + 4.0)
    box(ax, cx, cy, cw, ch, "white", "", ec="#9bb0c6", lw=1.0, round=0.8)
    ax.text(cx + cw / 2, cy + ch * 0.63, big, ha="center", va="center", fontsize=13, weight="bold", color=C["exec"])
    ax.text(cx + cw / 2, cy + ch * 0.20, lab, ha="center", va="center", fontsize=5.9, color=INK, linespacing=1.05)
    if i < len(facs) - 1:
        ax.text(cx + cw + 2.0, cy + ch / 2, "x", ha="center", va="center", fontsize=12, color="#888888")
ax.text(gx + gw / 2, cy - 1.8, "=  8,910 evaluations per condition  (each scored on phenotype A1, drug recommendation A2, lethal-class safety A3)",
        ha="center", va="center", fontsize=6.6, weight="bold", color=INK)
# per-condition breakdown row (the five conditions of Panel A, each summing to 44,550)
conds_b = [(C["free"], "Free-\nprompt", False), (C["rag"], "Retrieval", False),
           (C["reason"], "Skill-\nreasoning", True), (C["exec"], "Skill-\nexecution", True),
           (C["ctrl"], "Answer-\nsupplied", False)]
bw = (gw - 12) / 5; bx0 = gx + 6; byy = gy0 + 4.4; bhh = 7.0
for i, (fc, nm, skill) in enumerate(conds_b):
    bxx = bx0 + i * bw
    box(ax, bxx + 0.8, byy, bw - 1.6, bhh, fc, "", ec="black", lw=0.5, round=0.45)
    ax.text(bxx + bw / 2, byy + bhh * 0.70, nm + (" *" if skill else ""), ha="center", va="center",
            fontsize=5.8, color="white", weight="bold", linespacing=1.0)
    ax.text(bxx + bw / 2, byy + bhh * 0.24, "8,910", ha="center", va="center", fontsize=8.0, color="white", weight="bold")
    if i < 4:
        ax.text(bxx + bw - 0.4, byy + bhh / 2, "+", ha="center", va="center", fontsize=10, color="#888888")
ax.text(gx + gw / 2, gy0 + 1.4,
        "*Skill conditions exclude Mistral Large 2 from the headline (8 models, 7,920 evaluations; rate-limit artefact); the all-nine-model view is in the supplement.",
        ha="center", va="center", fontsize=5.8, style="italic", color="#7a6a6a")

# ---------------- Panel C: two arms ----------------
ptitle(40.5, "C   Two arms: controlled benchmark and real-genome validation")
# left arm
lx, ly, lw_, lh = 1.5, 3.0, 46, 33.0
ax.add_patch(FancyBboxPatch((lx, ly), lw_, lh, boxstyle="round,pad=0,rounding_size=1.2",
                            facecolor="#f6f8fa", edgecolor="#c4d0db", lw=0.9))
ax.text(lx + lw_ / 2, ly + lh - 3.0, "Curated analytical benchmark", ha="center", va="center",
        fontsize=8.6, weight="bold", color=INK)
ax.text(lx + lw_ / 2, ly + lh - 7.2, "known CPIC ground truth; balanced\nlethal-class coverage; the 5 conditions above",
        ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.25)
# mini gradient of 5 dots
for i, key in enumerate(["free", "rag", "reason", "exec", "ctrl"]):
    ax.add_patch(plt.Circle((lx + 9 + i * 7.0, ly + 16.5), 1.9, color=C[key], ec="black", lw=0.4))
ax.text(lx + lw_ / 2, ly + 9.5, "isolates WHERE correctness must reside\n(curated accuracy reaches ~96%)",
        ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.25)
# right arm
rx = lx + lw_ + 5
ax.add_patch(FancyBboxPatch((rx, ly), lw_, lh, boxstyle="round,pad=0,rounding_size=1.2",
                            facecolor="#f8f6f4", edgecolor="#dbcdc4", lw=0.9))
ax.text(rx + lw_ / 2, ly + lh - 3.0, "Real-genome validation", ha="center", va="center",
        fontsize=8.6, weight="bold", color=INK)
ax.text(rx + lw_ / 2, ly + lh - 7.2, "PyPGx-called diplotypes from three\nancestrally distinct cohorts (GRCh37)",
        ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.25)
cohorts = [("Corpas family", "European"), ("Peruvian Genome\nProject", "Latin American"),
           ("Uganda Genome\nResource", "East African")]
for i, (nm, anc) in enumerate(cohorts):
    cxx = rx + 4 + i * 13.0
    box(ax, cxx, ly + 13.0, 11.5, 7.5, "white", nm, ec="#cbb6a8", lw=0.8, fs=5.6, tc=INK, round=0.6)
    ax.text(cxx + 5.75, ly + 12.2, anc, ha="center", va="top", fontsize=5.6, color="#8a7a6a")
ax.text(rx + lw_ / 2, ly + 6.0, "tests whether curated accuracy transfers\n(it does not: 72% / 51% / 40% by ancestry)",
        ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.25)
# arrow between arms
arrow(ax, lx + lw_ + 0.5, ly + lh / 2, rx - 0.5, ly + lh / 2, color="#888888", lw=1.3, ms=10)

fig.savefig(PNG, dpi=300, bbox_inches="tight")
fig.savefig(TIFF, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
print("wrote", PNG.name, "and", TIFF.name)
