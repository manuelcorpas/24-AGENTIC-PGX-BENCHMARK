#!/usr/bin/env python3
"""Figure: the actual query sent to the model in each of the five matched cells.

Reads the prompts from data/v3_matched_factorial_prompts.json, which is written
by `60-matched-factorial.py --dry-run`, so the figure cannot drift from the text
that was actually issued. Nothing here is retyped by hand.

The point of the figure is that the scaffolding is identical and exactly one
block varies. That is what "matched" means in the design, and it is the thing a
reader otherwise has to take on trust.

Usage:
    python3 code/60-matched-factorial.py --dry-run --limit 1 \
        --cells free_generation rag_generation rag_execution \
                skill_generation skill_execution
    python3 code/84-figure-five-cell-queries.py
"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

BASE = Path(__file__).resolve().parent.parent
PROMPTS = BASE / "data" / "v3_matched_factorial_prompts.json"
OUTDIR = BASE / "figures"

# Blue/orange/grey. The red/green pair used in Figures 1-7 fails colour-vision
# separation at protan dE 3.5; this palette matches Figure 8's dE 24.7.
SHARED = "#5A5A5A"
KNOWLEDGE = {
    "none": "#9E9E9E",
    "retrieved_prose": "#E8873A",
    "structured_rules": "#2E6FB7",
}
# Row order follows Table 1, which is the headline ordering in the main text.
# The supplement previously used a different order (generation cells grouped
# first), so the same five cells appeared in two sequences across the package.
ROWS = [
    ("free_generation", "1. Free generation", "none", "no external knowledge"),
    ("rag_generation", "2. RAG generation", "retrieved_prose", "retrieved guideline prose"),
    ("rag_execution", "3. RAG-assisted call + authored execution", "retrieved_prose",
     "retrieved prose + authored vocabulary"),
    ("skill_generation", "4. Authored rules, generated", "structured_rules",
     "authored rule table"),
    ("skill_execution", "5. Authored rules, executed", "structured_rules",
     "authored rules + vocabulary"),
]


def split_prompt(p: str) -> dict:
    """Cut a prompt into its four structural blocks."""
    instr, _, rest = p.partition("## Knowledge provided")
    knowledge, _, tail = rest.partition("## Patient")
    patient, _, output = tail.partition("## Output")
    return {
        "instruction": instr.strip(),
        "knowledge": knowledge.strip(),
        "patient": ("## Patient" + patient).strip(),
        "output": ("## Output" + output).strip() if output else "(diplotype only)",
    }


def condense(block: str, max_out: int, width: int) -> str:
    """Wrap, then cap the number of RENDERED lines.

    Capping source lines before wrapping is what made the first version overflow:
    15 source lines became 25-plus rendered lines and spilled into the next
    panel. The cap has to apply after wrapping, which is the only count that
    corresponds to vertical space on the page.
    """
    src = [l for l in block.splitlines() if l.strip()]
    out: list[str] = []
    for i, l in enumerate(src):
        wrapped = textwrap.wrap(l, width) or [""]
        if len(out) + len(wrapped) > max_out - 1:
            out.append(f"... {len(src) - i} further lines ...")
            break
        out.extend(wrapped)
    return "\n".join(out)


def main() -> int:
    if not PROMPTS.exists():
        raise SystemExit(f"missing {PROMPTS}; run 60-matched-factorial.py --dry-run first")
    raw = {r["cell"]: r["prompt"] for r in json.loads(PROMPTS.read_text())}

    fig = plt.figure(figsize=(13.5, 16.5))
    gs = fig.add_gridspec(len(ROWS), 1, hspace=0.22, left=0.035, right=0.975,
                          top=0.935, bottom=0.045)

    for i, (cell, title, ktype, klabel) in enumerate(ROWS):
        parts = split_prompt(raw[cell])
        ax = fig.add_subplot(gs[i, 0])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

        colour = KNOWLEDGE[ktype]
        ax.add_patch(FancyBboxPatch((0.002, 0.02), 0.996, 0.96,
                                    boxstyle="round,pad=0.004",
                                    linewidth=1.1, edgecolor=colour,
                                    facecolor="white", zorder=0))
        ax.text(0.012, 0.925, title, fontsize=11.5, fontweight="bold",
                color=colour, va="top")
        n = len(parts["knowledge"])
        ax.text(0.988, 0.925, f"{klabel}   |   knowledge block {n:,} characters",
                fontsize=8.6, color=colour, va="top", ha="right")

        # Column 1: the shared instruction and patient blocks.
        instr = condense(parts["instruction"], 6, 62)
        ax.text(0.012, 0.79, "INSTRUCTION (shared)", fontsize=7.4,
                color=SHARED, fontweight="bold", va="top")
        ax.text(0.012, 0.735, instr, fontsize=7.3, family="monospace",
                color=SHARED, va="top", linespacing=1.35)
        ax.text(0.012, 0.36, "PATIENT (shared, identical in all five cells)",
                fontsize=7.4, color=SHARED, fontweight="bold", va="top")
        ax.text(0.012, 0.305, condense(parts["patient"], 5, 62),
                fontsize=7.3, family="monospace", color=SHARED, va="top",
                linespacing=1.35)

        # Column 2: the one block that varies.
        ax.add_patch(FancyBboxPatch((0.50, 0.05), 0.487, 0.80,
                                    boxstyle="round,pad=0.004",
                                    linewidth=0, facecolor=colour, alpha=0.09,
                                    zorder=0))
        ax.text(0.512, 0.79, "KNOWLEDGE PROVIDED (the only block that varies)",
                fontsize=7.4, color=colour, fontweight="bold", va="top")
        body = parts["knowledge"] or "(none: answer from your own knowledge)"
        ax.text(0.512, 0.735, condense(body, 14, 64),
                fontsize=7.3, family="monospace", color="#1A1A1A", va="top",
                linespacing=1.35)

    fig.suptitle("The query issued in each of the five matched cells",
                 fontsize=14, fontweight="bold", y=0.975)
    fig.text(0.5, 0.955,
             "Every cell receives the same instruction, the same patient text and the same "
             "output schema. Only the knowledge block differs.",
             fontsize=9.2, ha="center", color="#444444")

    OUTDIR.mkdir(exist_ok=True)
    for ext in ("png", "tiff"):
        out = OUTDIR / f"Figure2-five-cell-queries.{ext}"
        fig.savefig(out, dpi=300, facecolor="white")
        print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
