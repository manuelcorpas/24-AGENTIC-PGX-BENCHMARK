#!/usr/bin/env python3
"""Figure S8: the conditional-invariance control, rebuilt from the deposit.

The submitted FigureS8.tiff had no generator in any commit, so it could not be
checked, corrected or reproduced. This script replaces it. Everything it draws
is computed here from two files in the Zenodo deposit
(10.5281/zenodo.21969051, version 1.7.0); nothing is typed.

Three defects in the submitted figure, each of which this script makes
impossible rather than merely fixed:

  1. Its title read "Population-invariance: the EUR/AMR/AFR 95% CIs overlap;
     the gap vanishes under execution and control". Skill-execution carries the
     LARGEST gap of the five arms, so the title asserted the opposite of its
     own picture, which printed a delta of 1.7 directly above it.
  2. It printed 1.7 for skill-execution while its own bars give 1.63. The bars
     came from v3_armAB_fullgrid.json and the delta from
     v3_armA9_armBv2_POP.json, which is a different subset. Here the delta is
     the arithmetic of the bars and cannot disagree with them.
  3. It said "population-invariance" flat, where the legend and Minor 2 reserve
     invariance for the deterministic mapping conditional on identical valid
     inputs. That condition holds in the control arm and nowhere else on the
     panel.

Provenance of each arm, recovered by reproducing the submitted bar heights to
two decimals before anything was changed:

  Free-prompt       v3_raw_rescored_three_arm.json   cond=no_spec
  RAG               v3_raw_rescored_three_arm.json   cond=cpic_rag
  Skill-reasoning   v3_armAB_fullgrid.json           arm=A_reasoning
  Skill-execution   v3_armAB_fullgrid.json           arm=B_execution
  Control           v3_raw_rescored_three_arm.json   cond=with_spec

Metric is A1, phenotype accuracy against the reference. Panel is the
manuscript's eight models; Mistral Large 2 is excluded here as it is
everywhere else in the paper. Rows that failed to parse are counted as
failures in the three-arm source via its format_fail flag, matching
36-figure6-population-equity.py.

Run:  python3 code/test_figure_s8.py && python3 code/90-figure-s8-conditional-invariance.py
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
FIGDIR = BASE / "figures"

THREE_ARM = "v4_flash_population_three_arm_panel.json"
FULLGRID = "v4_flash_population_skill_arms_panel.json"

EXCLUDED_MODEL = "Mistral Large 2"
POPULATIONS = ["EUR", "AMR", "AFR"]
POPULATION_LABELS = {
    "EUR": "EUR (European)",
    "AMR": "AMR (Latin American)",
    "AFR": "AFR (East African)",
}
POPULATION_COLOURS = {
    "EUR": "#35b0a7",
    "AMR": "#3b6ea5",
    "AFR": "#3b2a55",
}

TITLE = "Conditional invariance is a property of the control arm, not of the pipeline"
SUBTITLE = (
    "Identical curated genotypes under European, Latin American and East African "
    "framing. The EUR minus AFR gap is zero only in the control, and largest "
    "under skill-execution."
)
YLABEL = "Phenotype concordance, A1 (%)"

# The supplement holds Figure S8 at a fixed wp:extent whose ratio is 1.6577.
# Drawing to any other ratio squashes the picture when the image bytes are
# swapped in place, so the canvas is derived from the extent rather than
# chosen by eye.
DOCX_ASPECT = 5760720 / 3475167
FIG_WIDTH_IN = 10.5

ARMS = [
    {"key": "free_prompt", "label": "Free-prompt", "source": THREE_ARM,
     "field": "cond", "value": "no_spec"},
    {"key": "rag", "label": "RAG", "source": THREE_ARM,
     "field": "cond", "value": "cpic_rag"},
    {"key": "skill_reasoning", "label": "Skill-\nreasoning", "source": FULLGRID,
     "field": "arm", "value": "A_reasoning"},
    {"key": "skill_execution", "label": "Skill-\nexecution", "source": FULLGRID,
     "field": "arm", "value": "B_execution"},
    {"key": "control", "label": "Control", "source": THREE_ARM,
     "field": "cond", "value": "with_spec"},
]

_CACHE: dict[str, list[dict]] = {}


def load(name: str) -> list[dict]:
    if name not in _CACHE:
        path = DATA / name
        if not path.exists():
            raise SystemExit(
                f"missing {path}\n"
                f"It is in the Zenodo deposit. Fetch it with:\n"
                f"  curl -sL -o {path} "
                f"'https://zenodo.org/records/21969051/files/{name}?download=1'"
            )
        _CACHE[name] = json.loads(path.read_text())
    return _CACHE[name]


def score(row: dict, drop_parse_failures: bool = False) -> float | None:
    """A1 for a row, from either file's shape, or None to drop it.

    The three-arm file nests scores and flags unparsed responses; the fullgrid
    file is flat and has no unparsed rows at all, so this choice only bites on
    the 292 flagged rows of the three-arm source.

    The default counts a parse failure as a failure rather than as a missing
    observation, which is what the manuscript says it does: "missing, empty and
    unparsed outputs count as failures". The submitted figure did the opposite
    and dropped them, following the older 36-figure6-population-equity.py, and
    that is why its free-prompt and RAG bars sit about 1.5 points higher than
    the ones drawn here. Pass drop_parse_failures to reproduce the old bars.
    """
    if "scores" in row:
        scores = row["scores"]
        if scores.get("format_fail"):
            return None if drop_parse_failures else 0.0
        return scores.get("A1")
    return row.get("A1")


def wilson(successes: float, n: int) -> tuple[float, float]:
    """Wilson 95% interval, in percentage points."""
    if n == 0:
        return (0.0, 0.0)
    z = 1.959963984540054
    p = successes / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (100 * (centre - half), 100 * (centre + half))


def compute(drop_parse_failures: bool = False) -> list[dict]:
    """Every number the figure prints, with the source that produced it."""
    out: list[dict] = []
    for arm in ARMS:
        rows = [r for r in load(arm["source"]) if r.get(arm["field"]) == arm["value"]]
        rows = [r for r in rows if r.get("model") != EXCLUDED_MODEL]
        models = {r["model"] for r in rows}

        pops: dict[str, dict] = {}
        for pop in POPULATIONS:
            values = [score(r, drop_parse_failures)
                      for r in rows if r.get("pop") == pop]
            values = [v for v in values if v is not None]
            n = len(values)
            successes = sum(values)
            lo, hi = wilson(successes, n)
            pops[pop] = {
                "pct": 100 * successes / n if n else float("nan"),
                "n": n,
                "lo": lo,
                "hi": hi,
            }

        out.append({
            "key": arm["key"],
            "label": arm["label"],
            "source": arm["source"],
            "selector": f"{arm['field']}={arm['value']}",
            "models": models,
            "n_models": len(models),
            "pops": pops,
            "delta": pops["EUR"]["pct"] - pops["AFR"]["pct"],
        })
    return out


def report(arms: list[dict]) -> None:
    print(f"{'arm':17s} {'source':34s} {'EUR':>6s} {'AMR':>6s} {'AFR':>6s} "
          f"{'EUR-AFR':>8s} {'n/pop':>7s}")
    for arm in arms:
        p = arm["pops"]
        print(f"{arm['key']:17s} {arm['source']:34s} "
              f"{p['EUR']['pct']:6.2f} {p['AMR']['pct']:6.2f} {p['AFR']['pct']:6.2f} "
              f"{arm['delta']:8.2f} {p['EUR']['n']:7d}")


def draw(arms: list[dict],save_output=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.edgecolor": "#333333",
    })

    fig, ax = plt.subplots(figsize=(FIG_WIDTH_IN, FIG_WIDTH_IN / DOCX_ASPECT))

    width = 0.26
    centres = np.arange(len(arms), dtype=float)

    for i, pop in enumerate(POPULATIONS):
        offset = (i - 1) * width
        heights = [a["pops"][pop]["pct"] for a in arms]
        # Wilson is not centred on p, so at p = 1 the arithmetic can land a
        # hair below zero. Clamp rather than let matplotlib refuse to draw.
        lower = [max(0.0, a["pops"][pop]["pct"] - a["pops"][pop]["lo"]) for a in arms]
        upper = [max(0.0, a["pops"][pop]["hi"] - a["pops"][pop]["pct"]) for a in arms]
        ax.bar(centres + offset, heights, width,
               color=POPULATION_COLOURS[pop], label=POPULATION_LABELS[pop],
               zorder=3)
        ax.errorbar(centres + offset, heights, yerr=[lower, upper],
                    fmt="none", ecolor="#222222", elinewidth=1.0,
                    capsize=2.5, capthick=1.0, zorder=4)

    for centre, arm in zip(centres, arms):
        top = max(arm["pops"][p]["hi"] for p in POPULATIONS)
        ax.annotate(
            f"EUR - AFR = {arm['delta']:.2f}",
            xy=(centre, top + 2.0),
            ha="center", va="bottom", fontsize=9.5,
            color="#222222" if arm["key"] != "control" else "#4a7a4a",
            fontweight="bold" if arm["key"] == "control" else "normal",
        )

    ax.set_xticks(centres)
    ax.set_xticklabels([a["label"] for a in arms], fontsize=11)
    ax.set_ylabel(YLABEL, fontsize=11.5)
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#e6e6e6", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_title(TITLE, fontsize=13.5, pad=26, loc="left")
    ax.text(0.0, 1.015, SUBTITLE, transform=ax.transAxes, fontsize=9.3,
            color="#555555", ha="left", va="bottom", wrap=True)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
              frameon=False, fontsize=10.5)

    fig.tight_layout()
    if not save_output:return fig
    FIGDIR.mkdir(exist_ok=True)
    png = FIGDIR / "FigureS8_conditional_invariance.png"
    tiff = FIGDIR / "FigureS8_conditional_invariance.tiff"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(tiff, dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"\nwrote {png}")
    print(f"wrote {tiff}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numbers-only", action="store_true",
                        help="print the computed table and do not draw")
    parser.add_argument("--drop-parse-failures", action="store_true",
                        help="drop unparsed rows instead of scoring them as "
                             "failures; reproduces the submitted figure's bars "
                             "and contradicts the manuscript's stated policy")
    args = parser.parse_args()

    arms = compute(args.drop_parse_failures)
    report(arms)

    if args.numbers_only:
        return 0
    draw(arms)
    return 0


if __name__ == "__main__":
    sys.exit(main())
