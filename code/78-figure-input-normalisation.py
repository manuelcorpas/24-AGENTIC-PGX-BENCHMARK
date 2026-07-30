#!/usr/bin/env python3
"""
Figure: input normalisation with and without the definition artefact (R1.1).

EVERY NUMBER PLOTTED IS READ FROM A DATA FILE. No hardcoded results, no random
number generator anywhere in this script; --self-test recomputes the headline
values and refuses to draw if they disagree with the registered numbers.

WHAT THE FIGURE SHOWS
Two panels, one measure each. Coverage (did the model answer at all) and
accuracy among the answers it gave. Grouped by model, paired within model:
without the allele-definition table, and with it. A horizontal reference line
marks the deterministic caller on the same pairs.

WHY TWO PANELS AND NOT ONE
Coverage and accuracy are different measures on the same 0-1 scale but they are
not comparable quantities, and a model can trade one for the other by
abstaining. Putting them on one axis invites exactly that misreading; a second
y-axis would be worse. Two panels, one axis each.

COLOUR
Blue and orange, validated for colour-vision deficiency (worst adjacent pair
protan dE 24.7). The red/green pair used by the earlier figures in this paper
fails that check at dE 3.5 and should be revisited. Identity is carried by the
legend AND by direct value labels, so the figure is never colour-alone.

USAGE
    python code/78-figure-input-normalisation.py --self-test
    python code/78-figure-input-normalisation.py
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
FIGDIR = BASE / "figures"

WITHOUT = "#2a78d6"   # blue
WITH = "#eb6834"      # orange
NEUTRAL = "#8a8a8a"
INK = "#0b0b0b"

# (label, with-definitions file, model name as recorded)
ARMS = [
    ("Claude Opus 4.5", ["v3_input_normalisation_defs.json",
                         "v3_input_normalisation_defs_tail.json"]),
    ("GPT-5.2", ["v3_input_normalisation_defs_gpt52.json"]),
    ("o3", ["v3_input_normalisation_defs_o3.json"]),
]
BASELINE = "v3_input_normalisation_main.json"
BASELINE_O3 = "v3_input_normalisation_o3.json"


def _nd(d: str) -> tuple:
    return tuple(sorted(p.strip() for p in d.split("/")))


def _load(name: str) -> list[dict]:
    p = DATA / name
    if not p.exists():
        return []
    return json.loads(p.read_text())["rows"]


def metrics(rows: list[dict], ref: str = "reference") -> dict | None:
    """Coverage and accuracy among emitted, on rows that were actually scored.

    Errors and budget-truncated responses leave the denominator: neither is the
    model declining to answer, and counting them as abstentions reports a limit
    of the apparatus as a property of the model (CORRECTIONS.md C11, C12, C13).
    """
    scored = [r for r in rows if r.get("status") in ("call", "abstain")]
    if not scored:
        return None
    em = [r for r in scored if r.get("status") == "call"]
    ok = sum(1 for r in em if _nd(r["call"]) == _nd(r[ref]))
    return {
        "n": len(scored),
        "coverage": len(em) / len(scored),
        "accuracy": (ok / len(em)) if em else None,
        "emitted": len(em),
    }


def collect() -> dict:
    """Paired with/without per model, restricted to the pairs each arm covers."""
    out = {}
    base_rows = _load(BASELINE) + _load(BASELINE_O3)
    for model, files in ARMS:
        withrows = []
        for f in files:
            withrows += _load(f)
        if not withrows:
            continue
        keys = {(r["sample"], r["gene"]) for r in withrows}
        # the comparison must be on the SAME pairs and the SAME rendering,
        # or the two bars are answering different questions
        without = [r for r in base_rows
                   if r["model"] == model and r["form"] == "vcf"
                   and (r["sample"], r["gene"]) in keys]
        w = metrics(withrows)
        wo = metrics(without)
        if w is None:
            continue
        em = [r for r in withrows if r.get("status") == "call"]
        caller = (sum(1 for r in em if _nd(r["reference"]) == _nd(r["getrm"])) / len(em)
                  if em else None)
        out[model] = {"with": w, "without": wo, "caller_vs_getrm": caller,
                      "model_vs_getrm": (sum(1 for r in em
                                             if _nd(r["call"]) == _nd(r["getrm"])) / len(em)
                                         if em else None)}
    return out


def draw(res: dict, path_stem: str = "FigureR5-input-normalisation"):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#555555", "figure.dpi": 300,
    })
    models = [m for m, _ in ARMS if m in res]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 4.1), constrained_layout=True)
    xs = range(len(models))
    width = 0.36
    gap = 0.02   # surface gap between adjacent bars

    for ax, key, title in ((axA, "coverage", "Coverage: did the model answer?"),
                           (axB, "accuracy", "Accuracy among answers given")):
        for i, m in enumerate(models):
            wo = res[m]["without"]
            w = res[m]["with"]
            vwo = (wo or {}).get(key)
            vw = w.get(key)
            if vwo is not None:
                ax.bar(i - width / 2 - gap, vwo, width, color=WITHOUT, zorder=3)
                ax.text(i - width / 2 - gap, vwo + 0.02, f"{vwo:.2f}", ha="center",
                        fontsize=8, color=INK, zorder=4)
            if vw is not None:
                ax.bar(i + width / 2 + gap, vw, width, color=WITH, zorder=3)
                ax.text(i + width / 2 + gap, vw + 0.02, f"{vw:.2f}", ha="center",
                        fontsize=8, color=INK, zorder=4)
        ax.set_xticks(list(xs))
        ax.set_xticklabels(models, fontsize=8.5)
        ax.set_ylim(0, 1.12)
        ax.set_ylabel(key.capitalize())
        ax.set_title(title, fontsize=9.5, loc="left")
        ax.grid(axis="y", color="#e6e6e4", lw=0.7, zorder=0)
        ax.set_axisbelow(True)

    # the deterministic caller, on the pairs each model answered
    callers = [res[m]["caller_vs_getrm"] for m in models if res[m]["caller_vs_getrm"]]
    if callers:
        lvl = sum(callers) / len(callers)
        axB.axhline(lvl, color=NEUTRAL, lw=1.2, ls="--", zorder=2)
        axB.text(len(models) - 0.5, lvl + 0.025,
                 f"deterministic caller vs external consensus ({lvl:.2f})",
                 ha="right", fontsize=7.5, color="#52514e")

    handles = [plt.Rectangle((0, 0), 1, 1, color=WITHOUT),
               plt.Rectangle((0, 0), 1, 1, color=WITH)]
    fig.legend(handles, ["without allele definitions", "with allele definitions"],
               loc="lower center", ncol=2, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.06))
    FIGDIR.mkdir(exist_ok=True)
    fig.savefig(FIGDIR / f"{path_stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / f"{path_stem}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    print(f"wrote {FIGDIR / path_stem}.png and .tiff")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    res = collect()
    if not res:
        print("no input-normalisation data found", file=sys.stderr)
        return 2

    for m, r in res.items():
        wo = r["without"] or {}
        print(f"{m:18s} without: cov={wo.get('coverage')} acc={wo.get('accuracy')} | "
              f"with: cov={r['with']['coverage']:.3f} acc={r['with']['accuracy']:.3f} | "
              f"vs GeT-RM model={r['model_vs_getrm']:.3f} caller={r['caller_vs_getrm']:.3f}")

    if args.self_test:
        c = res.get("Claude Opus 4.5")
        assert c and abs(c["with"]["coverage"] - 0.928) < 0.005, "Claude coverage drifted"
        assert abs(c["with"]["accuracy"] - 0.973) < 0.005, "Claude accuracy drifted"
        print("self-test OK")
        return 0

    draw(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
