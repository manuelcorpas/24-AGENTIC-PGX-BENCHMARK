#!/usr/bin/env python3
"""
Step 6 — decompose the agent's real-genome predictions and render the cross-population
figure. Reads the combined predictions (step 4 output for all cohorts), classifies each
determinate-cohort call as correct / abstained / wrong / parse-fail, and writes:
  - a summary table (aggregate per-population counts; safe to share)
  - the cross-population figure (PNG 300 dpi + TIFF 600 dpi)

Self-contained scorer with PGx phenotype-tier vocabulary (metaboliser, function,
acetylator, expresser, responder, susceptibility). Connection-error predictions are
dropped (API failures, not results).

Usage: 06_decompose_and_figure.py <predictions.tsv> <out_prefix>
  predictions.tsv columns: cohort, gene, diplotype, cohort_phenotype, model, pred
"""
import sys
import re
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# population display order = increasing genetic distance from the European-centric
# training literature; labels are generic so the script is cohort-agnostic.
ORDER = [
    ("CorpasFamily", "European\n(family)"),
    ("Peru", "Latin American"),
    ("UGR", "East African"),
]
CANONICAL = 96.0  # curated-benchmark phenotype accuracy, for reference

TIERS = [
    ("MH", r"malignant hyperthermia|mh[ -]?susceptib|uncertain susceptib"),
    ("ULTRA", r"ultra[ -]?rapid"), ("RAPID_AC", r"rapid acetylat"), ("RAPID", r"rapid metaboli"),
    ("SLOW_AC", r"slow acetylat"), ("INT_AC", r"intermediate acetylat"),
    ("POOR_M", r"poor metaboli"), ("IM", r"intermediate metaboli"),
    ("NORMAL_M", r"normal metaboli|extensive metaboli"),
    ("NONEXP", r"non[ -]?expresser"), ("EXP", r"expresser"),
    ("POOR_FN", r"poor function|low function"),
    ("DEC_FN", r"decreased function|reduced function|possible decreased|reduced response|reduced dpd"),
    ("INC_FN", r"increased function"),
    ("NORM_FN", r"normal function|normal dpd|normal responder|normal response"),
    ("FAV", r"favou?rable"), ("UNFAV", r"unfavou?rable|poor response"),
    ("POS", r"positive|carrier|at[ -]?risk"), ("NEG", r"negative|non[ -]?carrier|absent"),
    ("DEF", r"deficien"), ("IND", r"indetermin|uncertain|unknown|variable"),
]
FAMILY = {"NORMAL_M": "N", "NORM_FN": "N", "EXP": "N", "FAV": "N", "RAPID": "R", "RAPID_AC": "R",
          "ULTRA": "U", "IM": "I", "INT_AC": "I", "DEC_FN": "I", "POOR_M": "P", "POOR_FN": "P",
          "NONEXP": "P", "SLOW_AC": "P", "UNFAV": "P", "POS": "POS", "NEG": "NEG", "DEF": "DEF", "MH": "MH"}

def tier(s):
    s = (s or "").lower()
    for n, p in TIERS:
        if re.search(p, s):
            return n
    return None

def main():
    src, outpfx = sys.argv[1], sys.argv[2]
    rows = list(csv.DictReader(open(src), delimiter="\t"))
    cat = defaultdict(lambda: defaultdict(int))
    for r in rows:
        if "connection error" in r["pred"].lower():
            continue
        gt, pred = r["cohort_phenotype"], r["pred"]
        if tier(gt) == "IND":
            continue
        c = cat[r["cohort"]]; c["n"] += 1
        pt, gtt = tier(pred), tier(gt)
        if pt is None or not pred.strip():
            c["parsefail"] += 1
        elif pt == "IND":
            c["abstain"] += 1
        elif pt == gtt or (FAMILY.get(pt) and FAMILY.get(pt) == FAMILY.get(gtt)):
            c["correct"] += 1
        else:
            c["wrong"] += 1

    # summary table
    with open(f"{outpfx}_summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["population", "n", "correct_pct", "abstain_pct", "wrong_pct", "parsefail_pct"])
        for coh, lab in ORDER:
            c = cat[coh]; n = c["n"]
            if not n: continue
            w.writerow([lab.replace("\n", " "), n,
                        round(100*c["correct"]/n, 1), round(100*c["abstain"]/n, 1),
                        round(100*c["wrong"]/n, 1), round(100*c["parsefail"]/n, 1)])

    # figure
    plt.rcParams.update({"font.family": "Helvetica", "font.size": 11, "axes.titlesize": 12,
                         "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 300})
    mako = plt.get_cmap("mako") if "mako" in plt.colormaps() else None
    C = {"correct": "#2A7B7B", "abstain": "#E0A96D", "wrong": "#D1495B", "parsefail": "#CFCFCF"}
    labels = ["Curated\nbenchmark"] + [lab for _, lab in ORDER if cat[_]["n"]]
    cohs = [c for c, _ in ORDER if cat[c]["n"]]
    def pct(c, k): n = cat[c]["n"]; return 100*cat[c][k]/n if n else 0
    correct = [CANONICAL] + [pct(c, "correct") for c in cohs]
    abstain = [0] + [pct(c, "abstain") for c in cohs]
    wrong = [100-CANONICAL] + [pct(c, "wrong") for c in cohs]
    parse = [0] + [pct(c, "parsefail") for c in cohs]
    import numpy as np
    x = np.arange(len(labels)); w = 0.6
    fig, ax = plt.subplots(figsize=(8.4, 5))
    b = [0]*len(labels)
    for key, vals, lab in [("correct", correct, "Correct (matches CPIC)"),
                           ("abstain", abstain, "Abstained (uncertain allele)"),
                           ("wrong", wrong, "Wrong (clinically incorrect)"),
                           ("parsefail", parse, "Unparsed")]:
        ax.bar(x, vals, w, bottom=b, label=lab, color=C[key], edgecolor="white", lw=0.7, zorder=3)
        b = [bb+vv for bb, vv in zip(b, vals)]
    for i, cv in enumerate(correct):
        ax.text(i, cv/2, f"{cv:.0f}%", ha="center", va="center", color="white", fontweight="bold", fontsize=10.5)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("Agent responses (%)"); ax.set_ylim(0, 100)
    ax.yaxis.grid(True, color="#ECECEC"); ax.set_axisbelow(True); ax.tick_params(length=0)
    ax.legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.28), fontsize=9)
    ax.set_title("Curated accuracy does not transfer to real genomes, and degrades with ancestry",
                 loc="left", pad=10)
    fig.savefig(f"{outpfx}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{outpfx}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    print(f"wrote {outpfx}.png / .tiff and {outpfx}_summary.tsv")
    for coh, lab in ORDER:
        c = cat[coh]; n = c["n"]
        if n: print(f"  {lab.replace(chr(10),' '):22} n={n:5} correct={100*c['correct']/n:.0f}% abstain={100*c['abstain']/n:.0f}% wrong={100*c['wrong']/n:.0f}%")

if __name__ == "__main__":
    main()
