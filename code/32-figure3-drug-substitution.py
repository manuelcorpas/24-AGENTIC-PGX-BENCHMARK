#!/usr/bin/env python3
"""
Figure 3 for Cell Genomics Manuscript v8.

The drug-substitution finding: retrieval-augmented generation regresses on the A2
(drug-recommendation) dimension because CPIC guideline chunks are keyed on gene
rather than (gene, drug). A retrieved chunk for CYP2D6 contains recommendations
for codeine, tamoxifen, ondansetron, paroxetine and amitriptyline together;
models often echo the lead drug rather than the queried drug.

Two-panel design:
  Panel A - horizontal grouped bars: per-gene A2 mean under free-prompted,
            retrieval-augmented and specification-constrained conditions, sorted
            so that genes where RAG most underperforms no_spec appear at the top.
  Panel B - worked example callout: a real Claude Sonnet 4 response on the
            cyp2d6_tamox_im case under cpic_rag, illustrating the substitution
            mechanism (queried drug = tamoxifen; model answered for codeine).

Per-gene A2 means in Panel A are computed directly from the locked 3-of-3
replicate three-arm dataset. The drug-substitution rate annotation in Panel A
and the worked example in Panel B draw from v3_three_arm_a2_regression_classified.csv,
which is the locked 3-run classifier output (33-classify-a2-regressions.py).
Headline: 470/973 (48.3%) of A2 regressions are drug-substitution; combined
chunk/multi-drug structural confusion = 922/973 (94.8%).

Inputs:
  RESULTS/v3_raw_rescored_three_arm.json            (per-gene A2 means; 3-run locked)
  RESULTS/v3_three_arm_a2_regression_classified.csv (worked example + drug-sub rate)

Outputs (PNG 300 DPI + TIFF 600 DPI LZW):
  FIGURES/Figure3_drug_substitution.png
  FIGURES/Figure3_drug_substitution.tiff

Usage:
  python3 32-figure3-drug-substitution.py             # render
  python3 32-figure3-drug-substitution.py --self-test # verify inputs
"""
from __future__ import annotations

import argparse
import csv
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
MAIN = BASE / "RESULTS" / "v3_raw_rescored_three_arm.json"
REGRESSIONS = BASE / "RESULTS" / "v3_three_arm_a2_regression_classified.csv"
FIGDIR = BASE / "FIGURES"

CONDITIONS = ["no_spec", "cpic_rag", "with_spec"]
CONDITION_LABELS = {
    "no_spec": "Free-prompted",
    "cpic_rag": "Retrieval-augmented",
    "with_spec": "Specification-constrained",
}
CONDITION_COLOURS = {
    "no_spec": "#d4756f",
    "cpic_rag": "#e0a96d",
    "with_spec": "#5a8f5a",
}

EXPECTED = {
    "csv_n_rows": 973,
    "csv_drug_sub": 470,
    "csv_phrasing": 51,
    "csv_wrong_dir": 176,
    "csv_other": 276,
    "n_genes": 21,
    "cyp2d6_cpic_rag_a2_approx": 0.413,
    "cyp2d6_no_spec_a2_approx": 0.649,
    "worked_example_tc": "cyp2d6_tamox_im",
}


def load_main() -> list[dict]:
    return json.loads(MAIN.read_text())


def load_regressions() -> list[dict]:
    with REGRESSIONS.open() as fh:
        return list(csv.DictReader(fh))


def per_gene_a2(main: list[dict]) -> dict[str, dict[str, tuple[float, int]]]:
    by_gene: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in main:
        if r["scores"].get("format_fail"):
            continue
        a2 = r["scores"].get("A2")
        if a2 is None:
            continue
        by_gene[r["gene"]][r["cond"]].append(a2)
    out: dict[str, dict[str, tuple[float, int]]] = {}
    for gene, cond_map in by_gene.items():
        out[gene] = {c: (sum(vals) / len(vals) if vals else float("nan"), len(vals))
                     for c, vals in cond_map.items()}
    return out


def drug_sub_rate_per_gene(regressions: list[dict], main: list[dict]) -> dict[str, float]:
    tc2gene = {r["tc"]: r["gene"] for r in main}
    totals: Counter = Counter()
    subs: Counter = Counter()
    for r in regressions:
        gene = tc2gene.get(r["tc"])
        if gene is None:
            continue
        totals[gene] += 1
        if r["category"] == "a_drug_substitution":
            subs[gene] += 1
    return {g: (subs[g] / totals[g] if totals[g] else 0.0) for g in totals}


def worked_example(regressions: list[dict]) -> dict:
    candidates = [r for r in regressions
                  if r["tc"] == EXPECTED["worked_example_tc"]
                  and r["category"] == "a_drug_substitution"
                  and r["model"] == "Claude Sonnet 4"
                  and r["pop"] == "EUR"]
    if not candidates:
        raise RuntimeError("worked example not found")
    return candidates[0]


def self_test() -> int:
    failures: list[str] = []

    regs = load_regressions()
    if len(regs) != EXPECTED["csv_n_rows"]:
        failures.append(f"regression CSV row count: got {len(regs)}, want {EXPECTED['csv_n_rows']}")

    cats = Counter(r["category"] for r in regs)
    if cats.get("a_drug_substitution", 0) != EXPECTED["csv_drug_sub"]:
        failures.append(f"a_drug_substitution count: got {cats.get('a_drug_substitution', 0)}, "
                        f"want {EXPECTED['csv_drug_sub']}")
    if cats.get("b_phrasing_equiv", 0) != EXPECTED["csv_phrasing"]:
        failures.append(f"b_phrasing_equiv count: got {cats.get('b_phrasing_equiv', 0)}, "
                        f"want {EXPECTED['csv_phrasing']}")
    if cats.get("c_wrong_direction", 0) != EXPECTED["csv_wrong_dir"]:
        failures.append(f"c_wrong_direction count: got {cats.get('c_wrong_direction', 0)}, "
                        f"want {EXPECTED['csv_wrong_dir']}")
    if cats.get("d_other", 0) != EXPECTED["csv_other"]:
        failures.append(f"d_other count: got {cats.get('d_other', 0)}, "
                        f"want {EXPECTED['csv_other']}")

    main = load_main()
    pg = per_gene_a2(main)
    if len(pg) != EXPECTED["n_genes"]:
        failures.append(f"distinct gene count: got {len(pg)}, want {EXPECTED['n_genes']}")

    cyp2d6 = pg.get("CYP2D6", {})
    if abs(cyp2d6.get("cpic_rag", (0, 0))[0] - EXPECTED["cyp2d6_cpic_rag_a2_approx"]) > 0.005:
        failures.append(f"CYP2D6 cpic_rag A2: got {cyp2d6.get('cpic_rag', (0,0))[0]:.4f}, "
                        f"want ~{EXPECTED['cyp2d6_cpic_rag_a2_approx']}")
    if abs(cyp2d6.get("no_spec", (0, 0))[0] - EXPECTED["cyp2d6_no_spec_a2_approx"]) > 0.005:
        failures.append(f"CYP2D6 no_spec A2: got {cyp2d6.get('no_spec', (0,0))[0]:.4f}, "
                        f"want ~{EXPECTED['cyp2d6_no_spec_a2_approx']}")
    if cyp2d6.get("with_spec", (0, 0))[0] != 1.0:
        failures.append(f"CYP2D6 with_spec A2: got {cyp2d6.get('with_spec', (0,0))[0]}, want 1.0")

    try:
        worked_example(load_regressions())
    except RuntimeError as exc:
        failures.append(f"worked example lookup: {exc}")

    if failures:
        print("SELF-TEST FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"SELF-TEST PASSED ({EXPECTED['csv_n_rows']} regression rows; "
          f"{EXPECTED['n_genes']} genes; CYP2D6 A2 (no_spec/cpic_rag/with_spec) "
          f"= {cyp2d6['no_spec'][0]:.3f}/{cyp2d6['cpic_rag'][0]:.3f}/{cyp2d6['with_spec'][0]:.3f}).")
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
    main = load_main()
    regs = load_regressions()
    pg = per_gene_a2(main)
    drug_sub = drug_sub_rate_per_gene(regs, main)
    we = worked_example(regs)

    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 300,
    })

    genes = list(pg.keys())
    delta_rag_vs_nospec = {g: pg[g]["cpic_rag"][0] - pg[g]["no_spec"][0] for g in genes}
    genes_sorted = sorted(genes, key=lambda g: delta_rag_vs_nospec[g])

    fig = plt.figure(figsize=(10.5, 9.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.7, 1.0], hspace=0.36)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])

    y = np.arange(len(genes_sorted))
    height = 0.26

    for i, cond in enumerate(CONDITIONS):
        offsets = (1 - i) * height
        vals = [100 * pg[g][cond][0] for g in genes_sorted]
        axA.barh(y + offsets, vals, height,
                 label=CONDITION_LABELS[cond],
                 color=CONDITION_COLOURS[cond],
                 edgecolor="black", linewidth=0.3)

    axA.set_yticks(y)
    axA.set_yticklabels(genes_sorted)
    axA.set_xlabel("A2 mean accuracy (%)")
    axA.set_xlim(0, 134)
    axA.set_xticks([0, 25, 50, 75, 100])
    axA.axvline(100, color="grey", linewidth=0.5, linestyle=":")
    axA.set_title("A  Per-gene drug-recommendation accuracy by condition "
                  "(sorted by retrieval-augmented vs free-prompted delta; worst RAG regressions at top)",
                  loc="left", fontweight="bold", fontsize=9.5)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    axA.invert_yaxis()
    axA.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0), frameon=False, fontsize=8)

    axA.text(112, -0.7, "Drug-substitution rate\n(% of regressions per gene)",
             fontsize=7, color="#7a3a32", style="italic", va="center", ha="left", fontweight="bold")
    for i, g in enumerate(genes_sorted):
        rate = drug_sub.get(g)
        if rate is None or rate < 0.05:
            continue
        axA.text(112, i + 0.05,
                 f"{rate * 100:.0f}%",
                 fontsize=7.5, color="#7a3a32", va="center", fontweight="bold")

    axB.axis("off")
    axB.set_xlim(0, 1); axB.set_ylim(0, 1)

    box_xy = (0.01, 0.04)
    box_w, box_h = 0.98, 0.92
    axB.add_patch(mpatches.FancyBboxPatch(box_xy, box_w, box_h,
                                          boxstyle="round,pad=0.012,rounding_size=0.012",
                                          linewidth=1.0, edgecolor="#7a3a32",
                                          facecolor="#fbf2ee"))

    def wrap(text: str, width: int) -> str:
        import textwrap
        return "\n".join(textwrap.wrap(text, width=width)) or text

    axB.text(0.03, 0.92,
             "B  Worked example: CYP2D6 / tamoxifen (CYP2D6 IM, EUR), Claude Sonnet 4, cpic_rag",
             fontsize=9.5, fontweight="bold", va="top", color="#5a2820")

    axB.text(0.03, 0.78, "Queried drug:", fontsize=8.5, fontweight="bold", color="#333333", va="top")
    axB.text(0.21, 0.78, we["queried_drug"], fontsize=8.5, color="#333333", va="top")

    axB.text(0.03, 0.66, "CPIC truth:", fontsize=8.5, fontweight="bold", color="#2a5a2a", va="top")
    axB.text(0.21, 0.66, wrap(we["gt_drug"], 90),
             fontsize=8.5, color="#2a5a2a", va="top", family="monospace")

    axB.text(0.03, 0.48, "Model answered:", fontsize=8.5, fontweight="bold", color="#7a3a32", va="top")
    axB.text(0.21, 0.48, wrap(we["cpic_rag_DRUG"], 90),
             fontsize=8.5, color="#7a3a32", va="top", family="monospace")

    axB.text(0.03, 0.28, "Mechanism:", fontsize=8.5, fontweight="bold", color="#333333", va="top")
    axB.text(0.21, 0.28,
             wrap("The retrieved CPIC chunk for CYP2D6 bundles recommendations for codeine, "
                  "tamoxifen, ondansetron, paroxetine and amitriptyline. The model echoed "
                  "the lead drug (codeine) rather than the queried drug (tamoxifen). "
                  "A2 = 0 because the response is not about the queried drug.", 95),
             fontsize=8.5, color="#333333", va="top")

    fig.text(0.5, 0.01,
             "n = 8,738 / 8,790 / 8,910 parsed A2 evaluations per condition. "
             "Per-gene drug-substitution rate computed on 973 A2 regression combos "
             "(locked 3-run classifier; 33-classify-a2-regressions.py).",
             ha="center", fontsize=7, color="#555555")

    save_fig(fig, "Figure3_drug_substitution")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1] if __doc__ else None)
    parser.add_argument("--self-test", action="store_true",
                        help="Verify inputs; do not render.")
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
