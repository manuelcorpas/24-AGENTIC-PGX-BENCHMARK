#!/usr/bin/env python3
"""
Regenerate v3 manuscript figures from v3_raw_rescored.json.

Outputs (PNG + TIFF for journal submission):
  Figure1_v3_model_accuracy.png/tiff   - per-model A1 accuracy, no_spec vs with_spec
  Figure2_v3_gene_heatmap.png/tiff     - heatmap of A1 by gene-drug pair × model, two panels
  Figure3_v3_lethal_by_gene.png/tiff   - lethal-case errors by gene, no_spec
  Figure4_v3_population.png/tiff       - per-model A1 by population × model, no_spec

These figures replace v2 Figure 1 (consistency heatmap) and Figure 2 (population
accuracy). The deck of 4 panels is sized for BiB column widths (single ~85mm or
double ~178mm).
"""
from __future__ import annotations
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "RESULTS" / "v3_raw_rescored.json"
CASES = BASE / "SPECS" / "test_cases_v3.json"
FIGDIR = BASE / "FIGURES"
FIGDIR.mkdir(exist_ok=True)

# Style
plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
})


def save_fig(fig, name: str):
    fig.savefig(FIGDIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / f"{name}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  wrote {name}.png + .tiff")


def main():
    rows = json.loads(RAW.read_text())
    cases = json.loads(CASES.read_text())
    cbi = {c["id"]: c for c in cases}

    parsed = [r for r in rows if not r["scores"].get("format_fail")]
    models = sorted({r["model"] for r in rows})
    pops = ["EUR", "AMR", "AFR"]
    gene_drug_pairs = sorted({(c["gene"], c["drug"]) for c in cases}, key=lambda x: (x[0], x[1]))

    # =========================================================================
    # FIGURE 1: per-model A1 accuracy, no_spec vs with_spec
    # =========================================================================
    print("Building Figure 1 (per-model accuracy)...")
    fig, ax = plt.subplots(figsize=(7, 4.2))
    no_spec_pct = []
    with_spec_pct = []
    for m in models:
        ns = [r for r in parsed if r["model"] == m and r["cond"] == "no_spec"]
        ws = [r for r in parsed if r["model"] == m and r["cond"] == "with_spec"]
        no_spec_pct.append(100 * sum(1 for r in ns if r["scores"]["A1"] == 1.0) / len(ns))
        with_spec_pct.append(100 * sum(1 for r in ws if r["scores"]["A1"] == 1.0) / len(ws))

    x = np.arange(len(models))
    width = 0.4
    ax.bar(x - width/2, no_spec_pct, width, label="without skill", color="#d4756f")
    ax.bar(x + width/2, with_spec_pct, width, label="with ClawBio skill", color="#5a8f5a")
    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=35, ha="right")
    ax.set_ylabel("Phenotype accuracy A1 (%)")
    ax.set_ylim(0, 105)
    ax.axhline(100, linestyle=":", color="gray", linewidth=0.8)
    ax.legend(loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    save_fig(fig, "Figure1_v3_model_accuracy")

    # =========================================================================
    # FIGURE 2: heatmap A1 by gene-drug pair × model, two panels
    # =========================================================================
    print("Building Figure 2 (gene-drug × model heatmap)...")
    # Aggregate A1 per (model, gene, drug) for each condition (mean across pop and tier)
    def heatmap_data(cond):
        mat = np.full((len(models), len(gene_drug_pairs)), np.nan)
        for i, m in enumerate(models):
            for j, (gene, drug) in enumerate(gene_drug_pairs):
                rs = [r for r in parsed
                      if r["model"] == m and r["cond"] == cond
                      and cbi[r["tc"]]["gene"] == gene and cbi[r["tc"]]["drug"] == drug]
                if rs:
                    mat[i, j] = sum(r["scores"]["A1"] for r in rs) / len(rs)
        return mat

    mat_ns = heatmap_data("no_spec")
    mat_ws = heatmap_data("with_spec")

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.2), sharey=True)
    labels = [f"{g} / {d}" for g, d in gene_drug_pairs]
    for ax, mat, title in [(axes[0], mat_ns, "without skill"), (axes[1], mat_ws, "with ClawBio skill")]:
        im = ax.imshow(mat, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_yticks(np.arange(len(models)))
        ax.set_yticklabels(models)
        ax.set_title(title)
    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("A1 accuracy (mean)")
    save_fig(fig, "Figure2_v3_gene_heatmap")

    # =========================================================================
    # FIGURE 3: lethal-case errors by gene, no_spec
    # =========================================================================
    print("Building Figure 3 (lethal-case errors)...")
    by_gene = defaultdict(lambda: {"errors": 0, "n": 0})
    for r in parsed:
        if r["cond"] != "no_spec": continue
        c = cbi[r["tc"]]
        if "lethal" not in c["gt_drug"].lower(): continue
        by_gene[c["gene"]]["n"] += 1
        if r["scores"].get("A3", 1.0) < 1.0:
            by_gene[c["gene"]]["errors"] += 1
    sorted_genes = sorted(by_gene.items(), key=lambda x: -x[1]["errors"])
    genes = [g for g, _ in sorted_genes]
    err_counts = [v["errors"] for _, v in sorted_genes]
    err_pcts = [100 * v["errors"] / v["n"] for _, v in sorted_genes]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.barh(range(len(genes)), err_counts, color="#c44d4d")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes)
    ax.invert_yaxis()
    ax.set_xlabel("A3 errors (lethal-case clinical safety failures, no_spec)")
    for i, (n, p) in enumerate(zip(err_counts, err_pcts)):
        ax.text(n + 1, i, f"{n} ({p:.0f}%)", va="center", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, max(err_counts) * 1.2)
    save_fig(fig, "Figure3_v3_lethal_by_gene")

    # =========================================================================
    # FIGURE 4: A1 by population × model, no_spec
    # =========================================================================
    print("Building Figure 4 (population × model)...")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    pop_data = {pop: [] for pop in pops}
    for m in models:
        for pop in pops:
            rs = [r for r in parsed if r["model"] == m and r["cond"] == "no_spec" and r["pop"] == pop]
            pct = 100 * sum(1 for r in rs if r["scores"]["A1"] == 1.0) / len(rs) if rs else 0
            pop_data[pop].append(pct)

    x = np.arange(len(models))
    width = 0.27
    colors = {"EUR": "#3a6ea5", "AMR": "#c08a3a", "AFR": "#5a8f5a"}
    for i, pop in enumerate(pops):
        ax.bar(x + (i - 1) * width, pop_data[pop], width, label=pop, color=colors[pop])
    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=35, ha="right")
    ax.set_ylabel("Phenotype accuracy A1 (%, no_spec)")
    ax.set_ylim(0, 100)
    ax.legend(title="Population", loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    save_fig(fig, "Figure4_v3_population")

    print("\nAll figures written to:", FIGDIR)


if __name__ == "__main__":
    main()
