#!/usr/bin/env python3
"""Figure S9. One lethal-class case traced through all five configurations.

Every string shown in the figure is read from the deposited run records and the
case specification. Nothing is retyped: if a response changes in the deposit,
the figure changes. Usage:

    python3 code/93-figure-s9-worked-case.py --out FigureS9.tiff
"""
from __future__ import annotations

import argparse
import json
import textwrap
import re
from pathlib import Path

import importlib.util

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

BASE = Path(__file__).resolve().parent.parent
CASE_ID = "hlab1502_carba_pos"
MODEL = "GPT-5.2"
REP = 0

CELLS = [
    ("free_generation",   "Free generation",        "none",                   "model"),
    ("rag_generation",    "RAG generation",         "retrieved guideline prose", "model"),
    ("rag_execution",     "RAG-assisted execution", "retrieved prose (call only)", "code"),
    ("skill_generation",  "Authored-rule generation", "authored rule table",  "model"),
    ("skill_execution",   "Authored-rule execution",  "authored rule table",  "code"),
]

INK = "#111111"
MUTED = "#5b5b5b"
SAFE_FILL, SAFE_EDGE = "#e8f2ea", "#2f7d4f"
ERR_FILL, ERR_EDGE = "#fbeaea", "#b3261e"
BOX_FILL, BOX_EDGE = "#f6f6f4", "#c9c9c4"


def rule_engine(base: Path):
    """The authored mapping, loaded through the same module the runners use."""
    spec = importlib.util.spec_from_file_location(
        "_pgx_rules", base / "code" / "_pgx_rules.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    d2p, rec, _, _ = mod.build_rules(mod.load_cases(base / "specs" / "test_cases_v3.json"))
    return mod, d2p, rec


def load(base: Path):
    live = json.loads((base / "data" / "v4_flash_eight_model_five_cell.json").read_text())
    scored = json.loads((base / "data" / "v4_flash_matched_scored_rows_all5.json").read_text())
    cases = json.loads((base / "specs" / "test_cases_v3.json").read_text())
    return live, scored, cases


def collect(base: Path, case_id: str, model: str, rep: int):
    live, scored, cases = load(base)
    case = next(c for c in cases if c["id"] == case_id)
    key = lambda r: (r["case_id"], r["cell"], r["model"], r["rep"])
    sidx = {key(r): r for r in scored}
    mod, d2p, rec = rule_engine(base)
    out = []
    for cell, label, knowledge, agent in CELLS:
        raw = next(r for r in live if key(r) == (case_id, cell, model, rep))
        sc = sidx[(case_id, cell, model, rep)]
        assert sc["lethal_action"] is not None, f"{cell} is not scored as lethal-class"
        executed = None
        if agent == "code":
            nd = mod.norm_dip(raw["called_diplotype"], case["gene"])
            executed = (d2p.get((case["gene"], nd)), rec.get((case["gene"], case["drug"], nd)))
            assert executed[1], f"{cell}: authored mapping returned nothing for {nd!r}"
        out.append({
            "label": label, "knowledge": knowledge, "agent": agent,
            "text": (raw["raw"] or "").strip(),
            "executed": executed,
            "error": sc["lethal_action"] == 0.0,
            "widened_error": not bool(re.search(r'\bavoid\b|\bcontraindicat|\bdo not use\b',raw.get('parsed_drug') or '',re.I)) or sc['substituted'],
        })
    return case, out


def decision_lines(text: str) -> str:
    """The two fields that carry the clinical decision, as returned."""
    keep = []
    for line in text.splitlines():
        s = line.strip()
        if s.upper().startswith(("DIPLOTYPE:", "DRUG:")):
            keep.append(s)
    return "\n".join(keep) if keep else text.strip()


def panel_error_rates(base: Path, case_id: str):
    live,scored,_=load(base)
    raw={(r['cell'],r['case_id'],r['model'],r['rep']):r for r in live}
    gene = next(r["gene"] for r in scored if r["case_id"] == case_id)
    rows = [r for r in scored if r["gene"] == gene and r.get("lethal_action") is not None]
    out = {}
    for cell, *_ in CELLS:
        sub = [r for r in rows if r["cell"] == cell]
        widened=sum(not bool(re.search(r'\bavoid\b|\bcontraindicat|\bdo not use\b',raw[(r['cell'],r['case_id'],r['model'],r['rep'])].get('parsed_drug') or '',re.I)) or r['substituted'] for r in sub)
        out[cell] = (sum(1 for r in sub if r["lethal_action"] == 0.0), len(sub),widened)
    return gene, out


def render(case, rows, gene, rates, out_path: Path, model: str, rep: int, dpi: int = 600, save_output=True):
    """Lay out in inches from the top, then size the canvas to the content."""
    W = 7.2
    mono = {"family": "DejaVu Sans Mono"}
    LINE = 0.118          # inches per wrapped text line
    HEAD = 0.46           # inches for a row's title + meta line
    PAD = 0.10

    for r in rows:
        r["wrapped"] = "\n".join(
            textwrap.fill(ln, 104, subsequent_indent="            ")
            for ln in decision_lines(r["text"]).splitlines())
        r["nbody"] = r["wrapped"].count("\n") + 1
        r["exec_wrapped"] = None
        if r["executed"]:
            r["exec_wrapped"] = textwrap.fill(
                f"code applies the authored mapping  ->  {r['executed'][1]}",
                104, subsequent_indent="      ")
            r["nexec"] = r["exec_wrapped"].count("\n") + 1
        r["h"] = HEAD + LINE * r["nbody"] + (LINE * (r.get("nexec", 0) + 0.9) if r["exec_wrapped"] else 0) + PAD

    head_h = 0.86          # title block
    input_h = 0.86         # shared-input box
    FOOT = textwrap.fill(
        "Responses are reproduced verbatim from the deposited run records "
        "(v4_flash_eight_model_five_cell.json); only the DIPLOTYPE and DRUG fields are shown. "
        "Executed recommendations are computed by the authored mapping in _pgx_rules.py.",
        128)
    foot_h = 0.26 + 0.10 * (FOOT.count("\n") + 1)
    H = head_h + input_h + sum(r["h"] for r in rows) + foot_h

    fig = plt.figure(figsize=(W, H), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    Y = lambda inches_from_top: 1.0 - inches_from_top / H
    HH = lambda inches: inches / H

    ax.text(0.035, Y(0.30), "One high-severity case, two wording rules",
            fontsize=11.5, fontweight="bold", color=INK, va="top")
    ax.text(0.035, Y(0.53),
            f"Case {case['id']}  \u00b7  {model}, replicate {rep + 1}. "
            "Only the knowledge supplied and the decision agent differ.",
            fontsize=7.6, color=MUTED, va="top")

    y = head_h
    ax.add_patch(FancyBboxPatch((0.035, Y(y + input_h - 0.10)), 0.93, HH(input_h - 0.16),
                                boxstyle="round,pad=0.006,rounding_size=0.008",
                                facecolor=BOX_FILL, edgecolor=BOX_EDGE, linewidth=0.7))
    ax.text(0.05, Y(y + 0.16), "Identical input to all five configurations",
            fontsize=7.4, fontweight="bold", color=INK, va="top")
    ax.text(0.05, Y(y + 0.34),
            f"Gene: {case['gene']}        Genotype: {case['genotype']}        "
            f"Target drug: {case['drug']}",
            fontsize=7.0, color=INK, va="top", **mono)
    ax.text(0.05, Y(y + 0.50), f"Authored reference: {case['gt_drug']}",
            fontsize=7.0, color=INK, va="top", **mono)
    ax.text(0.05, Y(y + 0.66),
            "A single configuration-blind function scores every response.",
            fontsize=6.8, color=MUTED, va="top", style="italic")

    y = head_h + input_h
    for i, r in enumerate(rows):
        fill, edge = (ERR_FILL, ERR_EDGE) if r["error"] else (SAFE_FILL, SAFE_EDGE)
        ax.add_patch(FancyBboxPatch((0.035, Y(y + r["h"] - 0.06)), 0.93, HH(r["h"] - 0.12),
                                    boxstyle="round,pad=0.006,rounding_size=0.008",
                                    facecolor=fill, edgecolor=edge, linewidth=0.9))
        ax.text(0.05, Y(y + 0.16), f"{i + 1}.  {r['label']}",
                fontsize=8.2, fontweight="bold", color=INK, va="top")
        ax.text(0.05, Y(y + 0.34),
                f"knowledge: {r['knowledge']}    decision by: {r['agent']}",
                fontsize=6.9, color=MUTED, va="top")
        verdict = ('Original fail; widened '+('fail' if r['widened_error'] else 'pass')) if r['error'] else 'Both wording rules pass'
        ax.text(0.95, Y(y + 0.16), verdict, fontsize=7.4, fontweight="bold",
                color=edge, va="top", ha="right")
        n_err, n_tot, n_wide = rates[CELLS[i][0]]
        ax.text(0.95, Y(y + 0.34), f"Locus failures: original {n_err}/{n_tot}; widened {n_wide}/{n_tot}",
                fontsize=6.6, color=MUTED, va="top", ha="right")
        ax.text(0.05, Y(y + 0.50), r["wrapped"], fontsize=6.4, color=INK,
                va="top", linespacing=1.35, **mono)
        if r["exec_wrapped"]:
            ax.text(0.05, Y(y + 0.50 + LINE * (r["nbody"] + 0.9)), r["exec_wrapped"],
                    fontsize=6.4, color=edge, va="top", linespacing=1.35,
                    fontweight="bold", **mono)
        y += r["h"]

    ax.text(0.035, Y(H - foot_h + 0.10), FOOT,
            fontsize=6.2, color=MUTED, va="top", linespacing=1.4)
    if not save_output:return fig
    fig.savefig(out_path, dpi=dpi, format="tiff",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(out_path.with_suffix(".png"), dpi=200)
    plt.close(fig)
    print(f"  canvas {W:.2f} x {H:.2f} in")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="FigureS9.tiff")
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--model", default=MODEL)
    ap.add_argument("--rep", type=int, default=REP)
    ap.add_argument("--dpi", type=int, default=600)
    a = ap.parse_args()
    base = Path(a.base)
    case, rows = collect(base, CASE_ID, a.model, a.rep)
    gene, rates = panel_error_rates(base, CASE_ID)
    render(case, rows, gene, rates, Path(a.out), a.model, a.rep, a.dpi)
    print(f"wrote {a.out}")
    for r in rows:
        print(f"  {r['label']:26s} {'ERROR' if r['error'] else 'ok'}")


if __name__ == "__main__":
    main()
