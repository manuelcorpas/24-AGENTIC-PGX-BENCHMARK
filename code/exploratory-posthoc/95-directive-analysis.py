#!/usr/bin/env python3
# EXPLORATORY, POST HOC, NOT CLINICALLY VALIDATED.
# Written after the primary results were seen. Not part of the frozen scoring
# pipeline and not an endpoint. See README.md in this directory.
"""Directive-strength analysis of the five matched configurations.

Reads the deposited run records and reports, per configuration:
  (1) the distribution of responses across the directive ladder
  (2) unresolved conditionals: prohibitions gated on information not supplied
  (3) branch commitment where the ground truth permits dose modification
  (4) transfer of directive register from the CPIC source text to the output
  (5) lethal-class errors under the lexical scorer against the speech-act scorer
"""
from __future__ import annotations

import argparse
import collections
import importlib.util
import json
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
CELLS = ["free_generation", "rag_generation", "rag_execution",
         "skill_generation", "skill_execution"]
LABEL = {"free_generation": "free generation", "rag_generation": "RAG generation",
         "rag_execution": "RAG-assisted execution",
         "skill_generation": "authored-rule generation",
         "skill_execution": "authored-rule execution"}
HLA = {"hlaa3101_carba_pos", "hlab1502_carba_pos", "hlab1502_oxcarba_pos",
       "hlab5701_aba_pos", "hlab5801_allo_pos"}


def load(base: Path):
    spec = importlib.util.spec_from_file_location("ds", base / "code" / "exploratory-posthoc" / "94-directive-strength.py")
    ds = importlib.util.module_from_spec(spec); spec.loader.exec_module(ds)
    live = json.loads((base / "data" / "v3_five_cell_live.json").read_text())
    scored = json.loads((base / "data" / "v3_matched_scored_rows_all5.json").read_text())
    cases = {c["id"]: c for c in json.loads((base / "specs" / "test_cases_v3.json").read_text())}
    corpus = json.loads((base / "specs" / "cpic_rag_corpus_v3.json").read_text())["genes"]
    return ds, live, scored, cases, corpus


def build(base: Path):
    ds, live, scored, cases, corpus = load(base)
    key = lambda r: (r["case_id"], r["cell"], r["model"], r["rep"])
    LIVE = {key(r): r for r in live}
    rows = []
    for s in scored:
        if s.get("lethal_action") is None:
            continue
        raw = LIVE.get(key(s))
        text = (raw or {}).get("parsed_drug") or ""
        case = cases[s["case_id"]]
        rows.append({
            "cell": s["cell"], "case_id": s["case_id"], "gene": s["gene"],
            "drug": s["drug"], "model": s["model"], "rep": s["rep"],
            "text": text,
            "level": ds.directive_strength(text),
            "outcome": ds.classify(case, text),
            "lexical_ok": s["lethal_action"] == 1.0,
            "hla": s["case_id"] in HLA,
            "permits_dose": ds.permits_dose_modification(case),
        })
    return ds, rows, cases, corpus


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--base", default=str(BASE))
    a = ap.parse_args(); base = Path(a.base)
    ds, rows, cases, corpus = build(base)
    N = {c: sum(1 for r in rows if r["cell"] == c) for c in CELLS}
    names = {ds.CATEGORICAL: "categorical", ds.CONDITIONAL: "conditional",
             ds.AGAINST: "against", ds.MODIFY: "modify", ds.NONE: "none"}

    print("=" * 78)
    print("(1) DIRECTIVE LADDER, % of lethal-class responses")
    print(f"{'configuration':26s}" + "".join(f"{names[l][:11]:>12s}" for l in [4, 3, 2, 1, 0]))
    for c in CELLS:
        sub = [r for r in rows if r["cell"] == c]
        print(f"{LABEL[c]:26s}" + "".join(
            f"{100*sum(1 for r in sub if r['level']==l)/len(sub):11.1f}%" for l in [4, 3, 2, 1, 0]))

    print()
    print("(2) UNRESOLVED CONDITIONALS  (prohibition gated on information the case never supplies)")
    print(f"{'configuration':26s}{'n':>7s}{'% of responses':>17s}{'  on HLA cases':>16s}")
    for c in CELLS:
        sub = [r for r in rows if r["cell"] == c]
        u = [r for r in sub if r["outcome"] == ds.UNRESOLVED]
        uh = [r for r in u if r["hla"]]
        print(f"{LABEL[c]:26s}{len(u):7d}{100*len(u)/len(sub):16.1f}%{len(uh):16d}")

    print()
    print("(3) BRANCH COMMITMENT, cases whose ground truth permits dose modification")
    for c in CELLS:
        sub = [r for r in rows if r["cell"] == c and r["permits_dose"]]
        if not sub: continue
        dose = sum(1 for r in sub if r["level"] == ds.MODIFY)
        print(f"   {LABEL[c]:26s} n={len(sub):3d}  committed to dose branch {dose:3d} "
              f"({100*dose/len(sub):4.1f}%)   lexical scorer marks these wrong: "
              f"{sum(1 for r in sub if r['level']==ds.MODIFY and not r['lexical_ok'])}")

    print()
    print("(4) REGISTER TRANSFER: directive strength of the CPIC source vs the output")
    print(f"{'locus':16s}{'source':>9s}" + "".join(f"{LABEL[c][:11]:>12s}" for c in CELLS))
    loci = sorted({r["gene"] for r in rows})
    for g in loci:
        sub0 = [r for r in rows if r["gene"] == g]
        drug = sub0[0]["drug"]
        ex = (corpus.get(g) or {}).get("guideline_excerpt", "")
        src = ds.source_strength(ex, drug)
        line = f"{g:16s}{names[src][:8]:>9s}"
        for c in CELLS:
            sub = [r for r in sub0 if r["cell"] == c]
            m = sum(r["level"] for r in sub) / len(sub) if sub else 0
            line += f"{m:12.2f}"
        print(line)
    print("   (source column is the ladder level of the guideline text itself; "
          "cell columns are mean output level)")

    print()
    print("(5) LETHAL-CLASS ERRORS: lexical scorer vs speech-act scorer")
    print(f"{'configuration':26s}{'lexical':>9s}{'non-avoid':>11s}{'unresolved':>12s}{'executable':>12s}")
    for c in CELLS:
        sub = [r for r in rows if r["cell"] == c]
        lex = sum(1 for r in sub if not r["lexical_ok"])
        na = sum(1 for r in sub if r["outcome"] == ds.NON_AVOIDANCE)
        un = sum(1 for r in sub if r["outcome"] == ds.UNRESOLVED)
        ex = sum(1 for r in sub if r["outcome"] == ds.EXECUTABLE)
        print(f"{LABEL[c]:26s}{lex:9d}{na:11d}{un:12d}{ex:12d}")
    print()
    print("   HLA subset (5 cases, n=120 per configuration):")
    for c in CELLS:
        sub = [r for r in rows if r["cell"] == c and r["hla"]]
        print(f"      {LABEL[c]:26s} non-avoid {sum(1 for r in sub if r['outcome']==ds.NON_AVOIDANCE):3d}"
              f"   unresolved {sum(1 for r in sub if r['outcome']==ds.UNRESOLVED):3d}"
              f"   executable {sum(1 for r in sub if r['outcome']==ds.EXECUTABLE):3d}")


if __name__ == "__main__":
    main()
