#!/usr/bin/env python3
"""
Diff and execute model-extracted rule tables (reviewer point 2).

WHAT THIS ANSWERS
The manuscript claimed that extracting decision logic from guideline prose
reaches the same place as a hand-authored rule table. It never tested that: both
execution cells were handed the authored controlled vocabulary and both ran the
same executor, so their agreement was guaranteed by construction. The reviewer
asked for the real experiment: extract a table, freeze it, diff it against the
authored table, and execute the extracted table. This is the second half.

THE ONE RULE THAT MAKES THE MEASUREMENT HONEST
If a called diplotype is absent from the extracted table, the pipeline ABSTAINS.
It never falls back to the authored mapping. A fallback would measure the
authored table twice and rebuild exactly the circularity being corrected, while
looking like a higher score.

Abstention is therefore reported beside accuracy, not folded into it. An
extracted table that covers a third of the vocabulary perfectly is not as good
as one that covers all of it, and a single accuracy figure would hide that.

NO NEW MODEL CALLS
Execution reuses the diplotypes the panel already called in the matched
factorial's skill_execution cell. Only the mapping is swapped, which is the
variable under test.

USAGE
    python code/68-evaluate-extracted-rules.py
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
EXTRACTED = BASE / "data" / "v3_extracted_rules.json"
LIVE = BASE / "data" / "v3_five_cell_live.json"
OUT = BASE / "data" / "v3_extracted_rules_eval.json"
REPORT = BASE / "data" / "v3_extracted_rules_eval.txt"


def _load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rules = _load(CODE / "_pgx_rules.py", "_pgx_rules")

# Phenotype wording differs between the guideline and the authored table
# ("CYP2C19 normal metabolizer" against "Normal Metabolizer"). Comparing raw
# strings would measure capitalisation and gene prefixes, not agreement.
_TIERS = [
    ("ultrarapid", "ultrarapid"), ("rapid", "rapid"),
    ("normal", "normal"), ("extensive", "normal"),
    ("intermediate", "intermediate"), ("poor", "poor"),
    ("indeterminate", "indeterminate"), ("uncertain", "indeterminate"),
    ("increased", "increased"), ("decreased", "decreased"),
    ("normal function", "normal"), ("no function", "poor"),
    ("positive", "positive"), ("negative", "negative"),
]


def phen_key(text: str) -> str:
    t = (text or "").strip().lower()
    if not t:
        return ""
    likely = "likely" in t
    for needle, tier in _TIERS:
        if needle in t:
            return ("likely " if likely else "") + tier
    return re.sub(r"[^a-z]+", " ", t).strip()


def dip_key(gene: str, dip: str) -> str:
    return rules.norm_dip((dip or "").strip(), gene)


def diff_gene(gene: str, authored: dict, extracted: dict) -> dict:
    ext = {dip_key(gene, k): v for k, v in
           (extracted.get("diplotype_to_phenotype") or {}).items()}
    auth = {dip_key(gene, k): v for k, v in authored.items()}
    shared = sorted(set(auth) & set(ext))
    agree, examples = 0, []
    for d in shared:
        if phen_key(auth[d]) == phen_key(ext[d]):
            agree += 1
        elif len(examples) < 5:
            examples.append({"diplotype": d, "authored": auth[d], "extracted": ext[d]})
    return {
        "gene": gene,
        "authored_n": len(auth),
        "extracted_n": len(ext),
        "shared": len(shared),
        "agree": agree,
        "missing": sorted(set(auth) - set(ext)),
        "extra": sorted(set(ext) - set(auth)),
        "disagree_examples": examples,
    }


def execute_with(gene: str, called: str, extracted: dict) -> str | None:
    """Map a called diplotype using ONLY the extracted table, or abstain."""
    if not called:
        return None
    ext = {dip_key(gene, k): v for k, v in
           (extracted.get("diplotype_to_phenotype") or {}).items()}
    return ext.get(dip_key(gene, called))


def score_rows(rows: list[dict], tables: dict) -> dict:
    n = correct = wrong = abstained = 0
    per_gene = defaultdict(lambda: {"n": 0, "correct": 0, "abstained": 0})
    for r in rows:
        gene = r["gene"]
        n += 1
        per_gene[gene]["n"] += 1
        table = tables.get(gene)
        got = execute_with(gene, r.get("called_diplotype", ""), table) if table else None
        if got is None:
            abstained += 1
            per_gene[gene]["abstained"] += 1
            continue
        if phen_key(got) == phen_key(r["gt_phenotype"]):
            correct += 1
            per_gene[gene]["correct"] += 1
        else:
            wrong += 1
    emitted = correct + wrong
    return {
        "n": n, "correct": correct, "wrong": wrong, "abstained": abstained,
        "coverage": round(emitted / n, 4) if n else None,
        "accuracy_among_emitted": round(correct / emitted, 4) if emitted else None,
        "accuracy_overall": round(correct / n, 4) if n else None,
        "per_gene": {g: dict(v) for g, v in sorted(per_gene.items())},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--extracted", type=Path, default=EXTRACTED)
    ap.add_argument("--live", type=Path, default=LIVE)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--report", type=Path, default=REPORT)
    a = ap.parse_args(argv)

    if not a.extracted.exists():
        sys.stderr.write(f"missing {a.extracted}; run code/67-extract-rules-from-prose.py\n")
        return 1

    cases = rules.load_cases()
    _, _, rules_dip, _ = rules.build_rules(cases)
    case_by_id = {c["id"]: c for c in cases}

    extracted_rows = json.loads(a.extracted.read_text())
    by_model = defaultdict(dict)
    for r in extracted_rows:
        if not r.get("error"):
            by_model[r["model"]][r["gene"]] = r["table"]

    live = json.loads(a.live.read_text())
    exec_rows = []
    for r in live:
        if r.get("cell") != "skill_execution":
            continue
        c = case_by_id[r["case_id"]]
        exec_rows.append({"gene": c["gene"],
                          "called_diplotype": r.get("called_diplotype") or "",
                          "gt_phenotype": c["gt_phenotype"]})

    result = {"models": {}, "authored_totals": {g: len(v) for g, v in rules_dip.items()}}
    for model, tables in sorted(by_model.items()):
        diffs = [diff_gene(g, rules_dip.get(g, {}), t) for g, t in sorted(tables.items())]
        result["models"][model] = {
            "n_genes": len(tables),
            "diff": diffs,
            "diff_totals": {
                "authored": sum(d["authored_n"] for d in diffs),
                "extracted": sum(d["extracted_n"] for d in diffs),
                "shared": sum(d["shared"] for d in diffs),
                "agree": sum(d["agree"] for d in diffs),
                "missing": sum(len(d["missing"]) for d in diffs),
                "extra": sum(len(d["extra"]) for d in diffs),
            },
            "execution": score_rows(exec_rows, tables),
        }

    # The authored table executed on the same calls: the comparator.
    authored_tables = {g: {"diplotype_to_phenotype": v} for g, v in rules_dip.items()}
    result["authored_execution"] = score_rows(exec_rows, authored_tables)

    a.out.write_text(json.dumps(result, indent=2))

    L = ["EXTRACTED RULE TABLES: DIFF AGAINST THE AUTHORED TABLE, AND EXECUTION", ""]
    ae = result["authored_execution"]
    L.append(f"  authored table executed on the same {ae['n']} calls:")
    L.append(f"    coverage {ae['coverage']}   accuracy among emitted "
             f"{ae['accuracy_among_emitted']}   overall {ae['accuracy_overall']}")
    L.append("")
    for model, m in result["models"].items():
        t, e = m["diff_totals"], m["execution"]
        L.append(f"  {model}  ({m['n_genes']} genes)")
        L.append(f"    diff vs authored: {t['shared']} shared diplotypes, "
                 f"{t['agree']} phenotype agreements "
                 f"({t['agree']/t['shared']:.3f})" if t["shared"] else
                 f"    diff vs authored: no shared diplotypes")
        L.append(f"                      {t['missing']} authored diplotypes absent, "
                 f"{t['extra']} extracted not in the authored table")
        L.append(f"    executed:         coverage {e['coverage']}, "
                 f"accuracy among emitted {e['accuracy_among_emitted']}, "
                 f"overall {e['accuracy_overall']}")
        L.append(f"                      {e['correct']} correct, {e['wrong']} wrong, "
                 f"{e['abstained']} abstained of {e['n']}")
        L.append("")
    L.append("  Abstention is reported beside accuracy, never folded into it: a table")
    L.append("  covering part of the vocabulary perfectly is not equivalent to one")
    L.append("  covering all of it. A called diplotype absent from an extracted table")
    L.append("  abstains and never falls back to the authored mapping.")
    text = "\n".join(L)
    a.report.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
