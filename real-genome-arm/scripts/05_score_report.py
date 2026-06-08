#!/usr/bin/env python3
"""
Step 5 — score the agent's real-genome predictions against the cohort's CPIC
phenotype calls, and decompose the result the way the manuscript reports it:

  - DETERMINATE diplotypes (cohort gave a definite phenotype): exact-match accuracy
  - UNCERTAIN/Indeterminate diplotypes (abstaining is the correct answer):
    abstention rate (did the agent output Indeterminate rather than guess?)
  - GENUINE errors by gene (after removing terminology/format-only mismatches)

Self-contained phenotype-tier matcher (no external scorer needed).

Usage: 05_score_report.py <predictions.tsv> [out_report.txt]
"""
import sys
import re
import csv
from collections import defaultdict

TIERS = [
    ("MH", r"malignant hyperthermia|mh[\s-]?susceptible|\bmhs\b"),
    ("ULTRARAPID", r"ultra[\s-]?rapid"),
    ("RAPID", r"\brapid metaboli"),
    ("POOR", r"\bpoor metaboli"),
    ("INTERMEDIATE", r"intermediate metaboli"),
    ("NORMAL_M", r"normal metaboli|extensive metaboli"),
    ("POOR_FN", r"poor function|low function"),
    ("DECREASED_FN", r"decreased function|reduced function|intermediate function|possible decreased"),
    ("INCREASED_FN", r"increased function|possible increased"),
    ("NORMAL_FN", r"normal function"),
    ("POS", r"\b(positive|carrier|at[\s-]?risk|present|susceptib)"),
    ("NEG", r"\b(negative|non[\s-]?carrier|absent|not at risk)"),
    ("DEFICIENT", r"deficien"),
    ("FAVORABLE", r"favou?rable response|favou?rable"),
    ("UNFAVORABLE", r"unfavou?rable response|poor response"),
    ("INDETERMINATE", r"indetermin|uncertain|unknown|variable|n/a"),
]

def tier(s):
    s = (s or "").lower()
    for name, pat in TIERS:
        if re.search(pat, s):
            return name
    return None

def is_uncertain(s):
    return tier(s) == "INDETERMINATE"

def main():
    src = sys.argv[1]
    rows = list(csv.DictReader(open(src), delimiter="\t"))
    # per (cohort, model)
    agg = defaultdict(lambda: {"det_n": 0, "det_ok": 0, "unc_n": 0, "abst": 0})
    err_by_gene = defaultdict(int)
    for r in rows:
        coh, model = r["cohort"], r["model"]
        gt, pred = r["cohort_phenotype"], r["pred"]
        k = (coh, model)
        if is_uncertain(gt):
            agg[k]["unc_n"] += 1
            if is_uncertain(pred):
                agg[k]["abst"] += 1
        else:
            agg[k]["det_n"] += 1
            match = tier(pred) is not None and tier(pred) == tier(gt)
            if match:
                agg[k]["det_ok"] += 1
            else:
                a, b = (pred or "").lower(), (gt or "").lower()
                if not (a[:12] and (a[:12] in b or b[:12] in a)):  # not a format-only artefact
                    err_by_gene[(coh, r["gene"])] += 1

    out = ["REAL-GENOME ARM — agent vs cohort CPIC phenotype (PyPGx-called diplotypes)", ""]
    out.append(f"{'cohort':10}{'model':18}{'determinate':>14}{'abstain-on-uncertain':>22}")
    by_coh = defaultdict(lambda: {"det_n": 0, "det_ok": 0, "unc_n": 0, "abst": 0})
    for (coh, model), s in sorted(agg.items()):
        da = f"{s['det_ok']}/{s['det_n']}" + (f" ({100*s['det_ok']/s['det_n']:.0f}%)" if s['det_n'] else "")
        ab = f"{s['abst']}/{s['unc_n']}" + (f" ({100*s['abst']/s['unc_n']:.0f}%)" if s['unc_n'] else "")
        out.append(f"{coh:10}{model:18}{da:>14}{ab:>22}")
        for kk in s:
            by_coh[coh][kk] += s[kk]
    out.append("")
    out.append("=== PER-COHORT (panel pooled) ===")
    for coh, s in sorted(by_coh.items()):
        det = f"{100*s['det_ok']/s['det_n']:.1f}%" if s['det_n'] else "n/a"
        ab = f"{100*s['abst']/s['unc_n']:.1f}%" if s['unc_n'] else "n/a"
        out.append(f"  {coh:8} determinate accuracy {det}   abstention {ab}   (det n={s['det_n']}, uncertain n={s['unc_n']})")
    out.append("")
    out.append("=== GENUINE ERRORS BY (cohort, gene) ===")
    for (coh, gene), n in sorted(err_by_gene.items(), key=lambda x: -x[1])[:25]:
        out.append(f"  {coh:8} {gene:9} {n}")
    text = "\n".join(out)
    print(text)
    if len(sys.argv) > 2:
        open(sys.argv[2], "w").write(text)
        print(f"\nwrote {sys.argv[2]}")

if __name__ == "__main__":
    main()
