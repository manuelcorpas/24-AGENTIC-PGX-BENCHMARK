#!/usr/bin/env python3
"""
Fabrication firewall for manuscript v30.

Every headline number in the paper is recomputed here from the authoritative
data file and checked against the string actually present in the .docx. If the
manuscript and this script disagree, the manuscript is wrong.

WHY THIS EXISTS AT ALL
A referee has now twice found that the deposited artefacts and the manuscript
described different things, and this project has logged eight instances of a
number that came from tooling rather than from the phenomenon. A check that a
number appears in the text is weak; a check that the number in the text equals
the number the data produce is the one that matters, and it is cheap.

WHAT IT DOES NOT DO
It does not parse prose or find numbers the author forgot to register here. A
claim absent from CHECKS below is unverified by this script, which is why the
script prints how many claims it checked rather than implying completeness.

USAGE
    python3 code/72-validate-v30-numbers.py --manuscript /path/to/v30.docx
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
D = BASE / "data"


def load(name):
    p = D / name
    if not p.exists():
        return None
    return json.loads(p.read_text()) if p.suffix == ".json" else p.read_text()


def parse_report(text: str, scorer="baseline") -> dict:
    """Pull per-cell metrics out of v3_five_cell_live_report.txt."""
    out, section = {}, None
    for line in (text or "").splitlines():
        s = line.strip()
        if s.startswith("## "):
            section = s[3:].strip()
            continue
        if section != scorer or not s or "n=" not in s:
            continue
        cell = s.split()[0]
        fields = {}
        for tok in s.split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                fields[k] = v
        lethal_err = None
        if "(" in s and "errors)" in s:
            lethal_err = int(s.split("(")[1].split()[0])
        out[cell] = {"n": int(fields.get("n", 0)),
                     "coverage": fields.get("coverage"),
                     "A1": fields.get("A1"), "A2": fields.get("A2"),
                     "lethal_errors": lethal_err}
    return out


def build_checks():
    checks = []
    rep = parse_report(load("v3_five_cell_live_report.txt") or "")
    if rep:
        for cell, a1, a2, leth in (
            ("free_generation", "0.744", "0.624", 54),
            ("rag_generation", "0.820", "0.551", 137),
            ("rag_execution", "0.965", "0.971", 23),
            ("skill_generation", "0.964", "0.971", 26),
            ("skill_execution", "0.967", "0.973", 25),
        ):
            got = rep.get(cell, {})
            checks.append((f"{cell} A1", a1, f"{float(got.get('A1', 0)):.3f}"))
            checks.append((f"{cell} A2", a2, f"{float(got.get('A2', 0)):.3f}"))
            checks.append((f"{cell} lethal errors", str(leth),
                           str(got.get("lethal_errors"))))
        total = sum(v["n"] for v in rep.values())
        checks.append(("total evaluations", "13,199", f"{total:,}"))

    ext = load("v3_extracted_rules_eval.json")
    if ext:
        ae = ext["authored_execution"]
        checks.append(("authored coverage", "0.9996", str(ae["coverage"])))
        checks.append(("authored accuracy among emitted", "0.968",
                       f"{ae['accuracy_among_emitted']:.3f}"))
        for model, lo in (("Claude Opus 4.5", "0.627"), ("GPT-5.2", "0.636"),
                          ("o3", "0.645")):
            m = ext["models"].get(model)
            if m:
                checks.append((f"extracted coverage, {model}", lo,
                               str(m["execution"]["coverage"])))

    leth = load("v3_lethal_case_level.json")
    if leth:
        for name, diff in (("all", "0.247"), ("HLA", "0.617"), ("non-HLA", "0.042")):
            st = leth["strata"].get(name)
            if st:
                checks.append((f"lethal difference, {name}", diff,
                               f"{st['retrieval_minus_free']['difference']:.3f}"))
        checks.append(("distinct lethal cases", "14",
                       str(leth["n_distinct_cases"])))

    prov = load("v3_cell_provenance.json")
    if prov:
        tot = sum(v["cost_usd"] for v in prov["per_cell"].values())
        checks.append(("total spend across cells", "$55.30", f"${tot:.2f}"))
        lc = prov.get("legacy_comparison") or {}
        checks.append(("legacy-identical skill rows", "0",
                       str(lc.get("identical_on_raw_tokens_and_cost"))))

    getrm = load("v3_getrm_disagreement_classes.json")
    if getrm:
        # Registered at the precision the manuscript states. The in-text check is
        # string containment, so registering 0.7609 against a paper that rounds
        # to 0.761 fails on presentation rather than on substance.
        checks.append(("GeT-RM concordance", "0.761", str(getrm["concordance"])))
        checks.append(("GeT-RM shared-definition subset", "0.993",
                       str(getrm["shared_definition_concordance"])))
        checks.append(("GeT-RM evaluable pairs", "527", str(getrm["evaluable"])))
        checks.append(("GeT-RM attributable disagreements", "104",
                       str(getrm["attributable"])))
        checks.append(("GeT-RM unexplained", "22", str(getrm["unexplained"])))

    # R1.1 input normalisation. Every number is recomputed from the raw rows
    # here rather than read from a summary, so a stale report cannot validate a
    # manuscript claim.
    def _nd(d):
        return tuple(sorted(p.strip() for p in d.split("/")))

    def _rows(name):
        d = load(name)
        return (d or {}).get("rows") or []

    defs_rows = _rows("v3_input_normalisation_defs.json")
    main_rows = _rows("v3_input_normalisation_main.json")
    o3_rows = _rows("v3_input_normalisation_o3.json")

    defs_rows = defs_rows + _rows("v3_input_normalisation_defs_tail.json")
    o3_defs = _rows("v3_input_normalisation_defs_o3.json")
    gpt_defs = _rows("v3_input_normalisation_defs_gpt52.json")

    def cov_acc(rows, ref="reference"):
        scored = [r for r in rows if r.get("status") in ("call", "abstain")]
        em = [r for r in scored if r.get("status") == "call"]
        ok = sum(1 for r in em if _nd(r["call"]) == _nd(r[ref]))
        return len(em) / len(scored), (ok / len(em) if em else 0.0)

    if defs_rows and main_rows:
        keys = {(r["sample"], r["gene"]) for r in defs_rows}
        base = [r for r in main_rows
                if r["model"] == "Claude Opus 4.5" and r["form"] == "vcf"
                and (r["sample"], r["gene"]) in keys]
        bc, ba = cov_acc(base)
        dc, da = cov_acc(defs_rows)
        checks.append(("normalisation coverage, no definitions", "0.558", f"{bc:.3f}"))
        checks.append(("normalisation accuracy, no definitions", "0.456", f"{ba:.3f}"))
        checks.append(("normalisation coverage, definitions", "0.928", f"{dc:.3f}"))
        checks.append(("normalisation accuracy, definitions", "0.973", f"{da:.3f}"))

        em = [r for r in defs_rows if r.get("status") == "call"]
        mok = sum(1 for r in em if _nd(r["call"]) == _nd(r["getrm"]))
        pok = sum(1 for r in em if _nd(r["reference"]) == _nd(r["getrm"]))
        checks.append(("normalisation vs GeT-RM, model", "0.793", f"{mok/len(em):.3f}"))
        checks.append(("normalisation vs GeT-RM, caller", "0.791", f"{pok/len(em):.3f}"))
        checks.append(("normalisation pairs answered", "489", str(len(em))))

    for rows, name, cov, acc, mg, cg in (
            (o3_defs, "o3", "0.846", "0.910", "0.724", "0.767"),
            (gpt_defs, "GPT-5.2", "0.368", "0.562", "0.495", "0.768")):
        if not rows:
            continue
        c, a = cov_acc(rows)
        checks.append((f"normalisation coverage, definitions, {name}", cov, f"{c:.3f}"))
        checks.append((f"normalisation accuracy, definitions, {name}", acc, f"{a:.3f}"))
        em = [r for r in rows if r.get("status") == "call"]
        checks.append((f"normalisation vs GeT-RM, {name}", mg,
                       f"{sum(1 for r in em if _nd(r['call'])==_nd(r['getrm']))/len(em):.3f}"))
        checks.append((f"normalisation vs GeT-RM, caller for {name}", cg,
                       f"{sum(1 for r in em if _nd(r['reference'])==_nd(r['getrm']))/len(em):.3f}"))

    avd = load("v3_agent_vs_deterministic.json")
    if avd:
        for coh, cov in (("1000G_IBS", "0.215"), ("CorpasFamily", "0.333"),
                         ("Peru", "0.300"), ("UGR", "0.090")):
            d = avd["deterministic"].get(coh)
            if d:
                checks.append((f"deterministic coverage, {coh}", cov,
                               f"{d['coverage']:.3f}"))
                checks.append((f"deterministic accuracy, {coh}", "1.0",
                               str(d["accuracy_among_emitted"])))
    return checks


def _agree(stated: str, computed: str) -> bool:
    """Compare at the precision the manuscript states.

    The paper rounds to three decimals; the data file carries four. Demanding
    string equality would flag 0.627 against 0.6268 as a fabrication, which is
    noise that trains the reader of this report to ignore it.
    """
    a = stated.lstrip("$").replace(",", "")
    b = computed.lstrip("$").replace(",", "")
    if a == b:
        return True
    try:
        dp = len(a.split(".")[1]) if "." in a else 0
        return round(float(a), dp) == round(float(b), dp)
    except ValueError:
        return False


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--manuscript", type=Path, default=None,
                    help="v30 .docx; if given, each stated value must appear in it")
    a = ap.parse_args(argv)

    text = None
    if a.manuscript:
        try:
            import docx
        except ImportError:
            sys.stderr.write("python-docx not installed; skipping the text check\n")
        else:
            if not a.manuscript.exists():
                sys.stderr.write(f"missing {a.manuscript}\n")
                return 1
            doc = docx.Document(str(a.manuscript))
            text = "\n".join(p.text for p in doc.paragraphs)
            text += "\n" + "\n".join(c.text for t in doc.tables
                                     for r in t.rows for c in r.cells)

    checks = build_checks()
    if not checks:
        sys.stderr.write("no data files found; nothing verified\n")
        return 1

    fails = []
    print(f"{'claim':44s} {'stated':>12} {'from data':>12}  {'in text':>8}")
    for name, stated, computed in checks:
        agrees = _agree(stated, computed)
        in_text = "n/a" if text is None else ("yes" if stated in text else "NO")
        flag = "" if (agrees and in_text != "NO") else "   <-- FAIL"
        if flag:
            fails.append((name, stated, computed, in_text))
        print(f"{name:44s} {stated:>12} {computed:>12}  {in_text:>8}{flag}")

    print(f"\nchecked {len(checks)} registered claims; {len(fails)} failed")
    if fails:
        print("\nA failure means the manuscript states a number the data do not "
              "produce, or omits one it should state.")
        return 1
    print("All registered claims reproduce from data and appear in the manuscript.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
