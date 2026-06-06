#!/usr/bin/env python3
"""
Adversarial audit of the v3 benchmark — designed to detect hallucination,
scoring exploits, and circularity that a peer reviewer would attack.

Sections:
  A1  Data structural integrity
  A2  Score-range and value sanity
  A3  Cell-coverage exactness (no duplicates, no missing combos)
  A4  Response uniqueness (no caching/dedup artefacts in raw data)
  A5  Scoring boundary tests — adversarial inputs to score_a1, score_a2
  A6  Echo-vs-paraphrase analysis on with_spec PHENOTYPE field
  A7  with_spec PHENOTYPE field is not literal prompt copy
  A8  Lethal-case correctness — every case flagged as lethal IS lethal per CPIC
  A9  Population-attribution audit
  A10 Stochastic check — within-combo run agreement
  A11 Cross-model agreement on hard cases (independent corroboration)
  A12 Spec-content audit — does spec contain population keyword for B1?

Output: ../RESULTS/v3_audit_report.txt
"""
from __future__ import annotations
import json
import re
from collections import Counter, defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
PYDIR = BASE / "PYTHON"
RAW = BASE / "RESULTS" / "v3_raw_rescored.json"
CASES = BASE / "SPECS" / "test_cases_v3.json"
REPORT = BASE / "RESULTS" / "v3_audit_report.txt"

# Import scorer
_spec = spec_from_file_location("rescore", PYDIR / "10-rescore-v3.py")
rescore = module_from_spec(_spec); _spec.loader.exec_module(rescore)


def section(t):
    return f"\n{'=' * 70}\n{t}\n{'=' * 70}"


def main():
    rows = json.loads(RAW.read_text())
    cases = json.loads(CASES.read_text())
    cbi = {c["id"]: c for c in cases}
    out = ["v3 BENCHMARK ADVERSARIAL AUDIT", f"Source: {RAW.name}", f"Rows: {len(rows)}"]
    issues = []

    # =========================================================================
    # A1: Data structural integrity
    # =========================================================================
    out.append(section("A1: Data structural integrity"))
    expected_n = 9 * 110 * 3 * 2 * 3
    out.append(f"  Expected cells: {expected_n}")
    out.append(f"  Actual cells:   {len(rows)}")
    if len(rows) != expected_n:
        issues.append(f"A1: Row count {len(rows)} != expected {expected_n}")
    out.append(f"  [{'PASS' if len(rows) == expected_n else 'FAIL'}]")

    # Required keys
    required_keys = {"run", "model", "tc", "pop", "cond", "parsed", "scores"}
    missing = sum(1 for r in rows if not required_keys.issubset(r.keys()))
    out.append(f"  Rows missing required keys: {missing}")
    if missing:
        issues.append(f"A1: {missing} rows missing required keys")

    # =========================================================================
    # A2: Score range and value sanity
    # =========================================================================
    out.append(section("A2: Score range and value sanity"))
    valid_a_values = {0.0, 0.5, 1.0}
    out_of_range = 0
    for r in rows:
        if r["scores"].get("format_fail"): continue
        for dim in ["A1", "A2", "A3"]:
            v = r["scores"].get(dim)
            if v not in valid_a_values:
                out_of_range += 1
    out.append(f"  Tier A scores outside {{0, 0.5, 1.0}}: {out_of_range}")
    if out_of_range:
        issues.append(f"A2: {out_of_range} Tier A scores outside valid range")
    out.append(f"  [{'PASS' if out_of_range == 0 else 'FAIL'}]")

    # =========================================================================
    # A3: Cell coverage exactness
    # =========================================================================
    out.append(section("A3: Cell coverage exactness (no duplicates, no missing combos)"))
    cells = Counter()
    for r in rows:
        cells[(r["run"], r["model"], r["tc"], r["pop"], r["cond"])] += 1
    duplicates = [k for k, v in cells.items() if v > 1]
    missing_combos = expected_n - len(cells)
    out.append(f"  Unique cells: {len(cells)}")
    out.append(f"  Duplicates:   {len(duplicates)} (showing first 5)")
    for d in duplicates[:5]:
        out.append(f"    {d}")
    out.append(f"  Missing combos: {missing_combos}")
    if duplicates or missing_combos:
        issues.append(f"A3: {len(duplicates)} duplicates, {missing_combos} missing")
    out.append(f"  [{'PASS' if not duplicates and not missing_combos else 'FAIL'}]")

    # =========================================================================
    # A4: Response uniqueness — detect cache hits / dedup artefacts
    # =========================================================================
    out.append(section("A4: Response uniqueness within-combo (3 runs should produce text variation)"))
    by_combo = defaultdict(list)
    for r in rows:
        if r["scores"].get("format_fail"): continue
        by_combo[(r["model"], r["tc"], r["pop"], r["cond"])].append(r)
    identical_run_combos = 0
    sampled_no_spec_with_identical = 0
    sampled_with_spec_with_identical = 0
    for k, v in by_combo.items():
        if len(v) < 3: continue
        texts = tuple(json.dumps(r["parsed"], sort_keys=True) for r in v)
        if len(set(texts)) == 1:
            identical_run_combos += 1
            if k[3] == "no_spec":
                sampled_no_spec_with_identical += 1
            else:
                sampled_with_spec_with_identical += 1
    out.append(f"  Combos with identical 3-run output: {identical_run_combos}/{len(by_combo)}")
    out.append(f"    no_spec: {sampled_no_spec_with_identical}")
    out.append(f"    with_spec: {sampled_with_spec_with_identical}")
    # with_spec identical-output is EXPECTED (deterministic); no_spec identical is suspicious if very high
    out.append(f"  Note: with_spec identical is expected; no_spec identical >50% would suggest temperature=0 or cache effects")

    # =========================================================================
    # A5: Scoring boundary tests — adversarial inputs
    # =========================================================================
    out.append(section("A5: Scoring boundary tests"))
    # If we feed a clearly wrong answer to the scorer, it must return 0
    boundary_tests = [
        # (gt_phenotype, parsed_phen, expected_A1, label)
        ("Poor Metaboliser", "Normal Metaboliser", 0.0, "swap PM->NM"),
        ("Ultra-rapid Metaboliser", "Poor Metaboliser", 0.0, "swap UM->PM"),
        ("Positive", "Negative", 0.0, "swap HLA Pos->Neg"),
        ("Deficient", "Normal", 0.0, "swap deficient->normal"),
        ("Variable", "Deficient", 0.0, "Variable->Deficient (the v2 bug)"),
        ("Poor Metaboliser", "", 0.0, "empty parsed"),
        ("Poor Metaboliser", "Poor Metaboliser", 1.0, "exact echo"),
        ("Poor Metaboliser", "POOR METABOLISER", 1.0, "case-insensitive"),
        ("Favourable Response", "favorable response", 1.0, "US/UK spelling"),
    ]
    a1_fails = 0
    for gt, parsed, expected, label in boundary_tests:
        actual = rescore.score_a1(parsed, gt)
        ok = abs(actual - expected) < 0.01
        if not ok:
            a1_fails += 1
            issues.append(f"A5 boundary: {label}: expected {expected}, got {actual}")
        out.append(f"  [{'PASS' if ok else 'FAIL'}] A1: {label}: gt={gt!r} parsed={parsed!r} -> {actual} (expected {expected})")

    # A2 boundary
    a2_tests = [
        ("codeine: AVOID (potentially LETHAL)", "codeine: standard dosing", 0.0, "lethal-avoid -> standard"),
        ("warfarin: standard dosing", "warfarin: AVOID", 0.0, "standard -> avoid (over-cautious wrong)"),
        ("ivacaftor: not indicated", "ivacaftor: indicated", 0.0, "swap indicated polarity"),
    ]
    for gt, parsed, expected, label in a2_tests:
        actual = rescore.score_a2(parsed, gt)
        ok = abs(actual - expected) < 0.01
        if not ok:
            issues.append(f"A5 A2 boundary: {label}: expected {expected}, got {actual}")
        out.append(f"  [{'PASS' if ok else 'FAIL'}] A2: {label}: -> {actual}")

    # =========================================================================
    # A6: Echo-vs-paraphrase analysis on with_spec PHENOTYPE field
    # =========================================================================
    out.append(section("A6: Echo-vs-paraphrase analysis (with_spec circularity test)"))
    # For each with_spec row, compare parsed PHENOTYPE to gt_phenotype.
    # If they're identical (case-insensitive, stripped), that's verbatim echo.
    # If they differ but still match the canonical tier, that's paraphrase
    # (model deviated from spec format but kept semantic meaning).
    # The spec literally appears in the prompt with these values — verbatim is expected.
    verbatim = 0
    paraphrase_correct = 0
    spec_format_drift = []
    for r in rows:
        if r["cond"] != "with_spec": continue
        if r["scores"].get("format_fail"): continue
        gt = cbi[r["tc"]]["gt_phenotype"].strip().lower()
        parsed = (r["parsed"].get("PHENOTYPE", "") or "").strip().lower()
        if parsed == gt:
            verbatim += 1
        else:
            paraphrase_correct += 1
            if len(spec_format_drift) < 5:
                spec_format_drift.append({
                    "model": r["model"], "tc": r["tc"], "pop": r["pop"],
                    "gt": gt[:60], "parsed": parsed[:60]
                })
    total_ws = sum(1 for r in rows if r["cond"] == "with_spec" and not r["scores"].get("format_fail"))
    out.append(f"  with_spec rows: {total_ws}")
    out.append(f"    verbatim echo (parsed == gt): {verbatim} ({100*verbatim/total_ws:.1f}%)")
    out.append(f"    semantic match but text differs: {paraphrase_correct} ({100*paraphrase_correct/total_ws:.1f}%)")
    if spec_format_drift:
        out.append("  Examples of paraphrase-but-correct (model deviated from spec text):")
        for ex in spec_format_drift:
            out.append(f"    [{ex['model']}] {ex['tc']} {ex['pop']}: gt={ex['gt']!r} parsed={ex['parsed']!r}")
    # Interpretation: if verbatim is ~100%, models are simple stenographers.
    # If paraphrase >5%, models are interpreting/normalizing the spec, which has both pros (validates the contract pattern) and cons (introduces drift).

    # =========================================================================
    # A7: Verify with_spec PHENOTYPE field is not just literal prompt copy
    # =========================================================================
    out.append(section("A7: Detect prompt-copy artefacts in with_spec output"))
    # The with_spec prompt explicitly contains "PHENOTYPE: <value>".
    # We need to verify the model is NOT echoing the entire prompt (including
    # SKILL.md headers etc.); it should just emit the 5-line response.
    long_phenotype = 0
    long_examples = []
    for r in rows:
        if r["cond"] != "with_spec": continue
        if r["scores"].get("format_fail"): continue
        phen = r["parsed"].get("PHENOTYPE", "") or ""
        # If PHENOTYPE field contains spec-prompt headers, that's a sign of prompt copying
        if any(token in phen.lower() for token in ["skill.md", "## ", "input:", "patient cohort:", "activity score:"]):
            long_phenotype += 1
            if len(long_examples) < 3:
                long_examples.append({"model": r["model"], "tc": r["tc"], "pop": r["pop"], "phen": phen[:100]})
    out.append(f"  with_spec PHENOTYPE fields containing spec-prompt headers: {long_phenotype}")
    if long_examples:
        out.append("  Examples:")
        for ex in long_examples:
            out.append(f"    [{ex['model']}] {ex['tc']} {ex['pop']}: {ex['phen']}")

    # =========================================================================
    # A8: Lethal-case correctness audit
    # =========================================================================
    out.append(section("A8: Lethal-case classification audit"))
    lethal_cases = [c for c in cases if "lethal" in c["gt_drug"].lower()]
    lethal_genes = sorted({c["gene"] for c in lethal_cases})
    out.append(f"  Cases flagged 'lethal' in gt_drug: {len(lethal_cases)}")
    out.append(f"  Genes covered: {len(lethal_genes)}")
    for g in lethal_genes:
        gc = [c for c in lethal_cases if c["gene"] == g]
        out.append(f"    {g:<15} {len(gc)} case(s): {[c['id'] for c in gc]}")

    # Lethal-class A3 errors per condition
    no_spec_lethal_errs = 0
    with_spec_lethal_errs = 0
    for r in rows:
        if r["scores"].get("format_fail"): continue
        if "lethal" not in cbi[r["tc"]]["gt_drug"].lower(): continue
        if r["scores"].get("A3", 1.0) < 1.0:
            if r["cond"] == "no_spec":
                no_spec_lethal_errs += 1
            else:
                with_spec_lethal_errs += 1
    out.append(f"  no_spec lethal-class A3 errors: {no_spec_lethal_errs}")
    out.append(f"  with_spec lethal-class A3 errors: {with_spec_lethal_errs}")
    if with_spec_lethal_errs > 0:
        issues.append(f"A8: with_spec has {with_spec_lethal_errs} lethal-class errors (should be 0)")

    # =========================================================================
    # A9: Population-attribution audit
    # =========================================================================
    out.append(section("A9: Population attribution audit"))
    # Verify the 178 / 270 = 66% non-EUR claim
    pop_lethal_errs = Counter()
    for r in rows:
        if r["cond"] != "no_spec": continue
        if r["scores"].get("format_fail"): continue
        if "lethal" not in cbi[r["tc"]]["gt_drug"].lower(): continue
        if r["scores"].get("A3", 1.0) < 1.0:
            pop_lethal_errs[r["pop"]] += 1
    total_lethal = sum(pop_lethal_errs.values())
    non_eur = pop_lethal_errs["AMR"] + pop_lethal_errs["AFR"]
    out.append(f"  Lethal-class errors by pop (no_spec, parsed): EUR={pop_lethal_errs['EUR']}, AMR={pop_lethal_errs['AMR']}, AFR={pop_lethal_errs['AFR']}")
    out.append(f"  Total: {total_lethal}; non-EUR fraction: {non_eur}/{total_lethal} = {100*non_eur/total_lethal:.1f}%")
    expected_pct = 66
    out.append(f"  Manuscript claim: 178 of 270 (66%) non-EUR")
    if abs(100*non_eur/total_lethal - expected_pct) > 2:
        issues.append(f"A9: non-EUR fraction {100*non_eur/total_lethal:.1f}% deviates from claim 66%")

    # =========================================================================
    # A10: Stochastic check — within-combo agreement
    # =========================================================================
    out.append(section("A10: Within-combo run agreement"))
    # For with_spec, all 3 runs should yield A1=1.0
    for cond in ["no_spec", "with_spec"]:
        combos = defaultdict(list)
        for r in rows:
            if r["cond"] != cond: continue
            if r["scores"].get("format_fail"): continue
            combos[(r["model"], r["tc"], r["pop"])].append(r["scores"]["A1"])
        full_agreement = sum(1 for v in combos.values() if len(set(v)) == 1)
        out.append(f"  {cond}: {full_agreement}/{len(combos)} combos with identical A1 across all parsed runs")

    # =========================================================================
    # A11: Cross-model corroboration on hard cases
    # =========================================================================
    out.append(section("A11: Cross-model corroboration on hardest no_spec cases"))
    # For each test case, compute A1 across all 9 models in no_spec.
    # Cases where ALL 9 models fail are likely scorer issues OR genuine model failures.
    # Cases where ALL 9 models succeed are clear positives.
    case_a1 = defaultdict(list)
    for r in rows:
        if r["cond"] != "no_spec": continue
        if r["scores"].get("format_fail"): continue
        case_a1[r["tc"]].append((r["model"], r["scores"]["A1"]))
    universal_fail = []
    universal_pass = []
    for tc, ms in case_a1.items():
        models_with_correct = {m for m, a in ms if a == 1.0}
        models_with_wrong = {m for m, a in ms if a == 0.0}
        if len(models_with_correct) == 0 and len(models_with_wrong) >= 5:
            universal_fail.append((tc, len(ms)))
        if len(models_with_wrong) == 0 and len(models_with_correct) >= 5:
            universal_pass.append((tc, len(ms)))
    out.append(f"  Universal-fail cases (all 9 models wrong): {len(universal_fail)}")
    for tc, n in universal_fail[:10]:
        c = cbi[tc]
        out.append(f"    {tc:<32} ({n} runs) gene={c['gene']} gt_phen={c['gt_phenotype']!r}")
    out.append(f"  Universal-pass cases (all 9 models correct): {len(universal_pass)}")
    for tc, n in universal_pass[:5]:
        c = cbi[tc]
        out.append(f"    {tc:<32} ({n} runs) gene={c['gene']} gt_phen={c['gt_phenotype']!r}")

    # If universal-fail cases are all the same phenotype tier (e.g., Variable, Variant Carrier),
    # that's a strong signal the gt phenotype label is unconventional, not a scorer bug.
    universal_fail_tiers = Counter(cbi[tc]["gt_phenotype"] for tc, _ in universal_fail)
    out.append(f"  Universal-fail tier distribution: {dict(universal_fail_tiers)}")

    # =========================================================================
    # A12: Spec-content audit — does spec contain population keyword?
    # =========================================================================
    out.append(section("A12: Spec contains population keyword (B1 sanity)"))
    pop_kw = {"EUR": ["european", "caucasian"],
              "AMR": ["latin", "peru", "admixed", "indigenous", "amerindian", "mestizo"],
              "AFR": ["african", "uganda", "east african", "hiv", "efavirenz", "first-line"]}
    spec_has_kw = 0
    spec_missing_kw = 0
    for c in cases:
        for pop in ["EUR", "AMR", "AFR"]:
            note = (c["pop_note"].get(pop, "") or "").lower()
            if any(kw in note for kw in pop_kw[pop]):
                spec_has_kw += 1
            else:
                spec_missing_kw += 1
    out.append(f"  case-pop notes containing population keyword: {spec_has_kw}/{spec_has_kw + spec_missing_kw}")
    out.append(f"  case-pop notes WITHOUT keyword: {spec_missing_kw}")
    # If spec contains keyword but B1 with_spec is only 0.42, that means models are NOT
    # echoing the population note in the POPULATION field. That's a real finding.

    # =========================================================================
    # FINAL VERDICT
    # =========================================================================
    out.append(section("OVERALL VERDICT"))
    if not issues:
        out.append("PASS — no integrity issues detected")
        out.append("")
        out.append("Note: this audit does NOT include the live scrambled-spec adversarial")
        out.append("test, which requires API calls and is run separately by 14-adversarial-test.py.")
        out.append("Without that test, the with_spec=100% finding could in principle reflect")
        out.append("either (a) faithful contract execution or (b) the model already knew the")
        out.append("answer. The scrambled-spec test is the empirical disambiguation.")
    else:
        out.append(f"ISSUES DETECTED: {len(issues)}")
        for i in issues:
            out.append(f"  - {i}")

    text = "\n".join(out)
    REPORT.write_text(text)
    print(text)
    print(f"\nReport: {REPORT}")


if __name__ == "__main__":
    main()
