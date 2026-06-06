#!/usr/bin/env python3
"""
Clinical-equivalence rescorer extension for the v3 benchmark.

Wraps `10-rescore-v3.py` (which is locked per the v3 benchmark spec) with a
narrow `(gene, phrase) -> tier` table that recognises clinically-equivalent
natural-language phrasings as members of their canonical CPIC tier. Applied
uniformly to all conditions and runs.

Seed table:

  HLA-B*57:01-specific (abacavir hypersensitivity phrasing):
    "high risk of [drug?] hypersensitivity"      -> HLA_POSITIVE
    "increased risk of [drug?] hypersensitivity" -> HLA_POSITIVE
    "carrier" + locus mention                    -> HLA_POSITIVE
    "no increased risk"                          -> HLA_NEGATIVE
    "low risk of hypersensitivity"               -> HLA_NEGATIVE
    "non-carrier"                                -> HLA_NEGATIVE

  Shared HLA risk-allele patterns (apply across HLA-A*31:01, HLA-B*15:02,
  HLA-B*57:01, HLA-B*58:01 — the four CPIC carbamazepine/abacavir/allopurinol
  risk-allele loci):
    "(significantly|substantially)? increased risk of <reaction>" -> HLA_POSITIVE
    "(normal|low) (or|/)? reduced risk of <reaction>"             -> HLA_NEGATIVE
    "normal risk of <reaction>"                                   -> HLA_NEGATIVE
    "homozygous for alleles other than HLA-<locus>"               -> HLA_NEGATIVE
    "(heterozygote|carrier) for HLA-<locus>"                      -> HLA_POSITIVE

  Where <reaction> ∈ {sjs, ten, dress, mpe, scar, hypersensitivity,
                       cutaneous, adverse}.

The optional drug-name word inside "high/increased risk of [drug?] hypersensitivity"
is matched as ([A-Za-z]+\s+)? rather than (\w+\s+)? on purpose: \w admits digits
and underscores, which would let pathological tokens like "HLA_B5701" or numeric
strings sneak in. Clinical drug names are alphabetic; restrict accordingly.

Safety property — by construction:
  The wrapper only PROMOTES A1 from <1.0 to 1.0 when an equivalent phrase
  matches AND its tier == target_tier. It never demotes. 1->0 transitions
  cannot occur. The original score is the floor.

Reads:  ../RESULTS/v3_raw_rescored.json   (rigorous 10- baseline, locked)
Writes: ../RESULTS/v3_raw_rescored_clinical_eq.json
        ../RESULTS/v3_rescore_clinical_eq_report.txt

The baseline for comparison is the RIGOROUS 10- score in v3_raw_rescored.json,
NOT the legacy v2 substring score in v3_raw.json. The clinical-equivalence
wrapper extends 10-, so the only legitimate diff is the 0->1 promotions where
an HLA-B*57:01 equivalence phrase matched.

Run with --self-test for sanity checks (TDD-lite) before scoring real data.
"""
from __future__ import annotations
import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
PYDIR = BASE / "PYTHON"
BASELINE = BASE / "RESULTS" / "v3_raw_rescored.json"   # rigorous 10- output
CASES = BASE / "SPECS" / "test_cases_v3.json"
OUT = BASE / "RESULTS" / "v3_raw_rescored_clinical_eq.json"
REPORT = BASE / "RESULTS" / "v3_rescore_clinical_eq_report.txt"

# Import original rescorer (file name starts with digit, must use spec loader)
_spec = spec_from_file_location("rescore_v3", PYDIR / "10-rescore-v3.py")
rescore = module_from_spec(_spec)
_spec.loader.exec_module(rescore)


# =============================================================================
# Clinical equivalence table
# =============================================================================
# Each entry: (gene, list of (regex, canonical_tier)).
# Patterns operate on `normalize_text(parsed_PHENOTYPE)` (lowercase, collapsed
# whitespace, en/em-dashes normalised, US/UK spelling unified).
#
# Discipline: every entry must be justified by an audited dry-run example.
# Do not anticipate. Expand only after observed cases motivate the addition.

# Locus mention required for the bare "carrier" pattern: the response must
# reference the HLA-B*57:01 locus somewhere (so generic "carrier" wording in
# non-HLA contexts cannot accidentally promote).
# Per-locus regexes — for "carrier" disambiguation (so bare "carrier" inside the
# parsed PHENOTYPE only promotes when the response also mentions the right HLA).
LOCUS_PATTERNS = {
    "HLA-B*57:01": re.compile(
        r"hla[\s-]?b[\s\*:_-]*57[\s:_-]*01|hla[\s-]?b\*?5701|b\*?5701",
        re.IGNORECASE),
    "HLA-A*31:01": re.compile(
        r"hla[\s-]?a[\s\*:_-]*31[\s:_-]*01|hla[\s-]?a\*?3101|a\*?3101",
        re.IGNORECASE),
    "HLA-B*15:02": re.compile(
        r"hla[\s-]?b[\s\*:_-]*15[\s:_-]*02|hla[\s-]?b\*?1502|b\*?1502",
        re.IGNORECASE),
    "HLA-B*58:01": re.compile(
        r"hla[\s-]?b[\s\*:_-]*58[\s:_-]*01|hla[\s-]?b\*?5801|b\*?5801",
        re.IGNORECASE),
}

# Reaction keywords for the generic HLA-risk-allele patterns.
# Required so the "increased risk of <X>" template only fires on clinically
# relevant adverse reactions, not on unrelated "increased risk of ..." text.
_HLA_REACTION = (
    r"(?:sjs|ten|dress|mpe|scar|hypersensitivity|cutaneous|adverse)"
)
# Up to ~50 alphabetic chars and separators between "of" and reaction noun,
# admitting drug/condition words like "carbamazepine-induced", "allopurinol".
_HLA_INFIX = r"(?:[A-Za-z][\w\s\-/,]{0,60}?\s)?"

# Generic patterns for any HLA risk-allele locus. Applied IN ADDITION to each
# locus's specific entries. ORDER MATTERS — Negatives first.
GENERIC_HLA_RISK_ALLELE_PATTERNS: list[tuple[re.Pattern, str, str]] = [
    # Negative: "normal/low [or/]reduced risk of <reaction>"
    (re.compile(
        rf"\b(?:normal|low)\s*(?:or|/)?\s*reduced\s*risk\s+(?:of|for)\s+{_HLA_INFIX}{_HLA_REACTION}\b",
        re.IGNORECASE),
     "HLA_NEGATIVE", ""),
    # Negative: bare "normal risk of <reaction>"
    (re.compile(
        rf"\bnormal\s+risk\s+(?:of|for)\s+{_HLA_INFIX}{_HLA_REACTION}\b",
        re.IGNORECASE),
     "HLA_NEGATIVE", ""),
    # Negative: "homozygous for alleles other than HLA-..."
    (re.compile(
        r"\bhomozygous\s+for\s+alleles?\s+other\s+than\s+hla[\s-]",
        re.IGNORECASE),
     "HLA_NEGATIVE", ""),
    # Positive: "(significantly|substantially)? increased risk of <reaction>"
    (re.compile(
        rf"\b(?:significantly|substantially)?\s*increased\s+risk\s+(?:of|for)\s+{_HLA_INFIX}{_HLA_REACTION}\b",
        re.IGNORECASE),
     "HLA_POSITIVE", ""),
    # Positive: "heterozygote/heterozygous for HLA-<locus>" — locus check ensures
    # the right allele is named (caller passes the gt gene).
    (re.compile(
        r"\b(?:heterozygote|heterozygous)\b[^.]{0,40}\bfor\s+hla[\s-]",
        re.IGNORECASE),
     "HLA_POSITIVE", "LOCUS"),
]

CLINICAL_EQUIVALENCE: dict[str, list[tuple[re.Pattern, str, str]]] = {
    # gene -> [(compiled regex, canonical tier, locus_filter?)]
    # locus_filter: "" = no extra check; "LOCUS" = require the gt-gene's locus
    # regex to also match somewhere in the parsed text.
    "HLA-B*57:01": [
        # Negative patterns first (first-match-wins resolves "no increased risk"
        # before the broader Positive template would claim it).
        (re.compile(r"\bno increased risk\b", re.IGNORECASE),
         "HLA_NEGATIVE", ""),
        (re.compile(r"\blow risk of hypersensitivity\b", re.IGNORECASE),
         "HLA_NEGATIVE", ""),
        (re.compile(r"\bnon[\s-]?carrier\b", re.IGNORECASE),
         "HLA_NEGATIVE", ""),
        (re.compile(r"\bhigh risk of\s+(?:[A-Za-z]+\s+)?hypersensitivity\b", re.IGNORECASE),
         "HLA_POSITIVE", ""),
        (re.compile(r"\bincreased risk of\s+(?:[A-Za-z]+\s+)?hypersensitivity\b", re.IGNORECASE),
         "HLA_POSITIVE", ""),
        (re.compile(r"\bcarrier\b", re.IGNORECASE),
         "HLA_POSITIVE", "LOCUS"),
    ],
    # The three new HLA risk-allele loci share the same generic patterns
    # (carbamazepine/oxcarbazepine SJS/TEN/DRESS/MPE for HLA-A*31:01 +
    # HLA-B*15:02; allopurinol SCAR for HLA-B*58:01). Per-locus entries are
    # empty — the generic patterns plus locus filter handle all observed cases.
    "HLA-A*31:01": [],
    "HLA-B*15:02": [],
    "HLA-B*58:01": [],
}

# Genes that should additionally check the generic HLA-risk-allele patterns
HLA_RISK_ALLELE_GENES = frozenset({
    "HLA-A*31:01", "HLA-B*15:02", "HLA-B*57:01", "HLA-B*58:01"
})


def _applicable_patterns(gene: str) -> list[tuple[re.Pattern, str, str]]:
    """Return the ordered list of (pattern, tier, locus_filter) tuples to try
    for `gene`. Per-locus entries first, then generic HLA-risk-allele patterns
    if the gene is in HLA_RISK_ALLELE_GENES."""
    specific = CLINICAL_EQUIVALENCE.get(gene, [])
    if gene in HLA_RISK_ALLELE_GENES:
        return list(specific) + list(GENERIC_HLA_RISK_ALLELE_PATTERNS)
    return list(specific)


def _locus_filter_passes(filter_tag: str, gene: str, text: str) -> bool:
    """Check the locus-filter sentinel.
    - "" : no extra check.
    - "LOCUS" : require the gt-gene's locus regex to also match in the text.
    - explicit gene name : require that named locus's regex to match.
    """
    if not filter_tag:
        return True
    locus_gene = gene if filter_tag == "LOCUS" else filter_tag
    locus_re = LOCUS_PATTERNS.get(locus_gene)
    if locus_re is None:
        return False
    return bool(locus_re.search(text))


def equivalence_tier(gene: str, parsed_phen: str) -> str | None:
    """Return canonical tier if a clinical-equivalent phrase matches the gene's
    table, else None. Returns the FIRST matching tier; ambiguity (positive +
    negative phrases in the same response) is resolved as no-promotion at the
    score_a1 layer (only promote if matched_tier == target_tier)."""
    if not parsed_phen:
        return None
    text = rescore.normalize_text(parsed_phen)
    for pattern, tier, locus_filter in _applicable_patterns(gene):
        if not pattern.search(text):
            continue
        if not _locus_filter_passes(locus_filter, gene, text):
            continue
        return tier
    return None


def equivalence_tiers_all(gene: str, parsed_phen: str) -> list[str]:
    """Return all distinct equivalence tiers found (used by self-test to detect
    ambiguous responses)."""
    if not parsed_phen:
        return []
    text = rescore.normalize_text(parsed_phen)
    seen = []
    for pattern, tier, locus_filter in _applicable_patterns(gene):
        if not pattern.search(text):
            continue
        if not _locus_filter_passes(locus_filter, gene, text):
            continue
        if tier not in seen:
            seen.append(tier)
    return seen


def score_a1_clinical_eq(parsed_phen: str, gt_phen: str, gene: str) -> float:
    """Extended A1 scorer.

    1. Compute baseline = original rescore.score_a1(parsed, gt).
    2. If baseline == 1.0, return 1.0 (no override needed).
    3. Else, look up (gene, parsed_phen) in the equivalence table.
       If a matching equivalence tier == target tier, return 1.0.
    4. Else, return baseline.
    """
    baseline = rescore.score_a1(parsed_phen, gt_phen)
    if baseline >= 1.0:
        return baseline
    target = rescore.gt_tier(gt_phen)
    eq_tier = equivalence_tier(gene, parsed_phen)
    if eq_tier is not None and eq_tier == target:
        return 1.0
    return baseline


def rescore_row_clinical_eq(r: dict, tc: dict) -> dict:
    """Rescore one row using the clinical-equivalence A1 wrapper, leaving A2,
    A3, B1, B2, B3 unchanged (they don't depend on phenotype-tier matching)."""
    parsed = r.get("parsed", {}) or {}
    if r["scores"].get("format_fail"):
        return r["scores"]
    new = {}
    new["A1"] = score_a1_clinical_eq(
        parsed.get("PHENOTYPE", ""),
        tc["gt_phenotype"],
        tc["gene"],
    )
    new["A2"] = rescore.score_a2(parsed.get("DRUG", ""), tc["gt_drug"])
    new["A3"] = rescore.score_a3(parsed.get("DRUG", ""), tc["gt_drug"])
    new["B1"] = rescore.score_b1(parsed.get("POPULATION", ""), r["pop"])
    new["B2"] = rescore.score_b2(parsed)
    new["B3"] = rescore.score_b3(parsed)
    new["tier_a"] = (new["A1"] + new["A2"] + new["A3"]) / 3.0
    new["tier_b"] = (new["B1"] + new["B2"] + new["B3"]) / 3.0
    new["overall"] = (new["tier_a"] + new["tier_b"]) / 2.0
    new["format_fail"] = False
    return new


# =============================================================================
# Self-test (TDD red/green)
# =============================================================================

def run_self_test() -> int:
    """Returns 0 if all tests pass, 1 otherwise."""
    cases = [
        # (gene, gt_phen, parsed_phen, expected_a1, label)
        # 1. Bare canonical labels — parity with original
        ("HLA-B*57:01", "Positive", "Positive", 1.0, "canonical Positive parity"),
        ("HLA-B*57:01", "Negative", "Negative", 1.0, "canonical Negative parity"),
        # 2. Clinical-equivalent phrases promote correctly
        ("HLA-B*57:01", "Positive",
         "High risk of hypersensitivity reaction", 1.0,
         "high risk -> Positive (gt Positive)"),
        ("HLA-B*57:01", "Positive",
         "Increased risk of hypersensitivity", 1.0,
         "increased risk -> Positive (gt Positive)"),
        # 2b. Drug-word inserted between "of" and "hypersensitivity"
        ("HLA-B*57:01", "Positive",
         "High risk of abacavir hypersensitivity", 1.0,
         "high risk of abacavir hypersensitivity -> Positive"),
        ("HLA-B*57:01", "Positive",
         "Increased risk of abacavir hypersensitivity reaction", 1.0,
         "increased risk of abacavir hypersensitivity reaction -> Positive"),
        ("HLA-B*57:01", "Positive",
         "Increased risk of Abacavir Hypersensitivity Reaction", 1.0,
         "TitleCase variant -> Positive"),
        # 2c. Drug-word slot does NOT permit digits/underscores
        # (digits and underscores are not [A-Za-z], so the regex fails to
        # admit pathological tokens; original anchor still requires "of <word>
        # hypersensitivity")
        ("HLA-B*57:01", "Positive",
         "High risk of HLA_B5701 hypersensitivity", 0.0,
         "underscore-bearing token NOT matched (regex tightening rationale)"),
        ("HLA-B*57:01", "Negative",
         "No increased risk of abacavir hypersensitivity", 1.0,
         "no increased risk -> Negative (gt Negative)"),
        ("HLA-B*57:01", "Negative",
         "Low risk of hypersensitivity", 1.0,
         "low risk -> Negative (gt Negative)"),
        ("HLA-B*57:01", "Negative",
         "Non-carrier", 1.0,
         "non-carrier -> Negative (gt Negative)"),
        # 3. Equivalence does NOT promote when target mismatches
        ("HLA-B*57:01", "Positive",
         "Low risk of hypersensitivity", 0.0,
         "low risk does NOT promote when gt is Positive"),
        ("HLA-B*57:01", "Negative",
         "High risk of hypersensitivity reaction", 0.0,
         "high risk does NOT promote when gt is Negative"),
        # 4. Gene context required: same phrase, different gene -> no promotion
        ("CYP2D6", "Poor Metaboliser",
         "High risk of hypersensitivity", 0.0,
         "non-HLA gene: equivalence does not fire"),
        # 5. Locus mention required for "carrier" alone
        ("HLA-B*57:01", "Positive",
         "Carrier of HLA-B*57:01", 1.0,
         "carrier with locus -> Positive"),
        # NB: bare "carrier" without locus is already caught by the original
        # HLA_POSITIVE pattern (\bcarrier\b), so it scores 1.0 via baseline,
        # not via equivalence. We confirm parity:
        ("HLA-B*57:01", "Positive", "Carrier", 1.0,
         "bare carrier handled by original (baseline)"),
        # 6. Format-fail-style empty phenotype
        ("HLA-B*57:01", "Positive", "", 0.0,
         "empty parsed -> 0 (no override)"),
        # 7. HLA-A*31:01 (carbamazepine; SJS/TEN/DRESS/MPE)
        ("HLA-A*31:01", "Positive",
         "Increased risk of carbamazepine-induced SJS/TEN, DRESS and MPE",
         1.0, "A*31:01 Pos: increased risk of carbamazepine SJS"),
        ("HLA-A*31:01", "Negative",
         "Normal or reduced risk of carbamazepine-induced SJS/TEN, DRESS and MPE",
         1.0, "A*31:01 Neg: normal or reduced risk"),
        ("HLA-A*31:01", "Negative",
         "Homozygous for alleles other than HLA-B*15:02 and HLA-A*31:01",
         1.0, "A*31:01 Neg: homozygous for alleles other than"),
        # 8. HLA-B*15:02 (carbamazepine/oxcarbazepine; SJS/TEN)
        ("HLA-B*15:02", "Positive",
         "Increased risk of carbamazepine/oxcarbazepine-induced SJS/TEN",
         1.0, "B*15:02 Pos: increased risk of carbamazepine/oxcarbazepine"),
        ("HLA-B*15:02", "Negative",
         "Normal or reduced risk of carbamazepine/oxcarbazepine-induced SJS/TEN",
         1.0, "B*15:02 Neg: normal or reduced risk"),
        # 9. HLA-B*58:01 (allopurinol; SCAR)
        ("HLA-B*58:01", "Positive",
         "Significantly increased risk of allopurinol SCAR",
         1.0, "B*58:01 Pos: significantly increased risk of allopurinol SCAR"),
        ("HLA-B*58:01", "Positive",
         "Significantly increased risk of allopurinol-induced SCAR",
         1.0, "B*58:01 Pos: hyphenated drug variant"),
        ("HLA-B*58:01", "Negative",
         "Low or reduced risk of allopurinol SCAR",
         1.0, "B*58:01 Neg: low or reduced risk"),
        # 10. Cross-gene safety: generic patterns DON'T fire on non-HLA-risk genes
        ("CYP2D6", "Poor Metaboliser",
         "Increased risk of hypersensitivity",
         0.0, "non-risk-allele gene: generic HLA patterns do NOT fire"),
        # 11. Locus filter: "heterozygote for HLA-..." promotes only when locus
        # matches the gt gene
        ("HLA-A*31:01", "Positive",
         "Heterozygote or homozygous for HLA-A*31:01",
         1.0, "A*31:01 Pos: locus matches gt gene"),
    ]
    fails = []
    for gene, gt, parsed, expected, label in cases:
        actual = score_a1_clinical_eq(parsed, gt, gene)
        ok = abs(actual - expected) < 0.01
        marker = "PASS" if ok else "FAIL"
        print(f"  [{marker}] {label}")
        print(f"         gene={gene!r} gt={gt!r} parsed={parsed!r}")
        print(f"         expected={expected}  actual={actual}")
        if not ok:
            fails.append(label)
    print()
    if fails:
        print(f"SELF-TEST FAILED: {len(fails)} cases")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("SELF-TEST PASSED: all cases green")
    return 0


# =============================================================================
# Main rescore loop
# =============================================================================

def main_rescore() -> int:
    # Baseline is the rigorous 10- output (v3_raw_rescored.json). Each row has
    # `scores` = rigorous 10- A1/A2/A3, plus the parsed response. We extend the
    # A1 only via the clinical-equivalence wrapper.
    rows = json.loads(BASELINE.read_text())
    cases = json.loads(CASES.read_text())
    cases_by_id = {c["id"]: c for c in cases}

    rescored = []
    n_promoted = 0           # rows where new A1 > rigorous baseline A1
    n_demoted = 0            # rows where new A1 < rigorous baseline A1 (must be 0)
    promotion_examples = []  # for the report

    for r in rows:
        tc = cases_by_id[r["tc"]]
        baseline_scores = r["scores"]                    # rigorous 10-
        new = rescore_row_clinical_eq(r, tc)
        old_a1 = baseline_scores.get("A1", 0)
        new_a1 = new.get("A1", 0)
        if new_a1 > old_a1 + 0.01:
            n_promoted += 1
            if len(promotion_examples) < 25:
                promotion_examples.append({
                    "model": r["model"], "tc": r["tc"], "pop": r["pop"],
                    "cond": r["cond"], "run": r["run"],
                    "gene": tc["gene"],
                    "gt_phenotype": tc["gt_phenotype"],
                    "parsed_PHENOTYPE": (r["parsed"] or {}).get("PHENOTYPE", "")[:200],
                    "old_A1": old_a1, "new_A1": new_a1,
                })
        if new_a1 < old_a1 - 0.01:
            n_demoted += 1

        new_row = dict(r)
        new_row["scores_rigorous"] = baseline_scores  # preserve 10- baseline
        new_row["scores"] = new                         # clinical-eq scores
        rescored.append(new_row)

    OUT.write_text(json.dumps(rescored, indent=2))

    # Summary report
    lines = [
        f"# v3 clinical-equivalence rescore report",
        f"Baseline:   {BASELINE.name}  (rigorous 10- output)",
        f"Output:     {OUT.name}",
        f"Total rows: {len(rows)}",
        "",
        f"Promoted (A1 baseline<1 -> 1.0):  {n_promoted}",
        f"Demoted  (A1 baseline=1 -> <1):   {n_demoted}  (must be 0 by construction)",
        "",
    ]
    if promotion_examples:
        lines.append("## Sample promotions (up to 25)")
        for ex in promotion_examples:
            lines.append(
                f"  {ex['model']:<22} {ex['tc']:<28} {ex['pop']:<3} "
                f"{ex['cond']:<10} run={ex['run']} gene={ex['gene']}"
            )
            lines.append(
                f"    gt_phenotype: {ex['gt_phenotype']}"
            )
            lines.append(
                f"    parsed:       {ex['parsed_PHENOTYPE']}"
            )
            lines.append(
                f"    A1 {ex['old_A1']} -> {ex['new_A1']}"
            )
            lines.append("")
    REPORT.write_text("\n".join(lines))
    print("\n".join(lines))
    print()
    print(f"Wrote {OUT.name}, {REPORT.name}")
    if n_demoted > 0:
        print(f"WARNING: {n_demoted} rows demoted; this should be impossible by construction.")
        return 1
    return 0


# =============================================================================
# Entry point
# =============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true",
                        help="Run TDD self-test and exit")
    args = parser.parse_args()
    if args.self_test:
        return run_self_test()
    return main_rescore()


if __name__ == "__main__":
    sys.exit(main())
