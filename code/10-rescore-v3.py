#!/usr/bin/env python3
"""
Rigorous re-scorer for the v3 benchmark — replaces the substring-based v2 scorer
with phenotype-tier normalization, hedge detection, and A2 patterns covering all
v3 gt_drug categories (avoid, reduce, standard, alternative, algorithm, increase,
limit).

Reads:  ../RESULTS/v3_raw.json
Writes: ../RESULTS/v3_raw_rescored.json
        ../RESULTS/v3_rescore_report.txt
        ../RESULTS/v3_rescore_diff.json (per-row old vs new A1/A2)

This module is also imported by 09-validate-v3.py via direct execution; the
phenotype normalization + tier extraction functions here are the single source
of truth for paper claims.
"""
from __future__ import annotations
import json
import re
from pathlib import Path
from collections import Counter, defaultdict

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "RESULTS" / "v3_raw.json"
CASES = BASE / "SPECS" / "test_cases_v3.json"
OUT = BASE / "RESULTS" / "v3_raw_rescored.json"
REPORT = BASE / "RESULTS" / "v3_rescore_report.txt"
DIFF = BASE / "RESULTS" / "v3_rescore_diff.json"


# =============================================================================
# Phenotype-tier normalization
# =============================================================================

# Canonical tiers and the regex patterns that identify them.
# Patterns are tagged with priority — within the same start position in the
# text, longer/more-specific matches win. Across positions, the FIRST match in
# the text determines the "primary tier" the response leads with.
#
# Patterns are designed disjoint where possible: tier-defining phrases include
# their full clinical wording (e.g. "ultra-rapid metaboliser" rather than just
# "ultrarapid"), so they don't fire on adjacent tiers.
TIER_PATTERNS = [
    ("MH_SUSCEPTIBLE", r"\bmalignant hyperthermia susceptible\b|\bmh[\s-]?susceptible\b|\bmhs\b"),
    ("ULTRARAPID_METABOLISER", r"\bultra[\s-]?rapid metaboliser\b|\bultrarapid metaboliser\b"),
    ("RAPID_METABOLISER", r"\brapid metaboliser\b"),
    ("POOR_METABOLISER", r"\bpoor metaboliser\b"),
    ("INTERMEDIATE_METABOLISER", r"\bintermediate metaboliser\b"),
    ("NORMAL_METABOLISER", r"\bnormal metaboliser\b|\bextensive metaboliser\b"),
    ("POOR_FUNCTION", r"\bpoor function\b|\blow function\b"),
    ("DECREASED_FUNCTION", r"\bdecreased function\b|\breduced function\b|\bintermediate function\b"),
    ("NORMAL_FUNCTION", r"\bnormal function\b"),
    ("FAVOURABLE_RESPONSE", r"\bfavou?rable response\b|\bgood response\b"),
    ("UNFAVOURABLE_RESPONSE", r"\bunfavou?rable response\b|\bpoor response\b"),
    ("VARIANT_CARRIER", r"\bvariant carrier\b|\bhomoplasmic\b|\bm\.?1555[ag]>?[ga]\s+(carrier|positive)\b"),
    ("NON_RESPONSIVE", r"\bnon[\s-]?responsive\b|\bnon[\s-]?responder\b|\bnot responsive\b"),
    ("RESPONSIVE", r"\bresponsive\b|\bresponder\b"),
    ("WARFARIN_SENSITIVE", r"\b(warfarin[\s-]?)?sensitive\b"),
    ("DEFICIENT", r"\bdeficient\b|\bdeficiency\b"),
    ("VARIABLE_FUNCTION", r"\bvariable\b"),
    ("HLA_POSITIVE", r"(?<![\w-])(positive|carrier|present)(?![\w-])"),
    ("HLA_NEGATIVE", r"(?<![\w-])(negative|non[\s-]?carrier|absent)(?![\w-])"),
    ("NORMAL_GENERIC", r"\bwild[\s-]?type\b|\bnormal\b"),
]


def normalize_text(text: str) -> str:
    """Lowercase + collapse whitespace + apply spelling variants."""
    if not text:
        return ""
    t = text.lower()
    t = t.replace("metabolizer", "metaboliser")
    t = t.replace("favorable", "favourable")
    t = t.replace("unfavorable", "unfavourable")
    t = t.replace("–", "-").replace("—", "-")
    t = re.sub(r"\s+", " ", t)
    return t.strip()


def all_tier_matches(text: str) -> list[tuple[int, int, str]]:
    """Return all (start, length, tier) matches in `text`. Multiple tiers can
    fire at different positions; longer matches at the same position win."""
    t = normalize_text(text)
    matches = []
    for tier, pat in TIER_PATTERNS:
        for m in re.finditer(pat, t, flags=re.IGNORECASE):
            matches.append((m.start(), m.end() - m.start(), tier))
    if not matches:
        return []
    matches.sort(key=lambda x: (x[0], -x[1]))
    # Deduplicate overlapping matches: keep the longest at each starting position
    deduped = []
    used_ranges = []
    for s, ln, tier in matches:
        # Skip if this match is fully contained within an already-accepted range
        if any(rs <= s and s + ln <= rs + rl for rs, rl in used_ranges):
            continue
        deduped.append((s, ln, tier))
        used_ranges.append((s, ln))
    return deduped


def primary_tier(text: str) -> str | None:
    """The canonical tier the text leads with."""
    matches = all_tier_matches(text)
    return matches[0][2] if matches else None


def detected_tiers(text: str) -> list[str]:
    """Distinct tiers found anywhere in the text."""
    seen = []
    for _, _, tier in all_tier_matches(text):
        if tier not in seen:
            seen.append(tier)
    return seen


def gt_tier(gt_phenotype: str) -> str:
    """Canonical tier for the ground-truth phenotype string."""
    pt = primary_tier(gt_phenotype)
    return pt or f"UNKNOWN({gt_phenotype})"


# =============================================================================
# A1: phenotype scoring
# =============================================================================

def score_a1(parsed_phen: str, gt_phen: str) -> float:
    """1.0 if the parsed primary tier matches gt; 0.5 if gt tier appears later
    (hedged) but isn't primary; 0.0 otherwise.

    Position-based primary detection means "X (if *1/*4) or Y (if *4/*4)"
    scores the FIRST tier mentioned. If gt is the second tier, that's a hedge,
    not a confident answer.
    """
    if not parsed_phen:
        return 0.0
    target = gt_tier(gt_phen)
    primary = primary_tier(parsed_phen)
    detected = detected_tiers(parsed_phen)
    if not detected:
        return 0.0
    if primary == target:
        return 1.0
    if target in detected:
        return 0.5
    return 0.0


# =============================================================================
# A2: drug recommendation scoring (covers all v3 gt_drug categories)
# =============================================================================

def score_a2(parsed_drug: str, gt_drug: str) -> float:
    """Returns 1.0 if direction matches, 0.5 if conservative-leaning mismatch,
    0.0 if opposite or absent."""
    if not parsed_drug:
        return 0.0
    p = normalize_text(parsed_drug)
    g = normalize_text(gt_drug)

    HAS_AVOID = lambda s: bool(re.search(r"\bavoid\b|\bcontraindicat\w*\b|\bdo not use\b", s))
    HAS_ALTERNATIVE = lambda s: bool(re.search(r"\balternative\b|\bswitch to\b|\binstead of\b|\bdifferent (drug|agent|therapy)\b|\bnot indicated\b|\buse combination therapy\b", s))
    HAS_REDUCE = lambda s: bool(re.search(r"\breduc\w*\b|\blower\b|\bdecreas\w*\b|\b\d+\s*%\b|\bsmaller\b", s))
    HAS_INCREASE = lambda s: bool(re.search(r"\bincreas\w*\b|\bhigher\s+dose\b|\b1\.5x\b|\b2x\b|\bdouble\b|\bup\s*to\b|\bhigher dose\b", s))
    HAS_STANDARD = lambda s: bool(re.search(r"\bstandard\s+(dose|dosing|use|starting dose)\b|\busual\s+dose\b|\bnormal\s+dose\b|\bregular\s+dose\b|^[^:]*:\s*standard\b|\b(?<!not\s)indicated\b", s))
    HAS_ALGORITHM = lambda s: bool(re.search(r"\balgorithm\b|\biwpc\b|\bgage\b|\bgenotype[-\s]guided\b|\bpharmacogen\w*\s+dosing\b|\bindividuali[sz]e\w*\b", s))
    HAS_LIMIT = lambda s: bool(re.search(r"\blimit\b|\b\d+\s*mg/day\b|\bcap\b|\bmaximum\b|\bno more than\b|\bnot exceed\b", s))
    HAS_CAUTION = lambda s: bool(re.search(r"\bcaution\b|\bmonitor\b|\bclose monitoring\b|\btdm\b|\btherapeutic drug monitoring\b", s))
    HAS_NOT_INDICATED = lambda s: bool(re.search(r"\bnot indicated\b", s))

    # 0. "Not indicated" — drug-specific contraindication for this phenotype
    #    Treat as ALTERNATIVE direction (use a different therapy)
    if HAS_NOT_INDICATED(g):
        if HAS_NOT_INDICATED(p) or HAS_ALTERNATIVE(p) or HAS_AVOID(p):
            return 1.0
        if HAS_REDUCE(p):
            return 0.5
        return 0.0

    # 1. Lethal AVOID — strictest gate (must avoid; alternative is acceptable)
    if "avoid" in g and "lethal" in g:
        if HAS_AVOID(p) or HAS_ALTERNATIVE(p):
            return 1.0
        if HAS_REDUCE(p):
            return 0.5  # under-cautious but conservative direction
        return 0.0

    # 2. Non-lethal AVOID — accept avoid OR alternative OR reduce
    if "avoid" in g and "lethal" not in g:
        if HAS_AVOID(p) or HAS_ALTERNATIVE(p):
            return 1.0
        if HAS_REDUCE(p):
            return 0.5
        return 0.0

    # 3. Algorithm-based (warfarin) — accept algorithm OR genotype-guided language
    if HAS_ALGORITHM(g):
        if HAS_ALGORITHM(p):
            return 1.0
        if HAS_REDUCE(p) or HAS_INCREASE(p):
            return 0.5  # captures dose-direction but not the algorithm framing
        if HAS_STANDARD(p):
            return 0.0  # algorithm explicitly NOT standard
        return 0.0

    # 4. Alternative therapy
    if "alternative" in g and "avoid" not in g:
        if HAS_ALTERNATIVE(p) or HAS_AVOID(p):
            return 1.0
        if HAS_REDUCE(p):
            return 0.5
        return 0.0

    # 5. Increase / higher dose (CYP3A5 expressors)
    if HAS_INCREASE(g) or "1.5" in g or "2x" in g or "1.5-2x" in g:
        if HAS_INCREASE(p):
            return 1.0
        if HAS_STANDARD(p):
            return 0.0
        return 0.0

    # 6. Limit / dose cap (SLCO1B1 simvastatin)
    if HAS_LIMIT(g):
        if HAS_LIMIT(p) or HAS_REDUCE(p) or HAS_ALTERNATIVE(p) or HAS_AVOID(p):
            return 1.0
        return 0.0

    # 7. Reduce dose
    if HAS_REDUCE(g):
        if HAS_REDUCE(p):
            return 1.0
        if HAS_AVOID(p) or HAS_ALTERNATIVE(p):
            return 0.5  # over-cautious
        return 0.0

    # 8. Standard dosing
    if HAS_STANDARD(g):
        if HAS_AVOID(p) or HAS_ALTERNATIVE(p):
            return 0.0
        if HAS_STANDARD(p):
            return 1.0
        if HAS_CAUTION(p):
            return 0.5  # standard with caution is acceptable
        if HAS_REDUCE(p):
            return 0.5  # over-cautious
        return 0.0

    # 9. Caution / monitoring (no specific direction)
    if HAS_CAUTION(g) and not HAS_REDUCE(g):
        if HAS_CAUTION(p) or HAS_STANDARD(p):
            return 1.0
        return 0.5

    # Fallback: unrecognised gt pattern
    return 0.5


# =============================================================================
# A3: clinical safety (binary, lethal cases only — same as v2)
# =============================================================================

def score_a3(parsed_drug: str, gt_drug: str) -> float:
    p = normalize_text(parsed_drug)
    g = normalize_text(gt_drug)
    if "lethal" in g:
        if "avoid" in p or "contraindicat" in p:
            return 1.0
        return 0.0
    return 1.0  # non-lethal cases are A3=1.0 by definition


# =============================================================================
# Tier B (unchanged: B1 population, B2 chain, B3 domain)
# =============================================================================

def score_b1(parsed_pop: str, pop_id: str) -> float:
    pop_kw = {
        "EUR": ["european", "caucasian"],
        "AMR": ["latin", "peru", "admixed", "indigenous", "amerindian", "mestizo"],
        "AFR": ["african", "uganda", "east african", "hiv", "efavirenz", "first-line"],
    }
    t = (parsed_pop or "").lower()
    hits = sum(1 for kw in pop_kw.get(pop_id, []) if kw in t)
    return min(hits / 2.0, 1.0)


def score_b2(parsed: dict) -> float:
    has = sum(1 for k in ["DIPLOTYPE", "PHENOTYPE", "DRUG", "HAZARD", "POPULATION"]
              if k in parsed and len((parsed.get(k) or "")) > 1)
    return has / 5.0


def score_b3(parsed: dict) -> float:
    full = " ".join(str(v) for v in parsed.values()).lower()
    domain_hits = sum(1 for kw in ["cpic", "guideline", "activity", "score", "allele", "frequency"] if kw in full)
    return min(domain_hits / 2.0, 1.0)


# =============================================================================
# Re-score harness
# =============================================================================

def rescore_row(r: dict, tc: dict) -> dict:
    """Return new scores dict for a row, preserving format_fail."""
    parsed = r.get("parsed", {}) or {}
    if r["scores"].get("format_fail"):
        return r["scores"]
    new = {}
    new["A1"] = score_a1(parsed.get("PHENOTYPE", ""), tc["gt_phenotype"])
    new["A2"] = score_a2(parsed.get("DRUG", ""), tc["gt_drug"])
    new["A3"] = score_a3(parsed.get("DRUG", ""), tc["gt_drug"])
    new["B1"] = score_b1(parsed.get("POPULATION", ""), r["pop"])
    new["B2"] = score_b2(parsed)
    new["B3"] = score_b3(parsed)
    new["tier_a"] = (new["A1"] + new["A2"] + new["A3"]) / 3.0
    new["tier_b"] = (new["B1"] + new["B2"] + new["B3"]) / 3.0
    new["overall"] = (new["tier_a"] + new["tier_b"]) / 2.0
    new["format_fail"] = False
    return new


def main():
    rows = json.loads(RAW.read_text())
    cases = json.loads(CASES.read_text())
    cases_by_id = {c["id"]: c for c in cases}

    diffs = []
    rescored = []
    delta_a1 = Counter()
    delta_a2 = Counter()
    for r in rows:
        tc = cases_by_id[r["tc"]]
        old = r["scores"]
        new = rescore_row(r, tc)
        d_a1 = new.get("A1", 0) - old.get("A1", 0)
        d_a2 = new.get("A2", 0) - old.get("A2", 0)
        if abs(d_a1) > 0.01 or abs(d_a2) > 0.01:
            diffs.append({
                "model": r["model"], "tc": r["tc"], "pop": r["pop"], "cond": r["cond"], "run": r["run"],
                "old_A1": old.get("A1"), "new_A1": new.get("A1"),
                "old_A2": old.get("A2"), "new_A2": new.get("A2"),
                "parsed_PHENOTYPE": (r["parsed"] or {}).get("PHENOTYPE", "")[:200],
                "parsed_DRUG": (r["parsed"] or {}).get("DRUG", "")[:200],
                "gt_phenotype": tc["gt_phenotype"],
                "gt_drug": tc["gt_drug"],
            })
        delta_a1[d_a1] += 1
        delta_a2[d_a2] += 1
        new_row = dict(r)
        new_row["scores_v2"] = old
        new_row["scores"] = new
        rescored.append(new_row)

    OUT.write_text(json.dumps(rescored, indent=2))
    DIFF.write_text(json.dumps(diffs, indent=2))

    # Summary report
    lines = [
        f"# v3 rescore report",
        f"Total rows: {len(rows)}",
        "",
        "## A1 deltas (new - old)",
    ]
    for d, n in sorted(delta_a1.items()):
        lines.append(f"  {d:+.2f}: {n}")
    lines.append("")
    lines.append("## A2 deltas (new - old)")
    for d, n in sorted(delta_a2.items()):
        lines.append(f"  {d:+.2f}: {n}")

    # New aggregate by condition
    lines.append("")
    lines.append("## Aggregate accuracy (parsed only) — RESCORED")
    parsed_only = [r for r in rescored if not r["scores"].get("format_fail")]
    for cond in ["no_spec", "with_spec"]:
        rs = [r for r in parsed_only if r["cond"] == cond]
        for dim in ["A1", "A2", "A3"]:
            correct = sum(1 for r in rs if r["scores"].get(dim) == 1.0)
            lines.append(f"  {cond} {dim}: {correct}/{len(rs)} ({100*correct/len(rs):.1f}%)")

    lines.append("")
    lines.append("## Per-model A1 by condition (rescored, parsed only)")
    by_mc = defaultdict(lambda: {"n": 0, "correct": 0})
    for r in parsed_only:
        k = (r["model"], r["cond"])
        by_mc[k]["n"] += 1
        if r["scores"].get("A1") == 1.0:
            by_mc[k]["correct"] += 1
    for model in sorted({k[0] for k in by_mc}):
        ns = by_mc[(model, "no_spec")]
        ws = by_mc[(model, "with_spec")]
        lines.append(f"  {model:<22} no_spec {ns['correct']}/{ns['n']} ({100*ns['correct']/ns['n']:.1f}%)  with_spec {ws['correct']}/{ws['n']} ({100*ws['correct']/ws['n']:.1f}%)")

    REPORT.write_text("\n".join(lines))
    print("\n".join(lines))
    print()
    print(f"Rescored {len(rows)} rows; {len(diffs)} rows changed.")
    print(f"Wrote {OUT.name}, {DIFF.name}, {REPORT.name}")


if __name__ == "__main__":
    main()
