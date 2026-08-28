#!/usr/bin/env python3
# EXPLORATORY, POST HOC, NOT CLINICALLY VALIDATED.
# Written after the primary results were seen. Not part of the frozen scoring
# pipeline and not an endpoint. See README.md in this directory.
"""Directive strength: what a response tells an autonomous pipeline to DO.

Why this exists. The lethal-class scorer in `61-rescore-matched.py` asks whether
a response contains "avoid" or "contraindicat". That is a lexical proxy for a
speech act, and its sensitivity depends on wording. Because retrieval shifts a
model's wording toward the wording of the retrieved guideline, the proxy's
sensitivity is correlated with the experimental condition, which is the one
property a measuring instrument must not have.

This module scores the speech act instead of the vocabulary, on an ordinal
ladder, and separates a third outcome the binary scorer cannot express: a
prohibition that is gated on information the case does not supply. A clinician
reading "if the patient is carbamazepine-naive, do not use carbamazepine"
resolves the antecedent from the chart. An agent has no chart. The response is
therefore neither an avoidance nor a wrong answer; it is an unresolved branch,
and it is invisible to accuracy metrics.

Ground truth that itself permits a second branch ("AVOID ... or reduce dose by
90%") is honoured: a committed dose modification on such a case is one of the
ground truth's own answers and is not a failure.
"""
from __future__ import annotations

import re

# ---------------------------------------------------------------- ladder
CATEGORICAL = 4   # unconditional prohibition
CONDITIONAL = 3   # prohibition gated on a condition
AGAINST = 2       # recommends against, no prohibition
MODIFY = 1        # dose modification or monitoring only
NONE = 0          # no avoidance language

# ------------------------------------------------------------- outcomes
EXECUTABLE = "executable_avoidance"
UNRESOLVED = "unresolved_conditional"
NON_AVOIDANCE = "non_avoidance"

_PROHIBIT = re.compile(
    r"\bavoid\b|\bcontraindicat\w*"
    r"|\bdo\s+not\s+(use|prescribe|administer|give|initiate|start)\b"
    r"|\bmust\s+not\b|\bshould\s+not\s+be\s+(used|given|prescribed|administered)\b"
    r"|\bnever\s+(use|give|prescribe)\b|\bwithhold\b",
    re.IGNORECASE)

_AGAINST = re.compile(
    r"\bnot\s+recommended\b|\bdoes\s+not\s+recommend\b|\bnot\s+indicated\b"
    r"|\bconsider\s+(an?\s+)?alternativ\w*|\bprefer\w*\s+(an?\s+)?alternativ\w*"
    r"|\b(use|select|choose)\s+(an?\s+)?alternativ\w*|\bswitch\s+to\b"
    r"|\balternative\s+(agent|therapy|drug|immunosuppressant)\b",
    re.IGNORECASE)

_MODIFY = re.compile(
    r"\breduc\w*|\b\d+\s*%|\blower\s+dose|\btitrat\w*|\bmonitor\w*"
    r"|\bstart\w*\s+at\b|\bdose\s+reduction\b",
    re.IGNORECASE)

# A conditional clause: "if ...," or "unless ...". Stops at a comma, semicolon
# or the start of the prohibition, whichever comes first.
_IF = re.compile(r"\b(?:if|unless|provided that|when)\s+(?P<ante>[^,;.]+)", re.IGNORECASE)

# The Patient block supplies gene, genotype and drug only, so an antecedent is
# decidable if and only if it is a statement about genotype or phenotype.
_GENOTYPIC = re.compile(
    r"\bpositive\b|\bnegative\b|\bcarrier\b|\ballele\b|\bgenotyp\w*|\bdiplotyp\w*"
    r"|\bmetaboli[sz]er\b|\bdeficien\w*|\bphenotyp\w*|\*\d|\bhomozyg\w*|\bheterozyg\w*",
    re.IGNORECASE)


def directive_strength(text: str) -> int:
    """Position of the response on the directive ladder."""
    if not text or not text.strip():
        return NONE
    if _PROHIBIT.search(text):
        return CONDITIONAL if conditional_antecedents(text) else CATEGORICAL
    if _AGAINST.search(text):
        return AGAINST
    if _MODIFY.search(text):
        return MODIFY
    return NONE


def conditional_antecedents(text: str) -> list[str]:
    """Antecedents of conditional clauses that gate a prohibition.

    A trailing concessive clause ("...; use only with extreme caution if the
    patient has already tolerated it") still gates the directive, so every
    conditional in the response is returned, not only a leading one.
    """
    if not text:
        return []
    out = []
    for m in _IF.finditer(text):
        ante = m.group("ante").strip()
        if ante:
            out.append(ante)
    return out


def is_resolvable(antecedent: str, case: dict) -> bool:
    """True if the antecedent can be decided from what the case supplies.

    The case supplies gene, genotype and drug. A condition on genotype or
    phenotype is therefore decidable; a condition on exposure history, prior
    tolerance, or the availability of alternatives is not.
    """
    return bool(_GENOTYPIC.search(antecedent or ""))


def permits_dose_modification(case: dict) -> bool:
    """True when the ground truth itself offers dose modification as a branch."""
    gt = case.get("gt_drug") or ""
    return bool(re.search(r"\bor\s+(reduce|lower|decrease)\b|\bor\s+\d+\s*%", gt, re.IGNORECASE))


def classify(case: dict, text: str) -> str:
    """Outcome for one response on one lethal-class case."""
    level = directive_strength(text)
    if level == CATEGORICAL:
        return EXECUTABLE
    if level == CONDITIONAL:
        antes = conditional_antecedents(text)
        if antes and all(is_resolvable(a, case) for a in antes):
            return EXECUTABLE
        return UNRESOLVED
    if level == MODIFY and permits_dose_modification(case):
        return EXECUTABLE
    return NON_AVOIDANCE


def source_strength(excerpt: str, drug: str) -> int:
    """Directive strength of the guideline text itself, for the named drug.

    Scores only sentences that mention the drug, so a multi-drug excerpt is not
    credited with a prohibition issued about a different agent.
    """
    if not excerpt:
        return NONE
    sents = [s for s in re.split(r"(?<=[.])\s+", excerpt)
             if drug and re.search(re.escape(drug), s, re.IGNORECASE)]
    return max((directive_strength(s) for s in sents), default=NONE)
