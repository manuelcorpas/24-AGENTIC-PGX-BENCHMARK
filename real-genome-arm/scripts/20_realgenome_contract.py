#!/usr/bin/env python3
"""Strict output contract and phenotype vocabulary for the real-genome arm.

The original parser converted any unlabelled response into an 80-character
"prediction". This module makes output validity explicit. Clinical accuracy is
only evaluated after a response satisfies the supplied one-line contract and
contains a recognised phenotype or functional-status label.
"""
from __future__ import annotations

import re


CANONICAL_PHENOTYPES = {
    "indeterminate",
    "normal metabolizer",
    "intermediate metabolizer",
    "poor metabolizer",
    "rapid metabolizer",
    "ultrarapid metabolizer",
    "likely intermediate metabolizer",
    "likely poor metabolizer",
    "possible intermediate metabolizer",
    "normal function",
    "decreased function",
    "possible decreased function",
    "poor function",
    "increased function",
    "favorable response",
    "unfavorable response",
    "uncertain susceptibility",
    "increased susceptibility",
    "decreased susceptibility",
    "slow acetylator",
    "intermediate acetylator",
    "rapid acetylator",
    "malignant hyperthermia susceptible",
    "not susceptible",
}


ALIASES = {
    "extensive metabolizer": "normal metabolizer",
    "normal activity": "normal function",
    "reduced function": "decreased function",
    "no function": "poor function",
    "favorable response genotype": "favorable response",
    "unfavorable response genotype": "unfavorable response",
    "favorable genotype": "favorable response",
    "unfavorable genotype": "unfavorable response",
    "malignant hyperthermia susceptibility": "malignant hyperthermia susceptible",
    "extensive (normal) metabolizer": "normal metabolizer",
    "non-expresser": "poor metabolizer",
    "non-expressor": "poor metabolizer",
    "normal responder": "favorable response",
    "unfavorable ifnl3 genotype": "unfavorable response",
    "unfavorable treatment response": "unfavorable response",
    "normal dpd activity": "normal metabolizer",
    "normal dpd activity (likely normal metabolizer)": "normal metabolizer",
    "partial dpd deficiency": "intermediate metabolizer",
    "intermediate dpd activity": "intermediate metabolizer",
    "dihydropyrimidine dehydrogenase deficient": "poor metabolizer",
    "deficient metabolizer": "poor metabolizer",
    "low oatp1b1 function": "poor function",
    "uncertain function": "indeterminate",
    "uncertain/unknown significance": "indeterminate",
}


TRUNCATION_REASONS = {
    "length",
    "max_tokens",
    "max token",
    "max tokens",
    "max_output_tokens",
}


def _normalise_text(text: str) -> str:
    value = (text or "").strip().casefold()
    value = value.replace("metaboliser", "metabolizer")
    value = value.replace("favourable", "favorable")
    value = value.replace("unfavourable", "unfavorable")
    value = value.replace("ultra-rapid", "ultrarapid")
    value = value.replace("ultra rapid", "ultrarapid")
    value = value.replace("no-function", "no function")
    value = re.sub(r"\s+", " ", value)
    return value.strip(" .")


def normalise_phenotype(text: str | None, gene: str | None = None) -> str | None:
    """Return one exact canonical label, or None for an unsupported value."""
    value = _normalise_text(text or "")
    if gene:
        gene_value = _normalise_text(gene)
        value = re.sub(rf"^{re.escape(gene_value)}\s*[:\-]?\s+", "", value)
    value = re.sub(
        r"\s*\((?:nm|im|pm|rm|um|expressor|non[- ]?expressor)\)$",
        "",
        value,
    )
    if (
        value == "undefined"
        or value.startswith("not applicable")
        or value.startswith("not defined by cpic")
        or value.startswith("not assigned by cpic")
    ):
        return "indeterminate"
    if re.fullmatch(r"uncertain susceptibility(?: to malignant hyperthermia)?", value):
        return "uncertain susceptibility"
    value = ALIASES.get(value, value)
    return value if value in CANONICAL_PHENOTYPES else None


def stop_indicates_truncation(stop_reason: str | None) -> bool:
    value = _normalise_text(stop_reason or "").replace("-", "_")
    return value in TRUNCATION_REASONS or "max_token" in value


def _labelled_lines(lines: list[str]) -> list[tuple[str, str]]:
    found = []
    for line in lines:
        match = re.fullmatch(r"(?i)PHENOTYPE\s*:\s*(.*?)\s*", line)
        if match:
            found.append((line, match.group(1)))
    return found


def parse_contract_output(
    raw_text: str | None,
    stop_reason: str | None = None,
    error=None,
    gene: str | None = None,
) -> dict:
    """Parse a response without converting prose or errors into predictions.

    The supplied prompt requests exactly one line. A unique phenotype line is
    still extracted from a longer response for audit, but the response remains
    invalid because it did not follow the output contract.
    """
    base = {
        "parse_status": "invalid",
        "invalid_reason": None,
        "contract_status": None,
        "value_status": None,
        "extracted_phenotype": None,
        "phenotype_key": None,
    }
    if error:
        return {**base, "invalid_reason": "provider_error", "contract_status": "provider_error"}

    text = (raw_text or "").strip()
    if not text:
        reason = "truncated_empty_output" if stop_indicates_truncation(stop_reason) else "empty_output"
        return {**base, "invalid_reason": reason, "contract_status": "empty"}

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    labelled = _labelled_lines(lines)
    extracted = labelled[0][1] if len(labelled) == 1 else None

    if len(labelled) > 1:
        return {
            **base,
            "invalid_reason": "multiple_phenotype_lines",
            "contract_status": "multiple_phenotype_lines",
        }

    if not labelled:
        reason = (
            "truncated_without_phenotype"
            if stop_indicates_truncation(stop_reason)
            else "missing_phenotype_line"
        )
        return {**base, "invalid_reason": reason, "contract_status": "missing_phenotype_line"}

    value = labelled[0][1].strip()
    if not value:
        return {
            **base,
            "invalid_reason": "empty_phenotype",
            "contract_status": "empty_phenotype",
        }
    key = normalise_phenotype(value, gene=gene)
    return {
        "parse_status": "valid",
        "invalid_reason": None,
        "contract_status": "conformant" if len(lines) == 1 else "extra_text",
        "value_status": "recognised" if key is not None else "unrecognised",
        "extracted_phenotype": value,
        "phenotype_key": key,
    }


def reference_is_scorable(phenotype: str | None) -> bool:
    """A caller reference is scorable unless it explicitly abstained."""
    key = normalise_phenotype(phenotype)
    if key is None:
        raise ValueError(f"unsupported caller phenotype: {phenotype!r}")
    return key != "indeterminate"
