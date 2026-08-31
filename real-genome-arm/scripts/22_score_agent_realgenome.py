#!/usr/bin/env python3
"""Score provenance-preserving real-genome responses in four categories."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_RAW = REPO / "data" / "v3_realgenome_provider_responses_4cohorts.jsonl"
DEFAULT_STATES = REPO / "data" / "n5_four_cohorts.tsv"
DEFAULT_SCORED = REPO / "data" / "v3_realgenome_scored_4cohorts.jsonl"
DEFAULT_OUT = REPO / "data" / "v3_agent_vs_deterministic.json"
DEFAULT_AUDIT = REPO / "data" / "v3_realgenome_model_audit.json"
DEFAULT_REPORT = REPO / "data" / "v3_realgenome_model_audit.txt"
OUTCOMES = ("correct", "abstained", "clinically_wrong", "invalid_output")
RETURNED_MODEL_LABELS = {"deepseek-v4-flash": "DeepSeek V4 Flash"}


def _load_module(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


contract = _load_module(HERE / "20_realgenome_contract.py", "realgenome_contract_scorer")
rules = _load_module(REPO / "code" / "_pgx_rules.py", "realgenome_rules_scorer")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_jsonl(path: Path) -> list[dict]:
    rows = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f"invalid JSONL at line {line_number}: {exc}") from exc
    return rows


def read_states(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def score_record(record: dict) -> dict:
    reference = record.get("caller_phenotype")
    reference_key = contract.normalise_phenotype(reference)
    if reference_key is None:
        raise ValueError(f"unsupported caller phenotype: {reference!r}")

    parsed = contract.parse_contract_output(
        record.get("raw_text"),
        record.get("provider_stop_reason"),
        record.get("error"),
        gene=record.get("gene"),
    )
    if reference_key == "indeterminate":
        outcome = "excluded_reference"
    elif parsed["parse_status"] != "valid":
        outcome = "invalid_output"
    elif parsed["phenotype_key"] == "indeterminate":
        outcome = "abstained"
    elif parsed["phenotype_key"] == reference_key:
        outcome = "correct"
    else:
        outcome = "clinically_wrong"

    return {
        **record,
        "reference_key": reference_key,
        **parsed,
        "outcome": outcome,
    }


def summarise(records: list[dict]) -> dict:
    included = [row for row in records if row.get("outcome") in OUTCOMES]
    counts = Counter(row["outcome"] for row in included)
    n = len(included)
    valid_non_abstaining = counts["correct"] + counts["clinically_wrong"]
    return {
        "n": n,
        **{outcome: counts[outcome] for outcome in OUTCOMES},
        "excluded_reference": sum(row.get("outcome") == "excluded_reference" for row in records),
        "valid_output_rate": round((n - counts["invalid_output"]) / n, 6) if n else None,
        "accuracy_overall": round(counts["correct"] / n, 6) if n else None,
        "accuracy_among_valid_non_abstaining": (
            round(counts["correct"] / valid_non_abstaining, 6)
            if valid_non_abstaining
            else None
        ),
    }


def summarise_deterministic(counts: dict) -> dict:
    n = counts["n"]
    emitted = counts["emitted"]
    correct = counts["correct"]
    return {
        "n": n,
        "emitted": emitted,
        "correct": correct,
        "abstained": n - emitted,
        "coverage": emitted / n if n else None,
        "accuracy_among_emitted": correct / emitted if emitted else None,
        "accuracy_overall": correct / n if n else None,
    }


def deterministic_arm(states: list[dict], dip2phen: dict) -> dict:
    """Score the validated mapping on the same corrected reference states."""
    per_cohort = defaultdict(lambda: {"n": 0, "emitted": 0, "correct": 0})
    for state in states:
        gene = state["gene"]
        reference = contract.normalise_phenotype(state.get("phenotype"), gene=gene)
        if reference in {None, "indeterminate"}:
            continue
        cohort = state["cohort"]
        counts = per_cohort[cohort]
        counts["n"] += 1
        key = (gene, rules.norm_dip(state["diplotype"], gene))
        predicted = dip2phen.get(key)
        if predicted is None:
            continue
        counts["emitted"] += 1
        if contract.normalise_phenotype(predicted, gene=gene) == reference:
            counts["correct"] += 1
    return {
        cohort: summarise_deterministic(counts)
        for cohort, counts in sorted(per_cohort.items())
    }


def aggregate(scored: list[dict]) -> tuple[dict, dict]:
    by_model_cohort = defaultdict(list)
    by_cohort = defaultdict(list)
    for row in scored:
        by_model_cohort[(row["model"], row["cohort"])].append(row)
        by_cohort[row["cohort"]].append(row)

    models = sorted({row["model"] for row in scored})
    cohorts = sorted({row["cohort"] for row in scored})
    agent_by_model = {
        model: {
            cohort: summarise(by_model_cohort[(model, cohort)])
            for cohort in cohorts
        }
        for model in models
    }
    agent_pooled = {cohort: summarise(by_cohort[cohort]) for cohort in cohorts}
    return agent_pooled, agent_by_model


def validate_complete(rows: list[dict]) -> dict:
    keys = [
        (row["model"], row["cohort"], row["gene"], row["diplotype"])
        for row in rows
    ]
    duplicates = [key for key, count in Counter(keys).items() if count > 1]
    if duplicates:
        raise ValueError(f"duplicate model-state rows: {duplicates[:3]}")
    models = sorted({row["model"] for row in rows})
    states = sorted({(row["cohort"], row["gene"], row["diplotype"]) for row in rows})
    expected = len(models) * len(states)
    if len(rows) != expected:
        raise ValueError(f"incomplete panel: {len(rows)} rows, expected {expected}")
    missing_raw = [row for row in rows if row.get("error") is None and row.get("raw_text") is None]
    if missing_raw:
        raise ValueError(f"{len(missing_raw)} successful rows have no raw_text field")
    missing_stop = [
        row
        for row in rows
        if row.get("error") is None and row.get("provider_stop_reason") is None
    ]
    if missing_stop:
        raise ValueError(f"{len(missing_stop)} successful rows have no provider stop reason")
    return {"rows": len(rows), "models": len(models), "states": len(states)}


def returned_model_matches(requested: str | None, returned: str | None) -> bool:
    """Accept an exact ID or a dated snapshot of the requested model alias."""
    if not requested or not returned:
        return False
    return returned == requested or returned.startswith(requested + "-")


def audit(scored: list[dict]) -> dict:
    result = {}
    for model in sorted({row["model"] for row in scored}):
        rows = [row for row in scored if row["model"] == model]
        included = [row for row in rows if row["outcome"] in OUTCOMES]
        requested_ids = sorted({row["model_id"] for row in rows})
        returned_ids = sorted({str(row.get("model_version_returned")) for row in rows})
        version_match = all(
            returned_model_matches(row.get("model_id"), row.get("model_version_returned"))
            for row in rows
        )
        analysis_label = model
        if not version_match and len(returned_ids) == 1:
            analysis_label = RETURNED_MODEL_LABELS.get(returned_ids[0], model)
        result[model] = {
            "model_id": requested_ids,
            "returned_model_id": returned_ids,
            "version_match": version_match,
            "analysis_label": analysis_label,
            "provider": sorted({row["provider"] for row in rows}),
            "rows": len(rows),
            "scorable_rows": len(included),
            "outcomes": dict(Counter(row["outcome"] for row in included)),
            "invalid_reasons": dict(
                Counter(
                    row["invalid_reason"]
                    for row in included
                    if row["outcome"] == "invalid_output"
                )
            ),
            "contract_status": dict(Counter(row.get("contract_status") for row in included)),
            "unrecognised_values": sum(
                row.get("value_status") == "unrecognised" for row in included
            ),
            "stop_reasons": dict(Counter(str(row.get("provider_stop_reason")) for row in rows)),
            "provider_errors": sum(bool(row.get("error")) for row in rows),
            "retried_rows": sum((row.get("attempt_count") or 1) > 1 for row in rows),
            "input_tokens": sum(row.get("input_tokens") or 0 for row in rows),
            "output_tokens": sum(row.get("output_tokens") or 0 for row in rows),
        }
    return result


def report_text(audit_data: dict, agent_by_model: dict) -> str:
    lines = ["REAL-GENOME MODEL AUDIT", ""]
    for model, details in audit_data.items():
        lines.append(model)
        lines.append(
            f"  rows={details['rows']} scorable={details['scorable_rows']} "
            f"provider_errors={details['provider_errors']} retried={details['retried_rows']}"
        )
        lines.append(
            f"  requested={details['model_id']} returned={details['returned_model_id']} "
            f"version_match={details['version_match']} analysis_label={details['analysis_label']}"
        )
        lines.append(f"  stop_reasons={details['stop_reasons']}")
        lines.append(f"  invalid_reasons={details['invalid_reasons']}")
        for cohort, summary in agent_by_model[model].items():
            lines.append(
                f"  {cohort}: n={summary['n']} correct={summary['correct']} "
                f"abstained={summary['abstained']} wrong={summary['clinically_wrong']} "
                f"invalid={summary['invalid_output']}"
            )
        lines.append("")
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--states", type=Path, default=DEFAULT_STATES)
    parser.add_argument("--scored", type=Path, default=DEFAULT_SCORED)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args(argv)

    raw = read_jsonl(args.raw)
    states = read_states(args.states)
    dimensions = validate_complete(raw)
    scored = [score_record(row) for row in raw]
    agent_pooled, agent_by_model = aggregate(scored)
    dip2phen, _, _, _ = rules.build_rules(rules.load_cases())
    deterministic = deterministic_arm(states, dip2phen)
    audit_data = audit(scored)
    metadata = {
        "schema_version": "realgenome-four-outcome-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "raw_source": str(args.raw.relative_to(REPO)),
        "raw_sha256": file_sha256(args.raw),
        "state_source": str(args.states.relative_to(REPO)),
        "state_sha256": file_sha256(args.states),
        **dimensions,
        "outcomes": list(OUTCOMES),
        "version_matched_models": [
            model for model, details in audit_data.items() if details["version_match"]
        ],
        "version_mismatches": {
            model: {
                "requested": details["model_id"],
                "returned": details["returned_model_id"],
                "analysis_label": details["analysis_label"],
            }
            for model, details in audit_data.items()
            if not details["version_match"]
        },
    }
    aggregate_data = {
        "metadata": metadata,
        "model_provenance": {
            model: {
                "requested_model_id": details["model_id"],
                "returned_model_id": details["returned_model_id"],
                "version_match": details["version_match"],
                "analysis_label": details["analysis_label"],
            }
            for model, details in audit_data.items()
        },
        "deterministic": deterministic,
        "agent_pooled": agent_pooled,
        "agent_by_model": agent_by_model,
    }

    args.scored.write_text(
        "".join(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n" for row in scored)
    )
    args.out.write_text(json.dumps(aggregate_data, indent=2) + "\n")
    args.audit.write_text(json.dumps({"metadata": metadata, "models": audit_data}, indent=2) + "\n")
    args.report.write_text(report_text(audit_data, agent_by_model) + "\n")
    print(report_text(audit_data, agent_by_model))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
