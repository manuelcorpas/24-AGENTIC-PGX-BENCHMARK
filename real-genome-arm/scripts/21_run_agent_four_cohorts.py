#!/usr/bin/env python3
"""Run all real-genome states with raw provider-response provenance.

Every completed request is appended to JSONL immediately. Successful rows are
resumed rather than repeated, and the file is compacted to one row per
model-state pair at the end. The exact request payload, full provider response,
raw text, stop reason, usage, response ID, timing and retry history are kept.

Usage:
  21_run_agent_four_cohorts.py data/n5_four_cohorts.tsv \
      data/v3_realgenome_provider_responses_4cohorts.jsonl
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


HERE = Path(__file__).resolve().parent
ARM = HERE.parent
DEFAULT_PANEL = ARM / "config" / "models.txt"
PROVIDERS = {"anthropic", "openai", "openai-reasoning", "google", "deepseek"}
SCHEMA_VERSION = "realgenome-provider-response-v1"
DEEPSEEK_V4_ID = "deepseek-v4-flash"
DEEPSEEK_V4_OUTPUT = "data/v4_flash_realgenome_provider_responses.jsonl"


def _load_module(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


contract = _load_module(HERE / "20_realgenome_contract.py", "realgenome_contract_runner")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def build_prompt(gene: str, diplotype: str) -> str:
    """Preserve the exact prompt used by the original real-genome runner."""
    return (
        "You are executing a ClawBio pharmacogenomics skill. Apply CPIC Level A guidelines. "
        "Given the gene and the patient's star-allele diplotype, output the metaboliser phenotype "
        "(or CPIC functional status). If the diplotype contains an allele of uncertain or unknown "
        "function, output 'Indeterminate'.\n\n"
        f"Gene: {gene}\nDiplotype: {diplotype}\n\n"
        "Output one line:\nPHENOTYPE: [CPIC phenotype]"
    )


def load_panel(path: Path) -> list[tuple[str, str, str]]:
    panel = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"\s{2,}", line)
        if len(parts) < 3:
            raise ValueError(f"malformed model-panel line: {line!r}")
        name, model_id, provider = parts[:3]
        if provider not in PROVIDERS:
            raise ValueError(f"unknown provider {provider!r} for {name!r}")
        panel.append((name, model_id, provider))
    if not panel:
        raise ValueError(f"no models parsed from {path}")
    if len({name for name, _, _ in panel}) != len(panel):
        raise ValueError("duplicate display name in model panel")
    return panel


def read_states(path: Path) -> list[dict]:
    states, seen = [], set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["cohort"], row["gene"], row["diplotype"])
            if key in seen:
                continue
            seen.add(key)
            states.append(
                {
                    "cohort": row["cohort"],
                    "gene": row["gene"],
                    "diplotype": row["diplotype"],
                    "caller_phenotype": row["phenotype"],
                }
            )
    return states


def extract_provider_result(provider: str, body: dict) -> dict:
    """Extract common fields while retaining the caller's full raw body."""
    if provider == "anthropic":
        texts = [part.get("text", "") for part in body.get("content", []) if part.get("type") == "text"]
        usage = body.get("usage") or {}
        return {
            "raw_text": "".join(texts),
            "provider_stop_reason": body.get("stop_reason"),
            "response_id": body.get("id"),
            "model_version_returned": body.get("model"),
            "input_tokens": usage.get("input_tokens", 0),
            "output_tokens": usage.get("output_tokens", 0),
        }

    if provider in {"openai", "openai-reasoning", "deepseek"}:
        choices = body.get("choices") or []
        choice = choices[0] if choices else {}
        message = choice.get("message") or {}
        content = message.get("content")
        if isinstance(content, list):
            content = "".join(
                part.get("text", "") for part in content if isinstance(part, dict)
            )
        usage = body.get("usage") or {}
        return {
            "raw_text": content or "",
            "provider_stop_reason": choice.get("finish_reason"),
            "response_id": body.get("id"),
            "model_version_returned": body.get("model"),
            "input_tokens": usage.get("prompt_tokens", 0),
            "output_tokens": usage.get("completion_tokens", 0),
        }

    if provider == "google":
        candidates = body.get("candidates") or []
        candidate = candidates[0] if candidates else {}
        parts = (candidate.get("content") or {}).get("parts") or []
        usage = body.get("usageMetadata") or {}
        return {
            "raw_text": "".join(part.get("text", "") for part in parts),
            "provider_stop_reason": candidate.get("finishReason"),
            "response_id": body.get("responseId"),
            "model_version_returned": body.get("modelVersion"),
            "input_tokens": usage.get("promptTokenCount", 0),
            "output_tokens": usage.get("candidatesTokenCount", 0),
        }

    raise ValueError(f"unsupported provider: {provider}")


def request_payload(
    model_id: str,
    provider: str,
    prompt: str,
    standard_tokens: int,
    reasoning_tokens: int,
    google_tokens: int,
) -> dict:
    if provider == "anthropic":
        return {
            "model": model_id,
            "max_tokens": standard_tokens,
            "messages": [{"role": "user", "content": prompt}],
        }
    if provider in {"openai", "deepseek"}:
        payload = {
            "model": model_id,
            "max_tokens": standard_tokens,
            "messages": [{"role": "user", "content": prompt}],
        }
        if provider == "deepseek":
            if model_id != DEEPSEEK_V4_ID:
                raise ValueError(f"DeepSeek must be pinned to {DEEPSEEK_V4_ID}, got {model_id}")
            payload["thinking"] = {"type": "disabled"}
        return payload
    if provider == "openai-reasoning":
        return {
            "model": model_id,
            "max_completion_tokens": reasoning_tokens,
            "messages": [{"role": "user", "content": prompt}],
        }
    if provider == "google":
        return {
            "contents": [{"parts": [{"text": prompt}]}],
            "generationConfig": {"maxOutputTokens": google_tokens},
        }
    raise ValueError(provider)


class ProviderCallError(Exception):
    def __init__(self, kind: str, message: str, status=None, body=None, retry_after=None):
        super().__init__(message)
        self.details = {
            "type": kind,
            "message": message[:500],
            "status": status,
            "provider_body": body,
            "retry_after_seconds": retry_after,
        }


def _decode_body(raw: bytes):
    text = raw.decode(errors="replace")
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return {"unparsed_body": text[:2000]}


def _post(url: str, payload: dict, headers: dict, timeout: int) -> dict:
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", **headers},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read())
    except urllib.error.HTTPError as exc:
        body = _decode_body(exc.read())
        retry_after = exc.headers.get("Retry-After") if exc.headers else None
        try:
            retry_after = float(retry_after) if retry_after is not None else None
        except ValueError:
            retry_after = None
        raise ProviderCallError(
            "HTTPError",
            f"provider returned HTTP {exc.code}",
            status=exc.code,
            body=body,
            retry_after=retry_after,
        ) from None
    except (urllib.error.URLError, TimeoutError) as exc:
        raise ProviderCallError(type(exc).__name__, str(getattr(exc, "reason", exc))[:500]) from None


def call_provider(model_id: str, provider: str, payload: dict, timeout: int) -> dict:
    if provider == "anthropic":
        return _post(
            "https://api.anthropic.com/v1/messages",
            payload,
            {
                "x-api-key": os.environ["ANTHROPIC_API_KEY"],
                "anthropic-version": "2023-06-01",
            },
            timeout,
        )
    if provider in {"openai", "openai-reasoning"}:
        return _post(
            "https://api.openai.com/v1/chat/completions",
            payload,
            {"Authorization": "Bearer " + os.environ["OPENAI_API_KEY"]},
            timeout,
        )
    if provider == "deepseek":
        return _post(
            "https://api.deepseek.com/chat/completions",
            payload,
            {"Authorization": "Bearer " + os.environ["DEEPSEEK_API_KEY"]},
            timeout,
        )
    if provider == "google":
        model = payload.pop("_model_id")
        key = os.environ["GOOGLE_API_KEY"]
        return _post(
            f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={key}",
            payload,
            {},
            timeout,
        )
    raise ValueError(provider)


def is_retryable(error: dict) -> bool:
    status = error.get("status")
    return status in {408, 409, 429, 500, 502, 503, 504, 529} or status is None


def provider_group(provider: str) -> str:
    """Providers that share credentials also share one concurrency budget."""
    return "openai" if provider in {"openai", "openai-reasoning"} else provider


def job_key(row: dict) -> tuple[str, str, str, str]:
    return row["model"], row["cohort"], row["gene"], row["diplotype"]


def run_one(
    job: tuple,
    standard_tokens: int,
    reasoning_tokens: int,
    google_tokens: int,
    timeout: int,
    max_attempts: int,
    retry_base: float,
) -> dict:
    name, model_id, provider, state = job
    prompt = build_prompt(state["gene"], state["diplotype"])
    payload = request_payload(
        model_id,
        provider,
        prompt,
        standard_tokens,
        reasoning_tokens,
        google_tokens,
    )
    request_for_storage = json.loads(json.dumps(payload))
    call_payload = json.loads(json.dumps(payload))
    if provider == "google":
        call_payload["_model_id"] = model_id

    started = utc_now()
    started_clock = time.monotonic()
    attempt_log = []
    body = None
    extracted = {
        "raw_text": None,
        "provider_stop_reason": None,
        "response_id": None,
        "model_version_returned": None,
        "input_tokens": 0,
        "output_tokens": 0,
    }
    final_error = None

    for attempt in range(1, max_attempts + 1):
        attempt_started = utc_now()
        try:
            body = call_provider(model_id, provider, json.loads(json.dumps(call_payload)), timeout)
            extracted = extract_provider_result(provider, body)
            if provider == "deepseek" and extracted["model_version_returned"] != DEEPSEEK_V4_ID:
                raise RuntimeError(
                    f"requested {DEEPSEEK_V4_ID}, provider returned "
                    f"{extracted['model_version_returned']!r}"
                )
            attempt_log.append(
                {"attempt": attempt, "started_at_utc": attempt_started, "status": "success"}
            )
            final_error = None
            break
        except ProviderCallError as exc:
            final_error = exc.details
        except Exception as exc:  # noqa: BLE001
            final_error = {"type": type(exc).__name__, "message": str(exc)[:500], "status": None}
        attempt_log.append(
            {
                "attempt": attempt,
                "started_at_utc": attempt_started,
                "status": "error",
                "error": final_error,
            }
        )
        if attempt >= max_attempts or not is_retryable(final_error):
            break
        delay = final_error.get("retry_after_seconds") or retry_base * (2 ** (attempt - 1))
        time.sleep(min(float(delay), 60.0))

    completed = utc_now()
    parsed = contract.parse_contract_output(
        extracted["raw_text"],
        extracted["provider_stop_reason"],
        final_error,
        gene=state["gene"],
    )
    return {
        "schema_version": SCHEMA_VERSION,
        **state,
        "model": name,
        "model_id": model_id,
        "provider": provider,
        "prompt": prompt,
        "prompt_sha256": sha256_text(prompt),
        "request_payload": request_for_storage,
        "requested_at_utc": started,
        "completed_at_utc": completed,
        "latency_ms": round((time.monotonic() - started_clock) * 1000),
        "attempt_count": len(attempt_log),
        "attempt_log": attempt_log,
        "provider_response": body,
        **extracted,
        "error": final_error,
        "contract_parse": parsed,
    }


def load_jsonl(path: Path) -> dict[tuple[str, str, str, str], dict]:
    latest = {}
    if not path.exists():
        return latest
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"invalid JSONL at line {line_number}: {exc}") from exc
            latest[job_key(row)] = row
    return latest


def append_jsonl(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def compact_jsonl(path: Path, rows: list[dict]) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
    temp.replace(path)


def check_keys(panel: list[tuple[str, str, str]]) -> None:
    env_by_provider = {
        "anthropic": "ANTHROPIC_API_KEY",
        "openai": "OPENAI_API_KEY",
        "openai-reasoning": "OPENAI_API_KEY",
        "google": "GOOGLE_API_KEY",
        "deepseek": "DEEPSEEK_API_KEY",
    }
    missing = sorted({env_by_provider[provider] for _, _, provider in panel if not os.environ.get(env_by_provider[provider])})
    if missing:
        raise SystemExit("missing API keys: " + ", ".join(missing))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("states_tsv", type=Path)
    parser.add_argument("out_jsonl", type=Path)
    parser.add_argument("--models", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--only-model", action="append", default=[])
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--standard-max-tokens", type=int, default=1024)
    parser.add_argument("--reasoning-max-tokens", type=int, default=4096)
    parser.add_argument("--google-max-tokens", type=int, default=4096)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--max-attempts", type=int, default=5)
    parser.add_argument("--retry-base-seconds", type=float, default=2.0)
    args = parser.parse_args(argv)

    states = read_states(args.states_tsv)
    if args.limit is not None:
        states = states[: args.limit]
    panel = load_panel(args.models)
    if args.only_model:
        wanted = set(args.only_model)
        panel = [entry for entry in panel if entry[0] in wanted]
        missing_names = wanted - {entry[0] for entry in panel}
        if missing_names:
            raise SystemExit("unknown model names: " + ", ".join(sorted(missing_names)))
    check_keys(panel)

    expected_order = [
        (name, model_id, provider, state)
        for name, model_id, provider in panel
        for state in states
    ]
    latest = load_jsonl(args.out_jsonl)
    expected = len(expected_order)
    successful_existing = sum(
        1
        for name, _, _, state in expected_order
        if (name, state["cohort"], state["gene"], state["diplotype"]) in latest
        and latest[(name, state["cohort"], state["gene"], state["diplotype"])].get("error") is None
    )
    print(
        f"{len(states)} states x {len(panel)} models = {expected} rows; "
        f"resuming {successful_existing} successful rows",
        flush=True,
    )

    jobs = []
    pending_by_model = Counter()
    for state in states:
        for name, model_id, provider in panel:
            key = (name, state["cohort"], state["gene"], state["diplotype"])
            if key in latest and latest[key].get("error") is None:
                continue
            jobs.append((name, model_id, provider, state))
            pending_by_model[name] += 1
    if not jobs:
        print("  all model-state rows are already complete", flush=True)
    else:
        for name, _, _ in panel:
            print(f"  {name}: {pending_by_model[name]} rows pending", flush=True)
        groups = sorted({provider_group(provider) for _, _, provider in panel})
        executors = {
            group: ThreadPoolExecutor(max_workers=args.workers)
            for group in groups
        }
        futures = {}
        try:
            for job in jobs:
                provider = job[2]
                future = executors[provider_group(provider)].submit(
                    run_one,
                    job,
                    args.standard_max_tokens,
                    args.reasoning_max_tokens,
                    args.google_max_tokens,
                    args.timeout,
                    args.max_attempts,
                    args.retry_base_seconds,
                )
                futures[future] = job[0]
            completed_by_model = Counter()
            for index, future in enumerate(as_completed(futures), 1):
                row = future.result()
                append_jsonl(args.out_jsonl, row)
                latest[job_key(row)] = row
                completed_by_model[row["model"]] += 1
                if index % 50 == 0 or index == len(jobs):
                    failures = sum(1 for value in latest.values() if value.get("error"))
                    progress = ", ".join(
                        f"{name} {completed_by_model[name]}/{pending_by_model[name]}"
                        for name, _, _ in panel
                    )
                    print(
                        f"    {index}/{len(jobs)} new rows; retained errors={failures}; {progress}",
                        flush=True,
                    )
        finally:
            for executor in executors.values():
                executor.shutdown(wait=True)

    ordered_rows = []
    for name, _, _, state in expected_order:
        key = (name, state["cohort"], state["gene"], state["diplotype"])
        if key in latest:
            row = latest[key]
            row["contract_parse"] = contract.parse_contract_output(
                row.get("raw_text"),
                row.get("provider_stop_reason"),
                row.get("error"),
                gene=row.get("gene"),
            )
            ordered_rows.append(row)
    compact_jsonl(args.out_jsonl, ordered_rows)

    errors = [row for row in ordered_rows if row.get("error")]
    print(f"wrote {len(ordered_rows)} unique rows to {args.out_jsonl}", flush=True)
    if len(ordered_rows) != expected:
        print(f"FAIL: expected {expected} rows", file=sys.stderr)
        return 1
    if errors:
        print(f"FAIL: {len(errors)} provider errors remain", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
