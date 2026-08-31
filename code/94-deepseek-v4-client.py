#!/usr/bin/env python3
"""Pinned DeepSeek V4 Flash caller with auditable request provenance."""
from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from typing import Any


MODEL_ID = "deepseek-v4-flash"
MODEL_LABEL = "DeepSeek V4 Flash"
THINKING = {"thinking": {"type": "disabled"}}


class ReturnedModelMismatch(RuntimeError):
    """Raised when the provider does not return the explicitly requested model."""

    def __init__(self, message: str, provenance: dict[str, Any] | None = None):
        super().__init__(message)
        self.provenance = provenance


class DeepSeekCallError(RuntimeError):
    """A failed provider request with its complete attempt record attached."""

    def __init__(self, message: str, provenance: dict[str, Any]):
        super().__init__(message)
        self.provenance = provenance


class ProvenanceText(str):
    """String-compatible model output carrying the provider audit record."""

    def __new__(cls, value: str, provenance: dict[str, Any]):
        obj = super().__new__(cls, value)
        obj.provenance = provenance
        return obj


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _serialise_response(response: Any) -> dict[str, Any]:
    if hasattr(response, "model_dump"):
        try:
            return response.model_dump(mode="json")
        except TypeError:
            return response.model_dump()
    if hasattr(response, "to_dict"):
        return response.to_dict()
    return json.loads(json.dumps(response, default=lambda obj: getattr(obj, "__dict__", str(obj))))


def _usage(response: Any) -> dict[str, int]:
    usage = getattr(response, "usage", None)
    return {
        "input_tokens": int(getattr(usage, "prompt_tokens", 0) or 0),
        "output_tokens": int(getattr(usage, "completion_tokens", 0) or 0),
        "total_tokens": int(getattr(usage, "total_tokens", 0) or 0),
    }


def _error_record(exc: Exception) -> dict[str, Any]:
    record = {"type": type(exc).__name__, "message": str(exc)[:2000]}
    for attr in ("status_code", "code", "request_id"):
        value = getattr(exc, attr, None)
        if value is not None:
            record[attr] = value
    response = getattr(exc, "response", None)
    if response is not None:
        record["response"] = {
            "status_code": getattr(response, "status_code", None),
            "text": getattr(response, "text", "")[:5000],
        }
    return record


def call_deepseek_v4(
    client: Any,
    prompt: str,
    *,
    max_tokens: int,
    max_attempts: int = 3,
    retry_base_seconds: float = 2.0,
) -> dict[str, Any]:
    """Call pinned V4 Flash in non-thinking mode and retain full provenance."""
    request_payload = {
        "model": MODEL_ID,
        "max_tokens": max_tokens,
        "messages": [{"role": "user", "content": prompt}],
        "extra_body": THINKING,
    }
    requested_at = utc_now()
    start = time.monotonic()
    attempt_log = []

    for attempt in range(1, max_attempts + 1):
        attempt_started = utc_now()
        try:
            response = client.chat.completions.create(**request_payload)
            returned_model = getattr(response, "model", None)
            attempt_log.append({
                "attempt": attempt,
                "started_at_utc": attempt_started,
                "completed_at_utc": utc_now(),
                "status": "success",
            })
            if returned_model != MODEL_ID:
                choice = response.choices[0]
                mismatch_provenance = {
                    "text": choice.message.content or "",
                    "requested_model": MODEL_ID,
                    "returned_model": returned_model,
                    "thinking_mode": "disabled",
                    "request_payload": request_payload,
                    "provider_response": _serialise_response(response),
                    "response_id": getattr(response, "id", None),
                    "provider_created": getattr(response, "created", None),
                    "stop_reason": getattr(choice, "finish_reason", None),
                    "usage": _usage(response),
                    "requested_at_utc": requested_at,
                    "completed_at_utc": utc_now(),
                    "latency_ms": round((time.monotonic() - start) * 1000),
                    "attempt_count": len(attempt_log),
                    "attempt_log": attempt_log,
                    "error": {
                        "type": "ReturnedModelMismatch",
                        "message": f"requested {MODEL_ID}, provider returned {returned_model!r}",
                    },
                }
                raise ReturnedModelMismatch(
                    f"requested {MODEL_ID}, provider returned {returned_model!r}",
                    mismatch_provenance,
                )
            choice = response.choices[0]
            usage = _usage(response)
            completed_at = utc_now()
            return {
                "text": choice.message.content or "",
                "requested_model": MODEL_ID,
                "returned_model": returned_model,
                "thinking_mode": "disabled",
                "request_payload": request_payload,
                "provider_response": _serialise_response(response),
                "response_id": getattr(response, "id", None),
                "provider_created": getattr(response, "created", None),
                "stop_reason": getattr(choice, "finish_reason", None),
                "usage": usage,
                "requested_at_utc": requested_at,
                "completed_at_utc": completed_at,
                "latency_ms": round((time.monotonic() - start) * 1000),
                "attempt_count": len(attempt_log),
                "attempt_log": attempt_log,
                "error": None,
            }
        except ReturnedModelMismatch:
            raise
        except Exception as exc:  # noqa: BLE001
            attempt_log.append({
                "attempt": attempt,
                "started_at_utc": attempt_started,
                "completed_at_utc": utc_now(),
                "status": "error",
                "error": _error_record(exc),
            })
            if attempt < max_attempts:
                time.sleep(retry_base_seconds * (2 ** (attempt - 1)))
                continue
            provenance = {
                "requested_model": MODEL_ID,
                "returned_model": None,
                "thinking_mode": "disabled",
                "request_payload": request_payload,
                "provider_response": None,
                "response_id": None,
                "stop_reason": None,
                "usage": {"input_tokens": 0, "output_tokens": 0, "total_tokens": 0},
                "requested_at_utc": requested_at,
                "completed_at_utc": utc_now(),
                "latency_ms": round((time.monotonic() - start) * 1000),
                "attempt_count": len(attempt_log),
                "attempt_log": attempt_log,
                "error": _error_record(exc),
            }
            raise DeepSeekCallError(str(exc), provenance) from exc

    raise AssertionError("unreachable")


def as_text(result: dict[str, Any]) -> ProvenanceText:
    """Return a string-compatible response without discarding provenance."""
    return ProvenanceText(result["text"], result)
