#!/usr/bin/env python3
"""Concurrent, resumable DeepSeek V4 Flash runner for the three population arms."""
from __future__ import annotations

import argparse
import json
import os
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
CODE = BASE / "code"
OUT_JSONL = BASE / "data" / "v4_flash_population_three_arm.jsonl"
OUT_JSON = BASE / "data" / "v4_flash_population_three_arm.json"
MODEL = "DeepSeek V4 Flash"
HARD_CAP_USD = 15.0


def load(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_tasks() -> list[dict]:
    core = load(CODE / "02-run-benchmark-v3.py", "v4_population_core")
    rag = load(CODE / "16-run-rag-condition.py", "v4_population_rag")
    cases = json.loads(core.SPECS_FILE.read_text())
    corpus = json.loads(rag.CORPUS_FILE.read_text())["genes"]
    tasks = []
    for run in range(3):
        for case in cases:
            for pop in core.POPULATIONS:
                for cond in ("no_spec", "with_spec", "cpic_rag"):
                    if cond == "no_spec":
                        prompt = core.no_spec_prompt(case, pop)
                    elif cond == "with_spec":
                        prompt = core.with_spec_prompt(case, pop)
                    else:
                        prompt = rag.cpic_rag_prompt(
                            case, pop, corpus[case["gene"]]["guideline_excerpt"]
                        )
                    tasks.append({
                        "run": run,
                        "model": MODEL,
                        "tc": case["id"],
                        "gene": case["gene"],
                        "drug": case.get("drug", ""),
                        "pop": pop["id"],
                        "cond": cond,
                        "prompt": prompt,
                        "case": case,
                        "population": pop,
                    })
    return tasks


def unit_key(row: dict) -> tuple:
    return row["run"], row["model"], row["tc"], row["pop"], row["cond"]


def load_latest(path: Path) -> dict[tuple, dict]:
    latest = {}
    if not path.exists():
        return latest
    for line in path.read_text().splitlines():
        if line.strip():
            row = json.loads(line)
            latest[unit_key(row)] = row
    return latest


def append_row(path: Path, row: dict, lock: threading.Lock) -> None:
    with lock:
        with path.open("a") as sink:
            sink.write(json.dumps(row, ensure_ascii=False) + "\n")
            sink.flush()
            os.fsync(sink.fileno())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed-json", type=Path)
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args(argv)

    import openai
    ds_v4 = load(CODE / "94-deepseek-v4-client.py", "v4_population_client")
    key = os.environ.get("DEEPSEEK_API_KEY")
    if not key:
        raise SystemExit("DEEPSEEK_API_KEY is not set")
    client = openai.OpenAI(api_key=key, base_url="https://api.deepseek.com")
    core = load(CODE / "02-run-benchmark-v3.py", "v4_population_scoring_core")
    rag = load(CODE / "16-run-rag-condition.py", "v4_population_scoring_rag")
    tasks = build_tasks()
    lock = threading.Lock()

    latest = load_latest(OUT_JSONL)
    if args.seed_json and args.seed_json.exists() and not latest:
        seed = json.loads(args.seed_json.read_text())
        for row in seed:
            append_row(OUT_JSONL, row, lock)
            latest[unit_key(row)] = row

    pending = [task for task in tasks
               if unit_key(task) not in latest or latest[unit_key(task)].get("error")]
    print(f"population three-arm: {len(tasks)} expected, {len(latest)} retained, {len(pending)} pending")

    def run_one(task: dict) -> dict:
        result = ds_v4.call_deepseek_v4(client, task["prompt"], max_tokens=600)
        parsed = core.parse(result["text"])
        scores = (rag.score_preliminary(parsed, task["case"], task["population"])
                  if task["cond"] == "cpic_rag"
                  else core.score(parsed, task["case"], task["population"]))
        usage = result["usage"]
        return {
            "run": task["run"], "model": MODEL, "tc": task["tc"],
            "gene": task["gene"], "drug": task["drug"], "pop": task["pop"],
            "cond": task["cond"], "parsed": parsed, "scores": scores,
            "in_tok": usage["input_tokens"], "out_tok": usage["output_tokens"],
            "deepseek_provenance": result,
        }

    completed = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, task): task for task in pending}
        for future in as_completed(futures):
            task = futures[future]
            try:
                row = future.result()
            except Exception as exc:  # noqa: BLE001
                row = {
                    "run": task["run"], "model": MODEL, "tc": task["tc"],
                    "gene": task["gene"], "drug": task["drug"], "pop": task["pop"],
                    "cond": task["cond"], "parsed": {}, "scores": {},
                    "in_tok": 0, "out_tok": 0,
                    "error": f"{type(exc).__name__}: {exc}",
                }
                if hasattr(exc, "provenance"):
                    row["deepseek_provenance"] = exc.provenance
            append_row(OUT_JSONL, row, lock)
            latest[unit_key(row)] = row
            completed += 1
            if completed % 100 == 0 or completed == len(pending):
                errors = sum(bool(value.get("error")) for value in latest.values())
                print(f"  {completed}/{len(pending)} new, total {len(latest)}/2970, errors {errors}", flush=True)

    ordered = [latest[unit_key(task)] for task in tasks if unit_key(task) in latest]
    OUT_JSON.write_text(json.dumps(ordered, indent=2))
    failures = [row for row in ordered if row.get("error")]
    cost = sum((row.get("in_tok", 0) * 0.27 + row.get("out_tok", 0) * 1.10) / 1_000_000
               for row in ordered)
    print(f"wrote {len(ordered)} rows, errors {len(failures)}, recorded cost ${cost:.4f}")
    if cost > HARD_CAP_USD:
        raise SystemExit(f"hard spend cap exceeded: ${cost:.2f} > ${HARD_CAP_USD:.2f}")
    return 0 if len(ordered) == 2970 and not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
