from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "code" / "96-run-deepseek-v4-population.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("deepseek_v4_population", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_population_runner_builds_exact_three_arm_manifest():
    mod = load_runner()
    tasks = mod.build_tasks()
    assert len(tasks) == 2970
    keys = {mod.unit_key(task) for task in tasks}
    assert len(keys) == 2970
    assert {task["cond"] for task in tasks} == {"no_spec", "with_spec", "cpic_rag"}
    assert all(task["model"] == "DeepSeek V4 Flash" for task in tasks)
