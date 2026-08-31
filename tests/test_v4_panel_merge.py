from __future__ import annotations

import importlib.util
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "code" / "97-build-v4-panel-data.py"


def load_module():
    spec = importlib.util.spec_from_file_location("v4_panel_merge", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_merge_replaces_only_deepseek_and_preserves_panel_denominators():
    mod = load_module()
    outputs = mod.build_outputs()
    expected = {
        "five_cell": 13199,
        "extracted_rules": 168,
        "adversarial_forward": 45,
        "adversarial_reverse": 45,
        "drug_keyed": 5280,
        "population_three_arm": 26730,
        "population_skill_arms": 17820,
        "real_genome": 4288,
    }
    assert {name: len(rows) for name, rows in outputs.items()} == expected
    for rows in outputs.values():
        models = Counter(row["model"] for row in rows)
        assert "DeepSeek V3" not in models
        assert models["DeepSeek V4 Flash"] > 0
