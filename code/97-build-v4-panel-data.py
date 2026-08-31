#!/usr/bin/env python3
"""Replace only the DeepSeek arm in each pooled paper dataset."""
from __future__ import annotations

import json
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
OLD_LABEL = "DeepSeek V3"
NEW_LABEL = "DeepSeek V4 Flash"

SPECS = {
    "five_cell": (
        "v3_five_cell_live.json", "v4_flash_five_cell.json",
        "v4_flash_eight_model_five_cell.json", False,
    ),
    "extracted_rules": (
        "v3_extracted_rules.json", "v4_flash_extracted_rules.json",
        "v4_flash_eight_model_extracted_rules.json", False,
    ),
    "adversarial_forward": (
        "v3_adversarial_scrambled.json", "v4_flash_adversarial_forward.json",
        "v4_flash_adversarial_forward_panel.json", False,
    ),
    "adversarial_reverse": (
        "v3_adversarial_reverse.json", "v4_flash_adversarial_reverse.json",
        "v4_flash_adversarial_reverse_panel.json", False,
    ),
    "drug_keyed": (
        "v3_rag_drugkeyed_fullgrid.json", "v4_flash_rag_drugkeyed.json",
        "v4_flash_eight_model_rag_drugkeyed.json", False,
    ),
    "population_three_arm": (
        "v3_raw_rescored_three_arm.json", "v4_flash_population_three_arm.json",
        "v4_flash_population_three_arm_panel.json", False,
    ),
    "population_skill_arms": (
        "v3_armAB_fullgrid.json", "v4_flash_population_skill_arms.json",
        "v4_flash_population_skill_arms_panel.json", False,
    ),
    "real_genome": (
        "v3_realgenome_provider_responses_4cohorts.jsonl",
        "v4_flash_realgenome_provider_responses.jsonl",
        "v4_flash_eight_model_realgenome_provider_responses.jsonl", True,
    ),
}


def read_rows(path: Path, jsonl: bool = False) -> list[dict]:
    if jsonl:
        return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    return json.loads(path.read_text())


def merge(old_path: Path, new_path: Path, jsonl: bool = False) -> list[dict]:
    old = read_rows(old_path, jsonl)
    new = read_rows(new_path, jsonl)
    if any(row.get("model") != NEW_LABEL for row in new):
        raise AssertionError(f"non-V4 row in replacement source {new_path}")
    retained = [row for row in old if row.get("model") != OLD_LABEL]
    if len(retained) + len(new) != len(old):
        raise AssertionError(
            f"replacement cardinality changed for {old_path}: "
            f"old={len(old)}, retained={len(retained)}, new={len(new)}"
        )
    rows = retained + new
    if any(row.get("model") == OLD_LABEL for row in rows):
        raise AssertionError(f"retired label survived merge for {old_path}")
    return rows


def build_outputs() -> dict[str, list[dict]]:
    return {
        name: merge(DATA / old, DATA / new, jsonl)
        for name, (old, new, _, jsonl) in SPECS.items()
    }


def main() -> int:
    outputs = build_outputs()
    for name, rows in outputs.items():
        _, _, out_name, jsonl = SPECS[name]
        out = DATA / out_name
        if out.exists():
            raise SystemExit(f"refusing to overwrite existing panel data: {out}")
        if jsonl:
            with out.open("x") as sink:
                for row in rows:
                    sink.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
        else:
            out.write_text(json.dumps(rows, ensure_ascii=False, separators=(",", ":")))
        print(name, len(rows), out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
