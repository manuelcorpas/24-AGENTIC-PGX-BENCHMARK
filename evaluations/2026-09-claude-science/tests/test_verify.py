"""verify.py is the one command an outsider runs. These tests pin what it must
reproduce from the shipped CSVs, so a silent change to the scorer, the data or
the id map cannot leave the headline numbers standing.
"""
from __future__ import annotations

import csv
import json
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]


def load():
    spec = spec_from_file_location("verify", HERE / "verify.py")
    mod = module_from_spec(spec)
    sys.modules["verify"] = mod
    spec.loader.exec_module(mod)
    return mod


v = load()


def _write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case_id", "diplotype", "phenotype", "recommendation", "abstained"])
        w.writeheader()
        w.writerows(rows)


def test_hermetic_gold_answers_score_one_and_wrong_answers_show_spread(tmp_path):
    cases = v.load_cases()[:12]
    id_map = {f"h{i:02d}": c["id"] for i, c in enumerate(cases)}
    rev = {real: h for h, real in id_map.items()}
    gold = [{"case_id": rev[c["id"]], "diplotype": c["gt_diplotype"], "phenotype": c["gt_phenotype"],
             "recommendation": c["gt_drug"], "abstained": "false"} for c in cases]
    wrong = [dict(r, phenotype="Normal Metabolizer") for r in gold]
    _write(tmp_path / "cellA" / "rep1" / "claude_science_answers.csv", gold)
    _write(tmp_path / "cellA" / "rep2" / "claude_science_answers.csv", wrong)
    out = v.aggregate(tmp_path, cases, id_map)
    assert out["cellA"]["n_replicates"] == 2
    assert out["cellA"]["phenotype"]["per_replicate"][0] == 1.0
    assert out["cellA"]["phenotype"]["sd"] > 0


def test_shipped_runs_reproduce_the_reported_numbers():
    cases = v.load_cases()
    id_map = json.loads((HERE / "keys" / "id_map.json").read_text())
    out = v.aggregate(HERE / "runs", cases, id_map)
    assert out["api_corpus"]["n_replicates"] == 3
    assert out["corpus_only"]["n_replicates"] == 3
    assert round(out["api_corpus"]["phenotype"]["mean"], 3) == 0.942
    assert round(out["corpus_only"]["phenotype"]["mean"], 3) == 0.918
    assert out["api_corpus"]["abstained"] == [0, 5, 0]
    assert out["corpus_only"]["abstained"] == [2, 2, 2]
    assert out["api_corpus"]["lethal_errors"] == [6, 6, 6]
    assert out["corpus_only"]["lethal_errors"] == [5, 5, 0]


def test_leak_split_reproduces_the_identifier_finding():
    cases = v.load_cases()
    id_map = json.loads((HERE / "keys" / "id_map.json").read_text())
    res = v.leak_split(HERE / "pilot" / "claude_science_answers.csv",
                       HERE / "runs" / "api_corpus" / "rep1" / "claude_science_answers.csv",
                       cases, id_map)
    assert res["leaked"]["n"] == 104 and res["control"]["n"] == 6
    assert round(res["leaked"]["pilot"], 3) == 0.990
    assert round(res["leaked"]["hashed"], 3) == 0.962
    assert res["control"]["pilot"] == res["control"]["hashed"]


def test_expected_numbers_file_matches_a_fresh_recompute():
    expected = json.loads((HERE / "results" / "expected.json").read_text())
    assert v.check(expected) == []
