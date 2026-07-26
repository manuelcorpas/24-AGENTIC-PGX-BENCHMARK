#!/usr/bin/env python3
"""
Tests for the matched factorial runner (revision item N1; R1.3, R2.1).

The reviewers' central objection: the original conditions differed in more than
the variable of interest. The skill arms were told the target drug and scored
drug-specifically; the free-prompted and retrieval arms were not. So the
accuracy gap confounds "executed rules help" with "an easier task".

"Matched" is not a claim to make in prose, it is a property to enforce. These
tests define it: across every generation cell the prompt must be byte-identical
apart from the knowledge block, must name the same target drug, and must demand
the same output schema. If someone later edits one prompt and not the others,
these fail.
"""
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
RUNNER = REPO / "code" / "60-matched-factorial.py"


@pytest.fixture(scope="module")
def mf():
    spec = spec_from_file_location("matched_factorial", RUNNER)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def case(mf):
    """A CYP2D6 codeine case: has a lethal tier, a real drug, real rules."""
    return next(c for c in mf.CASES if c["id"] == "cyp2d6_codeine_pm")


# ---------------------------------------------------------------- factorial shape

def test_the_five_realisable_cells_are_present(mf):
    assert set(mf.CELLS) == {
        "free_generation",
        "rag_generation",
        "rag_execution",
        "skill_generation",
        "skill_execution",
    }


def test_no_knowledge_execution_cell_is_absent(mf):
    """Execution is undefined without a structured input call. Do not fake the cell."""
    assert "free_execution" not in mf.CELLS
    assert "no_knowledge_execution" not in mf.CELLS


def test_cells_declare_their_factorial_coordinates(mf):
    """Each cell states its knowledge representation and decision mechanism."""
    for name, cell in mf.CELLS.items():
        assert cell["knowledge"] in {"none", "retrieved_prose", "structured_rules"}
        assert cell["mechanism"] in {"model_generation", "deterministic_execution"}
    assert mf.CELLS["rag_execution"]["knowledge"] == "retrieved_prose"
    assert mf.CELLS["rag_execution"]["mechanism"] == "deterministic_execution"


def test_reused_cells_are_marked_as_not_rerun(mf):
    """Skill arms are reused unchanged; the plan declares their schema canonical."""
    assert mf.CELLS["skill_generation"]["rerun"] is False
    assert mf.CELLS["skill_execution"]["rerun"] is False
    assert mf.CELLS["free_generation"]["rerun"] is True
    assert mf.CELLS["rag_generation"]["rerun"] is True
    assert mf.CELLS["rag_execution"]["rerun"] is True


# ---------------------------------------------------------------- the matched property

GENERATION_CELLS = ["free_generation", "rag_generation", "skill_generation"]


def test_every_prompt_names_the_target_drug(mf, case):
    """The confound, killed: no condition has to guess which drug is being asked about."""
    for name in mf.CELLS:
        prompt = mf.build_prompt(name, case)
        assert f"Drug: {case['drug']}" in prompt, f"{name} does not name the target drug"


def test_generation_cells_share_one_output_schema(mf, case):
    for name in GENERATION_CELLS:
        assert mf.SCHEMA_GENERATION in mf.build_prompt(name, case)


def test_execution_cells_share_one_output_schema(mf, case):
    for name in ["rag_execution", "skill_execution"]:
        assert mf.SCHEMA_EXECUTION in mf.build_prompt(name, case)


def test_generation_prompts_differ_only_in_the_knowledge_block(mf, case):
    """The heart of the matched design, enforced mechanically.

    Strip each prompt's knowledge block and the remainder must be identical
    across cells: same patient text, same drug, same schema, same instruction.
    """
    remainders = {}
    for name in GENERATION_CELLS:
        prompt = mf.build_prompt(name, case)
        remainders[name] = prompt.replace(mf.knowledge_block(name, case), "<KNOWLEDGE>")
    distinct = set(remainders.values())
    assert len(distinct) == 1, (
        "generation prompts differ outside the knowledge block; that is a confound: "
        + repr({k: v[:200] for k, v in remainders.items()})
    )


def test_knowledge_blocks_actually_differ(mf, case):
    """Guard against the previous test passing because everything is empty."""
    blocks = {mf.knowledge_block(n, case) for n in GENERATION_CELLS}
    assert len(blocks) == 3


def test_single_framing_only(mf, case):
    """Population framings are identical genotypes with different wording, so the
    matched grid runs one framing and invariance is checked separately (R2.1)."""
    prompt = mf.build_prompt("free_generation", case)
    for wording in ("Peruvian", "Ugandan", "European", "cohort"):
        assert wording not in prompt


# ---------------------------------------------------------------- no ground-truth leakage

def test_free_cell_receives_nothing_at_all(mf, case):
    """The no-knowledge cell must contain no part of the answer, in any form."""
    prompt = mf.build_prompt("free_generation", case)
    assert case["gt_diplotype"] not in prompt
    assert case["gt_phenotype"] not in prompt
    assert case["gt_drug"] not in prompt


@pytest.mark.parametrize("cell", ["free_generation", "rag_generation", "rag_execution",
                                  "skill_generation", "skill_execution"])
def test_knowledge_never_keys_the_answer_to_the_patient_genotype(mf, case, cell):
    """The real leakage test, and the reason the naive one is wrong.

    The retrieved-prose and structured-rules cells legitimately contain
    diplotypes, phenotypes and recommendations: a CPIC guideline excerpt and a
    validated rule table are exactly what those conditions are DEFINED as
    providing. Asserting their absence would assert the conditions away.

    What must never appear is the bridge from THIS patient to the answer. Both
    knowledge representations are keyed by diplotype (*4/*4); the patient is
    described only by a raw genotype (rs3892097 T/T). The model still has to
    perform the genotype-to-diplotype call, which is the step the paper claims
    all residual error localises to. If the patient's genotype string ever
    appeared inside the knowledge block, that step could be skipped by lookup
    and every downstream claim about error localisation would be void.
    """
    block = mf.knowledge_block(cell, case)
    assert case["genotype"] not in block


# ---------------------------------------------------------------- deterministic execution

def test_execute_skill_is_deterministic_and_correct(mf, case):
    phen, rec = mf.execute_skill(case["gene"], case["drug"], case["gt_diplotype"])
    assert phen == case["gt_phenotype"]
    assert rec == case["gt_drug"]
    assert (phen, rec) == mf.execute_skill(case["gene"], case["drug"], case["gt_diplotype"])


def test_execute_skill_abstains_on_unknown_diplotype(mf, case):
    """Out-of-vocabulary input must abstain, not guess. Same contract as N0."""
    assert mf.execute_skill(case["gene"], case["drug"], "*99/*99") == (None, None)
    assert mf.execute_skill(case["gene"], case["drug"], "") == (None, None)


# ---------------------------------------------------------------- cost safety

def test_dry_run_makes_no_api_calls(mf, monkeypatch):
    """A dry run must be free. Guards a US$500-1,100 mistake."""
    called = []
    monkeypatch.setattr(mf, "MODELS", {"fake": lambda p: called.append(p) or "DIPLOTYPE: x"})
    mf.main(["--dry-run", "--limit", "3", "--models", "fake"])
    assert called == []


def test_planned_call_count_matches_the_registered_budget(mf):
    """The number quoted to the editor is computed here, not typed by hand."""
    n = mf.planned_calls(n_models=9, n_reps=3)
    assert n == 3 * 110 * 9 * 3          # 3 rerun cells x 110 cases x 9 models x 3 reps
    assert n == 8910


def test_norm_dip_still_matches_the_original_runner():
    """The shared module copied norm_dip out of 41-armA9-armBv2.py so the runner
    could be imported without API keys. If either copy drifts, the reused
    skill-arm data and the new matched cells would normalise diplotypes
    differently and the comparison would silently stop being matched.
    """
    import inspect
    import re as _re

    spec = spec_from_file_location("_pgx_rules", REPO / "code" / "_pgx_rules.py")
    shared = module_from_spec(spec)
    spec.loader.exec_module(shared)

    original_src = (REPO / "code" / "41-armA9-armBv2.py").read_text()
    match = _re.search(r"def norm_dip\(s, gene\):\n(?:.*\n)*?(?=\ncases = )", original_src)
    assert match, "could not locate norm_dip in 41-armA9-armBv2.py"

    def body(text):
        return [ln.rstrip() for ln in text.strip().splitlines() if ln.strip()]

    assert body(match.group(0)) == body(inspect.getsource(shared.norm_dip))
