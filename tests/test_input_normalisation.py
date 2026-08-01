"""Tests for 74-input-normalisation.py (reviewer point R1.1).

The experiment asks a model to do the one job the paper argues it should do:
turn heterogeneous variant-level input into a canonical star-allele diplotype,
which validated code then executes. Everything the paper has measured so far is
the model doing the mapping we argue it should NOT do.

Two failure modes are worth more than the rest, and both are tested here.

1. THE ANSWER LEAKING INTO THE PROMPT. The extraction claim in an earlier
   revision was vacuous because both cells were handed the authored vocabulary.
   The same trap is available here: if the rendered genotype block, the
   instructions or the gene name carry the star-allele call, the model is
   copying, not normalising. `test_prompt_never_contains_the_answer` and
   `test_render_emits_no_star_alleles` exist to keep that impossible.

2. AN ABSTENTION SCORED AS A CALL. Coverage is the headline of this experiment,
   so a refusal that is silently parsed into a diplotype, or a blank response
   counted as a correct abstention, would move the number that matters most.
"""
import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent


def _load(name, path):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


norm = _load("input_normalisation", BASE / "code" / "74-input-normalisation.py")


VARIANTS = [
    {"chrom": "22", "pos": 42522613, "ref": "G", "alt": "C", "gt": "1|0"},
    {"chrom": "22", "pos": 42524947, "ref": "C", "alt": "T", "gt": "1|1"},
    {"chrom": "22", "pos": 42526694, "ref": "G", "alt": "A", "gt": "0|1"},
]


# --- rendering -------------------------------------------------------------

@pytest.mark.parametrize("form", ["vcf", "hgvs", "prose"])
def test_render_includes_every_variant(form):
    text = norm.render(VARIANTS, form)
    for v in VARIANTS:
        assert str(v["pos"]) in text, f"{form} rendering dropped position {v['pos']}"


@pytest.mark.parametrize("form", ["vcf", "hgvs", "prose"])
def test_render_emits_no_star_alleles(form):
    """The input must not contain the thing being asked for."""
    assert "*" not in norm.render(VARIANTS, form)


def test_renderings_are_genuinely_different():
    forms = {f: norm.render(VARIANTS, f) for f in ("vcf", "hgvs", "prose")}
    assert len(set(forms.values())) == 3


def test_hgvs_uses_genomic_notation():
    text = norm.render(VARIANTS, "hgvs")
    assert "g.42522613G>C" in text


def test_prose_claims_completeness_only_when_complete():
    """The prose form closes with "no other non-reference genotype was observed".
    That sentence is a claim about the input, and it is false whenever the
    variant list has been truncated. Asserting it anyway would put a fabricated
    statement into the prompt and then score the model on it."""
    assert "No other non-reference genotype" in norm.render(VARIANTS, "prose")
    truncated = norm.render(VARIANTS, "prose", truncated=True)
    assert "No other non-reference genotype" not in truncated
    assert "truncated" in truncated.lower()


@pytest.mark.parametrize("form", ["vcf", "hgvs", "prose"])
def test_truncation_is_disclosed_in_every_form(form):
    assert "truncat" in norm.render(VARIANTS, form, truncated=True).lower()


def test_render_rejects_unknown_form():
    with pytest.raises(ValueError):
        norm.render(VARIANTS, "telepathy")


def test_render_refuses_empty_variant_set():
    """A sample with no non-reference genotype is a real state, but it must be
    an explicit decision upstream, not an empty prompt that quietly asks a model
    to call a diplotype from nothing."""
    with pytest.raises(ValueError):
        norm.render([], "vcf")


# --- prompt construction ---------------------------------------------------

def test_prompt_never_contains_the_answer():
    p = norm.build_prompt("CYP2D6", VARIANTS, "vcf")
    assert "*4" not in p and "*1/*1" not in p
    assert "/" not in p.split("## Output")[0].replace("0|1", "").replace("1|0", "")[-200:]


def test_prompt_offers_abstention():
    p = norm.build_prompt("CYP2D6", VARIANTS, "vcf")
    assert "ABSTAIN" in p


def test_prompt_names_the_gene():
    assert "CYP2D6" in norm.build_prompt("CYP2D6", VARIANTS, "vcf")


def test_prompt_does_not_supply_a_candidate_diplotype_list():
    """Supplying the valid diplotype list would make this a multiple-choice
    question and reintroduce the vocabulary hand-off that made the extraction
    claim vacuous."""
    p = norm.build_prompt("CYP2D6", VARIANTS, "vcf")
    assert "valid diplotypes" not in p.lower()
    assert "choose from" not in p.lower()


# --- parsing ---------------------------------------------------------------

def test_load_models_covers_the_whole_common_panel(monkeypatch):
    """The five-cell comparison used eight models; this experiment must be able
    to use the same eight, or the input-step result speaks for a different panel
    than the mapping-step result it is compared against."""
    for var in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "GOOGLE_API_KEY",
                "DEEPSEEK_API_KEY"):
        monkeypatch.setenv(var, "test-key-not-used")
    panel = ["Claude Opus 4.5", "Claude Sonnet 4.5", "GPT-5.2", "GPT-4.1",
             "o3", "o4-mini", "Gemini 2.5 Flash", "DeepSeek V3"]
    clients = norm.load_models(panel)
    assert sorted(clients) == sorted(panel)
    assert all(callable(c) for c in clients.values())


def test_gemini_gets_an_output_budget_that_survives_its_reasoning(monkeypatch):
    """Gemini 2.5 Flash spends output tokens on reasoning before it emits text.
    At the default budget it returns empty text, which the parser would record as
    a refusal to answer. That is the C3 class of artefact: a harness limit read as
    model behaviour. Its budget must therefore be strictly larger than the shared
    default."""
    assert norm.GEMINI_OUT_TOKENS > 320


def test_parse_reads_a_diplotype():
    assert norm.parse_call("DIPLOTYPE: *4/*41") == "*4/*41"


def test_parse_tolerates_surrounding_text():
    assert norm.parse_call("Reasoning here.\nDIPLOTYPE: *1/*2\nDone.") == "*1/*2"


def test_parse_accepts_markdown_wrapped_marker():
    """C15. Presentation markup must not convert a complete marked call into
    an abstention."""
    assert norm.parse_call("**DIPLOTYPE: *1/*2**") == "*1/*2"


def test_parse_accepts_latex_boxed_marker():
    """The policy is model-neutral: a marked value inside a LaTeX box is the
    same schema value, not a Gemini-specific salvage rule."""
    text = r"The final answer is $\boxed{\text{DIPLOTYPE: *1/*9}}$"
    assert norm.parse_call(text) == "*1/*9"


def test_parse_preserves_trailing_explanation_policy():
    """Coverage records whether a marked answer was given; the first value
    token remains the answer even if the model violates the one-line request by
    adding an explanation after it."""
    assert norm.parse_call("DIPLOTYPE: *1/*2 followed by explanation") == "*1/*2"


def test_parse_rejects_unsupported_tandem_value_without_partial_salvage():
    """A marker-anywhere rule must validate the whole value, not silently take
    the first two star alleles from a longer unsupported expression."""
    assert norm.parse_call("DIPLOTYPE: *1/*80+*28") is None


def test_last_marked_value_is_the_final_answer():
    text = "DIPLOTYPE: *1/*2\nOn reflection, the data are insufficient.\nDIPLOTYPE: ABSTAIN"
    assert norm.parse_call(text) is None


def test_nonvalue_heading_does_not_hide_final_marked_call():
    text = "**Diplotype:**\n- haplotype one: *1\n- haplotype two: *2\nDIPLOTYPE: *1/*2"
    assert norm.parse_call(text) == "*1/*2"


def test_bold_heading_with_value_is_a_marked_call():
    assert norm.parse_call("**Full diplotype:** *1/*6") == "*1/*6"


def test_parse_strips_gene_prefix():
    assert norm.parse_call("DIPLOTYPE: CYP2D6 *4/*4") == "*4/*4"


def test_parse_strips_gene_prefix_without_a_separator():
    """C14. The published run required a space or hyphen after the gene symbol,
    so the standard PharmVar rendering with no separator was scored as an
    abstention. 744 well-formed diplotypes were lost to it, and because coverage
    is this experiment's headline the artefact moved the headline."""
    assert norm.parse_call("DIPLOTYPE: CYP2D6*4/*4", "CYP2D6") == "*4/*4"
    assert norm.parse_call("DIPLOTYPE: CYP2D6*4/*4") == "*4/*4"


def test_parse_strips_gene_prefix_case_insensitively():
    assert norm.parse_call("DIPLOTYPE: cyp2d6*4/*4", "CYP2D6") == "*4/*4"


def test_parse_strips_a_hyphenated_gene_prefix():
    assert norm.parse_call("DIPLOTYPE: CYP2D6-*4/*4", "CYP2D6") == "*4/*4"


def test_parse_strips_a_prefix_naming_another_gene():
    """The prefix is stripped because it is a gene symbol, not because it
    matches the gene asked about. Scoring against the reference decides
    correctness; the parser only decides whether an answer was given."""
    assert norm.parse_call("DIPLOTYPE: CYP2C19*1/*2", "CYP2D6") == "*1/*2"


def test_parse_does_not_strip_a_hedge_word():
    """The fix must not become C1 in reverse. A leading word that is not a gene
    symbol means the model qualified its answer, and a qualified answer is not a
    call."""
    assert norm.parse_call("DIPLOTYPE: unknown*1/*2", "CYP2D6") is None
    assert norm.parse_call("DIPLOTYPE: possibly*1/*2", "CYP2D6") is None
    assert norm.parse_call("DIPLOTYPE: indeterminate*1/*2") is None


def test_parse_returns_none_on_abstention():
    assert norm.parse_call("DIPLOTYPE: ABSTAIN") is None


def test_parse_returns_none_on_prose_refusal():
    """A refusal in prose must not be parsed into a call."""
    assert norm.parse_call("I cannot determine the diplotype from these data.") is None


def test_parse_returns_none_on_empty():
    assert norm.parse_call("") is None
    assert norm.parse_call(None) is None


def test_parse_does_not_invent_from_a_bare_star():
    assert norm.parse_call("DIPLOTYPE: *") is None


# --- abstention against truncation -----------------------------------------

def test_explicit_abstention_is_an_abstention():
    assert norm.classify("DIPLOTYPE: ABSTAIN", out_tokens=12, cap=2000) == "abstain"


def test_a_call_is_a_call():
    assert norm.classify("DIPLOTYPE: *4/*4", out_tokens=12, cap=2000) == "call"


def test_response_cut_off_at_the_cap_is_not_an_abstention():
    """A model reasoning past its output budget never reaches the DIPLOTYPE
    line. Reading that as "the model declined" publishes a harness cap as a
    finding, which is what a 320-token cap did on the first run of this
    experiment."""
    assert norm.classify("Let me work through the variants...",
                         out_tokens=2000, cap=2000) == "truncated_output"


def test_short_non_answer_is_an_abstention_not_a_truncation():
    assert norm.classify("I cannot say.", out_tokens=5, cap=2000) == "abstain"


def test_a_failed_call_is_never_an_abstention():
    """A 429, a timeout or any exhausted retry leaves empty text and zero
    tokens. Reading that as "the model declined to answer" is C3 in this
    project's corrections log, and it recurred here: a quota error on 94 o3
    calls was first recorded as 94 abstentions."""
    assert norm.classify("", out_tokens=0, cap=6000,
                         error="RateLimitError: 429") == "error"


def test_error_wins_even_if_text_looks_like_a_call():
    """A partial response banked alongside an error is not a result."""
    assert norm.classify("DIPLOTYPE: *1/*1", out_tokens=9, cap=6000,
                         error="APIError: 500") == "error"


# --- scoring ---------------------------------------------------------------

ROWS = [
    {"sample": "A", "gene": "CYP2D6", "call": "*1/*1", "reference": "*1/*1"},
    {"sample": "B", "gene": "CYP2D6", "call": "*4/*4", "reference": "*1/*4"},
    {"sample": "C", "gene": "CYP2D6", "call": None, "reference": "*1/*1"},
    {"sample": "D", "gene": "CYP2D6", "call": "*2/*1", "reference": "*1/*2"},
]


def test_coverage_counts_emitted_calls():
    s = norm.score(ROWS)
    assert s["n"] == 4
    assert s["emitted"] == 3
    assert s["coverage"] == pytest.approx(0.75)


def test_abstention_is_the_complement_of_coverage():
    s = norm.score(ROWS)
    assert s["abstention"] == pytest.approx(0.25)


def test_accuracy_is_among_emitted_not_among_all():
    """Dividing correct answers by the total rather than by what was emitted
    lets a model improve its score by abstaining, which is the exact error the
    paper accuses other benchmarks of."""
    s = norm.score(ROWS)
    assert s["accuracy_among_emitted"] == pytest.approx(2 / 3)


def test_diplotype_order_does_not_change_correctness():
    """*2/*1 and *1/*2 are the same diplotype."""
    s = norm.score([ROWS[3]])
    assert s["accuracy_among_emitted"] == pytest.approx(1.0)


def test_abstention_is_never_scored_as_correct():
    s = norm.score([{"sample": "C", "gene": "CYP2D6", "call": None, "reference": None}])
    assert s["emitted"] == 0
    assert s["accuracy_among_emitted"] is None


def test_score_refuses_rows_without_a_reference():
    with pytest.raises(ValueError):
        norm.score([{"sample": "X", "gene": "CYP2D6", "call": "*1/*1"}])


# --- definition-supplied arm -----------------------------------------------

DEFS = [
    {"StarAllele": "*1", "Function": "Normal Function", "GRCh37Core": "N/A"},
    {"StarAllele": "*3", "Function": "No Function", "GRCh37Core": "22-42524243-CT-C"},
    {"StarAllele": "*4", "Function": "No Function",
     "GRCh37Core": "22-42524947-C-T,22-42526694-G-A"},
]


def test_definition_block_lists_every_allele():
    b = norm.definition_block(DEFS)
    for d in DEFS:
        assert d["StarAllele"] in b


def test_definition_block_carries_defining_variants():
    assert "22-42524243-CT-C" in norm.definition_block(DEFS)


def test_definition_block_carries_function():
    assert "No Function" in norm.definition_block(DEFS)


def test_definition_prompt_contains_no_diplotype():
    """Supplying single-allele definitions is fair: PyPGx has them too. Supplying
    a DIPLOTYPE would hand over the answer and rebuild the tautology that voided
    the extraction claim."""
    p = norm.build_prompt("CYP2D6", VARIANTS, "vcf", definitions=DEFS)
    assert "*3/*3" not in p and "*1/*1" not in p and "*4/*4" not in p


def test_genotype_rendering_still_carries_no_star_allele():
    """The definitions block may name alleles. The rendered patient genotypes
    must still not, or the model is reading the answer off its own input."""
    assert "*" not in norm.render(VARIANTS, "vcf")


def test_definition_prompt_still_offers_abstention():
    assert "ABSTAIN" in norm.build_prompt("CYP2D6", VARIANTS, "vcf", definitions=DEFS)


def test_prompt_without_definitions_is_unchanged():
    assert norm.build_prompt("CYP2D6", VARIANTS, "vcf") == \
           norm.build_prompt("CYP2D6", VARIANTS, "vcf", definitions=None)


# --- stratified subsampling ------------------------------------------------

KEYS = {f"S{i}|{g}": {"gene": g}
        for g, n in (("CYP2D6", 30), ("TPMT", 20), ("NUDT15", 4))
        for i in range(n)}


def test_subsample_is_deterministic():
    a = norm.stratified_subsample(KEYS, 20)
    b = norm.stratified_subsample(KEYS, 20)
    assert a == b, "a subsample that changes between runs cannot be republished"


def test_subsample_returns_the_requested_size():
    assert len(norm.stratified_subsample(KEYS, 20)) == 20


def test_subsample_covers_every_gene():
    """A proportional draw must not drop the small genes entirely; NUDT15 has
    only 11 pairs in the real data and would vanish under naive sampling."""
    got = norm.stratified_subsample(KEYS, 20)
    assert {KEYS[k]["gene"] for k in got} == {"CYP2D6", "TPMT", "NUDT15"}


def test_subsample_larger_than_population_returns_all():
    assert len(norm.stratified_subsample(KEYS, 999)) == len(KEYS)


# --- input freezing --------------------------------------------------------

def test_inputs_are_hash_frozen():
    a = norm.freeze({"CYP2D6": VARIANTS})
    b = norm.freeze({"CYP2D6": list(VARIANTS)})
    assert a == b and len(a) == 64


def test_freeze_changes_when_a_genotype_changes():
    mutated = [dict(VARIANTS[0], gt="1|1")] + VARIANTS[1:]
    assert norm.freeze({"CYP2D6": VARIANTS}) != norm.freeze({"CYP2D6": mutated})
