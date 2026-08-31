"""Regression tests for the model-level real-genome Figure 3 data."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import json
import numpy as np


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "code" / "82-figure7-real-genome.py"
DATA = REPO / "data" / "v4_flash_agent_vs_deterministic.json"


def load_module():
    spec = spec_from_file_location("figure7_real_genome", SCRIPT)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_model_level_matrices_preserve_models_and_four_outcome_types():
    module = load_module()
    data = json.loads(DATA.read_text())

    labels, correct, abstained, wrong, invalid = module.model_level_matrices(data)

    assert labels == module.MODELS
    assert correct.shape == abstained.shape == wrong.shape == invalid.shape == (8, 5)
    total = correct + abstained + wrong + invalid
    np.testing.assert_allclose(total[~np.isnan(total)], 100.0)


def test_matched_iberian_cohort_is_adjacent_to_corpas_family_in_both_panels():
    module = load_module()
    data = json.loads(DATA.read_text())

    assert module.COLUMN_KEYS == [
        "benchmark",
        "CorpasFamily",
        "1000G_IBS",
        "Peru",
        "UGR",
    ]
    assert module.COLUMN_LABELS == [
        "Curated\ncases",
        "Family",
        "IBS",
        "Peru",
        "Uganda",
    ]

    correct, _, _, _ = module.pooled_percentages(data)
    np.testing.assert_allclose(
        correct,
        [96.5530303030, 63.9423076923, 61.1842105263, 55.3571428571, 38.3552631579],
    )

    fig = module.make_figure(data)
    try:
        expected = ["Curated\ncases", "Family", "IBS", "Peru", "Uganda"]
        for axis in fig.axes:
            assert [label.get_text() for label in axis.get_xticklabels()] == expected
    finally:
        module.plt.close(fig)


def test_model_level_matrices_do_not_conflate_invalid_and_wrong():
    module = load_module()
    data = json.loads(DATA.read_text())

    labels, correct, abstained, wrong, invalid = module.model_level_matrices(data)
    col = module.COLUMN_KEYS.index("UGR")

    opus = labels.index("Claude Opus 4.5")
    gemini = labels.index("Gemini 2.5 Flash")
    np.testing.assert_allclose(
        [correct[opus, col], abstained[opus, col], wrong[opus, col], invalid[opus, col]],
        [44.2105263158, 39.4736842105, 16.3157894737, 0.0],
    )
    np.testing.assert_allclose(
        [correct[gemini, col], abstained[gemini, col], wrong[gemini, col], invalid[gemini, col]],
        [34.2105263158, 14.7368421053, 50.5263157895, 0.5263157895],
    )
    np.testing.assert_allclose(
        correct[:, col] + abstained[:, col] + wrong[:, col] + invalid[:, col],
        100.0,
    )


def test_pooled_panel_reproduces_the_equal_weight_summary():
    module = load_module()
    data = json.loads(DATA.read_text())

    correct, abstained, wrong, invalid = module.pooled_percentages(data)

    assert correct.shape == abstained.shape == wrong.shape == invalid.shape == (5,)
    np.testing.assert_allclose(correct + abstained + wrong + invalid, 100.0)
    np.testing.assert_allclose(
        [correct[-1], abstained[-1], wrong[-1], invalid[-1]],
        [38.3552631579, 39.6710526316, 21.8421052632, 0.1315789474],
    )


def test_curated_benchmark_is_computed_separately_for_each_model():
    module = load_module()
    by_model = module.benchmark_by_model()

    assert set(by_model) == set(module.MODELS)
    assert all(by_model[model][1] == 0.0 for model in module.MODELS)
    assert all(len(values) == 4 for values in by_model.values())
    assert not any(np.isnan(value) for value in by_model["DeepSeek V4 Flash"])


def test_deepseek_v4_flash_is_version_matched_and_included_in_pooling():
    module = load_module()
    data = json.loads(DATA.read_text())

    assert module.MODELS[-1] == "DeepSeek V4 Flash"
    assert module.REAL_MODEL_KEYS["DeepSeek V4 Flash"] == "DeepSeek V4 Flash"
    assert module.BENCHMARK_MODEL_KEYS["DeepSeek V4 Flash"] == "DeepSeek V4 Flash"
    assert module.POOLED_MODELS == module.MODELS
    assert data["model_provenance"]["DeepSeek V4 Flash"]["returned_model_id"] == [
        "deepseek-v4-flash"
    ]


def test_figure_has_model_level_panel_and_pooled_bar_panel():
    module = load_module()
    data = json.loads(DATA.read_text())

    fig = module.make_figure(data)
    try:
        assert len(fig.axes) == 5
        assert [axis.get_title() for axis in fig.axes[:4]] == [
            "Correct",
            "Abstained",
            "Clinically wrong",
            "Invalid output",
        ]
        assert fig.axes[4].get_title(loc="left") == "B  Equal-weight pooled summary"
        assert len(fig.axes[4].patches) == 20
        assert any(text.get_text() == "A" for text in fig.texts)
    finally:
        module.plt.close(fig)
