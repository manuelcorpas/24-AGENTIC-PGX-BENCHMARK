"""Layout gates for the manuscript's matched-factorial Figure 2."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "code" / "70-figures-revision.py"


def load_module():
    spec = spec_from_file_location("figures_revision", SCRIPT)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_figure_r1_preserves_vertical_bars_with_readable_external_labels():
    module = load_module()
    summary, _, _ = module.scored_cells()
    figure = module.figure_r1(summary, save_output=False)
    try:
        assert len(figure.axes) == 2
        assert len(figure.legends) == 1
        assert all(axis.get_legend() is None for axis in figure.axes)
        assert all(patch.get_y() == 0 for patch in figure.axes[0].patches)
        assert all(round(patch.get_width(), 2) == 0.68 for patch in figure.axes[0].patches)
        expected_accuracy = [summary[cell]["a2_recommendation_mean"] for cell in module.CELLS]
        assert [patch.get_height() for patch in figure.axes[0].patches] == pytest.approx(
            expected_accuracy
        )
        assert figure.axes[0].get_ylim()[1] >= 1.08
        assert [round(patch.get_height()) for patch in figure.axes[1].patches] == [
            46, 137, 22, 24, 23,
        ]
        expected = [
            "No\nknowledge",
            "RAG\ngeneration",
            "RAG call +\nauthored\nexecution",
            "Authored\nrules,\ngenerated",
            "Authored\nrules,\nexecuted",
        ]
        assert [label.get_text() for label in figure.axes[0].get_xticklabels()] == expected
        assert [label.get_text() for label in figure.axes[1].get_xticklabels()] == expected
        assert all(label.get_fontsize() >= 8.5 for label in figure.axes[0].get_xticklabels())
        assert all(label.get_fontsize() >= 8.5 for label in figure.axes[1].get_xticklabels())
        assert figure.get_figheight() <= 4.0
        assert figure.subplotpars.bottom <= 0.24
        assert all(text.get_position()[1] > patch.get_height()
                   for text, patch in zip(figure.axes[0].texts, figure.axes[0].patches))
    finally:
        module.plt.close(figure)
