"""Two-rule Figure 2: data parity and every-pair text geometry at placement size."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]

def load():
 s=spec_from_file_location('figures_revision',ROOT/'code/70-figures-revision.py');m=module_from_spec(s);s.loader.exec_module(m);return m

def test_stacked_horizontal_panels_and_all_label_pairs():
 m=load();a=m._load('98-third-revision.py','r3');rows=m.json.loads(m.FIVE_CELL.read_text());summary=a.attempt_summary(rows,a.cases(),a.models(rows))
 fig=m.figure_r1(summary,save_output=False)
 try:
  top,bottom=fig.axes
  assert top.get_position().y0 > bottom.get_position().y0
  assert all(p.get_x()==0 for p in top.patches)
  assert [p.get_width() for p in top.patches]==pytest.approx([summary[c]['a2_recommendation_mean']*100 for c in m.CELLS])
  assert [p.get_width() for p in bottom.patches]==[46,137,22,24,23,43,75,22,24,23]
  assert any(p.get_hatch() for p in bottom.patches)
  assert {p.get_facecolor() for p in bottom.patches}.__len__()==2
  fig.canvas.draw();renderer=fig.canvas.get_renderer()
  labels=[t for ax in fig.axes for t in ax.texts]
  for i,left in enumerate(labels):
   for right in labels[i+1:]:assert not left.get_window_extent(renderer).overlaps(right.get_window_extent(renderer))
  assert all(t.get_fontsize()>=8 for t in labels)
  assert fig.get_figwidth()<=7.2
 finally:m.plt.close(fig)
