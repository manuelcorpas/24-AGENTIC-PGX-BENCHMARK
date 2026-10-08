from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
import json
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def load(name):
    s=spec_from_file_location(name,ROOT/'code'/name);m=module_from_spec(s);s.loader.exec_module(m);return m

def test_supplementary_labels_describe_measured_endpoints():
    rows=json.loads((ROOT/'docs/R3-2026-09-16/DATA/attempted-scored-rows.json').read_text())
    m=load('86-figures-s2-s3-five-cell.py')
    high=[r for r in rows if r['lethal_action'] is not None]
    for fig in [m.figure_s2(high,save_output=False),m.figure_s3(high,save_output=False)]:
        labels=' '.join(t.get_text() for t in fig.findobj(plt.Text)).lower()
        assert not any(s in labels for s in ['lethal-class','information without action','correctness by coincidence','accuracy'])
        assert 'wording' in labels
        plt.close(fig)
    m=load('70-figures-revision.py')
    # Full attempted grid matters for the one missing execution response.
    fig=m.figure_r2({},rows=rows,already_scored=True,save_output=False)
    assert fig.axes[0].get_ylabel()=='Phenotype concordance with authored references'
    plt.close(fig)
    by={c:{model:sum(r['a1_phenotype'] for r in rows if r['cell']==c and r['model']==model)/330 for model in {r['model'] for r in rows}} for c in m.CELLS}
    fig=m.figure_r3(by,save_output=False)
    assert 'concordance' in fig.axes[0].get_ylabel()
    plt.close(fig)
    fig=m.figure_r4(save_output=False)
    assert 'authored' in fig.axes[0].get_ylabel()
    plt.close(fig)


def test_worked_example_shows_wording_sensitivity(tmp_path):
    m=load('93-figure-s9-worked-case.py')
    case,rows=m.collect(ROOT,m.CASE_ID,m.MODEL,m.REP)
    assert rows[1]['error'] and not rows[1]['widened_error']
    gene,rates=m.panel_error_rates(ROOT,m.CASE_ID)
    assert rates['rag_generation']==(45,48,0)
    fig=m.render(case,rows,gene,rates,tmp_path/'unused.tiff',m.MODEL,m.REP,save_output=False)
    texts=' '.join(t.get_text() for t in fig.findobj(plt.Text))
    assert 'LETHAL-CLASS ERROR' not in texts
    assert 'Original fail; widened pass' in texts
    assert not (tmp_path/'unused.tiff').exists()
    plt.close(fig)
    assert load('90-figure-s8-conditional-invariance.py').YLABEL=='Phenotype concordance, A1 (%)'

def test_locus_panels_fit_readable_page_width():
    m=load('86-figures-s2-s3-five-cell.py');fig=m.figure_s2(m.load(),save_output=False)
    assert fig.get_figwidth()<=7.2
    assert len(fig.axes)==6
    plt.close(fig)

def test_population_plot_uses_measured_label():
    m=load('90-figure-s8-conditional-invariance.py')
    fig=m.draw(m.compute(),save_output=False)
    assert fig.axes[0].get_ylabel()=='Phenotype concordance, A1 (%)'
    plt.close(fig)

def test_supplementary_title_and_gene_labels_do_not_overlap():
    import itertools
    m=load('86-figures-s2-s3-five-cell.py')
    fig=m.figure_s2(m.load(),save_output=False)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    title=fig._suptitle.get_window_extent(renderer)
    subtitle=next(t for t in fig.texts if t.get_text().startswith('Each point'))
    assert not title.overlaps(subtitle.get_window_extent(renderer))
    plt.close(fig)
    m=load('70-figures-revision.py');m.style()
    rows=json.loads((ROOT/'docs/R3-2026-09-16/DATA/attempted-scored-rows.json').read_text())
    caller=json.loads((ROOT/'docs/R3-2026-09-16/DATA/caller-comparison.json').read_text())
    fig=m.figure_r2({},rows=rows,already_scored=True,save_output=False,caller=caller)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    labels=[t for t in fig.axes[0].texts if isinstance(t,plt.Annotation)]
    for a,b in itertools.combinations(labels,2):
        assert not a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)),(a.get_text(),b.get_text())
    plt.close(fig)
    # A gene label must not cover a plotted gene, even if text pairs are separate.
    cyp=next(t for t in labels if t.get_text()=='CYP2C19')
    box=cyp.get_window_extent(renderer).expanded(1.01,1.1)
    assert not any(box.contains(*fig.axes[0].transData.transform(xy)) for xy in fig.axes[0].collections[0].get_offsets())

def test_locus_annotation_does_not_cover_other_locus():
 m=load('86-figures-s2-s3-five-cell.py');fig=m.figure_s2(m.load(),save_output=False)
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 for ax in fig.axes:
  for t in ax.texts:
   if not isinstance(t,plt.Annotation):continue
   box=t.get_window_extent(renderer)
   assert not any(box.contains(*ax.transData.transform(xy)) for c in ax.collections if hasattr(c,'get_offsets') for xy in c.get_offsets()),t.get_text()
 plt.close(fig)
