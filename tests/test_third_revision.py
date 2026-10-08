"""Acceptance gates for the separately derived third-revision analyses."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import json
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load():
    spec = spec_from_file_location('revision3', ROOT / 'code/98-third-revision.py')
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope='module')
def analysis():
    m = load()
    return m, m.analyse(json.loads((ROOT / 'data/v4_flash_eight_model_five_cell.json').read_text()))


def test_attempts_and_exact_sensitivity(analysis):
    m, a = analysis
    s = a['baseline']['skill_execution']
    assert s['n_attempted'] == 2640 and s['n_returned'] == 2639
    assert s['missing_attempts'] == 1
    assert s['coverage'] == pytest.approx(2638 / 2640)
    assert s['a2_recommendation_mean'] == pytest.approx(2566.5 / 2640)
    assert a['ryr1_zero_credit']['n_zeroed'] == 24
    assert a['ryr1_zero_credit']['recommendation_mean'] == pytest.approx(2542.5 / 2640)
    assert round(a['ryr1_zero_credit']['recommendation_mean'] * 100, 2) == 96.31
    assert a['excluded_ryr1']['baseline']['skill_execution']['n_attempted'] == 2616
    assert round(a['excluded_ryr1']['baseline']['skill_execution']['a2_recommendation_mean'] * 100, 2) == 97.19


def test_both_panels_and_intervals(analysis):
    m, a = analysis
    for dataset, expected in [('v4_flash_eight_model_five_cell.json', [46,137,43,75]), ('v3_five_cell_live.json',[54,137,50,75])]:
        rows = json.loads((ROOT / 'data' / dataset).read_text())
        s = m.attempt_summary(rows, m.cases(), m.models(rows))
        assert [s[c][k] for k in ['strict_failures','widened_failures'] for c in ['free_generation','rag_generation']] == expected
        assert all(v['widened_failures'] <= v['strict_failures'] for v in s.values())
        assert all(v['high_severity_n'] == 336 for v in s.values())
    d = a['wording_intervals']['widened']['all']['retrieval_minus_free']
    assert [round(d[k]*100,1) for k in ['difference','ci_low','ci_high']] == [9.5,-2.7,25.3]
    assert a['recommendation_wording_disagreement'] == {'free_generation':15,'rag_generation':96,'rag_execution':0,'skill_generation':0,'skill_execution':0}
    m.scorer().load_equivalence_scorer().verify_frozen()


def test_screen_and_equivalence(analysis):
    m,a=analysis
    candidates={c['case_id']:c for c in a['candidates']}
    assert len(candidates)==8
    for cid, numerator in {'hlab5701_aba_pos':24,'hlaa3101_carba_neg':24,'hlab5801_allo_pos':23,'hlab5801_allo_neg':22,'hlab5701_aba_neg':19,'mtrnr1_amino_var':0}.items():
        assert candidates[cid]['rag_equivalence_sum']==numerator
    assert candidates['ryr1_anaes_norm']['phenotype_means']['rag_generation']==0
    assert candidates['cyp2d6_amitrip_im']['phenotype_means']['skill_execution']==1
    assert round(a['equivalence']['rag_generation']['a1_phenotype_mean']*100,1)==88.7


def test_invalid_inputs_and_missing_accounting():
    m=load(); rows=json.loads((ROOT/'data/v4_flash_eight_model_five_cell.json').read_text())
    with pytest.raises(ValueError,match='empty'):m.attempt_summary([],m.cases(),m.models(rows))
    with pytest.raises(ValueError,match='duplicate'):m.attempt_summary(rows+[rows[0]],m.cases(),m.models(rows))
    s=m.attempt_summary(rows[1:],m.cases(),m.models(rows))['free_generation']
    assert s['n_attempted']==2640 and s['missing_attempts']==1
    assert s['a1_phenotype_mean'] < m.attempt_summary(rows,m.cases(),m.models(rows))['free_generation']['a1_phenotype_mean']


def test_statistics_export_retains_full_grid():
    m=load();rows=json.loads((ROOT/'data/v4_flash_eight_model_five_cell.json').read_text())
    scored=m.statistics_rows(rows)
    assert len(scored)==13200
    assert {r['framing'] for r in scored}=={'single'}
    missing=[r for r in scored if r['missing_attempt']]
    assert len(missing)==1 and missing[0]['a1_phenotype']==0


def test_caller_comparison_uses_v4_calls():
    m=load();rows=json.loads((ROOT/'data/v4_flash_eight_model_five_cell.json').read_text())
    out=m.caller_comparison(rows)
    assert out['evaluable']==5272 and out['call_errors']==225
    assert out['call_concordance']==0.9573
    assert not any('DeepSeek V3' in e['sample'] for v in out['per_gene'].values() for e in v['errors'])


def test_table2_and_opaque_inputs(analysis):
    m,a=analysis
    t=a['table2']['authored_execution']
    assert t['n']==2639 and t['coverage']==0.9996
    assert round(t['accuracy_among_emitted']*100,1)==96.7
    assert 'DeepSeek V4 Flash' in a['table2']['models']
    inputs,refs,mapping,proof=m.separate_inputs(m.cases())
    assert len(inputs)==len(refs)==len(mapping)==110
    assert all(set(x)=={'id','gene','drug','genotype'} for x in inputs)
    assert len({x['id'] for x in inputs})==110
    assert proof['compared_prompts']==550 and proof['identical']
    corrected=m.corrected_references(m.cases())
    assert [c['id'] for c,o in zip(corrected,m.cases()) if c!=o]==['ryr1_anaes_norm']
    assert next(c for c in corrected if c['id']=='ryr1_anaes_norm')['gt_phenotype']=='Uncertain susceptibility'
