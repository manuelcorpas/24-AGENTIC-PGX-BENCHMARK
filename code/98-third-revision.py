#!/usr/bin/env python3
"""R3 derived analyses. Frozen scorers, references and responses are never edited.

Run: python code/98-third-revision.py
The zero-credit RYR1 calculation is arithmetic sensitivity, not a new scorer.
"""
from __future__ import annotations
import argparse
from collections import defaultdict, Counter
from copy import deepcopy
from functools import lru_cache
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
import itertools
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
CELLS = ('free_generation','rag_generation','rag_execution','skill_generation','skill_execution')
AUTHORED = CELLS[2:]
RYR1 = 'ryr1_anaes_norm'
DATA = ROOT / 'data/v4_flash_eight_model_five_cell.json'
OUT = ROOT / 'docs/R3-2026-09-16/DATA'


@lru_cache(None)
def module(filename):
    s = spec_from_file_location(filename.replace('-','_'),ROOT/'code'/filename)
    m = module_from_spec(s); s.loader.exec_module(m)
    return m


def scorer(): return module('61-rescore-matched.py')
def cases(): return json.loads((ROOT/'specs/test_cases_v3.json').read_text())
def models(rows): return sorted({r['model'] for r in rows})
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def complete_rows(rows, references, panel, reps=(0,1,2)):
    if not rows: raise ValueError('empty input')
    if not references or not panel or not reps: raise ValueError('empty design')
    by_id = {c['id']:c for c in references}
    if len(by_id)!=len(references): raise ValueError('duplicate case')
    expected = set(itertools.product(CELLS,by_id,panel,reps))
    indexed={}
    for r in rows:
        key=tuple(r[k] for k in ('cell','case_id','model','rep'))
        if key in indexed: raise ValueError(f'duplicate evaluation: {key}')
        if key not in expected: raise ValueError(f'unexpected evaluation: {key}')
        indexed[key]=r
    out=list(rows)
    for cell,cid,model,rep in sorted(expected-indexed.keys()):
        c=by_id[cid]
        out.append(dict(cell=cell,case_id=cid,model=model,rep=rep,gene=c['gene'],drug=c['drug'],raw='',error='missing attempt',missing_attempt=True))
    return out


def score_attempts(rows, references, panel, equivalence=False):
    rs=scorer(); by_id={c['id']:c for c in references}
    fn=None
    if equivalence:
        eq=rs.load_equivalence_scorer()
        fn=eq.score_a1_clinical_eq
    out=[]
    for r in complete_rows(rows,references,panel):
        s=rs.score_row(r,by_id[r['case_id']],fn)
        s['missing_attempt']=bool(r.get('missing_attempt'))
        # Wording sensitivity is separate from the frozen phenotype equivalence layer.
        if s['lethal_action'] is not None:
            widened=bool(re.search(r'\bavoid\b|\bcontraindicat|\bdo not use\b',r.get('parsed_drug') or '',re.I)) and not s['substituted']
            s['widened_action']=float(widened)
        if not s['parsed_ok'] or s['abstained']:
            s['a1_phenotype']=s['a2_recommendation']=0.0
        out.append(s)
    return out


def attempt_summary(rows,references,panel,equivalence=False):
    scored=score_attempts(rows,references,panel,equivalence)
    out={}
    for cell in CELLS:
        ss=[s for s in scored if s['cell']==cell]; n=len(ss)
        hs=[s for s in ss if s['lethal_action'] is not None]
        out[cell]={
            'n_attempted':n,'n_returned':sum(not s['missing_attempt'] for s in ss),
            'missing_attempts':sum(s['missing_attempt'] for s in ss),
            'parse_failures':sum(not s['parsed_ok'] for s in ss),
            'abstentions':sum(s['abstained'] for s in ss),
            'coverage':sum(s['parsed_ok'] and not s['abstained'] for s in ss)/n,
            'a1_phenotype_mean':sum(s['a1_phenotype'] for s in ss)/n,
            'a2_recommendation_mean':sum(s['a2_recommendation'] for s in ss)/n,
            'high_severity_n':len(hs),
            'strict_failures':sum(s['lethal_action']==0 for s in hs),
            'widened_failures':sum(s['widened_action']==0 for s in hs),
            'drug_substitutions':sum(s['substituted'] for s in ss),
        }
    return out


def wording_intervals(scored):
    bs=module('71-lethal-case-level-stats.py');out={}
    for rule,key in [('strict','lethal_action'),('widened','widened_action')]:
        by=defaultdict(lambda:defaultdict(list)); genes={}
        for r in scored:
            if r['lethal_action'] is None:continue
            by[r['case_id']][r['cell']].append(1-r[key]);genes[r['case_id']]=r['gene']
        out[rule]={}
        for name,ids in [('all',sorted(by)),('HLA',sorted(c for c in by if genes[c].startswith('HLA'))),('non-HLA',sorted(c for c in by if not genes[c].startswith('HLA')))]:
            out[rule][name]={'cases':ids,'retrieval_minus_free':bs.bootstrap(by,ids,'free_generation','rag_generation')}
    return out


def screen(rows,scored,eq_scored):
    result=[]
    for c in cases():
        by={cell:[s for s in scored if s['case_id']==c['id'] and s['cell']==cell] for cell in CELLS}
        means={cell:sum(s['a1_phenotype'] for s in ss)/len(ss) for cell,ss in by.items()}
        if not (all(means[cell]>=0.95 for cell in AUTHORED) and means['rag_generation']<=0.10):continue
        rr=[r for r in rows if r['case_id']==c['id']]
        eq=[s for s in eq_scored if s['case_id']==c['id'] and s['cell']=='rag_generation']
        result.append({'case_id':c['id'],'reference':c,'phenotype_means':means,
            'rag_equivalence_sum':sum(s['a1_phenotype'] for s in eq),'rag_n':len(eq),
            'rag_recommendation_mean':sum(s['a2_recommendation'] for s in by['rag_generation'])/len(by['rag_generation']),
            'observed_phenotypes':{cell:dict(Counter(r.get('parsed_phenotype','') for r in rr if r['cell']==cell)) for cell in CELLS},
            'responses':[{k:r.get(k) for k in ('cell','model','rep','raw','parsed_phenotype','parsed_drug')} for r in rr]})
    return result


def extraction_table(rows):
    ev=module('68-evaluate-extracted-rules.py');rs=scorer();cc=cases();by={c['id']:c for c in cc}
    _,_,rd,_=rs.rules.build_rules(cc)
    calls=[dict(gene=r['gene'],called_diplotype=r.get('called_diplotype') or '',gt_phenotype=by[r['case_id']]['gt_phenotype']) for r in rows if r['cell']=='skill_execution']
    extracted=json.loads((ROOT/'data/v4_flash_eight_model_extracted_rules.json').read_text()); tables=defaultdict(dict)
    for r in extracted:
        if not r.get('error'):tables[r['model']][r['gene']]=r['table']
    out={'authored_execution':ev.score_rows(calls,{g:{'diplotype_to_phenotype':v} for g,v in rd.items()}),'models':{}}
    for model,tt in sorted(tables.items()):
        diffs=[ev.diff_gene(g,rd.get(g,{}),t) for g,t in sorted(tt.items())]
        out['models'][model]={'execution':ev.score_rows(calls,tt),'shared':sum(d['shared'] for d in diffs),'agree':sum(d['agree'] for d in diffs),'missing':sum(len(d['missing']) for d in diffs),'extra':sum(len(d['extra']) for d in diffs)}
    return out


def corrected_references(original):
    corrected=deepcopy(original)
    for c in corrected:
        if c['id']==RYR1:
            c['gt_phenotype']='Uncertain susceptibility'
            c['gt_drug']='volatile-anaesthetics: a negative RYR1 result does not establish normal susceptibility; assess clinical findings, family history and other laboratory data to guide anaesthetic choice and further evaluation'
    return corrected


def separate_inputs(original):
    runner=module('60-matched-factorial.py'); inputs=[]; refs=[];mapping=[]; hashes=[]
    for i,c in enumerate(original,1):
        opaque=f'PGX-{i:04d}'
        inp={k:c[k] for k in ('gene','drug','genotype')};inp['id']=opaque;inputs.append(inp)
        refs.append(dict(id=opaque,**{k:v for k,v in c.items() if k.startswith('gt_')}))
        mapping.append(dict(opaque_id=opaque,historical_id=c['id']))
        for cell in CELLS:
            old=runner.build_prompt(cell,c);new=runner.build_prompt(cell,inp)
            if old!=new:raise ValueError(f'prompt changed for {c["id"]} {cell}')
            hashes.append(dict(id=opaque,cell=cell,sha256=hashlib.sha256(old.encode()).hexdigest()))
    return inputs,refs,mapping,dict(compared_prompts=len(hashes),identical=True,prompts=hashes)


def analyse(rows):
    cc=cases();panel=models(rows);rs=scorer()
    frozen=rs.score_both(rows)
    sc=score_attempts(rows,cc,panel);eq=score_attempts(rows,cc,panel,True)
    baseline=attempt_summary(rows,cc,panel)
    excluded=[r for r in rows if r['case_id']!=RYR1]; excases=[c for c in cc if c['id']!=RYR1]
    out={'frozen_returned_record_results':frozen,'baseline':baseline,'equivalence':attempt_summary(rows,cc,panel,True),
        'excluded_ryr1':{'baseline':attempt_summary(excluded,excases,panel),'equivalence':attempt_summary(excluded,excases,panel,True)},
        'wording_intervals':wording_intervals(sc),'candidates':screen(rows,sc,eq),'table2':extraction_table(rows),
        'recommendation_wording_disagreement':{cell:sum(s['lethal_action']==0 and s['a2_recommendation']==1 for s in sc if s['cell']==cell) for cell in CELLS}}
    se=[s for s in sc if s['cell']=='skill_execution'];zero=[s for s in se if s['case_id']==RYR1]
    out['ryr1_zero_credit']={'definition':'Arithmetic sensitivity only; not corrected-reference output of the frozen scorer','n_zeroed':len(zero),'denominator':len(se),'original_sum':sum(s['a2_recommendation'] for s in se),'removed_credit':sum(s['a2_recommendation'] for s in zero),'recommendation_mean':sum(s['a2_recommendation'] for s in se if s['case_id']!=RYR1)/len(se)}
    out['exclusion_max_change_pp']=max(abs(baseline[c][metric]-out['excluded_ryr1']['baseline'][c][metric])*100 for c in CELLS for metric in ['a1_phenotype_mean','a2_recommendation_mean'])
    out['exclusion_authored_max_change_pp']=max(abs(baseline[c][metric]-out['excluded_ryr1']['baseline'][c][metric])*100 for c in AUTHORED for metric in ['a1_phenotype_mean','a2_recommendation_mean'])
    return out


def statistics_rows(rows):
    return [dict(s,framing='single') for s in score_attempts(rows,cases(),models(rows))]


def caller_comparison(rows):
    path=ROOT/'real-genome-arm/scripts/08_caller_truth_eval.py'
    s=spec_from_file_location('caller_truth',path);ev=module_from_spec(s);s.loader.exec_module(ev)
    by={c['id']:c for c in cases()};truth=[];calls=[]
    for r in rows:
        if r['cell'].endswith('execution') and r.get('called_diplotype') and not r.get('abstained'):
            key=dict(sample='|'.join(str(r[k]) for k in ('case_id','model','cell','rep')),gene=r['gene'])
            truth.append(dict(key,diplotype=by[r['case_id']]['gt_diplotype']))
            calls.append(dict(key,diplotype=r['called_diplotype']))
    return ev.evaluate(truth,calls)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=OUT);args=ap.parse_args()
    rows=json.loads(DATA.read_text());result=analyse(rows)
    inp,refs,mapping,proof=separate_inputs(cases());corrected=corrected_references(cases())
    corrected_opaque=[dict(id=m['opaque_id'],**{k:v for k,v in c.items() if k.startswith('gt_')}) for c,m in zip(corrected,mapping)]
    source_paths=[DATA,ROOT/'data/v3_five_cell_live.json',ROOT/'data/v4_flash_eight_model_extracted_rules.json',ROOT/'specs/test_cases_v3.json',ROOT/'specs/cpic_rag_corpus_v3.json',ROOT/'SCORING-PREREG.md',ROOT/'real-genome-arm/scripts/08_caller_truth_eval.py',ROOT/'real-genome-arm/scripts/07_executed_pipeline_n0.py']+[ROOT/'code'/f for f in ['10-rescore-v3.py','10b-rescore-v3-clinical-equivalence.py','60-matched-factorial.py','61-rescore-matched.py','65-hierarchical-stats.py','68-evaluate-extracted-rules.py','70-figures-revision.py','71-lethal-case-level-stats.py','73-crossed-mixed-model.py','86-figures-s2-s3-five-cell.py','_pgx_rules.py','98-third-revision.py']]
    args.out.mkdir(parents=True,exist_ok=True)
    outputs={'analysis.json':result,'attempted-scored-rows.json':statistics_rows(rows),'model-inputs-v1.json':inp,'evaluator-references-original.json':refs,'evaluator-references-corrected-v1.json':corrected_opaque,'identifier-map.json':mapping,'prompt-identity-proof.json':proof,'test-cases-corrected-v1.json':corrected,'source-hashes.json':{str(p.relative_to(ROOT)):digest(p) for p in source_paths}}
    outputs['caller-comparison.json']=caller_comparison(rows)
    for name,value in outputs.items():(args.out/name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:result[k] for k in ['ryr1_zero_credit','exclusion_max_change_pp','exclusion_authored_max_change_pp','recommendation_wording_disagreement']},indent=2))


if __name__=='__main__':main()
