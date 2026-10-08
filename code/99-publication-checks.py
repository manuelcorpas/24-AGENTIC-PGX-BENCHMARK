#!/usr/bin/env python3
"""Expose screen boundaries, frozen RYR1 fallback and all registered raw rows.

Run from any directory. No scorer, reference or historical response is changed.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from copy import deepcopy
from importlib.util import module_from_spec, spec_from_file_location
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'docs/R3-2026-09-17/DATA'

def load(name):
    s=spec_from_file_location(name, ROOT/'code'/name)
    m=module_from_spec(s);s.loader.exec_module(m);return m

def analysis():
    m=load('98-third-revision.py'); rows=json.loads(m.DATA.read_text())
    scored=m.score_attempts(rows,m.cases(),m.models(rows))
    by=defaultdict(lambda:defaultdict(list))
    for r in scored:by[r['case_id']][r['cell']].append(r['a1_phenotype'])
    partition=dict(selected=0,authored_only=0,rag_only=0,neither=0)
    all_cases=[];near=[]
    for cid,cs in sorted(by.items()):
        means={c:sum(v)/len(v) for c,v in cs.items()}
        authored=min(means[c] for c in m.AUTHORED);rag=means['rag_generation']
        a=authored>=.95;r=rag<=.10
        group='selected' if a and r else 'authored_only' if a else 'rag_only' if r else 'neither'
        partition[group]+=1
        if group!='selected' and authored>=.90 and rag<=.15:near.append(cid)
        all_cases.append(dict(case_id=cid,phenotype_means=means,authored_min=authored,rag_mean=rag,group=group))
    corrected=next(c for c in m.corrected_references(m.cases()) if c['id']==m.RYR1)
    baseline=m.scorer().baseline;rr=[]
    for r in rows:
        if r['case_id']!=m.RYR1:continue
        s=m.scorer().score_row(r,corrected)
        rr.append(dict(cell=r['cell'],model=r['model'],rep=r['rep'],parsed_drug=r.get('parsed_drug'),a1_corrected=s['a1_phenotype'],a2_corrected=s['a2_recommendation']))
    remaining=sum(s['a2_recommendation'] for s in scored if s['cell']=='skill_execution' and s['case_id']!=m.RYR1)
    denominator=sum(s['cell']=='skill_execution' for s in scored)
    diagnostic=dict(corrected_reference=corrected,n_responses=len(rr),responses=rr,
        correct_reference_against_itself_a1=baseline.score_a1(corrected['gt_phenotype'],corrected['gt_phenotype']),
        correct_reference_against_itself_a2=baseline.score_a2(corrected['gt_drug'],corrected['gt_drug']),
        authored_execution_full_grid_fallback_mean=(remaining+sum(r['a2_corrected'] for r in rr if r['cell']=='skill_execution'))/denominator,
        zero_credit_mean=remaining/denominator,
        interpretation='The corrected recommendation matches no reference-action branch. Every nonempty RYR1 recommendation reaches the generic 0.5 fallback, including the corrected reference scored against itself. The resulting 96.76% is a software diagnostic, not clinically meaningful corrected-reference concordance. Zero-credit arithmetic sensitivity remains 96.31%.')
    return dict(screen=dict(partition=partition,outside_count=110-partition['selected'],near_miss_margin=.05,near_miss_thresholds=dict(authored_min=.90,rag_max=.15),near_miss_ids=near,cases=all_cases),ryr1_diagnostic=diagnostic)

def reparser_inventory(root):
    parser=load('79-reparse-normalisation.py');norm=load('74-input-normalisation.py');files=[]
    for name in parser.FILES:
        path=root/'data'/name;blob=json.loads(path.read_text());rows=blob['rows'] if isinstance(blob,dict) else blob
        copy=deepcopy(rows);tally=parser.reparse(copy,norm)
        files.append(dict(path=str(path.relative_to(root)),rows=len(rows),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            idempotent=rows==copy,raw_unchanged=[r.get('raw') for r in rows]==[r.get('raw') for r in copy],reparse_tally=tally))
    return dict(total_rows=sum(f['rows'] for f in files),files=files)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=OUT);args=ap.parse_args()
    result=analysis();result['reparser_inventory']=reparser_inventory(ROOT)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'publication-checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(screen=result['screen']['partition'],near_misses=result['screen']['near_miss_ids'],reparsed_rows=result['reparser_inventory']['total_rows']),indent=2))

if __name__=='__main__':main()
