#!/usr/bin/env python3
"""Render the candidate evidence ledger from the derived data and reviewed dispositions."""
from pathlib import Path
import json

REV=Path(__file__).resolve().parents[1]

def candidate_audit():
    a=json.loads((REV/'DATA/analysis.json').read_text())
    evidence=json.loads((REV/'DATA/candidate-dispositions.json').read_text())
    lines=['# Post hoc candidate audit, 16 September 2026','',
      'Screen: all 110 original-reference cases. Select baseline phenotype means >=0.95 in each of RAG-assisted execution, authored-rule generation and authored-rule execution, and <=0.10 in RAG generation. Eight cases meet the rule. This is a full-case screen followed by candidate review, not independent clinical adjudication of the benchmark. Dispositions concern the observed phenotype discrepancy; they are not blanket approval of every recommendation.',
      '', 'Outcome: one confirmed reference error (RYR1), one outdated retrieval interpretation (CYP2D6 activity score 1.0), and six terminology mismatches. No scoring patterns were changed. All 120 observed responses per candidate, including 24 from each configuration, are retained in DATA/analysis.json. The following examples are the fixed GPT-5.2 replicate 0 pair; complete frequency counts prevent those examples being mistaken for exhaustive response coverage.',
      '', 'Evidence: primary guideline pages were read in the browser, with DOI resolution attempted directly. The CYP2D6 consensus was checked at its Wiley DOI page after resolver access failed. The source-linked dispositions below are reasoned candidate review. RYR1, CYP2D6-amitriptyline and HLA-B*57:01-positive were also hand-checked against the preserved raw responses and references. No new model calls were made.',
      '', 'The terminology cases have frozen equivalence totals 24, 24, 23, 22, 19 and 0 of 24 in the order reported in the manuscript. MT-RNR1 receives zero because the frozen equivalence layer does not recognise that risk phrasing. A terminology disposition does not imply that every response parses or earns equivalence credit.']
    for c in a['candidates']:
        e=evidence[c['case_id']];ref=c['reference']
        lines.extend(['',f"## Case: {c['case_id']}",'',f"Gene/drug: {ref['gene']} / {ref['drug']}. Input: {ref['genotype']}. Authored diplotype: {ref['gt_diplotype']}. Original phenotype: {ref['gt_phenotype']}.",'',f"Disposition: **{e['disposition']}**. {e['evidence']}",'',f"Source: {e['source']}",'','| Configuration | Baseline phenotype mean |','|---|---:|'])
        for cell,value in c['phenotype_means'].items():lines.append(f'| {cell} | {value:.4f} |')
        lines.extend(['',f"Frozen RAG equivalence: {int(c['rag_equivalence_sum'])}/24 ({c['rag_equivalence_sum']/24:.2f}).",'', 'Observed RAG phenotype strings and frequencies:',''])
        for text,n in c['observed_phenotypes']['rag_generation'].items():lines.append(f'- {n}/24: {text or "[empty]"}')
        for cell in ('rag_generation','skill_execution'):
            r=next(r for r in c['responses'] if r['cell']==cell and r['model']=='GPT-5.2' and r['rep']==0)
            lines.extend(['',f"Example {cell}, GPT-5.2, replicate 0:",'','```text',r['raw'],'```',f"Parsed or executed phenotype: {r['parsed_phenotype']}",f"Parsed or executed recommendation: {r['parsed_drug']}"])
    return '\n'.join(lines)+'\n'

def main():
    (REV/'CANDIDATE-AUDIT.md').write_text(candidate_audit())

if __name__=='__main__':main()
