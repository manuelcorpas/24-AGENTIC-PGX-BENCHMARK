#!/usr/bin/env python3
"""Recompute the R3 analyses and figures from preserved data, without model calls."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys
from importlib.util import module_from_spec,spec_from_file_location

ROOT=Path(__file__).resolve().parents[3]
REV=Path(__file__).resolve().parents[1]

def commands(python):
    data=REV/'DATA';review=REV/'REVIEW'
    return [[python,str(ROOT/'code/98-third-revision.py')],
        [python,str(ROOT/'code/65-hierarchical-stats.py'),'--input',str(data/'attempted-scored-rows.json'),'--out',str(data/'hierarchical-stats.json'),'--report',str(review/'hierarchical-stats.txt')],
        [python,str(ROOT/'code/73-crossed-mixed-model.py'),'--rows',str(data/'attempted-scored-rows.json'),'--stats',str(data/'hierarchical-stats.json'),'--out',str(data/'crossed-mixed-model.json'),'--report',str(review/'crossed-mixed-model.txt')],
        [python,str(ROOT/'code/88-supp-normalisation-by-gene.py'),'--out',str(review/'normalisation-by-gene.txt'),'--json-out',str(data/'normalisation-by-gene.json')]]

def load(name):
    s=spec_from_file_location(name,ROOT/'code'/name);m=module_from_spec(s);s.loader.exec_module(m);return m

def figures():
    import matplotlib.pyplot as plt
    from PIL import Image
    out=REV/'FIGURES';out.mkdir(parents=True,exist_ok=True)
    a=json.loads((REV/'DATA/analysis.json').read_text())
    rows=json.loads((REV/'DATA/attempted-scored-rows.json').read_text())
    m=load('70-figures-revision.py');m.style()
    models=sorted({r['model'] for r in rows})
    by={c:{model:sum(r['a1_phenotype'] for r in rows if r['cell']==c and r['model']==model)/330 for model in models} for c in m.CELLS}
    h=load('86-figures-s2-s3-five-cell.py');high=[r for r in rows if r['lethal_action'] is not None]
    pop=load('90-figure-s8-conditional-invariance.py')
    worked=load('93-figure-s9-worked-case.py');case,examples=worked.collect(ROOT,worked.CASE_ID,worked.MODEL,worked.REP);gene,rates=worked.panel_error_rates(ROOT,worked.CASE_ID)
    fs={'Figure2':m.figure_r1(a['baseline'],save_output=False),
        'FigureS2':h.figure_s2(high,save_output=False),
        'FigureS3':h.figure_s3(high,save_output=False),
        'FigureS4':m.figure_r2({},rows=rows,already_scored=True,save_output=False,caller=json.loads((REV/'DATA/caller-comparison.json').read_text())),
        'FigureS5':m.figure_r3(by,save_output=False),
        'FigureS6':m.figure_r4(save_output=False),
        'FigureS8':pop.draw(pop.compute(),save_output=False),
        'FigureS9':worked.render(case,examples,gene,rates,out/'FigureS9.tiff',worked.MODEL,worked.REP,save_output=False)}
    for name,fig in fs.items():
        for ext in ('png','pdf','tiff'):
            kwargs={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
            fig.savefig(out/(name+'.'+ext),dpi=600 if ext=='tiff' else 300,facecolor='white',**kwargs)
        plt.close(fig)
    with Image.open(out/'Figure2.png') as im:im.convert('L').save(out/'Figure2-grayscale.png')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--figures-only',action='store_true');args=ap.parse_args()
    if not args.figures_only:
        (REV/'REVIEW').mkdir(parents=True,exist_ok=True)
        for cmd in commands(sys.executable):subprocess.run(cmd,cwd=ROOT,check=True)
    figures()

if __name__=='__main__':main()
