#!/usr/bin/env python3
# EXPLORATORY, POST HOC, NOT CLINICALLY VALIDATED.
# Written after the primary results were seen. Not part of the frozen scoring
# pipeline and not an endpoint. See README.md in this directory.
"""Figure S10. Directive strength of high-severity responses, by configuration.

Every value is computed from the deposited run records by 94-directive-strength.py.
"""
from __future__ import annotations
import argparse, importlib.util, json
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE = Path(__file__).resolve().parents[2]
CELLS = [("free_generation","free\ngeneration"),("rag_generation","RAG\ngeneration"),
         ("rag_execution","RAG-assisted\nexecution"),("skill_generation","authored-rule\ngeneration"),
         ("skill_execution","authored-rule\nexecution")]
COL = {4:"#1b6b3a", 3:"#e0a busy", 2:"#c8842a", 1:"#8a6fb0", 0:"#b3261e"}
COL = {4:"#1b6b3a", 3:"#e8b04b", 2:"#c8842a", 1:"#8a6fb0", 0:"#b3261e"}
NAME = {4:"categorical prohibition", 3:"gated on a condition", 2:"recommends against",
        1:"dose modification", 0:"no avoidance language"}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--out",default="FigureS10.tiff")
    ap.add_argument("--base",default=str(BASE)); ap.add_argument("--dpi",type=int,default=600)
    a=ap.parse_args(); base=Path(a.base)
    spec=importlib.util.spec_from_file_location("ds",base/"code"/"exploratory-posthoc"/"94-directive-strength.py")
    ds=importlib.util.module_from_spec(spec); spec.loader.exec_module(ds)
    live=json.loads((base/"data"/"v3_five_cell_live.json").read_text())
    scored=json.loads((base/"data"/"v3_matched_scored_rows_all5.json").read_text())
    cases={c["id"]:c for c in json.loads((base/"specs"/"test_cases_v3.json").read_text())}
    K=lambda r:(r["case_id"],r["cell"],r["model"],r["rep"]); LIVE={K(r):r for r in live}
    data={}; unres={}
    for cell,_ in CELLS:
        rows=[s for s in scored if s.get("lethal_action") is not None and s["cell"]==cell]
        lv=[ds.directive_strength((LIVE.get(K(s)) or {}).get("parsed_drug") or "") for s in rows]
        data[cell]=[100*lv.count(l)/len(lv) for l in [4,3,2,1,0]]
        unres[cell]=sum(1 for s in rows
                        if ds.classify(cases[s["case_id"]],(LIVE.get(K(s)) or {}).get("parsed_drug") or "")==ds.UNRESOLVED)
    fig,ax=plt.subplots(figsize=(7.2,3.9),dpi=a.dpi)
    xs=range(len(CELLS)); bottom=[0]*len(CELLS)
    for idx,l in enumerate([4,3,2,1,0]):
        vals=[data[c][idx] for c,_ in CELLS]
        ax.bar(xs,vals,bottom=bottom,color=COL[l],edgecolor="white",linewidth=.8,width=.68)
        for x,(v,b) in enumerate(zip(vals,bottom)):
            if v>=4.5: ax.text(x,b+v/2,f"{v:.0f}%",ha="center",va="center",fontsize=7.4,
                               color="white",fontweight="bold")
        bottom=[b+v for b,v in zip(bottom,vals)]
    ax.set_xticks(list(xs)); ax.set_xticklabels([lab for _,lab in CELLS],fontsize=8)
    ax.set_ylabel("% of high-severity responses",fontsize=8.5); ax.set_ylim(0,100)
    ax.tick_params(axis="y",labelsize=8)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    for x,(c,_) in enumerate(CELLS):
        ax.text(x,101.5,f"unresolved: {unres[c]}",ha="center",fontsize=7,color="#5b5b5b")
    ax.legend(handles=[Patch(facecolor=COL[l],label=NAME[l]) for l in [4,3,2,1,0]],
              fontsize=7.2,frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,-.16))
    fig.tight_layout(rect=(0,.09,1,.97))
    fig.savefig(a.out,dpi=a.dpi,format="tiff",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(Path(a.out).with_suffix(".png"),dpi=200)
    print("wrote",a.out)
    for c,_ in CELLS: print(f"   {c:20s} "+" ".join(f"{v:5.1f}" for v in data[c])+f"   unresolved {unres[c]}")

if __name__=="__main__": main()
