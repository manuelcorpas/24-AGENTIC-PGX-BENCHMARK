#!/usr/bin/env python3
"""
Case-level cluster-bootstrap statistics for the v26 manuscript.

Reproduces the clustered confidence intervals and contrasts reported in
"Trustworthy agentic genomics through versioned skill libraries" (Corpas et al.,
submitted to Cell Genomics):

  - headline phenotype accuracy (mean A1) per condition, with naive binomial
    Wilson intervals AND case-level cluster-bootstrap 95% intervals + design effect;
  - the free-prompted -> retrieval lethal-class contrast, clustered over the 14
    distinct lethal-class cases (bootstrap CI and P; the nominal two-proportion
    z-test is a pseudoreplication artefact and is not used in the manuscript);
  - the mechanism-specific (HLA risk-allele vs non-HLA) lethal-class breakdown
    behind Table S5;
  - the all-nine-model skill-arm aggregate and skill-arm lethal-class counts.

Why clustering: the 110 cases are each evaluated across 9 models, 3 population
framings of identical genotypes and 3 replicates, so the per-cell observations
are not independent. Naive binomial intervals understate uncertainty by a design
factor of ~23-43; the case-level cluster bootstrap is the honest measure.

Inputs (place the Zenodo archive contents in ../data/; see data/README.md):
  ../data/v3_raw_rescored_three_arm.json   three core conditions (no_spec, cpic_rag, with_spec)
  ../data/v3_three_arm_per_case_a1.csv      per-case A1 and the lethal-class flag
  ../data/v3_armAB_fullgrid.json            skill arm (A_reasoning, B_execution), all 9 models
                                            [needed for the skill-arm rows; if absent the
                                            three-arm statistics are still produced]

Outputs:
  ../data/v26_cluster_stats.json   machine-readable
  ../data/v26_cluster_stats.txt    human-readable

Run:  python3 code/11b-cluster-bootstrap-stats.py
Deterministic: seeded bootstrap (numpy). scipy is optional (used only for the
Wilcoxon signed-rank check); the bootstrap runs without it.
"""
import os
import json
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
SEED = 20260611
B = 10000  # bootstrap resamples


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def load_json(name):
    return json.loads((DATA / name).read_text())


def mean_A1(cells):
    v = [c["A1"] for c in cells if c["parsed"] and c["A1"] is not None]
    return (sum(v) / len(v) if v else float("nan")), len(v)


def lethal_A3_errors(cells):
    leth = [c for c in cells if c["lethal"] and c["parsed"] and c["A3"] is not None]
    return sum(1 for c in leth if c["A3"] < 1.0), len(leth)


def by_cluster(cells, key):
    d = defaultdict(list)
    for c in cells:
        d[c[key]].append(c)
    return d


def boot_mean_A1(cells, key, rng):
    groups = list(by_cluster(cells, key).values())
    g = len(groups)
    ests = []
    for _ in range(B):
        idx = rng.integers(0, g, g)
        s = 0.0
        n = 0
        for i in idx:
            for c in groups[i]:
                if c["parsed"] and c["A1"] is not None:
                    s += c["A1"]
                    n += 1
        if n:
            ests.append(s / n)
    return np.percentile(ests, [2.5, 97.5])


def boot_lethal_diff(a, b, key, rng):
    ga = by_cluster([c for c in a if c["lethal"]], key)
    gb = by_cluster([c for c in b if c["lethal"]], key)
    keys = sorted(set(ga) | set(gb))
    k = len(keys)
    ea, na = lethal_A3_errors(a)
    eb, nb = lethal_A3_errors(b)
    diffs = []
    for _ in range(B):
        idx = rng.integers(0, k, k)
        ea_ = na_ = eb_ = nb_ = 0
        for i in idx:
            kk = keys[i]
            for c in ga.get(kk, []):
                if c["parsed"] and c["A3"] is not None:
                    na_ += 1
                    ea_ += c["A3"] < 1.0
            for c in gb.get(kk, []):
                if c["parsed"] and c["A3"] is not None:
                    nb_ += 1
                    eb_ += c["A3"] < 1.0
        if na_ and nb_:
            diffs.append(eb_ / nb_ - ea_ / na_)
    diffs = np.array(diffs)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p_two = 2 * min(float(np.mean(diffs <= 0)), float(np.mean(diffs >= 0)))
    return dict(rate_a=round(100 * ea / na, 1), rate_b=round(100 * eb / nb, 1),
                err_a=ea, n_a=na, err_b=eb, n_b=nb,
                diff_pp=round(100 * (eb / nb - ea / na), 1),
                ci_pp=[round(100 * lo, 1), round(100 * hi, 1)],
                boot_p_two_sided=round(p_two, 4), n_clusters=k)


def main():
    rng = np.random.default_rng(SEED)
    if not DATA.exists():
        raise SystemExit(f"Data directory not found: {DATA}. "
                         "Unpack the Zenodo archive into ../data/ (see data/README.md).")

    # lethal flag per case from the per-case CSV
    lethal_case = {}
    with open(DATA / "v3_three_arm_per_case_a1.csv") as f:
        for row in csv.DictReader(f):
            lethal_case[row["case_id"]] = int(row["lethal"])
    n_lethal_cases = sum(lethal_case.values())

    # ---- three core conditions (mean-A1 headline; parsed = not format_fail) ----
    def core_cells(cond):
        out = []
        for r in load_json("v3_raw_rescored_three_arm.json"):
            if r["cond"] != cond:
                continue
            sc = r.get("scores") or {}
            out.append(dict(cond=cond, tc=r["tc"], model=r["model"], pop=r["pop"],
                            gene=r.get("gene", ""),
                            lethal=lethal_case.get(r["tc"], 0),
                            parsed=(not sc.get("format_fail", False)),
                            A1=sc.get("A1"), A3=sc.get("A3")))
        return out

    free_ce = core_cells("no_spec")
    rag = core_cells("cpic_rag")
    control_ce = core_cells("with_spec")

    # ---- skill arms (optional file) ----
    have_skill = (DATA / "v3_armAB_fullgrid.json").exists()
    reason8 = exec8 = reason9 = exec9 = []
    if have_skill:
        fg = load_json("v3_armAB_fullgrid.json")

        def skill_cells(arm, include_mistral):
            out = []
            for r in fg:
                if r["arm"] != arm:
                    continue
                if (not include_mistral) and "istral" in r["model"]:
                    continue
                raw = (r.get("raw_text") or "")
                out.append(dict(cond=arm, tc=r["tc"], model=r["model"], pop=r["pop"],
                                gene=r.get("gene", ""), lethal=bool(r.get("lethal")),
                                parsed=(raw.strip() not in ("", "<error: 'choices'>")),
                                A1=r.get("A1"), A3=r.get("A3"), raw=raw))
            return out
        reason8 = skill_cells("A_reasoning", False)
        exec8 = skill_cells("B_execution", False)
        reason9 = skill_cells("A_reasoning", True)
        exec9 = skill_cells("B_execution", True)

    # ---- clustered accuracy ----
    clustered = {}
    rows = [("free_prompted", free_ce), ("retrieval", rag), ("control", control_ce)]
    if have_skill:
        rows += [("skill_reasoning_8mdl", reason8), ("skill_execution_8mdl", exec8)]
    for name, cells in rows:
        m, n = mean_A1(cells)
        k = round(m * n)
        nw = wilson(k, n)
        lo_c, hi_c = boot_mean_A1(cells, "tc", rng)
        lo_m, hi_m = boot_mean_A1(cells, "model", rng)
        nw_w, cl_w = (nw[1] - nw[0]), (hi_c - lo_c)
        clustered[name] = dict(point=round(100 * m, 1), n_cells=n,
                               n_cases=len(set(c["tc"] for c in cells)),
                               naive_wilson=[round(100 * nw[0], 1), round(100 * nw[1], 1)],
                               cluster_by_case=[round(100 * lo_c, 1), round(100 * hi_c, 1)],
                               cluster_by_model=[round(100 * lo_m, 1), round(100 * hi_m, 1)],
                               design_effect=round((cl_w / nw_w) ** 2, 1) if nw_w > 0 else None)

    # ---- lethal-class contrast + mechanism breakdown (Table S5) ----
    lethal_contrast = boot_lethal_diff(free_ce, rag, "tc", rng)

    def pooled_lethal(cells, pred):
        sub = [c for c in cells if c["lethal"] and pred(c)]
        e, n = lethal_A3_errors(sub)
        return dict(errors=e, n=n, rate=round(100 * e / n, 1) if n else None)
    is_hla = lambda c: str(c["gene"]).upper().startswith("HLA")
    mechanism = dict(
        hla=dict(free=pooled_lethal(free_ce, is_hla), retrieval=pooled_lethal(rag, is_hla)),
        non_hla=dict(free=pooled_lethal(free_ce, lambda c: not is_hla(c)),
                     retrieval=pooled_lethal(rag, lambda c: not is_hla(c))),
    )
    try:
        from scipy import stats as _st
        gf = by_cluster([c for c in free_ce if c["lethal"]], "tc")
        gr = by_cluster([c for c in rag if c["lethal"]], "tc")
        cases = sorted(set(gf) | set(gr))

        def rate(g):
            e, n = lethal_A3_errors(g)
            return 100 * e / n if n else float("nan")
        fr = np.array([rate(gf.get(c, [])) for c in cases])
        rr = np.array([rate(gr.get(c, [])) for c in cases])
        ok = ~(np.isnan(fr) | np.isnan(rr))
        w = _st.wilcoxon(rr[ok], fr[ok], zero_method="wilcox")
        mechanism["wilcoxon_p"] = round(float(w.pvalue), 3)
    except Exception:
        mechanism["wilcoxon_p"] = None

    # ---- skill-arm counts (optional) ----
    skill = None
    if have_skill:
        def agg(cells, parsed_only):
            v = [c["A1"] for c in cells if c["A1"] is not None and (c["parsed"] or not parsed_only)]
            return round(100 * sum(v) / len(v), 1), len(v)

        def per_model(cells):
            d = defaultdict(lambda: dict(s=0.0, np_=0, nt=0))
            for c in cells:
                if c["A1"] is None:
                    continue
                d[c["model"]]["nt"] += 1
                if c["parsed"]:
                    d[c["model"]]["np_"] += 1
                    d[c["model"]]["s"] += c["A1"]
            return {m: dict(overall_pct=round(100 * v["s"] / v["nt"], 1) if v["nt"] else None,
                            responding_pct=round(100 * v["s"] / v["np_"], 1) if v["np_"] else None,
                            response_rate_pct=round(100 * v["np_"] / v["nt"], 1) if v["nt"] else None)
                    for m, v in sorted(d.items())}
        skill = dict(
            reasoning_8mdl=mean_A1(reason8)[0],
            execution_8mdl=mean_A1(exec8)[0],
            reasoning_9mdl_error_as_wrong=agg(reason9, False)[0],
            reasoning_9mdl_responding_only=agg(reason9, True)[0],
            execution_9mdl_error_as_wrong=agg(exec9, False)[0],
            execution_9mdl_responding_only=agg(exec9, True)[0],
            lethal_reasoning_8mdl=dict(zip(("errors", "n"), lethal_A3_errors(reason8))),
            lethal_execution_8mdl=dict(zip(("errors", "n"), lethal_A3_errors(exec8))),
            per_model_reasoning=per_model(reason9),
            per_model_execution=per_model(exec9),
        )
        skill["reasoning_8mdl"] = round(100 * skill["reasoning_8mdl"], 1)
        skill["execution_8mdl"] = round(100 * skill["execution_8mdl"], 1)

    out = dict(meta=dict(B=B, seed=SEED, n_total_cases=len(lethal_case),
                         n_lethal_cases=n_lethal_cases, skill_arm_available=have_skill),
               clustered_accuracy=clustered,
               lethal_contrast_free_to_retrieval=lethal_contrast,
               mechanism_table_s5=mechanism,
               skill_arm=skill)
    (DATA / "v26_cluster_stats.json").write_text(json.dumps(out, indent=2))

    # human-readable
    L = [f"V26 CLUSTER-BOOTSTRAP STATISTICS  (B={B}, seed={SEED})", "=" * 64]
    if not have_skill:
        L.append("NOTE: v3_armAB_fullgrid.json not found in ../data/; skill-arm rows skipped.")
    L.append("\nClustered accuracy (naive Wilson vs case-level cluster bootstrap):")
    for name, v in clustered.items():
        L.append(f"  {name}: {v['point']}%  Wilson {v['naive_wilson']}  "
                 f"cluster-by-case {v['cluster_by_case']}  design effect x{v['design_effect']}")
    lc = lethal_contrast
    L.append("\nLethal-class contrast free -> retrieval (clustered by case):")
    L.append(f"  free {lc['err_a']}/{lc['n_a']}={lc['rate_a']}%  retrieval {lc['err_b']}/{lc['n_b']}={lc['rate_b']}%"
             f"  diff +{lc['diff_pp']}pp  95% CI {lc['ci_pp']}pp  bootstrap P={lc['boot_p_two_sided']}"
             f"  ({lc['n_clusters']} lethal clusters)")
    mh, mn = mechanism["hla"], mechanism["non_hla"]
    L.append("\nMechanism (Table S5): lethal-class A3 error rate by locus class")
    L.append(f"  HLA risk-allele:   free {mh['free']['rate']}%  retrieval {mh['retrieval']['rate']}%")
    L.append(f"  non-HLA:           free {mn['free']['rate']}%  retrieval {mn['retrieval']['rate']}%")
    L.append(f"  Wilcoxon signed-rank P (per case): {mechanism['wilcoxon_p']}")
    if skill:
        L.append("\nSkill-arm aggregate:")
        L.append(f"  reasoning  8-model {skill['reasoning_8mdl']}% | 9-model responding-only "
                 f"{skill['reasoning_9mdl_responding_only']}% | 9-model error-as-wrong {skill['reasoning_9mdl_error_as_wrong']}%")
        L.append(f"  execution  8-model {skill['execution_8mdl']}% | 9-model responding-only "
                 f"{skill['execution_9mdl_responding_only']}% | 9-model error-as-wrong {skill['execution_9mdl_error_as_wrong']}%")
        L.append(f"  skill-arm lethal-class: reasoning {skill['lethal_reasoning_8mdl']['errors']}/"
                 f"{skill['lethal_reasoning_8mdl']['n']}, execution {skill['lethal_execution_8mdl']['errors']}/"
                 f"{skill['lethal_execution_8mdl']['n']}")
    (DATA / "v26_cluster_stats.txt").write_text("\n".join(L))
    print("\n".join(L))
    print(f"\nWrote {DATA/'v26_cluster_stats.json'} and .txt")


if __name__ == "__main__":
    main()
