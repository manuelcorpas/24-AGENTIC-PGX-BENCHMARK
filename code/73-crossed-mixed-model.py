#!/usr/bin/env python3
"""
Crossed random-effects model for phenotype accuracy (Reviewer 2, point 5).

WHY THIS EXISTS
R2.5 asked for multiway clustering OR a hierarchical model with crossed effects
for case, model, framing and replicate. We supplied the clustering and a
method-of-moments variance decomposition, and declined to fit a mixed model on
the grounds that statsmodels was not a dependency and that calling a moments
estimate a fitted model would be a misdescription. That reasoning was sound but
it left the referee able to say we answered a request with a refusal. This fits
the model.

WHAT IS FITTED
A mixed model per cell, with the scored outcome as response, a random intercept
for case and a random intercept for model. statsmodels' MixedLM fits one grouping
factor at a time, so the crossed structure is obtained with `vc_formula` on a
constant group: case and model both enter as variance components, which is the
standard way to express crossed effects in this library.

Framing and replicate are omitted from the fitted model deliberately. This design
has a single framing, so that component is zero by construction, and the
replicate component in the moments decomposition is of order 1e-5, three orders
below the smallest other term. Fitting terms known to be zero invites a singular
fit and reports a precision the data do not carry.

WHAT IT DOES NOT CLAIM
The response is a 0/1 correctness score fitted with a linear mixed model, not a
GLMM with a binomial link. On proportions away from the boundary the variance
components are interpretable on the probability scale and comparable with the
moments decomposition, which is the comparison being asked for. Cells whose
accuracy sits near 1.0 are close enough to the boundary that the linear
approximation is strained, and the report says so per cell rather than hiding it.

USAGE
    .venv-stats/bin/python code/73-crossed-mixed-model.py
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
ROWS = BASE / "data" / "v3_five_cell_live_rows.json"
STATS = BASE / "data" / "v3_five_cell_live_stats.json"
OUT = BASE / "data" / "v3_crossed_mixed_model.json"
REPORT = BASE / "data" / "v3_crossed_mixed_model.txt"


def fit_cell(df, sm, smf):
    """Random intercepts for case and model, crossed."""
    df = df.copy()
    df["grp"] = 1
    vc = {"case": "0 + C(case_id)", "model": "0 + C(model)"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        md = smf.mixedlm("y ~ 1", df, groups="grp", vc_formula=vc)
        res = md.fit(reml=True, method="lbfgs")
    vcomp = dict(zip(md.exog_vc.names, res.vcomp)) if hasattr(md, "exog_vc") else {}
    return {
        "n": int(len(df)),
        "intercept": float(res.fe_params.iloc[0]),
        "var_case": float(vcomp.get("case", float("nan"))),
        "var_model": float(vcomp.get("model", float("nan"))),
        "var_residual": float(res.scale),
        "converged": bool(res.converged),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--rows", type=Path, default=ROWS)
    ap.add_argument("--stats", type=Path, default=STATS)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--report", type=Path, default=REPORT)
    a = ap.parse_args(argv)

    try:
        import pandas as pd
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
    except ImportError:
        sys.stderr.write(
            "statsmodels and pandas are required. This script is run from the "
            "dedicated environment:\n  .venv-stats/bin/python code/73-crossed-mixed-model.py\n")
        return 1

    if not a.rows.exists():
        sys.stderr.write(f"missing {a.rows}; run 61-rescore-matched.py --rows\n")
        return 1
    rows = json.loads(a.rows.read_text())
    df = pd.DataFrame([{"cell": r["cell"], "case_id": r["case_id"],
                        "model": r["model"], "y": float(r["a1_phenotype"])}
                       for r in rows])

    moments = {}
    if a.stats.exists():
        st = json.loads(a.stats.read_text())
        for cell, v in (st.items() if isinstance(st, dict) else []):
            comp = (v or {}).get("variance_components") or {}
            if comp:
                moments[cell] = comp

    result = {}
    for cell, sub in df.groupby("cell"):
        try:
            result[cell] = fit_cell(sub, sm, smf)
        except Exception as exc:                                   # noqa: BLE001
            result[cell] = {"error": f"{type(exc).__name__}: {exc}"}
        m = moments.get(cell) or {}
        result[cell]["moments_case"] = m.get("case")
        result[cell]["moments_model"] = m.get("model")

    a.out.write_text(json.dumps(result, indent=2))

    L = ["CROSSED RANDOM-EFFECTS MODEL, PHENOTYPE ACCURACY", ""]
    L.append("  Random intercepts for case and model, fitted by REML. Compared with the")
    L.append("  method-of-moments decomposition reported in Table 3.")
    L.append("")
    L.append(f"  {'cell':18s} {'mean':>7} {'var(case)':>11} {'var(model)':>11} "
             f"{'moments case':>13} {'moments model':>14}")
    for cell, v in sorted(result.items()):
        if v.get("error"):
            L.append(f"  {cell:18s} FIT FAILED: {v['error']}")
            continue
        L.append(f"  {cell:18s} {v['intercept']:7.3f} {v['var_case']:11.5f} "
                 f"{v['var_model']:11.5f} {str(v['moments_case'] or '')[:11]:>13} "
                 f"{str(v['moments_model'] or '')[:12]:>14}")
    L.append("")
    near = [c for c, v in result.items()
            if not v.get("error") and v.get("intercept", 0) > 0.95]
    if near:
        L.append("  Cells with mean accuracy above 0.95 sit close to the boundary, where a")
        L.append("  linear mixed model understates variance on the probability scale: "
                 + ", ".join(sorted(near)) + ".")
        L.append("  Their components should be read as indicative and compared with the")
        L.append("  clustered intervals rather than relied on alone.")
        L.append("")
    L.append("  Framing and replicate are not fitted: this design has a single framing, so")
    L.append("  that component is zero by construction, and the replicate component in the")
    L.append("  moments decomposition is of order 1e-5. Fitting terms known to be zero")
    L.append("  invites a singular fit and reports precision the data do not carry.")
    text = "\n".join(L)
    a.report.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
