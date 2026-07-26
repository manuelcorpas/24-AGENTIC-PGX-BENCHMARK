#!/usr/bin/env python3
"""
Hierarchical and multiway-clustered statistics (revision items N6 and N3;
Reviewer 2 point 5).

WHAT THE REVIEWER SAID
R2.5: a single case-level bootstrap understates uncertainty in a design where
the same nine models are evaluated on every case. If whole models are
systematically strong or weak, the model dimension carries variance that
case-clustering cannot see. The reviewer asked for a hierarchical model with
crossed effects for case, model, framing and replicate, and for model-clustered
intervals reported alongside case-clustered ones.

WHAT THIS DOES, AND WHAT IT HONESTLY DOES NOT
It reports three interval families for any headline proportion:
    case_clustered      resample cases (the original analysis)
    model_clustered     resample models
    two_way_clustered   resample both dimensions independently
and a variance decomposition across the four crossed factors.

It does NOT fit a generalised linear mixed model. statsmodels is not a
dependency of this package, and claiming a fitted mixed model while running a
method-of-moments decomposition would be a misdescription. The decomposition
below is a crossed one-way ANOVA-style estimate of between-level variance for
each factor, computed on factor-level means, and it is described as that in
STAR Methods. It answers the question the reviewer actually posed (how much
variance sits on each dimension) without overstating the machinery.

THE MEASUREMENT THAT MATTERS
Under deterministic execution the answer does not depend on which model
supplied the input call, so the between-model variance component must be
exactly zero. That turns the paper's model-invariance claim from a statement
into a measurement with a falsifier: any nonzero model component refutes it.

WHAT IS REPORTED ALONGSIDE
Coverage, abstention and parse failure, with missing, empty and unparsed
outputs counted as failures rather than dropped (R2.5 point 4). Aggregate A3 is
not produced: it is 1.0 by definition for non-lethal cases, so its aggregate
measures case composition. Lethal-class accuracy and error counts are reported
directly instead.

USAGE
    python code/65-hierarchical-stats.py --input data/v3_matched_factorial_scored.json
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent.parent
OUT_DEFAULT = BASE / "data" / "v3_hierarchical_stats.json"
REPORT_DEFAULT = BASE / "data" / "v3_hierarchical_stats.txt"

FACTORS = ("case_id", "model", "framing", "rep")
FACTOR_LABELS = {"case_id": "case", "model": "model", "framing": "framing", "rep": "replicate"}
N_BOOT = 2000


def _values(rows: list[dict], field: str = "correct") -> np.ndarray:
    return np.asarray([float(r.get(field) or 0.0) for r in rows], dtype=float)


def _groups(rows: list[dict], key: str) -> dict:
    out: dict = defaultdict(list)
    for r in rows:
        out[r[key]].append(r)
    return out


def _cluster_bootstrap(rows: list[dict], key: str, seed: int, field: str) -> dict:
    """Resample whole clusters with replacement; the mean is the statistic."""
    groups = list(_groups(rows, key).values())
    rng = np.random.default_rng(seed)
    point = float(np.mean(_values(rows, field))) if rows else 0.0
    if len(groups) < 2:
        return {"mean": round(point, 6), "ci95": (round(point, 6), round(point, 6)),
                "n_clusters": len(groups)}
    draws = np.empty(N_BOOT)
    idx = np.arange(len(groups))
    for b in range(N_BOOT):
        pick = rng.choice(idx, size=len(groups), replace=True)
        vals = np.concatenate([_values(groups[i], field) for i in pick])
        draws[b] = vals.mean()
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return {"mean": round(point, 6), "ci95": (round(float(lo), 6), round(float(hi), 6)),
            "n_clusters": len(groups)}


def _two_way_bootstrap(rows: list[dict], seed: int, field: str) -> dict:
    """Resample cases and models independently, keeping their intersection.

    The multiway analogue of a cluster bootstrap: dependence within a case and
    within a model are both respected, which a one-way bootstrap cannot do.
    """
    by_case = _groups(rows, "case_id")
    by_model = _groups(rows, "model")
    cases, models = list(by_case), list(by_model)
    point = float(np.mean(_values(rows, field))) if rows else 0.0
    if len(cases) < 2 or len(models) < 2:
        return {"mean": round(point, 6), "ci95": (round(point, 6), round(point, 6)),
                "n_cases": len(cases), "n_models": len(models)}
    cell: dict = defaultdict(list)
    for r in rows:
        cell[(r["case_id"], r["model"])].append(r)
    rng = np.random.default_rng(seed)
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        cs = rng.choice(len(cases), size=len(cases), replace=True)
        ms = rng.choice(len(models), size=len(models), replace=True)
        picked = [v for ci in cs for mi in ms
                  for v in _values(cell.get((cases[ci], models[mi]), []), field)]
        draws[b] = np.mean(picked) if picked else point
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return {"mean": round(point, 6), "ci95": (round(float(lo), 6), round(float(hi), 6)),
            "n_cases": len(cases), "n_models": len(models)}


def clustered_intervals(rows: list[dict], seed: int = 20260726,
                        field: str = "correct") -> dict:
    """Case-, model- and two-way-clustered 95% intervals for the mean."""
    return {
        "case_clustered": _cluster_bootstrap(rows, "case_id", seed, field),
        "model_clustered": _cluster_bootstrap(rows, "model", seed + 1, field),
        "two_way_clustered": _two_way_bootstrap(rows, seed + 2, field),
    }


def variance_components(rows: list[dict], field: str = "correct") -> dict:
    """Between-level variance for each crossed factor, plus the residual.

    Method of moments on factor-level means, not a fitted mixed model; see the
    module docstring. A factor whose levels all share the same mean contributes
    exactly zero, which is the property the model-invariance claim rests on.
    """
    out: dict[str, float] = {}
    values = _values(rows, field)
    grand = float(values.mean()) if len(values) else 0.0
    for key in FACTORS:
        if not rows or key not in rows[0]:
            out[FACTOR_LABELS.get(key, key)] = 0.0
            continue
        means = [float(_values(g, field).mean()) for g in _groups(rows, key).values()]
        out[FACTOR_LABELS[key]] = (round(float(np.var(means, ddof=0)), 12)
                                   if len(means) > 1 else 0.0)
    explained = sum(out.values())
    total = float(np.var(values, ddof=0)) if len(values) else 0.0
    out["residual"] = round(max(0.0, total - explained), 12)
    out["total"] = round(total, 12)
    out["grand_mean"] = round(grand, 6)
    return out


def report(rows: list[dict], seed: int = 20260726) -> dict:
    """Accuracy with clustered intervals, beside coverage and failure counts."""
    n = len(rows)
    emitted = [r for r in rows if r.get("parsed_ok") and not r.get("abstained")]
    lethal = [r for r in rows if r.get("lethal_action") is not None]
    out = {
        "n": n,
        "parse_failures": sum(1 for r in rows if not r.get("parsed_ok")),
        "abstentions": sum(1 for r in rows if r.get("abstained")),
        "coverage": round(len(emitted) / n, 6) if n else 0.0,
        "intervals": clustered_intervals(rows, seed),
        "variance_components": variance_components(rows),
        "lethal_n": len(lethal),
        "lethal_errors": sum(1 for r in lethal if not r["lethal_action"]),
        "lethal_action_accuracy": (round(sum(r["lethal_action"] for r in lethal) / len(lethal), 6)
                                   if lethal else None),
    }
    return out


def _rows_from_scored(payload) -> list[dict]:
    """Accept either raw scored rows or the 61- output, and normalise field names."""
    raw = payload["rows"] if isinstance(payload, dict) and "rows" in payload else payload
    rows = []
    for r in raw:
        rows.append({
            "case_id": r.get("case_id"), "model": r.get("model"),
            "rep": r.get("rep", 0), "framing": r.get("framing", "single"),
            "cell": r.get("cell", "all"),
            "correct": r.get("a1_phenotype", r.get("correct", 0.0)),
            "parsed_ok": r.get("parsed_ok", True),
            "abstained": r.get("abstained", False),
            "lethal_action": r.get("lethal_action"),
        })
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True,
                    help="scored rows (61-rescore-matched output, or any row list)")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--report", type=Path, default=REPORT_DEFAULT)
    args = ap.parse_args(argv)

    if not args.input.exists():
        sys.stderr.write(f"no input at {args.input}\n")
        return 1

    rows = _rows_from_scored(json.loads(args.input.read_text()))
    by_cell: dict = defaultdict(list)
    for r in rows:
        by_cell[r["cell"]].append(r)

    result = {cell: report(rs) for cell, rs in sorted(by_cell.items())}
    args.out.write_text(json.dumps(result, indent=2))

    lines = ["HIERARCHICAL / MULTIWAY-CLUSTERED STATISTICS", ""]
    for cell, s in result.items():
        i, v = s["intervals"], s["variance_components"]
        lines += [
            f"## {cell}   n={s['n']}  coverage={s['coverage']}  "
            f"parse_failures={s['parse_failures']}  abstentions={s['abstentions']}",
            f"   mean accuracy               {i['case_clustered']['mean']}",
            f"   case-clustered 95% CI       {i['case_clustered']['ci95']}",
            f"   model-clustered 95% CI      {i['model_clustered']['ci95']}",
            f"   two-way clustered 95% CI    {i['two_way_clustered']['ci95']}",
            f"   variance: case={v['case']}  model={v['model']}  "
            f"framing={v['framing']}  replicate={v['replicate']}  residual={v['residual']}",
            f"   lethal-class                {s['lethal_action_accuracy']} "
            f"({s['lethal_errors']} errors of {s['lethal_n']})",
            "",
        ]
    text = "\n".join(lines)
    args.report.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
