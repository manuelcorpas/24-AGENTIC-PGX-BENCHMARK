#!/usr/bin/env python3
"""Gates for Figure S8, the conditional-invariance control.

Written before the generator. Each test states a property the figure must
have, and each one failed against the submitted FigureS8.tiff, which is why
the figure is being rebuilt:

  1. Every bar height and every delta must be computed from a deposited file,
     never typed. The submitted figure printed a delta of 1.7 for
     skill-execution while its own bars give 1.6, because the bars and the
     delta came from different files.
  2. The title must not claim the gap vanishes under execution. Skill-execution
     carries the LARGEST gap of the five arms, so the submitted title asserted
     the opposite of its own picture.
  3. The title must not assert population-invariance unconditionally. The
     legend, and Minor 2, reserve invariance for the deterministic mapping
     conditional on identical valid inputs, which is the control arm alone.
  4. The panel must be the manuscript's eight models. Mistral Large 2 is
     excluded everywhere else in the paper and must be excluded here.

Run:  python3 code/test_figure_s8.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
GEN = Path(__file__).resolve().parent / "90-figure-s8-conditional-invariance.py"


def load_generator():
    spec = importlib.util.spec_from_file_location("figs8", GEN)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    failures: list[str] = []

    def check(name: str, condition: bool, detail: str = "") -> None:
        if condition:
            print(f"  pass  {name}")
        else:
            print(f"  FAIL  {name}{('  ' + detail) if detail else ''}")
            failures.append(name)

    if not GEN.exists():
        print(f"FAIL  generator not found: {GEN}")
        return 1

    fig = load_generator()
    arms = fig.compute()

    print("Figure S8 gates")

    # 1. Provenance. Every arm names the deposited file it came from, and that
    #    file is on disk.
    for arm in arms:
        check(
            f"{arm['key']} names an on-disk source",
            (BASE / "data" / arm["source"]).exists(),
            arm["source"],
        )

    # 2. The delta is the arithmetic of the bars, to the precision printed.
    for arm in arms:
        computed = arm["pops"]["EUR"]["pct"] - arm["pops"]["AFR"]["pct"]
        check(
            f"{arm['key']} delta equals EUR minus AFR",
            abs(arm["delta"] - computed) < 1e-9,
            f"printed {arm['delta']:.2f} vs bars {computed:.2f}",
        )

    # 3. The eight-model panel, Mistral excluded, on every arm.
    for arm in arms:
        check(
            f"{arm['key']} uses the eight-model panel",
            arm["n_models"] == 8,
            f"{arm['n_models']} models: {sorted(arm['models'])}",
        )
        check(
            f"{arm['key']} excludes Mistral Large 2",
            "Mistral Large 2" not in arm["models"],
        )

    # 4. The control is the only arm that is flat, which is the whole point of
    #    calling it a control. If another arm were flat too, the figure would
    #    not distinguish construction from measurement.
    control = [a for a in arms if a["key"] == "control"]
    check("there is exactly one control arm", len(control) == 1)
    if control:
        check(
            "the control arm is exactly flat",
            control[0]["delta"] == 0.0,
            f"delta {control[0]['delta']}",
        )
    non_control = [a for a in arms if a["key"] != "control"]
    check(
        "every non-control arm carries a non-zero gap",
        all(abs(a["delta"]) > 0.0 for a in non_control),
        str([(a["key"], round(a["delta"], 2)) for a in non_control]),
    )

    # 5. The title. The three claims the submitted artwork made and could not
    #    support.
    title = fig.TITLE + " " + fig.SUBTITLE
    lower = title.lower()
    check(
        "the title does not claim the gap vanishes under execution",
        "vanish" not in lower,
        fig.TITLE,
    )
    check(
        "the title does not assert unconditional population-invariance",
        "population-invariance" not in lower and "population invariance" not in lower,
        fig.TITLE,
    )
    check(
        "the title says this is a control",
        "control" in lower,
        fig.TITLE,
    )
    check(
        "the title carries the conditional",
        "conditional" in lower or "identical valid" in lower,
        fig.TITLE,
    )

    # 6. No long dashes anywhere in the artwork strings, per house style. The
    #    minus sign in "EUR - AFR" is a real minus, not an en dash.
    strings = [fig.TITLE, fig.SUBTITLE, fig.YLABEL] + [a["label"] for a in arms]
    for s in strings:
        check(
            f"no long dash in {s[:40]!r}",
            not any(ch in s for ch in "—–"),
        )

    # 7. The largest gap must not be in the control, and the figure must be
    #    honest about where it is. This is the assertion the submitted title
    #    inverted.
    worst = max(non_control, key=lambda a: abs(a["delta"]))
    check(
        "skill-execution carries the largest gap, as the data says",
        worst["key"] == "skill_execution",
        f"largest is {worst['key']} at {worst['delta']:.2f}",
    )

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print(f"all {len(arms)} arms pass every gate")
    return 0


if __name__ == "__main__":
    sys.exit(main())
