#!/usr/bin/env python3
"""Generate the final bidirectional adversarial-contract figure.

The plot is derived directly from the stored forward and reverse response rows.
It reports whether models followed the deliberately corrupted contract, hedged
while retaining its unsafe fields, or reverted to the canonical CPIC answer.

Usage:
    python3 code/80-figure6-adversarial-contract.py --self-test
    python3 code/80-figure6-adversarial-contract.py
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
FIGURES = BASE / "figures"
FORWARD = DATA / "v3_adversarial_scrambled.json"
REVERSE = DATA / "v3_adversarial_reverse.json"

ECHO = "#d1495b"
HEDGE = "#e79aa8"
REVERT = "#2d7f7c"


def load(path: Path) -> list[dict]:
    if not path.exists():
        raise FileNotFoundError(
            f"missing {path}; download the deposited raw rows into data/ first")
    return json.loads(path.read_text())


def counts(rows: list[dict]) -> tuple[int, int, int]:
    c = Counter(r.get("classification") for r in rows)
    echo = c.get("ECHO_SCRAMBLED", 0)
    hedge = c.get("HEDGE", 0)
    reverted = len(rows) - echo - hedge
    return echo, hedge, reverted


def self_test() -> int:
    failures: list[str] = []
    try:
        forward = load(FORWARD)
        reverse = load(REVERSE)
    except FileNotFoundError as exc:
        print(f"SELF-TEST FAIL: {exc}")
        return 1

    expected = {
        "forward": (45, (43, 2, 0)),
        "reverse": (45, (45, 0, 0)),
    }
    for label, rows in (("forward", forward), ("reverse", reverse)):
        n, want = expected[label]
        got = counts(rows)
        if len(rows) != n:
            failures.append(f"{label} row count: got {len(rows)}, want {n}")
        if got != want:
            failures.append(f"{label} classifications: got {got}, want {want}")

    if failures:
        for failure in failures:
            print(f"SELF-TEST FAIL: {failure}")
        return 1
    print("SELF-TEST PASSED: forward 43 echo + 2 hedge; reverse 45 echo; zero reversions")
    return 0


def render() -> None:
    forward = load(FORWARD)
    reverse = load(REVERSE)
    vals = [counts(forward), counts(reverse)]
    echoed = [v[0] for v in vals]
    hedged = [v[1] for v in vals]
    reverted = [v[2] for v in vals]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
    })
    fig, ax = plt.subplots(figsize=(7.4, 4.45))
    x = [0, 1]
    width = 0.55
    ax.bar(x, echoed, width, color=ECHO,
           label="Followed corrupted contract (verbatim)", zorder=3)
    ax.bar(x, hedged, width, bottom=echoed, color=HEDGE,
           label="Hedged; unsafe fields retained", zorder=3)
    bottoms = [a + b for a, b in zip(echoed, hedged)]
    ax.bar(x, reverted, width, bottom=bottoms, color=REVERT,
           label="Reverted to canonical CPIC answer", zorder=3)

    for i, rows in enumerate((forward, reverse)):
        ax.text(i, len(rows) + 0.7,
                f"{echoed[i] + hedged[i]}/{len(rows)} followed\n{reverted[i]} reverted",
                ha="center", va="bottom", fontsize=8.5, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([
        "Forward\n(lethal contract called safe)",
        "Reverse\n(safe contract called dangerous)",
    ])
    ax.set_ylabel("Responses (of 45)")
    ax.set_ylim(0, 55)
    ax.set_title("Models follow the supplied contract, including when it is corrupted",
                 loc="left", fontweight="bold")
    ax.grid(axis="y", color="#e7e7e7", linewidth=0.7, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              frameon=False, ncol=1)
    fig.tight_layout()

    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "Figure6_adversarial_contract.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / "Figure6_adversarial_contract.tiff", dpi=600,
                bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("wrote figures/Figure6_adversarial_contract.png and .tiff")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    result = self_test()
    if result or args.self_test:
        return result
    render()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
