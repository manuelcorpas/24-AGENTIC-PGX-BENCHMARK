#!/usr/bin/env python3
"""Draw Figure 8 from the frozen seven-model normalisation analysis.

Both panels use the same paired variant-call units within each model. Provider
errors and output-budget truncations remain visible in the operational counts
but leave behavioural denominators. Error bars are percentile 95% intervals
from the frozen sample-cluster bootstrap. No model calls are issued here.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
FIGDIR = BASE / "figures"
FREEZE = DATA / "v3_input_normalisation_seven_model_freeze.json"

WITHOUT = "#2a78d6"
WITH = "#eb6834"
CALLER = "#5c5c5c"
INK = "#111111"
MODEL_LABELS = {
    "Claude Opus 4.5": "Claude\nOpus 4.5",
    "Claude Sonnet 4.5": "Claude\nSonnet 4.5",
    "GPT-5.2": "GPT-5.2",
    "GPT-4.1": "GPT-4.1",
    "o3": "o3",
    "o4-mini": "o4-mini",
    "DeepSeek V3": "DeepSeek\nV3",
}
MIN_CALLS_FOR_ACCURACY = 10


def load_freeze() -> dict:
    if not FREEZE.exists():
        raise FileNotFoundError(
            f"missing {FREEZE}; run code/83-freeze-seven-model-normalisation.py"
        )
    return json.loads(FREEZE.read_text())


FONT_REGULAR = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"


def font(size: int, bold: bool = False):
    return ImageFont.truetype(FONT_BOLD if bold else FONT_REGULAR, size)


def centered(draw, xy, text, fnt, fill=INK, anchor="mm"):
    draw.text(xy, text, font=fnt, fill=fill, anchor=anchor, align="center")


def marker(draw, x, y, colour, filled=True, diamond=False, radius=17):
    if diamond:
        pts = [(x, y-radius), (x+radius, y), (x, y+radius), (x-radius, y)]
        draw.polygon(pts, fill="white" if not filled else colour,
                     outline=colour, width=5)
    else:
        draw.ellipse((x-radius, y-radius, x+radius, y+radius),
                     fill=colour if filled else "white", outline=colour, width=5)


def errorbar(draw, x, y, lo, hi, colour, ymap):
    ylo, yhi = ymap(lo), ymap(hi)
    draw.line((x, ylo, x, yhi), fill=colour, width=5)
    draw.line((x-12, ylo, x+12, ylo), fill=colour, width=5)
    draw.line((x-12, yhi, x+12, yhi), fill=colour, width=5)


def draw(data: dict, stem: str = "Figure8-input-normalisation") -> None:
    models = data["models"]
    analysis = data["models_analysis"]
    # Draw at 600 dpi and downsample the PNG to 300 dpi. This keeps text and
    # confidence intervals crisp without relying on a GUI plotting backend.
    W, H = 6600, 3000
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    centered(d, (W//2, 120),
             "Validated allele definitions improve every tested model, but model choice still matters",
             font(70, True))

    panels = [(420, 390, 3150, 2250), (3510, 390, 6240, 2250)]
    titles = [("A  Coverage", "Fraction emitting a diplotype"),
              ("B  Accuracy against GeT-RM", "Accuracy among emitted calls")]

    for (left, top, right, bottom), (title, ylabel) in zip(panels, titles):
        d.text((left, top-100), title, font=font(54, True), fill=INK, anchor="la")
        for tick in [0, .2, .4, .6, .8, 1.0]:
            y = bottom - int((tick / 1.09) * (bottom-top))
            d.line((left, y, right, y), fill="#e4e4e4", width=3)
            d.text((left-35, y), f"{tick:.1f}", font=font(34), fill="#333333",
                   anchor="rm")
        d.line((left, top, left, bottom), fill="#555555", width=5)
        d.line((left, bottom, right, bottom), fill="#555555", width=5)
        # Rotated y label.
        layer = Image.new("RGBA", (900, 80), (255, 255, 255, 0))
        ld = ImageDraw.Draw(layer)
        centered(ld, (450, 40), ylabel, font(38), fill="#222222")
        layer = layer.rotate(90, expand=True)
        img.paste(layer, (left-310, (top+bottom-layer.height)//2), layer)

    def coordinates(panel, i, offset, value):
        left, top, right, bottom = panel
        step = (right-left) / len(models)
        x = left + step * (i + .5 + offset)
        y = bottom - (value / 1.09) * (bottom-top)
        return int(x), int(y)

    for i, model in enumerate(models):
        cell = analysis[model]
        for side, key, colour, offset in (
            ("without", "without_definitions", WITHOUT, -0.16),
            ("with", "with_definitions", WITH, 0.10),
        ):
            arm = cell[key]
            metric = arm["coverage"]
            x, y = coordinates(panels[0], i, offset, metric["estimate"])
            ymap = lambda v, p=panels[0]: coordinates(p, i, offset, v)[1]
            errorbar(d, x, y, metric["ci95"][0], metric["ci95"][1], colour, ymap)
            marker(d, x, y, colour, filled=side == "with")
            centered(d, (x, y-52), f"{metric['estimate']:.2f}", font(29), anchor="ms")

            calls = arm["operational"]["call"]
            if calls < MIN_CALLS_FOR_ACCURACY:
                x2, _ = coordinates(panels[1], i, offset, 0)
                centered(d, (x2, panels[1][3]-45), f"n={calls}", font(29),
                         fill="#555555", anchor="ms")
            else:
                metric = arm["accuracy_vs_getrm_among_calls"]
                x2, y2 = coordinates(panels[1], i, offset, metric["estimate"])
                ymap = lambda v, p=panels[1]: coordinates(p, i, offset, v)[1]
                errorbar(d, x2, y2, metric["ci95"][0], metric["ci95"][1], colour, ymap)
                marker(d, x2, y2, colour, filled=side == "with")
                if side == "with":
                    centered(d, (x2, y2-52), f"{metric['estimate']:.2f}",
                             font(29), anchor="ms")

        caller = cell["with_definitions"]["caller_vs_getrm_on_model_calls"]
        x, y = coordinates(panels[1], i, 0.31, caller["estimate"])
        ymap = lambda v, p=panels[1]: coordinates(p, i, 0.31, v)[1]
        errorbar(d, x, y, caller["ci95"][0], caller["ci95"][1], CALLER, ymap)
        marker(d, x, y, CALLER, filled=False, diamond=True, radius=15)

        for panel in panels:
            xlab, _ = coordinates(panel, i, 0, 0)
            label = MODEL_LABELS[model]
            centered(d, (xlab, panel[3]+95), label, font(33), anchor="ma")

    # Legend.
    legend_y = 2540
    items = [(WITHOUT, False, False, "without allele definitions"),
             (WITH, True, False, "with allele definitions"),
             (CALLER, False, True, "deterministic caller on model-emitted pairs")]
    starts = [720, 2450, 4030]
    for start, (colour, filled, diamond, label) in zip(starts, items):
        marker(d, start, legend_y, colour, filled=filled, diamond=diamond, radius=16)
        d.text((start+38, legend_y), label, font=font(34), fill="#222222", anchor="lm")

    centered(
        d, (W//2, 2790),
        "Paired variant-call inputs: n=527 per model except o3 (prespecified n=150 subsample).",
        font(31), fill="#444444"
    )
    centered(
        d, (W//2, 2850),
        "Error bars: sample-cluster bootstrap 95% intervals. Accuracy is omitted when fewer than 10 calls were emitted.",
        font(31), fill="#444444"
    )
    FIGDIR.mkdir(exist_ok=True)
    img.save(FIGDIR / f"{stem}.tiff", dpi=(600, 600), compression="tiff_lzw")
    png = img.resize((W//2, H//2), Image.Resampling.LANCZOS)
    png.save(FIGDIR / f"{stem}.png", dpi=(300, 300), optimize=True)
    print(f"wrote {FIGDIR / stem}.png and .tiff")


def self_test(data: dict) -> None:
    assert len(data["models"]) == 7
    assert data["definition_arm_operational_total"] == {
        "attempted": 3689, "call": 2905, "abstain": 780,
        "error": 0, "truncated_output": 4, "scored": 3685,
    }
    opus = data["models_analysis"]["Claude Opus 4.5"]["with_definitions_full_527"]
    assert round(opus["accuracy_vs_pypgx_among_calls"]["estimate"], 3) == 0.967
    diff = opus["model_minus_caller_vs_getrm"]
    assert diff["ci95"][0] < 0 < diff["ci95"][1]
    assert data["models_analysis"]["o3"]["paired_units"] == 150
    print("self-test OK")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    try:
        data = load_freeze()
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 2
    self_test(data)
    if not args.self_test:
        draw(data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
