"""Digitize the bar charts of mean + SD of a figure of a study.

The images are those of the PDF of the study (`pdfimages -png`), which the
repository does not hold. A configuration (JSON) gives per panel the bottom
of its y axis (the row of the baseline), the top of the axis line (`top`)
and its ticks `[row, value]`, read
from the labels of the figure; a value between the ticks is interpolated
piecewise linearly, so an uneven scan of the page is followed. The bars are
found on the row just above the baseline as pairs of vertical edges; the mean
is the top line of a bar, the SD the top of its whisker above it.

    uv run python scripts/database/digitize_barchart.py CONFIG.json IMAGE_DIR

prints the mean and SD of every bar of every panel in the order of the
`groups` of the configuration. `scripts/database/leithe1984` holds the
configuration of Leithe 1984 (Figure 1, page 3 of the PDF).
"""

import itertools
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image

DARK = 128
"""The gray value below which a pixel is dark."""


def value(row: float, ticks: list[list[float]]) -> float:
    """The value of a row, piecewise linear between the ticks, extrapolated outside."""
    points = sorted((float(r), float(v)) for r, v in ticks)
    for (r0, v0), (r1, v1) in itertools.pairwise(points):
        if row <= r1:
            return v0 + (v1 - v0) * (row - r0) / (r1 - r0)
    (r0, v0), (r1, v1) = points[-2], points[-1]
    return v0 + (v1 - v0) * (row - r0) / (r1 - r0)


def runs(column: np.ndarray) -> list[tuple[int, int]]:
    """The (start, end) of the runs of dark pixels of a vector."""
    d = np.diff(np.concatenate([[0], column.astype(int), [0]]))
    return [
        (int(s), int(e))
        for s, e in zip(np.where(d == 1)[0], np.where(d == -1)[0], strict=True)
    ]


def bars(img: np.ndarray, axis_x: int, bottom: int, top: int) -> list[dict[str, Any]]:
    """The top rows of the bars of a panel and of their whiskers."""
    columns = [
        x
        for x in range(axis_x + 8, axis_x + 800)
        if img[bottom - 40 : bottom - 8, x].mean() > 0.9
    ]
    edges: list[list[int]] = []
    for x in columns:
        if edges and x - edges[-1][-1] <= 3:
            edges[-1].append(x)
        else:
            edges.append([x])
    found = []
    for left, right in zip(edges[0::2], edges[1::2], strict=False):
        lo, hi = min(left), max(right)
        bar_top = next(
            y for y in range(top - 20, bottom) if img[y, lo + 3 : hi - 3].mean() > 0.8
        )
        whisker = None
        for x in range(lo, hi + 1):
            for s, e in runs(img[top - 40 : bar_top, x]):
                if e + top - 40 >= bar_top - 3 and (
                    whisker is None or e - s > whisker[1]
                ):
                    whisker = (s + top - 40, e - s)
        found.append(
            {"top": bar_top + 2, "whisker": None if whisker is None else whisker[0]}
        )
    return found


def main() -> None:
    """Print the mean and SD of every bar of a configuration."""
    config_path, image_dir = Path(sys.argv[1]), Path(sys.argv[2])
    config = json.loads(config_path.read_text(encoding="utf-8"))
    img = np.array(Image.open(image_dir / config["image"]).convert("L")) < DARK
    for name, panel in config["panels"].items():
        ticks = panel["ticks"]
        # the top of the axis line, above the highest tick a whisker may reach
        top = int(panel.get("top", min(r for r, _ in ticks)))
        for group, bar in zip(
            config["groups"],
            bars(img, panel["axis_x"], panel["bottom"], top),
            strict=False,
        ):
            mean = value(bar["top"], ticks)
            sd = None if bar["whisker"] is None else value(bar["whisker"], ticks) - mean
            print(f"{name}\t{group}\t{mean:.4g}\t{'' if sd is None else f'{sd:.3g}'}")


if __name__ == "__main__":
    main()
