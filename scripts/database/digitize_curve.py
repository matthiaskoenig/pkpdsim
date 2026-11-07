"""Digitize a solid curve of a panel of a figure of a study.

The images are those of the PDF of the study (`pdfimages -png`), which the
repository does not hold. A configuration (JSON) gives per panel the box of
the plot interior (`box`, the rows and columns of the image which hold the
curve), the calibration of the axes (`x` for a decadic logarithmic axis, the
column of its first decade and the columns per decade; `y` for a linear axis,
the row of zero and the rows per unit), the `grid` of the x values at which
the curve is read and an optional `close` point which the curve is taken to
reach at the right (the adult value of an ontogeny).

The curve is traced column by column: a pixel is on the curve when it is dark
and unsaturated (so a colored curve of the same panel is left out) and when
its run of dark pixels is no thicker than `max_thick` (so the frame and a
filled marker are left out). The run of each column which continues the curve
is chosen by dynamic programming over the columns, minimising the sum of the
steps between the runs of neighbouring columns, which follows the curve
through the data markers it crosses. The trace is smoothed over `smooth`
columns against the noise of the scan and read at the `grid`.

    uv run python scripts/database/digitize_curve.py CONFIG.json IMAGE

prints the points of every panel of the configuration as the `age` and
`value` arrays of a `points` entry of the database.
`scripts/database/upreti2016` holds the configuration of Upreti 2016
(Figure 1, page 5 of the PDF: the in vivo ontogeny of CYP1A2 in panel A and
of CYP3A in panel I).
"""

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image

DARK = 110
"""The value below which the brightest channel of a pixel is dark."""

SATURATION = 45
"""The spread of the channels above which a pixel is colored, not black."""


def _runs(column: np.ndarray, max_thick: int) -> list[float]:
    """The centres of the dark runs of a column which are no thicker than `max_thick`."""
    out: list[float] = []
    start: int | None = None
    for i, value in enumerate(column):
        if value and start is None:
            start = i
        elif not value and start is not None:
            if i - start <= max_thick:
                out.append((start + i - 1) / 2.0)
            start = None
    if start is not None and len(column) - start <= max_thick:
        out.append((start + len(column) - 1) / 2.0)
    return out


def trace(
    image: Path, box: list[int], max_thick: int = 14, jump: float = 26.0
) -> list[tuple[int, float]]:
    """Trace the solid dark curve of the box of an image.

    Args:
        image: the image of the figure.
        box: the plot interior, `[top, bottom, left, right]` in pixels.
        max_thick: the thickest run of dark pixels which is still the curve.
        jump: the step between two columns above which the curve is taken to
            be interrupted (a marker drawn over it), in pixels.

    Returns:
        The column and the row of the curve, left to right.
    """
    pixels = np.asarray(Image.open(image).convert("RGB")).astype(int)
    black = (pixels.max(axis=2) < DARK) & (
        pixels.max(axis=2) - pixels.min(axis=2) < SATURATION
    )
    top, bottom, left, right = box
    columns: list[tuple[int, dict[float, tuple[float, float | None]]]] = []
    best: dict[float, tuple[float, float | None]] = {}
    for x in range(left, right + 1):
        current: dict[float, tuple[float, float | None]] = {}
        for y in _runs(black[top:bottom, x], max_thick):
            if not best:
                current[y] = (0.0, None)
                continue
            cost, previous = min((best[q][0] + min(abs(y - q), jump), q) for q in best)
            current[y] = (cost, previous)
        if current:
            best = current
        columns.append((x, dict(current)))
    last = next(i for i in range(len(columns) - 1, -1, -1) if columns[i][1])
    y = min(columns[last][1], key=lambda q: columns[last][1][q][0])
    points: list[tuple[int, float]] = []
    for i in range(last, -1, -1):
        if not columns[i][1]:
            continue
        points.append((columns[i][0], y + top))
        nxt = columns[i][1][y][1]
        if nxt is None:
            break
        y = nxt
    return list(reversed(points))


def digitize(image: Path, panel: dict[str, Any]) -> list[tuple[float, float]]:
    """The points of a panel at the grid of its configuration."""
    points = trace(image, panel["box"], panel.get("max_thick", 14))
    x0, per_decade = panel["x"]
    y0, per_unit = panel["y"]
    ages = np.array([10.0 ** (panel["x0"] + (x - x0) / per_decade) for x, _y in points])
    values = np.array([(y0 - y) / per_unit for _x, y in points])
    keep = ages <= panel["cut"]
    ages, values = ages[keep], values[keep]
    k = int(panel.get("smooth", 9))
    smooth = np.convolve(values, np.ones(k) / k, mode="same")
    smooth[: k // 2] = values[: k // 2]
    smooth[-(k // 2) :] = values[-(k // 2) :]
    out = [(float(g), float(np.interp(g, ages, smooth))) for g in panel["grid"]]
    if "close" in panel:
        out.append((float(panel["close"][0]), float(panel["close"][1])))
    return out


def main(config: Path, image: Path) -> None:
    """Print the points of every panel of a configuration."""
    panels = json.loads(config.read_text(encoding="utf-8"))["panels"]
    for name, panel in panels.items():
        points = digitize(image, panel)
        ages = ", ".join(f"{a:g}" for a, _v in points)
        values = ", ".join(f"{v:.3f}" for _a, v in points)
        print(f"# {name}: {panel.get('what', '')}")
        print(f"age = [{ages}]")
        print(f"value = [{values}]")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
