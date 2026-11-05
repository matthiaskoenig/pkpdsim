"""Digitize the box plots of R (`boxplot` with the mean as "+") of a figure of a study.

The images are those of the PDF of the study (`pdfimages -png PDF fig`), which
the repository does not hold. A configuration (JSON) gives per panel its y
axis (the x of the axis line) and its axis segments (a broken axis has two),
each calibrated by two ticks of known value, read from the labels of the
figure; the group centers are the ticks of the x axis below the panel (or
those of the panel above). Per group the rows through the center are
classified by the length of their horizontal dark run: the box edges and the
median span the box (the median thicker), the whisker caps half the box, the
horizontal stroke of the mean "+" less. A cap hidden by a box edge (no
whisker) is the edge; a "+" hidden by the median or an edge has no mean.

    uv run python scripts/database/digitize_boxplot.py CONFIG.json IMAGE_DIR

writes `<config>.out.json` next to the configuration (the values of every
group of every panel and the raw rows). `scripts/database/drozdzik2020` holds
the configurations of Drozdzik 2020 (Figures 1, 2, 5 and 6 of the accepted
manuscript, images `fig-000.png`, `fig-001.png`, `fig-004.png` and
`fig-005.png` of `pdfimages`).
"""

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image

DARK = 200
"""The gray value below which a pixel is dark (the lines are anti-aliased)."""


@dataclass(frozen=True)
class Segment:
    """A segment of a y axis with two calibration ticks and its extent."""

    y_top: float
    v_top: float
    y_bottom: float
    v_bottom: float
    y_min: float
    y_max: float

    def value(self, y: float) -> float:
        """The value of a row of the segment."""
        scale = (self.v_top - self.v_bottom) / (self.y_bottom - self.y_top)
        return self.v_bottom + (self.y_bottom - y) * scale


def load(path: Path) -> np.ndarray:
    """The dark pixels of an image."""
    return np.array(Image.open(path).convert("L")) < DARK


def runs(arr: np.ndarray) -> list[tuple[int, int]]:
    """The (start, end) of the runs of true values of a vector."""
    d = np.diff(np.concatenate([[0], arr.astype(int), [0]]))
    return [
        (int(s), int(e))
        for s, e in zip(np.where(d == 1)[0], np.where(d == -1)[0], strict=True)
    ]


def axis_ticks(img: np.ndarray, x: int, y0: float, y1: float) -> list[float]:
    """The rows of the ticks of a y axis at x (marks to its left)."""
    rows = [
        y
        for y in range(int(y0), int(y1) + 1)
        if img[y, x - 6 : x].all() and not img[y, x - 14 : x - 9].any()
    ]
    merged: list[list[int]] = []
    for y in rows:
        if merged and y - merged[-1][-1] <= 2:
            merged[-1].append(y)
        else:
            merged.append([y])
    return [float(np.mean(m)) for m in merged]


def x_ticks(
    img: np.ndarray, x_axis: int, y_bottom: float, x_end: int
) -> tuple[int, list[float]]:
    """The row of the x axis under a panel and its ticks, the group centers."""
    for y in range(int(y_bottom) + 3, int(y_bottom) + 60):
        long = [(s, e) for s, e in runs(img[y, x_axis:x_end]) if e - s > 150]
        if not long:
            continue
        s, e = long[0]
        ticks = [
            x_axis + x
            for x in range(max(s - 2, 0), e + 2)
            if img[y + 2 : y + 8, x_axis + x].all()
        ]
        merged: list[list[int]] = []
        for x in ticks:
            if merged and x - merged[-1][-1] <= 2:
                merged[-1].append(x)
            else:
                merged.append([x])
        if len(merged) >= 2:
            return y, [float(np.mean(m)) for m in merged]
    raise ValueError(f"no x axis under the panel at x {x_axis}, y {y_bottom}")


def run_through(img: np.ndarray, y: int, cx: int, reach: int = 3) -> int:
    """The length of the horizontal dark run through (y, cx), or the nearest within reach."""
    row = img[y]
    for dx in range(reach + 1):
        for x in (cx - dx, cx + dx):
            if row[x]:
                left = x
                while left > 0 and row[left - 1]:
                    left -= 1
                right = x
                while right < len(row) - 1 and row[right + 1]:
                    right += 1
                return right - left + 1
    return 0


def line_segments(
    img: np.ndarray, cx: float, y0: float, y1: float
) -> list[tuple[float, int, int]]:
    """The horizontal lines through a center: (mean row, thickness, length)."""
    out: list[list[tuple[int, int]]] = []
    current: list[tuple[int, int]] = []
    for y in range(int(y0), int(y1) + 1):
        n = run_through(img, y, round(cx))
        if n >= 8:
            if current and (y - current[-1][0] > 1 or abs(n - current[-1][1]) > 6):
                out.append(current)
                current = []
            current.append((y, n))
        elif current:
            out.append(current)
            current = []
    if current:
        out.append(current)
    return [
        (float(np.mean([r[0] for r in s])), len(s), max(r[1] for r in s)) for s in out
    ]


def digitize_group(
    img: np.ndarray, cx: float, box_width: float, y0: float, y1: float
) -> dict[str, Any]:
    """The rows of the median, the quartiles, the whisker ends and the mean of a group."""
    lines = line_segments(img, cx, y0, y1)
    wide = [s for s in lines if s[2] >= 0.7 * box_width]
    caps = [s for s in lines if 0.3 * box_width <= s[2] < 0.7 * box_width]
    plus = [s for s in lines if 10 <= s[2] < 0.3 * box_width]
    result: dict[str, Any] = {"notes": [], "lines": lines}
    thick = [s for s in wide if s[1] >= 4]
    edges = sorted(s[0] for s in wide if s[1] < 4)
    if thick:
        result["median"] = max(thick, key=lambda s: s[1])[0]
        edges = sorted([*edges, result["median"]])
    if edges:
        result["q3"], result["q1"] = edges[0], edges[-1]
    q3, q1 = result.get("q3"), result.get("q1")
    above = [c[0] for c in caps if q3 is not None and c[0] < q3 - 2]
    below = [c[0] for c in caps if q1 is not None and c[0] > q1 + 2]
    result["max"] = min(above) if above else q3
    result["min"] = max(below) if below else q1
    lo, hi = result["max"], result["min"]
    inside = [
        p[0]
        for p in plus
        if lo is not None
        and hi is not None
        and lo - 3 <= p[0] <= hi + 3
        and not any(abs(p[0] - c[0]) <= 3 for c in caps)
    ]
    if len(inside) == 1:
        result["mean"] = inside[0]
    elif len(inside) > 1:
        result["notes"].append(f"{len(inside)} candidates for the mean")
    else:
        result["notes"].append(
            "the mean is not separable from the median or a box edge"
        )
    return result


def digitize(config: dict[str, Any], image_dir: Path) -> list[dict[str, Any]]:
    """Digitize every panel of a configuration."""
    img = load(image_dir / config["image"])
    out = []
    known: dict[int, list[float]] = {}
    for panel in config["panels"]:
        segments = [Segment(**s) for s in panel["segments"]]
        y_bottom = max(s.y_max for s in segments)
        y_top = min(s.y_min for s in segments)
        try:
            y_axis_line, centers = x_ticks(
                img, panel["x_axis"], y_bottom, panel["x_end"]
            )
        except ValueError:
            y_axis_line, centers = round(y_bottom) + 15, []
        if len(centers) != len(panel["groups"]):
            # the groups of a panel are at the centers of the panel above
            centers = known[panel["x_axis"] // 10]
            y_axis_line = round(y_bottom) + 15
        known[panel["x_axis"] // 10] = centers
        box_width = 0.8 * (
            float(np.median(np.diff(centers))) if len(centers) > 1 else 100.0
        )
        lower = segments[-1]
        pixel = abs((lower.v_top - lower.v_bottom) / (lower.y_bottom - lower.y_top))
        for group, cx in zip(panel["groups"], centers, strict=True):
            res = digitize_group(img, cx, box_width, y_top - 5, y_axis_line - 2)
            row = {
                "panel": panel["name"],
                "entity": panel["entity"],
                "kind": panel["kind"],
                "group": group,
                "pixel": pixel,
                "notes": res["notes"],
                "lines": res["lines"],
            }
            for key in ("median", "q1", "q3", "min", "max", "mean"):
                y = res.get(key)
                row[key] = None if y is None else _value(segments, y)
            out.append(row)
    return out


def _value(segments: list[Segment], y: float) -> float | None:
    """The value of a row, `None` between the segments of a broken axis."""
    for s in segments:
        if s.y_min - 3 <= y <= s.y_max + 3:
            return round(s.value(y), 4)
    return None


def main() -> None:
    """Digitize a configuration into `<config>.out.json`."""
    config_path, image_dir = Path(sys.argv[1]), Path(sys.argv[2])
    config = json.loads(config_path.read_text(encoding="utf-8"))
    result = digitize(config, image_dir)
    out = config_path.with_suffix(".out.json")
    out.write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(f"written {out}")


if __name__ == "__main__":
    main()
