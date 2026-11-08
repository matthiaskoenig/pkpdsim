"""The deformation of a human drawing into the body of an individual.

The outline of the drawing is rasterized and segmented into head, trunk, arms
and legs (`Drawing`). The body of an individual is the drawing mapped by a
smooth field (`Field`): a monotone vertical map through the landmarks gives
the proportions of the age (the head, the trunk, the legs), and every part is
stretched horizontally about its own axis by a sum of broad Gaussians over the
height, blended between the parts by blurred masks, so the outline changes
continuously and mainly by stretching in x and y. The field maps points of the
drawing (`Field.__call__`) and is inverted by a fixed point iteration
(`Field.inverse`), which the raster of the body needs.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from functools import cache

import numpy as np
from PIL import Image, ImageDraw
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter, map_coordinates

from pkpdsim.pbpk.visual.resources import Anatomy, anatomy

__all__ = ["RES", "Drawing", "Field", "Gaussian", "drawing", "raster"]

RES = 0.25
"""The raster of the drawing: units of the drawing per pixel."""

Gaussian = tuple[float, float, float]
"""A bump of a profile over the height: centre, sigma (units of the drawing) and amplitude."""


def raster(
    polygons: Sequence[np.ndarray],
    shape: tuple[int, int],
    origin: tuple[float, float] = (0.0, 0.0),
    res: float = RES,
) -> np.ndarray:
    """The pixels inside polygons.

    Args:
        polygons: closed polylines, arrays of shape (n, 2).
        shape: the shape of the raster (rows, columns).
        origin: the point of the pixel (0, 0).
        res: the units per pixel.

    Returns:
        A boolean array of the shape.
    """
    im = Image.new("L", (shape[1], shape[0]), 0)
    draw = ImageDraw.Draw(im)
    for p in polygons:
        if len(p) > 2:
            draw.polygon(
                [((x - origin[0]) / res, (y - origin[1]) / res) for x, y in p], fill=1
            )
    return np.asarray(im) > 0


def _runs(row: np.ndarray) -> list[tuple[float, float]]:
    """The runs of inside pixels of a row, as intervals of the drawing."""
    out, i, n = [], 0, len(row)
    while i < n:
        if row[i]:
            j = i
            while j < n and row[j]:
                j += 1
            out.append((i * RES, j * RES))
            i = j
        else:
            i += 1
    return out


def _filled(values: np.ndarray, sigma: float = 6.0) -> np.ndarray:
    """A profile defined on some rows, interpolated over all rows and smoothed."""
    ok = ~np.isnan(values)
    if not ok.any():
        return np.zeros_like(values)
    rows = np.arange(len(values))
    return gaussian_filter(np.interp(rows, rows[ok], values[ok]), sigma)


@dataclass(frozen=True, eq=False)
class Drawing:
    """The parts of a human drawing on its raster.

    Attributes:
        anatomy: the resource.
        mask: the inside of the outline.
        ys: the y of every row.
        xs: the x of every column.
        cx: the midline.
        labels: per pixel 1 head, 2 trunk, 3 arm, 4 leg, 0 outside.
        arm_axis: the x of the axis of the arm per row, by side (-1, 1).
        leg_axis: the x of the axis of the leg per row, by side.
        half: the half width of the trunk per row.
        weights: the blurred masks of the parts, summing to 1.
    """

    anatomy: Anatomy
    mask: np.ndarray
    ys: np.ndarray
    xs: np.ndarray
    cx: float
    labels: np.ndarray
    arm_axis: dict[int, np.ndarray]
    leg_axis: dict[int, np.ndarray]
    half: np.ndarray
    weights: dict[str, np.ndarray]

    @property
    def landmarks(self) -> dict[str, float]:
        """The landmarks of the drawing."""
        return dict(self.anatomy.landmarks)


@cache
def drawing(sex: str) -> Drawing:
    """The segmented human drawing of a sex (cached).

    Args:
        sex: `male` or `female`.

    Returns:
        The drawing with its parts.
    """
    a = anatomy("human", sex)
    lm = a.landmarks
    x0, y0, w, h = a.viewbox
    shape = (int((y0 + h + 2) / RES), int((x0 + w + 2) / RES))
    mask = raster([a.outline], shape)
    rows, cols = mask.shape
    ys, xs = np.arange(rows) * RES, np.arange(cols) * RES
    cx = a.midline
    labels = np.zeros(mask.shape, np.int8)
    arm = {-1: np.full(rows, np.nan), 1: np.full(rows, np.nan)}
    leg = {-1: np.full(rows, np.nan), 1: np.full(rows, np.nan)}
    half = np.full(rows, np.nan)
    hand_end = lm["crotch"] + 0.1 * (lm["bottom"] - lm["crotch"])
    for i, y in enumerate(ys):
        runs = _runs(mask[i])
        if not runs:
            continue
        if y < lm["neck"]:
            for a0, b0 in runs:
                labels[i, int(a0 / RES) : int(b0 / RES)] = 1
            continue
        trunk = next(((a0, b0) for a0, b0 in runs if a0 - 0.5 <= cx <= b0 + 0.5), None)
        if y <= lm["crotch"]:
            for a0, b0 in runs:
                if (a0, b0) == trunk:
                    labels[i, int(a0 / RES) : int(b0 / RES)] = 2
                    half[i] = (b0 - a0) / 2
                else:
                    labels[i, int(a0 / RES) : int(b0 / RES)] = 3
                    arm[1 if (a0 + b0) / 2 > cx else -1][i] = (a0 + b0) / 2
            continue
        legs = [
            (a0, b0) for a0, b0 in runs if a0 < cx < b0 or abs((a0 + b0) / 2 - cx) < 22
        ]
        for a0, b0 in runs:
            if (a0, b0) in legs:
                labels[i, int(a0 / RES) : int(b0 / RES)] = 4
                for s, lo, hi in ((-1, a0, min(b0, cx)), (1, max(a0, cx), b0)):
                    if hi > lo + 0.5 and np.isnan(leg[s][i]):
                        leg[s][i] = (lo + hi) / 2
            elif y < hand_end:
                labels[i, int(a0 / RES) : int(b0 / RES)] = 3
                arm[1 if (a0 + b0) / 2 > cx else -1][i] = (a0 + b0) / 2
    blurred = [gaussian_filter((labels == k).astype(float), 6) for k in (1, 2, 3, 4)]
    total = sum(blurred) + 1e-9
    return Drawing(
        anatomy=a,
        mask=mask,
        ys=ys,
        xs=xs,
        cx=cx,
        labels=labels,
        arm_axis={s: _filled(v) for s, v in arm.items()},
        leg_axis={s: _filled(v) for s, v in leg.items()},
        half=_filled(half),
        weights={
            k: b / total
            for k, b in zip(("head", "trunk", "arm", "leg"), blurred, strict=True)
        },
    )


def _profile(ys: np.ndarray, gaussians: Sequence[Gaussian]) -> np.ndarray:
    p = np.zeros(len(ys))
    for c, s, amplitude in gaussians:
        p += amplitude * np.exp(-0.5 * ((ys - c) / s) ** 2)
    return p


class Field:
    """The map of a drawing onto the body of an individual (units of the drawing to cm).

    Args:
        d: the drawing.
        height: the height of the individual [cm].
        head_frac: the head (top to neck) as a fraction of the height.
        leg_frac: the legs (crotch to floor) as a fraction of the height.
        trunk: the relative widening of the trunk over the height.
        arm: the relative widening of the arms about their axes.
        leg: the relative widening of the legs about their axes.
        head_w: the width of the head relative to its height.
        belly: a pregnant belly, a Gaussian of the trunk (centre, sigma,
            amplitude).
    """

    def __init__(
        self,
        d: Drawing,
        height: float,
        head_frac: float,
        leg_frac: float,
        trunk: Sequence[Gaussian] = (),
        arm: Sequence[Gaussian] = (),
        leg: Sequence[Gaussian] = (),
        head_w: float = 1.0,
        belly: Gaussian | None = None,
    ) -> None:
        """Build the field of an individual."""
        self.d = d
        lm = d.landmarks
        ys = d.ys
        src = [lm["top"], lm["neck"], lm["crotch"], lm["bottom"]]
        dst = [0.0, head_frac * height, height * (1 - leg_frac), height]
        self.Y = PchipInterpolator(src, dst)
        self.kh = (dst[1] - dst[0]) / (src[1] - src[0])
        self.kt = (dst[2] - dst[1]) / (src[2] - src[1])
        self.kx = gaussian_filter(
            np.where(ys < lm["neck"], self.kh * head_w, self.kt), 12
        )
        st, sa, sl = _profile(ys, trunk), _profile(ys, arm), _profile(ys, leg)
        x = d.xs[None, :] - d.cx
        w = d.weights
        dx = w["trunk"] * x * st[:, None]
        edge = gaussian_filter(d.half * st, 4)
        side = np.sign(x) * np.ones((len(ys), 1))
        crotch = lm["crotch"]
        hip_profile = np.interp(
            ys,
            [crotch, crotch + 0.35 * (lm["bottom"] - crotch), lm["bottom"]],
            [edge[int(crotch / RES)] * 0.8, 0.0, 0.0],
        )
        hip = gaussian_filter(np.where(ys > crotch, hip_profile, 0.0), 4)
        for s in (-1, 1):
            m = side == s
            arm_x = d.xs[None, :] - d.arm_axis[s][:, None]
            dx += np.where(m, w["arm"] * (s * edge[:, None] + arm_x * sa[:, None]), 0)
            leg_x = d.xs[None, :] - d.leg_axis[s][:, None]
            dx += np.where(m, w["leg"] * (s * hip[:, None] + leg_x * sl[:, None]), 0)
        if belly is not None:
            dx += w["trunk"] * x * _profile(ys, [belly])[:, None]
        self.dx = gaussian_filter(dx, 2)
        self.dy = np.zeros_like(self.dx)

    def __call__(self, points: np.ndarray) -> np.ndarray:
        """Map points of the drawing into the body [cm], the midline at x 0, the top at y 0."""
        p = np.asarray(points, float)
        r, c = p[:, 1] / RES, p[:, 0] / RES
        dx = map_coordinates(self.dx, [r, c], order=1, mode="nearest")
        dy = map_coordinates(self.dy, [r, c], order=1, mode="nearest")
        y = p[:, 1] + dy
        kx = np.interp(y, self.d.ys, self.kx)
        return np.c_[(p[:, 0] + dx - self.d.cx) * kx, self.Y(y)]

    def inverse(
        self, x: np.ndarray, y: np.ndarray, iterations: int = 8
    ) -> tuple[np.ndarray, np.ndarray]:
        """The points of the drawing which the field maps onto (x, y) [cm]."""
        ys = self.d.ys
        yp = np.interp(y, self.Y(ys), ys)
        u = x / np.interp(yp, ys, self.kx) + self.d.cx
        px, py = u.copy(), yp.copy()
        for _ in range(iterations):
            r, c = py / RES, px / RES
            px = u - map_coordinates(self.dx, [r, c], order=1, mode="nearest")
            py = yp - map_coordinates(self.dy, [r, c], order=1, mode="nearest")
        return px, py

    @property
    def height(self) -> float:
        """The height of the body [cm]."""
        return float(self.Y(self.d.landmarks["bottom"]))
