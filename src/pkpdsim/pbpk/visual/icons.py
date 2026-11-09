"""The icons of the figures, one graphical language for all of them.

Every icon sits on a white circle: the organs (their shapes in the anatomy
resources; the bowel for the intestine and its lumen, the body for the skin,
fat cells for the adipose tissue, whose drawn depot has no recognisable form at
this size), the oral dose (a tablet), the intravenous dose (a syringe) and the
blood sample (a drop). An icon is an SVG fragment centred at a point. The
pictogram of an impaired organ (`impaired_organ`) is the shape of the liver, a
kidney or the heart in the color of an impairment class, its lost part pale.
"""

import html
import math
import re
from functools import cache

import numpy as np

from pkpdsim.pbpk.visual import palette
from pkpdsim.pbpk.visual.resources import anatomy

__all__ = ["drop", "impaired_organ", "organ", "organ_keys", "ring", "syringe", "tablet"]

_FAT_CELLS = (
    (3.2, 3.4, 2.6),
    (8.4, 3.0, 2.4),
    (5.8, 7.6, 2.7),
    (10.6, 7.8, 2.2),
    (2.6, 8.6, 1.9),
)


def _circle(x: float, y: float, r: float) -> str:
    steps = [i * math.pi / 12 for i in range(24)]
    return (
        "M"
        + " ".join(
            f"{x + r * math.cos(t):.2f},{y + r * math.sin(t):.2f}" for t in steps
        )
        + "Z"
    )


@cache
def _icons(sex: str) -> dict[str, tuple[str, tuple[tuple[str, str], ...], float]]:
    """The organ icons of a sex: the view box, the coloured paths and the extent."""
    a = anatomy("human", sex)
    parts: dict[str, tuple[tuple[str, str], ...]] = {
        k: ((k, s.d),) for k, s in a.organs.items()
    }
    outline = "M" + " ".join(f"{x:.2f},{y:.2f}" for x, y in a.outline) + "Z"
    parts["sk"] = (("sk", outline),)
    bowel = (("la", a.organs["la"].d), ("si", a.organs["si"].d))
    parts["gu"] = bowel
    parts["lumen"] = bowel
    parts["ad"] = (("ad", "".join(_circle(*c) for c in _FAT_CELLS)),)
    out = {}
    for k, items in parts.items():
        nums = [float(v) for _, d in items for v in re.findall(r"-?\d+\.?\d*", d)]
        xs, ys = nums[0::2], nums[1::2]
        x, y = min(xs), min(ys)
        w, h = max(xs) - x, max(ys) - y
        m = max(w, h) * 0.04
        out[k] = (
            f"{x - m:.2f} {y - m:.2f} {w + 2 * m:.2f} {h + 2 * m:.2f}",
            items,
            max(w, h),
        )
    return out


def organ_keys(sex: str = "male") -> frozenset[str]:
    """The keys which have an organ icon: the organs, `gu` and `lumen`."""
    return frozenset(_icons(sex))


def ring(cx: float, cy: float, r: float, color: str = palette.RING) -> str:
    """The white circle of an icon."""
    return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="#fff" stroke="{color}" stroke-width="1.3"/>'


def organ(key: str, cx: float, cy: float, d: float = 24.0, sex: str = "male") -> str:
    """The icon of an organ on a circle of diameter `d`, centred at (cx, cy); empty without an icon.

    The icon is in the greys of `palette.ORGAN_GREYS`, so that the colours of a
    figure are those of the substances and the tissues.
    """
    icons = _icons(sex)
    if key not in icons:
        return ""
    box, items, extent = icons[key]
    s = d * 0.62
    width = extent / s * 0.9
    paths = "".join(
        f'<path d="{p}" fill="{palette.ORGAN_GREYS.get(k, "#ccc")}" stroke="{palette.ICON_INK}" '
        f'stroke-width="{width:.3f}" fill-rule="evenodd"/>'
        for k, p in items
    )
    return (
        ring(cx, cy, d / 2)
        + f'<svg x="{cx - s / 2:.1f}" y="{cy - s / 2:.1f}" width="{s:.1f}" height="{s:.1f}" '
        f'viewBox="{box}" preserveAspectRatio="xMidYMid meet">{paths}</svg>'
    )


_PICTOGRAMS: dict[str, tuple[int, str]] = {
    "li": (0, "right"),
    "ki": (0, "top"),
    "he": (0, "top"),
}
"""The organs of a pictogram of an impairment: the subpath of the drawing (one
kidney of the two) and the side from which the lost part grows (the liver from
the right as its left lobe, the kidney and the heart from the top)."""


def _clip(poly: np.ndarray, axis: int, cut: float, keep_above: bool) -> np.ndarray:
    """The part of a polygon on one side of the line `x[axis] = cut` (Sutherland-Hodgman)."""
    out: list[np.ndarray] = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        p_in = (p[axis] >= cut) if keep_above else (p[axis] <= cut)
        q_in = (q[axis] >= cut) if keep_above else (q[axis] <= cut)
        if p_in:
            out.append(p)
        if p_in != q_in:
            t = (cut - p[axis]) / (q[axis] - p[axis])
            out.append(p + t * (q - p))
    return np.array(out) if out else np.empty((0, 2))


def _area(poly: np.ndarray) -> float:
    """The area of a polygon (shoelace)."""
    if len(poly) < 3:
        return 0.0
    x, y = poly[:, 0], poly[:, 1]
    return abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))) / 2


def _split(poly: np.ndarray, side: str, loss: float) -> tuple[np.ndarray, np.ndarray]:
    """The kept and the lost part of a polygon, the lost part `loss` of its area.

    The lost part grows from `side` (`right` or `top`); the line is found by
    bisection on the area.
    """
    axis = 0 if side == "right" else 1
    lost_above = side == "right"  # the right is the large x, the top the small y
    lo, hi = float(poly[:, axis].min()), float(poly[:, axis].max())
    total = _area(poly)
    a, b = lo, hi
    for _ in range(60):
        cut = (a + b) / 2
        lost = _area(_clip(poly, axis, cut, lost_above)) / total
        # moving the cut towards the lost side shrinks the lost part
        if (lost > loss) == lost_above:
            a = cut
        else:
            b = cut
    cut = (a + b) / 2
    return _clip(poly, axis, cut, not lost_above), _clip(poly, axis, cut, lost_above)


def impaired_organ(
    key: str,
    cx: float,
    cy: float,
    width: float,
    height: float,
    loss: float,
    color: str,
    *,
    title: str | None = None,
    sex: str = "male",
) -> str:
    """The pictogram of an impaired organ: its shape, the lost part pale.

    The shape of the organ (`li`, `ki` or `he`, the anatomy resources) fills the
    box `width` x `height` centred at (cx, cy); the functional part is filled with
    `color`, the lost part, `loss` of the area, with `color` at a low opacity, and
    the outline is `currentColor`, so the pictogram reads in a light and a dark
    theme.

    Args:
        key: the organ.
        cx: the centre of the box, x.
        cy: the centre of the box, y.
        width: the width of the box.
        height: the height of the box.
        loss: the lost fraction of the organ, 0 to 1.
        color: the color of the class.
        title: the hover text, plain text (escaped here).
        sex: the sex of the drawing.

    Returns:
        The SVG fragment, a group.

    Raises:
        ValueError: if the organ has no pictogram or `loss` is outside [0, 1].
    """
    if key not in _PICTOGRAMS:
        raise ValueError(
            f"No pictogram of the organ '{key}', only {sorted(_PICTOGRAMS)}."
        )
    if not 0.0 <= loss <= 1.0:
        raise ValueError(f"The lost fraction {loss} is outside [0, 1].")
    index, side = _PICTOGRAMS[key]
    poly = anatomy("human", sex).organs[key].paths[index]
    (x0, y0), (x1, y1) = poly.min(axis=0), poly.max(axis=0)
    s = min(width / (x1 - x0), height / (y1 - y0))
    poly = (poly - [(x0 + x1) / 2, (y0 + y1) / 2]) * s + [cx, cy]
    kept, lost = _split(poly, side, loss)

    def path(p: np.ndarray) -> str:
        return "M" + " ".join(f"{x:.2f},{y:.2f}" for x, y in p) + "Z"

    out = ["<g>"]
    if title:
        out.append(f"<title>{html.escape(title)}</title>")
    if _area(kept) > 1e-6 * _area(poly):
        out.append(f'<path d="{path(kept)}" fill="{color}"/>')
    if _area(lost) > 1e-6 * _area(poly):
        out.append(f'<path d="{path(lost)}" fill="{color}" fill-opacity="0.15"/>')
    out.append(
        f'<path d="{path(poly)}" fill="none" stroke="currentColor" '
        'stroke-opacity="0.7" stroke-width="1"/>'
    )
    out.append("</g>")
    return "".join(out)


def tablet(cx: float, cy: float, d: float = 24.0) -> str:
    """The oral dose: a tablet on a circle."""
    k = d / 22
    c = palette.ORAL
    return (
        ring(cx, cy, d / 2, c)
        + f'<g transform="translate({cx:.1f},{cy:.1f}) scale({k:.3f}) rotate(-30)">'
        f'<rect x="-7" y="-3.2" width="14" height="6.4" rx="3.2" fill="#fff" stroke="{c}" stroke-width="1.2"/>'
        f'<path d="M0,-3.2 v6.4" stroke="{c}" stroke-width="1"/>'
        f'<rect x="-6.4" y="-2.6" width="6.4" height="5.2" rx="2.6" fill="{c}" fill-opacity="0.45"/></g>'
    )


def syringe(cx: float, cy: float, d: float = 24.0, rot: float = 35.0) -> str:
    """The intravenous dose: a syringe on a circle."""
    k = d / 22
    c = palette.IV
    return (
        ring(cx, cy, d / 2, c)
        + f'<g transform="translate({cx:.1f},{cy:.1f}) scale({k:.3f}) rotate({rot:g}) translate(0,8)">'
        f'<path d="M0,0 v-4" stroke="{c}" stroke-width="0.9"/>'
        f'<rect x="-2.2" y="-13" width="4.4" height="9" rx="0.8" fill="#fff" stroke="{c}" stroke-width="1"/>'
        f'<rect x="-2.2" y="-8" width="4.4" height="4" fill="{palette.BLOOD}" fill-opacity="0.5"/>'
        f'<path d="M-3.4,-13 h6.8 M0,-13 v-3 M-2.2,-16 h4.4" stroke="{c}" stroke-width="1"/></g>'
    )


def drop(cx: float, cy: float, d: float = 24.0) -> str:
    """The blood sample: a drop on a circle."""
    k = d / 22
    c = palette.BLOOD
    return (
        ring(cx, cy, d / 2, c)
        + f'<path transform="translate({cx:.1f},{cy:.1f}) scale({k:.3f})" '
        f'd="M0,-6.5 C3,-2.5 4.4,0.5 4.4,2.2 A4.4,4.4 0 0 1 -4.4,2.2 C-4.4,0.5 -3,-2.5 0,-6.5z" fill="{c}"/>'
    )
