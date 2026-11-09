"""The icons of the figures, one graphical language for all of them.

Every icon sits on a white circle: the organs (their shapes in the anatomy
resources; the bowel for the intestine and its lumen, the body for the skin,
fat cells for the adipose tissue, whose drawn depot has no recognisable form at
this size), the oral dose (a tablet), the intravenous dose (a syringe) and the
blood sample (a drop). An icon is an SVG fragment centred at a point. The
pictogram of an impaired organ (`impaired_organ`) is the liver, a kidney or
the heart with its details on a card framed in the color of an impairment
class, its lost part faded.
"""

import hashlib
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


def _hull(points: np.ndarray) -> np.ndarray:
    """The convex hull of points, counter-clockwise (monotone chain)."""
    pts = sorted(map(tuple, points))

    def half(seq: list[tuple[float, float]]) -> list[tuple[float, float]]:
        out: list[tuple[float, float]] = []
        for p in seq:
            while len(out) >= 2:
                (ax, ay), (bx, by) = out[-2], out[-1]
                if (bx - ax) * (p[1] - ay) - (by - ay) * (p[0] - ax) > 0:
                    break
                out.pop()
            out.append(p)
        return out

    lower, upper = half(pts), half(pts[::-1])
    return np.array(lower[:-1] + upper[:-1])


def _hilum(poly: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The deepest point of the concavity of a shape and the direction out of it."""
    hull = _hull(poly)
    best, depth, out = poly[0], -1.0, np.array([1.0, 0.0])
    for p in poly:
        dists = []
        for a, b in zip(hull, np.roll(hull, -1, axis=0), strict=True):
            ab = b - a
            t = float(np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0.0, 1.0))
            q = a + t * ab
            dists.append((float(np.linalg.norm(p - q)), q))
        d, q = min(dists, key=lambda e: e[0])
        if d > depth:
            best, depth = p, d
            out = (q - p) / (np.linalg.norm(q - p) or 1.0)
    return best, out


def _vertical(poly: np.ndarray, x: float) -> tuple[float, float]:
    """The top and the bottom of a shape at `x` (the crossings of its outline)."""
    ys = []
    for p, q in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        if (p[0] - x) * (q[0] - x) <= 0 and p[0] != q[0]:
            ys.append(p[1] + (x - p[0]) / (q[0] - p[0]) * (q[1] - p[1]))
    return (
        (min(ys), max(ys)) if ys else (float(poly[:, 1].min()), float(poly[:, 1].max()))
    )


def _path(p: np.ndarray) -> str:
    return "M" + " ".join(f"{x:.2f},{y:.2f}" for x, y in p) + "Z"


def _details(key: str, poly: np.ndarray) -> list[str]:
    """The drawn details of an organ.

    The lobes and the gallbladder of the liver, the cortex and the vessels and
    the ureter at the hilum of a kidney (the cortex first), the
    interventricular groove of the heart.
    """
    (x0, y0), (x1, y1) = poly.min(axis=0), poly.max(axis=0)
    w, h = x1 - x0, y1 - y0
    ink = palette.ICON_INK
    out = []
    if key == "li":
        x = x0 + 0.57 * w
        top, bottom = _vertical(poly, x)
        out.append(
            f'<path d="M{x:.2f},{top + 0.02 * h:.2f} Q{x - 0.06 * w:.2f},{(top + bottom) / 2:.2f} '
            f'{x - 0.02 * w:.2f},{bottom - 0.02 * h:.2f}" fill="none" stroke="{ink}" '
            f'stroke-opacity="0.55" stroke-width="{0.012 * w:.2f}"/>'
        )
        gx = x0 + 0.40 * w
        _, gb = _vertical(poly, gx)
        out.append(
            f'<ellipse cx="{gx:.2f}" cy="{gb - 0.04 * h:.2f}" rx="{0.06 * w:.2f}" '
            f'ry="{0.1 * h:.2f}" fill="{palette.BILE}" stroke="{ink}" '
            f'stroke-width="{0.01 * w:.2f}" transform="rotate(-25 {gx:.2f} {gb - 0.04 * h:.2f})"/>'
        )
    elif key == "ki":
        c = poly.mean(axis=0)
        inner = (poly - c) * 0.68 + c
        out.append(
            f'<path d="{_path(inner)}" fill="#fff" fill-opacity="0.22" stroke="{ink}" '
            f'stroke-opacity="0.45" stroke-width="{0.025 * w:.2f}"/>'
        )
        hilum, d = _hilum(poly)
        n = np.array([-d[1], d[0]])
        length = 0.7 * w
        for color, offset, bend, width in (
            (palette.BLOOD, -0.11 * h, -0.05 * h, 0.1),
            ("#3d6fb6", 0.0, 0.0, 0.11),
            (palette.URINE, 0.11 * h, 0.35 * h, 0.09),
        ):
            a = hilum + n * offset * 0.4
            b = hilum + d * length + n * offset
            m = (a + b) / 2 + np.array([0.0, bend])
            out.append(
                f'<path d="M{a[0]:.2f},{a[1]:.2f} Q{m[0]:.2f},{m[1]:.2f} {b[0]:.2f},{b[1] + bend:.2f}" '
                f'fill="none" stroke="{color}" stroke-width="{width * w:.2f}" stroke-linecap="round"/>'
            )
    elif key == "he":
        apex = poly[poly[:, 1].argmax()]
        x = x0 + 0.52 * w
        top, _ = _vertical(poly, x)
        start = np.array([x, top + 0.28 * h])
        m = (start + apex) / 2 + np.array([0.06 * w, 0.0])
        out.append(
            f'<path d="M{start[0]:.2f},{start[1]:.2f} Q{m[0]:.2f},{m[1]:.2f} '
            f'{apex[0] - 0.02 * w:.2f},{apex[1] - 0.05 * h:.2f}" fill="none" stroke="{ink}" '
            f'stroke-opacity="0.5" stroke-width="{0.014 * w:.2f}"/>'
        )
    return out


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
    """The pictogram of an impaired organ on a card framed in the color of its class.

    The organ (`li`, `ki` or `he`, its shape in the anatomy resources) is drawn in
    its color of `palette.ORGAN_COLORS`, shaded, with its details (the lobes and the
    gallbladder of the liver, the cortex and the vessels and the ureter at the hilum
    of a kidney, the groove of the heart); its lost part, `loss` of the area, is
    faded. The card is white with a frame in `color`, so the pictogram reads in a
    light and a dark theme alike. The ids of the gradient are a hash of the
    arguments, unique on a page.

    Args:
        key: the organ.
        cx: the centre of the card, x.
        cy: the centre of the card, y.
        width: the width of the card.
        height: the height of the card.
        loss: the lost fraction of the organ, 0 to 1.
        color: the color of the class, the frame of the card.
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
    # the organ in the card, room for the vessels of the kidney
    room = (0.62 if key == "ki" else 0.84) * width, 0.8 * height
    s = min(room[0] / (x1 - x0), room[1] / (y1 - y0))
    shift = -0.12 * width if key == "ki" else 0.0
    poly = (poly - [(x0 + x1) / 2, (y0 + y1) / 2]) * s + [cx + shift, cy]
    _, lost = _split(poly, side, loss)
    uid = hashlib.sha1(
        f"{key}|{cx}|{cy}|{width}|{height}|{loss}|{color}|{sex}".encode()
    ).hexdigest()[:10]
    frame = "#9e9e9e" if color == "black" else color
    fill = palette.ORGAN_COLORS[key]
    shape = _path(poly)
    out = ["<g>"]
    if title:
        out.append(f"<title>{html.escape(title)}</title>")
    out += [
        f'<defs><linearGradient id="io-{uid}" x1="0" y1="0" x2="1" y2="1">'
        '<stop offset="0" stop-color="#fff" stop-opacity="0.35"/>'
        '<stop offset="0.55" stop-color="#fff" stop-opacity="0"/>'
        '<stop offset="1" stop-color="#000" stop-opacity="0.25"/>'
        "</linearGradient></defs>",
        f'<rect x="{cx - width / 2:.2f}" y="{cy - height / 2:.2f}" width="{width:.2f}" '
        f'height="{height:.2f}" rx="{0.08 * min(width, height):.2f}" fill="#fff" '
        f'stroke="{frame}" stroke-width="2"/>',
    ]
    vessels = _details(key, poly) if key == "ki" else []
    # the vessels and the ureter behind the kidney, the rest on the organ
    out += vessels[1:]
    out += [
        f'<path d="{shape}" fill="{fill}" data-part="organ"/>',
        f'<path d="{shape}" fill="url(#io-{uid})"/>',
    ]
    out += vessels[:1] if key == "ki" else _details(key, poly)
    if _area(lost) > 1e-6 * _area(poly):
        out.append(
            f'<path d="{_path(lost)}" fill="#fff" fill-opacity="0.72" data-part="lost"/>'
        )
    out.append(
        f'<path d="{shape}" fill="none" stroke="{palette.ICON_INK}" '
        f'stroke-width="{max(0.8, 0.012 * room[0]):.2f}" stroke-linejoin="round" data-part="outline"/>'
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
