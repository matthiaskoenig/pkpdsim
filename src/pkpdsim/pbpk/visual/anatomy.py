"""The anatomical figure of a model: the male and the female build on their reference individuals.

The man and the woman stand next to each other (their midlines `MIDLINES`
apart, the arms overlapping, the man drawn over the woman), each captioned
with its build. The bodies carry what the definition has:

- the organs in their drawn forms, sized by their volumes (`body.Body`);
- the connections of the model: the portal vein into the liver, the bile from
  the liver over the gallbladder into the duodenum, the ureters to the bladder
  and out, the feces, the oral route from the mouth to the stomach;
- circled icons: the tablet of an oral dose at the mouths, the drop of the
  blood sample at the arm vein of the man, the syringe of an intravenous dose
  at the arm of the woman.

Every organ has a label in one of two columns, the left one (right-aligned)
pointing into the man, the right one into the woman; the builds share the
structure of the definition, so one label serves both. A label shows the
circled icon of the organ, its substances (filled: enters the tissue; dashed:
plasma only) and, for an organ with processes, its submodel: a row per
substance, the compartments plasma, interstitium, cells and the apical side
(bile, urine, lumen) as columns, a chip per species the build has, arrows for
the exchange, the metabolism (with its enzyme), the apical transport (with its
transporter) and the filtration; below it the process lines of the schematic
unless the figure is compact. The credit (the anatomy, the license, pkpdsim)
belongs in the legend of a figure, never inside it.
"""

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np

from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.visual import icons, palette
from pkpdsim.pbpk.visual.body import Body, reference_body
from pkpdsim.pbpk.visual.content import ModelContent, Process

__all__ = [
    "LEFT",
    "MIDLINES",
    "RIGHT",
    "TITLES",
    "AnatomyFigure",
    "anatomy_figure",
]

FONT = "'DejaVu Sans', 'Noto Sans', sans-serif"
MONO = "'DejaVu Sans Mono', monospace"
_CW = 0.6  # the advance of a character of the monospace font, em

MIDLINES = 70.0
"""The distance of the midlines of the man and the woman [cm]."""

LEFT = ("br", "lu", "he", "li", "sk", "ki", "la", "blood")
"""The labels of the left column, pointing into the man."""

RIGHT = ("st", "sp", "pa", "si", "iv", "lumen", "ad", "bo", "mu", "go", "urine")
"""The labels of the right column, pointing into the woman."""

TITLES = {
    "lu": "lung",
    "he": "heart",
    "br": "brain",
    "li": "liver",
    "ki": "kidney",
    "st": "stomach",
    "si": "small intestine",
    "la": "large intestine",
    "sp": "spleen",
    "pa": "pancreas",
    "go": "gonads",
    "mu": "muscle",
    "sk": "skin",
    "ad": "adipose tissue",
    "bo": "bone",
    "lumen": "intestinal lumen",
    "urine": "urine",
    "blood": "blood sample",
    "iv": "intravenous dose",
}
"""The headings of the labels."""

_BADGES = {
    "MET": "#7b5ea7",
    "TRA": "#2f7f9c",
    "CUS": "#777777",
    "URI": "#b08a00",
    "SLW": "#777777",
    "ORL": palette.ORAL,
    "ABS": "#4d8a3a",
    "BIL": palette.BILE,
    "CNV": "#777777",
    "FEC": palette.FECES,
    "OBS": palette.BLOOD,
    "IV": palette.IV,
}
_APICAL = {"li": "bile", "ki": "urine", "si": "lumen", "la": "lumen", "st": "lumen"}
_COLUMNS = (
    ("pl", "plasma", "#f6d5d2"),
    ("is", "interstitium", "#f3ece2"),
    ("ic", "cells", "#ead8c6"),
)
_APICAL_BG = {"bile": "#d8ecd9", "urine": "#f7eeb8", "lumen": "#f6e3c8"}
_CWID, _ROWH, _HEAD, _LANE = 44.0, 17.0, 15.0, 78.0
_ICON = 22.0


def _short(text: str) -> str:
    """A process line without its id and rate law, arrows as arrows."""
    t = re.sub(r"^[A-Za-z0-9_]+: ", "", text)
    t = re.sub(r" \((first order|Michaelis-Menten|Hill)\)$", "", t)
    return t.replace(" -> ", " → ").replace(" by ", " · ")


def _wrap(text: str, n: int) -> list[str]:
    lines: list[str] = []
    cur = ""
    for w in text.split(" "):
        if cur and len(cur) + 1 + len(w) > n:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}" if cur else w
    return [*lines, cur] if cur else lines


def _curve(
    a: np.ndarray,
    b: np.ndarray,
    color: str,
    width: float,
    bend: float = 0.0,
    dash: str = "",
    arrow: bool = False,
) -> str:
    """A quadratic curve from a to b, bent to the side by `bend` [cm]."""
    m = (a + b) / 2
    d = b - a
    n = np.array([-d[1], d[0]]) / (float(np.hypot(*d)) + 1e-9)
    q = m + n * bend
    da = f' stroke-dasharray="{dash}"' if dash else ""
    mk = f' marker-end="url(#arr{color[1:]})"' if arrow else ""
    return (
        f'<path d="M{a[0]:.2f},{a[1]:.2f} Q{q[0]:.2f},{q[1]:.2f} {b[0]:.2f},{b[1]:.2f}" fill="none" '
        f'stroke="{color}" stroke-width="{width}" stroke-linecap="round" stroke-opacity="0.95"{da}{mk}/>'
    )


@dataclass(frozen=True)
class _Glyph:
    """An icon of the overlays: a function of the icons at a point of the body [cm]."""

    kind: str
    at: np.ndarray


class _View:
    """A body with the content of the model: anchors and overlays [cm]."""

    def __init__(
        self, body: Body, content: ModelContent, glyphs: Sequence[str]
    ) -> None:
        self.b = body
        self.c = content
        self.glyphs = tuple(glyphs)
        lm = body.landmarks
        cx = body.anatomy.midline
        self.mouth = body.at(cx, lm["top"] + 0.86 * (lm["neck"] - lm["top"]))
        y_fore = lm["armpit"] + 0.72 * (lm["crotch"] - lm["armpit"])
        self.arm = body.edge_point(y_fore, -1, inset=1.2)
        self.iv = body.edge_point(y_fore - 4, 1, inset=1.2)
        bladder = body.pieces("bladder")
        self.bladder = bladder[0] if bladder else body.centroid("bladder")
        self.urine_end = self.bladder + np.array([1.2, 9.0])

    def anchor(self, key: str, side: int) -> np.ndarray:
        """The point of the body a label of a key points to."""
        b, lm, cx = self.b, self.b.landmarks, self.b.anatomy.midline
        leg = lm["bottom"] - lm["crotch"]
        if key in ("lu", "he", "li", "ki", "la", "st", "sp", "pa", "si"):
            return b.toward(key, side)
        if key == "br":
            return b.centroid("br") + np.array([side * 2.0, 0.0])
        if key == "mu":
            return b.at(cx + side * 8.5, lm["crotch"] + 0.2 * leg)
        if key == "bo":
            femur = np.vstack(b.anatomy.organs["bo"].paths).mean(0)
            x = cx + side * abs(femur[0] - cx)
            return b.at(x, lm["crotch"] + 0.12 * leg)
        if key == "sk":
            return b.edge_point(
                lm["armpit"] + 0.2 * (lm["crotch"] - lm["armpit"]), side
            )
        if key == "ad":
            return b.edge_point(
                lm["waist"] + 0.55 * (lm["crotch"] - lm["waist"]), side, inset=1.0
            )
        if key == "go":
            if b.sex == "male" and "go" in b.anatomy.organs:
                return b.field(
                    np.vstack(b.anatomy.organs["go"].paths).mean(0)[None, :]
                )[0]
            return b.toward("uterus" if "uterus" in b.organs else "bladder", side)
        if key == "lumen":
            return b.toward("si", side) + np.array([0.0, 4.0])
        if key == "urine":
            return self.urine_end
        if key == "blood":
            return self.arm
        if key == "iv":
            return self.iv
        return b.centroid(key)

    def overlays(self, icon_scale: float) -> str:
        """The connections of the model and the circled icons [cm]."""
        b, c = self.b, self.c
        o = []
        liver = b.toward("li", 1)
        liver = liver + (b.centroid("li") - liver) * 0.5
        hub = liver + np.array([3.0, 3.5])
        o.append(_curve(hub, liver, palette.PORTAL, 1.0))
        for k in ("st", "sp", "pa", "si", "la"):
            if k in b.organs:
                pieces = b.pieces(k) if k in ("si", "la") else []
                p = pieces[0] if pieces else b.centroid(k)
                o.append(_curve(p, hub, palette.PORTAL, 0.6, bend=1.0))
                o.append(
                    f'<circle cx="{p[0]:.2f}" cy="{p[1]:.2f}" r="0.55" fill="{palette.PORTAL}"/>'
                )
        if c.bile and "gallbladder" in b.organs:
            gb = b.centroid("gallbladder")
            du = b.centroid("pa") + np.array([-3.5, 2.0])
            o.append(_curve(liver, gb, palette.BILE, 0.45, bend=-0.8))
            o.append(_curve(gb, du, palette.BILE, 0.45, bend=0.8, arrow=True))
        if c.urine:
            bl = self.bladder
            for p in b.pieces("ki"):
                s = float(np.sign(p[0] - bl[0]))
                o.append(
                    _curve(
                        p,
                        bl + np.array([s * 1.5, -1.0]),
                        palette.URINE,
                        0.4,
                        bend=s * 1.2,
                    )
                )
            o.append(
                _curve(
                    bl + np.array([0.6, 2.0]),
                    self.urine_end,
                    palette.URINE,
                    0.5,
                    bend=-0.6,
                    dash="1.0,0.7",
                    arrow=True,
                )
            )
        if c.luminal and "la" in b.organs:
            pieces = b.pieces("la") or [b.centroid("la")]
            low = max(pieces, key=lambda p: p[1])
            end = np.array(
                [self.bladder[0] - 1.2, max(low[1], self.urine_end[1] - 9.0) + 9.0]
            )
            o.append(
                _curve(
                    low, end, palette.FECES, 0.5, bend=1.0, dash="1.0,0.7", arrow=True
                )
            )
        if c.has_route("oral"):
            m = self.mouth
            o.append(
                _curve(
                    m + np.array([0, 1.0]),
                    m + np.array([0, 10.0]),
                    palette.ORAL,
                    0.45,
                    dash="0.9,0.7",
                )
            )
            o.append(
                _curve(
                    m + np.array([0, 10.0]),
                    b.centroid("st"),
                    palette.ORAL,
                    0.45,
                    bend=-2.5,
                    dash="0.9,0.7",
                    arrow=True,
                )
            )
        k = icon_scale
        if "oral" in self.glyphs and c.has_route("oral"):
            p = self.mouth + np.array([4.4, 0.4])
            o.append(
                f'<g data-k="route:oral" transform="translate({p[0]:.2f},{p[1]:.2f}) scale({k:.4f})">{icons.tablet(0, 0, 22)}</g>'
            )
        if "sample" in self.glyphs:
            p = self.arm
            o.append(
                f'<g data-k="sample" transform="translate({p[0]:.2f},{p[1]:.2f}) scale({k:.4f})">{icons.drop(0, 0, 22)}</g>'
            )
        if "iv" in self.glyphs and c.has_route("iv"):
            p = self.iv
            o.append(
                f'<g data-k="route:iv" transform="translate({p[0]:.2f},{p[1]:.2f}) scale({k:.4f})">{icons.syringe(0, 0, 22)}</g>'
            )
        return "".join(o)


class AnatomyFigure:
    """The anatomical figure of the content of a model on the reference man and woman.

    Args:
        content: the content of the definition.
        compact: labels without the process lines (the panel of a publication
            figure, whose schematic shows them).
        man: the body of the man, the reference man by default.
        woman: the body of the woman, the reference woman by default.
    """

    def __init__(
        self,
        content: ModelContent,
        *,
        compact: bool = False,
        man: Body | None = None,
        woman: Body | None = None,
    ) -> None:
        """Prepare the two views."""
        self.c = content
        self.compact = compact
        self.man = _View(man or reference_body("male"), content, ("oral", "sample"))
        self.woman = _View(woman or reference_body("female"), content, ("oral", "iv"))
        self.size = (0, 0)

    # --- labels ---------------------------------------------------------------------
    def shown(self, key: str) -> bool:
        """Whether a label of a key applies to the model."""
        if key == "lumen":
            return bool(self.c.lumen_lines)
        if key == "urine":
            return bool(self.c.urine)
        if key == "iv":
            return self.c.has_route("iv")
        return True

    def has_submodel(self, key: str) -> bool:
        """Whether an organ has processes (or the urine), so its label shows its submodel."""
        return any(p.organ == key and p.kind != "CUS" for p in self.c.processes) or (
            key == "ki" and bool(self.c.urine)
        )

    def submodel(self, key: str, x: float, y: float) -> tuple[str, float, float]:
        """The submodel grid of an organ at (x, y): the SVG, its width and height."""
        org = self.c.organs[key]
        subs = [s.sid for s in self.c.substances]
        apical = _APICAL.get(key)
        chem = [p for p in self.c.processes if p.organ == key]
        to_apical = {
            p.substrate
            for p in chem
            if p.kind == "TRA" and p.membrane == "apical" and p.direction == "efflux"
        }
        if key == "ki":
            to_apical |= set(self.c.urine)
        xs = {"pl": x, "is": x + _CWID + 4, "ic": x + 2 * (_CWID + 4)}
        mets = [p for p in chem if p.kind == "MET"]
        tag_w = max(
            [len(p.protein) * 8.5 * 0.7 + 6 for p in mets if p.protein] or [0.0]
        )
        lane_w = max(_LANE, 8 + 12 * len(mets) + (tag_w + 12 if tag_w else 0))
        xl = xs["ic"] + _CWID + 4
        xa = xl + lane_w
        h = _HEAD + len(subs) * _ROWH + 6
        w = (xa + _CWID - x) if apical else (xl - x)
        o = []
        for c, label, bg in _COLUMNS:
            o.append(
                f'<rect x="{xs[c]:.1f}" y="{y:.1f}" width="{_CWID}" height="{h}" rx="4" fill="{bg}"/>'
                f'<text x="{xs[c] + _CWID / 2:.1f}" y="{y + 10:.1f}" class="sh">{label}</text>'
            )
        if apical:
            o.append(
                f'<rect x="{xa:.1f}" y="{y:.1f}" width="{_CWID}" height="{h}" rx="4" fill="{_APICAL_BG[apical]}"/>'
                f'<text x="{xa + _CWID / 2:.1f}" y="{y + 10:.1f}" class="sh">{apical}</text>'
            )
        ry = {s: y + _HEAD + _ROWH * i + _ROWH / 2 for i, s in enumerate(subs)}
        colors = {s.sid: s.color for s in self.c.substances}

        def chip(cx: float, cy: float, s: str, width: float = 0.0) -> str:
            if width:
                return f'<rect x="{cx - width / 2:.1f}" y="{cy - 5:.1f}" width="{width:.1f}" height="10" rx="5" fill="{colors[s]}"/>'
            return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="5" fill="{colors[s]}"/>'

        def cen(c: str) -> float:
            return xs[c] + _CWID / 2

        def exchange(x0: float, x1: float, yy: float) -> str:
            return f'<path d="M{x0:.1f},{yy:.1f} H{x1:.1f}" stroke="#666" stroke-width="0.9" marker-end="url(#sa)" marker-start="url(#sa)"/>'

        for s in subs:
            form, yy = org.parts[s], ry[s]
            o.append(chip(cen("pl"), yy, s))
            if form == "merged":
                o.append(chip((cen("is") + cen("ic")) / 2, yy, s, width=_CWID + 14))
                o.append(exchange(cen("pl") + 6, cen("is") - 9, yy))
            elif form == "full":
                o.append(chip(cen("is"), yy, s))
                o.append(chip(cen("ic"), yy, s))
                o.append(exchange(cen("pl") + 6, cen("is") - 6, yy))
                o.append(exchange(cen("is") + 6, cen("ic") - 6, yy))
            elif form == "organ":
                o.append(chip((cen("pl") + cen("ic")) / 2, yy, s, width=2 * _CWID + 30))
            if apical and s in to_apical:
                o.append(chip(xa + _CWID / 2, yy, s))
        for p in chem:
            if (
                p.kind == "TRA"
                and p.membrane == "apical"
                and apical
                and p.substrate in ry
            ):
                yy = ry[p.substrate]
                o.append(
                    f'<path d="M{cen("ic") + 6:.1f},{yy - 2.5:.1f} H{xa + _CWID / 2 - 7:.1f}" stroke="#2f7f9c" '
                    f'stroke-width="1.1" marker-end="url(#st)"/>'
                )
                if p.protein:
                    o.append(
                        f'<text x="{xa - 22:.1f}" y="{yy - 5:.1f}" class="st2">{escape(p.protein)}</text>'
                    )

        # the enzyme of a metabolism as a tag in a gap between two rows within its span,
        # the narrowest spans first, so no tag lies on another or on a transport
        def span(p: Process) -> tuple[int, int]:
            rows = [subs.index(p.substrate)] + [
                subs.index(q) for q in p.products if q in subs
            ]
            return min(rows), max(rows)

        used: set[float] = set()
        slot_of: dict[str, float] = {}
        for p in sorted(
            (q for q in mets if q.protein), key=lambda q: span(q)[1] - span(q)[0]
        ):
            lo, hi = span(p)
            ideal = (lo + hi) / 2
            inside = [i + 0.5 for i in range(lo, hi)]
            every = [i + 0.5 for i in range(-1, len(subs))]
            free = [s for s in inside if s not in used] or [
                s for s in every if s not in used
            ]
            slot_of[p.sid] = min(free, key=lambda s: (abs(s - ideal), s))
            used.add(slot_of[p.sid])
        tx = xl + 8 + 12 * len(mets) + 4
        for lane, p in enumerate(mets):
            if p.substrate not in ry:
                continue
            y0, lx = ry[p.substrate], xl + 8 + 12 * lane
            for prod in p.products:
                if prod in ry:
                    o.append(
                        f'<path d="M{cen("ic") + 6:.1f},{y0:.1f} H{lx:.1f} V{ry[prod]:.1f} H{cen("ic") + 7:.1f}" '
                        f'fill="none" stroke="#7b5ea7" stroke-width="1.1" marker-end="url(#sm)"/>'
                    )
            if p.protein:
                ym = y + _HEAD + _ROWH * slot_of[p.sid] + _ROWH / 2
                tw = len(p.protein) * 8.5 * 0.7 + 6
                o.append(
                    f'<path d="M{lx:.1f},{ym:.1f} H{tx:.1f}" stroke="#7b5ea7" stroke-width="0.7" stroke-dasharray="1.5,1.5"/>'
                    f'<rect x="{tx:.1f}" y="{ym - 5.5:.1f}" width="{tw:.1f}" height="11" rx="2" fill="#fff" '
                    f'stroke="#7b5ea7" stroke-width="0.7"/><text x="{tx + 3:.1f}" y="{ym + 3:.1f}" class="se">{escape(p.protein)}</text>'
                )
        if key == "ki" and self.c.urine and apical:
            yb = y + h + 6
            o.append(
                f'<path d="M{cen("pl"):.1f},{y + h - 2:.1f} V{yb:.1f} H{xa + _CWID / 2:.1f} V{y + h - 2:.1f}" fill="none" '
                f'stroke="{palette.URINE}" stroke-width="1.2" marker-end="url(#su)"/>'
                f'<text x="{(cen("pl") + xa) / 2 + 20:.1f}" y="{yb + 10:.1f}" class="sh">filtration</text>'
            )
            h += 18
        return "".join(o), w, h

    def label(
        self, key: str, x: float, y: float, width: float, side: int
    ) -> tuple[str, float]:
        """A label at (x, y) of a width; the left column (`side` -1) right-aligned; the SVG and its height."""
        o = []
        right = side < 0
        xe = x + width
        org = self.c.organs.get(key)
        title = TITLES.get(key, key)
        variant = f" ({org.variant})" if org is not None and org.variant else ""
        tail = f"{key if org is not None else ''}{variant}"
        if key == "iv":
            icon = icons.syringe(0, 0, _ICON)
        elif key == "blood":
            icon = icons.drop(0, 0, _ICON)
        else:
            icon = icons.organ(key, 0, 0, _ICON)
        icon_x = (xe - _ICON / 2) if right else (x + _ICON / 2)
        if icon:
            o.append(f'<g transform="translate({icon_x:.1f},{y + 8:.1f})">{icon}</g>')
        off = _ICON + 6 if icon else 0
        if right:
            tw = len(tail) * 11 * 0.64
            xe2 = xe - off
            o.append(
                f'<text x="{xe2:.1f}" y="{y + 12:.1f}" class="lid" text-anchor="end">{escape(tail)}</text>'
                f'<text x="{xe2 - tw - (7 if tail else 0):.1f}" y="{y + 12:.1f}" class="lt" text-anchor="end">{escape(title)}</text>'
            )
        else:
            o.append(
                f'<text x="{x + off:.1f}" y="{y + 12:.1f}" class="lt">{escape(title)}'
                f'<tspan class="lid"> {escape(tail)}</tspan></text>'
            )
        h = 16.0
        chips: list[tuple[str, str, bool]] = (
            [(ch.sid, ch.color, ch.tissue) for ch in org.chips]
            if org is not None
            else []
        )
        if key == "lumen":
            chips = [(s, self.c.color(s), True) for s in self.c.luminal]
        elif key == "urine":
            chips = [(s, self.c.color(s), True) for s in self.c.urine]
        elif key == "blood":
            chips = [(s.sid, s.color, True) for s in self.c.substances]
        r, step = 4.2, 11.4
        cx = (xe - off - len(chips) * step + 3) if right else (x + off)
        for sid, color, tissue in chips:
            if tissue:
                o.append(
                    f'<circle data-k="sub:{escape(sid)}" cx="{cx + r:.1f}" cy="{y + h + 5:.1f}" r="{r}" fill="{color}"/>'
                )
            else:
                o.append(
                    f'<circle data-k="sub:{escape(sid)}" cx="{cx + r:.1f}" cy="{y + h + 5:.1f}" r="{r - 0.7:.1f}" fill="#fff" '
                    f'stroke="{color}" stroke-width="1.4" stroke-dasharray="2,1.4"/>'
                )
            cx += step
        if chips:
            h += 14
        if org is not None and self.has_submodel(key):
            g, sw, sh = self.submodel(key, 0, 0)
            ox = (xe - sw) if right else x
            o.append(f'<g transform="translate({ox:.1f},{y + h + 2:.1f})">{g}</g>')
            h += sh + 8
        if self.compact:
            return "".join(o), h + 2
        lines: list[tuple[str, str]] = list(org.lines) if org is not None else []
        if key == "lumen":
            lines = list(self.c.lumen_lines)
        elif key == "urine":
            lines = []
        elif key == "blood":
            lines = [("OBS", "arm vein: the observer Carm of every substance")]
        elif key == "iv":
            lines = [
                ("IV", f"{s} into the venous blood")
                for k, s in self.c.routes
                if k == "iv"
            ]
        n = int(width / (10 * _CW)) - 5
        for badge, text in lines:
            t = _short(text) if badge in ("MET", "TRA", "CUS") else text
            parts = _wrap(t, n)
            tw = max(len(p) for p in parts) * 10 * _CW
            bx = (xe - tw - 31) if right else x
            o.append(
                f'<g data-k="txt:{escape(badge)}: {escape(text, {chr(34): "&quot;"})}">'
                f'<rect x="{bx:.1f}" y="{y + h + 1:.1f}" width="26" height="12" rx="2" fill="{_BADGES.get(badge, "#777")}"/>'
                f'<text x="{bx + 13:.1f}" y="{y + h + 10.2:.1f}" class="bd">{escape(badge)}</text>'
            )
            for j, p in enumerate(parts):
                o.append(
                    f'<text x="{bx + 31:.1f}" y="{y + h + 10.5 + 12 * j:.1f}" class="pl">{escape(p)}</text>'
                )
            o.append("</g>")
            h += 12 * len(parts) + 3
        return "".join(o), h + 4

    @staticmethod
    def place(
        items: Sequence[tuple[str, float, float]],
        top: float,
        bottom: float,
        gap: float = 8.0,
    ) -> list[tuple[str, float, float]]:
        """The tops of labels near the heights of their anchors.

        Overlapping labels merge into a cluster centred on the anchors of its
        members (pool adjacent violators); the clusters stay inside [top, bottom].
        """
        ordered = sorted(items, key=lambda t: t[1])

        def start(members: list[tuple[str, float, float]]) -> tuple[float, float]:
            off, starts = 0.0, []
            for _, want, h in members:
                starts.append(want - 8 - off)
                off += h + gap
            return float(np.mean(starts)), off - gap

        clusters: list[list] = []
        for item in ordered:
            clusters.append([[item], item[1] - 8])
            while True:
                st, _ = start(clusters[-1][0])
                clusters[-1][1] = max(st, top)
                if len(clusters) > 1:
                    _, length = start(clusters[-2][0])
                    if clusters[-2][1] + length + gap > clusters[-1][1]:
                        last = clusters.pop()
                        clusters[-1][0] += last[0]
                        continue
                break
        end = clusters[-1][1] + start(clusters[-1][0])[1] if clusters else top
        if end > bottom:
            for c in clusters:
                c[1] -= end - bottom
        out = []
        prev = top - gap
        for members, st in clusters:
            yy = max(st, prev + gap)
            for k, _, h in members:
                out.append((k, yy, h))
                yy += h + gap
            prev = yy - gap
        return out

    # --- the figure ------------------------------------------------------------------
    def svg(
        self,
        *,
        title: bool = True,
        substances: bool = True,
        scale: float = 3.3,
        text: float = 1.3,
        left: Sequence[str] = LEFT,
        right: Sequence[str] = RIGHT,
        tags: bool = False,
    ) -> str:
        """The figure as a standalone SVG.

        Args:
            title: the name of the model and its builds at the top.
            substances: the substances with their colours at the top.
            scale: the bodies [px per cm].
            text: the scale of the labels, icons and the key.
            left: the labels of the left column.
            right: the labels of the right column.
            tags: mark the organs of the bodies for hover cards.

        Returns:
            The SVG, deterministic.
        """
        man, woman = self.man, self.woman
        wl = 330.0
        wcol, gap = wl * text, 34.0
        xm0, xm1, hm = man.b.bounds()
        xf0, xf1, hf = woman.b.bounds()
        lo, hi = min(xm0, MIDLINES + xf0), max(xm1, MIDLINES + xf1)
        top = 112.0 if title else (70.0 if substances else 12.0)
        hmax = max(hm, hf)
        oy = top + (10 if (title or substances) else 0)
        oxm = wcol + gap - lo * scale
        oym = oy + (hmax - hm) * scale
        oxf = oxm + MIDLINES * scale
        oyf = oy + (hmax - hf) * scale
        width = round(2 * wcol + 2 * gap + (hi - lo) * scale)
        foot = oy + hmax * scale
        k = text / scale
        bodies = (
            f'<g transform="translate({oxf:.1f},{oyf:.1f}) scale({scale})">{woman.b.svg(tags=tags)}{woman.overlays(k)}</g>'
            f'<g transform="translate({oxm:.1f},{oym:.1f}) scale({scale})">{man.b.svg(tags=tags)}{man.overlays(k)}</g>'
        )

        def px_m(p: np.ndarray) -> np.ndarray:
            return np.array([oxm + p[0] * scale, oym + p[1] * scale])

        def px_f(p: np.ndarray) -> np.ndarray:
            return np.array([oxf + p[0] * scale, oyf + p[1] * scale])

        col_top = max(4.0, top - (20 if (title or substances) else 6))
        key_h = 50 * text
        need = 0.0
        for side, keys in ((-1, left), (1, right)):
            hs = [
                self.label(k2, 0, 0, wl - 20, side)[1] * text
                for k2 in keys
                if self.shown(k2)
            ]
            need = max(need, sum(hs) + 8 * len(hs))
        height = round(max(foot + 30 * text + key_h + 8, col_top + need + key_h + 8))
        out = []
        for side, keys, view, px in ((-1, left, man, px_m), (1, right, woman, px_f)):
            items = [
                (
                    k2,
                    float(px(view.anchor(k2, side))[1]),
                    self.label(k2, 0, 0, wl - 20, side)[1] * text,
                )
                for k2 in keys
                if self.shown(k2)
            ]
            for k2, yy, _ in self.place(items, col_top, height - key_h - 8):
                lx = 0.0 if side < 0 else width - wcol
                g, _ = self.label(k2, 6 if side < 0 else 14, 0, wl - 20, side)
                out.append(
                    f'<g data-k="{k2}" transform="translate({lx:.1f},{yy:.1f}) scale({text})">{g}</g>'
                )
                a = px(view.anchor(k2, side))
                ex = (wcol + 2) if side < 0 else (width - wcol + 6)
                ly = yy + 8 * text
                out.append(
                    f'<path d="M{a[0]:.1f},{a[1]:.1f} L{ex - side * 14:.1f},{ly:.1f} L{ex:.1f},{ly:.1f}" fill="none" '
                    f'stroke="#444" stroke-width="{0.8 * text:.2f}"/>'
                    f'<circle cx="{a[0]:.1f}" cy="{a[1]:.1f}" r="{2.4 * text:.1f}" fill="#fff" stroke="#222" stroke-width="1"/>'
                )
        captions = []
        for view, cx in ((man, oxm), (woman, oxf)):
            q = view.b.quantities
            yb = foot + 14 * text
            captions.append(
                f'<text x="{cx:.1f}" y="{yb:.1f}" class="bc" text-anchor="middle">{view.b.sex} build</text>'
                f'<text x="{cx:.1f}" y="{yb + 14 * text:.1f}" class="bs" text-anchor="middle">{q.height:.0f} cm, {q.bw:.0f} kg</text>'
            )
        if self.c.has_route("oral"):
            for view, ox, oyb in ((man, oxm, oym), (woman, oxf, oyf)):
                p = (
                    ox + (view.mouth[0] + 4.4) * scale + 13 * text,
                    oyb + (view.mouth[1] + 0.4) * scale + 4 * text,
                )
                captions.append(
                    f'<text x="{p[0]:.1f}" y="{p[1]:.1f}" class="gl" fill="{palette.ORAL}">oral</text>'
                )
        header = []
        if title:
            header += [
                f'<text x="6" y="22" class="tt">{escape(self.c.name)}</text>',
                f'<text x="6" y="40" class="st">{escape(self.c.sid)} · male and female builds · reference individuals</text>',
            ]
        if title or substances:
            x, y = 6.0, 62.0 if title else 22.0
            for s in self.c.substances:
                t = f"{s.sid}  {s.name}"
                w = 14 + len(t) * 6.4 + 18
                if x + w > width:
                    x, y = 6.0, y + 18
                header.append(
                    f'<g data-k="sub:{escape(s.sid)}"><circle cx="{x + 5:.1f}" cy="{y - 4:.1f}" r="5" fill="{s.color}"/>'
                    f'<text x="{x + 14:.1f}" y="{y:.1f}" class="sub"><tspan font-weight="600">{escape(s.sid)}</tspan> {escape(s.name)}</text></g>'
                )
                x += w
        key = f'<g transform="translate(0,{height - key_h + 10 * text:.1f}) scale({text})">{self._key()}</g>'
        self.size = (width, height)
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
            f'font-family="{FONT}">{self._style(text)}{self._defs()}<rect width="{width}" height="{height}" fill="#fff"/>'
            f"{''.join(header)}{bodies}{''.join(out)}{''.join(captions)}{key}</svg>"
        )

    def _key(self) -> str:
        """The key of the lines, chips, layers and icons (two rows)."""
        k = []
        kx = 6.0
        for color, label, dash in (
            (palette.PORTAL, "portal vein", ""),
            (palette.BILE, "bile", ""),
            (palette.URINE, "urine", ""),
            (palette.FECES, "feces", "4,3"),
            (palette.ORAL, "oral route", "4,3"),
            (palette.BLOOD, "arteries", ""),
        ):
            da = f' stroke-dasharray="{dash}"' if dash else ""
            k.append(
                f'<path d="M{kx:.1f},-4 h22" stroke="{color}" stroke-width="3"{da}/><text x="{kx + 27:.1f}" y="0" class="ky">{label}</text>'
            )
            kx += 34 + len(label) * 6.6 + 14
        for fill, label in ((True, "enters the tissue"), (False, "plasma only")):
            style = (
                'fill="#555"'
                if fill
                else 'fill="#fff" stroke="#555" stroke-width="1.4" stroke-dasharray="2,1.4"'
            )
            k.append(
                f'<circle cx="{kx + 5:.1f}" cy="-4" r="4.2" {style}/><text x="{kx + 14:.1f}" y="0" class="ky">{label}</text>'
            )
            kx += 22 + len(label) * 6.6 + 14
        kx, ky = 6.0, 22.0
        for color, label in (
            (palette.FAT, "subcutaneous adipose"),
            (palette.VAT, "visceral adipose"),
            (palette.MUSCLE, "muscle"),
        ):
            k.append(
                f'<rect x="{kx:.1f}" y="{ky - 10:.1f}" width="14" height="10" fill="{color}" stroke="#999" stroke-width="0.5"/>'
                f'<text x="{kx + 19:.1f}" y="{ky:.1f}" class="ky">{label}</text>'
            )
            kx += 26 + len(label) * 6.6 + 14
        for kind, label in (
            ("tablet", "oral dose"),
            ("syringe", "intravenous dose"),
            ("drop", "blood sample"),
        ):
            icon = getattr(icons, kind)(kx + 9, ky - 4, 18)
            k.append(
                f'{icon}<text x="{kx + 22:.1f}" y="{ky:.1f}" class="ky">{label}</text>'
            )
            kx += 34 + len(label) * 6.6 + 14
        return "".join(k)

    @staticmethod
    def _style(text: float) -> str:
        f = text
        return (
            f"<style>.tt{{font:700 18px {FONT};fill:#111}}.st{{font:12px {FONT};fill:#555}}.sub{{font:12px {FONT};fill:#222}}"
            f".lt{{font:600 12.5px {FONT};fill:#111}}.lid{{font:400 11px {MONO};fill:#777}}"
            f".bd{{font:700 8.5px {MONO};fill:#fff;text-anchor:middle}}.pl{{font:10px {MONO};fill:#222}}"
            f".gl{{font:600 {11 * f:.1f}px {FONT}}}.ky{{font:11px {FONT};fill:#333}}"
            f".bc{{font:600 {12.5 * f:.1f}px {FONT};fill:#111}}.bs{{font:{11 * f:.1f}px {FONT};fill:#555}}"
            f".sh{{font:8.5px {FONT};fill:#555;text-anchor:middle}}.se{{font:600 8.5px {FONT};fill:#7b5ea7}}"
            f".st2{{font:600 8.5px {FONT};fill:#2f7f9c;text-anchor:middle}}</style>"
        )

    @staticmethod
    def _defs() -> str:
        def marker(mid: str, color: str, size: int = 4, reverse: bool = False) -> str:
            orient = "auto-start-reverse" if reverse else "auto"
            return (
                f'<marker id="{mid}" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="{size}" markerHeight="{size}" '
                f'orient="{orient}"><path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>'
            )

        colored = (
            palette.PORTAL,
            palette.BILE,
            palette.URINE,
            palette.FECES,
            palette.ORAL,
        )
        return (
            "<defs>"
            + "".join(marker(f"arr{c[1:]}", c) for c in colored)
            + marker("sa", "#666", 4, reverse=True)
            + marker("sm", "#7b5ea7", 5)
            + marker("st", "#2f7f9c", 5)
            + marker("su", palette.URINE, 5)
            + "</defs>"
        )


def anatomy_figure(
    definition: PBPKDefinition,
    path: Path,
    *,
    compact: bool = False,
    title: bool = True,
    substances: bool = True,
    scale: float = 3.3,
    text: float = 1.3,
) -> Path:
    """Write the anatomical figure of a definition, the male and the female build.

    Args:
        definition: the definition.
        path: the SVG file to write.
        compact: labels without the process lines.
        title: the name of the model at the top.
        substances: the substances at the top.
        scale: the bodies [px per cm].
        text: the scale of the labels, icons and the key.

    Returns:
        The path of the figure.
    """
    figure = AnatomyFigure(ModelContent.of(definition), compact=compact)
    svg = figure.svg(title=title, substances=substances, scale=scale, text=text)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(svg + "\n", encoding="utf-8")
    return path
