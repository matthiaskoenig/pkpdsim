"""The anatomical figure of a machine perfusion: the organ in its circuit (panel A).

`circuit_anatomy_figure(definition, path)` draws the mounted organ and the
apparatus around it, the companion of the schematic of
`pkpdsim.perfusion.figure` (panel B), in the language of the anatomical
figure of the PBPK library (`pkpdsim.pbpk.visual.anatomy`: the organ in its
drawn form and tissue colour, circled grey icons, label columns with leader
lines, the chips of the quantities):

- in the middle the organ, large, its shape from the anatomy resources of
  the PBPK library (`resources.anatomy("human", "male")`, the human shape
  for every species as the PBPK figures of the animals have it; one kidney of
  the two, the largest part of a shape of several parts) in its colour of
  `palette.ORGAN_COLORS`, an organ without a shape as a neutral rounded
  shape;
- the circuit as apparatus: the reservoir as a vessel with its perfusate,
  the pump, the oxygenator with the O2 of the sweep gas in and the CO2 out,
  the inflow tubing into the inflow of the organ (the portal vein of the
  liver, the renal artery of the kidney), the outflow (the hepatic vein, the
  renal vein) through the outflow tubing back into the reservoir, the flow
  as arrows, the bile duct into the bile collector of the liver or the ureter
  into the urine collector of the kidney, and the medium change as dashed
  arrows into and out of the reservoir;
- two label columns, the apparatus on the left (right-aligned), the organ
  and its collector on the right, each label with its circled grey icon and
  a leader line: the reservoir with the chips of the perfusate, the pump, the
  oxygenator with its gases, the temperature, the organ with its chips
  (filled where the cells hold the quantity as well, dashed where it stays in
  the perfusate) and the inset of its zones, a strip in the order of the
  flow (from periportal to pericentral along the sinusoid of the liver; the
  zones of another organ when it has more than one), each zone a hover key
  `zone:<zone>`; the key of the lines and chips at the bottom.

Every text is monospace and start-anchored with its font size (the
`_Canvas.text` of the PBPK schematic), so the label checks of the figures
apply; the hover cards are `figure.circuit_cards` and the card of the medium
change (`anatomy_cards`). Pure python and deterministic.
"""

from collections.abc import Callable, Sequence
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np

from pkpdsim.pbpk.figure import (
    _ICON,
    _PAD,
    _SIZE,
    _SMALL,
    _STYLE,
    _TITLE,
    _Canvas,
    _f,
    _flow,
    _tag_width,
    _width,
    _wrap,
)
from pkpdsim.pbpk.visual import icons, palette
from pkpdsim.pbpk.visual.anatomy import AnatomyFigure
from pkpdsim.pbpk.visual.resources import anatomy
from pkpdsim.perfusion import ids
from pkpdsim.perfusion.circuit import quantities
from pkpdsim.perfusion.definition import (
    KidneyOrgan,
    LiverOrgan,
    PerfusionDefinition,
)
from pkpdsim.perfusion.figure import (
    BADGES,
    _check,
    _colors,
    _ends,
    _gas_pressure,
    _n,
    _organ_name,
    _oxygenator,
    _pump,
    _thermometer,
    _tissue,
    _urine,
    _zones,
    circuit_cards,
)

__all__ = [
    "NEUTRAL",
    "anatomy_cards",
    "circuit_anatomy_figure",
    "circuit_anatomy_svg",
]

NEUTRAL = "#d9d4cf"
"""The colour of an organ without a colour of the PBPK figures."""

_COLUMN = 250.0
"""The width of a label column (px)."""
_GAP = 36.0
"""The gap between a label column and the drawing, the bend of the leaders (px)."""
_WD, _HD = 440.0, 420.0
"""The size of the drawing (px)."""
_TOP = 70.0
"""The top of the drawing and of the label columns (px), below the title."""
_PIPE = "#8c959f"
"""The tubing of the perfusate."""
_PERFUSATE = "#e3edf6"
"""The perfusate in the reservoir."""
_FIT = (150.0, 110.0, 270.0, 170.0)
"""The box of the organ in the drawing: x, y, width, height (px)."""
_OXY, _PUMP_AT = (60.0, 60.0), (60.0, 175.0)
"""The oxygenator and the pump in the drawing (px)."""
_RES = (25.0, 270.0, 70.0, 90.0)
"""The reservoir in the drawing: x, y, width, height (px)."""
_RETURN = 335.0
"""The height of the outflow tubing into the reservoir in the drawing (px)."""
_COLLECTOR = (318.0, 44.0, 36.0)
"""The bile or urine collector: its top, height and width (px)."""
_VEINS = {
    "li": ("portal vein", "hepatic vein"),
    "ki": ("renal artery", "renal vein"),
}
"""The inflow and the outflow of an organ, the artery and the vein of another."""
_LINE = 18.0
"""The height of a line of a label (px)."""


def _polyline(points: np.ndarray) -> str:
    """The path data of a closed polyline."""
    return "M" + " ".join(f"{x:.2f},{y:.2f}" for x, y in points) + "Z"


def _organ_shape(key: str) -> np.ndarray | None:
    """The outline of an organ in the drawing of the human man, `None` without a shape.

    A shape of several parts gives its largest part (a muscle of a limb),
    the organ of one perfusion; the kidney is the one drawn on the right
    (the left kidney of the body), whose hilum faces the circuit.
    """
    shape = anatomy("human", "male").shape(key)
    if shape is None:
        return None
    if key == "ki":
        return max(shape.paths, key=lambda p: float(p[:, 0].mean()))

    def area(p: np.ndarray) -> float:
        w, h = p.max(axis=0) - p.min(axis=0)
        return float(w * h)

    return max(shape.paths, key=area)


class _Anatomy:
    """The layout of the anatomical figure of a definition."""

    def __init__(self, definition: PerfusionDefinition) -> None:
        """Collect what the figure shows."""
        self.d = definition
        self.organ = definition.organ.resolved(definition.species)
        self.quantities = list(quantities(definition.perfusate, definition.released))
        self.colors = _colors(self.quantities, definition.perfusate.tracers)
        self.canvas = _Canvas(badges=BADGES)
        self.dx = _PAD + _COLUMN + _GAP
        self.rx = self.dx + _WD + _GAP
        self.width = self.rx + _COLUMN + _PAD
        self.anchors: dict[str, tuple[float, float]] = {}
        self.hilum = False
        """Whether the vessels enter the organ from the side (the hilum of
        the kidney), else the inflow from above and the outflow below."""
        self.left = 0.0
        """The left edge of the organ in the drawing."""

    def p(self, x: float, y: float) -> tuple[float, float]:
        """A point of the drawing in the figure."""
        return self.dx + x, _TOP + y

    # --- the drawing -------------------------------------------------------

    def _pipe(self, points: Sequence[tuple[float, float]], color: str = _PIPE) -> None:
        """Draw a tube along the points of the drawing, its flow an arrow at the end."""
        xy = " ".join(f"{_f(x)},{_f(y)}" for x, y in (self.p(*q) for q in points))
        marker = "pka-flow" if color == _PIPE else f"pka-{color[1:]}"
        width = 3 if color == _PIPE else 2
        self.canvas.elements.append(
            f'<polyline points="{xy}" fill="none" stroke="{color}" '
            f'stroke-width="{width}" stroke-linejoin="round" marker-end="url(#{marker})"/>'
        )

    def _text(self, x: float, y: float, text: str, cls: str = "pkf-s pkf-dim") -> None:
        """Draw a small text of the drawing at its baseline (drawing coordinates)."""
        gx, gy = self.p(x, y)
        self.canvas.text(gx, gy, text, cls, _SMALL)

    def _organ(
        self,
    ) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
        """Draw the organ; get its inflow, outflow and excretion points (drawing)."""
        c = self.canvas
        key = self.organ.key
        fx, fy, fw, fh = _FIT
        outline = _organ_shape(key)
        if outline is None:
            # a neutral rounded shape of the proportions of a small organ
            w, h = fw * 0.6, fh * 0.75
            x0, y0 = fx + (fw - w) / 2, fy + (fh - h) / 2
            gx, gy = self.p(x0, y0)
            c.open("pka-organ", key)
            c.rect(
                gx,
                gy,
                w,
                h,
                style=f"fill:{NEUTRAL};stroke:{palette.ICON_INK};stroke-width:1.2",
                rx=h / 3,
            )
            c.close()
            points = np.array(
                [[x0 + t * w, y] for t in np.linspace(0, 1, 41) for y in (y0, y0 + h)]
                + [[x, y0 + t * h] for t in np.linspace(0, 1, 41) for x in (x0, x0 + w)]
            )
        else:
            lo, hi = outline.min(axis=0), outline.max(axis=0)
            size = hi - lo
            s = min(fw / size[0], fh / size[1])
            offset = np.array(
                [fx + (fw - s * size[0]) / 2, fy + (fh - s * size[1]) / 2]
            )
            points = (outline - lo) * s + offset
            gpoints = points + np.array([self.dx, _TOP])
            color = palette.ORGAN_COLORS.get(key, NEUTRAL)
            c.open("pka-organ", key)
            c.elements.append(
                f'<path d="{_polyline(gpoints)}" fill="{color}" '
                f'stroke="{palette.ICON_INK}" stroke-width="1.2" stroke-linejoin="round"/>'
            )
            c.close()
        x0, y0 = points.min(axis=0)
        x1, y1 = points.max(axis=0)
        w = x1 - x0
        self.left = float(x0)

        def nearest(tx: float, ty: float) -> tuple[float, float]:
            i = int(np.argmin(np.hypot(points[:, 0] - tx, points[:, 1] - ty)))
            return float(points[i, 0]), float(points[i, 1])

        right = nearest(x1, (y0 + y1) / 2)
        self.anchors[key] = self.p(*right)
        if isinstance(self.organ, KidneyOrgan):
            # the renal artery and vein at the hilum, the deepest point of the
            # concave side which faces the circuit
            h, ym = y1 - y0, (y0 + y1) / 2
            side = points[
                (points[:, 0] < x0 + w / 2) & (np.abs(points[:, 1] - ym) < 0.2 * h)
            ]
            hx, hy = side[int(np.argmax(side[:, 0]))]

            def at(ty: float) -> tuple[float, float]:
                i = int(np.argmin(np.hypot(side[:, 0] - hx, side[:, 1] - ty)))
                return float(side[i, 0]), float(side[i, 1])

            self.hilum = True
            return at(hy - 0.1 * h), at(hy + 0.1 * h), nearest(x0 + 0.6 * w, y1)
        return (
            nearest(x0 + 0.3 * w, y0),
            nearest(x0 + 0.3 * w, y1),
            nearest(x0 + 0.7 * w, y1),
        )

    def _drawing(self) -> None:
        """Draw the organ and the circuit."""
        c = self.canvas
        organ = self.organ
        inflow, outflow, excretion = self._organ()
        vein_in, vein_out = _VEINS.get(organ.key, ("artery", "vein"))
        circuit = self.d.circuit
        ox, oy = _OXY
        px, py = _PUMP_AT
        rx, ry, rw, rh = _RES
        r = _ICON / 2 + 2
        # the reservoir: a vessel with its perfusate, the dip tube to the pump
        c.open("pka-res", "reservoir")
        gx, gy = self.p(rx, ry)
        c.rect(gx, gy + 18, rw, rh - 18, style=f"fill:{_PERFUSATE};stroke:none", rx=4)
        c.elements.append(
            f'<path d="M{_f(gx)},{_f(gy)} V{_f(gy + rh - 6)} Q{_f(gx)},{_f(gy + rh)} '
            f"{_f(gx + 6)},{_f(gy + rh)} H{_f(gx + rw - 6)} Q{_f(gx + rw)},{_f(gy + rh)} "
            f'{_f(gx + rw)},{_f(gy + rh - 6)} V{_f(gy)}" fill="none" '
            f'stroke="{palette.ICON_INK}" stroke-width="1.6"/>'
        )
        self._text(rx + 8, ry + rh - 10, f"{_n(circuit.reservoir.value)} ml")
        c.close()
        self.anchors["reservoir"] = self.p(rx, ry + 30)
        # the medium change into and out of the reservoir
        c.open("pka-medium", "medium")
        for x, (y_from, y_to) in (
            (rx + 18, (ry + rh + 34, ry + rh + 2)),
            (rx + rw - 18, (ry + rh + 2, ry + rh + 34)),
        ):
            a, b = self.p(x, y_from), self.p(x, y_to)
            c.elements.append(
                f'<line x1="{_f(a[0])}" y1="{_f(a[1])}" x2="{_f(b[0])}" y2="{_f(b[1])}" '
                f'stroke="{_PIPE}" stroke-width="2" stroke-dasharray="4 3" '
                'marker-end="url(#pka-medium)"/>'
            )
        self._text(rx - 8, ry + rh + 48, "medium change")
        c.close()
        # the pump and the oxygenator, the tubing between them and the reservoir
        c.open("pka-pump", "pump")
        self._pipe([(px, ry + 26), (px, py + r)])
        c.elements.append(_pump(*self.p(px, py), _ICON + 4))
        c.close()
        self.anchors["pump"] = self.p(px - r, py)
        c.open("pka-gas", "gas")
        for x, (y_from, y_to) in (
            (ox - 8, (oy - r - 30, oy - r - 2)),
            (ox + 8, (oy - r - 2, oy - r - 30)),
        ):
            a, b = self.p(x, y_from), self.p(x, y_to)
            c.elements.append(
                f'<line x1="{_f(a[0])}" y1="{_f(a[1])}" x2="{_f(b[0])}" y2="{_f(b[1])}" '
                f'stroke="{_PIPE}" stroke-width="1.6" marker-end="url(#pka-gas)"/>'
            )
        self._text(ox - 14 - _width("O2", _SMALL), oy - r - 14, "O2")
        self._text(ox + 14, oy - r - 14, "CO2")
        c.close()
        c.open("pka-inflow", "inflow")
        self._pipe([(px, py - r), (ox, oy + r)])
        c.elements.append(_oxygenator(*self.p(ox, oy), _ICON + 4))
        x_in = self.left - 34 if self.hilum else inflow[0]
        self._pipe([(ox + r, oy), (x_in, oy), (x_in, inflow[1]), inflow])
        self._text(
            ox + r + 12, oy - 7, f"inflow tubing {_n(circuit.tubing_in.value)} ml"
        )
        mid = (oy + inflow[1]) / 2 + 4 if not self.hilum else inflow[1] - 8
        self._text(x_in - 6 - _width(vein_in, _SMALL), mid, vein_in)
        c.close()
        self.anchors["inflow"] = self.p(ox - r, oy)
        # the outflow back into the reservoir
        c.open("pka-outflow", "outflow")
        x_out = self.left - 18 if self.hilum else outflow[0]
        self._pipe([outflow, (x_out, outflow[1]), (x_out, _RETURN), (rx + rw, _RETURN)])
        mid = (outflow[1] + _RETURN) / 2 + 4
        self._text(x_out - 6 - _width(vein_out, _SMALL), mid, vein_out)
        self._text(
            rx + rw + 10,
            _RETURN + 17,
            f"outflow tubing {_n(circuit.tubing_out.value)} ml",
        )
        c.close()
        # the bile or the urine
        sink = (
            "bile"
            if isinstance(organ, LiverOrgan)
            else "urine"
            if isinstance(organ, KidneyOrgan)
            else ""
        )
        if sink:
            color = palette.BILE if sink == "bile" else palette.URINE
            fill = "#e8f3eb" if sink == "bile" else "#fff8d6"
            top, h, w = _COLLECTOR
            x = excretion[0]
            c.open("pka-sink", sink)
            self._pipe([excretion, (x, top + 16)], color)
            gx, gy = self.p(x - w / 2, top)
            c.rect(gx, gy + 12, w, h - 12, style=f"fill:{fill};stroke:none", rx=3)
            c.elements.append(
                f'<path d="M{_f(gx)},{_f(gy)} V{_f(gy + h)} H{_f(gx + w)} V{_f(gy)}" '
                f'fill="none" stroke="{palette.ICON_INK}" stroke-width="1.4"/>'
            )
            c.close()
            self.anchors[sink] = self.p(x + w / 2, top + h / 2)

    # --- the labels --------------------------------------------------------

    def _label(
        self,
        key: str,
        x: float,
        y: float,
        side: int,
        icon: Callable[[float, float, float], str],
        lines: Sequence[str],
        chips: Sequence[tuple[str, str, bool]] = (),
        zones: bool = False,
    ) -> float:
        """Draw a label in a column from `x`; get its height.

        The icon sits at the side of the drawing: a label of the left
        column (`side` -1) is right-aligned, one of the right column
        left-aligned. The first line is the title.
        """
        c = self.canvas
        right = side < 0
        xe = x + _COLUMN
        t0, t1 = (x, xe - _ICON - 6) if right else (x + _ICON + 6, xe)
        c.open("pka-label", key)
        c.elements.append(
            icon(xe - _ICON / 2 if right else x + _ICON / 2, y + _ICON / 2, _ICON)
        )
        cy = y + 2
        for i, line in enumerate(lines):
            cls = "pkf-b" if i == 0 else "pkf-s pkf-dim"
            size = _SIZE if i == 0 else _SMALL
            for part in _wrap(line, int((t1 - t0) // (size * 0.62))):
                w = _width(part, size)
                c.text(
                    t1 - w if right else t0,
                    cy + _LINE / 2 + 0.35 * size,
                    part,
                    cls,
                    size,
                )
                cy += _LINE
        widths = [_tag_width(q) for q, _, _ in chips]
        for row in _flow(widths, t1 - t0, 5.0) if chips else []:
            total = sum(widths[i] for i in row) + 5.0 * (len(row) - 1)
            cx = t1 - total if right else t0
            for i in row:
                q, color, filled = chips[i]
                c.chip(cx, cy - 1, q, color, tissue=filled, substance=True)
                cx += widths[i] + 5.0
            cy += _LINE
        if zones:
            cy = self._inset(t0, t1, cy + 4)
        c.close()
        return max(cy - y, _ICON + 4)

    def _inset(self, x0: float, x1: float, y: float) -> float:
        """Draw the zones of the organ as a strip in the order of the flow; get the next top."""
        c = self.canvas
        n = self.organ.zones
        first, last = _ends(self.organ)
        what = "sinusoid" if isinstance(self.organ, LiverOrgan) else "zones"
        c.text(
            x0,
            y + _LINE / 2 + 0.35 * _SMALL,
            f"{what}, {_zones(n)}",
            "pkf-s pkf-dim",
            _SMALL,
        )
        y += _LINE
        cell = (x1 - x0) / n
        labelled = cell >= _tag_width(ids.zone_id(n))
        for i in range(1, n + 1):
            zone = ids.zone_id(i)
            grey = round(246 - 40 * (i - 1) / max(n - 1, 1))
            c.open("pkf-zone", f"zone:{zone}")
            c.rect(
                x0 + (i - 1) * cell,
                y,
                cell,
                _LINE - 2,
                style=f"fill:rgb({grey},{grey},{grey});stroke:#8c959f;stroke-width:1",
                rx=0,
            )
            if labelled:
                tx = x0 + (i - 1) * cell + (cell - _width(zone, _SMALL)) / 2
                c.text(tx, y + (_LINE - 2) / 2 + 0.35 * _SMALL, zone, "pkf-s", _SMALL)
            c.close()
        y += _LINE + 2
        if not labelled:
            first, last = f"{first} {ids.zone_id(1)}", f"{ids.zone_id(n)} {last}"
        base = y + _LINE / 2 + 0.35 * _SMALL
        c.text(x0, base, first, "pkf-s pkf-dim", _SMALL)
        lx = x1 - _width(last, _SMALL)
        c.text(lx, base, last, "pkf-s pkf-dim", _SMALL)
        a0, a1 = x0 + _width(first, _SMALL) + 6, lx - 6
        if a1 - a0 >= 16:
            c.arrow(a0, y + _LINE / 2, a1, y + _LINE / 2)
        return y + _LINE

    def _labels(self) -> float:
        """Draw the label columns with their leaders; get the bottom of the columns."""
        d, organ = self.d, self.organ
        pressure, _ = _gas_pressure(d)
        circuit = d.circuit
        c = self.canvas
        left: list[
            tuple[
                str,
                Callable[[float, float, float], str],
                list[str],
                list[tuple[str, str, bool]],
                bool,
            ]
        ] = [
            (
                "inflow",
                _oxygenator,
                [
                    "oxygenator",
                    f"O2 in towards {_n(d.perfusate.o2_target_kpa.value)} kPa, "
                    f"CO2 out towards {_n(pressure)} kPa",
                ],
                [(q, self.colors[q], True) for q in ("o2", "co2")],
                False,
            ),
            (
                "pump",
                _pump,
                ["pump", f"{_n(circuit.flow.value)} ml/min"],
                [],
                False,
            ),
            (
                "reservoir",
                _reservoir_icon,
                [
                    "reservoir",
                    f"{_n(circuit.reservoir.value)} ml of perfusate",
                ],
                [(q, self.colors[q], True) for q in self.quantities],
                False,
            ),
        ]
        key = organ.key

        def organ_icon(cx: float, cy: float, size: float) -> str:
            return icons.organ(key, cx, cy, size) or _tissue(cx, cy, size)

        vein_in, vein_out = _VEINS.get(key, ("artery", "vein"))
        right: list[
            tuple[
                str,
                Callable[[float, float, float], str],
                list[str],
                list[tuple[str, str, bool]],
                bool,
            ]
        ] = [
            (
                key,
                organ_icon,
                [
                    _organ_name(organ),
                    f"{_n(organ.mass.value)} g, {_zones(organ.zones)}, from the "
                    f"{vein_in} to the {vein_out}",
                ],
                [(q, self.colors[q], q in ids.CELL) for q in self.quantities],
                isinstance(organ, LiverOrgan) or organ.zones > 1,
            )
        ]
        if isinstance(organ, LiverOrgan):
            right.append(
                (
                    "bile",
                    _collector,
                    ["bile", "secreted in every zone"],
                    [],
                    False,
                )
            )
        elif isinstance(organ, KidneyOrgan):
            right.append(
                (
                    "urine",
                    _collector,
                    ["urine", "the filtrate less the reabsorbed"],
                    [(q, self.colors[q], True) for q in _urine(d)],
                    False,
                )
            )
        bottom = _TOP + _HD
        out_bottom = bottom
        for side, labels, x in ((-1, left, _PAD), (1, right, self.rx)):
            # measure, place near the anchors, draw
            heights = {}
            for lkey, icon, lines, chips, zones in labels:
                mark = len(c.elements)
                heights[lkey] = self._label(lkey, x, 0, side, icon, lines, chips, zones)
                del c.elements[mark:]
            items = [
                (lkey, self.anchors[lkey][1], heights[lkey]) for lkey, *_ in labels
            ]
            placed = {
                k: top for k, top, _ in AnatomyFigure.place(items, _TOP, bottom, 12.0)
            }
            for lkey, icon, lines, chips, zones in labels:
                top = placed[lkey]
                self._label(lkey, x, top, side, icon, lines, chips, zones)
                ax, ay = self.anchors[lkey]
                edge = x + _COLUMN + 4 if side < 0 else x - 4
                ly = top + _ICON / 2
                c.elements.append(
                    f'<path d="M{_f(ax)},{_f(ay)} L{_f(edge - side * 14)},{_f(ly)} '
                    f'L{_f(edge)},{_f(ly)}" fill="none" stroke="#444" stroke-width="0.8"/>'
                    f'<circle cx="{_f(ax)}" cy="{_f(ay)}" r="2.4" fill="#fff" '
                    'stroke="#222" stroke-width="1"/>'
                )
                out_bottom = max(out_bottom, top + heights[lkey])
            if side < 0:
                # the temperature below the apparatus, without a leader
                top = max(placed[k] + heights[k] for k, *_ in labels) + 16
                text = f"temperature {_n(d.temperature.value)} °C"
                h = self._label(
                    "temperature",
                    x,
                    top,
                    side,
                    _thermometer,
                    [text, "rates scaled by f_T"],
                )
                out_bottom = max(out_bottom, top + h)
        return out_bottom

    def _key(self, y: float) -> float:
        """Draw the key of the lines and the chips from `y`; get the next top."""
        c = self.canvas
        x = _PAD
        items: list[tuple[str, str, str]] = [("line", _PIPE, "perfusate")]
        if isinstance(self.organ, LiverOrgan):
            items.append(("line", palette.BILE, "bile"))
        if isinstance(self.organ, KidneyOrgan):
            items.append(("line", palette.URINE, "urine"))
        items += [
            ("dash", _PIPE, "medium change"),
            ("cells", "", "in the perfusate and the cells"),
            ("perfusate", "", "in the perfusate only"),
        ]
        c.text(x, y + _LINE / 2 + 0.35 * _SIZE, "key", "pkf-b")
        x += _width("key ", _SIZE)
        for kind, color, label in items:
            if kind in ("line", "dash"):
                dash = ' stroke-dasharray="4 3"' if kind == "dash" else ""
                c.elements.append(
                    f'<path d="M{_f(x)},{_f(y + _LINE / 2)} h22" stroke="{color}" '
                    f'stroke-width="3"{dash}/>'
                )
                w = 22.0
            else:
                w = c.chip(x, y - 1, "q", "#57606a", tissue=kind == "cells")
            c.text(x + w + 6, y + _LINE / 2 + 0.35 * _SIZE, label)
            x += w + 6 + _width(label, _SIZE) + 18
        return y + _LINE

    def draw(self) -> str:
        """Draw the figure; get the SVG."""
        c, d = self.canvas, self.d
        organ = self.organ
        c.text(_PAD, _PAD + _TITLE, d.name, "pkf-h", _TITLE)
        c.text(
            _PAD,
            _PAD + _TITLE + 22,
            f"machine perfusion model {d.sid}: {d.species} {_organ_name(organ)} in "
            "its circuit",
            "pkf-dim",
        )
        self._drawing()
        bottom = self._labels()
        height = self._key(bottom + 14) + _PAD
        width = self.width

        def marker(mid: str, color: str, size: float) -> str:
            return (
                f'<marker id="{mid}" viewBox="0 0 10 10" refX="9" refY="5" '
                f'markerUnits="userSpaceOnUse" markerWidth="{size}" markerHeight="{size}" '
                f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>'
            )

        return "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{_f(width)}" '
                f'height="{_f(height)}" viewBox="0 0 {_f(width)} {_f(height)}" '
                'role="img" aria-labelledby="pka-title">',
                f'<title id="pka-title">{escape(d.name)}</title>',
                f"<style>{_STYLE}</style>",
                "<defs>",
                '<marker id="pkf-arrow" viewBox="0 0 10 10" refX="9" refY="5" '
                'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
                '<path d="M0,0 L10,5 L0,10 z" fill="#57606a"/></marker>',
                marker("pka-flow", _PIPE, 11),
                marker("pka-gas", _PIPE, 8),
                marker("pka-medium", _PIPE, 8),
                marker(f"pka-{palette.BILE[1:]}", palette.BILE, 9),
                marker(f"pka-{palette.URINE[1:]}", palette.URINE, 9),
                "</defs>",
                f'<rect x="0" y="0" width="{_f(width)}" height="{_f(height)}" '
                'fill="#ffffff"/>',
                *c.elements,
                "</svg>",
                "",
            ]
        )


def _reservoir_icon(cx: float, cy: float, d: float = _ICON) -> str:
    """The reservoir: a vessel with its perfusate."""
    w, h = d * 0.42, d * 0.48
    x0, y0 = cx - w / 2, cy - h / 2
    return (
        icons.ring(cx, cy, d / 2)
        + f'<rect x="{x0:.1f}" y="{y0 + h * 0.35:.1f}" width="{w:.1f}" '
        f'height="{h * 0.65:.1f}" fill="#dbdbdb"/>'
        f'<path d="M{x0:.1f},{y0:.1f} V{y0 + h:.1f} H{x0 + w:.1f} V{y0:.1f}" '
        f'fill="none" stroke="{palette.ICON_INK}" stroke-width="1"/>'
    )


def _collector(cx: float, cy: float, d: float = _ICON) -> str:
    """A collector of bile or urine: a drop on a circle, in grey."""
    k = d / 22
    return icons.ring(cx, cy, d / 2) + (
        f'<path transform="translate({cx:.1f},{cy:.1f}) scale({k:.3f})" '
        'd="M0,-6.5 C3,-2.5 4.4,0.5 4.4,2.2 A4.4,4.4 0 0 1 -4.4,2.2 '
        f'C-4.4,0.5 -3,-2.5 0,-6.5z" fill="#a3a3a3" stroke="{palette.ICON_INK}" '
        'stroke-width="0.8"/>'
    )


def anatomy_cards(definition: PerfusionDefinition) -> dict[str, list[str]]:
    """Get the hover cards of the anatomical figure by key, each a list of HTML lines.

    Args:
        definition: the definition of the circuit.

    Returns:
        The cards of `figure.circuit_cards` and the card of the medium change
        (`medium`).

    Raises:
        ValueError: for an invalid definition (`PerfusionDefinition.validate`).
    """
    out = circuit_cards(definition)
    out["medium"] = [
        "<b>medium change</b>",
        "an experiment replaces a fraction of the reservoir with fresh perfusate "
        "between the parts of a simulation",
    ]
    return out


def circuit_anatomy_svg(definition: PerfusionDefinition) -> str:
    """Get the SVG of the anatomical figure of a definition (see the module).

    Args:
        definition: the definition of the circuit.

    Returns:
        The SVG, deterministic.

    Raises:
        ValueError: for an invalid definition (`PerfusionDefinition.validate`).
    """
    _check(definition)
    return _Anatomy(definition).draw()


def circuit_anatomy_figure(definition: PerfusionDefinition, path: Path) -> Path:
    """Write the anatomical figure of a definition, the organ in its circuit.

    Args:
        definition: the definition of the circuit.
        path: the SVG file to write.

    Returns:
        The path of the figure.

    Raises:
        ValueError: for an invalid definition (`PerfusionDefinition.validate`).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(circuit_anatomy_svg(definition), encoding="utf-8")
    return path
