"""The figure of a machine perfusion: an SVG schematic of a circuit and its organ.

`circuit_figure(definition, path)` draws a `PerfusionDefinition` in the
language of the figures of the PBPK library (`pkpdsim.pbpk.figure`: rows
with an icon, a title, chips and process badges, monospace text placed by
its width at `CHAR_WIDTH`, `data-k` keys for the hover cards of an
interactive page), as the COMBINE archive of a model carries it
(`pkpdsim.perfusion.archive.create_archive`):

- the reservoir as a column on the left with a chip of every quantity of
  the perfusate (`circuit.quantities`, the tracers included) and its volume;
- the circuit as rows in series, from the reservoir to the pump (its flow
  `Q_pump`), the inflow tubing with the oxygenator (its volume, the oxygen
  towards its target and the CO2 towards the pCO2 of the sweep gas, the gas
  as a compartment with O2 in and CO2 out), the organ and the outflow tubing
  back to the reservoir, the flow as arrows between them;
- the organ with its grey icon (`icons.organ`, the generic tissue on a
  circle for an organ without an icon), its token, mass and zones; the zones
  as a strip in the order of the flow, from periportal to pericentral in the
  liver, from arterial to venous in the kidney and a generic organ, each
  zone with its hover card (`zone:<zone>`); a chip per quantity of the
  perfusate, filled where the cells hold it as well (`ids.CELL`), dashed
  where it stays in the perfusate (CO2, bicarbonate, the protons, the
  tracers); the processes of the cellular core of every zone as badges, read
  from the composed model (`process_groups`): transport (`TRA`), metabolism
  (`MET`), the redox (`ROS`), the necrosis and the release (`NEC`), the
  zonal profiles (`ZON`), the bile (`BIL`) or the filtration into the urine
  (`URI`); the bile of the liver or the urine of the kidney (with the
  filtered quantities) as a compartment at the right end of the row;
- the temperature in the band above the rows, the legend (the title, the
  quantities of the perfusate with their fresh composition, the key of the
  chips and badges) in a column right of the diagram.

The figure is pure python and deterministic: the same definition gives the
same text. The anatomical figure of the organ in its circuit (panel A) is
`pkpdsim.perfusion.anatomy`. `interactive_circuit` writes the
self-contained page of either figure with the hover cards of `circuit_cards` (volumes, flow, temperature, the processes of
the core per zone, the released enzymes, the tracers) through
`pkpdsim.pbpk.visual.interactive.html_page`.
"""

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from xml.sax.saxutils import escape

from pkpdsim.pbpk.figure import (
    _BADGES,
    _COLORS,
    _ICON,
    _PAD,
    _SIZE,
    _SMALL,
    _STYLE,
    _TITLE,
    LEGEND_WIDTH,
    MAX_TEXT,
    _Canvas,
    _chars,
    _f,
    _flow,
    _tag_width,
    _width,
    _wrap,
)
from pkpdsim.pbpk.visual import icons, palette
from pkpdsim.pbpk.visual.anatomy import TITLES
from pkpdsim.perfusion import core, ids
from pkpdsim.perfusion.builder import compose, zone_factors
from pkpdsim.perfusion.circuit import fresh_co2_pressure, quantities, quantity_label
from pkpdsim.perfusion.definition import (
    KidneyOrgan,
    LiverOrgan,
    Organ,
    PerfusionDefinition,
    table_value,
)

__all__ = [
    "BADGES",
    "circuit_cards",
    "circuit_figure",
    "circuit_svg",
    "interactive_circuit",
    "process_groups",
]

BADGES: dict[str, tuple[str, str, str]] = {
    "FLO": ("pump flow", "#eef6ff", "#0969da"),
    "OXY": ("oxygenator", "#e7f7f9", "#0e7c86"),
    "TRA": _BADGES["TRA"],
    "MET": _BADGES["MET"],
    "ROS": ("redox", "#fff0f6", "#bf3989"),
    "NEC": ("necrosis and release", "#fbe9e9", "#a40e26"),
    "ZON": ("zonal profile", "#f3f0ff", "#6639ba"),
    "BIL": _BADGES["BIL"],
    "URI": _BADGES["URI"],
}
"""The badges of the circuit and of the core, in the order of the key: the
label of the key, fill and stroke (those of the PBPK library for its codes)."""

_GROUPS: dict[str, tuple[str, ...]] = {
    "TRA": ("GLCIM", "O2IM", "LACEX"),
    "MET": ("GLGP", "GLC2PYR", "LDH", "ATPU", "TCA", "OXPHOS", "PPP"),
    "ROS": ("ROSP", "ROSDETOX", "GR"),
}
"""The reactions of the cells of a zone by the badge of their line."""

_COLUMN = 76
"""The smallest width of the reservoir column (px)."""
_ARROW = 40
"""The length of the arrows between the reservoir and the rows (px)."""
_GAP = 24
"""The gap between two rows, the arrows of the flow (px)."""
_STRIP = 240
"""The width of the strip of the zones (px), wider for many zones."""
_GREY = "#a3a3a3"
"""The fill of the drawn icons, a grey of `palette.ORGAN_GREYS`."""
_LIGHT = "#dbdbdb"
"""The light fill of the drawn icons."""
_SINKS = {
    "gas": "fill:#eef6fb;stroke:#4a90a4;stroke-width:1.4",
    "bile": f"fill:#e8f3eb;stroke:{palette.BILE};stroke-width:1.4",
    "urine": "fill:#fff8d6;stroke:#b08a00;stroke-width:1.4",
}
"""The style of a compartment of a row."""


_TRACER_COLORS = ("#17becf", "#bcbd22", "#393b79", "#e377c2", "#8c564b")
"""The colours of the tracers, after the ten colours of the quantities of the
perfusate (`_COLORS` of the PBPK figure), so a tracer has a colour of its own."""


def _colors(quantities: Sequence[str], tracers: Sequence[str]) -> dict[str, str]:
    """The colour of every quantity of the perfusate and of every tracer.

    The quantities take the colours of the substances of the PBPK figure in
    their order, the tracers those of `_TRACER_COLORS`; both repeat beyond
    their number.
    """
    plain = [q for q in quantities if q not in tracers]
    out = {q: _COLORS[i % len(_COLORS)] for i, q in enumerate(plain)}
    out |= {t: _TRACER_COLORS[i % len(_TRACER_COLORS)] for i, t in enumerate(tracers)}
    return out


def _n(x: float) -> str:
    """A number of the figure, four significant digits."""
    return f"{x:.4g}"


# --- the icons of the circuit, on the circle of the organ icons ------------


def _pump(cx: float, cy: float, d: float = _ICON) -> str:
    """The pump: a circle with the triangle of its flow."""
    r = d * 0.3
    tri = (
        f"M{cx - r * 0.45:.1f},{cy - r * 0.7:.1f} L{cx + r * 0.8:.1f},{cy:.1f} "
        f"L{cx - r * 0.45:.1f},{cy + r * 0.7:.1f} Z"
    )
    return (
        icons.ring(cx, cy, d / 2)
        + f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{_LIGHT}" '
        f'stroke="{palette.ICON_INK}" stroke-width="1"/>'
        f'<path d="{tri}" fill="{_GREY}" stroke="{palette.ICON_INK}" stroke-width="0.8"/>'
    )


def _oxygenator(cx: float, cy: float, d: float = _ICON) -> str:
    """The membrane oxygenator: a module with its hollow fibres."""
    w, h = d * 0.56, d * 0.38
    fibres = "".join(
        f'<path d="M{cx - w / 2 + 2:.1f},{cy + dy:.1f} h{w - 4:.1f}" '
        f'stroke="{palette.ICON_INK}" stroke-width="0.7"/>'
        for dy in (-h / 4, 0.0, h / 4)
    )
    return (
        icons.ring(cx, cy, d / 2)
        + f'<rect x="{cx - w / 2:.1f}" y="{cy - h / 2:.1f}" width="{w:.1f}" '
        f'height="{h:.1f}" rx="1.5" fill="{_LIGHT}" stroke="{palette.ICON_INK}" '
        'stroke-width="1"/>' + fibres
    )


def _tubing(cx: float, cy: float, d: float = _ICON) -> str:
    """The tubing: a bent tube."""
    a = d * 0.3
    tube = (
        f"M{cx - a:.1f},{cy + a * 0.5:.1f} C{cx - a * 0.2:.1f},{cy + a * 0.5:.1f} "
        f"{cx - a * 0.2:.1f},{cy - a * 0.5:.1f} {cx + a * 0.1:.1f},{cy - a * 0.5:.1f} "
        f"L{cx + a:.1f},{cy - a * 0.5:.1f}"
    )
    return (
        icons.ring(cx, cy, d / 2)
        + f'<path d="{tube}" fill="none" stroke="{palette.ICON_INK}" stroke-width="4" '
        'stroke-linecap="round"/>'
        f'<path d="{tube}" fill="none" stroke="{_LIGHT}" stroke-width="2.2" '
        'stroke-linecap="round"/>'
    )


def _thermometer(cx: float, cy: float, d: float = _ICON) -> str:
    """The temperature: a thermometer."""
    h = d * 0.5
    return (
        icons.ring(cx, cy, d / 2)
        + f'<rect x="{cx - 1.6:.1f}" y="{cy - h / 2 - 1:.1f}" width="3.2" '
        f'height="{h:.1f}" rx="1.6" fill="{_LIGHT}" stroke="{palette.ICON_INK}" '
        'stroke-width="0.9"/>'
        f'<circle cx="{cx:.1f}" cy="{cy + h / 2:.1f}" r="2.9" fill="{_GREY}" '
        f'stroke="{palette.ICON_INK}" stroke-width="0.9"/>'
    )


def _tissue(cx: float, cy: float, d: float = _ICON) -> str:
    """An organ without an icon of the anatomy: its cells."""
    r = d * 0.13
    cells = "".join(
        f'<circle cx="{cx + dx * d:.1f}" cy="{cy + dy * d:.1f}" r="{r:.1f}" '
        f'fill="{_GREY}" stroke="{palette.ICON_INK}" stroke-width="0.8"/>'
        for dx, dy in ((-0.14, -0.1), (0.14, -0.1), (0.0, 0.14))
    )
    return icons.ring(cx, cy, d / 2) + cells


# --- the content -----------------------------------------------------------


def process_groups(definition: PerfusionDefinition) -> dict[str, list[tuple[str, str]]]:
    """Get the reactions of the cells of a zone by the badge of their line.

    The reactions are those of the first zone of the composed model
    (`builder.compose`), so the figure shows what the model has.

    Args:
        definition: the definition.

    Returns:
        The id (without the zone) and the name (without the zone) of every
        reaction of the cells by badge: `TRA`, `MET`, `ROS`, `NEC` (the
        death of the cells and the release of every released quantity) and,
        for a kidney, `URI` (the filtration of every filtered quantity).
    """
    zone = ids.zone_id(1)
    names: dict[str, str] = {}
    for r in compose(definition).reactions or []:
        sid = r.sid or ""
        if sid.endswith(f"_{zone}"):
            name = (r.name or sid).removesuffix(f" ({zone})")
            names[sid.removesuffix(f"_{zone}")] = name
    out: dict[str, list[tuple[str, str]]] = {
        code: [(sid, names[sid]) for sid in sids if sid in names]
        for code, sids in _GROUPS.items()
    }
    out["NEC"] = [
        (sid, names[sid])
        for q in definition.released
        for sid in (f"{q.upper()}D", f"{q.upper()}R")
        if sid in names
    ]
    filtration = [
        ids.filtration(q, zone).removesuffix(f"_{zone}") for q in _urine(definition)
    ]
    out["URI"] = [(sid, names[sid]) for sid in filtration if sid in names]
    return {code: rows for code, rows in out.items() if rows}


def _urine(definition: PerfusionDefinition) -> list[str]:
    """The quantities which a kidney filters into the urine, the tracers last."""
    if not isinstance(definition.organ, KidneyOrgan):
        return []
    return [*ids.FILTERED, *definition.perfusate.tracers]


def _organ_name(organ: Organ) -> str:
    """The title of an organ: its name and token, `organ <token>` without a name."""
    name = TITLES.get(organ.key)
    return f"{name} {organ.key}" if name else f"organ {organ.key}"


def _ends(organ: Organ) -> tuple[str, str]:
    """The ends of the zones in the order of the flow."""
    if isinstance(organ, LiverOrgan):
        return "periportal", "pericentral"
    return "arterial", "venous"


def _zones(n: int) -> str:
    """The number of zones in words."""
    return "1 zone" if n == 1 else f"{n} zones"


def _fresh(definition: PerfusionDefinition, q: str) -> str:
    """The fresh perfusate of a quantity in words."""
    p = definition.perfusate
    plain = {"glc": p.glucose, "lac": p.lactate, "hco3": p.hco3, "k": p.k}
    if q in plain:
        text = f"fresh {_n(plain[q].value)} mM"
    elif q == "o2":
        text = f"fresh at pO2 {_n(p.o2_target_kpa.value)} kPa"
    elif q == "co2":
        text = f"fresh at the buffer equilibrium, pCO2 {_n(fresh_co2_pressure(p).value)} kPa"
    elif q == "h":
        text = f"fresh at pH {_n(p.ph.value)}"
    elif q in ids.ENZYMES:
        text = "none fresh"
    else:
        text = "inert, none fresh"
        if isinstance(definition.organ, KidneyOrgan):
            text += ", filtered and not reabsorbed"
    if q in definition.released:
        text += ", released by necrosis"
    return text


def _gas_pressure(definition: PerfusionDefinition) -> tuple[float, str]:
    """The pCO2 of the sweep gas (kPa) and where it comes from."""
    gas = definition.circuit.pco2_gas
    if gas is not None:
        return gas.value, "the sweep gas"
    return fresh_co2_pressure(definition.perfusate).value, "the fresh perfusate"


def _co2_rate(definition: PerfusionDefinition) -> float:
    """The rate constant of the exchange of CO2 (1/min)."""
    c = definition.circuit
    return (c.oxygenator_k_co2 or c.oxygenator_k).value


@dataclass
class _Sink:
    """A compartment at the right end of a row: the gas, the bile or the urine.

    Attributes:
        key: the key of its hover card and of its style.
        name: its title.
        chips: its chips (label, colour, filled).
        arrows: the arrows between the row and the compartment, by the line
            of the compartment (0 its title, `j + 1` its chip `j`), inward
            (from the compartment into the row) or outward.
    """

    key: str
    name: str
    chips: list[tuple[str, str, bool]]
    arrows: list[tuple[int, bool]]


@dataclass
class _Row:
    """A row of the circuit: its icon, title, zones, chips, processes and compartment."""

    key: str
    icon: Callable[[float, float, float], str]
    title: list[str]
    lines: list[tuple[str, str]] = field(default_factory=list)
    chips: list[tuple[str, str, bool]] = field(default_factory=list)
    zones: int = 0
    ends: tuple[str, str] = ("", "")
    sink: _Sink | None = None


class _Circuit:
    """The layout of the figure of a definition."""

    def __init__(self, definition: PerfusionDefinition, width: float) -> None:
        """Collect what the figure shows; `width` is the width of the diagram."""
        self.d = definition
        self.organ = definition.organ.resolved(definition.species)
        self.width = width
        self.quantities = list(quantities(definition.perfusate, definition.released))
        self.colors = _colors(self.quantities, definition.perfusate.tracers)
        self.canvas = _Canvas(badges=BADGES)
        self.k = self.canvas.k
        self.maxx = 0.0
        self.codes: set[str] = set()

    # --- the rows ----------------------------------------------------------

    def _rows(self) -> list[_Row]:
        """The rows in the order of the flow: pump, inflow tubing, organ, outflow tubing."""
        d, c, organ = self.d, self.d.circuit, self.organ
        flow = c.flow.value
        pressure, of = _gas_pressure(d)
        pump = _Row(
            "pump",
            _pump,
            ["pump", f"{_n(flow)} ml/min"],
            [
                (
                    "FLO",
                    f"Q_pump {_n(flow)} ml/min, {_n(flow / organ.mass.value)} ml/min/g "
                    "of the organ",
                )
            ],
        )
        inflow = _Row(
            "inflow",
            _oxygenator,
            ["inflow tubing", f"{_n(c.tubing_in.value)} ml", "oxygenator"],
            [
                (
                    "OXY",
                    f"O2 towards pO2 {_n(d.perfusate.o2_target_kpa.value)} kPa, k_oxy "
                    f"{_n(c.oxygenator_k.value)} 1/min",
                ),
                (
                    "OXY",
                    f"CO2 towards pCO2 {_n(pressure)} kPa of {of}, k_oxy_co2 "
                    f"{_n(_co2_rate(d))} 1/min",
                ),
            ],
            sink=_Sink(
                "gas",
                "sweep gas",
                [(q, self.colors[q], True) for q in ("o2", "co2")],
                [(1, True), (2, False)],
            ),
        )
        outflow = _Row(
            "outflow", _tubing, ["outflow tubing", f"{_n(c.tubing_out.value)} ml"]
        )
        if isinstance(organ, KidneyOrgan):
            outflow.lines.append(("FLO", "the pump flow less the urine"))
        return [pump, inflow, self._organ_row(), outflow]

    def _organ_row(self) -> _Row:
        """The row of the organ: its zones, chips, processes and bile or urine."""
        d, organ = self.d, self.organ
        groups = process_groups(d)
        lines = [
            (code, ", ".join(name for _, name in groups[code]))
            for code in ("TRA", "MET", "ROS")
            if code in groups
        ]
        released = [core.quantity_name(q) for q in d.released]
        if released:
            lines.append(("NEC", f"necrosis releases {_and(released)}"))
        factors = zone_factors(organ)
        n = organ.zones
        first, last = ids.zone_id(1), ids.zone_id(n)
        profiles = [
            f"{p} {_n(factors[first][p])} ({first}) to {_n(factors[last][p])} ({last})"
            for p in organ.profiles
        ]
        if profiles:
            lines.append(("ZON", "; ".join(profiles)))
        sink = None
        if isinstance(organ, LiverOrgan):
            lines.append(("BIL", "bile secretion of the hepatocytes of every zone"))
            sink = _Sink("bile", "bile", [], [(0, False)])
        elif isinstance(organ, KidneyOrgan):
            assert organ.gfr is not None and organ.ff_max is not None
            assert organ.reabsorption is not None
            lines.append(
                (
                    "URI",
                    f"filtration {_n(organ.gfr.value)} ml/min/g, at most "
                    f"{_n(organ.ff_max.value)} of the pump flow; "
                    f"{_n(organ.reabsorption.value)} of the water and small solutes "
                    "reabsorbed, no tracer",
                )
            )
            sink = _Sink(
                "urine",
                "urine",
                [(q, self.colors[q], True) for q in _urine(d)],
                [(0, False)],
            )
        key = organ.key

        def icon(cx: float, cy: float, size: float) -> str:
            return icons.organ(key, cx, cy, size) or _tissue(cx, cy, size)

        return _Row(
            key,
            icon,
            [_organ_name(organ), f"{_n(organ.mass.value)} g", _zones(n)],
            lines,
            [(q, self.colors[q], q in ids.CELL) for q in self.quantities],
            n,
            _ends(organ),
            sink,
        )

    # --- the drawing -------------------------------------------------------

    def _strip(self, row: _Row, x: float, top: float) -> float:
        """Draw the zones of the organ as a strip in the order of the flow; get its right edge."""
        c, k = self.canvas, self.k
        first, last = row.ends
        n = row.zones
        cell = max(_STRIP / n, 8.0)
        labelled = cell >= _tag_width(ids.zone_id(n))
        if not labelled:
            first, last = f"{first} {ids.zone_id(1)}", f"{ids.zone_id(n)} {last}"
        c.text(x, k.baseline(top, _SIZE), first, "pkf-dim")
        sx = x + _width(first, _SIZE) + 8
        for i in range(1, n + 1):
            zone = ids.zone_id(i)
            c.open("pkf-zone", f"zone:{zone}")
            # lighter at the inflow, darker at the outflow of the organ
            grey = round(246 - 40 * (i - 1) / max(n - 1, 1))
            c.rect(
                sx + (i - 1) * cell,
                top + 2,
                cell,
                k.line - 4,
                style=f"fill:rgb({grey},{grey},{grey});stroke:#8c959f;stroke-width:1",
                rx=0,
            )
            if labelled:
                tx = sx + (i - 1) * cell + (cell - _width(zone, _SMALL)) / 2
                c.text(tx, k.baseline(top, _SMALL), zone, "pkf-s", _SMALL)
            c.close()
        end = sx + n * cell
        c.arrow(end + 4, top + k.line / 2, end + 20, top + k.line / 2)
        c.text(end + 26, k.baseline(top, _SIZE), last, "pkf-dim")
        return end + 26 + _width(last, _SIZE)

    def _draw_row(self, row: _Row, x0: float, y: float, x1: float, cx: float) -> float:
        """Draw a row from `x0` to `x1`, its content from `cx`; get its height.

        The box of the row is drawn first, behind its labels, with the
        height of its content, as a row of the PBPK figure.
        """
        c, k = self.canvas, self.k
        top = y + k.vpad
        c.open("pkf-box", row.key)
        box = len(c.elements)
        c.elements.append(row.icon(x0 + 7 + _ICON / 2, top + _ICON / 2 - 1, _ICON))
        ix = _ICON + 6
        t = top
        for i, line in enumerate(row.title):
            c.text(
                x0 + 10 + ix,
                k.baseline(t, _SIZE),
                line,
                "pkf-b" if i == 0 else "pkf-dim",
            )
            t += k.line
        right = min(x1 - 10, cx + MAX_TEXT)
        sink: tuple[float, float] | None = None
        if row.sink is not None:
            sw = (
                max(
                    [_width(row.sink.name, _SIZE)]
                    + [_tag_width(label) for label, _, _ in row.sink.chips]
                )
                + 16
            )
            sink = (x1 - 10 - sw, sw)
            right = min(right, x1 - 10 - sw - 26)
        rowmax = cx
        cy = top
        if row.zones:
            rowmax = max(rowmax, self._strip(row, cx, cy))
            cy += k.line
        widths = [_tag_width(label) for label, _, _ in row.chips]
        for line in _flow(widths, right - cx, 6.0) if row.chips else []:
            x = cx
            for i in line:
                label, color, filled = row.chips[i]
                c.chip(x, cy, label, color, tissue=filled, substance=True)
                rowmax = max(rowmax, x + widths[i])
                x += widths[i] + 6
            cy += k.line
        for code, text in row.lines:
            self.codes.add(code)
            tx = cx + _tag_width(code) + 6
            c.open("pkf-proc", f"txt:{code}: {text}")
            c.badge(cx, cy, code)
            for part in _wrap(text, _chars(right - tx, _SIZE)):
                c.text(tx, k.baseline(cy, _SIZE), part)
                rowmax = max(rowmax, tx + _width(part, _SIZE))
                cy += k.line
            c.close()
        if sink is not None and row.sink is not None:
            sx, sw = sink
            s = row.sink
            sh = (1 + len(s.chips)) * k.line + 6
            c.open("pkf-sink", s.key)
            c.rect(sx, top, sw, sh, style=_SINKS[s.key], rx=6)
            c.text(sx + 8, k.baseline(top + 2, _SIZE), s.name, "pkf-b")
            for j, (label, color, filled) in enumerate(s.chips):
                c.chip(
                    sx + 8,
                    top + 2 + (j + 1) * k.line,
                    label,
                    color,
                    tissue=filled,
                    substance=True,
                )
            c.close()
            for j, inward in s.arrows:
                ay = top + 2 + j * k.line + k.line / 2
                if inward:
                    c.arrow(sx - 2, ay, sx - 22, ay)
                else:
                    c.arrow(sx - 22, ay, sx - 2, ay)
            cy = max(cy, top + sh)
            rowmax += sw + 26
        self.maxx = max(self.maxx, rowmax)
        height = max(t, cy, top + _ICON) - y + k.vpad
        content = c.elements[box:]
        del c.elements[box:]
        c.rect(x0, y, x1 - x0, height, "pkf-org")
        c.elements += content
        c.close()
        return height

    def _paragraph(
        self, x: float, y: float, right: float, text: str, cls: str = ""
    ) -> float:
        """Draw a wrapped text from `x` to `right`; get the next top."""
        for line in _wrap(text, _chars(right - x, _SIZE)):
            self.canvas.text(x, self.k.baseline(y, _SIZE), line, cls)
            y += self.k.line
        return y

    def _legend(self, x: float, y: float) -> float:
        """Draw the title, the quantities of the perfusate and the key from `x`; get the next top."""
        c, d, k = self.canvas, self.d, self.k
        right = x + LEGEND_WIDTH - 2 * _PAD
        # a long name wraps within the legend
        for line in _wrap(d.name, _chars(right - x, _TITLE)):
            c.text(x, y + _TITLE, line, "pkf-h", _TITLE)
            y += _TITLE + 4
        y += 6
        organ = self.organ
        y = self._paragraph(x, y, right, f"machine perfusion model {d.sid}", "pkf-dim")
        y = self._paragraph(
            x,
            y,
            right,
            f"{d.species} {_organ_name(organ)}, {_zones(organ.zones)}, "
            f"{_n(d.temperature.value)} °C",
            "pkf-dim",
        )
        y += 6
        c.text(x, k.baseline(y, _SIZE), "perfusate", "pkf-b")
        y += k.line
        token_w = max(_width(q, _SIZE) for q in self.quantities)
        x_text = x + 16 + token_w + 12
        for q in self.quantities:
            c.dot(x + 5, y + k.line / 2, self.colors[q])
            c.text(x + 16, k.baseline(y, _SIZE), q, "pkf-b")
            y = self._paragraph(
                x_text, y, right, f"{quantity_label(q)}, {_fresh(d, q)}"
            )
        y = self._paragraph(
            x,
            y + 4,
            right,
            "The bicarbonate buffer acts at every site; every quantity is "
            "observable at every site of the circuit and in every zone.",
        )
        y += 6
        items: list[tuple[str, str]] = [
            ("cells", "in the perfusate and the cells"),
            ("perfusate", "in the perfusate only"),
        ]
        items += [(code, BADGES[code][0]) for code in BADGES if code in self.codes]
        widths = []
        for kind, label in items:
            head = (
                _tag_width("q") if kind in ("cells", "perfusate") else _tag_width(kind)
            )
            widths.append(head + 6 + _width(label, _SIZE))
        lines = _flow(widths, right - x - _width("key ", _SIZE), 18)
        c.text(x, k.baseline(y, _SIZE), "key", "pkf-b")
        x0 = x + _width("key ", _SIZE)
        for line in lines:
            xi = x0
            for i in line:
                kind, label = items[i]
                if kind in ("cells", "perfusate"):
                    w = c.chip(xi, y, "q", "#57606a", tissue=kind == "cells")
                else:
                    w = c.badge(xi, y, kind)
                c.text(xi + w + 6, k.baseline(y, _SIZE), label)
                xi += widths[i] + 18
            y += k.line
        return y

    def draw(self) -> str:
        """Draw the figure; get the SVG."""
        c, k, d = self.canvas, self.k, self.d
        rows = self._rows()
        res_title = ["reservoir", f"{_n(d.circuit.reservoir.value)} ml"]
        column = (
            max(
                [_COLUMN - 16]
                + [_width(t, _SIZE) for t in res_title]
                + [_tag_width(q) for q in self.quantities]
            )
            + 16
        )
        ve0, ve1 = _PAD, _PAD + column
        mx0, mx1 = ve1 + _ARROW, self.width - _PAD
        cx = mx0 + max(_width(t, _SIZE) for r in rows for t in r.title) + 30 + _ICON + 6
        body_top = _PAD
        rows_top = body_top + 2 * k.line + 14
        # the temperature in the band above the rows
        c.open("pkf-temp", "temperature")
        mid = body_top + 4 + k.line / 2
        c.elements.append(_thermometer(mx0 + 7 + _ICON / 2, mid, 22))
        text = f"temperature {_n(d.temperature.value)} °C, the rates scaled by f_T"
        c.text(mx0 + 10 + _ICON + 6, k.baseline(body_top + 4, _SIZE), text)
        c.close()
        self.maxx = max(self.maxx, mx0 + 10 + _ICON + 6 + _width(text, _SIZE))
        behind = len(c.elements)
        top = rows_top
        icon_x = mx0 + 7 + _ICON / 2
        bottom = top
        for i, row in enumerate(rows):
            if i:
                c.arrow(icon_x, bottom, icon_x, top)
            height = self._draw_row(row, mx0, top, mx1, cx)
            bottom = top + height
            mid = top + k.vpad + k.line / 2
            if i == 0:
                c.arrow(ve1, mid, mx0, mid)
            elif i == len(rows) - 1:
                c.arrow(mx0, mid, ve1, mid)
            top = bottom + _GAP
        # the reservoir column, behind the arrows
        chips_bottom = body_top + 10 + (2 + len(self.quantities)) * k.line + 10
        res = _Canvas(k, badges=BADGES)
        res.open("pkf-box", "reservoir")
        res.rect(ve0, body_top, column, max(bottom, chips_bottom) - body_top, "pkf-res")
        for j, t in enumerate(res_title):
            res.text(
                ve0 + 8,
                k.baseline(body_top + 4 + j * k.line, _SIZE),
                t,
                "pkf-b" if j == 0 else "pkf-dim",
            )
        for j, q in enumerate(self.quantities):
            res.chip(
                ve0 + 8,
                body_top + 10 + (2 + j) * k.line,
                q,
                self.colors[q],
                tissue=True,
                substance=True,
            )
        res.close()
        c.elements[behind:behind] = res.elements
        diagram_bottom = max(bottom, chips_bottom)
        width = self.width - 4 + LEGEND_WIDTH
        legend_bottom = self._legend(self.width - 4 + _PAD, _PAD)
        height = max(diagram_bottom, legend_bottom) + _PAD
        style = _STYLE + ".pkf-res{fill:#eef3f8;stroke:#57606a;stroke-width:1.2}\n"
        return "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{_f(width)}" '
                f'height="{_f(height)}" viewBox="0 0 {_f(width)} {_f(height)}" '
                'role="img" aria-labelledby="pkf-title">',
                f'<title id="pkf-title">{escape(d.name)}</title>',
                f"<style>{style}</style>",
                "<defs>",
                '<marker id="pkf-arrow" viewBox="0 0 10 10" refX="9" refY="5" '
                'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
                '<path d="M0,0 L10,5 L0,10 z" fill="#57606a"/></marker>',
                "</defs>",
                f'<rect x="0" y="0" width="{_f(width)}" height="{_f(height)}" '
                'fill="#ffffff"/>',
                *c.elements,
                "</svg>",
                "",
            ]
        )


def _check(definition: PerfusionDefinition) -> None:
    """Check that a definition builds, the figure shows only what builds.

    Raises:
        ValueError: with the problems of an invalid definition.
    """
    problems = definition.validate()
    if problems:
        raise ValueError(
            f"The definition '{definition.sid}' is invalid, no figure: "
            + "; ".join(problems)
        )


def _and(items: Sequence[str]) -> str:
    """A list in words: `a, b and c`."""
    if len(items) <= 1:
        return "".join(items)
    return f"{', '.join(items[:-1])} and {items[-1]}"


def circuit_svg(definition: PerfusionDefinition) -> str:
    """Get the SVG of the figure of a definition (see the module).

    The figure is drawn twice: once wide to measure its content, then at
    the width the content needs.

    Args:
        definition: the definition of the circuit.

    Returns:
        The SVG, deterministic.

    Raises:
        ValueError: for an invalid definition (`PerfusionDefinition.validate`).
    """
    _check(definition)
    probe = _Circuit(definition, 2400)
    probe.draw()
    return _Circuit(definition, math.ceil(probe.maxx + 10 + _PAD)).draw()


def circuit_figure(definition: PerfusionDefinition, path: Path) -> Path:
    """Write the figure of a definition, an SVG schematic of the circuit.

    Args:
        definition: the definition of the circuit.
        path: the SVG file to write.

    Returns:
        The path of the figure.

    Raises:
        ValueError: for an invalid definition (`PerfusionDefinition.validate`).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(circuit_svg(definition), encoding="utf-8")
    return path


def circuit_cards(definition: PerfusionDefinition) -> dict[str, list[str]]:
    """Get the hover cards of the figure of a definition by key, each a list of HTML lines.

    Args:
        definition: the definition of the circuit.

    Returns:
        The cards of the reservoir, the pump, the inflow tubing and its
        sweep gas, the organ, every zone (`zone:<zone>`), the outflow
        tubing, every quantity of the perfusate (`sub:<q>`), the bile or the
        urine and the temperature.

    Raises:
        ValueError: for an invalid definition (`PerfusionDefinition.validate`).
    """
    _check(definition)
    d, c = definition, definition.circuit
    organ = d.organ.resolved(d.species)
    assert organ.density is not None and organ.vascular_fraction is not None
    assert organ.cell_fraction is not None
    qs = list(quantities(d.perfusate, d.released))
    kidney = isinstance(organ, KidneyOrgan)
    flow, mass = c.flow.value, organ.mass.value
    volume = mass / organ.density.value
    f_vas, f_cell = organ.vascular_fraction.value, organ.cell_fraction.value
    n = organ.zones
    out: dict[str, list[str]] = {}
    out["reservoir"] = [
        f"<b>reservoir</b> <code>{ids.volume(ids.RESERVOIR)}</code>",
        f"volume {_n(c.reservoir.value)} ml",
        "the perfusate: " + ", ".join(escape(q) for q in qs),
        "a medium change of an experiment replaces a fraction of it with fresh "
        "perfusate",
    ]
    if kidney:
        out["reservoir"].append("loses the urine, its volume a rate rule")
    out["pump"] = [
        f"<b>pump</b> <code>{ids.Q_PUMP}</code>",
        f"flow {_n(flow)} ml/min, {_n(flow / mass)} ml/min/g of the organ",
        "the total flow of the circuit; an experiment sets the flow of its group",
    ]
    pressure, of = _gas_pressure(d)
    out["inflow"] = [
        f"<b>inflow tubing with the oxygenator</b> <code>{ids.volume(ids.INFLOW)}</code>",
        f"volume {_n(c.tubing_in.value)} ml, residence time "
        f"{_n(c.tubing_in.value / flow) if flow else 'infinite'} min",
        f"O2 towards pO2 {_n(d.perfusate.o2_target_kpa.value)} kPa "
        f"(<code>{ids.O2_TARGET}</code>), <code>{ids.OXYGENATOR_K}</code> "
        f"{_n(c.oxygenator_k.value)} 1/min",
        f"CO2 towards pCO2 {_n(pressure)} kPa of {of} (<code>{ids.CO2_GAS}</code>), "
        f"<code>{ids.OXYGENATOR_K_CO2}</code> {_n(_co2_rate(d))} 1/min",
    ]
    out["gas"] = [
        "<b>sweep gas</b> of the oxygenator",
        f"O2 into the perfusate, towards {_n(d.perfusate.o2_target_kpa.value)} kPa",
        f"CO2 out of the perfusate, towards {_n(pressure)} kPa",
    ]
    out["outflow"] = [
        f"<b>outflow tubing</b> <code>{ids.volume(ids.OUTFLOW)}</code>",
        f"volume {_n(c.tubing_out.value)} ml",
    ]
    if kidney:
        out["outflow"].append("flow: the pump flow less the urine")
    groups = process_groups(d)
    core_names = [
        name for code in ("TRA", "MET", "ROS") for _, name in groups.get(code, [])
    ]
    first, last = _ends(organ)
    rows = [
        f"<b>{escape(_organ_name(organ))}</b> <code>{ids.mass(organ.key)}</code>",
        f"mass {_n(mass)} g, density {_n(organ.density.value)} kg/l, volume "
        f"{_n(volume)} ml",
        f"vascular space {_n(100 * f_vas)} % ({_n(f_vas * volume)} ml), cells "
        f"{_n(100 * f_cell)} % ({_n(f_cell * volume)} ml)",
        f"{_zones(n)} in series from {first} to {last}, each 1/{n} of the volumes",
        "the core in the cells of every zone: " + escape(", ".join(core_names)),
        "released by necrosis: "
        + (", ".join(escape(core.quantity_name(q)) for q in d.released) or "nothing"),
    ]
    if isinstance(organ, LiverOrgan):
        rows.append("secretes bile")
    if isinstance(organ, KidneyOrgan):
        assert organ.gfr is not None and organ.ff_max is not None
        assert organ.reabsorption is not None
        rows.append(
            f"GFR {_n(organ.gfr.value)} ml/min/g ({_n(organ.gfr.value * mass)} ml/min), "
            f"at most {_n(organ.ff_max.value)} of the pump flow; reabsorbed "
            f"{_n(organ.reabsorption.value)}"
        )
    out[organ.key] = rows
    factors = zone_factors(organ)
    for i in range(1, n + 1):
        zone = ids.zone_id(i)
        where = (
            "the whole organ"
            if n == 1
            else first
            if i == 1
            else last
            if i == n
            else f"between {first} and {last}"
        )
        zone_rows = [
            f"<b>zone {zone}</b> of {n}, {where}",
            f"vascular space <code>{ids.vascular(zone)}</code> "
            f"{_n(f_vas * volume / n)} ml, cells <code>{ids.cells(zone)}</code> "
            f"{_n(f_cell * volume / n)} ml",
        ]
        zone_rows += [
            f"<code>{ids.profile(p, zone)}</code> {_n(v)}"
            for p, v in factors[zone].items()
        ]
        if not factors[zone]:
            zone_rows.append("uniform: no zonal profile")
        out[f"zone:{zone}"] = zone_rows
    for q in qs:
        out[f"sub:{q}"] = [
            f"<b>{escape(q)}</b> {escape(quantity_label(q))}",
            escape(_fresh(d, q)),
            f"at every site <code>{escape(ids.species(q, '&lt;site&gt;'))}</code>",
        ]
    if isinstance(organ, LiverOrgan):
        out["bile"] = [
            "<b>bile</b>",
            "secreted by the hepatocytes of every zone",
            f"observers <code>{ids.BILE}</code> (µl) and <code>{ids.BILE_RATE}</code> "
            "(µl/min)",
        ]
    if kidney:
        out["urine"] = [
            f"<b>urine</b> <code>{ids.volume(ids.URINE_SITE)}</code>",
            "filtered: " + ", ".join(escape(q) for q in _urine(d)),
            f"observers <code>{ids.URINE}</code> (µl) and <code>{ids.URINE_RATE}</code> "
            "(µl/min)",
        ]
    reference = table_value("temperature", "reference").value
    out["temperature"] = [
        f"<b>temperature</b> <code>{ids.TEMPERATURE}</code> {_n(d.temperature.value)} °C",
        f"the rates of the core scale with the Arrhenius factor <code>{ids.F_T}</code>, "
        f"1 at {_n(reference)} °C",
        "the solubility of oxygen rises in the cold",
    ]
    return out


def interactive_circuit(
    definition: PerfusionDefinition, path: Path, *, figure: str = "anatomy"
) -> Path:
    """Write the interactive page of a figure of a definition.

    A figure with its hover cards, one self-contained HTML file
    (`pkpdsim.pbpk.visual.interactive.html_page`), which posts its height to
    the docs page which embeds it; the keyword and its default are those of
    `pkpdsim.pbpk.visual.interactive.interactive_figure`.

    Args:
        definition: the definition of the circuit.
        path: the HTML file to write.
        figure: `anatomy` (panel A, the organ in its circuit,
            `pkpdsim.perfusion.anatomy.circuit_anatomy_svg`, the default) or
            `schematic` (panel B, `circuit_svg`).

    Returns:
        The path of the page.

    Raises:
        ValueError: for another figure or an invalid definition
            (`PerfusionDefinition.validate`).
    """
    from pkpdsim.pbpk.visual.interactive import CREDIT, ICON_CREDIT, html_page
    from pkpdsim.perfusion.anatomy import anatomy_cards, circuit_anatomy_svg

    if figure == "anatomy":
        svg = circuit_anatomy_svg(definition)
        cards = anatomy_cards(definition)
        legend = (
            "Hover the organ, its zones, the apparatus, the labels and the chips "
            "for details."
        )
        credit = CREDIT
    elif figure == "schematic":
        svg = circuit_svg(definition)
        cards = circuit_cards(definition)
        legend = (
            "Hover the circuit, the organ, its zones, the chips and the process "
            "lines for details."
        )
        credit = ICON_CREDIT
    else:
        raise ValueError(
            f"No interactive figure '{figure}'; the figures are 'anatomy' and "
            "'schematic'."
        )
    width = int(float(svg.split('width="', 1)[1].split('"', 1)[0]))
    page = html_page(svg, cards, f"{definition.name} ({figure})", width, legend, credit)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(page, encoding="utf-8")
    return path
