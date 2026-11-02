"""The figure of a model: an SVG schematic of a definition.

`model_figure(definition, path)` draws the body of a definition as the
COMBINE archive of a model carries it (`create_archive(..., figure=True)`):

- the title, the substances (a chip in the color of the substance, its name
  and its mechanisms: binding, erythrocytes, distribution, baseline (an
  endogenous substance: with a given baseline, or produced with the baseline
  solved by the library), inhibitions and inductions), the routes, the PD components and the key;
- the circulation: the venous blood on the left, the arterial blood on the
  right, the lung in series between them, the organs which the arterial
  blood perfuses and which drain into the venous blood, the liver with the
  hepatic artery and the portal vein, and the organs which drain into the
  portal vein; the variant of an organ is named under it;
- in every organ a chip per substance, filled where the substance enters
  the tissue of the organ, dashed where it stays in the plasma, and the
  processes of the organ as badges: metabolism (`MET`, substrate, products
  and enzyme), transport (`TRA`, the side and the transporter), custom
  reaction (`CUS`, its equation), the excretion into the urine (`URI`, with
  the reabsorption), the slow pool (`SLW`, the exchange with its source, the
  tissue or the plasma);
- the lumen of the intestine: its segments in the order of the transit to
  the feces, the oral routes (`ORL`), the absorption (`ABS`), the bile
  (`BIL`), the conversions in the lumen (`CNV`) and the feces (`FEC`).

The figure is pure python and deterministic: the same definition gives the
same bytes. Its text is monospace and every label is placed by its width at
`CHAR_WIDTH` per character, an upper bound of the advance of the monospace
fonts it names (DejaVu Sans Mono 0.602 em), so a label stays inside its box
and labels do not overlap, also with many substances: a row grows by lines
instead of running out of its box.
"""

import math
import textwrap
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from xml.sax.saxutils import escape

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.chemistry import (
    MM,
    CustomReaction,
    Direction,
    FirstOrder,
    Hill,
    Membrane,
    Metabolism,
    Transport,
)
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.organ import OrganComponent
from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.pbpk.parts import parts_for
from pkpdsim.pbpk.routes import (
    DoseDependentDissolution,
    FirstOrderDissolution,
    Route,
    WeibullDissolution,
)
from pkpdsim.pbpk.stationary import steady_plan
from pkpdsim.pbpk.substance import Kp, KpMethod, Mode, Substance
from pkpdsim.pbpk.values import Formula, Q, Value, number
from pkpdsim.pbpk.variants import DEFAULTS, lookup
from pkpdsim.pbpk.variants.gut import bile_substances, luminal_substances, wall_parts
from pkpdsim.pbpk.variants.intestine import SegmentedIntestine, SingleLumen
from pkpdsim.pbpk.variants.kidney_nephron import NephronKidney
from pkpdsim.pbpk.variants.liver_zonated import ZonatedLiver
from pkpdsim.pbpk.visual import icons
from pkpdsim.sex import Sex

__all__ = [
    "CHAR_WIDTH",
    "COMPACT",
    "FONT",
    "LEGEND_WIDTH",
    "MAX_TEXT",
    "model_figure",
    "schematic_svg",
]

CHAR_WIDTH = 0.62
"""The width of a character in units of the font size: an upper bound of the
advance of the monospace fonts of `FONT` (DejaVu Sans Mono and Liberation
Mono 0.6 em, Menlo 0.602 em, Consolas 0.55 em)."""

FONT = "'DejaVu Sans Mono', Menlo, Consolas, 'Liberation Mono', monospace"
"""The fonts of the figure, monospace, so that `CHAR_WIDTH` bounds every label."""

WIDTH = 1200
"""The width of the figure (px)."""

_PAD = 16
_SIZE = 12
_SMALL = 11
_TITLE = 16
_LINE = 20
_VPAD = 5
_GAP = 8
_COLUMN = 76
_ARROW = 40
_BAR = 72
_ICON = 24

MAX_TEXT = 420
"""The longest line of the processes of a row (px); longer texts wrap, so the
figure is as wide as its content needs (`model_figure` trims the width)."""

LEGEND_WIDTH = 430
"""The width of the legend column right of the diagram (px)."""


@dataclass(frozen=True)
class _Dims:
    """The dimensions of a layout of the figure (px).

    Attributes:
        width: the width of the diagram (the circulation, without a legend
            column).
        line: the height of a line of text.
        vpad: the padding above and below the content of a box.
        gap: the gap between two rows.
        arrow: the length of the arrows between the blood and the rows.
    """

    width: float = WIDTH
    line: float = _LINE
    vpad: float = _VPAD
    gap: float = _GAP
    arrow: float = _ARROW

    @property
    def mx0(self) -> float:
        """The left edge of the rows of the circulation."""
        return _PAD + _COLUMN + self.arrow

    @property
    def mx1(self) -> float:
        """The right edge of the rows."""
        return self.width - _PAD - _COLUMN - self.arrow

    @property
    def px0(self) -> float:
        """The left edge of the rows of the portal organs."""
        return self.mx0 + _BAR + self.arrow

    def baseline(self, top: float, size: float) -> float:
        """The baseline of a text centred in a line which starts at `top`."""
        return top + self.line / 2 + 0.35 * size


COMPACT = _Dims(line=18, vpad=3, gap=4, arrow=28)
"""The dimensions of the compact figure, the panel of a publication figure."""

_COLORS = (
    "#4269d0",
    "#efb118",
    "#ff725c",
    "#6cc5b0",
    "#3ca951",
    "#ff8ab7",
    "#a463f2",
    "#97bbf5",
    "#9c6b4e",
    "#9498a0",
)
"""The colors of the substances, in the order of the definition (Observable 10)."""

_BADGES = {
    "MET": ("metabolism", "#e8f0fe", "#1a56db"),
    "TRA": ("transport", "#e6f6ec", "#1a7f37"),
    "CUS": ("custom reaction", "#fdf0e6", "#bc4c00"),
    "URI": ("urine", "#fff8c5", "#9a6700"),
    "SLW": ("slow pool", "#eef1f4", "#57606a"),
    "ORL": ("oral dose", "#fbefff", "#8250df"),
    "ABS": ("absorption", "#e6f6ec", "#1a7f37"),
    "BIL": ("bile", "#fff1e5", "#953800"),
    "CNV": ("lumen conversion", "#ffebe9", "#cf222e"),
    "FEC": ("feces", "#f6f1ea", "#6e5634"),
}
"""The badges of the processes: the label of the key, fill and stroke."""

_SYSTEMIC = (
    Organ.HEART,
    Organ.BRAIN,
    Organ.MUSCLE,
    Organ.SKIN,
    Organ.ADIPOSE,
    Organ.BONE,
    Organ.GONADS,
    Organ.KIDNEY,
)
"""The organs between the arterial and the venous blood, in the figure's order."""

_PORTAL = (
    Organ.SPLEEN,
    Organ.PANCREAS,
    Organ.STOMACH,
    Organ.SMALL_INTESTINE,
    Organ.LARGE_INTESTINE,
)
"""The organs which drain into the portal vein; the intestine last, above its lumen."""

_STYLE = f"""
.pkf-t{{font-family:{FONT};font-size:{_SIZE}px;fill:#1f2328}}
.pkf-s{{font-size:{_SMALL}px}}
.pkf-h{{font-size:{_TITLE}px;font-weight:700}}
.pkf-b{{font-weight:700}}
.pkf-dim{{fill:#57606a}}
.pkf-org{{fill:#f6f8fa;stroke:#8c959f;stroke-width:1.2}}
.pkf-ven{{fill:#dbe7fb;stroke:#4269d0;stroke-width:1.2}}
.pkf-art{{fill:#fbe0e0;stroke:#d0413e;stroke-width:1.2}}
.pkf-por{{fill:#efe3fa;stroke:#8a4fc0;stroke-width:1.2}}
.pkf-lum{{fill:#fff7e0;stroke:#b7891d;stroke-width:1.2}}
.pkf-seg{{fill:#ffffff;stroke:#b7891d;stroke-width:1}}
.pkf-line{{stroke:#57606a;stroke-width:1.4;fill:none}}
"""


def _width(text: str, size: float) -> float:
    """The width of a text in the figure (px)."""
    return len(text) * size * CHAR_WIDTH


def _chars(width: float, size: float) -> int:
    """The number of characters which fit into a width, at least 8."""
    return max(8, int(width // (size * CHAR_WIDTH)))


def _wrap(text: str, chars: int) -> list[str]:
    """The lines of a text of at most `chars` characters."""
    return textwrap.wrap(
        text, width=chars, break_long_words=True, break_on_hyphens=False
    ) or [""]


def _f(x: float) -> str:
    """A coordinate."""
    return f"{x:.1f}".removesuffix(".0")


@dataclass
class _Canvas:
    """The elements of the figure, in the order they are drawn.

    Attributes:
        k: the dimensions of the layout.
        elements: the SVG elements drawn so far.
        badges: the badges of the processes by code (`_BADGES`), the label
            of the key, fill and stroke; another figure (the circuit of a
            machine perfusion) draws its own codes on the same canvas.
    """

    k: _Dims = field(default_factory=_Dims)
    elements: list[str] = field(default_factory=list)
    badges: Mapping[str, tuple[str, str, str]] = field(default_factory=lambda: _BADGES)

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        cls: str = "",
        style: str = "",
        rx: float = 5.0,
    ) -> None:
        """Draw a rectangle."""
        c = f' class="{cls}"' if cls else ""
        s = f' style="{style}"' if style else ""
        self.elements.append(
            f'<rect{c}{s} x="{_f(x)}" y="{_f(y)}" width="{_f(w)}" '
            f'height="{_f(h)}" rx="{_f(rx)}"/>'
        )

    def text(
        self,
        x: float,
        y: float,
        text: str,
        cls: str = "",
        size: float = _SIZE,
    ) -> None:
        """Draw a text at its baseline `y`, its font size in `data-size`."""
        classes = f"pkf-t {cls}".strip()
        self.elements.append(
            f'<text class="{classes}" x="{_f(x)}" y="{_f(y)}" '
            f'data-size="{_f(size)}">{escape(text)}</text>'
        )

    def arrow(self, x1: float, y1: float, x2: float, y2: float) -> None:
        """Draw an arrow from the first to the second point."""
        self.elements.append(
            f'<line class="pkf-line" x1="{_f(x1)}" y1="{_f(y1)}" x2="{_f(x2)}" '
            f'y2="{_f(y2)}" marker-end="url(#pkf-arrow)"/>'
        )

    def dot(self, x: float, y: float, color: str) -> None:
        """Draw a dot of a color, centred at the point."""
        self.elements.append(
            f'<circle cx="{_f(x)}" cy="{_f(y)}" r="5" style="fill:{color}"/>'
        )

    def bend(self, points: Sequence[tuple[float, float]]) -> None:
        """Draw an arrow along the points, ending at the last one."""
        coordinates = " ".join(f"{_f(x)},{_f(y)}" for x, y in points)
        self.elements.append(
            f'<polyline class="pkf-line" points="{coordinates}" '
            'marker-end="url(#pkf-arrow)"/>'
        )

    def open(self, cls: str, key: str = "") -> None:
        """Open a group of a box and its labels, `key` its hover card in an interactive page."""
        k = f' data-k="{escape(key, {chr(34): "&quot;"})}"' if key else ""
        self.elements.append(f'<g class="{cls}"{k}>')

    def close(self) -> None:
        """Close a group."""
        self.elements.append("</g>")

    def chip(
        self,
        x: float,
        top: float,
        label: str,
        color: str,
        *,
        tissue: bool,
        substance: bool = False,
    ) -> float:
        """Draw the chip of a substance (or of the key) in a line; get its width."""
        w = _tag_width(label)
        if substance:
            self.elements.append(f'<g data-k="sub:{escape(label)}">')
        if tissue:
            style = f"fill:{color};fill-opacity:.28;stroke:{color};stroke-width:1.2"
        else:
            style = f"fill:#ffffff;stroke:{color};stroke-width:1.2;stroke-dasharray:3 2"
        self.rect(x, top + 2, w, self.k.line - 4, style=style, rx=4)
        self.text(x + 5, self.k.baseline(top, _SMALL), label, "pkf-s", _SMALL)
        if substance:
            self.elements.append("</g>")
        return w

    def badge(self, x: float, top: float, code: str) -> float:
        """Draw the badge of a process in a line; get its width."""
        _, fill, stroke = self.badges[code]
        w = _tag_width(code)
        self.rect(
            x, top + 2, w, self.k.line - 4, style=f"fill:{fill};stroke:{stroke}", rx=4
        )
        self.text(x + 5, self.k.baseline(top, _SMALL), code, "pkf-s pkf-b", _SMALL)
        return w


@dataclass
class _Row:
    """An organ (or the lumen) of the figure: its title, chips and processes."""

    title: list[str]
    chips: list[tuple[str, str, bool]]
    lines: list[tuple[str, str]]
    chain: bool = False
    key: str = ""
    """The organ key of the row (its icon and hover card), `gu` or `lumen`."""
    sink: tuple[str, list[tuple[str, str, bool]]] | None = None
    """A compartment at the right end of the row: its name and chips (urine, feces)."""


def _flow(widths: Sequence[float], width: float, gap: float) -> list[list[int]]:
    """The items of a flow by line, each line at most `width` wide."""
    lines: list[list[int]] = [[]]
    used = 0.0
    for i, w in enumerate(widths):
        if lines[-1] and used + gap + w > width:
            lines.append([])
            used = 0.0
        used += (gap if lines[-1] else 0.0) + w
        lines[-1].append(i)
    return lines


def _tag_width(label: str) -> float:
    """The width of a chip or a badge."""
    return _width(label, _SMALL) + 10


class _Figure:
    """The layout of the figure of a definition."""

    def __init__(
        self,
        definition: PBPKDefinition,
        sexes: Sequence[Sex],
        k: _Dims | None = None,
        *,
        compact: bool = False,
    ) -> None:
        """Collect what the figure shows."""
        self.k = k or _Dims()
        self.compact = compact
        self.maxx = 0.0
        self.text_width = self.k.width
        self.lx = 0.0
        self.d = definition
        self.sex = Sex(definition.sex).value
        self.sexes = tuple(sexes)
        self.colors = {
            s.sid: _COLORS[i % len(_COLORS)]
            for i, s in enumerate(definition.substances)
        }
        self.components: dict[Organ, OrganComponent] = {}
        for organ in Organ:
            choice = definition.variants.organs.get(
                organ, DEFAULTS.get(organ, "standard")
            )
            self.components[organ] = (
                lookup(organ, choice) if isinstance(choice, str) else choice
            )
        gut = self.components[Organ.SMALL_INTESTINE]
        self.gut = gut if isinstance(gut, SegmentedIntestine | SingleLumen) else None
        order = [s.sid for s in definition.substances]
        if self.gut is not None:
            reach = luminal_substances(
                definition.substances,
                definition.chemistry,
                definition.routes,
                self.gut.conversions,
            )
            bile = bile_substances(definition.substances, definition.chemistry)
            self.luminal = [sid for sid in order if sid in reach]
            self.bile = [sid for sid in order if sid in bile]
            # an absorption fixed at 0 cannot happen (`ka={"cr": 0.0}`, a
            # substance which a transport takes up instead): not drawn
            self.absorbed = {
                sid: options
                for sid, options in self.gut.absorption(definition).items()
                if not _fixed_zero(_ka_of(self.gut, sid))
            }
        else:
            self.luminal, self.bile, self.absorbed = [], [], {}
        # a substance without baseline whose steady state the library solves
        # and which a declared rate produces (`pkpdsim.pbpk.stationary`)
        plan = steady_plan(definition)
        self.produced = {
            sid
            for (_, sid), steady in plan.organs.items()
            if sid in plan.solved and steady.has_rate
        }
        self.canvas = _Canvas(self.k)
        self.codes: set[str] = set()

    # --- content -----------------------------------------------------------

    def _variant(self, organ: Organ) -> str | None:
        """The name of the variant of an organ, `None` for the standard organ."""
        c = self.components[organ]
        if isinstance(c, ZonatedLiver):
            return f"zonated, {c.n_zones} zones"
        if isinstance(c, NephronKidney):
            return "nephron"
        if isinstance(c, SegmentedIntestine):
            return "segmented"
        if isinstance(c, SingleLumen):
            return "single lumen"
        name = getattr(c, "name", "standard")
        return None if name == "standard" else str(name)

    def _chips(self, organ: Organ) -> list[tuple[str, str, bool]]:
        """The chips of the substances in an organ: id, color, enters the tissue."""
        wall = self.gut is not None and organ in self.gut.replaces
        return [
            (
                s.sid,
                self.colors[s.sid],
                (
                    wall_parts(s, self.d.chemistry)
                    if wall
                    else parts_for(s, organ, self.d.chemistry)
                ).tissue
                != "none",
            )
            for s in self.d.substances
        ]

    @staticmethod
    def _kinetics(k: FirstOrder | MM | Hill) -> str:
        """The rate law in words."""
        if isinstance(k, MM):
            return "Michaelis-Menten"
        if isinstance(k, Hill):
            return "Hill"
        return "first order"

    def _apical(self, organ: Organ) -> str:
        """The side of the apical membrane of an organ."""
        if organ is Organ.LIVER:
            return "bile"
        if organ is Organ.KIDNEY:
            return "urine"
        if organ in (Organ.SMALL_INTESTINE, Organ.LARGE_INTESTINE):
            return "lumen"
        return "apical side"

    def _process(self, p: Metabolism | Transport | CustomReaction) -> tuple[str, str]:
        """The badge and the text of a process."""
        if isinstance(p, Metabolism):
            products = ", ".join(p.products) if p.products else "eliminated"
            enzyme = f" by {p.protein}" if p.protein else ""
            return (
                "MET",
                f"{p.sid}: {p.substrate} -> {products}{enzyme} "
                f"({self._kinetics(p.kinetics)})",
            )
        if isinstance(p, Transport):
            side = (
                "blood" if p.membrane is Membrane.BASOLATERAL else self._apical(p.organ)
            )
            way = (
                f"{side} -> cell"
                if p.direction is Direction.INFLUX
                else f"cell -> {side}"
            )
            carrier = f" by {p.protein}" if p.protein else ""
            return (
                "TRA",
                f"{p.sid}: {p.substance} {way}{carrier} ({self._kinetics(p.kinetics)})",
            )
        return "CUS", f"{p.sid}: {' '.join(p.equation.split())}"

    def _organ_lines(self, organs: Sequence[Organ]) -> list[tuple[str, str]]:
        """The badges and texts of the processes of organs."""
        lines = [self._process(p) for p in self.d.chemistry if p.organ in organs]
        if Organ.KIDNEY in organs:
            urine = []
            for s in self.d.substances:
                e = s.excretion
                if not e.urine:
                    continue
                if e.reabsorption is not None:
                    urine.append(f"{s.sid} (threshold reabsorption)")
                elif number(e.fraction_reabsorbed) != 0.0:
                    urine.append(f"{s.sid} (reabsorbed)")
                else:
                    urine.append(s.sid)
            if urine:
                lines.append(("URI", "filtered into the urine: " + ", ".join(urine)))
        for pool in self.d.slow_pools:
            if any(o in pool.organs for o in organs):
                lines.append(
                    (
                        "SLW",
                        f"slow pool of {pool.substance}: "
                        f"{Part(pool.source).name.lower()} <-> pool",
                    )
                )
        return lines

    def _title(self, organs: Sequence[Organ], name: str) -> list[str]:
        """The title of a row: the name and ids, the variant."""
        ids = " ".join(o.value for o in organs)
        title = [f"{name} {ids}"]
        variant = self._variant(organs[0])
        if variant is not None:
            title.append(variant)
        return title

    def _row(self, organs: Sequence[Organ], name: str) -> _Row:
        """The row of an organ (or of the organs the intestine replaces)."""
        return _Row(
            self._title(organs, name),
            self._chips(organs[0]),
            self._organ_lines(organs),
            key=organs[0].value if len(organs) == 1 else "gu",
        )

    def _lumen_row(self, gut: SegmentedIntestine | SingleLumen) -> _Row:
        """The row of the lumen: the segments, the routes and the fluxes.

        Which substance reaches the lumen, is excreted into the bile and is
        absorbed are the rules of the build (`luminal_substances`,
        `bile_substances`, the `absorption` of the intestine variant).
        """
        if isinstance(gut, SegmentedIntestine):
            segments = [*ids.GUT_SEGMENTS]
        else:
            segments = ["lumen"]
        chips = [(seg, "", True) for seg in segments]
        lines: list[tuple[str, str]] = []
        for r in self.d.routes:
            if r.kind == "oral":
                lines.append(("ORL", f"tablet {r.substance} ({_dissolution(r)})"))
            elif r.kind == "oral_solution":
                lines.append(("ORL", f"solution {r.substance}"))
        if self.absorbed:
            absorbed = [
                f"{sid} ({', '.join(options)})" if options else sid
                for sid, options in self.absorbed.items()
            ]
            lines.append(("ABS", "into the intestinal wall: " + "; ".join(absorbed)))
        if self.bile:
            entry = segments[1] if len(segments) > 1 else segments[0]
            lines.append(("BIL", f"bile into {entry}: " + ", ".join(self.bile)))
        for c in gut.conversions:
            lines.append(
                ("CNV", f"{c.substrate} -> {c.product} in {', '.join(c.segments)}")
            )
        lines.append(("FEC", "from the last segment: " + ", ".join(self.luminal)))
        title = f"{len(segments)} segments" if len(segments) > 1 else "one segment"
        sink = ("feces", [(sid, self.colors[sid], True) for sid in self.luminal])
        return _Row(["lumen", title], chips, lines, chain=True, key="lumen", sink=sink)

    # --- legend ------------------------------------------------------------

    def _details(self, s: Substance) -> str:
        """The mechanisms of a substance in words."""
        parts = [s.name]
        if s.binding is not None:
            protein = (
                "" if s.binding.protein == "unknown" else f" ({s.binding.protein})"
            )
            parts.append(f"bound in plasma{protein}")
        if s.erythrocytes is not None:
            parts.append(f"erythrocytes ({s.erythrocytes.mode})")
        dist = s.distribution
        if dist is not None:
            if isinstance(dist, KpMethod):
                kind = f"Kp of {dist.method}"
            elif dist.interstitial is not None:
                kind = "extracellular (interstitial share)"
            else:
                kind = "Kp"
            words = [kind, _mode(dist.mode)]
            if dist.organs is not None:
                words.append(f"in {_organs(dist.organs)} only")
            if isinstance(dist, Kp) and dist.overrides:
                words.append(f"own Kp in {_organs(dist.overrides)}")
            if dist.ps:
                words.append(f"own exchange in {_organs(dist.ps)}")
            parts.append(", ".join(words))
        if s.baseline is not None:
            parts.append("endogenous (baseline)")
        elif s.sid in self.produced:
            parts.append("endogenous (production, baseline solved)")
        for i in s.inhibits:
            parts.append(f"inhibits {i.protein} ({i.kind})")
        for i in s.induces:
            parts.append(f"induces {i.protein}")
        return "; ".join(parts)

    def _route_text(self) -> str:
        """The routes in words."""
        words = []
        for r in self.d.routes:
            if r.kind == "iv":
                words.append(f"iv {r.substance}")
            elif r.kind == "oral":
                words.append(f"oral tablet {r.substance}")
            elif r.kind == "oral_solution":
                words.append(f"oral solution {r.substance}")
            else:
                words.append(f"{r.kind} {r.substance}")
        return "routes: " + ("; ".join(words) if words else "none")

    # --- drawing -----------------------------------------------------------

    def _paragraph(self, x: float, y: float, text: str, cls: str = "") -> float:
        """Draw a wrapped text from `x` to the right margin; get the next top."""
        for line in _wrap(text, _chars(self.text_width - _PAD - (x - self.lx), _SIZE)):
            self.canvas.text(x, self.k.baseline(y, _SIZE), line, cls)
            y += self.k.line
        return y

    def _legend(self, y: float) -> float:
        """Draw the title, the substances, routes, PD and the key; get the next top."""
        c, d = self.canvas, self.d
        # a long name wraps within the legend
        width = self.text_width - 2 * _PAD
        for line in _wrap(d.name or d.sid, _chars(width, _TITLE)):
            c.text((self.lx + _PAD), y + _TITLE, line, "pkf-h", _TITLE)
            y += _TITLE + 4
        y += 6
        y = self._paragraph(
            (self.lx + _PAD),
            y,
            f"PBPK model {d.sid}: {_builds(self.sexes)}, population {d.population}, "
            f"arm vein {d.arm_vein}",
            "pkf-dim",
        )
        y += 6
        c.text((self.lx + _PAD), self.k.baseline(y, _SIZE), "substances", "pkf-b")
        y += self.k.line
        # a dot of the color of the substance, not a chip of the key
        sid_w = max((_width(s.sid, _SIZE) for s in d.substances), default=0.0)
        x_text = (self.lx + _PAD) + 16 + sid_w + 12
        for s in d.substances:
            c.dot((self.lx + _PAD) + 5, y + self.k.line / 2, self.colors[s.sid])
            c.text((self.lx + _PAD) + 16, self.k.baseline(y, _SIZE), s.sid, "pkf-b")
            y = self._paragraph(x_text, y, self._details(s))
        y = self._paragraph((self.lx + _PAD), y + 4, self._route_text())
        for pd in d.pd:
            reads = ", ".join(pd.inputs.values()) or "nothing"
            sets = ", ".join(pd.outputs.values()) or "nothing"
            y = self._paragraph(
                (self.lx + _PAD), y, f"PD {pd.model.sid}: reads {reads}; sets {sets}"
            )
        return y + 6

    def _key(self, y: float) -> float:
        """Draw the key of the chips and the badges used; get the next top."""
        c = self.canvas
        items: list[tuple[str, str]] = [
            ("tissue", "enters the tissue"),
            ("plasma", "plasma only"),
        ]
        items += [(code, _BADGES[code][0]) for code in _BADGES if code in self.codes]
        widths = []
        for kind, label in items:
            head = _tag_width("s") if kind in ("tissue", "plasma") else _tag_width(kind)
            widths.append(head + 6 + _width(label, _SIZE))
        lines = _flow(widths, self.text_width - 2 * _PAD - _width("key ", _SIZE), 18)
        c.text((self.lx + _PAD), self.k.baseline(y, _SIZE), "key", "pkf-b")
        x0 = (self.lx + _PAD) + _width("key ", _SIZE)
        for line in lines:
            x = x0
            for i in line:
                kind, label = items[i]
                if kind in ("tissue", "plasma"):
                    w = c.chip(x, y, "s", "#57606a", tissue=kind == "tissue")
                else:
                    w = c.badge(x, y, kind)
                c.text(x + w + 6, self.k.baseline(y, _SIZE), label)
                x += widths[i] + 18
            y += self.k.line
        return y + 10

    def _draw_row(
        self,
        row: _Row,
        x0: float,
        y: float,
        x1: float,
        cx: float,
        cls: str = "pkf-org",
    ) -> float:
        """Draw a row from `x0` to `x1` with its chips and badges from `cx`; get its height.

        The box of the row is drawn first, behind its labels, with the height
        of its content: the icon of the organ, the title, the chips, the
        processes (each line a group with its full text, the hover card of an
        interactive page) and a compartment at the right end (urine, feces).
        """
        c = self.canvas
        top = y + self.k.vpad
        c.open("pkf-box", row.key)
        box = len(c.elements)
        icon = (
            icons.organ(
                row.key, x0 + 7 + _ICON / 2, top + _ICON / 2 - 1, _ICON, self.sex
            )
            if row.key
            else ""
        )
        ix = _ICON + 6 if icon else 0
        if icon:
            c.elements.append(icon)
        title_chars = _chars(cx - x0 - 20 - ix, _SIZE)
        t = top
        for i, line in enumerate(row.title):
            for part in _wrap(line, title_chars):
                c.text(
                    x0 + 10 + ix,
                    self.k.baseline(t, _SIZE),
                    part,
                    "pkf-b" if i == 0 else "pkf-dim",
                )
                t += self.k.line
        right = min(x1 - 10, cx + MAX_TEXT)
        sink: tuple[float, float] | None = None
        if row.sink is not None:
            name, schips = row.sink
            sw = (
                max([_width(name, _SIZE)] + [_tag_width(lab) for lab, _, _ in schips])
                + 16
            )
            sink = (x1 - 10 - sw, sw)
            right = min(right, x1 - 10 - sw - 26)
        gap = 18.0 if row.chain else 6.0
        widths = [_tag_width(label) for label, _, _ in row.chips]
        lines = _flow(widths, right - cx, gap)
        rowmax = cx
        cy = top
        for line in lines:
            x = cx
            for k, i in enumerate(line):
                label, color, tissue = row.chips[i]
                if row.chain:
                    c.rect(
                        x,
                        cy + 2,
                        widths[i],
                        self.k.line - 4,
                        "pkf-seg" if tissue else "pkf-org",
                        rx=4,
                    )
                    c.text(x + 5, self.k.baseline(cy, _SMALL), label, "pkf-s", _SMALL)
                    if k + 1 < len(line):
                        c.arrow(
                            x + widths[i] + 2,
                            cy + self.k.line / 2,
                            x + widths[i] + gap - 1,
                            cy + self.k.line / 2,
                        )
                else:
                    c.chip(x, cy, label, color, tissue=tissue, substance=True)
                rowmax = max(rowmax, x + widths[i])
                x += widths[i] + gap
            cy += self.k.line
        for code, text in row.lines:
            bw = _tag_width(code)
            tx = cx + bw + 6
            c.open("pkf-proc", f"txt:{code}: {text}")
            c.badge(cx, cy, code)
            for part in _wrap(text, _chars(right - tx, _SIZE)):
                c.text(tx, self.k.baseline(cy, _SIZE), part)
                rowmax = max(rowmax, tx + _width(part, _SIZE))
                cy += self.k.line
            c.close()
        if sink is not None and row.sink is not None:
            sx, sw = sink
            name, schips = row.sink
            sh = (1 + len(schips)) * self.k.line + 6
            style = (
                "fill:#fff8d6;stroke:#b08a00;stroke-width:1.4"
                if name == "urine"
                else "fill:#efe2d3;stroke:#7a5a3a;stroke-width:1.4"
            )
            c.open("pkf-sink", name)
            c.rect(sx, top, sw, sh, style=style, rx=6)
            c.text(sx + 8, self.k.baseline(top + 2, _SIZE), name, "pkf-b")
            for j, (lab, col, tis) in enumerate(schips):
                c.chip(
                    sx + 8,
                    top + 2 + (j + 1) * self.k.line,
                    lab,
                    col,
                    tissue=tis,
                    substance=True,
                )
            c.close()
            c.arrow(
                sx - 22, top + self.k.line / 2 + 2, sx - 2, top + self.k.line / 2 + 2
            )
            cy = max(cy, top + sh)
            rowmax += sw + 26
        self.maxx = max(self.maxx, rowmax)
        height = max(t, cy) - y + self.k.vpad
        content = c.elements[box:]
        del c.elements[box:]
        c.rect(x0, y, x1 - x0, height, cls)
        c.elements += content
        c.close()
        return height

    def _band(self, x0: float, x1: float, top: float) -> float:
        """The legend of the compact figure in the band above the first row.

        The substances and the key of the chips and badges, flowing in lines;
        get the next top.
        """
        c = self.canvas
        items: list[tuple[str, str, str, str]] = [
            ("dot", s.sid, s.sid, self.colors[s.sid]) for s in self.d.substances
        ]
        items += [
            ("tissue", "s", "enters the tissue", ""),
            ("plasma", "s", "plasma only", ""),
        ]
        items += [
            (code, code, _BADGES[code][0], "") for code in _BADGES if code in self.codes
        ]

        def width(item: tuple[str, str, str, str]) -> float:
            kind, head, label, _ = item
            if kind == "dot":
                return 14 + _width(label, _SIZE)
            return _tag_width(head) + 6 + _width(label, _SIZE)

        y = top
        x = x0
        for item in items:
            w = width(item)
            if x > x0 and x + w > x1:
                x, y = self.k.mx0, y + self.k.line
            kind, head, label, color = item
            if kind == "dot":
                c.dot(x + 5, y + self.k.line / 2, color)
                c.text(x + 14, self.k.baseline(y, _SIZE), label, "pkf-b")
            elif kind in ("tissue", "plasma"):
                w0 = c.chip(x, y, "s", "#57606a", tissue=kind == "tissue")
                c.text(x + w0 + 6, self.k.baseline(y, _SIZE), label)
            else:
                w0 = c.badge(x, y, head)
                c.text(x + w0 + 6, self.k.baseline(y, _SIZE), label)
            x += w + 16
        return y + self.k.line

    def draw(self) -> str:
        """Draw the figure; get the SVG."""
        c = self.canvas
        rows = self._rows()
        self.codes = {code for row in rows.values() for code, _ in row.lines}
        legend: list[str] = []
        legend_bottom = 0.0
        if not self.compact:
            # the legend in a column right of the diagram
            diagram, self.canvas = self.canvas, _Canvas(self.k)
            self.text_width, self.lx = LEGEND_WIDTH, self.k.width - 4
            legend_bottom = self._key(self._legend(_PAD))
            legend, self.canvas = self.canvas.elements, diagram
            self.text_width, self.lx = self.k.width, 0.0
            c = self.canvas
        y = _PAD
        circulation = {
            Organ.LUNG.value,
            *(o.value for o in _SYSTEMIC),
            Organ.LIVER.value,
        }

        def chips_x(x0: float, keys: Iterable[str]) -> float:
            """The column of the chips of rows: after their longest title."""
            titles = max(_width(t, _SIZE) for k in keys for t in rows[k].title)
            return x0 + titles + 30 + _ICON + 6

        portal = [k for k in rows if k not in circulation and k != "lumen"]
        # one column of chips for all organ rows, after the longest title
        column = max(chips_x(self.k.mx0, circulation), chips_x(self.k.px0, portal))
        columns_x = dict.fromkeys([*circulation, *portal], column)
        ve0, ve1 = _PAD, _PAD + _COLUMN
        ar0 = self.k.width - _PAD - _COLUMN

        body_top = y
        rows_top = body_top + 2 * self.k.line + 14
        # the intravenous and the depot routes enter the venous blood from the
        # band above the first row, so the blood columns start at the top
        entries = [
            f"{kind}: " + ", ".join(sids)
            for kind in ("iv", "im", "sc", "rectal")
            if (
                sids := list(
                    dict.fromkeys(r.substance for r in self.d.routes if r.kind == kind)
                )
            )
        ]
        band_x = self.k.mx0
        if entries:
            text = "; ".join(entries)
            mid = body_top + 4 + self.k.line / 2
            c.elements.append(icons.syringe(self.k.mx0 + 13, mid, 22, rot=90))
            c.text(self.k.mx0 + 30, self.k.baseline(body_top + 4, _SIZE), text)
            c.arrow(self.k.mx0 - 2, mid, ve1 + 2, mid)
            band_x = self.k.mx0 + 30 + _width(text, _SIZE) + 24
            self.maxx = max(self.maxx, band_x)
        if self.compact:
            band_bottom = self._band(band_x, ar0 - 10, body_top + 4)
            rows_top = max(rows_top, band_bottom + 8)
        columns = len(c.elements)

        def place(key: str, x0: float, top: float, cls: str = "pkf-org") -> float:
            return top + self._draw_row(
                rows[key], x0, top, self.k.mx1, columns_x[key], cls
            )

        # the lung between the venous and the arterial blood
        bottom = place(Organ.LUNG.value, self.k.mx0, rows_top)
        mid = (rows_top + bottom) / 2
        c.arrow(ve1, mid, self.k.mx0, mid)
        c.arrow(self.k.mx1, mid, ar0, mid)
        top = bottom + self.k.gap
        for organ in [*_SYSTEMIC, Organ.LIVER]:
            bottom = place(organ.value, self.k.mx0, top)
            mid = (top + bottom) / 2
            c.arrow(ar0, mid, self.k.mx1, mid)
            c.arrow(self.k.mx0, mid, ve1, mid)
            top = bottom + self.k.gap
        venous_bottom = bottom
        # the portal vein from the organs of the portal circulation into the liver
        bar_top = venous_bottom + 24
        top = bar_top
        last = top
        for key in portal:
            bottom = place(key, self.k.px0, top)
            mid = (top + bottom) / 2
            c.arrow(ar0, mid, self.k.mx1, mid)
            c.arrow(self.k.px0, mid, self.k.mx0 + _BAR, mid)
            last = bottom
            top = bottom + self.k.gap
        # the portal vein lists the circulating substances: as long as its chips need
        chips_h = (2 + len(self.d.substances)) * self.k.line + 14
        arterial_bottom = max(last, bar_top + chips_h, body_top + chips_h)
        c.arrow(self.k.mx0 + _BAR / 2, bar_top, self.k.mx0 + _BAR / 2, venous_bottom)
        # the blood columns and the portal vein, behind the arrows
        blood = _Canvas(self.k)
        blood.open("pkf-box")
        blood.rect(ve0, body_top, _COLUMN, venous_bottom - body_top, "pkf-ven")
        blood.text(ve0 + 8, self.k.baseline(body_top + 4, _SIZE), "venous", "pkf-b")
        blood.text(
            ve0 + 8,
            self.k.baseline(body_top + 4 + self.k.line, _SIZE),
            "blood",
            "pkf-b",
        )
        self._circulating(blood, ve0, body_top)
        blood.close()
        blood.open("pkf-box")
        blood.rect(ar0, body_top, _COLUMN, arterial_bottom - body_top, "pkf-art")
        blood.text(ar0 + 8, self.k.baseline(body_top + 4, _SIZE), "arterial", "pkf-b")
        blood.text(
            ar0 + 8,
            self.k.baseline(body_top + 4 + self.k.line, _SIZE),
            "blood",
            "pkf-b",
        )
        self._circulating(blood, ar0, body_top)
        blood.close()
        blood.open("pkf-box")
        blood.rect(self.k.mx0, bar_top, _BAR, arterial_bottom - bar_top, "pkf-por")
        blood.text(
            self.k.mx0 + 8, self.k.baseline(bar_top + 4, _SIZE), "portal", "pkf-b"
        )
        blood.text(
            self.k.mx0 + 8,
            self.k.baseline(bar_top + 4 + self.k.line, _SIZE),
            "vein",
            "pkf-b",
        )
        self._circulating(blood, self.k.mx0, bar_top)
        blood.close()
        c.elements[columns:columns] = blood.elements
        y = arterial_bottom
        # the lumen of the intestine below its wall, absorbing into it
        if "lumen" in rows:
            lumen_top = y + 30
            lumen_cx = (
                self.k.mx0
                + max(_width(t, _SIZE) for t in rows["lumen"].title)
                + 36
                + _ICON
            )
            bottom = lumen_top + self._draw_row(
                rows["lumen"], self.k.mx0, lumen_top, self.k.mx1, lumen_cx, "pkf-lum"
            )
            if self.absorbed:
                # the absorption of `ka_<s>`: the oral substances and those a
                # model reabsorbs after their biliary excretion
                x = self.k.px0 + 24
                c.arrow(x, lumen_top, x, y)
                c.text(
                    x + 10,
                    self.k.baseline(y + 5, _SMALL),
                    "absorption",
                    "pkf-s pkf-dim",
                    _SMALL,
                )
            if any(r.kind != "iv" for r in self.d.routes):
                c.text(_PAD, self.k.baseline(lumen_top + self.k.vpad, _SIZE), "oral")
                c.elements.append(
                    icons.tablet(
                        _PAD + 13, lumen_top + self.k.vpad + self.k.line + 12, 22
                    )
                )
                mid = lumen_top + self.k.vpad + self.k.line / 2
                c.arrow(_PAD + _width("oral", _SIZE) + 6, mid, self.k.mx0, mid)
            if self.bile:
                # the bile duct from the liver into the lumen, left of the
                # portal vein, crossing no other arrow
                channel = (ve1 + self.k.mx0) / 2
                below = y + 12
                c.bend(
                    [
                        (self.k.mx0 + 12, venous_bottom),
                        (self.k.mx0 + 12, venous_bottom + 10),
                        (channel, venous_bottom + 10),
                        (channel, below),
                        (self.k.mx0 + 40, below),
                        (self.k.mx0 + 40, lumen_top),
                    ]
                )
                c.text(
                    self.k.mx0 + 48,
                    self.k.baseline(y + 5, _SMALL),
                    "bile",
                    "pkf-s pkf-dim",
                    _SMALL,
                )
            y = bottom
        height = y + _PAD
        width = self.k.width
        if legend:
            height = max(height, legend_bottom + _PAD)
            c.elements += legend
            width = self.k.width - 4 + LEGEND_WIDTH
        return "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{_f(width)}" '
                f'height="{_f(height)}" viewBox="0 0 {_f(width)} {_f(height)}" '
                'role="img" aria-labelledby="pkf-title">',
                f'<title id="pkf-title">{escape(self.d.name or self.d.sid)}</title>',
                f"<style>{_STYLE}</style>",
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

    def _circulating(self, canvas: _Canvas, x: float, top: float) -> None:
        """The chips of the substances which circulate in a blood column, below its title."""
        for j, s in enumerate(self.d.substances):
            canvas.chip(
                x + 8,
                top + 10 + (2 + j) * self.k.line,
                s.sid,
                self.colors[s.sid],
                tissue=True,
                substance=True,
            )

    def _rows(self) -> dict[str, _Row]:
        """The rows of the figure by key: the organ id, `gu` or `lumen`."""
        names = {
            Organ.LUNG: "lung",
            Organ.HEART: "heart",
            Organ.BRAIN: "brain",
            Organ.MUSCLE: "muscle",
            Organ.SKIN: "skin",
            Organ.ADIPOSE: "adipose",
            Organ.BONE: "bone",
            Organ.GONADS: "gonads",
            Organ.KIDNEY: "kidney",
            Organ.LIVER: "liver",
            Organ.SPLEEN: "spleen",
            Organ.PANCREAS: "pancreas",
            Organ.STOMACH: "stomach",
            Organ.SMALL_INTESTINE: "small intestine",
            Organ.LARGE_INTESTINE: "large intestine",
        }
        rows = {
            o.value: self._row([o], names[o])
            for o in (Organ.LUNG, *_SYSTEMIC, Organ.LIVER)
        }
        for o in _PORTAL:
            if self.gut is not None and o in self.gut.replaces:
                if o is Organ.SMALL_INTESTINE:
                    rows["gu"] = self._row(list(self.gut.replaces), "intestine")
                continue
            rows[o.value] = self._row([o], names[o])
        urine = [s.sid for s in self.d.substances if s.excretion.urine]
        if urine:
            rows[Organ.KIDNEY.value].sink = (
                "urine",
                [(sid, self.colors[sid], True) for sid in urine],
            )
        if self.gut is not None and self.luminal:
            # a lumen which no substance reaches is not drawn
            rows["lumen"] = self._lumen_row(self.gut)
        return rows


def _mode(mode: Mode) -> str:
    """The exchange of a distribution in words."""
    return {
        Mode.TISSUE: "tissue exchange",
        Mode.PERFUSION: "perfusion-limited",
        Mode.PERMEABILITY: "permeability-limited",
    }.get(mode, str(mode))


def _ka_of(gut: SegmentedIntestine | SingleLumen, sid: str) -> Q | None:
    """The `ka` which the definition gives a substance, `None` for the default."""
    if isinstance(gut.ka, Mapping):
        return gut.ka.get(sid)
    return gut.ka


def _fixed_zero(q: Q | None) -> bool:
    """Whether a value is fixed at 0: a number or an unfitted `Value` of 0."""
    if q is None or isinstance(q, Formula):
        return False
    if isinstance(q, Value) and q.uncertainties:
        return False
    return number(q) == 0.0


def _organs(organs: Iterable[Organ]) -> str:
    """The ids of organs in the order of `Organ`."""
    chosen = set(organs)
    return ", ".join(o.value for o in Organ if o in chosen)


def _dissolution(route: Route) -> str:
    """The dissolution of a tablet in words."""
    f = route.formulation
    if isinstance(f, WeibullDissolution):
        return "Weibull dissolution"
    if isinstance(f, DoseDependentDissolution):
        return "dose-dependent dissolution"
    if isinstance(f, FirstOrderDissolution):
        factors = [
            name
            for name, q in (("fed factor", f.fed), ("release factor", f.release))
            if q is not None
        ]
        return ", ".join(["first order dissolution", *factors])
    return "first order dissolution"


def _builds(sexes: Sequence[Sex]) -> str:
    """The builds of the figure in words: `male build`, `male and female builds`."""
    if len(sexes) == 1:
        return f"{sexes[0]} build"
    return f"{', '.join(str(s) for s in sexes[:-1])} and {sexes[-1]} builds"


def model_figure(
    definition: PBPKDefinition,
    path: Path,
    sexes: Sequence[Sex] | None = None,
    *,
    compact: bool = False,
) -> Path:
    """Write the figure of a definition, an SVG schematic of the model.

    The circulation with the organs, the substances in each organ, the
    routes and the processes (see the module), deterministic and without a
    dependency beyond the library. The builds of the sexes have the same
    organs, substances and processes, so one figure shows them all; its
    title names the builds (`male build`, `male and female builds`).

    Args:
        definition: the definition of the model.
        path: the SVG file to write.
        sexes: the sexes which are built (`create_archive(..., sexes=...)`),
            the sex of the definition among them; `None` is the sex of the
            definition only.
        compact: the panel of a publication figure: no legend column (the
            substance dots and the key in the band above the first row, the
            names and routes in the caption of the figure), tighter rows and
            shorter arrows (`COMPACT`).

    Returns:
        The path of the figure.

    Raises:
        ValueError: if `sexes` leaves out the sex of the definition or names
            a sex twice.
    """
    sexes = (definition.sex,) if sexes is None else tuple(Sex(s) for s in sexes)
    if definition.sex not in sexes or len(set(sexes)) != len(sexes):
        raise ValueError(
            f"The figure of '{definition.sid}' names the builds "
            f"{[str(s) for s in sexes]}: each sex once, the sex of the definition "
            f"'{definition.sex}' among them."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(schematic_svg(definition, sexes, compact=compact), encoding="utf-8")
    return path


def schematic_svg(
    definition: PBPKDefinition, sexes: Sequence[Sex], *, compact: bool = False
) -> str:
    """The SVG of the figure of a definition (see `model_figure`).

    The figure is drawn twice: once wide to measure its content, then at the
    width the content needs.

    Args:
        definition: the definition of the model.
        sexes: the sexes which are built.
        compact: the panel of a publication figure.

    Returns:
        The SVG, deterministic.
    """
    dims = COMPACT if compact else _Dims()
    probe = _Figure(definition, sexes, replace(dims, width=2400), compact=compact)
    probe.draw()
    width = math.ceil(probe.maxx + 16 + _PAD + _COLUMN + dims.arrow)
    return _Figure(
        definition, sexes, replace(dims, width=width), compact=compact
    ).draw()
