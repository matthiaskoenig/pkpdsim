"""The anatomical figure of a model built on a laboratory animal (`PBPKDefinition.species`).

The drawing is the anatomogram of the male mouse of the EBI Expression Atlas, for
the mouse and for the rat (the rat drawing of the atlas is a side view which
lacks the stomach, the small intestine, the pancreas, the muscle, the bone and
the adipose tissue), and for the pig a silhouette of pkpdsim with the organs of
the human drawing of its sex in the trunk: the reference animal of the library
as drawn, without deformation and without layers, of the sex of the build (the
library has the male mouse and rat, the male and female pig); the figure fills
its height, so the species are not on one scale. The labels, submodels, icons
and key are those of the human figure (`pkpdsim.pbpk.visual.anatomy`); the oral
gavage enters at the snout, the intravenous dose and the blood sample are at
the tail vein of a rodent and the ear vein of the pig (`VEINS`).
"""

from collections.abc import Mapping, Sequence
from functools import cache
from types import MappingProxyType
from xml.sax.saxutils import escape

import numpy as np

from pkpdsim.pbpk.visual import icons, palette
from pkpdsim.pbpk.visual.anatomy import FONT, AnatomyFigure, _curve
from pkpdsim.pbpk.visual.body import path
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.resources import Anatomy, anatomy

__all__ = [
    "DRAWINGS",
    "LEFT",
    "RIGHT",
    "SKINS",
    "VEINS",
    "AnimalBody",
    "AnimalFigure",
    "figure_of",
    "reference_volumes",
]

DRAWINGS: Mapping[str, str] = MappingProxyType(
    {"mouse": "mouse", "rat": "mouse", "pig": "pig"}
)
"""The drawing of every laboratory animal: the resource of its species (the mouse for the rat)."""

VEINS: Mapping[str, str] = MappingProxyType(
    {"mouse": "tail vein", "rat": "tail vein", "pig": "ear vein"}
)
"""The vein of the intravenous dose and the blood sample of every animal."""

SKINS: Mapping[str, str] = MappingProxyType({"pig": "#f4c2c2"})
"""The skin of an animal, `palette.SKIN` by default."""

LEFT = ("br", "lu", "he", "li", "ki", "la", "sk", "mu", "blood")
"""The labels of the left column."""

RIGHT = ("st", "sp", "pa", "si", "iv", "lumen", "ad", "bo", "go", "urine")
"""The labels of the right column."""

_ORDER = (
    "mu", "ad", "lu", "la", "si", "ki", "sp", "pa", "st", "li",
    "gallbladder", "he", "bladder", "uterus", "go", "bo", "br",
)  # fmt: skip


@cache
def reference_volumes(species: str, sex: str = "male") -> Mapping[str, float]:
    """The organ masses of the reference animal of the library [kg], by organ key."""
    from pkpdsim.database.tables import load

    organs = load(f"physiology_{species}")[sex].get("organs", {})
    return MappingProxyType(
        {k: float(v["mass"]["value"]) for k, v in organs.items() if "mass" in v}
    )


@cache
def _bodyweight(species: str, sex: str) -> float:
    from pkpdsim.database.tables import load

    return float(load(f"physiology_{species}")[sex]["bodyweight"]["value"])


def _largest(pieces: Sequence[np.ndarray]) -> np.ndarray:
    return max(pieces, key=lambda s: float(np.ptp(s[:, 0]) * np.ptp(s[:, 1])))


class AnimalBody:
    """A drawn animal: its outline, organs and the points of the routes (drawing units).

    Args:
        species: `mouse` or `rat`.
        sex: the sex, of the physiology of the library.
        drawing: the resource to draw, `DRAWINGS` of the species by default.
    """

    def __init__(
        self, species: str, sex: str = "male", drawing: Anatomy | None = None
    ) -> None:
        """Read the drawing and find the snout and the tail vein."""
        if drawing is None:
            if species not in DRAWINGS:
                raise ValueError(
                    f"No drawing of the {species}; the animals drawn are {sorted(DRAWINGS)}."
                )
            resource = DRAWINGS[species]
            drawing = anatomy(resource, sex if resource == "pig" else "male")
        self.species, self.sex, self.a = species, sex, drawing
        self.organs: dict[str, list[np.ndarray]] = {
            k: list(s.paths) for k, s in drawing.organs.items()
        }
        self.organs.update(
            {
                k: list(s.paths)
                for k, s in drawing.anatomy.items()
                if k in ("gallbladder", "bladder", "uterus")
            }
        )
        pts = drawing.outline
        self.lo, self.hi = pts.min(0), pts.max(0)
        span = self.hi - self.lo
        mid_x = drawing.midline
        near = pts[np.abs(pts[:, 0] - mid_x) < 0.08 * span[0]]
        # head up: the snout at the top on the midline, the tail vein below the hind legs
        self.snout = near[np.argmin(near[:, 1])] + np.array([0.0, 2.0])
        base = near[near[:, 1] > self.lo[1] + 0.82 * span[1]]
        self.tail = (
            base[np.argmin(base[:, 1])]
            if len(base)
            else np.array([mid_x, self.hi[1] - 0.12 * span[1]])
        )
        self.mid = (self.lo + self.hi) / 2
        lm = drawing.landmarks
        if "snout_x" in lm:  # the points of the routes given by the drawing (the pig)
            self.snout = np.array([lm["snout_x"], lm["snout_y"]])
            self.tail = np.array([lm["vein_x"], lm["vein_y"]])
        self.vein = VEINS.get(species, "tail vein")

    def centroid(self, key: str) -> np.ndarray:
        """The centre of the largest piece of an organ."""
        return _largest(self.organs[key]).mean(0)

    def drawn(self, key: str) -> bool:
        """Whether the drawing has the organ (or the route or sample of the key)."""
        return key in self.organs or key in ("blood", "iv", "urine", "lumen", "sk")

    def anchor(self, key: str, side: int) -> np.ndarray:
        """The point a label of a key points to."""
        if key in self.organs:
            big = _largest(self.organs[key])
            c = big.mean(0)
            band = big[np.abs(big[:, 1] - c[1]) < 1.0]
            if len(band) == 0:
                return c
            x = band[:, 0].min() if side < 0 else band[:, 0].max()
            return np.array([x - side * min(1.2, abs(x - c[0]) * 0.5), c[1]])
        if key in ("blood", "iv"):
            return self.tail
        if key == "urine":
            return (
                self.centroid("bladder") + np.array([1.0, 5.0])
                if "bladder" in self.organs
                else self.tail
            )
        if key == "lumen":
            return self.centroid("la") if "la" in self.organs else self.mid
        if key == "sk":
            pts = self.a.outline
            row = pts[np.abs(pts[:, 1] - self.mid[1]) < 2.0]
            x = row[:, 0].min() if side < 0 else row[:, 0].max()
            return np.array([x - side * 1.0, self.mid[1]])
        return np.array(
            [self.mid[0] + side * 0.25 * (self.hi[0] - self.lo[0]), self.mid[1]]
        )

    def svg(self, *, tags: bool = False) -> str:
        """The drawing [drawing units]."""
        a = self.a
        o = [
            f'<path d="{path([a.outline])}" fill="{SKINS.get(self.species, palette.SKIN)}" fill-opacity="0.7" stroke="{palette.INK}" '
            f'stroke-width="0.8" stroke-linejoin="round"/>'
        ]
        for k in _ORDER:
            if k in self.organs:
                tag = f' data-k="{k}"' if tags else ""
                o.append(
                    f'<path{tag} d="{path(self.organs[k])}" fill="{palette.ORGAN_COLORS.get(k, "#ccc")}" '
                    f'stroke="#3b2f28" stroke-width="0.25" fill-rule="evenodd"/>'
                )
        if a.details:
            o.append(f'<path d="{path(a.details)}" fill="{palette.INK}"/>')
        return "".join(o)

    def overlays(self, content: ModelContent, icon_scale: float) -> str:
        """The portal vein, the ureters, the oral route and the circled icons [drawing units]."""
        o = []
        if "li" in self.organs:
            liver = self.centroid("li")
            for g in ("st", "sp", "pa", "si", "la"):
                if g in self.organs:
                    p = self.centroid(g)
                    o.append(_curve(p, liver, palette.PORTAL, 0.45, bend=0.8))
                    o.append(
                        f'<circle cx="{p[0]:.2f}" cy="{p[1]:.2f}" r="0.4" fill="{palette.PORTAL}"/>'
                    )
        if content.urine and "ki" in self.organs and "bladder" in self.organs:
            bladder = self.centroid("bladder")
            for s in self.organs["ki"][:2]:
                o.append(_curve(s.mean(0), bladder, palette.URINE, 0.35, bend=0.8))
        k = icon_scale
        if content.has_route("oral"):
            target = self.centroid("st") if "st" in self.organs else self.mid
            o.append(
                _curve(
                    self.snout,
                    target,
                    palette.ORAL,
                    0.4,
                    bend=-1.5,
                    dash="0.8,0.6",
                    arrow=True,
                )
            )
            p = self.snout
            o.append(
                f'<g data-k="route:oral" transform="translate({p[0]:.2f},{p[1]:.2f}) scale({k:.4f})">{icons.tablet(0, 0, 22)}</g>'
            )
        p = self.tail
        o.append(
            f'<g data-k="sample" transform="translate({p[0] - 3.2:.2f},{p[1]:.2f}) scale({k:.4f})">{icons.drop(0, 0, 22)}</g>'
        )
        if content.has_route("iv"):
            o.append(
                f'<g data-k="route:iv" transform="translate({p[0] + 3.2:.2f},{p[1]:.2f}) scale({k:.4f})">{icons.syringe(0, 0, 22)}</g>'
            )
        return "".join(o)

    def caption(self) -> str:
        """The caption: the species, the sex and the body weight of the reference animal."""
        bw = _bodyweight(self.species, self.sex)
        weight = f"{bw:g} kg" if bw >= 1 else f"{bw * 1000:g} g"
        return f"{self.species} {self.sex}, {weight}"


def _species(definition: object) -> str:
    return str(getattr(definition, "species", "human"))


class AnimalFigure(AnatomyFigure):
    """The anatomical figure of the content of a model on a laboratory animal.

    Args:
        content: the content of the definition.
        species: the species of the build.
        sex: the sex of the build.
        compact: labels without the process lines.
        body: the drawn animal, `AnimalBody(species)` by default.
    """

    def __init__(
        self,
        content: ModelContent,
        species: str,
        *,
        sex: str = "male",
        compact: bool = False,
        body: AnimalBody | None = None,
    ) -> None:
        """Prepare the animal."""
        self.c = content
        self.compact = compact
        self.animal = body or AnimalBody(species, sex)
        self.size = (0, 0)

    def svg(
        self,
        *,
        title: bool = True,
        substances: bool = True,
        height: float = 560.0,
        text: float = 1.3,
        left: Sequence[str] = LEFT,
        right: Sequence[str] = RIGHT,
        tags: bool = False,
        **unused: object,
    ) -> str:
        """The figure as a standalone SVG; the animal fills `height` [px]."""
        a = self.animal
        wl = 330.0
        wcol, gap = wl * text, 34.0
        span = a.hi - a.lo
        scale = height / span[1]
        top = 112.0 if title else (70.0 if substances else 12.0)
        oy = top + (10 if (title or substances) else 0)
        ox = wcol + gap - a.lo[0] * scale
        oyy = oy - a.lo[1] * scale
        width = round(2 * wcol + 2 * gap + span[0] * scale)
        foot = oy + height
        key_h = 50 * text
        body = (
            f'<g transform="translate({ox:.1f},{oyy:.1f}) scale({scale:.4f})">'
            f"{a.svg(tags=tags)}{a.overlays(self.c, text / scale)}</g>"
        )

        def px(p: np.ndarray) -> np.ndarray:
            return np.array([ox + p[0] * scale, oyy + p[1] * scale])

        col_top = max(4.0, top - (20 if (title or substances) else 6))
        need = 0.0
        for side, keys in ((-1, left), (1, right)):
            hs = [
                self.label(k, 0, 0, wl - 20, side)[1] * text
                for k in keys
                if self.shown(k)
            ]
            need = max(need, sum(hs) + 8 * len(hs))
        h = round(max(foot + 30 * text + key_h + 8, col_top + need + key_h + 8))
        out = []
        for side, keys in ((-1, left), (1, right)):
            items = [
                (
                    k,
                    float(px(a.anchor(k, side))[1]),
                    self.label(k, 0, 0, wl - 20, side)[1] * text,
                )
                for k in keys
                if self.shown(k)
            ]
            for k, yy, _ in self.place(items, col_top, h - key_h - 8):
                lx = 0.0 if side < 0 else width - wcol
                g, _ = self.label(k, 6 if side < 0 else 14, 0, wl - 20, side)
                g = g.replace(
                    "arm vein: the observer Carm", f"{a.vein}: the observer Carm"
                )
                out.append(
                    f'<g data-k="{k}" transform="translate({lx:.1f},{yy:.1f}) scale({text})">{g}</g>'
                )
                p = px(a.anchor(k, side))
                ex = (wcol + 2) if side < 0 else (width - wcol + 6)
                ly = yy + 8 * text
                drawn = a.drawn(k)
                dash = "" if drawn else ' stroke-dasharray="3,3"'
                out.append(
                    f'<path d="M{p[0]:.1f},{p[1]:.1f} L{ex - side * 14:.1f},{ly:.1f} L{ex:.1f},{ly:.1f}" fill="none" '
                    f'stroke="#444" stroke-width="{0.8 * text:.2f}"{dash}/>'
                )
                if drawn:
                    out.append(
                        f'<circle cx="{p[0]:.1f}" cy="{p[1]:.1f}" r="{2.4 * text:.1f}" fill="#fff" stroke="#222" stroke-width="1"/>'
                    )
        cx = ox + (a.lo[0] + a.hi[0]) / 2 * scale
        captions = [
            f'<text x="{cx:.1f}" y="{foot + 14 * text:.1f}" class="bc" text-anchor="middle">{escape(a.caption())}</text>',
            f'<text x="{cx:.1f}" y="{foot + 28 * text:.1f}" class="bs" text-anchor="middle">'
            "the reference animal of the library</text>",
        ]
        header = []
        if title:
            header += [
                f'<text x="6" y="22" class="tt">{escape(self.c.name)}</text>',
                f'<text x="6" y="40" class="st">{escape(self.c.sid)} · {escape(a.species)} build</text>',
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
        key = f'<g transform="translate(0,{h - key_h + 10 * text:.1f}) scale({text})">{self._key()}</g>'
        self.size = (width, h)
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {h}" width="{width}" height="{h}" '
            f'font-family="{FONT}">{self._style(text)}{self._defs()}<rect width="{width}" height="{h}" fill="#fff"/>'
            f"{''.join(header)}{body}{''.join(out)}{''.join(captions)}{key}</svg>"
        )

    def _key(self) -> str:
        """The key of the lines, chips and icons (the animals have no layers)."""
        k = []
        kx = 6.0
        for color, label, dash in (
            (palette.PORTAL, "portal vein", ""),
            (palette.URINE, "urine", ""),
            (palette.ORAL, "oral route", "4,3"),
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
        for kind, label in (
            ("tablet", "oral gavage"),
            ("syringe", "intravenous dose"),
            ("drop", "blood sample"),
        ):
            k.append(
                f'{getattr(icons, kind)(kx + 9, ky - 4, 18)}<text x="{kx + 22:.1f}" y="{ky:.1f}" class="ky">{label}</text>'
            )
            kx += 34 + len(label) * 6.6 + 14
        return "".join(k)


def figure_of(
    definition: object, content: ModelContent, *, compact: bool = False
) -> AnatomyFigure:
    """The anatomical figure of a definition: the human builds or the laboratory animal of its species."""
    species = _species(definition)
    if species == "human":
        return AnatomyFigure(content, compact=compact)
    sex = str(getattr(definition, "sex", "male"))
    return AnimalFigure(content, species, sex=sex, compact=compact)
