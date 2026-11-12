"""The body of an individual of the library, drawn from the anatomy resources.

An `Individual` is a sex, an age, a body weight and a height (the reference
individual of the library by default), a functional liver and a gestational
progress; its `Quantities` are the numbers of the library: the volumes of the
organs (`reference_individual` scaled as the build scales the lean tissues,
`body_composition`), the waist and the visceral and subcutaneous adipose tissue
(`visceral_composition`).

A `Body` draws it in cm, the midline at x 0 and the top of the head at y 0:

- the outline is the drawing of the sex mapped by the field of `field.py`:
  the proportions of the age in y, the trunk stretched by the waist, the limbs
  and the face by the body mass index, the pregnant belly;
- the layers on a raster of 0.25 cm: the subcutaneous adipose tissue as a band
  under the skin, its thickness the subcutaneous volume over the body surface
  (twice on the trunk), the muscle inside it, the thorax and the abdomen;
- the organs keep the forms of the drawing (a uniform scale and a shift only):
  their positions follow the box of their cavity, their size is the cube root
  of the volume over the volume of the reference individual of the sex, so the
  organs of the reference fill their cavity but a share (`FREE_REF`) and the
  space between them is visceral adipose tissue (abdomen) or mediastinum
  (thorax); the kidneys are moved out of the shadow of liver and spleen and
  the heart onto the diaphragm, so every organ shows in a frontal view.
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from types import MappingProxyType
from xml.sax.saxutils import escape

import contourpy
import numpy as np
from scipy.ndimage import (
    binary_erosion,
    distance_transform_edt,
    gaussian_filter,
    map_coordinates,
)

from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.visual import palette
from pkpdsim.pbpk.visual.field import RES, Drawing, Field, Gaussian, drawing, raster
from pkpdsim.sex import Sex

__all__ = [
    "CAVITY_ORGANS",
    "FREE_MIN",
    "FREE_REF",
    "Body",
    "Individual",
    "Quantities",
    "body_strip",
    "path",
    "reference_body",
]

RES_CM = 0.25
"""The raster of the layers of a body [cm per pixel]."""

CAVITY_ORGANS: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        "thorax": ("lu", "he"),
        "abdomen": (
            "la",
            "si",
            "li",
            "st",
            "sp",
            "ki",
            "pa",
            "gallbladder",
            "bladder",
            "uterus",
        ),
    }
)
"""The organs of each cavity, in the order they are drawn (later over earlier)."""

FREE_REF: Mapping[str, float] = MappingProxyType({"thorax": 0.10, "abdomen": 0.24})
"""The share of a cavity of the reference individual which no organ covers."""

FREE_MIN = 0.04
"""The least free share of a cavity: larger organs (a child) are scaled down together."""

_LAYOUT = {"ki": (1.4, 6.5), "he": (0.0, 4.0), "sp": (0.6, 0.0)}
_LAYOUT_SCALE = {"he": 1.2}

# the proportions over age, fractions of the height: the head (top to chin) and
# the legs (crotch to floor), and the waist over the height of the children
_AGES = np.array([0, 1, 3, 6, 10, 15, 20, 120])
_HEAD = np.array([0.25, 0.22, 0.195, 0.17, 0.155, 0.14, 0.13, 0.13])
_LEG = np.array([0.30, 0.36, 0.42, 0.45, 0.47, 0.475, 0.475, 0.475])
_WHTR = np.array([0.65, 0.58, 0.52, 0.47, 0.45, 0.43, 0.44, 0.44])


def path(subpaths: Sequence[np.ndarray], close: bool = True) -> str:
    """The SVG path data of polylines [cm], to a millimetre."""
    z = "Z" if close else ""
    return "".join(
        "M" + " ".join(f"{x:.1f},{y:.1f}" for x, y in s) + z for s in subpaths
    )


def _area(subpaths: Sequence[np.ndarray]) -> float:
    total = 0.0
    for s in subpaths:
        x, y = s[:, 0], s[:, 1]
        total += 0.5 * abs(float(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))
    return total


def _trace(
    mask: np.ndarray, origin: np.ndarray, sigma: float = 0.9
) -> list[np.ndarray]:
    """The smooth closed outlines of a mask, in the coordinates of the raster."""
    if not mask.any():
        return []
    m = gaussian_filter(np.pad(mask.astype(float), 2), sigma)
    lines = contourpy.contour_generator(z=m, line_type="Separate").lines(0.5)
    out = []
    for c in lines:
        c = np.asarray(c)
        if len(c) < 6:
            continue
        pts = np.c_[
            (c[:, 0] - 2) * RES_CM + origin[0], (c[:, 1] - 2) * RES_CM + origin[1]
        ]
        out.append(pts[::3] if len(pts) > 60 else pts)
    return out


@dataclass(frozen=True)
class Quantities:
    """The numbers of the library of an individual.

    Attributes:
        sex: the sex.
        age: the age [years].
        bw: the body weight [kg].
        height: the height [cm].
        bmi: the body mass index [kg/m2].
        bsa: the body surface area [m2].
        waist: the waist circumference [cm].
        vat: the visceral adipose tissue [kg].
        sat: the subcutaneous adipose tissue [kg].
        volumes: the volume of every organ by its key [l].
    """

    sex: Sex
    age: float
    bw: float
    height: float
    bmi: float
    bsa: float
    waist: float
    vat: float
    sat: float
    volumes: Mapping[str, float]


@dataclass(frozen=True)
class Individual:
    """An individual of the library to draw.

    Attributes:
        sex: the sex.
        age: the age [years].
        bw: the body weight [kg], the reference individual of the age by default.
        height: the height [cm], the reference individual of the age by default.
        liver: the functional liver relative to the healthy one (the factor of
            a Child-Pugh class).
        pregnant: the progress of a pregnancy, 0 none to 1 at term.
        label: the caption.
    """

    sex: Sex | str
    age: float = 35.0
    bw: float | None = None
    height: float | None = None
    liver: float = 1.0
    pregnant: float = 0.0
    label: str = ""

    def quantities(self) -> Quantities:
        """The numbers of the library: volumes, waist, adipose depots."""
        return _quantities(self)


@cache
def _quantities(ind: Individual) -> Quantities:
    from pkpdsim.pbpk.adipose import visceral_composition
    from pkpdsim.pbpk.physiology import body_composition, reference_individual

    sex = Sex(ind.sex)
    r = reference_individual(sex, age=ind.age)
    bw = ind.bw if ind.bw is not None else r.bw
    height = ind.height if ind.height is not None else r.height
    c = body_composition(sex, bw, height, ind.age)
    rc = body_composition(sex, r.bw, r.height, ind.age)
    volumes = {}
    for organ, v in r.volumes.items():
        if organ is Organ.ADIPOSE:
            continue
        if organ is Organ.BRAIN:
            k = 1.0
        elif organ is Organ.SKIN:
            k = c.bsa / rc.bsa
        else:
            k = c.sv_tis / rc.sv_tis
        volumes[organ.value] = v * k * (ind.liver if organ is Organ.LIVER else 1.0)
    vc = visceral_composition(sex, bw, height, ind.age)
    return Quantities(
        sex=sex,
        age=ind.age,
        bw=bw,
        height=height,
        bmi=c.bmi,
        bsa=c.bsa,
        waist=vc.waist,
        vat=vc.visceral,
        sat=vc.subcutaneous,
        volumes=MappingProxyType(volumes),
    )


def _levels(d: Drawing) -> dict[str, float]:
    lm = d.landmarks
    leg = lm["bottom"] - lm["crotch"]
    return {
        "neck": lm["neck"],
        "shoulder": lm["neck"] + 0.35 * (lm["armpit"] - lm["neck"]),
        "chest": lm["armpit"] + 2,
        "waist": lm["waist"],
        "belly": lm["waist"] + 0.3 * (lm["crotch"] - lm["waist"]),
        "hip": lm["crotch"] - 0.2 * (lm["crotch"] - lm["waist"]),
        "crotch": lm["crotch"],
        "thigh": lm["crotch"] + 0.13 * leg,
    }


def _field(ind: Individual, q: Quantities) -> Field:
    """The field of an individual: the proportions of its age, its stretches, its belly."""
    sex = Sex(ind.sex)
    d = drawing(sex.value)
    lm = d.landmarks
    y = _levels(d)
    ref = _quantities(Individual(sex))
    bmi_ref = ref.bmi
    h0 = lm["bottom"] - lm["top"]
    head = (
        float(np.interp(ind.age, _AGES, _HEAD)) * (lm["neck"] - lm["top"]) / h0 / 0.13
    )
    leg = (
        float(np.interp(ind.age, _AGES, _LEG))
        * (lm["bottom"] - lm["crotch"])
        / h0
        / 0.475
    )
    adult = ind.age >= 18
    rw = q.waist / ref.waist - 1 if adult else 0.0
    rb = q.bmi / bmi_ref - 1 if adult else 0.0
    female = sex is Sex.FEMALE
    uniform = 400.0  # a sigma far beyond the drawing: a uniform stretch
    trunk: list[Gaussian] = [
        (y["chest"], 16, 0.55 * rw),
        (y["waist"], 16, 0.75 * rw),
        (y["hip"], 16, (0.8 if female else 0.6) * rw),
        (y["neck"] + 3, 4, 0.4 * rb),
    ]
    arm: list[Gaussian] = [(0, uniform, 0.55 * rb)]
    legs: list[Gaussian] = [
        (y["thigh"], 30, (0.45 if female else 0.35) * rb),
        (0, uniform, 0.25 * rb),
    ]
    if ind.age < 18:
        w = float(np.interp(ind.age, _AGES, _WHTR) / _WHTR[-2]) - 1
        c = max(0.0, 1 - ind.age / 5)
        trunk += [(y["waist"], 16, w), (y["shoulder"], 8, -0.12 * (1 - ind.age / 18))]
        arm.append((0, uniform, 0.2 * c))
        legs.append((0, uniform, 0.2 * c))
    e = float(np.clip((ind.age - 60) / 30, 0, 1))
    if e > 0:
        trunk.append((y["shoulder"], 8, -0.05 * e))
        arm.append((0, uniform, -0.07 * e))
        legs.append((0, uniform, -0.07 * e))
    belly: Gaussian | None = None
    if ind.pregnant > 0:
        p = ind.pregnant
        belly = (
            y["waist"] + 0.18 * (y["crotch"] - y["waist"]),
            6.5 + 2.5 * p,
            0.32 * p,
        )
    return Field(
        d,
        q.height,
        head,
        leg,
        trunk=trunk,
        arm=arm,
        leg=legs,
        head_w=1 + 0.08 * max(rb, -0.5),
        belly=belly,
    )


@cache
def _cavity_raster(sex: str) -> np.ndarray:
    """The cavities on the raster of the drawing: 1 thorax, 2 abdomen."""
    d = drawing(sex)
    out = np.zeros(d.mask.shape, np.int8)
    out[raster(d.anatomy.cavities.get("abdomen", ()), d.mask.shape)] = 2
    out[raster(d.anatomy.cavities.get("thorax", ()), d.mask.shape)] = 1
    return out


@dataclass(frozen=True, eq=False)
class _Layout:
    """The organs of a cavity: their pieces (drawing), positions (cm) and the cavity mask."""

    pieces: dict[str, list[np.ndarray]]
    positions: dict[str, list[np.ndarray]]
    cavity: np.ndarray


class Body:
    """The drawing of an individual [cm].

    Args:
        individual: the individual.
        reference: the body of the reference individual of the sex, whose
            organs set the scale; `None` for the reference itself and, for
            every other individual, the cached `reference_body` of its sex.
    """

    def __init__(self, individual: Individual, reference: Body | None = None) -> None:
        """Compute the outline, the layers and the organs."""
        self.individual = individual
        self.sex = Sex(individual.sex).value
        if reference is None and not _is_reference(individual):
            reference = reference_body(self.sex)
        self.reference = reference
        self.quantities = individual.quantities()
        self.drawing = drawing(self.sex)
        self.anatomy = self.drawing.anatomy
        self.field = _field(individual, self.quantities)
        self.outline = self.field(self.anatomy.outline)
        self._layers()
        self._place()

    # --- the layers on a raster -----------------------------------------------
    def _layers(self) -> None:
        q, f, d = self.quantities, self.field, self.drawing
        lo = self.outline.min(0) - 1
        hi = self.outline.max(0) + 1
        self.origin = lo
        shape = (int((hi[1] - lo[1]) / RES_CM) + 1, int((hi[0] - lo[0]) / RES_CM) + 1)
        self.shape = shape
        inside = raster([self.outline], shape, (lo[0], lo[1]), RES_CM)
        dist = distance_transform_edt(inside) * RES_CM
        gy, gx = np.mgrid[0 : shape[0], 0 : shape[1]]
        xr, yr = f.inverse((gx * RES_CM + lo[0]).ravel(), (gy * RES_CM + lo[1]).ravel())
        r, c = yr / RES, xr / RES
        w_trunk = map_coordinates(
            d.weights["trunk"], [r, c], order=1, mode="nearest"
        ).reshape(shape)
        w_head = map_coordinates(
            d.weights["head"], [r, c], order=1, mode="nearest"
        ).reshape(shape)
        cav = map_coordinates(
            _cavity_raster(self.sex).astype(float), [r, c], order=0, mode="constant"
        ).reshape(shape)
        mean = q.sat / 0.92 / q.bsa / 10.0  # l/m2 -> cm
        t_sat = (
            np.minimum(mean * (0.8 + 1.2 * w_trunk), 3.2) * (1 - w_head) + 0.25 * w_head
        )
        t_wall = 1.1 * q.height / 176.0
        self.sat = inside & (dist <= t_sat)
        deep = inside & (dist > t_sat + t_wall)
        self.cavities = {"thorax": deep & (cav == 1), "abdomen": deep & (cav == 2)}
        self.muscle = inside & ~self.sat & (w_head < 0.5)
        self.inner = binary_erosion(inside & ~self.sat, iterations=2)

    # --- the organs ---------------------------------------------------------------
    def _components(self, key: str) -> list[np.ndarray]:
        shape = self.anatomy.shape(key)
        if shape is None:
            return []
        subs = [s for s in shape.paths if len(s) > 2]
        cx = self.anatomy.midline
        if key in _LAYOUT:
            dx, dy = _LAYOUT[key]
            subs = [s + np.array([dx * np.sign(s[:, 0].mean() - cx), dy]) for s in subs]
        if key in _LAYOUT_SCALE:
            c = np.vstack(subs).mean(0)
            subs = [(s - c) * _LAYOUT_SCALE[key] + c for s in subs]
        areas = [_area([s]) for s in subs]
        return [s for s, a in zip(subs, areas, strict=True) if a > 0.04 * max(areas)]

    def size_ratio(self, key: str) -> float:
        """The linear size of an organ relative to the reference individual, (V/V_ref)^(1/3)."""
        ref = self.reference
        if ref is None:
            return 1.0
        v, v_ref = self.quantities.volumes, ref.quantities.volumes
        if key in v and key in v_ref:
            return (v[key] / v_ref[key]) ** (1 / 3)
        r = self.quantities.height / ref.quantities.height
        if key == "uterus":
            r *= (1.0 + 40.0 * self.individual.pregnant**1.6) ** (1 / 3)
        return r

    def _place(self) -> None:
        self.organs = {}
        self.cavity_paths: dict[str, list[np.ndarray]] = {}
        self.k_cav: dict[str, float] = {}
        self.free: dict[str, float] = {}
        for cav, keys in CAVITY_ORGANS.items():
            m = self.cavities[cav]
            if not m.any():
                continue
            a_cav = m.sum() * RES_CM**2
            comps = {k: self._components(k) for k in keys}
            comps = {k: v for k, v in comps.items() if v}
            centres = {
                k: [self.field(s.mean(0)[None, :])[0] for s in v]
                for k, v in comps.items()
            }
            warped = np.vstack([self.field(s) for v in comps.values() for s in v])
            bx0, by0 = warped.min(0)
            bx1, by1 = warped.max(0)
            ys, xs = np.nonzero(m)
            cx0, cx1 = (
                xs.min() * RES_CM + self.origin[0],
                xs.max() * RES_CM + self.origin[0],
            )
            cy0, cy1 = (
                ys.min() * RES_CM + self.origin[1],
                ys.max() * RES_CM + self.origin[1],
            )
            sx, sy = (cx1 - cx0) / (bx1 - bx0), (cy1 - cy0) / (by1 - by0)
            sx, sy = min(sx, 1.12 * sy), min(sy, 1.12 * sx)
            cb = np.array([(bx0 + bx1) / 2, (by0 + by1) / 2])
            cc = np.array([(cx0 + cx1) / 2, (cy0 + cy1) / 2])
            pos = {k: [cc + (p - cb) * [sx, sy] for p in v] for k, v in centres.items()}

            layout = _Layout(comps, pos, m)
            if self.reference is None:
                scale = self._solve(layout, (1 - FREE_REF[cav]) * a_cav, 0.1, 10.0)
            else:
                scale = self.reference.k_cav[cav]
                cap = (1 - FREE_MIN) * a_cav
                if self._covered(layout, scale) > cap:
                    scale = self._solve(layout, cap, 0.1 * scale, scale)
            self.k_cav[cav] = scale
            for k, v in self._shapes(layout, scale).items():
                self.organs[k] = [self._inside(piece) for piece in v]
            self.cavity_paths[cav] = _trace(m, self.origin)
            self.free[cav] = a_cav - self._covered(layout, scale)

    def _shapes(self, layout: _Layout, scale: float) -> dict[str, list[np.ndarray]]:
        """The organs of a cavity at a scale: every piece scaled about its centre and moved to its position."""
        return {
            k: [
                (c - c.mean(0)) * scale * self.size_ratio(k) + p
                for c, p in zip(v, layout.positions[k], strict=True)
            ]
            for k, v in layout.pieces.items()
        }

    def _covered(self, layout: _Layout, scale: float) -> float:
        """The area of a cavity which its organs cover at a scale [cm2]."""
        polys = [c for v in self._shapes(layout, scale).values() for c in v]
        im = raster(polys, self.shape, (self.origin[0], self.origin[1]), RES_CM)
        return float((im & layout.cavity).sum() * RES_CM**2)

    def _solve(self, layout: _Layout, target: float, lo: float, hi: float) -> float:
        """The scale at which the organs cover a target area (bisection on a log scale)."""
        for _ in range(30):
            mid = math.sqrt(lo * hi)
            if self._covered(layout, mid) < target:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)

    def _inside(self, piece: np.ndarray) -> np.ndarray:
        """A piece of an organ shrunk about its centre until it lies under the subcutaneous fat."""
        c = piece.mean(0)
        for _ in range(20):
            col = np.clip(
                ((piece[:, 0] - self.origin[0]) / RES_CM).round().astype(int),
                0,
                self.shape[1] - 1,
            )
            row = np.clip(
                ((piece[:, 1] - self.origin[1]) / RES_CM).round().astype(int),
                0,
                self.shape[0] - 1,
            )
            if self.inner[row, col].all():
                break
            piece = (piece - c) * 0.97 + c
        return piece

    # --- geometry -----------------------------------------------------------------
    def bounds(self) -> tuple[float, float, float]:
        """The left and right edge and the height of the body [cm]."""
        return (
            float(self.outline[:, 0].min()),
            float(self.outline[:, 0].max()),
            self.field.height,
        )

    def at(self, x: float, y: float) -> np.ndarray:
        """A point of the drawing in the body [cm]."""
        return self.field(np.array([[x, y]], float))[0]

    def centroid(self, key: str) -> np.ndarray:
        """The centre of the largest piece of an organ [cm]."""
        if key in self.organs:
            return max(self.organs[key], key=lambda s: _area([s])).mean(0)
        if key == "br":
            return self.at(*self.anatomy.organs["br"].paths[0].mean(0))
        shape = self.anatomy.shape(key)
        if shape is None:
            raise KeyError(key)
        return self.field(np.vstack(shape.paths).mean(0)[None, :])[0]

    def pieces(self, key: str, n: int = 2) -> list[np.ndarray]:
        """The centres of the largest pieces of an organ (the two kidneys) [cm]."""
        big = sorted(
            (s for s in self.organs.get(key, []) if _area([s]) > 4.0),
            key=lambda s: -_area([s]),
        )
        return [s.mean(0) for s in big[:n]]

    def toward(self, key: str, side: int) -> np.ndarray:
        """A point inside an organ towards a side (-1 left, 1 right of the figure) [cm]."""
        if key not in self.organs:
            return self.centroid(key)
        big = max(self.organs[key], key=lambda s: _area([s]))
        c = big.mean(0)
        band = big[np.abs(big[:, 1] - c[1]) < 1.5]
        if len(band) == 0:
            return c
        x = band[:, 0].min() if side < 0 else band[:, 0].max()
        return np.array([x - side * min(2.0, abs(x - c[0]) * 0.5), c[1]])

    def edge_point(self, y: float, side: int, inset: float = 0.6) -> np.ndarray:
        """A point inside the outline at the height y of the drawing, at the outer edge of a side [cm]."""
        yy = self.at(self.anatomy.midline, y)[1]
        near = self.outline[np.abs(self.outline[:, 1] - yy) < 0.8]
        x = near[:, 0].max() if side > 0 else near[:, 0].min()
        return np.array([x - side * inset, yy])

    @property
    def landmarks(self) -> dict[str, float]:
        """The landmarks of the drawing (units of the drawing)."""
        return self.drawing.landmarks

    # --- drawing -------------------------------------------------------------------
    def _arteries(self) -> list[np.ndarray]:
        """The arteries: aortic arch, carotids, arms and legs along their axes [cm]."""
        d, lm, cx = self.drawing, self.landmarks, self.anatomy.midline
        ys = d.ys
        out = []
        top_aorta = lm["armpit"] - 0.4 * (lm["armpit"] - lm["neck"])
        bottom_aorta = lm["waist"] + 0.35 * (lm["crotch"] - lm["waist"])
        out.append(
            np.array(
                [
                    [cx, top_aorta],
                    [cx + 0.6, (top_aorta + bottom_aorta) / 2],
                    [cx, bottom_aorta],
                ]
            )
        )
        jaw = lm["top"] + 0.8 * (lm["neck"] - lm["top"])
        hand = lm["crotch"] + 0.06 * (lm["bottom"] - lm["crotch"])
        foot = lm["bottom"] - 0.03 * (lm["bottom"] - lm["top"])
        for s in (-1, 1):
            out.append(
                np.array(
                    [
                        [cx + s * 0.8, top_aorta],
                        [cx + s * 2.2, lm["neck"]],
                        [cx + s * 2.0, jaw],
                    ]
                )
            )
            rows = (ys >= lm["armpit"]) & (ys <= hand)
            arm = np.c_[d.arm_axis[s][rows], ys[rows]][::8]
            start = np.array(
                [[cx + s * 1.0, top_aorta], [cx + s * 6.0, lm["shoulder"] + 3]]
            )
            out.append(np.vstack([start, arm]))
            rows = (ys >= lm["crotch"]) & (ys <= foot)
            leg = np.c_[d.leg_axis[s][rows], ys[rows]][::8]
            out.append(np.vstack([[[cx, bottom_aorta]], leg]))
        return [self.field(a) for a in out]

    def _pregnancy_line(self) -> str:
        p = self.individual.pregnant
        if p <= 0.05:
            return ""
        d, y = self.drawing, _levels(self.drawing)
        half = float(np.interp(y["belly"] - 2, d.ys, d.half))
        yc = y["belly"] - 2
        hw = half * (0.75 + 0.25 * p)
        ry = (y["crotch"] - y["waist"]) * (0.55 + 0.35 * p)
        t = np.linspace(math.radians(18), math.radians(162), 50)
        pts = self.field(np.c_[d.cx + hw * np.cos(t), yc - ry * 0.55 + ry * np.sin(t)])
        return f'<path d="{path([pts], close=False)}" fill="none" stroke="{palette.INK}" stroke-width="0.34" stroke-linecap="round"/>'

    def svg(
        self,
        *,
        fills: Mapping[str, str] | None = None,
        skin: str = palette.SKIN,
        layers: bool = True,
        arteries: bool = True,
        tags: bool = False,
        stroke: float = 0.22,
    ) -> str:
        """The SVG elements of the body [cm], for a group which places and scales them.

        Args:
            fills: colours by organ key over `palette.ORGAN_COLORS`.
            skin: the colour of the skin.
            layers: draw the subcutaneous adipose band, the muscle and the bones.
            arteries: draw the arteries.
            tags: mark every organ with `data-k="<key>"` (the interactive figure).
            stroke: the width of the outlines of the organs [cm].

        Returns:
            The elements, deterministic.
        """
        colors = {**palette.ORGAN_COLORS, **(fills or {})}
        a = self.anatomy

        def tag(k: str) -> str:
            return f' data-k="{escape(k)}"' if tags else ""

        o = [f'<path d="{path([self.outline])}" fill="{skin}" fill-opacity="0.6"/>']
        if layers:
            o.append(
                f'<path d="{path(_trace(self.muscle, self.origin, 1.0))}" fill="{palette.MUSCLE}" '
                f'fill-opacity="0.25" fill-rule="evenodd"/>'
            )
            o.append(
                f'<path d="{path(_trace(self.sat, self.origin, 1.0))}" fill="{palette.FAT}" '
                f'fill-opacity="0.75" fill-rule="evenodd"/>'
            )
            if "bo" in a.organs:
                cx = a.midline
                subs = list(a.organs["bo"].paths)
                both = subs + [np.c_[2 * cx - s[:, 0], s[:, 1]] for s in subs]
                o.append(
                    f'<path{tag("bo")} d="{path([self.field(s) for s in both])}" fill="{palette.BONE}" '
                    f'fill-opacity="0.9" stroke="#b9ad93" stroke-width="0.15"/>'
                )
        if arteries:
            o.append(
                f'<path d="{path(self._arteries(), close=False)}" fill="none" stroke="{palette.BLOOD}" '
                f'stroke-width="0.35" stroke-opacity="0.7" stroke-linecap="round" stroke-linejoin="round"/>'
            )
        o.append(self._pregnancy_line())
        for cav, bg in (("thorax", palette.THORAX), ("abdomen", palette.VAT)):
            if cav in self.cavity_paths:
                o.append(
                    f'<path d="{path(self.cavity_paths[cav])}" fill="{bg}" fill-rule="evenodd"/>'
                )
        for k in (*CAVITY_ORGANS["thorax"], *CAVITY_ORGANS["abdomen"]):
            if k in self.organs:
                o.append(
                    f'<path{tag(k)} d="{path(self.organs[k])}" fill="{colors.get(k, "#ccc")}" '
                    f'stroke="#3b2f28" stroke-width="{stroke}" stroke-opacity="0.85" stroke-linejoin="round"/>'
                )
        brain = a.organs["br"].paths
        c = np.vstack(brain).mean(0)
        cc = self.field(c[None, :])[0]
        kb = (
            self.reference.field.kh if self.reference is not None else self.field.kh
        ) * self.size_ratio("br")
        kb = min(kb, self.field.kh * 1.05)
        o.append(
            f'<path{tag("br")} d="{path([(s - c) * kb + cc for s in brain])}" fill="{colors["br"]}" '
            f'stroke="#3b2f28" stroke-width="{stroke}" fill-rule="evenodd"/>'
        )
        if self.sex == "male" and "go" in a.organs:
            o.append(
                f'<path{tag("go")} d="{path([self.field(s) for s in a.organs["go"].paths])}" '
                f'fill="{colors["go"]}" stroke="#3b2f28" stroke-width="{stroke}"/>'
            )
        if self.sex == "female" and "breast" in a.anatomy and self.individual.age >= 13:
            o.append(
                f'<path d="{path([self.field(s) for s in a.anatomy["breast"].paths])}" fill="none" '
                f'stroke="{palette.INK}" stroke-width="0.25" stroke-opacity="0.6"/>'
            )
        o.append(
            f'<path d="{path([self.outline])}" fill="none" stroke="{palette.INK}" stroke-width="0.6" '
            f'stroke-linejoin="round"/>'
        )
        if a.details:
            o.append(
                f'<path d="{path([self.field(s) for s in a.details])}" fill="{palette.INK}" fill-opacity="0.9"/>'
            )
        return "".join(o)


def _is_reference(ind: Individual) -> bool:
    return (
        ind.age == 35.0
        and ind.bw is None
        and ind.height is None
        and ind.liver == 1.0
        and ind.pregnant == 0.0
    )


@cache
def reference_body(sex: str) -> Body:
    """The body of the reference individual of a sex (cached), the scale of the organs."""
    return Body(Individual(Sex(sex)))


def body_strip(
    individuals: Sequence[Individual], scale: float = 2.4, font: str = "sans-serif"
) -> str:
    """Bodies side by side on one scale, the feet aligned, each with its caption.

    Args:
        individuals: the individuals, their `label` the first line of the caption.
        scale: px per cm.
        font: the font family of the captions.

    Returns:
        A standalone SVG.
    """
    items = []
    for ind in individuals:
        b = reference_body(Sex(ind.sex).value) if _is_reference(ind) else Body(ind)
        q = b.quantities
        x0, x1, h = b.bounds()
        lines = (
            ind.label or f"{Sex(ind.sex).value}, {ind.age:g} years",
            f"{q.height:.0f} cm, {q.bw:.0f} kg, liver {q.volumes['li']:.2f} l",
            f"VAT {q.vat:.1f} kg, SAT {q.sat:.0f} kg",
        )
        items.append((b.svg(), x0, x1, h, lines))
    top = max(h for _, _, _, h, _ in items)
    x = 0.0
    parts = []
    for g, x0, x1, h, lines in items:
        slot = max((x1 - x0) * scale, 170.0)
        ox = x + (slot - (x1 - x0) * scale) / 2 - x0 * scale
        parts.append(
            f'<g transform="translate({ox:.1f},{(top - h) * scale + 4:.1f}) scale({scale})">{g}</g>'
        )
        for j, t in enumerate(lines):
            cls = "a" if j == 0 else "b"
            parts.append(
                f'<text x="{x + slot / 2:.1f}" y="{top * scale + 22 + 14 * j:.1f}" text-anchor="middle" class="{cls}">{escape(t)}</text>'
            )
        x += slot + 14
    width, height = round(x), round(top * scale + 58)
    style = f"<style>.a{{font:600 12px {font};fill:#1f2a2e}}.b{{font:10.5px {font};fill:#5b6b70}}</style>"
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">'
        f'{style}<rect width="{width}" height="{height}" fill="#fff"/>{"".join(parts)}</svg>'
    )
