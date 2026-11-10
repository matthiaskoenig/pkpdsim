"""The gradients along the lobule: analytic profiles and the weight of a zone.

A protein, a transporter or any other quantity of a zonated liver varies along
the periportal to pericentral axis. A `Profile` declares that gradient as a
shape of the position, `0` periportal and `1` pericentral, and hands the engine
the level of every zone (`Profile.values`): a fit estimates the two parameters
of the shape instead of one level per zone, and the same declaration gives the
same gradient at any number of zones.

The levels are normalized: their mean over the volumes of the zones is 1, so
the total amount of the quantity in the liver does not change with its
distribution, whatever the shape. The volumes are those the engine builds, and
`pkpdsim.zonation.geometry.layer_volumes` splits every layer evenly over the
zones, so a profile without a weight normalizes to the arithmetic mean 1 over
the zones and follows the volumes if a later geometry makes them unequal. Only
the shape of a declaration survives the normalization, not the scale it is
written on: `Linear(1.0, 4.0)` is the four-fold pericentral gradient, and so is
`Linear(2.0, 8.0)`.

The area of the lobule is `lobular_weight`, which is passed in where it is the
weight: the figure of the lobule and the comparison with literature values
which are weighted by the area. It is the area of the idealized hexagonal
lobule at a position: the central vein in its centre, the six portal triads in
its corners, and the position of a point the relative distance
`d_pp / (d_pp + d_pc)` of the nearest triad and the central vein, as the
zonation coordinate of a histological section is measured. The geometry gives
the area in closed form: the points of one position lie on an Apollonius circle
of a triad and the centre, so the area up to a position is a single integral
over one twelfth of the hexagon (`lobular_area`, `_WEDGE_NODES`, 32 nodes of
Gauss-Legendre are exact to 1e-11) and `lobular_weight` is its derivative. The
area of a zone is the difference of the areas of its ends, so the lobular
weighting of the zones is the sequence `zone_areas(n)`, which is exact; a
density handed over as a callable is taken at the centre of every zone, which
for `lobular_weight` overstates the zone of its peak by 6.7 % as computed
(5.6 % after the levels are divided by their weighted mean, which is what a
profile does to them) and shifts a sharp profile by up to 1.5 % at 10 zones. The weight is 0 at both vessels,
peaks at
`1 / (1 + sqrt(3))` where the isoline leaves the hexagon through the edges, and
cuts the lobule in the ratio 2 to 1 at the bisectors, its 1/3 beyond the
position 0.5 being the inner hexagon of the perpendicular bisectors.

The area of the hexagonal lobule reproduces the normalization constants of the
linear patterns of the spt cell model (pkdb_models
`tpm/spt/models/spt_liver.py`, `protein = 0.5 * <constant> * <shape>`) to 1 %,
where a uniform weight is 14 % and the weight of a disc 24 % off: the constants
of spt were computed from a hexagonal lobule. Its two sharp patterns it misses
by 23.78 % and 11.08 %, and no weight of the lobule reproduces them (no
distance power, vessel radius, Hill exponent or threshold of this geometry
meets both); the `TODO.md` of spt records those patterns as too sharp and not
symmetrical, see `tests/zonation/test_profiles.py`.
"""

import math
import random
from abc import ABC, abstractmethod
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import pairwise
from typing import Protocol, runtime_checkable

import numpy as np

__all__ = [
    "Constant",
    "Exponential",
    "Linear",
    "Profile",
    "Random",
    "Sigmoid",
    "Weight",
    "Weights",
    "face_positions",
    "lobular_area",
    "lobular_edge",
    "lobular_radius",
    "lobular_weight",
    "normalized",
    "zone_areas",
    "zone_positions",
]

type Weight = Callable[[float], float] | Sequence[float] | None
"""The weight of the zones: a density of the position, the weights of the zones
themselves (`zone_areas`) or none, the zones of equal volume."""

_WEDGE = math.pi / 6.0
"""The angle of one twelfth of the hexagon, the wedge of one portal triad."""

_APOTHEM = math.sqrt(3.0) / 2.0
"""The apothem of the hexagonal lobule of the circumradius 1."""

_WEDGE_AREA = _APOTHEM**2 / math.sqrt(3.0)
"""The integral of the squared radius of the hexagon over the wedge."""

_WEDGE_NODES: tuple[tuple[float, float], ...] = tuple(
    zip(*(a.tolist() for a in np.polynomial.legendre.leggauss(32)), strict=True)
)
"""The nodes and weights of the Gauss-Legendre quadrature over the wedge."""


def _radius_of_the_edge(theta: float) -> float:
    """The radius of the hexagon of the circumradius 1 in the wedge direction."""
    return _APOTHEM / math.cos(_WEDGE - theta)


def _radius_of_the_isoline(ratio: float, theta: float) -> float:
    """The radius of a position in the wedge direction.

    The points of the position `x` lie on the Apollonius circle of the triad
    `(1, 0)` and the centre with the distance ratio `ratio = x / (1 - x)`,
    `r = 1 / (cos(theta) + sqrt(cos(theta)**2 + ratio**2 - 1))`; the root is
    the one which encloses the centre.

    Args:
        ratio: the ratio of the distances of the position.
        theta: the direction in the wedge of the triad (rad).

    Returns:
        The radius, `math.inf` where the circle misses the direction and the
        whole ray lies beyond the position.
    """
    under = math.cos(theta) ** 2 + ratio * ratio - 1.0
    if under < 0.0:
        return math.inf
    return 1.0 / (math.cos(theta) + math.sqrt(under))


def _direction_of_the_edge(ratio: float) -> float:
    """The direction in which a position leaves the hexagon through its edge.

    The Apollonius circle of the ratio meets the edge where
    `sin(theta) = sqrt(3) / 2 * ratio`; beyond that direction the position is
    outside the lobule and the whole ray lies beyond it.
    """
    return min(_WEDGE, math.asin(min(1.0, _APOTHEM * ratio)))


def lobular_edge(x: float) -> float:
    """Get the direction in which the isoline of a position leaves the hexagon.

    Up to this direction the isoline is the arc of the Apollonius circle of the
    position (`lobular_radius`); beyond it the position lies outside the lobule
    and the edge of the hexagon bounds it, so the isoline follows the straight
    edge. The isoline of 0 leaves at once, it is the edge over the whole wedge,
    and the isoline of 1, the central vein, never leaves. A renderer samples
    the arc and takes the rest straight (`pkpdsim.zonation.render`).

    Args:
        x: the lobular position, 0 periportal (the portal triads) and 1
            pericentral (the central vein).

    Returns:
        The direction in the wedge of a triad, in `[0, pi / 6]` (rad).

    Raises:
        ValueError: for a position outside the lobule.
    """
    if not 0.0 <= x <= 1.0:
        raise ValueError(f"The lobular position {x} is outside [0, 1].")
    if x == 1.0:
        return _WEDGE
    return _direction_of_the_edge(x / (1.0 - x))


def lobular_area(x: float) -> float:
    """Get the area fraction of the idealized hexagonal lobule up to a position.

    The periportal share of the lobule, 0 at the portal triads and 1 at the
    central vein. The area beyond the position is bounded by the isoline of the
    position up to the direction in which it leaves the hexagon and by the edge
    of the hexagon beyond it; `lobular_weight` is the derivative of the area
    and `zone_areas` the areas of the zones.

    Args:
        x: the lobular position, 0 periportal (the portal triads) and 1
            pericentral (the central vein).

    Returns:
        The fraction of the area of the lobule with a position below `x`.

    Raises:
        ValueError: for a position outside the lobule.
    """
    if not 0.0 <= x <= 1.0:
        raise ValueError(f"The lobular position {x} is outside [0, 1].")
    if x == 0.0:
        return 0.0
    if x == 1.0:
        return 1.0
    ratio = x / (1.0 - x)
    edge = _direction_of_the_edge(ratio)
    half = 0.5 * edge
    inside = 0.0
    for node, weight in _WEDGE_NODES:
        theta = half * (node + 1.0)
        inside += weight * half * _radius_of_the_isoline(ratio, theta) ** 2
    # beyond `edge` the hexagon bounds the area, and the integral of its
    # squared radius over the rest of the wedge is `apothem**2 * tan`
    beyond = inside + _APOTHEM**2 * math.tan(_WEDGE - edge)
    return 1.0 - beyond / _WEDGE_AREA


def lobular_radius(x: float, theta: float) -> float:
    """Get the radius at which the lobule has a position in a direction.

    The idealized hexagon has its centre (the central vein, position 1) at the
    origin and the circumradius 1; `theta` is the direction measured from the
    nearest portal triad, in the wedge `[0, pi / 6]` (a twelfth of the
    hexagon). The points of a position lie on an Apollonius circle of the
    triad and the centre, and the circle is cut off by the edge of the
    hexagon: the area up to this radius and beyond is what `lobular_area`
    integrates. A point at a smaller radius has a larger position.

    Args:
        x: the lobular position, 0 periportal and 1 pericentral.
        theta: the direction in the wedge of the triad (rad).

    Returns:
        The radius of the position in the direction, the radius of the edge of
        the hexagon where the isoline has left the lobule (always for the
        position 0), 0 for the position 1.

    Raises:
        ValueError: for a position outside the lobule, or a direction outside
            the wedge.
    """
    if not 0.0 <= x <= 1.0:
        raise ValueError(f"The lobular position {x} is outside [0, 1].")
    if not 0.0 <= theta <= _WEDGE:
        raise ValueError(f"The direction {theta} is outside the wedge [0, pi/6].")
    if x == 1.0:
        return 0.0
    edge = _radius_of_the_edge(theta)
    if x == 0.0:
        return edge
    return min(edge, _radius_of_the_isoline(x / (1.0 - x), theta))


def lobular_weight(x: float) -> float:
    """Get the area of the idealized hexagonal lobule at a lobular position.

    The density of the area over the position (the derivative of
    `lobular_area`): the area of the lobule between two positions is its
    integral, and over the whole lobule it is 1.

    Args:
        x: the lobular position, 0 periportal (the portal triads) and 1
            pericentral (the central vein).

    Returns:
        The area of the lobule at the position, 0 at both vessels.

    Raises:
        ValueError: for a position outside the lobule.
    """
    if not 0.0 <= x <= 1.0:
        raise ValueError(f"The lobular position {x} is outside [0, 1].")
    if x in (0.0, 1.0):
        return 0.0
    ratio = x / (1.0 - x)
    edge = _direction_of_the_edge(ratio)
    # the weight is minus the derivative of the area beyond the position: that
    # area is the integral of the squared radius of the isoline over the wedge,
    # the radius moves with the ratio by `dr / dratio = ratio * r**3 /
    # (r cos(theta) - 1)` and the ratio with the position by `1 / (1 - x)**2`;
    # beyond `edge` the isoline has left the lobule and does not move
    half = 0.5 * edge
    total = 0.0
    for node, weight in _WEDGE_NODES:
        theta = half * (node + 1.0)
        r = _radius_of_the_isoline(ratio, theta)
        total += weight * half * r**4 / (r * math.cos(theta) - 1.0)
    return -2.0 * ratio * total / (_WEDGE_AREA * (1.0 - x) ** 2)


def zone_positions(n: int) -> tuple[float, ...]:
    """Get the lobular positions of the centres of the zones.

    Args:
        n: the number of zones.

    Returns:
        The position `(k - 0.5) / n` of every zone, numbered from 1
        (periportal).

    Raises:
        ValueError: for less than one zone.
    """
    if n < 1:
        raise ValueError(f"A lobule has at least one zone, not {n}.")
    return tuple((k - 0.5) / n for k in range(1, n + 1))


def face_positions(n: int) -> tuple[float, ...]:
    """Get the lobular positions of the faces of the zones.

    The `n + 1` faces of `n` zones, `j / n` from the periportal face of the
    first zone (0) to the pericentral face of the last (1): the positions at
    which the chain evaluates its flows, where `zone_positions` gives the
    centres at which it evaluates the volumes, the profiles and the
    concentrations.

    Args:
        n: the number of zones.

    Returns:
        The position `j / n` of every face, `j` from 0 (periportal).

    Raises:
        ValueError: for less than one zone.
    """
    if n < 1:
        raise ValueError(f"A lobule has at least one zone, not {n}.")
    return tuple(j / n for j in range(n + 1))


def zone_areas(n: int) -> tuple[float, ...]:
    """Get the area fraction of the lobule of every zone.

    The exact area of a zone, the difference of the areas of its ends
    (`lobular_area`), and not the weight at its centre: at 10 zones the density
    at the centre overstates the zone which holds the peak of the weight by
    6.7 % as computed, 5.6 % after the normalisation of a profile.

    Args:
        n: the number of zones.

    Returns:
        The fraction of the area of the lobule of every zone, numbered from 1
        (periportal), which sum to 1.

    Raises:
        ValueError: for less than one zone.
    """
    if n < 1:
        raise ValueError(f"A lobule has at least one zone, not {n}.")
    edges = [lobular_area(k / n) for k in range(n + 1)]
    return tuple(b - a for a, b in pairwise(edges))


def _zone_weights(n: int, weight: Weight) -> tuple[float, ...]:
    """The weight of every zone, its share of the volume of the liver.

    Args:
        n: the number of zones.
        weight: `None` for zones of equal volume, a sequence of the weights of
            the `n` zones, or a density of the lobular position, which is taken
            at the centre of every zone.

    Returns:
        The weight of every zone.

    Raises:
        ValueError: for a sequence which is not of the number of zones.
    """
    if weight is None:
        return (1.0,) * n
    if isinstance(weight, Sequence):
        weights = tuple(weight)
        if len(weights) != n:
            raise ValueError(f"{len(weights)} weights for the {n} zones of the lobule.")
        return weights
    return tuple(weight(x) for x in zone_positions(n))


def normalized(
    raw: Sequence[float], n: int, weight: Weight = None
) -> tuple[float, ...]:
    """Normalize the levels of the zones to the volume-weighted mean 1.

    Args:
        raw: the level of every zone on the scale of its declaration.
        n: the number of zones.
        weight: the weight of the zones, their share of the volume of the
            liver. `None` weights the zones of equal volume the engine builds
            (the arithmetic mean, see
            `pkpdsim.zonation.geometry.layer_volumes`). A sequence is used as
            it is given, one weight per zone: `zone_areas(n)` is the exact area
            of the lobule and the lobular weighting. A callable is a density of
            the lobular position and is taken at the centre of every zone,
            always, so `lobular_weight` as a callable samples the density and
            overstates the zone which holds its peak by 6.7 % at 10 zones as
            computed, 5.6 % after the normalisation below.

    Returns:
        The levels divided by their weighted mean, so that a liver with the
        gradient holds the amount of a liver without one.

    Raises:
        ValueError: for a number of levels or of weights which is not the
            number of zones, or for levels whose weighted mean is not positive.
    """
    if len(raw) != n:
        raise ValueError(f"{len(raw)} levels for the {n} zones of the lobule.")
    weights = _zone_weights(n, weight)
    mean = sum(v * w for v, w in zip(raw, weights, strict=True)) / sum(weights)
    if not mean > 0.0:
        raise ValueError(f"The levels {tuple(raw)} have the weighted mean {mean}.")
    return tuple(v / mean for v in raw)


@runtime_checkable
class Profile(Protocol):
    """A gradient of a quantity along the lobule."""

    def values(self, n: int, weight: Weight = None) -> tuple[float, ...]:
        """The level of every zone, normalized to the weighted mean 1.

        Args:
            n: the number of zones.
            weight: the weight of the zones, as of `normalized`: `None` for
                the zones of equal volume the engine builds, a sequence used
                as it is given (`zone_areas(n)` for the exact area of the
                lobule) or a density taken at the centre of every zone
                (`lobular_weight`, which samples it).
        """
        ...


class _Pointwise(ABC):
    """A profile which is a shape of the lobular position."""

    @abstractmethod
    def _shape(self, x: float) -> float:
        """The level at a lobular position, on the scale of the declaration."""

    def values(self, n: int, weight: Weight = None) -> tuple[float, ...]:
        """The shape at the centre of every zone, normalized.

        Args:
            n: the number of zones.
            weight: the weight of the zones, as of `normalized`: `None` for
                the zones of equal volume the engine builds, a sequence used
                as it is given (`zone_areas(n)` for the exact area of the
                lobule) or a density taken at the centre of every zone
                (`lobular_weight`, which samples it).
        """
        return normalized([self._shape(x) for x in zone_positions(n)], n, weight)


@dataclass(frozen=True)
class Constant(_Pointwise):
    """No gradient, the same level in every zone.

    Attributes:
        level: the level, which the normalization divides out.
    """

    level: float = 1.0

    def __post_init__(self) -> None:
        """Check the level.

        Raises:
            ValueError: if it is not positive.
        """
        if not self.level > 0.0:
            raise ValueError(f"Constant: the level {self.level} is not positive.")

    def _shape(self, x: float) -> float:
        """The level."""
        return self.level


@dataclass(frozen=True)
class Linear(_Pointwise):
    """A linear gradient between the two ends of the lobule.

    Attributes:
        pp: the periportal level (position 0).
        pc: the pericentral level (position 1).
    """

    pp: float
    pc: float

    def __post_init__(self) -> None:
        """Check the ends.

        Raises:
            ValueError: for a negative level or for two levels of 0.
        """
        if self.pp < 0.0 or self.pc < 0.0 or self.pp + self.pc <= 0.0:
            raise ValueError(f"Linear: the levels {self.pp} and {self.pc} are no ends.")

    def _shape(self, x: float) -> float:
        """The level at a position."""
        return self.pp + (self.pc - self.pp) * x


@dataclass(frozen=True)
class Exponential(_Pointwise):
    """An exponential gradient from the periportal end.

    Attributes:
        level: the periportal level (position 0), which the normalization
            divides out.
        decay: the decay over the lobule, `level * exp(-decay * x)`; a negative
            decay rises pericentrally.
    """

    level: float
    decay: float

    def __post_init__(self) -> None:
        """Check the level.

        Raises:
            ValueError: if it is not positive.
        """
        if not self.level > 0.0:
            raise ValueError(f"Exponential: the level {self.level} is not positive.")

    def _shape(self, x: float) -> float:
        """The level at a position."""
        return self.level * math.exp(-self.decay * x)


@dataclass(frozen=True)
class Sigmoid(_Pointwise):
    """A sharp gradient, high in one half of the lobule and low in the other.

    The Hill shape `1 / (1 + (x / midpoint)**steepness)` of the position, the
    sharp patterns of the spt cell model: `Sigmoid(0.25, 10.0)` is its
    periportal pattern and `Sigmoid(0.25, 10.0, reverse=True)` its pericentral
    one.

    Attributes:
        midpoint: the position of the half level, measured from the periportal
            end, or from the pericentral end of a reversed shape.
        steepness: the Hill exponent of the step.
        reverse: mirror the shape, low periportal and high pericentral.
    """

    midpoint: float
    steepness: float
    reverse: bool = False

    def __post_init__(self) -> None:
        """Check the midpoint and the steepness.

        Raises:
            ValueError: for a midpoint outside the lobule or a steepness which
                is not positive.
        """
        if not 0.0 < self.midpoint < 1.0:
            raise ValueError(f"Sigmoid: the midpoint {self.midpoint} is no position.")
        if not self.steepness > 0.0:
            raise ValueError(
                f"Sigmoid: the steepness {self.steepness} is not positive."
            )

    def _shape(self, x: float) -> float:
        """The level at a position."""
        position = 1.0 - x if self.reverse else x
        return 1.0 / (1.0 + (position / self.midpoint) ** self.steepness)


@dataclass(frozen=True)
class Weights:
    """A measured gradient, one level per zone.

    The levels are those of a lobule of as many zones: a profile of a shape
    serves any number of zones, these levels only the zones they were measured
    in. The field is `levels`, since `values` is the method of a profile.

    Attributes:
        levels: the level of every zone, from periportal to pericentral.
    """

    levels: Sequence[float]

    def __post_init__(self) -> None:
        """Check the levels and keep them as a tuple.

        Raises:
            ValueError: for no level, a negative level or levels which are all 0.
        """
        levels = tuple(self.levels)
        if not levels or any(v < 0.0 for v in levels) or sum(levels) <= 0.0:
            raise ValueError(f"Weights: the levels {levels} are no profile.")
        object.__setattr__(self, "levels", levels)

    def values(self, n: int, weight: Weight = None) -> tuple[float, ...]:
        """The measured levels, normalized.

        Args:
            n: the number of zones, which the levels are those of.
            weight: the weight of the zones, as of `normalized`: `None` for
                the zones of equal volume the engine builds, a sequence used
                as it is given (`zone_areas(n)` for the exact area of the
                lobule) or a density taken at the centre of every zone
                (`lobular_weight`, which samples it).

        Returns:
            The levels divided by their weighted mean.

        Raises:
            ValueError: for a lobule of another number of zones.
        """
        return normalized(self.levels, n, weight)


@dataclass(frozen=True)
class Random:
    """A random gradient, the variation of a protein between the zones.

    The levels are lognormal around 1 with the coefficient of variation of the
    declaration, drawn from the seed and the number of zones, so a build
    repeats them.

    Attributes:
        seed: the seed of the draw.
        cv: the coefficient of variation of the levels.
    """

    seed: int
    cv: float

    def __post_init__(self) -> None:
        """Check the coefficient of variation.

        Raises:
            ValueError: if it is negative.
        """
        if self.cv < 0.0:
            raise ValueError(f"Random: the cv {self.cv} is negative.")

    def values(self, n: int, weight: Weight = None) -> tuple[float, ...]:
        """The drawn levels, normalized.

        Args:
            n: the number of zones.
            weight: the weight of the zones, as of `normalized`: `None` for
                the zones of equal volume the engine builds, a sequence used
                as it is given (`zone_areas(n)` for the exact area of the
                lobule) or a density taken at the centre of every zone
                (`lobular_weight`, which samples it).
        """
        rng = random.Random(self.seed)
        sigma = math.sqrt(math.log(1.0 + self.cv**2))
        raw = [math.exp(sigma * rng.gauss(0.0, 1.0)) for _ in range(n)]
        return normalized(raw, n, weight)
