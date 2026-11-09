"""The gradients along the lobule: analytic profiles and the weight of a zone.

A protein, a transporter or any other quantity of a zonated liver varies along
the periportal to pericentral axis. A `Profile` declares that gradient as a
shape of the position, `0` periportal and `1` pericentral, and hands the engine
the level of every zone (`Profile.values`): a fit estimates the two parameters
of the shape instead of one level per zone, and the same declaration gives the
same gradient at any number of zones.

The levels are normalized: their mean over the weight of the zones is 1, so the
total amount of the quantity in the liver is the one of a liver without a
gradient, whatever the shape. Only the shape of a declaration survives the
normalization, not the scale it is written on: `Linear(1.0, 4.0)` is the
four-fold pericentral gradient, and so is `Linear(2.0, 8.0)`.

The weight of a position is `lobular_weight`, the area of the idealized
hexagonal lobule at that position: the central vein in its centre, the six
portal triads in its corners, and the position of a point the relative distance
`d_pp / (d_pp + d_pc)` of the nearest triad and the central vein, as the
zonation coordinate of a histological section is measured. The weight is the
area of the lobule between two positions, which the geometry gives in closed
form: the points of one position lie on an Apollonius circle of a triad and the
centre, so the area beyond a position is a single integral over one twelfth of
the hexagon, and the weight is its derivative (`_WEDGE_NODES`, 32 nodes of
Gauss-Legendre are exact to 1e-11). The weight is 0 at both vessels, peaks at
`1 / (1 + sqrt(3))` where the isoline leaves the hexagon through the edges, and
cuts the lobule in the ratio 2 to 1 at the bisectors, its 1/3 beyond the
position 0.5 being the inner hexagon of the perpendicular bisectors.

The hexagonal lobule reproduces the normalization constants of the linear
patterns of the spt cell model (pkdb_models `tpm/spt/models/spt_liver.py`,
`protein = 0.5 * <constant> * <shape>`) to 1 %, where a uniform weight is 14 %
and the weight of a disc 24 % off. Its two sharp patterns it misses by 24 % and
11 %, and no weight of this geometry carries all four constants: the sharp
constants of spt are not the moments of the weight of its linear ones (the
`TODO.md` of spt asks for the sharp patterns to be fixed), see
`tests/zonation/test_profiles.py`.
"""

import math
import random
from abc import ABC, abstractmethod
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

__all__ = [
    "Constant",
    "Exponential",
    "Linear",
    "Profile",
    "Random",
    "Sigmoid",
    "Weights",
    "lobular_weight",
    "normalized",
    "zone_positions",
]

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


def lobular_weight(x: float) -> float:
    """Get the area of the idealized hexagonal lobule at a lobular position.

    The density of the area over the position: the area of the lobule between
    two positions is its integral, and over the whole lobule it is 1.

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


def normalized(
    raw: Sequence[float], n: int, weight: Callable[[float], float] = lobular_weight
) -> tuple[float, ...]:
    """Normalize the levels of the zones to the volume-weighted mean 1.

    Args:
        raw: the level of every zone on the scale of its declaration.
        n: the number of zones.
        weight: the weight of a lobular position, the area of the lobule.

    Returns:
        The levels divided by their weighted mean, so that a liver with the
        gradient holds the amount of a liver without one.

    Raises:
        ValueError: for a number of levels which is not the number of zones,
            or for levels whose weighted mean is not positive.
    """
    if len(raw) != n:
        raise ValueError(f"{len(raw)} levels for the {n} zones of the lobule.")
    weights = [weight(x) for x in zone_positions(n)]
    mean = sum(v * w for v, w in zip(raw, weights, strict=True)) / sum(weights)
    if not mean > 0.0:
        raise ValueError(f"The levels {tuple(raw)} have the weighted mean {mean}.")
    return tuple(v / mean for v in raw)


@runtime_checkable
class Profile(Protocol):
    """A gradient of a quantity along the lobule."""

    def values(self, n: int, weight: Callable[[float], float]) -> tuple[float, ...]:
        """The level of every zone, normalized to the weighted mean 1.

        Args:
            n: the number of zones.
            weight: the weight of a lobular position.
        """
        ...


class _Pointwise(ABC):
    """A profile which is a shape of the lobular position."""

    @abstractmethod
    def _shape(self, x: float) -> float:
        """The level at a lobular position, on the scale of the declaration."""

    def values(self, n: int, weight: Callable[[float], float]) -> tuple[float, ...]:
        """The shape at the centre of every zone, normalized.

        Args:
            n: the number of zones.
            weight: the weight of a lobular position.
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

    def values(self, n: int, weight: Callable[[float], float]) -> tuple[float, ...]:
        """The measured levels, normalized.

        Args:
            n: the number of zones, which the levels are those of.
            weight: the weight of a lobular position.

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

    def values(self, n: int, weight: Callable[[float], float]) -> tuple[float, ...]:
        """The drawn levels, normalized.

        Args:
            n: the number of zones.
            weight: the weight of a lobular position.
        """
        rng = random.Random(self.seed)
        sigma = math.sqrt(math.log(1.0 + self.cv**2))
        raw = [math.exp(sigma * rng.gauss(0.0, 1.0)) for _ in range(n)]
        return normalized(raw, n, weight)
