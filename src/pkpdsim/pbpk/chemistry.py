"""The processes of the chemistry: metabolism, transport and custom reactions."""

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum

from sbmlutils.factory import Parameter, Species

from pkpdsim.pbpk.organs import Organ, Part
from pkpdsim.pbpk.values import Q

__all__ = [
    "MM",
    "CustomReaction",
    "Direction",
    "FirstOrder",
    "Hill",
    "Membrane",
    "Metabolism",
    "Transport",
]


@dataclass(frozen=True)
class FirstOrder:
    """First order kinetics on the unbound concentration.

    Attributes:
        k: the clearance of the whole organ (l/min), with `per_volume` the
            clearance per litre of organ (1/min).
        per_volume: `k` is per litre of organ, so the whole organ has the
            clearance `k * V<o>` and follows the organ of the individual, as a
            Vmax does; `False` keeps the clearance of the whole organ, which
            does not scale with the organ volume.
    """

    k: Q
    per_volume: bool = False


@dataclass(frozen=True)
class MM:
    """Michaelis-Menten kinetics on the unbound concentration.

    Attributes:
        vmax: the maximal rate of the whole organ (mmole/min/l of organ).
        km: the Michaelis constant (mmole/l).
    """

    vmax: Q
    km: Q


@dataclass(frozen=True)
class Hill:
    """Hill kinetics on the unbound concentration.

    Attributes:
        vmax: the maximal rate of the whole organ (mmole/min/l of organ).
        km: the concentration of the half-maximal rate (mmole/l).
        n: the Hill coefficient.
    """

    vmax: Q
    km: Q
    n: Q


class Direction(StrEnum):
    """The direction of a transport relative to the cell."""

    INFLUX = "influx"
    EFFLUX = "efflux"


class Membrane(StrEnum):
    """The membrane a transporter sits in."""

    BASOLATERAL = "basolateral"
    """The blood side."""
    APICAL = "apical"
    """The side of bile, urine or gut lumen."""


@dataclass(frozen=True)
class Metabolism:
    """The conversion of a substrate in an organ.

    Attributes:
        sid: the id of the process, unique in the definition.
        organ: the organ.
        substrate: the id of the substrate.
        products: the fractions of the substrate converted into each product
            (sum 1); no products is an elimination.
        kinetics: the rate law.
        protein: the id of the enzyme, its activity factor multiplies the rate.
        profile: the name of the abundance profile of the enzyme.
    """

    sid: str
    organ: Organ
    substrate: str
    products: Mapping[str, float]
    kinetics: FirstOrder | MM | Hill
    protein: str | None = None
    profile: str | None = None


@dataclass(frozen=True)
class Transport:
    """The active transport of a substance in an organ.

    Attributes:
        sid: the id of the process, unique in the definition.
        organ: the organ.
        substance: the id of the transported substance.
        direction: the direction relative to the cell.
        membrane: the membrane.
        kinetics: the rate law.
        protein: the id of the transporter.
    """

    sid: str
    organ: Organ
    substance: str
    direction: Direction
    membrane: Membrane
    kinetics: FirstOrder | MM | Hill
    protein: str | None = None


@dataclass(frozen=True)
class CustomReaction:
    """A mechanism the library has no class for, placed like any process.

    Attributes:
        sid: the id of the process, unique in the definition.
        organ: the organ.
        part: the part of the organ it acts in.
        equation: the reaction equation.
        formula: the rate law.
        parameters: the parameters of the rate law.
        species: additional species of the reaction.
    """

    sid: str
    organ: Organ
    part: Part
    equation: str
    formula: str
    parameters: list[Parameter]
    species: list[Species] = field(default_factory=list)
