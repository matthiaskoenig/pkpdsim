"""The processes of the chemistry: metabolism, transport and custom reactions.

A process states its own constants, or takes its capacity as a ratio of
another process of the organ (`RateOf`); a metabolism converts its substrate
into its products by fractions, numbers or values with the remainder `REST`;
a custom reaction has the constants of its mechanism as `ModelParameter`s, its
own constants with provenance.

The enum fields of a process (`organ`, `direction`, `membrane`, `part`) may be
given as member, value or name (`"li"`, `"liver"`, `"basolateral"`): the
process holds the member, anything else raises a `ValueError` which names the
field and the allowed values.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum, StrEnum

from sbmlutils.factory import Parameter, Species

from pkpdsim.pbpk.organs import Organ, Part, as_member
from pkpdsim.pbpk.values import ModelParameter, Q, Value

__all__ = [
    "MM",
    "REST",
    "CustomReaction",
    "Direction",
    "FirstOrder",
    "Hill",
    "Membrane",
    "Metabolism",
    "RateOf",
    "Rest",
    "Transport",
]


@dataclass(frozen=True)
class RateOf:
    """The constant of a process as a ratio of the constant of another process.

    The `vmax` of `MM` or `Hill` or the `k` of `FirstOrder` of a metabolism or
    transport, given as a factor times the same constant of another
    metabolism or transport of the same organ (two transporters of one
    substrate with a fixed ratio of their capacities): the constant
    `<process>_Vmax` (or `_k`) is the initial assignment
    `<process>_Vmax_ratio * <other>_Vmax` of the organ, so it follows a change
    of the other constant before the simulation (a fit of one capacity moves
    both), and the factor is the constant `<process>_Vmax_ratio` with its
    provenance (`ids.rate_ratio`).

    Attributes:
        process: the id of the other process, a metabolism or transport of
            the organ which states the same constant itself (a `vmax` for a
            `vmax`, a `k` with the same `per_volume` for a `k`).
        factor: the ratio (dimensionless), a number, a `Value` or a `Formula`.
    """

    process: str
    factor: Q


class Rest(Enum):
    """The remainder of the product fractions of a metabolism, `REST`."""

    REST = "rest"


REST = Rest.REST
"""The product of a metabolism which takes the remainder of the fractions:
`1 - sum` of the fractions of the other products."""

type Fraction = float | Value | Rest
"""The fraction of a product of a metabolism: a number, a value (fitted) or the
remainder `REST`."""


@dataclass(frozen=True)
class FirstOrder:
    """First order kinetics on the unbound concentration.

    Attributes:
        k: the clearance per litre of organ (1/min), with `per_volume=False`
            the clearance of the whole organ (l/min); `RateOf` a ratio of the
            `k` of another process of the organ.
        per_volume: `k` is per litre of organ (the default), so the whole organ
            has the clearance `k * V<o>` and follows the organ of the
            individual, as the `vmax` of `MM` and `Hill` does; `False` keeps
            the clearance of the whole organ, which does not scale with the
            organ volume.
    """

    k: Q | RateOf
    per_volume: bool = True


@dataclass(frozen=True)
class MM:
    """Michaelis-Menten kinetics on the unbound concentration.

    Attributes:
        vmax: the maximal rate per litre of organ (mmole/min/l), the whole
            organ has the maximal rate `vmax * V<o>`; `RateOf` a ratio of the
            `vmax` of another process of the organ.
        km: the Michaelis constant (mmole/l).
    """

    vmax: Q | RateOf
    km: Q


@dataclass(frozen=True)
class Hill:
    """Hill kinetics on the unbound concentration.

    Attributes:
        vmax: the maximal rate per litre of organ (mmole/min/l), the whole
            organ has the maximal rate `vmax * V<o>`; `RateOf` a ratio of the
            `vmax` of another process of the organ.
        km: the concentration of the half-maximal rate (mmole/l).
        n: the Hill coefficient.
    """

    vmax: Q | RateOf
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
        products: the fractions of the substrate converted into each product;
            no products is an elimination. Numbers which sum to 1 are the
            stoichiometry of one reaction. A `Value` (a fitted fraction) or
            the remainder `REST` (one product, `1 - sum` of the others) makes
            every fraction a constant of the organ, `f_<process>_<product>`
            (`ids.product_fraction`), the reaction consumes the substrate and
            one reaction per product forms it with its fraction of the rate,
            `<process>_to_<product>` (`ids.product_reaction`); the fractions
            can never sum above 1: the range of a fitted fraction is its
            value and the support of the distribution of its draws (the
            bounds of the fit), a fitted fraction needs the remainder.
        kinetics: the rate law.
        protein: the id of the enzyme, its activity factor multiplies the rate.
        profile: the name of the abundance profile of the enzyme.
    """

    sid: str
    organ: Organ
    substrate: str
    products: Mapping[str, Fraction]
    kinetics: FirstOrder | MM | Hill
    protein: str | None = None
    profile: str | None = None

    def __post_init__(self) -> None:
        """Normalise the organ to an `Organ` member."""
        object.__setattr__(
            self, "organ", as_member(Organ, self.organ, "Metabolism.organ")
        )


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

    def __post_init__(self) -> None:
        """Normalise organ, direction and membrane to their enum members."""
        object.__setattr__(
            self, "organ", as_member(Organ, self.organ, "Transport.organ")
        )
        object.__setattr__(
            self,
            "direction",
            as_member(Direction, self.direction, "Transport.direction"),
        )
        object.__setattr__(
            self, "membrane", as_member(Membrane, self.membrane, "Transport.membrane")
        )


@dataclass(frozen=True)
class CustomReaction:
    """A mechanism the library has no class for, placed like any process.

    Attributes:
        sid: the id of the process, unique in the definition.
        organ: the organ.
        part: the part of the organ it acts in.
        equation: the reaction equation.
        formula: the rate law.
        parameters: the constants of the rate law, local ids of the organ
            (`<ORGAN>__<id>` in the flat model), none by default: a
            `ModelParameter` (a number, a `Value` with its provenance and the
            uncertainty of a fit, a `Formula` over ids of the top model), or a
            sbmlutils `Parameter` as it is. An id is defined once in an organ,
            the other reactions of the organ read it.
        species: additional species of the reaction.
    """

    sid: str
    organ: Organ
    part: Part
    equation: str
    formula: str
    parameters: Sequence[ModelParameter | Parameter] = field(default_factory=list)
    species: list[Species] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Normalise organ and part to their enum members."""
        object.__setattr__(
            self, "organ", as_member(Organ, self.organ, "CustomReaction.organ")
        )
        object.__setattr__(
            self, "part", as_member(Part, self.part, "CustomReaction.part")
        )
