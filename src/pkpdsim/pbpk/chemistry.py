"""The processes of the chemistry: metabolism, transport and custom reactions.

A process states its own constants, or takes its capacity as a ratio of
another process of the organ (`RateOf`); a metabolism converts its substrate
into its products by fractions, numbers or values with the remainder `REST`;
a custom reaction has the constants of its mechanism as `ModelParameter`s, its
own constants with provenance, and may declare its rate at the steady state
(`stationary`): a `Formula` of the rate, `"linear"` for a rate law which is
linear in its substrate, or `Saturable(vmax, km)` for a Michaelis-Menten rate
of a plasma species; the builder then writes the starts of its substances at
their steady state (`pkpdsim.pbpk.stationary`).

The enum fields of a process (`organ`, `direction`, `membrane`, `part`) may be
given as member, value or name (`"li"`, `"liver"`, `"basolateral"`): the
process holds the member, anything else raises a `ValueError` which names the
field and the allowed values.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum, StrEnum
from typing import Literal

import libsbml
from sbmlutils.factory import Parameter, Species

from pkpdsim import distrib
from pkpdsim.pbpk.organs import Organ, Part, as_member
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import (
    SID,
    Formula,
    ModelParameter,
    Q,
    Value,
    check_positive,
    check_unit,
    number,
)

__all__ = [
    "LINEAR",
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
    "Saturable",
    "Stationary",
    "Transport",
    "check_custom_parameters",
    "check_kinetics",
    "check_products",
    "check_rate_of",
    "check_stationary",
    "fraction_range",
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
class Saturable:
    """A Michaelis-Menten rate of a plasma species, `vmax*C/(km + C)`.

    As the steady state of a custom reaction (`CustomReaction.stationary`):
    its rate law is the saturable rate of the plasma species of its
    substrate `C`, so its organ starts at the root of a quadratic; the
    formulas are formulas of the organ, as the rate law (the elements of the
    organ, the parameters of the reaction, constants of the top model). As an
    argument of the formulas of `pkpdsim.pbpk.steady_state` (a saturable
    secretion of the kidney plasma) the formulas read ids of the top model.

    Attributes:
        vmax: the capacity of the whole organ (mmole/min), a formula.
        km: the Michaelis constant (mM), a formula.
    """

    vmax: str
    km: str


LINEAR = "linear"
"""The steady state of a custom reaction whose rate law is linear in its
substrate (`CustomReaction.stationary`): a first order clearance of a plasma
or tissue species."""

type Stationary = Formula | Literal["linear"] | Saturable
"""The steady state of a custom reaction (`CustomReaction.stationary`)."""


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
        stationary: the rate of the reaction at the steady state of its
            substances, which makes the builder start every organ part of
            them at their steady state (`pkpdsim.pbpk.stationary`), `None`
            for none: a `Formula` of the rate (mmole/min), a formula of the
            organ as the rate law (it reads the elements of the organ, the
            parameters of the reaction, the starts of species of other
            substances and constants of the top model, never a species of a
            substance it moves); `"linear"` (`LINEAR`), a rate law linear in
            its substrate, a plasma or tissue species of one substance, whose
            clearance the builder takes from the rate law; `Saturable(vmax,
            km)`, a Michaelis-Menten rate of the plasma species of its
            substrate. Every custom reaction which moves a substance with a
            steady state declares its own.
    """

    sid: str
    organ: Organ
    part: Part
    equation: str
    formula: str
    parameters: Sequence[ModelParameter | Parameter] = field(default_factory=list)
    species: list[Species] = field(default_factory=list)
    stationary: Stationary | None = None

    def __post_init__(self) -> None:
        """Normalise organ and part to their enum members."""
        object.__setattr__(
            self, "organ", as_member(Organ, self.organ, "CustomReaction.organ")
        )
        object.__setattr__(
            self, "part", as_member(Part, self.part, "CustomReaction.part")
        )


# --- the checks of the processes (`PBPKDefinition.validate`) -----------------


def check_kinetics(where: str, k: FirstOrder | MM | Hill) -> list[str]:
    """Check that the parameters of a rate law are positive.

    A `RateOf` is checked by `check_rate_of`.
    """
    problems: list[str] = []
    for name in ("k", "vmax", "km", "n"):
        q = getattr(k, name, None)
        if not isinstance(q, RateOf):
            problems += check_positive(where, name, q)
    return problems


def check_rate_of(
    p: Metabolism | Transport,
    chemistry: Mapping[str, Metabolism | Transport | CustomReaction],
) -> list[str]:
    """Check a constant of a process which is a ratio of another process, `RateOf`.

    The other process is a metabolism or transport of the same organ which
    states the same constant itself (a Vmax for a Vmax, a first order `k`
    with the same `per_volume` for a `k`, so that both have one unit), the
    factor is positive and dimensionless; a chain of ratios (and so a cycle)
    is refused.

    Args:
        p: the process.
        chemistry: the processes of the definition by id.

    Returns:
        The problems, naming the processes.
    """
    problems: list[str] = []
    for constant in ("vmax", "k"):
        q = getattr(p.kinetics, constant, None)
        if not isinstance(q, RateOf):
            continue
        problems += check_positive(p.sid, "RateOf factor", q.factor)
        problems += check_unit(p.sid, "RateOf factor", q.factor, U.dimensionless)
        other = chemistry.get(q.process)
        where = f"{p.sid}: RateOf '{q.process}'"
        if other is None:
            problems.append(f"{where} names no process of the chemistry")
            continue
        if not isinstance(other, Metabolism | Transport):
            problems.append(
                f"{where} is no metabolism or transport, a custom reaction has no "
                "constant of the library"
            )
            continue
        if other.organ is not p.organ:
            problems.append(
                f"{where} is a process of organ '{other.organ.value}', not of "
                f"'{p.organ.value}'"
            )
            continue
        theirs = getattr(other.kinetics, constant, None)
        if constant == "vmax":
            what = "Vmax"
        else:
            per_volume = isinstance(p.kinetics, FirstOrder) and p.kinetics.per_volume
            what = "first order constant k " + (
                "per litre of organ" if per_volume else "of the whole organ"
            )
            if (
                not isinstance(other.kinetics, FirstOrder)
                or other.kinetics.per_volume != per_volume
            ):
                theirs = None
        if theirs is None:
            problems.append(f"{where} has no {what}")
        elif isinstance(theirs, RateOf):
            problems.append(
                f"{where} is itself a ratio: the other process states its own constant"
            )
    return problems


def check_products(p: Metabolism) -> list[str]:
    """Check the product fractions of a metabolism: their sum is never above 1.

    Numbers sum to 1 (the stoichiometry of one reaction). A `Value` (fixed or
    fitted) is a constant of the organ which a fit, a scan or an experiment
    may change, so the products need the remainder `REST`, which keeps the
    sum at 1 for every value; the remainder is never negative: the sum of
    the upper ends of the ranges of the other fractions (`fraction_range`)
    is at most 1, and not always 1 (a remainder of 0). Every fraction is
    positive and its range not negative; a `Formula` cannot be bounded.

    Args:
        p: the metabolism.

    Returns:
        The problems, naming the metabolism and the products.
    """
    problems: list[str] = []
    rest = [k for k, f in p.products.items() if isinstance(f, Rest)]
    if len(rest) > 1:
        problems.append(
            f"{p.sid}: more than one product is the remainder REST: {sorted(rest)}"
        )
    elif rest and len(p.products) == 1:
        problems.append(f"{p.sid}: the remainder REST needs another product")
    low_sum = high_sum = total = 0.0
    values: list[str] = []
    ranges: list[str] = []
    for product, f in p.products.items():
        if isinstance(f, Rest):
            continue
        where = f"product fraction of '{product}'"
        if not isinstance(f, float | int | Value):
            problems.append(
                f"{p.sid}: {where} is a {type(f).__name__}, which cannot be "
                "bounded: a fraction is a number, a Value or REST"
            )
            continue
        problems += check_unit(p.sid, where, f, U.dimensionless)
        x = number(f)
        low, high = fraction_range(f)
        if x is None or not x > 0:
            problems.append(f"{p.sid}: {where} must be positive, got {x}")
        elif low < 0:
            ranges.append(
                f"{p.sid}: {where} can be {low:g} (the support of the "
                "distribution of its draws), a fraction is never negative"
            )
        total += x or 0.0
        low_sum += low
        high_sum += high
        if isinstance(f, Value):
            values.append(product)
    if problems:
        return problems
    problems += ranges
    if rest:
        if high_sum > 1.0 + 1e-9:
            problems.append(
                f"{p.sid}: the product fractions can sum to {high_sum:g} > 1 (the "
                "upper ends of the ranges of the fitted fractions), the remainder "
                f"REST of '{rest[0]}' would be negative"
            )
        elif low_sum >= 1.0 - 1e-9:
            problems.append(
                f"{p.sid}: the other fractions sum to {low_sum:g}, the remainder "
                f"REST of '{rest[0]}' is 0"
            )
    elif values:
        problems.append(
            f"{p.sid}: the product fraction of '{values[0]}' is a Value, a "
            "constant which a fit or an experiment may change: give one product "
            "the remainder REST (the sum stays 1), or the fixed fractions as "
            "numbers (the stoichiometry, the provenance in the notes)"
        )
    elif p.products and abs(total - 1.0) > 1e-9:
        problems.append(f"{p.sid}: product fractions sum to {total:g}, not 1")
    return problems


def fraction_range(f: float | Value) -> tuple[float, float]:
    """Get the range of a fraction: its value and the support of the draws of a fit.

    The distribution of the draws of a fitted value (the uncertainty of
    `fitted`, truncated to the bounds of the fit) is the range a fit or an
    uncertainty run gives the fraction (`pkpdsim.distrib.support`).
    """
    if not isinstance(f, Value):
        return float(f), float(f)
    low = high = f.value
    for u in f.uncertainties or []:
        for child in u.uncertParameters or []:
            if child.type != libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION or not child.math:
                continue
            lo, hi = distrib.support(child.math)
            low, high = min(low, lo), max(high, hi)
    return low, high


def check_stationary(p: CustomReaction) -> list[str]:
    """Check the declared steady state of a custom reaction.

    A `Formula` in mmole/min, the text `"linear"` or a `Saturable` with two
    formulas; the build checks what the declaration reads and whether the
    organs have a closed steady state (`pkpdsim.pbpk.stationary`).

    Args:
        p: the custom reaction.

    Returns:
        The problems, naming the reaction.
    """
    s = p.stationary
    if s is None or s == LINEAR:
        return []
    if isinstance(s, Formula):
        return check_unit(p.sid, "stationary", s, U.mmole_per_min)
    if isinstance(s, Saturable):
        return [
            f"{p.sid}: the stationary Saturable has an empty {name}"
            for name in ("vmax", "km")
            if not str(getattr(s, name)).strip()
        ]
    return [
        f"{p.sid}: stationary is a {type(s).__name__}, not a Formula of the rate "
        "(mmole/min), 'linear' or a Saturable"
    ]


def check_custom_parameters(p: CustomReaction, defined: set[str]) -> list[str]:
    """Check the parameters of a custom reaction, local ids of its organ.

    An SId without `__` (the separator of the flat ids), defined once in the
    organ (`defined`, the ids of the other custom reactions of the organ: they
    read it), a `ModelParameter` with a value in its unit. That an id is no
    other element of the organ is checked by the build
    (`OrganBuilder.model_definition`).

    Args:
        p: the custom reaction.
        defined: the parameter ids of the custom reactions of the organ so
            far, to which the ids of `p` are added.

    Returns:
        The problems, naming the reaction and the ids.
    """
    problems: list[str] = []
    for par in p.parameters:
        sid = str(par.sid)
        if not SID.fullmatch(sid):
            problems.append(f"{p.sid}: parameter id '{sid}' is not an SBML SId")
        elif "__" in sid:
            problems.append(
                f"{p.sid}: parameter id '{sid}' contains '__', the separator of "
                "the ids of the organs in the flat model"
            )
        elif sid in defined:
            problems.append(
                f"{p.sid}: parameter '{sid}' of organ '{p.organ.value}' is defined "
                "twice: one custom reaction defines it, the others of the organ "
                "read it"
            )
        defined.add(sid)
        if isinstance(par, ModelParameter):
            problems += check_unit(
                f"{p.sid}: parameter '{sid}'", "value", par.value, par.unit
            )
    return problems
