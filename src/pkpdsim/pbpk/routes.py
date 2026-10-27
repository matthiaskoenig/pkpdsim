"""The routes of administration and their dose states.

A dose is a dose state in mg, which an experiment sets or adds to
(`PKPDSimulationExperiment.add_doses`: `"PODOSE_<s> + 10.0"`), and which is
applied by a first order process, the conventions of the hand-written models:

- intravenous: `IVDOSE_<s>` drains into the venous plasma `Cve_<s>` with the
  rate `Kiv_<s> = ln(2)/ti_<s>` (`ti_<s>` the half-life of the injection in s);
  a constant infusion `Ri_<s>` (mg/min) adds to the dose state;
- oral tablet: `PODOSE_<s>` is the undissolved drug in the stomach (the
  first segment of the intestine), which dissolves into the stomach lumen
  with the dissolution rate of the formulation and empties undissolved with
  the stomach; in the next segments the undissolved drug is a species of
  its own, it dissolves there and leaves undissolved into the feces;
- oral solution: `SOLDOSE_<s>` enters the stomach lumen dissolved with the
  fast rate `Ka_sol_<s>` (1/hr).

The dose in mg is converted to mmole with the molecular weight `Mr_<s>`. The
dose states and the parameters of the formulations are parameters of the
top model (`route_objects`); the oral doses are applied by the intestine,
which reads them through ports and owns their reactions and the rate rules
of the oral dose states (the top dose states replace its inputs, so the
rules act on them).

The dissolution rate `kdis` (1/min) acts on the undissolved amount of every
segment:

- `FirstOrderDissolution(k, fed, release)`: `kdis = k` (`Ka_dis_<s>`, 1/hr),
  times the factor of a modified-release formulation `f_dis_release_<s>`
  and `(1 + FED*(f_dis_fed_<s> - 1))` of the fed state, of those which are set;
- `WeibullDissolution(t50, shape, lag)`: the fraction dissolved
  `F = 1 - exp(-ln2*((t - lag)/t50)^shape)` after the lag time, applied as
  its hazard `kdis = ln2*shape/t50*((t - lag)/t50)^(shape - 1)` on the
  undissolved amount, with the time since the dose `t` the mean age of the
  undissolved drug (`age_dis_<s>`, the state `W_dis_<s>` of the
  age-weighted amount divided by the undissolved amount): exact for one
  dose; doses which overlap share one clock, the amount-weighted mean age:
  both dissolve with the hazard of the mean age, the rest of the earlier
  dose as if it were younger (not at all while the mean age is below the
  lag time), the new dose as if it were older; a dose given after the
  earlier one has dissolved is exact again;
- `DoseDependentDissolution(k, dose50)`: `kdis = k*D50/(D50 + D)` with the
  undissolved amount `D` (mg) of all segments, a dissolution which saturates
  with the undissolved dose (solubility-limited), at most `k*D50`.
"""

import math
from dataclasses import dataclass
from typing import Literal

from sbmlutils.factory import (
    AssignmentRule,
    Parameter,
    RateRule,
    Reaction,
    Sbase,
    UnitDefinition,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate
from pkpdsim.pbpk.data import value
from pkpdsim.pbpk.organs import Pool
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, check_literals

__all__ = [
    "ORAL",
    "DoseDependentDissolution",
    "FirstOrderDissolution",
    "Formulation",
    "Route",
    "WeibullDissolution",
    "dissolution_inputs",
    "dissolution_rate",
    "dissolution_reads_fed",
    "is_weibull",
    "route_objects",
]


@dataclass(frozen=True)
class FirstOrderDissolution:
    """The first order dissolution of a tablet.

    Attributes:
        k: the rate of the dissolution (1/hr).
        fed: the factor of the dissolution rate in the fed state (`FED = 1`),
            the body parameter `f_dis_fed_<s>`; the rate is
            `k*(1 + FED*(fed - 1))`. No parameter if `None`. The single lumen
            has no fed state and raises.
        release: the factor of a modified-release formulation, the body
            parameter `f_dis_release_<s>` (the given value, set by the
            experiments of such a formulation; 1 for the immediate-release
            tablet), which multiplies the rate. No parameter if `None`.
    """

    k: Q
    fed: Q | None = None
    release: Q | None = None


@dataclass(frozen=True)
class WeibullDissolution:
    """The Weibull dissolution of a tablet.

    The fraction dissolved `1 - exp(-ln2*((t - lag)/t50)^shape)` after the
    lag time.

    Attributes:
        t50: the time of half dissolution after the lag time (min).
        shape: the shape of the Weibull function (1: first order).
        lag: the lag time (min).
    """

    t50: Q
    shape: Q
    lag: Q = 0.0


@dataclass(frozen=True)
class DoseDependentDissolution:
    """The dissolution of a tablet which saturates with the undissolved dose.

    `kdis = k*dose50/(dose50 + D)` with the undissolved amount `D`.

    Attributes:
        k: the rate of the dissolution of a small dose (1/hr).
        dose50: the undissolved amount of half the rate (mg).
    """

    k: Q
    dose50: Q


type Formulation = FirstOrderDissolution | WeibullDissolution | DoseDependentDissolution
"""The dissolution of a tablet."""

ORAL = frozenset({"oral", "oral_solution"})
"""The routes into the intestine."""


@dataclass(frozen=True)
class Route:
    """The administration of a substance.

    Attributes:
        kind: intravenous, oral tablet (with dissolution) or oral solution.
        substance: the id of the administered substance.
        formulation: the dissolution of the tablet, the first order
            dissolution of the data (`Ka_dis`) if `None`.
    """

    kind: Literal["iv", "oral", "oral_solution"]
    substance: str
    formulation: Formulation | None = None

    def __post_init__(self) -> None:
        """Check the kind against its choices."""
        check_literals(self)

    @classmethod
    def iv(cls, s: str) -> "Route":
        """Get the intravenous route of a substance: bolus and infusion."""
        return cls("iv", s)

    @classmethod
    def oral(cls, s: str, formulation: Formulation | None = None) -> "Route":
        """Get the oral route of a substance as tablet, with a formulation."""
        return cls("oral", s, formulation)

    @classmethod
    def oral_solution(cls, s: str) -> "Route":
        """Get the oral route of a substance given as solution."""
        return cls("oral_solution", s)


_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


def _formulation(route: Route) -> Formulation:
    """The formulation of a tablet, the first order dissolution of the data by default."""
    if route.formulation is not None:
        return route.formulation
    return FirstOrderDissolution(k=value("routes", "oral", "Ka_dis"))


def dissolution_inputs(route: Route) -> list[tuple[str, UnitDefinition, str]]:
    """Get the parameters of the top model which the dissolution of a tablet reads.

    Args:
        route: an oral route (tablet).

    Returns:
        (id, unit, name) of every parameter of the formulation.
    """
    s = route.substance
    f = _formulation(route)
    if isinstance(f, WeibullDissolution):
        return [
            (ids.dissolution_t50(s), U.min, f"time of half dissolution {s} [min]"),
            (ids.dissolution_shape(s), U.dimensionless, f"dissolution shape {s}"),
            (ids.dissolution_lag(s), U.min, f"dissolution lag time {s} [min]"),
        ]
    pars = [(ids.dissolution(s), U.per_hr, f"dissolution rate {s} [1/hr]")]
    if isinstance(f, DoseDependentDissolution):
        pars.append(
            (
                ids.dissolution_dose50(s),
                U.mg,
                f"undissolved amount of half the dissolution rate {s} [mg]",
            )
        )
    if isinstance(f, FirstOrderDissolution):
        if f.fed is not None:
            pars.append(_fed_input(s))
        if f.release is not None:
            pars.append(_release_input(s))
    return pars


def _fed_input(s: str) -> tuple[str, UnitDefinition, str]:
    """The input of the factor of the fed state."""
    return ids.dissolution_fed_factor(s), U.dimensionless, f"dissolution factor fed {s}"


def _release_input(s: str) -> tuple[str, UnitDefinition, str]:
    """The input of the factor of the modified-release formulation."""
    return (
        ids.dissolution_release_factor(s),
        U.dimensionless,
        f"dissolution factor release {s}",
    )


def dissolution_reads_fed(route: Route) -> bool:
    """Check whether the dissolution of a tablet depends on the fed state."""
    f = _formulation(route)
    return isinstance(f, FirstOrderDissolution) and f.fed is not None


def is_weibull(route: Route) -> bool:
    """Check whether a tablet dissolves by a Weibull function (it needs the age)."""
    return isinstance(route.formulation, WeibullDissolution)


def dissolution_rate(route: Route, undissolved: str, age: str) -> tuple[str, bool]:
    """Get the dissolution rate of a tablet (1/min) on its undissolved amount.

    Args:
        route: an oral route (tablet).
        undissolved: the undissolved amount of all segments (mg), formula.
        age: the mean time since the dose of the undissolved drug (min), id.

    Returns:
        The formula and whether it reads a species (then the intestine makes
        it an assignment rule).
    """
    s = route.substance
    f = _formulation(route)
    k = f"{ids.dissolution(s)} / 60 min_per_hr"
    if isinstance(f, FirstOrderDissolution):
        if f.release is not None:
            k += f" * {ids.dissolution_release_factor(s)}"
        if f.fed is not None:
            fed = ids.dissolution_fed_factor(s)
            k += f" * (1 dimensionless + {ids.FED} * ({fed} - 1 dimensionless))"
        return k, False
    if isinstance(f, DoseDependentDissolution):
        d50 = ids.dissolution_dose50(s)
        return f"{k} * {d50} / ({d50} + {undissolved})", True
    t50, shape = ids.dissolution_t50(s), ids.dissolution_shape(s)
    lag = ids.dissolution_lag(s)
    # the offset keeps the hazard finite at the end of the lag time (shape < 1)
    x = f"(({age} - {lag}) / {t50} + 1e-9 dimensionless)"
    hazard = f"{math.log(2)!r} dimensionless * {shape} / {t50} * {x}^({shape} - 1 dimensionless)"
    return f"piecewise(0 per_min, {age} <= {lag}, {hazard})", True


def route_objects(route: Route, substance: Substance) -> list[Sbase]:
    """Get the dose state of a route and its parameters in the top model.

    The intravenous route applies its dose into the venous plasma; the oral
    routes are applied by the intestine.

    Args:
        route: the route.
        substance: the administered substance.

    Returns:
        The parameters (and for the intravenous route the rules and the
        reaction) of the route.
    """
    s = substance.sid
    if route.kind == "iv":
        return _iv(substance)
    if route.kind == "oral":
        objects: list[Sbase] = [
            _dose(
                ids.po_dose(s),
                f"oral dose {s} (tablet) [mg]",
                f"Oral dose of {substance.name} as tablet: the undissolved drug in "
                "the stomach, which dissolves and empties into the intestine.",
            )
        ]
        return objects + _formulation_parameters(route)
    return [
        _dose(
            ids.sol_dose(s),
            f"oral dose {s} (solution) [mg]",
            f"Oral dose of {substance.name} as solution, which enters the stomach "
            "lumen dissolved.",
        ),
        annotate(
            value_parameter(
                ids.solution(s),
                value("routes", "oral", "Ka_sol"),
                U.per_hr,
                f"application rate {s} (solution) [1/hr]",
                _SBO_PARAMETER,
                "Swallowing of the solution into the stomach lumen.",
            ),
            "application_rate",
        ),
    ]


def _formulation_parameters(route: Route) -> list[Sbase]:
    """The parameters of the dissolution of a tablet."""
    return [annotate(p, "dissolution") for p in _dissolution_parameters(route)]


def _dissolution_parameters(route: Route) -> list[Parameter]:
    """The parameters of the dissolution of a tablet, not annotated."""
    s = route.substance
    f = _formulation(route)
    if isinstance(f, WeibullDissolution):
        return [
            value_parameter(
                ids.dissolution_t50(s),
                f.t50,
                U.min,
                f"time of half dissolution {s} [min]",
                _SBO_PARAMETER,
                "Weibull dissolution: the time of half dissolution after the lag time.",
            ),
            value_parameter(
                ids.dissolution_shape(s),
                f.shape,
                U.dimensionless,
                f"dissolution shape {s}",
                _SBO_PARAMETER,
                "Weibull dissolution: the shape (1 is a first order dissolution).",
            ),
            value_parameter(
                ids.dissolution_lag(s),
                f.lag,
                U.min,
                f"dissolution lag time {s} [min]",
                _SBO_PARAMETER,
                "Weibull dissolution: the time after the dose before the tablet "
                "dissolves.",
            ),
        ]
    pars: list[Parameter] = [
        value_parameter(
            ids.dissolution(s),
            f.k,
            U.per_hr,
            f"dissolution rate {s} [1/hr]",
            _SBO_PARAMETER,
            "First order dissolution of the undissolved drug in every segment of "
            "the intestine."
            if isinstance(f, FirstOrderDissolution)
            else "Dissolution rate of a small undissolved amount; it saturates "
            "with the undissolved amount.",
        )
    ]
    if isinstance(f, DoseDependentDissolution):
        pars.append(
            value_parameter(
                ids.dissolution_dose50(s),
                f.dose50,
                U.mg,
                f"undissolved amount of half the dissolution rate {s} [mg]",
                _SBO_PARAMETER,
                "The dissolution rate is `Ka_dis*D50/(D50 + D)` with the "
                "undissolved amount `D`.",
            )
        )
    if isinstance(f, FirstOrderDissolution):
        if f.fed is not None:
            pars.append(
                value_parameter(
                    ids.dissolution_fed_factor(s),
                    f.fed,
                    U.dimensionless,
                    f"dissolution factor fed {s}",
                    _SBO_PARAMETER,
                    "Factor of the dissolution rate in the fed state: the rate is "
                    "`Ka_dis*(1 + FED*(f_dis_fed - 1))`.",
                )
            )
        if f.release is not None:
            pars.append(
                value_parameter(
                    ids.dissolution_release_factor(s),
                    f.release,
                    U.dimensionless,
                    f"dissolution factor release {s}",
                    _SBO_PARAMETER,
                    "Factor of the dissolution rate of a modified-release "
                    "formulation (1 is the immediate-release tablet); the "
                    "experiments of such a formulation set it.",
                )
            )
    return pars


def _dose(pid: str, name: str, text: str) -> Parameter:
    """The dose state of a route (mg), set or added to by the experiments."""
    return annotate(
        Parameter(
            pid,
            0.0,
            U.mg,
            constant=False,
            name=name,
            sboTerm=_SBO_PARAMETER,
            notes=f"{text} An experiment sets the dose or adds to the rest of an "
            "earlier one (`add_doses`).",
        ),
        "dose",
    )


def _iv(s: Substance) -> list[Sbase]:
    """The intravenous bolus and infusion into the venous plasma."""
    sid = s.sid
    dose, ti = ids.iv_dose(sid), ids.injection_time(sid)
    kiv, ri = ids.injection_rate(sid), ids.infusion_rate(sid)
    reaction = ids.dose_reaction("iv", sid)
    kiv_formula = f"ln(2) / {ti} * 60 s_per_min"
    return [
        _dose(
            dose,
            f"intravenous dose {sid} [mg]",
            f"Intravenous dose of {s.name}, applied into the venous plasma.",
        ),
        annotate(
            value_parameter(
                ti,
                value("routes", "iv", "ti"),
                U.second,
                f"injection time {sid} [s]",
                _SBO_PARAMETER,
                "Half-life of the application of the intravenous dose.",
            ),
            "iv_time",
        ),
        annotate(
            Parameter(
                kiv,
                None,
                U.per_min,
                constant=False,
                name=f"injection rate {sid} [1/min]",
                sboTerm=_SBO_PARAMETER,
                notes="Rate of the application of the intravenous dose; a rule, so "
                f"that a dose protocol may change `{ti}` between doses.\n\n"
                f"`{kiv} = {kiv_formula}`",
            ),
            "iv_rate",
        ),
        AssignmentRule(kiv, kiv_formula, U.per_min),
        annotate(
            Parameter(
                ri,
                0.0,
                U.mg_per_min,
                constant=True,
                name=f"infusion rate {sid} [mg/min]",
                sboTerm=_SBO_PARAMETER,
                notes="Rate of a constant intravenous infusion, which adds to the "
                "intravenous dose state; an experiment sets it.",
            ),
            "infusion_rate",
        ),
        annotate(
            Reaction(
                sid=reaction,
                equation=f"-> {ids.pool_species(Pool.VENOUS, sid)}",
                formula=(f"{kiv} * {dose} / {ids.mr(sid)}", U.mmole_per_min),
                name=f"intravenous application {sid}",
                sboTerm=SBO.TRANSPORT_REACTION,
            ),
            "iv_application",
            pool=Pool.VENOUS,
        ),
        RateRule(
            dose,
            f"-{reaction} * {ids.mr(sid)} + {ri}",
            U.mg_per_min,
            name=f"intravenous dose {sid}",
        ),
    ]
