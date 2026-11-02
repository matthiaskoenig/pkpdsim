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
  its own, it dissolves there and leaves undissolved into the feces; an
  enteric-coated tablet (`enteric=True`) does not dissolve in the stomach,
  it empties undissolved with the gastric emptying of the condition and
  dissolves from the duodenum on;
- oral solution: `SOLDOSE_<s>` enters the stomach lumen dissolved with the
  fast rate `Ka_sol_<s>` (1/hr);
- depot routes, intramuscular, subcutaneous and rectal: `IMDOSE_<s>`,
  `SCDOSE_<s>`, `RECDOSE_<s>` drain with the half-life `ti_<k>_<s>` (s) of
  the application into the depot of the route, an amount of the top model
  (`Aim_<s>`, `Asc_<s>`, `Arec_<s>`, mmole, in the compartment `Vim`,
  `Vsc`, `Vrec`), which is absorbed with first order `ka_<k>_<s>` (1/hr):
  the fraction `f_<k>_<s>` into the venous plasma `Cve_<s>`, the rest into
  a sink (`Aloss_im_<s>`, `Aloss_sc_<s>`, degraded at the site of the
  injection; the feces `Afeces_<s>` from the rectum); of the absorbed
  rectal dose the fraction `f_rec_portal_<s>` enters the portal plasma
  `Cpo_<s>` and passes the liver (the superior rectal vein), the rest the
  venous plasma (the middle and inferior rectal veins).

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
- `DoseDependentDissolution(k, dose50, fed)`: `kdis = k*D50/(D50 + D)` with
  the undissolved amount `D` (mg) of all segments, a dissolution which
  saturates with the undissolved dose (solubility-limited), at most
  `k*D50`, times `(1 + FED*(f_dis_fed_<s> - 1))` of the fed state if set;
  of an enteric-coated tablet `D` leaves out the coated tablet in the
  stomach (`PODOSE_<s>`), which is not exposed to the intestinal fluid.
"""

import math
from dataclasses import dataclass
from typing import Literal, NamedTuple

from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    Parameter,
    RateRule,
    Reaction,
    Sbase,
    Species,
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
    "DEPOT",
    "ORAL",
    "Depot",
    "DoseDependentDissolution",
    "FirstOrderDissolution",
    "Formulation",
    "Route",
    "WeibullDissolution",
    "depot_site",
    "dissolution_inputs",
    "dissolution_rate",
    "is_enteric",
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
        enteric: an enteric-coated tablet, which does not dissolve in the
            stomach: it empties undissolved with the gastric emptying of the
            condition and dissolves from the duodenum on. The single lumen
            has no stomach and raises.
    """

    k: Q
    fed: Q | None = None
    release: Q | None = None
    enteric: bool = False


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
        fed: the factor of the dissolution rate in the fed state (`FED = 1`),
            the body parameter `f_dis_fed_<s>`; the rate is
            `k*(1 + FED*(fed - 1))*dose50/(dose50 + D)`. No parameter if
            `None`. The single lumen has no fed state and raises.
        enteric: an enteric-coated tablet, which does not dissolve in the
            stomach (see `FirstOrderDissolution`); `D` is the undissolved
            amount outside the stomach, without the coated tablet.
    """

    k: Q
    dose50: Q
    fed: Q | None = None
    enteric: bool = False


type Formulation = FirstOrderDissolution | WeibullDissolution | DoseDependentDissolution
"""The dissolution of a tablet."""

ORAL = frozenset({"oral", "oral_solution"})
"""The routes into the intestine."""

DEPOT = frozenset({"im", "sc", "rectal"})
"""The routes into a depot at the site of the administration."""


@dataclass(frozen=True)
class Depot:
    """The parameters of a depot route which a model gives instead of the defaults.

    Every value is a number, a `Value` with its source or a fitted value
    (`fitted(FITS.parameter(...))`); `None` takes the default of the data
    (`data/routes.toml`).

    Attributes:
        ka: the rate of the absorption from the depot (1/hr), `ka_<k>_<s>`.
        fraction: the fraction of the depot which is absorbed, `f_<k>_<s>`
            (1 by default), the rest leaves into the sink.
        portal: the fraction of the absorbed dose which enters the portal
            vein, `f_rec_portal_<s>`; the rectal route only.
        application: the half-life of the application of the dose into the
            depot (s), `ti_<k>_<s>`; a formulation which releases slowly (a
            suppository) sets it.
    """

    ka: Q | None = None
    fraction: Q | None = None
    portal: Q | None = None
    application: Q | None = None


@dataclass(frozen=True)
class Route:
    """The administration of a substance.

    Attributes:
        kind: intravenous, oral tablet (with dissolution), oral solution,
            intramuscular, subcutaneous or rectal.
        substance: the id of the administered substance.
        formulation: the dissolution of the tablet, the first order
            dissolution of the data (`Ka_dis`) if `None`.
        depot: the parameters of a depot route (intramuscular, subcutaneous,
            rectal), the defaults of the data if `None`.
    """

    kind: Literal["iv", "oral", "oral_solution", "im", "sc", "rectal"]
    substance: str
    formulation: Formulation | None = None
    depot: Depot | None = None

    def __post_init__(self) -> None:
        """Check the kind against its choices and the depot against the kind.

        Raises:
            ValueError: if the kind is none of its choices, a route which is
                no depot route has a depot, or a route which is not rectal a
                portal fraction.
        """
        check_literals(self)
        if self.depot is None:
            return
        if self.kind not in DEPOT:
            raise ValueError(
                f"Route '{self.kind}' of '{self.substance}': only the routes "
                f"{sorted(DEPOT)} have a depot."
            )
        if self.depot.portal is not None and self.kind != "rectal":
            raise ValueError(
                f"Route '{self.kind}' of '{self.substance}': only the rectal "
                "route has a portal fraction."
            )

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

    @classmethod
    def im(
        cls,
        s: str,
        *,
        ka: Q | None = None,
        fraction: Q | None = None,
        application: Q | None = None,
    ) -> "Route":
        """Get the intramuscular route of a substance, see `Depot`."""
        return cls("im", s, depot=_given(ka, fraction, None, application))

    @classmethod
    def sc(
        cls,
        s: str,
        *,
        ka: Q | None = None,
        fraction: Q | None = None,
        application: Q | None = None,
    ) -> "Route":
        """Get the subcutaneous route of a substance, see `Depot`."""
        return cls("sc", s, depot=_given(ka, fraction, None, application))

    @classmethod
    def rectal(
        cls,
        s: str,
        *,
        ka: Q | None = None,
        fraction: Q | None = None,
        portal: Q | None = None,
        application: Q | None = None,
    ) -> "Route":
        """Get the rectal route of a substance, see `Depot`."""
        return cls("rectal", s, depot=_given(ka, fraction, portal, application))


def _given(
    ka: Q | None, fraction: Q | None, portal: Q | None, application: Q | None
) -> Depot | None:
    """The depot of the given values, `None` if none is given (the route of the defaults)."""
    values = (ka, fraction, portal, application)
    return None if all(v is None for v in values) else Depot(*values)


_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


def _formulation(route: Route) -> Formulation:
    """The formulation of a tablet, the first order dissolution of the data by default."""
    if route.formulation is not None:
        return route.formulation
    return FirstOrderDissolution(k=value("routes", "oral", "Ka_dis"))


class _DissolutionParameter(NamedTuple):
    """A parameter of the dissolution of a tablet: the one place of its name."""

    pid: str
    value: Q
    unit: UnitDefinition
    name: str
    notes: str


def _dissolution_table(route: Route) -> list[_DissolutionParameter]:
    """The parameters of the dissolution of a tablet with their values and texts."""
    s = route.substance
    f = _formulation(route)
    P = _DissolutionParameter
    if isinstance(f, WeibullDissolution):
        return [
            P(
                ids.dissolution_t50(s),
                f.t50,
                U.min,
                f"time of half dissolution {s} [min]",
                "Weibull dissolution: the time of half dissolution after the lag time.",
            ),
            P(
                ids.dissolution_shape(s),
                f.shape,
                U.dimensionless,
                f"dissolution shape {s}",
                "Weibull dissolution: the shape (1 is a first order dissolution).",
            ),
            P(
                ids.dissolution_lag(s),
                f.lag,
                U.min,
                f"dissolution lag time {s} [min]",
                "Weibull dissolution: the time after the dose before the tablet "
                "dissolves.",
            ),
        ]
    table = [
        P(
            ids.dissolution(s),
            f.k,
            U.per_hr,
            f"dissolution rate {s} [1/hr]",
            "First order dissolution of the undissolved drug in every segment of "
            "the intestine."
            if isinstance(f, FirstOrderDissolution)
            else "Dissolution rate of a small undissolved amount; it saturates "
            "with the undissolved amount.",
        )
    ]
    if isinstance(f, DoseDependentDissolution):
        table.append(
            P(
                ids.dissolution_dose50(s),
                f.dose50,
                U.mg,
                f"undissolved amount of half the dissolution rate {s} [mg]",
                "The dissolution rate is `Ka_dis*D50/(D50 + D)` with the "
                "undissolved amount `D`.",
            )
        )
    if f.fed is not None:
        table.append(
            P(
                ids.dissolution_fed_factor(s),
                f.fed,
                U.dimensionless,
                f"dissolution factor fed {s}",
                "Factor of the dissolution rate in the fed state: the rate is "
                "`Ka_dis*(1 + FED*(f_dis_fed - 1))`"
                + (
                    " times the saturation `D50/(D50 + D)`."
                    if isinstance(f, DoseDependentDissolution)
                    else "."
                ),
            )
        )
    if isinstance(f, FirstOrderDissolution) and f.release is not None:
        table.append(
            P(
                ids.dissolution_release_factor(s),
                f.release,
                U.dimensionless,
                f"dissolution factor release {s}",
                "Factor of the dissolution rate of a modified-release "
                "formulation (1 is the immediate-release tablet); the "
                "experiments of such a formulation set it.",
            )
        )
    return table


def dissolution_inputs(route: Route) -> list[tuple[str, UnitDefinition, str]]:
    """Get the parameters of the top model which the dissolution of a tablet reads.

    Args:
        route: an oral route (tablet).

    Returns:
        (id, unit, name) of every parameter of the formulation.
    """
    return [(p.pid, p.unit, p.name) for p in _dissolution_table(route)]


def dissolution_reads_fed(route: Route) -> bool:
    """Check whether the dissolution of a tablet depends on the fed state."""
    f = _formulation(route)
    return not isinstance(f, WeibullDissolution) and f.fed is not None


def is_enteric(route: Route) -> bool:
    """Check whether a tablet is enteric-coated (it does not dissolve in the stomach)."""
    f = _formulation(route)
    return not isinstance(f, WeibullDissolution) and f.enteric


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
    if isinstance(f, FirstOrderDissolution) and f.release is not None:
        k += f" * {ids.dissolution_release_factor(s)}"
    if not isinstance(f, WeibullDissolution) and f.fed is not None:
        fed = ids.dissolution_fed_factor(s)
        k += f" * (1 dimensionless + {ids.FED} * ({fed} - 1 dimensionless))"
    if isinstance(f, FirstOrderDissolution):
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

    The intravenous route applies its dose into the venous plasma, a depot
    route into its depot, from which it is absorbed into the venous (and the
    portal) plasma; the oral routes are applied by the intestine. The
    compartment of the depots of a route is `depot_site`.

    Args:
        route: the route.
        substance: the administered substance.

    Returns:
        The parameters (and for the intravenous and the depot routes the
        species, rules and reactions) of the route.
    """
    s = substance.sid
    if route.kind == "iv":
        return _iv(substance)
    if route.kind in DEPOT:
        return _depot(route, substance)
    if route.kind == "oral":
        fate = (
            "which empties undissolved into the intestine (enteric-coated, it "
            "does not dissolve in the stomach)."
            if is_enteric(route)
            else "which dissolves and empties into the intestine."
        )
        objects: list[Sbase] = [
            _dose(
                ids.po_dose(s),
                f"oral dose {s} (tablet) [mg]",
                f"Oral dose of {substance.name} as tablet: the undissolved drug in "
                f"the stomach, {fate}",
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
    return [
        value_parameter(p.pid, p.value, p.unit, p.name, _SBO_PARAMETER, p.notes)
        for p in _dissolution_table(route)
    ]


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


class _Site(NamedTuple):
    """The words of a depot route: its name and the site of its depot."""

    name: str
    site: str


_SITES: dict[str, _Site] = {
    "im": _Site("intramuscular", "the muscle at the site of the injection"),
    "sc": _Site("subcutaneous", "the subcutis at the site of the injection"),
    "rectal": _Site("rectal", "the rectum"),
}
"""The depot routes by kind, the sections of `data/routes.toml` of their defaults."""


def _depot_kind(route: Route) -> Literal["im", "sc", "rectal"]:
    """The kind of a depot route, narrowed."""
    kind = route.kind
    if kind == "im" or kind == "sc" or kind == "rectal":
        return kind
    raise ValueError(f"Route '{kind}' is no depot route.")


def depot_site(route: Route) -> Compartment:
    """Get the compartment of the depots of the kind of a route, once per model.

    Args:
        route: a depot route.

    Returns:
        The compartment of amounts with the nominal volume 1 l, `Vim`,
        `Vsc`, `Vrec`.

    Raises:
        ValueError: if the route is no depot route.
    """
    kind = _depot_kind(route)
    site = _SITES[kind]
    return annotate(
        Compartment(
            ids.depot_compartment(kind),
            1.0,
            U.litre,
            constant=True,
            name=f"{site.name} depot",
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            notes=f"The depot of the {site.name} doses in {site.site}. Nominal "
            "volume of a compartment of amounts, which enters no rate.",
        ),
        f"depot_site_{kind}",
    )


def _depot(route: Route, s: Substance) -> list[Sbase]:
    """The dose state, the depot and the absorption of a depot route."""
    kind = _depot_kind(route)
    site = _SITES[kind]
    given = route.depot if route.depot is not None else Depot()
    sid, mr = s.sid, ids.mr(s.sid)
    dose, depot = ids.depot_dose(kind, sid), ids.depot(kind, sid)
    ti = ids.depot_application_time(kind, sid)
    ka, f = ids.depot_absorption_rate(kind, sid), ids.depot_fraction(kind, sid)
    application = ids.depot_application(kind, sid)
    loss, sink = ids.depot_loss_reaction(kind, sid), ids.depot_loss(kind, sid)

    def default(key: str) -> Q:
        return value("routes", kind, key)

    rate = f"{ka} / 60 min_per_hr * {depot}"
    objects: list[Sbase] = [
        _dose(
            dose,
            f"{site.name} dose {sid} [mg]",
            f"{site.name.capitalize()} dose of {s.name}, applied into its depot "
            f"`{depot}` in {site.site}.",
        ),
        annotate(
            value_parameter(
                ti,
                given.application
                if given.application is not None
                else value("routes", "depot", "ti"),
                U.second,
                f"application time {sid} ({site.name}) [s]",
                _SBO_PARAMETER,
                f"Half-life of the application of the {site.name} dose into its depot.",
            ),
            f"depot_time_{kind}",
        ),
        annotate(
            value_parameter(
                ka,
                given.ka if given.ka is not None else default("ka"),
                U.per_hr,
                f"absorption rate {sid} ({site.name}) [1/hr]",
                _SBO_PARAMETER,
                f"First order absorption from the {site.name} depot.",
            ),
            f"depot_rate_{kind}",
        ),
        annotate(
            value_parameter(
                f,
                given.fraction if given.fraction is not None else default("f"),
                U.dimensionless,
                f"fraction absorbed {sid} ({site.name})",
                _SBO_PARAMETER,
                f"Fraction of the {site.name} depot which is absorbed; the rest "
                f"leaves into `{sink}`.",
            ),
            f"depot_fraction_{kind}",
        ),
        annotate(
            Species(
                depot,
                ids.depot_compartment(kind),
                initialAmount=0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=True,
                name=f"{s.name} ({site.name} depot)",
                sboTerm=SBO.SIMPLE_CHEMICAL,
            ),
            f"depot_amount_{kind}",
            substance=s,
        ),
        annotate(
            Reaction(
                sid=application,
                equation=f"-> {depot}",
                formula=(
                    f"ln(2) / {ti} * 60 s_per_min * {dose} / {mr}",
                    U.mmole_per_min,
                ),
                name=f"{site.name} application {sid}",
                sboTerm=SBO.TRANSPORT_REACTION,
            ),
            f"depot_application_{kind}",
        ),
        RateRule(
            dose,
            f"-{application} * {mr}",
            U.mg_per_min,
            name=f"{site.name} dose {sid}",
        ),
    ]
    if kind != "rectal":
        objects.append(
            annotate(
                Species(
                    sink,
                    ids.depot_compartment(kind),
                    initialAmount=0.0,
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=True,
                    name=f"{s.name} ({site.name} depot, not absorbed)",
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                ),
                f"depot_loss_amount_{kind}",
                substance=s,
            )
        )
    venous = f"{f} * {rate}"
    if kind == "rectal":
        portal = ids.rectal_portal_fraction(sid)
        objects.append(
            annotate(
                value_parameter(
                    portal,
                    given.portal if given.portal is not None else default("f_portal"),
                    U.dimensionless,
                    f"portal fraction {sid} (rectal)",
                    _SBO_PARAMETER,
                    "Fraction of the absorbed rectal dose which enters the portal "
                    "vein (superior rectal vein) and passes the liver; the rest "
                    "enters the venous blood (middle and inferior rectal veins).",
                ),
                "rectal_portal_fraction",
            )
        )
        venous = f"{f} * (1 dimensionless - {portal}) * {rate}"
        objects.append(
            annotate(
                Reaction(
                    sid=ids.depot_absorption(kind, sid, portal=True),
                    equation=f"{depot} -> {ids.pool_species(Pool.PORTAL, sid)}",
                    formula=(f"{f} * {portal} * {rate}", U.mmole_per_min),
                    name=f"{site.name} absorption {sid} (portal)",
                    sboTerm=SBO.TRANSPORT_REACTION,
                ),
                f"depot_absorption_{kind}",
                pool=Pool.PORTAL,
            )
        )
    objects += [
        annotate(
            Reaction(
                sid=ids.depot_absorption(kind, sid),
                equation=f"{depot} -> {ids.pool_species(Pool.VENOUS, sid)}",
                formula=(venous, U.mmole_per_min),
                name=f"{site.name} absorption {sid}",
                sboTerm=SBO.TRANSPORT_REACTION,
            ),
            f"depot_absorption_{kind}",
            pool=Pool.VENOUS,
        ),
        annotate(
            Reaction(
                sid=loss,
                equation=f"{depot} -> {sink}",
                formula=(f"(1 dimensionless - {f}) * {rate}", U.mmole_per_min),
                name=f"{site.name} depot not absorbed {sid}",
                sboTerm=SBO.TRANSPORT_REACTION,
            ),
            f"depot_loss_{kind}",
        ),
    ]
    return objects
