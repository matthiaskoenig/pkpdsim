"""The routes of administration and their dose states.

A dose is a dose state in mg, which an experiment sets or adds to
(`PKPDSimulationExperiment.add_doses`: `"PODOSE_<s> + 10.0"`), and which is
applied by a first order process, the conventions of the hand-written models:

- intravenous: `IVDOSE_<s>` drains into the venous plasma `Cve_<s>` with the
  rate `Kiv_<s> = ln(2)/ti_<s>` (`ti_<s>` the half-life of the injection in s);
  a constant infusion `Ri_<s>` (mg/min) adds to the dose state;
- oral tablet: `PODOSE_<s>` dissolves into the intestinal lumen `Clumen_<s>`
  with the first order dissolution `Ka_dis_<s>` (1/hr);
- oral solution: `SOLDOSE_<s>` enters the intestinal lumen with the fast rate
  `Ka_sol_<s>` (1/hr).

The dose in mg is converted to mmole with the molecular weight `Mr_<s>`.
Until the segmented intestine, the lumen is the single lumen of the builder.
"""

from dataclasses import dataclass
from typing import Literal

from sbmlutils.factory import (
    AssignmentRule,
    Parameter,
    RateRule,
    Reaction,
    Sbase,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import value
from pkpdsim.pbpk.organs import Pool
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q

__all__ = ["FirstOrderDissolution", "Route", "route_objects"]


@dataclass(frozen=True)
class FirstOrderDissolution:
    """The first order dissolution of a tablet.

    Attributes:
        k: the rate of the dissolution (1/hr).
    """

    k: Q


@dataclass(frozen=True)
class Route:
    """The administration of a substance.

    Attributes:
        kind: intravenous, oral tablet (with dissolution) or oral solution.
        substance: the id of the administered substance.
        formulation: the formulation of the tablet, the default first order
            dissolution if `None`.
    """

    kind: Literal["iv", "oral", "oral_solution"]
    substance: str
    formulation: FirstOrderDissolution | None = None

    @classmethod
    def iv(cls, s: str) -> "Route":
        """Get the intravenous route of a substance: bolus and infusion."""
        return cls("iv", s)

    @classmethod
    def oral(cls, s: str, formulation: FirstOrderDissolution | None = None) -> "Route":
        """Get the oral route of a substance as tablet, with a formulation."""
        return cls("oral", s, formulation)

    @classmethod
    def oral_solution(cls, s: str) -> "Route":
        """Get the oral route of a substance given as solution."""
        return cls("oral_solution", s)


_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


def route_objects(route: Route, substance: Substance) -> list[Sbase]:
    """Get the dose state of a route with its application into the body.

    Args:
        route: the route.
        substance: the administered substance.

    Returns:
        The parameters, rules, reaction and rate rule of the route.
    """
    s = substance.sid
    if route.kind == "iv":
        return _iv(substance)
    if route.kind == "oral":
        rate = value("routes", "oral", "Ka_dis")
        if route.formulation is not None:
            rate = route.formulation.k
        return _oral(
            substance,
            "po",
            dose=ids.po_dose(s),
            rate=(ids.dissolution(s), rate),
            what="tablet",
            text="First order dissolution of the tablet, swallowing and transfer "
            "into the intestinal lumen.",
        )
    return _oral(
        substance,
        "sol",
        dose=ids.sol_dose(s),
        rate=(ids.solution(s), value("routes", "oral", "Ka_sol")),
        what="solution",
        text="Swallowing of the solution and transfer into the intestinal lumen.",
    )


def _dose(pid: str, name: str, text: str) -> Parameter:
    """The dose state of a route (mg), set or added to by the experiments."""
    return Parameter(
        pid,
        0.0,
        U.mg,
        constant=False,
        name=name,
        sboTerm=_SBO_PARAMETER,
        notes=f"{text} An experiment sets the dose or adds to the rest of an "
        "earlier one (`add_doses`).",
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
        value_parameter(
            ti,
            value("routes", "iv", "ti"),
            U.second,
            f"injection time {sid} [s]",
            _SBO_PARAMETER,
            "Half-life of the application of the intravenous dose.",
        ),
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
        AssignmentRule(kiv, kiv_formula, U.per_min),
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
        Reaction(
            sid=reaction,
            equation=f"-> {ids.pool_species(Pool.VENOUS, sid)}",
            formula=(f"{kiv} * {dose} / {ids.mr(sid)}", U.mmole_per_min),
            name=f"intravenous application {sid}",
            sboTerm=SBO.TRANSPORT_REACTION,
        ),
        RateRule(
            dose,
            f"-{reaction} * {ids.mr(sid)} + {ri}",
            U.mg_per_min,
            name=f"intravenous dose {sid}",
        ),
    ]


def _oral(
    s: Substance,
    kind: str,
    dose: str,
    rate: tuple[str, Q],
    what: str,
    text: str,
) -> list[Sbase]:
    """An oral dose state which enters the intestinal lumen with first order."""
    sid = s.sid
    reaction = ids.dose_reaction(kind, sid)
    rid, q = rate
    return [
        _dose(
            dose,
            f"oral dose {sid} ({what}) [mg]",
            f"Oral dose of {s.name} as {what}, which enters the intestinal lumen.",
        ),
        value_parameter(
            rid,
            q,
            U.per_hr,
            f"application rate {sid} ({what}) [1/hr]",
            _SBO_PARAMETER,
            text,
        ),
        Reaction(
            sid=reaction,
            equation=f"-> {ids.lumen(sid)}",
            formula=(
                f"{rid} / 60 min_per_hr * {dose} / {ids.mr(sid)}",
                U.mmole_per_min,
            ),
            name=f"oral application {sid} ({what})",
            sboTerm=SBO.TRANSPORT_REACTION,
        ),
        RateRule(
            dose,
            f"-{reaction} * {ids.mr(sid)}",
            U.mg_per_min,
            name=f"oral dose {sid} ({what})",
        ),
    ]
