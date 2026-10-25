"""The sinks of the body and the single intestinal lumen.

Every substance has the amounts in urine `Aurine_<s>` and feces `Afeces_<s>`
(compartments of amounts, `Vurine` and `Vfeces`), the bile `Cbile_<s>` in
`Vbile`, which drains into the lumen with first order `k_bile`, and the
intestinal lumen `Clumen_<s>` in `Vlumen`.

Until the segmented intestine the lumen is one compartment, the fluid of the
small intestine: the oral doses enter it, a substance is absorbed from it
with first order `ka_<s>` into the portal vein (the gut wall is passed
without distribution; a top reaction cannot reach the plasma of the small
intestine without renaming its id) and leaves it with the small intestinal
transit `k_feces = 1/t_lumen` into the feces. `ka_<s>` is the default of
the data for a substance with an oral route and 0 for any other, which is
then not reabsorbed after biliary excretion unless a model sets it.
"""

from sbmlutils.factory import (
    Compartment,
    InitialAssignment,
    Parameter,
    Reaction,
    Sbase,
    Species,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import cite, value
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.organs import Pool
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value

__all__ = ["sink_objects"]

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_NOMINAL = "Nominal volume of a compartment of amounts, which enters no rate."
_ORAL = {"oral", "oral_solution"}
"""The routes into the intestinal lumen."""


def sink_objects(definition: PBPKDefinition) -> list[Sbase]:
    """Get the sinks and the single intestinal lumen of a definition.

    Args:
        definition: the definition.

    Returns:
        The compartments, parameters, species and transits of the top model.
    """
    objects: list[Sbase] = [
        _nominal(ids.VURINE, "urine"),
        _nominal(ids.VFECES, "feces"),
        _value_compartment(
            ids.VBILE,
            value("intestine_human", "bile", "volume"),
            "bile",
            "Bile, which drains into the intestinal lumen.",
        ),
        _value_compartment(
            ids.VLUMEN,
            value("intestine_human", "lumen", "volume"),
            "intestinal lumen",
            "The single intestinal lumen: the fluid of the small intestine, into "
            "which the oral doses and the bile enter.",
        ),
        value_parameter(
            ids.K_BILE,
            value("intestine_human", "bile", "k_bile"),
            U.per_min,
            "transit bile to intestinal lumen [1/min]",
            _SBO_PARAMETER,
        ),
        value_parameter(
            ids.T_LUMEN,
            value("intestine_human", "lumen", "transit_time"),
            U.min,
            "transit time intestinal lumen [min]",
            _SBO_PARAMETER,
            "Transit time of the single lumen into the feces.",
        ),
        Parameter(
            ids.K_FECES,
            None,
            U.per_min,
            constant=True,
            name="transit intestinal lumen to feces [1/min]",
            sboTerm=_SBO_PARAMETER,
            notes=f"The reciprocal of the transit time.\n\n"
            f"`{ids.K_FECES} = 1 dimensionless / {ids.T_LUMEN}`",
        ),
        InitialAssignment(ids.K_FECES, f"1 dimensionless / {ids.T_LUMEN}", U.per_min),
    ]
    oral = {r.substance for r in definition.routes if r.kind in _ORAL}
    ka = value("intestine_human", "lumen", "ka")
    for s in definition.substances:
        sid = s.sid
        lumen, bile = ids.lumen(sid), ids.bile(sid)
        objects += [
            _species(ids.urine(sid), ids.VURINE, f"{s.name} (urine)", amount=True),
            _species(ids.feces(sid), ids.VFECES, f"{s.name} (feces)", amount=True),
            _species(bile, ids.VBILE, f"{s.name} (bile)"),
            _species(lumen, ids.VLUMEN, f"{s.name} (intestinal lumen)"),
            value_parameter(
                ids.absorption(sid),
                ka if sid in oral else 0.0,
                U.per_hr,
                f"absorption {sid} [1/hr]",
                _SBO_PARAMETER,
                "First order absorption from the intestinal lumen into the portal "
                "vein; a model fits it."
                if sid in oral
                else "No absorption: the substance has no oral route, it is not "
                "reabsorbed after biliary excretion unless a model sets it.",
            ),
            _transit(
                ids.bile_reaction(sid),
                f"{bile} -> {lumen}",
                f"{ids.K_BILE} * {ids.VBILE} * {bile}",
                f"transit {sid} bile to intestinal lumen",
            ),
            _transit(
                ids.absorption_reaction(sid),
                f"{lumen} -> {ids.pool_species(Pool.PORTAL, sid)}",
                f"{ids.absorption(sid)} / 60 min_per_hr * {ids.VLUMEN} * {lumen}",
                f"absorption {sid} intestinal lumen to portal vein",
            ),
            _transit(
                ids.feces_reaction(sid),
                f"{lumen} -> {ids.feces(sid)}",
                f"{ids.K_FECES} * {ids.VLUMEN} * {lumen}",
                f"transit {sid} intestinal lumen to feces",
            ),
        ]
    return objects


def _nominal(cid: str, name: str) -> Compartment:
    """A compartment of amounts with the nominal volume 1 l."""
    return Compartment(
        cid,
        1.0,
        U.litre,
        constant=True,
        name=name,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        notes=_NOMINAL,
    )


def _value_compartment(cid: str, v: Value, name: str, text: str) -> Compartment:
    """A constant compartment of a value of the data, its source in the notes."""
    return Compartment(
        cid,
        v.value,
        U.litre,
        constant=True,
        name=name,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        notes=f"{text} {v.notes}\n\nSource:\n\n- {cite(str(v.source))}",
    )


def _species(sid: str, compartment: str, name: str, amount: bool = False) -> Species:
    """A sink species, initially zero."""
    return Species(
        sid,
        compartment,
        initialAmount=0.0 if amount else None,
        initialConcentration=None if amount else 0.0,
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=amount,
        name=name,
        sboTerm=SBO.SIMPLE_CHEMICAL,
    )


def _transit(sid: str, equation: str, rate: str, name: str) -> Reaction:
    """A transit of the top model, the rate in mmole/min."""
    return Reaction(
        sid=sid,
        equation=equation,
        formula=(rate, U.mmole_per_min),
        name=name,
        sboTerm=SBO.TRANSPORT_REACTION,
    )
