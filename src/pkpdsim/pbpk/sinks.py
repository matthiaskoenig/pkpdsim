"""The sinks of the body: urine, feces and bile.

Every substance has the amounts in urine `Aurine_<s>` and feces `Afeces_<s>`
(compartments of amounts, `Vurine` and `Vfeces`) and the bile `Cbile_<s>` in
`Vbile`, which drains with first order `k_bile` into the duodenum of the
intestine (`pkpdsim.pbpk.variants.intestine`, which owns the lumen and that
transit); with the standard organs as intestine, which have no lumen, the
bile drains into the feces.
"""

from sbmlutils.factory import (
    Compartment,
    Reaction,
    Sbase,
    Species,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import cite, value
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value

__all__ = ["sink_objects"]

_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_NOMINAL = "Nominal volume of a compartment of amounts, which enters no rate."


def sink_objects(definition: PBPKDefinition, lumen: bool = True) -> list[Sbase]:
    """Get the sinks of a definition.

    Args:
        definition: the definition.
        lumen: the intestine has a lumen, which takes the bile; without one
            the bile drains into the feces here.

    Returns:
        The compartments, parameters and species of the top model, and the
        transit of the bile into the feces without lumen.
    """
    objects: list[Sbase] = [
        _nominal(ids.VURINE, "urine"),
        _nominal(ids.VFECES, "feces"),
        _value_compartment(
            ids.VBILE,
            value("intestine_human", "bile", "volume"),
            "bile",
            "Bile, which drains into the duodenum.",
        ),
        value_parameter(
            ids.K_BILE,
            value("intestine_human", "bile", "k_bile"),
            U.per_min,
            "transit bile to intestinal lumen [1/min]",
            _SBO_PARAMETER,
        ),
    ]
    for s in definition.substances:
        sid = s.sid
        bile = ids.bile(sid)
        objects += [
            _species(ids.urine(sid), ids.VURINE, f"{s.name} (urine)", amount=True),
            _species(ids.feces(sid), ids.VFECES, f"{s.name} (feces)", amount=True),
            _species(bile, ids.VBILE, f"{s.name} (bile)"),
        ]
        if not lumen:
            objects.append(
                _transit(
                    ids.bile_reaction(sid),
                    f"{bile} -> {ids.feces(sid)}",
                    f"{ids.K_BILE} * {ids.VBILE} * {bile}",
                    f"transit {sid} bile to feces (intestine without lumen)",
                )
            )
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
