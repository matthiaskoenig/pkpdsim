"""Shared helpers of the PBPK tests.

`simple_drug` and `definition` are plain functions (not fixtures), so a test
varies them with arguments: `definition(routes=[Route.iv("drg")])`.
"""

from typing import TYPE_CHECKING, Any

from pkpdsim.pbpk import (
    MM,
    Binding,
    Excretion,
    Kp,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    ids,
)
from pkpdsim.pbpk.builder import BuildResult

if TYPE_CHECKING:
    import roadrunner


def simple_drug() -> Substance:
    """Get the simple drug: Mr 300, fu 0.1 (albumin), Kp 2, merged tissue, urine."""
    return Substance(
        sid="drg",
        name="drug",
        mr=300.0,
        binding=Binding(fu=0.1, protein="albumin"),
        distribution=Kp(2.0, ftissue=1.0),
        excretion=Excretion(urine=True),
    )


def definition(**overrides: Any) -> PBPKDefinition:  # noqa: ANN401
    """Get the definition of the simple drug with a liver metabolism `DRG2M`.

    Args:
        **overrides: fields of `PBPKDefinition` which replace the defaults.

    Returns:
        The definition.
    """
    fields: dict[str, Any] = {
        "sid": "simple",
        "substances": [
            simple_drug(),
            Substance(sid="drg_m", name="metabolite", mr=316.0),
        ],
        "chemistry": [
            Metabolism(
                sid="DRG2M",
                organ=Organ.LIVER,
                substrate="drg",
                products={"drg_m": 1.0},
                kinetics=MM(vmax=1e-3, km=0.01),
                protein="CYP3A4",
            )
        ],
        "routes": [Route.iv("drg"), Route.oral("drg")],
    }
    fields.update(overrides)
    return PBPKDefinition(**fields)


def load(result: BuildResult) -> "roadrunner.RoadRunner":
    """Load the flat model of a build result with tight integrator tolerances.

    The model is compiled with the LLJIT backend of roadrunner: the MCJIT
    backend (the default) takes about ten times longer for the deep rules of
    the physiology.
    """
    import roadrunner

    config = roadrunner.Config
    backend = config.getValue(config.LLVM_BACKEND)
    config.setValue(config.LLVM_BACKEND, config.LLJIT)
    try:
        r = roadrunner.RoadRunner(str(result.flat_path))
    finally:
        config.setValue(config.LLVM_BACKEND, backend)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    return r


def total_amount(r: "roadrunner.RoadRunner", sid: str) -> float:
    """Get the total amount (mmole) of a substance: all species and the sinks.

    The amount of a species is decided from the SBML: a species with
    `hasOnlySubstanceUnits` is an amount, any other a concentration which is
    multiplied with the volume of its compartment. The sinks are the urine and
    feces amounts, the bile and the intestinal lumen; the dose states are not
    counted (`dose_amount`).
    """
    import libsbml

    doc = libsbml.readSBMLFromString(r.getCurrentSBML())
    model = doc.getModel()
    by_id = {sp.getId(): sp for sp in model.getListOfSpecies()}
    sinks = {ids.urine(sid), ids.feces(sid), ids.bile(sid), ids.lumen(sid)}
    selected = ids.species_of(sid, by_id) + [
        i for i in by_id if i.rsplit("__", 1)[-1] in sinks
    ]
    total = 0.0
    for i in selected:
        sp = by_id[i]
        if sp.getHasOnlySubstanceUnits():
            total += float(r[i])
        else:
            total += float(r[f"[{i}]"]) * float(r[sp.getCompartment()])
    return total


def dose_amount(r: "roadrunner.RoadRunner", sid: str) -> float:
    """Get the amount (mmole) of a substance in its dose states, not yet applied."""
    doses = (ids.iv_dose(sid), ids.po_dose(sid), ids.sol_dose(sid))
    mg = sum(float(r[d]) for d in doses if d in r.model.getGlobalParameterIds())
    return mg / float(r[ids.mr(sid)])
