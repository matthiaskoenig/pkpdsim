"""Shared helpers of the perfusion tests.

The tests compile their models at `JIT_OPTIMIZATION` (`O0`) as the PBPK tests
do, with the same autouse fixture.
"""

import dataclasses
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from pbpk.conftest import JIT_OPTIMIZATION, compiled, jit_optimization
from sbmlutils.factory import Compartment, FactoryResult, Model, Parameter, create_model
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion import core, ids
from pkpdsim.perfusion.definition import (
    Circuit,
    GenericOrgan,
    KidneyOrgan,
    LiverOrgan,
    Organ,
    Perfusate,
    PerfusionDefinition,
)
from pkpdsim.perfusion.units import PU
from pkpdsim.species import Species as BuildSpecies

if TYPE_CHECKING:
    import roadrunner

__all__ = [
    "CELLS",
    "JIT_OPTIMIZATION",
    "VASCULAR",
    "ZONE",
    "generic_definition",
    "jit_optimization",
    "kidney_definition",
    "liver_definition",
    "load_circuit",
    "load_zone",
    "zone_model",
]

ZONE = "z1"
CELLS = "Vcell_z1"
VASCULAR = "Vvas_z1"
CELL_VOLUME = 0.0069
"""The cells of a liver of 10 g (l), 0.69 of its volume."""
VASCULAR_VOLUME = 0.1
"""The closed vascular space of the test zone (l), as large as a reservoir."""


def zone_model(
    path: Path,
    params: core.CoreParameters | None = None,
    composition: core.Composition | None = None,
    bile: bool = True,
    zones: tuple[str, ...] = (ZONE,),
    profiles: Mapping[str, Mapping[str, float]] | None = None,
) -> FactoryResult:
    """Write a model of closed zones: the core with fixed vascular spaces.

    The temperature factor `f_T` is a constant 1, the constants of the core
    and of the buffer are defined once; nothing flows in or out.
    """
    params = params or core.CoreParameters.default()
    objects = core.CoreObjects()
    compartments = []
    for zone in zones:
        cells, vascular = f"Vcell_{zone}", f"Vvas_{zone}"
        objects += core.zone_objects(
            zone,
            cells,
            vascular,
            params,
            profile=(profiles or {}).get(zone),
            composition=composition,
            bile=bile,
        )
        compartments += [
            Compartment(
                cells,
                CELL_VOLUME,
                unit=U.litre,
                name=f"cells of {zone}",
                sboTerm=SBO.PHYSICAL_COMPARTMENT,
            ),
            Compartment(
                vascular,
                VASCULAR_VOLUME,
                unit=U.litre,
                name=f"vascular space of {zone}",
                sboTerm=SBO.PHYSICAL_COMPARTMENT,
            ),
        ]
    model = Model(
        "perfusion_zone",
        name="closed zones of a perfused organ",
        units=U,
        model_units=MODEL_UNITS,
        compartments=compartments,
        parameters=[
            Parameter(
                ids.F_T,
                1.0,
                U.dimensionless,
                constant=True,
                name="temperature factor",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            ),
            *core.buffer_parameters(),
            *core.core_parameters(params),
            *objects.parameters,
        ],
        species=objects.species,
        reactions=objects.reactions,
        rules=objects.rules,
        rate_rules=objects.rate_rules,
        assignments=objects.assignments,
    )
    return create_model(model, filepath=path, sbml_level=3, sbml_version=2)


def load_zone(path: Path) -> "roadrunner.RoadRunner":
    """Load a written zone model at `O0` with tight integrator tolerances."""
    r = compiled(path)
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    return r


def liver_definition(
    zones: int = 1,
    *,
    profiles: Mapping[str, Sequence[float]] | None = None,
    tracers: Sequence[str] = (),
    temperature: float = 37.0,
    sid: str = "perfusion_liver",
) -> PerfusionDefinition:
    """Get the rat liver in the circuit of the parameter database."""
    perfusate = Perfusate.default()
    if tracers:
        perfusate = dataclasses.replace(perfusate, tracers=tuple(tracers))
    return PerfusionDefinition(
        sid=sid,
        name="perfused rat liver",
        species=BuildSpecies.RAT,
        organ=LiverOrgan.default(zones=zones, profiles=profiles or {}),
        circuit=Circuit.default(),
        perfusate=perfusate,
        temperature=Value(temperature, PU.degC, "ModelChoice", "the test temperature"),
        core=core.CoreParameters.default(),
    )


def _definition(
    organ: Organ,
    core_section: str,
    *,
    tracers: Sequence[str],
    sid: str,
    name: str,
) -> PerfusionDefinition:
    """Get an organ in the circuit of the parameter database."""
    perfusate = Perfusate.default()
    if tracers:
        perfusate = dataclasses.replace(perfusate, tracers=tuple(tracers))
    return PerfusionDefinition(
        sid=sid,
        name=name,
        species=BuildSpecies.RAT,
        organ=organ,
        circuit=Circuit.default(),
        perfusate=perfusate,
        temperature=Value(37.0, PU.degC, "ModelChoice", "the test temperature"),
        core=core.CoreParameters.default(core_section),
    )


def kidney_definition(
    zones: int = 1, *, tracers: Sequence[str] = (), sid: str = "perfusion_kidney"
) -> PerfusionDefinition:
    """Get the rat kidney in the circuit of the parameter database."""
    return _definition(
        KidneyOrgan.default(zones=zones),
        "kidney_rat",
        tracers=tracers,
        sid=sid,
        name="perfused rat kidney",
    )


def generic_definition(
    key: str = "mu",
    mass: float = 5.0,
    *,
    tracers: Sequence[str] = (),
    sid: str = "perfusion_muscle",
) -> PerfusionDefinition:
    """Get a generic organ (the muscle by default) of `mass` g in the circuit."""
    organ = GenericOrgan(key, Value(mass, PU.g, "ModelChoice", "the test mass"))
    return _definition(
        organ, "generic", tracers=tracers, sid=sid, name=f"perfused rat {key}"
    )


def load_circuit(path: Path) -> "roadrunner.RoadRunner":
    """Load a built circuit at `O0` with tight integrator tolerances."""
    return load_zone(path)
