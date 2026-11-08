"""Shared helpers of the perfusion tests.

The tests compile their models at `JIT_OPTIMIZATION` (`O0`) as the PBPK tests
do, with the same autouse fixture.
"""

from pathlib import Path
from typing import TYPE_CHECKING

from pbpk.conftest import JIT_OPTIMIZATION, compiled, jit_optimization
from sbmlutils.factory import Compartment, FactoryResult, Model, Parameter, create_model
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.perfusion import core, ids

if TYPE_CHECKING:
    import roadrunner

__all__ = [
    "CELLS",
    "JIT_OPTIMIZATION",
    "VASCULAR",
    "ZONE",
    "jit_optimization",
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
) -> FactoryResult:
    """Write a model of one closed zone: the core with a fixed vascular space.

    The temperature factor `f_T` is a constant 1, the buffer constants are
    the defaults; nothing flows in or out.
    """
    objects = core.zone_objects(
        ZONE,
        CELLS,
        VASCULAR,
        params or core.CoreParameters.default(),
        composition=composition,
        bile=bile,
    )
    model = Model(
        "perfusion_zone",
        name="one zone of a perfused organ",
        units=U,
        model_units=MODEL_UNITS,
        compartments=[
            Compartment(
                CELLS,
                CELL_VOLUME,
                unit=U.litre,
                name="cells of zone 1",
                sboTerm=SBO.PHYSICAL_COMPARTMENT,
            ),
            Compartment(
                VASCULAR,
                VASCULAR_VOLUME,
                unit=U.litre,
                name="vascular space of zone 1",
                sboTerm=SBO.PHYSICAL_COMPARTMENT,
            ),
        ],
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
