"""Shared helpers of the PBPK tests.

`simple_drug` and `definition` are plain functions (not fixtures), so a test
varies them with arguments: `definition(routes=[Route.iv("drg")])`. `load` and
`total_amount` need the builder, which exists from the task of the builder on,
and import it inside the functions, so that this module imports before.
"""

from typing import TYPE_CHECKING, Any, Protocol

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

if TYPE_CHECKING:
    import roadrunner


class _Built(Protocol):
    """What `load` needs of a build result (`BuildResult` of the builder)."""

    flat_sbml: str


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


def load(result: "_Built") -> "roadrunner.RoadRunner":
    """Load the flat model of a build result with tight integrator tolerances."""
    import roadrunner

    r = roadrunner.RoadRunner(result.flat_sbml)
    r.integrator.absolute_tolerance = 1e-12
    r.integrator.relative_tolerance = 1e-10
    return r


def total_amount(r: "roadrunner.RoadRunner", sid: str) -> float:
    """Get the total amount (mmole) of a substance: all species and the sinks."""
    model = r.model
    all_ids = [*model.getFloatingSpeciesIds(), *model.getBoundarySpeciesIds()]
    species = ids.species_of(sid, all_ids)
    sinks = [i for i in (ids.urine(sid), ids.feces(sid), ids.bile(sid)) if i in all_ids]
    return float(sum(r[f"{i}"] for i in [*species, *sinks]))
