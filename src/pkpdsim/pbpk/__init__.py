"""The PBPK model library: the building blocks of the PBPK models of a substance.

The base of the experiments of a model built on the library
(`PBPKSimulationExperiment`, `ScaledRenalMap`, `healthy_gfr`, `REFERENCE_GFR`
of `pkpdsim.pbpk.experiment`) is exported lazily: it imports sbmlsim, which a
factory does not need, so `import pkpdsim.pbpk` imports no simulator.
"""

from typing import TYPE_CHECKING

from pkpdsim.pbpk.archive import create_archive
from pkpdsim.pbpk.builder import BuildResult, build, compose
from pkpdsim.pbpk.chemistry import (
    MM,
    REST,
    CustomReaction,
    Direction,
    FirstOrder,
    Hill,
    Membrane,
    Metabolism,
    RateOf,
    Rest,
    Transport,
)
from pkpdsim.pbpk.data import Source, cite, load, sources
from pkpdsim.pbpk.definition import PBPKDefinition, Variants
from pkpdsim.pbpk.distribution import interstitial_share
from pkpdsim.pbpk.observers import Observer
from pkpdsim.pbpk.organ import (
    Boundary,
    GutComponent,
    OrganComponent,
    OrganContext,
    StandardOrgan,
)
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.parts import Parts, parts_for
from pkpdsim.pbpk.pd import PDComponent
from pkpdsim.pbpk.physiology import (
    ReferenceIndividual,
    Sex,
    physiology_objects,
    reference_individual,
)
from pkpdsim.pbpk.proteins import PROTEINS, Induction, Inhibition, Protein
from pkpdsim.pbpk.routes import (
    DoseDependentDissolution,
    FirstOrderDissolution,
    Route,
    WeibullDissolution,
)
from pkpdsim.pbpk.slow_pool import SlowPool
from pkpdsim.pbpk.substance import (
    Binding,
    Erythrocytes,
    Excretion,
    Kp,
    KpMethod,
    Mode,
    Saturable,
    Substance,
)
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Formula, ModelParameter, Q, Value, as_value, fitted
from pkpdsim.pbpk.variants.intestine import SegmentedIntestine, SingleLumen
from pkpdsim.pbpk.variants.kidney_nephron import NephronKidney
from pkpdsim.pbpk.variants.liver_zonated import ZonatedLiver

if TYPE_CHECKING:
    from pkpdsim.pbpk.experiment import (
        REFERENCE_GFR,
        KidneyOutflow,
        PBPKSimulationExperiment,
        ScaledRenalMap,
        healthy_gfr,
        kidney_outflow,
    )

#: the names of `pkpdsim.pbpk.experiment`, imported on their first access
_EXPERIMENT_NAMES = frozenset(
    {
        "REFERENCE_GFR",
        "KidneyOutflow",
        "PBPKSimulationExperiment",
        "ScaledRenalMap",
        "healthy_gfr",
        "kidney_outflow",
    }
)


def __getattr__(name: str) -> object:
    """Import a name of `pkpdsim.pbpk.experiment` on its first access (PEP 562).

    Raises:
        AttributeError: for any other name.
    """
    if name in _EXPERIMENT_NAMES:
        from pkpdsim.pbpk import experiment

        return getattr(experiment, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """List the names of the module with the names imported on their first access."""
    return sorted({*globals(), *_EXPERIMENT_NAMES})


__all__ = [
    "MM",
    "MODEL_UNITS",
    "PROTEINS",
    "REFERENCE_GFR",
    "REST",
    "Binding",
    "Boundary",
    "BuildResult",
    "CustomReaction",
    "Direction",
    "DoseDependentDissolution",
    "Erythrocytes",
    "Excretion",
    "FirstOrder",
    "FirstOrderDissolution",
    "Formula",
    "GutComponent",
    "Hill",
    "Induction",
    "Inhibition",
    "KidneyOutflow",
    "Kp",
    "KpMethod",
    "Membrane",
    "Metabolism",
    "Mode",
    "ModelParameter",
    "NephronKidney",
    "Observer",
    "Organ",
    "OrganComponent",
    "OrganContext",
    "PBPKDefinition",
    "PBPKSimulationExperiment",
    "PDComponent",
    "Part",
    "Parts",
    "Pool",
    "Protein",
    "Q",
    "RateOf",
    "ReferenceIndividual",
    "Rest",
    "Route",
    "Saturable",
    "ScaledRenalMap",
    "SegmentedIntestine",
    "Sex",
    "SingleLumen",
    "SlowPool",
    "Source",
    "StandardOrgan",
    "Substance",
    "Transport",
    "U",
    "Value",
    "Variants",
    "WeibullDissolution",
    "ZonatedLiver",
    "as_value",
    "build",
    "cite",
    "compose",
    "create_archive",
    "fitted",
    "healthy_gfr",
    "interstitial_share",
    "kidney_outflow",
    "load",
    "parts_for",
    "physiology_objects",
    "reference_individual",
    "sources",
]
