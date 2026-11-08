"""The PBPK model library: the building blocks of the PBPK models of a substance.

The base of the experiments of a model built on the library
(`PBPKSimulationExperiment`, `ScaledRenalMap`, `healthy_gfr`, `REFERENCE_GFR`,
`molecular_weights` of `pkpdsim.pbpk.experiment`), the equilibration outside
an experiment (`equilibrated_simulation` of `pkpdsim.pbpk.equilibration`)
and the validity of the kidney outflow (`kidney_outflow`, `KidneyOutflow` of
`pkpdsim.pbpk.kidney_validity`) are
exported lazily: they import sbmlsim, which a factory does not need, so
`import pkpdsim.pbpk` imports no simulator.
"""

import importlib
from typing import TYPE_CHECKING

from pkpdsim.pbpk.adipose import (
    Adipose,
    VisceralComposition,
    visceral_composition,
    waist_factor,
)
from pkpdsim.pbpk.archive import create_archive, sex_definition
from pkpdsim.pbpk.builder import BuildResult, build, compose
from pkpdsim.pbpk.chemistry import (
    LINEAR,
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
    Saturable,
    Transport,
)
from pkpdsim.pbpk.circadian import Circadian
from pkpdsim.pbpk.data import Source, cite, load, sources
from pkpdsim.pbpk.definition import PBPKDefinition, Variants
from pkpdsim.pbpk.distribution import KCellFloorWarning, interstitial_share
from pkpdsim.pbpk.figure import model_figure
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
    BodyComposition,
    ReferenceIndividual,
    Sex,
    body_composition,
    physiology_objects,
    reference_individual,
)
from pkpdsim.pbpk.pregnancy import Pregnancy
from pkpdsim.pbpk.proteins import PROTEINS, Induction, Inhibition, Protein
from pkpdsim.pbpk.regional import (
    derive_regional_absorption,
    regional_key,
    release_auc_ratios,
)
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
    Substance,
    ThresholdReabsorption,
)
from pkpdsim.pbpk.units import MODEL_UNITS, UNIT_LABELS, U, axis, unit_label
from pkpdsim.pbpk.values import Formula, ModelParameter, Q, Value, as_value, fitted
from pkpdsim.pbpk.variants.intestine import SegmentedIntestine, SingleLumen
from pkpdsim.pbpk.variants.kidney_nephron import NephronKidney
from pkpdsim.pbpk.variants.liver_zonated import ZonatedLiver
from pkpdsim.pbpk.variants.lumen_options import LumenConversion, RegionalAbsorption
from pkpdsim.pbpk.virtual import (
    VirtualIndividual,
    VirtualPopulation,
    virtual_population,
)
from pkpdsim.species import Species

if TYPE_CHECKING:
    from pkpdsim.pbpk.equilibration import Equilibration, equilibrated_simulation
    from pkpdsim.pbpk.experiment import (
        REFERENCE_GFR,
        PBPKSimulationExperiment,
        ScaledRenalMap,
        healthy_gfr,
        molecular_weights,
    )
    from pkpdsim.pbpk.kidney_validity import KidneyOutflow, kidney_outflow

#: the names of the modules which import sbmlsim (or draw the figures), imported
#: on their first access: name -> module of `pkpdsim.pbpk`
_LAZY_NAMES = {
    "REFERENCE_GFR": "experiment",
    "PBPKSimulationExperiment": "experiment",
    "ScaledRenalMap": "experiment",
    "Equilibration": "equilibration",
    "equilibrated_simulation": "equilibration",
    "healthy_gfr": "experiment",
    "molecular_weights": "experiment",
    "KidneyOutflow": "kidney_validity",
    "kidney_outflow": "kidney_validity",
    # the figures: numpy, scipy and the raster of the bodies, not needed to build
    "anatomy_figure": "visual.anatomy",
    "body_strip": "visual.body",
    "interactive_figure": "visual.interactive",
    "publication_figure": "visual.publication",
}


def __getattr__(name: str) -> object:
    """Import a name of a module which imports sbmlsim on its first access (PEP 562).

    Raises:
        AttributeError: for any other name.
    """
    module = _LAZY_NAMES.get(name)
    if module is not None:
        return getattr(importlib.import_module(f"pkpdsim.pbpk.{module}"), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """List the names of the module with the names imported on their first access."""
    return sorted({*globals(), *_LAZY_NAMES})


__all__ = [
    "LINEAR",
    "MM",
    "MODEL_UNITS",
    "PROTEINS",
    "REFERENCE_GFR",
    "REST",
    "UNIT_LABELS",
    "Adipose",
    "Binding",
    "BodyComposition",
    "Boundary",
    "BuildResult",
    "Circadian",
    "CustomReaction",
    "Direction",
    "DoseDependentDissolution",
    "Equilibration",
    "Erythrocytes",
    "Excretion",
    "FirstOrder",
    "FirstOrderDissolution",
    "Formula",
    "GutComponent",
    "Hill",
    "Induction",
    "Inhibition",
    "KCellFloorWarning",
    "KidneyOutflow",
    "Kp",
    "KpMethod",
    "LumenConversion",
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
    "Pregnancy",
    "Protein",
    "Q",
    "RateOf",
    "ReferenceIndividual",
    "RegionalAbsorption",
    "Rest",
    "Route",
    "Saturable",
    "ScaledRenalMap",
    "SegmentedIntestine",
    "Sex",
    "SingleLumen",
    "SlowPool",
    "Source",
    "Species",
    "StandardOrgan",
    "Substance",
    "ThresholdReabsorption",
    "Transport",
    "U",
    "Value",
    "Variants",
    "VirtualIndividual",
    "VirtualPopulation",
    "VisceralComposition",
    "WeibullDissolution",
    "ZonatedLiver",
    "anatomy_figure",
    "as_value",
    "axis",
    "body_composition",
    "body_strip",
    "build",
    "cite",
    "compose",
    "create_archive",
    "derive_regional_absorption",
    "equilibrated_simulation",
    "fitted",
    "healthy_gfr",
    "interactive_figure",
    "interstitial_share",
    "kidney_outflow",
    "load",
    "model_figure",
    "molecular_weights",
    "parts_for",
    "physiology_objects",
    "publication_figure",
    "reference_individual",
    "regional_key",
    "release_auc_ratios",
    "sex_definition",
    "sources",
    "unit_label",
    "virtual_population",
    "visceral_composition",
    "waist_factor",
]
