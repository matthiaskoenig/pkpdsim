"""The machine perfusion library: a circuit with a mounted organ (pkpdsim#146).

A model is built from a definition (`pkpdsim.perfusion.definition`: the
circuit, the perfusate, the organ, the temperature) by
`pkpdsim.perfusion.builder.build`, one SBML file;
`pkpdsim.perfusion.archive.create_archive` writes it with its ODEs, its
antimony, the figures of the circuit (`pkpdsim.perfusion.figure`, the
schematic, and `pkpdsim.perfusion.anatomy`, the organ in its circuit) and
the COMBINE archive; `pkpdsim.perfusion.experiment` is the base of the
simulation experiments of a model, with the condition of a perfusion and its
medium changes. An organ in a machine perfusion is a chain of zones,
each a vascular space in series with the flow and its cells with the cellular core
(`pkpdsim.perfusion.core`: the energy metabolism, the redox and the damage of
the cells, the release of their enzymes and potassium, the bicarbonate buffer
of the perfusate). The ids are those of `pkpdsim.perfusion.ids`.

The names of the modules are exported lazily (PEP 562), so `import
pkpdsim.perfusion` imports neither sbmlutils nor the PBPK library.
"""

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pkpdsim.perfusion.anatomy import (
        circuit_anatomy_figure,
        circuit_anatomy_svg,
    )
    from pkpdsim.perfusion.archive import create_archive
    from pkpdsim.perfusion.builder import PerfusionExtension, build, compose
    from pkpdsim.perfusion.core import (
        PROCESSES,
        Buffer,
        Composition,
        CoreNames,
        CoreObjects,
        CoreParameters,
        buffer_objects,
        buffer_parameters,
        core_objects,
        core_parameters,
        perfusate_objects,
        quantity_name,
    )
    from pkpdsim.perfusion.definition import (
        Circuit,
        GenericOrgan,
        KidneyOrgan,
        LiverOrgan,
        Organ,
        Perfusate,
        PerfusionDefinition,
    )
    from pkpdsim.perfusion.experiment import (
        PerfusionSimulationExperiment,
        medium_change,
        observable_labels,
        observable_units,
        perfusion_changes,
        site_label,
        temperature_quantities,
    )
    from pkpdsim.perfusion.figure import (
        circuit_cards,
        circuit_figure,
        circuit_svg,
        interactive_circuit,
    )

#: the names imported on their first access: name -> module of `pkpdsim.perfusion`
_LAZY_NAMES = {
    "PROCESSES": "core",
    "Buffer": "core",
    "Composition": "core",
    "CoreNames": "core",
    "CoreObjects": "core",
    "CoreParameters": "core",
    "buffer_objects": "core",
    "buffer_parameters": "core",
    "core_objects": "core",
    "core_parameters": "core",
    "perfusate_objects": "core",
    "quantity_name": "core",
    "Circuit": "definition",
    "GenericOrgan": "definition",
    "KidneyOrgan": "definition",
    "LiverOrgan": "definition",
    "Organ": "definition",
    "Perfusate": "definition",
    "PerfusionDefinition": "definition",
    "circuit_anatomy_figure": "anatomy",
    "circuit_anatomy_svg": "anatomy",
    "create_archive": "archive",
    "PerfusionSimulationExperiment": "experiment",
    "medium_change": "experiment",
    "observable_labels": "experiment",
    "observable_units": "experiment",
    "perfusion_changes": "experiment",
    "site_label": "experiment",
    "temperature_quantities": "experiment",
    "circuit_cards": "figure",
    "circuit_figure": "figure",
    "circuit_svg": "figure",
    "interactive_circuit": "figure",
    "PerfusionExtension": "builder",
    "build": "builder",
    "compose": "builder",
}


def __getattr__(name: str) -> object:
    """Import a name of a module of the library on its first access (PEP 562).

    Raises:
        AttributeError: for any other name.
    """
    module = _LAZY_NAMES.get(name)
    if module is not None:
        return getattr(importlib.import_module(f"pkpdsim.perfusion.{module}"), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """List the names of the module with the names imported on their first access."""
    return sorted({*globals(), *_LAZY_NAMES})


__all__ = [
    "PROCESSES",
    "Buffer",
    "Circuit",
    "Composition",
    "CoreNames",
    "CoreObjects",
    "CoreParameters",
    "GenericOrgan",
    "KidneyOrgan",
    "LiverOrgan",
    "Organ",
    "Perfusate",
    "PerfusionDefinition",
    "PerfusionExtension",
    "PerfusionSimulationExperiment",
    "buffer_objects",
    "buffer_parameters",
    "build",
    "circuit_anatomy_figure",
    "circuit_anatomy_svg",
    "circuit_cards",
    "circuit_figure",
    "circuit_svg",
    "compose",
    "core_objects",
    "core_parameters",
    "create_archive",
    "interactive_circuit",
    "medium_change",
    "observable_labels",
    "observable_units",
    "perfusate_objects",
    "perfusion_changes",
    "quantity_name",
    "site_label",
    "temperature_quantities",
]
