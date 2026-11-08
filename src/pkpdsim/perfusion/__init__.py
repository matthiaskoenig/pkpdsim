"""The machine perfusion library: a circuit with a mounted organ (pkpdsim#146).

An organ in a machine perfusion is a chain of zones, each a vascular space
in series with the flow and its cells with the cellular core
(`pkpdsim.perfusion.core`: the energy metabolism, the redox and the damage of
the cells, the release of their enzymes and potassium, the bicarbonate buffer
of the perfusate). The ids are those of `pkpdsim.perfusion.ids`.

The names of the modules are exported lazily (PEP 562), so `import
pkpdsim.perfusion` imports neither sbmlutils nor the PBPK library.
"""

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pkpdsim.perfusion.core import (
        PROCESSES,
        Buffer,
        Composition,
        CoreObjects,
        CoreParameters,
        buffer_objects,
        buffer_parameters,
        core_parameters,
        perfusate_objects,
        zone_objects,
    )

#: the names imported on their first access: name -> module of `pkpdsim.perfusion`
_LAZY_NAMES = {
    "PROCESSES": "core",
    "Buffer": "core",
    "Composition": "core",
    "CoreObjects": "core",
    "CoreParameters": "core",
    "buffer_objects": "core",
    "buffer_parameters": "core",
    "core_parameters": "core",
    "perfusate_objects": "core",
    "zone_objects": "core",
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
    "Composition",
    "CoreObjects",
    "CoreParameters",
    "buffer_objects",
    "buffer_parameters",
    "core_parameters",
    "perfusate_objects",
    "zone_objects",
]
