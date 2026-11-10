"""The SBML writing kit: the names it holds are reachable at their old paths.

The kit (`pkpdsim.sbml`) is the layer below the PBPK library, the zonation
engine and the machine perfusion library, and it was filled by moving the
pieces the three shared out of them. A model of pkdb_models imports pkpdsim
by its public names, so every moved name keeps working where it was, as a
re-export, for this release; the layering itself is pinned by
`tests/zonation/test_limits.py`.

The checks are generated and not listed: the public names of the kit
(`KIT_MODULES`) are collected from their `__all__`, and every module which
re-exports them (`SHIMS`) is walked over its own `__all__`. A name added to a
shim later is checked the day it is added, which a hand-written table of
names would not be.
"""

from collections.abc import Mapping
from importlib import import_module
from typing import Any

import pytest

KIT_MODULES = (
    "pkpdsim.sbml.annotation",
    "pkpdsim.sbml.objects",
    "pkpdsim.sbml.parameters",
    "pkpdsim.sbml.units",
    "pkpdsim.sbml.values",
)
"""The modules of the kit, which own the moved names."""

SHIMS = (
    "pkpdsim.pbpk.annotation",
    "pkpdsim.pbpk.placement",
    "pkpdsim.pbpk.substance",
    "pkpdsim.pbpk.units",
    "pkpdsim.pbpk.values",
    "pkpdsim.perfusion.core",
    "pkpdsim.perfusion.observers",
    "pkpdsim.perfusion.units",
)
"""The modules a name moved out of, which re-export it at its old path."""

RENAMED: Mapping[str, str] = {"CoreObjects": "SbmlObjects"}
"""An old name which the kit holds under a new one."""

ALIASES = (("pkpdsim.sbml.units", "PU", "OrganUnits"),)
"""Two public names of one kit module which are the same object.

`PU` is the name the machine perfusion library has always called the units of
an organ by, `OrganUnits` the neutral name of the kit. An alias and not a
subclass: a subclass would be a second class, hence a second set of unit
definitions, hence a different model."""

MOVED_TO_ENGINE: tuple[tuple[str, str, str], ...] = (
    ("pkpdsim.perfusion.ids", "flow", "pkpdsim.zonation.ids"),
    ("pkpdsim.perfusion.ids", "is_zone", "pkpdsim.zonation.ids"),
)
"""What moved from the perfusion library into the engine, by name.

These two are listed and not walked: the two id modules share many public
names (`INFLOW`, `TEMPERATURE`, `volume`, `species`, ...) which are
deliberately different objects, the ids of one package and the ids of the
other, so a walk over the names they have in common would check the wrong
thing."""


def _kit_names() -> dict[str, Any]:
    """The public names of the kit, by name.

    Returns:
        The name of every entry of the `__all__` of a module of
        `KIT_MODULES`, mapped to its object.
    """
    out: dict[str, Any] = {}
    owner: dict[str, str] = {}
    for name in KIT_MODULES:
        module = import_module(name)
        for public in module.__all__:
            assert public not in out, (
                f"'{public}' is public in both {owner[public]} and {name}"
            )
            out[public] = getattr(module, public)
            owner[public] = name
    return out


KIT = _kit_names()
"""The public names of the kit, by name."""


def _reexports(shim: str) -> list[str]:
    """The names of the kit which a shim module re-exports, sorted."""
    public = import_module(shim).__all__
    return sorted(n for n in public if RENAMED.get(n, n) in KIT)


PAIRS = [(shim, name) for shim in SHIMS for name in _reexports(shim)]
"""Every shim module and every name of the kit it re-exports."""


@pytest.mark.parametrize(("shim", "name"), PAIRS, ids=[f"{m}.{n}" for m, n in PAIRS])
def test_a_name_of_the_kit_is_the_same_object_at_its_old_path(
    shim: str, name: str
) -> None:
    """The old path gives the object of the kit, not a copy of it."""
    assert getattr(import_module(shim), name) is KIT[RENAMED.get(name, name)]


@pytest.mark.parametrize("shim", SHIMS)
def test_a_shim_module_re_exports_a_name_of_the_kit(shim: str) -> None:
    """Every module of `SHIMS` is a shim, so the walk above cannot pass empty.

    A shim whose `__all__` loses its re-exports breaks the import path of a
    caller, and would otherwise contribute no check at all.
    """
    assert _reexports(shim) != []


@pytest.mark.parametrize(
    ("module", "name", "alias"), ALIASES, ids=[f"{m}.{n}" for m, n, _ in ALIASES]
)
def test_an_alias_of_the_kit_is_the_object_it_aliases(
    module: str, name: str, alias: str
) -> None:
    """Both names of one object are public and are that one object."""
    kit = import_module(module)
    assert name in kit.__all__
    assert alias in kit.__all__
    assert getattr(kit, name) is getattr(kit, alias)


@pytest.mark.parametrize(
    ("old_module", "name", "new_module"),
    MOVED_TO_ENGINE,
    ids=[f"{m}.{n}" for m, n, _ in MOVED_TO_ENGINE],
)
def test_a_name_moved_into_the_engine_is_the_same_object_at_its_old_path(
    old_module: str, name: str, new_module: str
) -> None:
    """The name is public at both paths and is one object."""
    old = import_module(old_module)
    assert name in old.__all__
    assert getattr(old, name) is getattr(import_module(new_module), name)
