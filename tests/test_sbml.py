"""The SBML writing kit: the names it holds are reachable at their old paths.

The kit (`pkpdsim.sbml`) is the layer below the PBPK library, the zonation
engine and the machine perfusion library, and it was filled by moving the
pieces the three shared out of them. A model of pkdb_models imports pkpdsim
by its public names, so every moved name keeps working where it was, as a
re-export, for this release; the layering itself is pinned by
`tests/zonation/test_limits.py`.

The table below is the list of what moved, so a name which is moved on and
loses its old path fails here.
"""

from importlib import import_module

import pytest

MOVED: list[tuple[str, str, str, str]] = [
    ("pkpdsim.pbpk.units", "U", "pkpdsim.sbml.units", "U"),
    ("pkpdsim.pbpk.units", "MODEL_UNITS", "pkpdsim.sbml.units", "MODEL_UNITS"),
    ("pkpdsim.perfusion.units", "PU", "pkpdsim.sbml.units", "PU"),
    ("pkpdsim.pbpk.values", "SID", "pkpdsim.sbml.values", "SID"),
    ("pkpdsim.pbpk.values", "Formula", "pkpdsim.sbml.values", "Formula"),
    ("pkpdsim.pbpk.values", "FrozenMap", "pkpdsim.sbml.values", "FrozenMap"),
    ("pkpdsim.pbpk.values", "ModelParameter", "pkpdsim.sbml.values", "ModelParameter"),
    ("pkpdsim.pbpk.values", "Q", "pkpdsim.sbml.values", "Q"),
    ("pkpdsim.pbpk.values", "Value", "pkpdsim.sbml.values", "Value"),
    ("pkpdsim.pbpk.values", "as_value", "pkpdsim.sbml.values", "as_value"),
    ("pkpdsim.pbpk.values", "check_literals", "pkpdsim.sbml.values", "check_literals"),
    ("pkpdsim.pbpk.values", "check_positive", "pkpdsim.sbml.values", "check_positive"),
    ("pkpdsim.pbpk.values", "check_unit", "pkpdsim.sbml.values", "check_unit"),
    ("pkpdsim.pbpk.values", "equal_by_value", "pkpdsim.sbml.values", "equal_by_value"),
    ("pkpdsim.pbpk.values", "fitted", "pkpdsim.sbml.values", "fitted"),
    ("pkpdsim.pbpk.values", "number", "pkpdsim.sbml.values", "number"),
    ("pkpdsim.pbpk.values", "unit_text", "pkpdsim.sbml.values", "unit_text"),
    (
        "pkpdsim.pbpk.substance",
        "DIFFUSION_UNITS",
        "pkpdsim.sbml.values",
        "DIFFUSION_UNITS",
    ),
    ("pkpdsim.pbpk.annotation", "TAXA", "pkpdsim.sbml.annotation", "TAXA"),
    (
        "pkpdsim.pbpk.annotation",
        "initial_assignment",
        "pkpdsim.sbml.annotation",
        "initial_assignment",
    ),
    (
        "pkpdsim.pbpk.placement",
        "value_parameter",
        "pkpdsim.sbml.parameters",
        "value_parameter",
    ),
    (
        "pkpdsim.pbpk.placement",
        "value_notes",
        "pkpdsim.sbml.parameters",
        "value_notes",
    ),
    (
        "pkpdsim.pbpk.placement",
        "FormulaParameter",
        "pkpdsim.sbml.parameters",
        "FormulaParameter",
    ),
    (
        "pkpdsim.pbpk.placement",
        "formula_assignment",
        "pkpdsim.sbml.parameters",
        "formula_assignment",
    ),
    ("pkpdsim.perfusion.observers", "RULE", "pkpdsim.sbml.parameters", "RULE"),
    ("pkpdsim.perfusion.observers", "observer", "pkpdsim.sbml.parameters", "observer"),
    ("pkpdsim.perfusion.core", "CoreObjects", "pkpdsim.sbml.objects", "SbmlObjects"),
    ("pkpdsim.perfusion.ids", "flow", "pkpdsim.zonation.ids", "flow"),
    ("pkpdsim.perfusion.ids", "is_zone", "pkpdsim.zonation.ids", "is_zone"),
]
"""What moved into the kit (and into the engine): the old module and name, the
new module and name."""


@pytest.mark.parametrize(
    ("old_module", "old_name", "new_module", "new_name"),
    MOVED,
    ids=[f"{m}.{n}" for m, n, _, _ in MOVED],
)
def test_a_moved_name_is_the_same_object_at_its_old_path(
    old_module: str, old_name: str, new_module: str, new_name: str
) -> None:
    """The old path gives the object of the new one, not a copy of it."""
    old = getattr(import_module(old_module), old_name)
    new = getattr(import_module(new_module), new_name)
    assert old is new


@pytest.mark.parametrize(
    ("old_module", "old_name"),
    [(m, n) for m, n, _, _ in MOVED],
    ids=[f"{m}.{n}" for m, n, _, _ in MOVED],
)
def test_a_moved_name_stays_public_at_its_old_path(
    old_module: str, old_name: str
) -> None:
    """The name is in the `__all__` of the module it moved out of.

    A re-export which is not public is one a caller is told not to use, which
    is not what a kept import path is for.
    """
    assert old_name in import_module(old_module).__all__
