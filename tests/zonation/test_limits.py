"""The limits of the engine build: valid SBML of a plain chain, and the imports.

The engine models of `test_chain.py`, `test_barrier.py`, `test_bile.py` and
`test_flow.py` are written with `validate=False` and simulated, so a model
roadrunner loads can still carry an invalid SBO term or annotation. This
module validates a plain build (the chain alone, no composition and no cell
model) with libsbml, as `test_observers.py` validates the composed model.

The engine and the PBPK library import each other (the engine builds on the
substances, the units and the placement of the library, the zonated liver of
the library builds on the geometry and the transport of the engine), and both
report the zonal means of the perfusion library's weighting, so the three
must import in any order: the tests of the engine import the library first,
which hides a cycle a user meets who imports the engine alone. The guard
walks the modules of the three packages instead of listing them
(`_entry_points`), so a new module is guarded the day it is written.
"""

import importlib
import pkgutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import libsbml
import pytest
from sbmlutils.factory import create_model

from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.zonation import ZonationBuilder
from pkpdsim.zonation.flow import ImposedFlow
from pkpdsim.zonation.settings import BileTransport, Geometry, ZonationSettings

from .conftest import mass_g

pytestmark = pytest.mark.filterwarnings("ignore:.*numerical dispersion:UserWarning")

DRG = Substance("drg", "drug", Value(180.0, U.g_per_mole))
"""A small solute, which every layer and the bile carry."""


@pytest.mark.parametrize(
    "geometry", [Geometry.SINUSOID, Geometry.WELL_STIRRED], ids=["sinusoid", "stirred"]
)
def test_a_plain_build_is_valid_sbml(tmp_path: Path, geometry: Geometry) -> None:
    """libsbml reports no error and no warning for the chain of the engine alone."""
    settings = ZonationSettings(
        n_zones=4,
        geometry=geometry,
        flow=ImposedFlow(1.0),
        bile=BileTransport.DIFFUSIVE,
        liver_mass=mass_g(1000.0),
    )
    path = tmp_path / "chain.xml"
    create_model(
        ZonationBuilder(settings, [DRG], None).build(),
        filepath=path,
        sbml_level=3,
        sbml_version=2,
        validate=False,
    )
    document = libsbml.readSBMLFromFile(str(path))
    document.checkConsistency()
    problems = [
        document.getError(i).getMessage()[:400] for i in range(document.getNumErrors())
    ]
    assert problems == []  # errors and warnings both count


PACKAGES = ("pkpdsim.pbpk", "pkpdsim.zonation", "pkpdsim.perfusion")
"""The three libraries whose modules import each other."""


def _entry_points() -> list[str]:
    """Every module of the three libraries, each one a possible first import.

    The modules are collected and not listed: a list only guards the modules
    somebody remembered to add to it, and the module which broke the imports
    last was one nobody had listed. A subpackage which does not import is no
    reason to leave it out of the guard, so the walk raises instead of
    skipping it.

    Returns:
        The dotted names of the three packages and of every module below
        them, sorted (the ids of the parametrization, so every worker of
        `pytest-xdist` sees the same ones).
    """

    def fail(name: str) -> None:
        """Refuse to skip a module which does not import."""
        raise ImportError(name)

    out: list[str] = []
    for name in PACKAGES:
        package: ModuleType = importlib.import_module(name)
        out.append(name)
        out += [
            m.name
            for m in pkgutil.walk_packages(package.__path__, f"{name}.", onerror=fail)
        ]
    return sorted(out)


@pytest.mark.parametrize("first", _entry_points())
def test_the_packages_import_in_any_order(first: str) -> None:
    """Every module of the three libraries imports first in a fresh interpreter.

    The engine builds on the PBPK library and the zonated liver of the PBPK
    library builds on the engine, and both report the means of the perfusion
    library's weighting, so the modules of the three import each other and a
    module below one of them can find a partially initialized module
    (`ImportError`). Every module is imported first in a subprocess of its
    own, the other two libraries after it: importing a module of a package
    imports the package, so a cycle through any of the three shows up
    wherever a user enters.
    """
    others = [p for p in PACKAGES if not first.startswith(p)]
    code = "; ".join(f"import {m}" for m in (first, *others))
    done = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert done.returncode == 0, done.stderr[-2000:]
