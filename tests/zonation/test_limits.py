"""The limits of the engine build: valid SBML of a plain chain, and the layering.

The engine models of `test_chain.py`, `test_barrier.py`, `test_bile.py` and
`test_flow.py` are written with `validate=False` and simulated, so a model
roadrunner loads can still carry an invalid SBO term or annotation. This
module validates a plain build (the chain alone, no composition and no cell
model) with libsbml, as `test_observers.py` validates the composed model.

The four packages which write SBML are layered: the writing kit
(`pkpdsim.sbml`: the units, the values with their provenance, the parameters
and observers they become, the initial assignment and the taxon, the
container of the objects), then the zonation engine, then the PBPK library,
then the machine perfusion library. A module imports its own package and the
packages below it, never one above: the kit is below all three, the engine
names the ids of a zone and a flow which the perfusion library takes from it,
and the zonated liver of the PBPK library and the organs of the perfusion
library both build on the engine. An order is stronger than the former guard
of this module, which started every module of three packages in a subprocess
and checked that the import went through: an order has no cycle to find, in
any entry order. The guard walks the packages instead of listing their
modules (`_modules`), so a new module is guarded the day it is written, and
it imports every module it reads, so a module which does not import at all
still fails here.
"""

import ast
import importlib
import pkgutil
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


LAYERS = ("pkpdsim.sbml", "pkpdsim.zonation", "pkpdsim.pbpk", "pkpdsim.perfusion")
"""The packages which write SBML, from the lowest layer up.

A module of `LAYERS[i]` imports its own package and the packages `LAYERS[:i]`,
never one of `LAYERS[i + 1:]`."""


def _modules() -> list[str]:
    """Every module of the layered packages.

    The modules are collected and not listed: a list only guards the modules
    somebody remembered to add to it, and the module which broke the layering
    last was one nobody had listed. A subpackage which does not import is no
    reason to leave it out of the guard, so the walk raises instead of
    skipping it.

    Returns:
        The dotted names of the packages and of every module below them,
        sorted (the ids of the parametrization, so every worker of
        `pytest-xdist` sees the same ones).
    """

    def fail(name: str) -> None:
        """Refuse to skip a module which does not import."""
        raise ImportError(name)

    out: list[str] = []
    for name in LAYERS:
        package: ModuleType = importlib.import_module(name)
        out.append(name)
        out += [
            m.name
            for m in pkgutil.walk_packages(package.__path__, f"{name}.", onerror=fail)
        ]
    return sorted(out)


def _imported(tree: ast.Module) -> set[str]:
    """The modules a parsed module imports, by their dotted name.

    Every `import a.b`, `from a.b import c` and `importlib.import_module("a.b")`
    of the file, those inside a function included: an import which is deferred
    to break a cycle is the layering violation this guard is about. A relative
    import stays inside the package of the module and is no violation.

    Args:
        tree: the parsed module.

    Returns:
        The dotted names.
    """
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            if node.level == 0 and node.module is not None:
                out.add(node.module)
        elif (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute | ast.Name)
            and _called(node.func) == "import_module"
        ):
            out |= {
                a.value
                for a in node.args
                if isinstance(a, ast.Constant) and isinstance(a.value, str)
            }
    return out


def _called(func: ast.Attribute | ast.Name) -> str:
    """The name of a called function, without the module it is called on."""
    return func.attr if isinstance(func, ast.Attribute) else func.id


@pytest.mark.parametrize("name", _modules())
def test_a_module_imports_only_lower_layers(name: str) -> None:
    """Every module of the four packages imports no package above its own layer.

    The kit is below the three libraries, the engine below the PBPK and the
    perfusion library. The module is imported (a module which does not import
    fails here) and its file is parsed: the imports are read and not run, so a
    deferred import inside a function counts as well.
    """
    module = importlib.import_module(name)
    layer = next(
        i for i, p in enumerate(LAYERS) if name == p or name.startswith(f"{p}.")
    )
    above = LAYERS[layer + 1 :]
    assert module.__file__ is not None
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    forbidden = sorted(
        m for m in _imported(tree) for p in above if m == p or m.startswith(f"{p}.")
    )
    assert forbidden == [], (
        f"{name} is in the layer '{LAYERS[layer]}' and imports {forbidden}, "
        f"which is above it ({' < '.join(LAYERS)})"
    )
