"""The equations of the page `docs/perfusion.md` are those the library generates.

The page shows the equations of a machine perfusion in fences with the info
string `equations`, in the syntax of the PBPK pages
(`tests/docs/test_pbpk_equations.py`, whose parser and comparison this test
reuses), with the placeholders `<q>` for a quantity of the perfusate, `<z>`
for a zone, `<m>` for the submodel of the zone (`M01`, the cells of a zone are
the cell model of its submodel) and `<o>` for the token of the organ. The test builds a rat liver
and a rat kidney of two zones with an inert tracer `trc` and finds every
equation, by its id, in one of them with the same math.
"""

import itertools
from pathlib import Path

import libsbml
import pytest
from perfusion.conftest import kidney_definition, liver_definition

from pkpdsim.perfusion import circuit
from pkpdsim.perfusion.builder import build
from pkpdsim.perfusion.definition import PerfusionDefinition

from .test_pbpk_equations import EQUATION, _lookup, equations

PAGE = Path(__file__).parent.parent.parent / "docs" / "perfusion.md"
ZONES = ("z1", "z2")


@pytest.fixture(scope="module")
def models(
    tmp_path_factory: pytest.TempPathFactory,
) -> list[tuple[libsbml.SBMLDocument, PerfusionDefinition]]:
    """The liver and the kidney of the page, two zones each (the documents keep
    their models alive)."""
    out = tmp_path_factory.mktemp("perfusion_equations")
    result = []
    for definition in (
        liver_definition(zones=2, tracers=("trc",), sid="doc_liver"),
        kidney_definition(zones=2, tracers=("trc",), sid="doc_kidney"),
    ):
        doc = libsbml.readSBMLFromFile(str(build(definition, out)))
        result.append((doc, definition))
    return result


def test_the_page_has_equations() -> None:
    assert len(equations(PAGE)) > 30


def test_every_equation_is_generated(
    models: list[tuple[libsbml.SBMLDocument, PerfusionDefinition]],
) -> None:
    problems: list[str] = []
    for line in equations(PAGE):
        match = EQUATION.match(line)
        if match is None:
            problems.append(f"not an equation: {line}")
            continue
        found: list[str] = []
        for doc, definition in models:
            model = doc.getModel()
            quantities = circuit.quantities(definition.perfusate, definition.released)
            organ = definition.organ.key
            for q, z in itertools.product(quantities, ZONES):
                difference = _lookup(
                    model,
                    match,
                    {"<q>": q, "<z>": z, "<o>": organ, "<m>": f"M{int(z[1:]):02}"},
                )
                if difference == "":
                    break
                if difference != "missing":
                    found.append(difference)
            else:
                continue
            break
        else:
            closest = f"\n    generated: {found[0]}" if found else " (no such id)"
            problems.append(f"{line}{closest}")
    assert not problems, "\n  ".join(["equations which no model generates:", *problems])
