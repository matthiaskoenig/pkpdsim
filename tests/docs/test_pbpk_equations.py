"""The equations of the PBPK pages are those the library generates.

A page of `docs/pbpk/` shows the equations of the library in fences with the
info string `equations` (rendered as text). Every line of such a fence is a
comment (`# ...`) or an equation in the infix syntax of SBML L3:

- `<id> = <formula>`: an assignment rule, an initial assignment or the rate
  of a reaction;
- `<id>: <reactants> -> <products> = <formula>`: a reaction, its species
  (joined by ` + `, empty for a source or a sink) and its rate;
- `d(<id>)/dt = <formula>`: a rate rule.

The ids are the local ids of an organ component (`Cpl_<s>`, `Q<o>`) or the
ids of the top model (`Carm_<s>`), with the placeholders `<s>` for a
substance and `<o>` for the token of an organ (`li`). The test builds the
toy definitions below and finds every equation, by its id, in the comp model
of a toy: in an organ component (its model definition, `<o>` the organ of
the component) or in the top model (`<o>` any organ), for some substance of
the toy. The formulas are compared as math: both are parsed, the units of
their numbers are dropped, and the printed formulas must be equal.
"""

import itertools
import re
from dataclasses import dataclass
from pathlib import Path

import libsbml
import pytest
from sbmlutils.factory import AssignmentRule, Model, Parameter, Port, PortType

from pkpdsim.pbpk import (
    MM,
    MODEL_UNITS,
    REST,
    Binding,
    Direction,
    Erythrocytes,
    Excretion,
    FirstOrder,
    FirstOrderDissolution,
    Hill,
    Induction,
    Inhibition,
    Kp,
    Membrane,
    Metabolism,
    Mode,
    Organ,
    PBPKDefinition,
    PDComponent,
    RateOf,
    Route,
    SegmentedIntestine,
    SlowPool,
    Substance,
    Transport,
    U,
    Value,
    Variants,
    WeibullDissolution,
    ZonatedLiver,
    build,
)

#: the pages of the PBPK library
PBPK_DIR = Path(__file__).parent.parent.parent / "docs" / "pbpk"

#: a fence of equations
FENCE = re.compile(r"^```equations[^\n]*\n(.*?)^```", re.S | re.M)

#: an equation: a reaction with its species, a rate rule or any other math
EQUATION = re.compile(
    r"^(?:d\((?P<rate>[\w<>]+)\)/dt"
    r"|(?P<id>[\w<>]+)(?::\s*(?P<reactants>.*?)\s*->\s*(?P<products>.*?))?)"
    r"\s+=\s+(?P<math>.+)$"
)


def standard() -> PBPKDefinition:
    """Get the toy with the standard organs and every collapsing mode.

    `drg` (merged tissue, urine), its metabolite `drg_m` (vascular), `prf`
    (flow-limited, erythrocytes in equilibrium, a Kp override of adipose),
    `prm` (full form, kinetic erythrocytes, transport), `hil` (Hill kinetics
    in the kidney) and `inh` (a mechanism-based inhibitor and inducer of the
    CYP2D6 which metabolises it), `spl` (split by a fraction into `spl_m1`
    and the remainder `spl_m2`); the OATP1B3 uptake `PRMIN3` of `prm` has a
    capacity as a ratio of `PRMIN`; `drg` inhibits CYP2C9 competitively;
    `drg`, `prf` and `prm` have an endogenous baseline; a PD component reads
    the liver cells of `prm`; `drg` has a slow pool in muscle.
    """
    return PBPKDefinition(
        sid="toy_standard",
        substances=[
            Substance(
                "drg",
                name="drug",
                mr=300.0,
                binding=Binding(fu=0.1, protein="albumin"),
                distribution=Kp(2.0, ftissue=1.0),
                excretion=Excretion(urine=True, fraction_reabsorbed=0.1),
                inhibits=(Inhibition("CYP2C9", "competitive", ki=0.01),),
                baseline=0.01,
            ),
            Substance("drg_m", name="metabolite", mr=316.0),
            Substance(
                "prf",
                name="flow-limited drug",
                mr=250.0,
                binding=Binding(fu=0.5),
                erythrocytes=Erythrocytes("equilibrium", bp=0.8),
                distribution=Kp(
                    3.0, mode=Mode.PERFUSION, overrides={Organ.ADIPOSE: 10.0}
                ),
                baseline=0.02,
            ),
            Substance(
                "prm",
                name="permeability-limited drug",
                mr=400.0,
                binding=Binding(fu=0.2),
                erythrocytes=Erythrocytes("kinetic", kbc=1.5, k_exchange=0.1),
                distribution=Kp(4.0, mode=Mode.PERMEABILITY, ftissue=0.5),
                baseline=0.03,
            ),
            Substance("hil", name="drug", mr=200.0, distribution=Kp(1.5, ftissue=1.0)),
            Substance(
                "inh",
                name="inhibitor",
                mr=350.0,
                binding=Binding(fu=0.3),
                distribution=Kp(1.0, ftissue=1.0),
                inhibits=(Inhibition("CYP2D6", "mechanism", ki=0.005, kinact=0.05),),
                induces=(Induction("CYP2D6", emax=2.0, ec50=0.1),),
            ),
            Substance(
                "spl",
                name="drug split into two metabolites",
                mr=300.0,
                binding=Binding(fu=0.5),
                distribution=Kp(1.0, ftissue=1.0),
            ),
            Substance("spl_m1", name="metabolite 1", mr=316.0),
            Substance("spl_m2", name="metabolite 2", mr=332.0),
        ],
        chemistry=[
            Metabolism(
                "DRG2M",
                Organ.LIVER,
                "drg",
                {"drg_m": 1.0},
                MM(vmax=1e-3, km=0.01),
                protein="CYP3A4",
            ),
            Metabolism(
                "PRFEL",
                Organ.LIVER,
                "prf",
                {},
                FirstOrder(k=0.2),
                protein="CYP2C9",
            ),
            Metabolism(
                "HILEL", Organ.KIDNEY, "hil", {}, Hill(vmax=1e-3, km=0.01, n=2.0)
            ),
            Metabolism(
                "INHEL",
                Organ.LIVER,
                "inh",
                {},
                MM(vmax=1e-3, km=0.02),
                protein="CYP2D6",
            ),
            Transport(
                "PRMIN",
                Organ.LIVER,
                "prm",
                Direction.INFLUX,
                Membrane.BASOLATERAL,
                MM(vmax=1e-2, km=0.05),
                protein="SLCO1B1",
            ),
            Transport(
                "PRMBIL",
                Organ.LIVER,
                "prm",
                Direction.EFFLUX,
                Membrane.APICAL,
                FirstOrder(k=0.1),
                protein="ABCC2",
            ),
            Transport(
                "PRMIN3",
                Organ.LIVER,
                "prm",
                Direction.INFLUX,
                Membrane.BASOLATERAL,
                MM(vmax=RateOf("PRMIN", 0.5), km=0.02),
                protein="SLCO1B3",
            ),
            Metabolism(
                "SPL2M",
                Organ.LIVER,
                "spl",
                {"spl_m1": Value(0.4, U.dimensionless), "spl_m2": REST},
                FirstOrder(k=0.3),
            ),
        ],
        routes=[Route.iv("drg"), Route.oral("drg"), Route.oral_solution("prf")],
        pd=[_cell_reader()],
        slow_pools=[SlowPool("drg", frozenset({Organ.MUSCLE}), k=1e-3, kp=1.5)],
    )


def _cell_reader() -> PDComponent:
    """Get a PD component which reads the concentration in the liver cells of `prm`."""
    model = Model(
        "cells",
        units=U,
        model_units=MODEL_UNITS,
        objects=[
            Parameter(
                "C_cell",
                0.0,
                U.mM,
                constant=False,
                port=Port("C_cell_port", portType=PortType.INPUT_PORT),
            ),
            Parameter("E", None, U.mM, constant=False),
            AssignmentRule("E", "C_cell", U.mM),
        ],
    )
    return PDComponent(model=model, inputs={"C_cell": "LI__Cic_prm"}, outputs={})


def variants() -> PBPKDefinition:
    """Get the toy with the zonated liver, the nephron kidney and a gut metabolism."""
    return PBPKDefinition(
        sid="toy_variants",
        substances=[
            Substance(
                "drg",
                name="drug",
                mr=300.0,
                binding=Binding(fu=0.1, protein="albumin"),
                distribution=Kp(2.0, ftissue=1.0),
                excretion=Excretion(urine=True, fraction_reabsorbed=0.1),
            ),
            Substance("drg_m", name="metabolite", mr=316.0),
        ],
        chemistry=[
            Metabolism(
                "DRG2M",
                Organ.LIVER,
                "drg",
                {"drg_m": 1.0},
                MM(vmax=1e-3, km=0.01),
                protein="CYP3A4",
            ),
            Metabolism(
                "DRGGU",
                Organ.SMALL_INTESTINE,
                "drg",
                {"drg_m": 1.0},
                MM(vmax=1e-3, km=0.01),
                protein="CYP3A4",
            ),
        ],
        variants=Variants(
            organs={
                Organ.LIVER: ZonatedLiver(
                    n_zones=3, profiles={"CYP3A4": (1.0, 2.0, 4.0)}
                ),
                Organ.KIDNEY: "nephron",
            }
        ),
        routes=[
            Route.iv("drg"),
            Route.oral("drg", WeibullDissolution(t50=30.0, shape=1.5)),
        ],
    )


def absorption() -> PBPKDefinition:
    """Get the toy with the fraction absorbed and the colon factor of the segmented intestine."""
    return PBPKDefinition(
        sid="toy_absorption",
        substances=[Substance("drg", name="drug", mr=300.0)],
        variants=Variants(
            organs={Organ.SMALL_INTESTINE: SegmentedIntestine(fabs=0.8, f_ka_colon=0.5)}
        ),
        routes=[
            Route.oral_solution("drg"),
            Route.oral("drg", FirstOrderDissolution(k=2.0, fed=0.5, release=0.2)),
        ],
    )


@dataclass(frozen=True)
class Context:
    """A model of a toy in which an equation is looked up.

    Attributes:
        document: the comp model of the toy, which owns the model.
        model: the top model or the model definition of an organ component.
        organs: the organ tokens `<o>` stands for.
        substances: the substance ids `<s>` stands for.
    """

    document: libsbml.SBMLDocument
    model: libsbml.Model
    organs: tuple[str, ...]
    substances: tuple[str, ...]


@pytest.fixture(scope="module")
def contexts(tmp_path_factory: pytest.TempPathFactory) -> list[Context]:
    """The models of the built toys: the top models and the organ components."""
    directory = tmp_path_factory.mktemp("pbpk_equations")
    result: list[Context] = []
    for definition in (standard(), variants(), absorption()):
        paths = build(definition, directory)
        doc = libsbml.readSBMLFromFile(str(paths.comp_path))
        substances = tuple(s.sid for s in definition.substances)
        organs = tuple(o.value for o in Organ)
        result.append(Context(doc, doc.getModel(), organs, substances))
        for md in doc.getPlugin("comp").getListOfModelDefinitions():
            organ = md.getId().rsplit("_", 1)[-1]
            result.append(Context(doc, md, (organ,), substances))
    return result


def _normalised(formula: str | libsbml.ASTNode) -> str:
    """Get a formula printed from its math, without the units of its numbers."""
    ast = libsbml.parseL3Formula(formula) if isinstance(formula, str) else formula
    if ast is None:
        raise ValueError(f"Cannot parse '{formula}': {libsbml.getLastParseL3Error()}")
    ast = ast.deepCopy()
    nodes = [ast]
    while nodes:
        node = nodes.pop()
        if node.isNumber():
            node.unsetUnits()
        nodes.extend(node.getChild(i) for i in range(node.getNumChildren()))
    return libsbml.formulaToL3String(ast)


def _species(text: str) -> list[str]:
    """Get the species of one side of a reaction equation."""
    return sorted(s.strip() for s in text.split("+") if s.strip())


def _lookup(model: libsbml.Model, match: re.Match[str], sub: dict[str, str]) -> str:
    """Compare an equation with a model.

    Args:
        model: the model.
        match: the equation.
        sub: the values of the placeholders.

    Returns:
        "" if the model has the element of the equation with the same math,
        "missing" if it has no element of the id, else what differs.
    """

    def fill(text: str) -> str:
        for key, value in sub.items():
            text = text.replace(key, value)
        return text

    expected = _normalised(fill(match["math"]))
    if match["rate"]:
        sid = fill(match["rate"])
        rule = model.getRateRule(sid)
        if rule is None:
            return "missing"
        actual = _normalised(rule.getMath())
        return "" if actual == expected else f"d({sid})/dt = {actual}"
    sid = fill(match["id"])
    reaction = model.getReaction(sid)
    if reaction is not None:
        actual = _normalised(reaction.getKineticLaw().getMath())
        if match["reactants"] is not None:
            reactants = sorted(r.getSpecies() for r in reaction.getListOfReactants())
            products = sorted(p.getSpecies() for p in reaction.getListOfProducts())
            if (reactants, products) != (
                _species(fill(match["reactants"])),
                _species(fill(match["products"])),
            ):
                return (
                    f"{sid}: {' + '.join(reactants)} -> {' + '.join(products)}"
                    f" = {actual}"
                )
        return "" if actual == expected else f"{sid} = {actual}"
    if match["reactants"] is not None:
        return "missing"
    element = model.getAssignmentRule(sid) or model.getInitialAssignment(sid)
    if element is None:
        return "missing"
    actual = _normalised(element.getMath())
    return "" if actual == expected else f"{sid} = {actual}"


def equations(page: Path) -> list[str]:
    """Get the equations of a page, the lines of its fences without comments."""
    return [
        line.strip()
        for block in FENCE.findall(page.read_text(encoding="utf-8"))
        for line in block.splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


#: the pages with equations
PAGES: list[str] = sorted(p.name for p in PBPK_DIR.glob("*.md") if equations(p))


def test_the_pages_with_equations_are_found() -> None:
    """The pages of the organ and the chemistry show equations."""
    assert {"organ.md", "chemistry.md"} <= set(PAGES)


def _problems(name: str, contexts: list[Context]) -> list[str]:
    """Get the equations of a page which no toy generates, with the closest match."""
    problems: list[str] = []
    for line in equations(PBPK_DIR / name):
        match = EQUATION.match(line)
        if match is None:
            problems.append(f"not an equation: {line}")
            continue
        found: list[str] = []
        for context in contexts:
            for s, o in itertools.product(context.substances, context.organs):
                difference = _lookup(context.model, match, {"<s>": s, "<o>": o})
                if difference == "":
                    break
                if difference != "missing":
                    found.append(difference)
            else:
                continue
            break
        else:
            closest = f"\n    generated: {found[0]}" if found else " (no such id)"
            problems.append(f"{name}: {line}{closest}")
    return problems


def test_every_equation_is_generated(contexts: list[Context]) -> None:
    """Every equation of every page is an element of a toy, with the same math.

    One test for all pages, so the toys are built once (the build is the cost
    of the test, and builds in several workers slow the parallel test run).
    """
    problems = [p for name in PAGES for p in _problems(name, contexts)]
    assert not problems, "\n  ".join(["equations which no toy generates:", *problems])
