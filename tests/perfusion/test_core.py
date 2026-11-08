"""Tests of the cellular core of a zone of a perfused organ."""

import dataclasses
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import pytest

from pkpdsim.pbpk.values import FrozenMap, Value
from pkpdsim.perfusion import core, ids

from .conftest import CELLS, VASCULAR, ZONE, load_zone, zone_model

if TYPE_CHECKING:
    import roadrunner

CELL = ids.cell_site(ZONE)
DEAD = ids.dead_site(ZONE)


def _species_ids(model: libsbml.Model) -> set[str]:
    return {model.getSpecies(k).getId() for k in range(model.getNumSpecies())}


def _has_species(node: libsbml.ASTNode, species: set[str]) -> bool:
    if node.getType() == libsbml.AST_NAME and node.getName() in species:
        return True
    return any(
        _has_species(node.getChild(k), species) for k in range(node.getNumChildren())
    )


def _flat(node: libsbml.ASTNode, kind: int) -> list[libsbml.ASTNode]:
    """The terms of a sum or the factors of a product, nested ones flattened."""
    if node.getType() != kind:
        return [node]
    return [
        t for k in range(node.getNumChildren()) for t in _flat(node.getChild(k), kind)
    ]


def _safe_denominator(node: libsbml.ASTNode, species: set[str]) -> bool:
    """A denominator which cannot reach 0: free of species, a sum with a term
    free of species (`Km + S`, `1 + S/Km`), a product or power of such."""
    if not _has_species(node, species):
        return True
    t = node.getType()
    children = [node.getChild(k) for k in range(node.getNumChildren())]
    if t == libsbml.AST_PLUS:
        return any(not _has_species(c, species) for c in _flat(node, t))
    if t == libsbml.AST_TIMES:
        return all(_safe_denominator(c, species) for c in _flat(node, t))
    if t in (libsbml.AST_POWER, libsbml.AST_FUNCTION_POWER):
        return _safe_denominator(children[0], species)
    return False


def _unsafe_divisions(node: libsbml.ASTNode, species: set[str]) -> list[str]:
    found: list[str] = []
    if node.getType() == libsbml.AST_DIVIDE and not _safe_denominator(
        node.getChild(1), species
    ):
        found.append(libsbml.formulaToL3String(node))
    for k in range(node.getNumChildren()):
        found += _unsafe_divisions(node.getChild(k), species)
    return found


def test_zone_objects_build_a_valid_model(tmp_path: Path) -> None:
    """A one-zone model of the objects validates in sbmlutils with units."""
    result = zone_model(tmp_path / "zone.xml")
    assert result.validation is not None
    assert result.validation.error_count == 0
    assert result.validation.warning_count == 0


def test_every_element_has_a_name_and_an_sbo_term(tmp_path: Path) -> None:
    objects = core.zone_objects(ZONE, CELLS, VASCULAR, core.CoreParameters.default())
    elements = [
        *objects.species,
        *objects.reactions,
        *objects.parameters,
        *objects.rules,
        *objects.rate_rules,
        *objects.assignments,
    ]
    for e in elements:
        assert e.name, e
        assert e.sboTerm, e


def test_the_ids_follow_the_scheme() -> None:
    objects = core.zone_objects(ZONE, CELLS, VASCULAR, core.CoreParameters.default())
    sids = {s.sid for s in objects.species}
    for q in ids.PERFUSATE:
        assert ids.species(q, ZONE) in sids
    for q in ids.CELL:
        assert ids.species(q, CELL) in sids
    for q in ids.RELEASED:
        assert ids.species(q, DEAD) in sids
    assert ids.zone_id(3) == "z3"
    assert ids.species("o2", "res") == "o2_res"
    with pytest.raises(ValueError, match="numbered from 1"):
        ids.zone_id(0)


def test_every_constant_is_a_value_with_a_source() -> None:
    params = core.CoreParameters.default()
    for field in dataclasses.fields(params):
        v = getattr(params, field.name)
        values = v.values() if isinstance(v, FrozenMap) else [v]
        for x in values:
            assert isinstance(x, Value), field.name
            assert x.source is not None, field.name
            if x.source == "ModelChoice":
                assert x.notes, field.name


def test_no_rate_divides_by_a_species(tmp_path: Path) -> None:
    """No kinetic law has a species in a denominator that can reach 0."""
    zone_model(tmp_path / "zone.xml")
    doc = libsbml.readSBMLFromFile(str(tmp_path / "zone.xml"))
    model = doc.getModel()
    species = _species_ids(model)
    problems: dict[str, list[str]] = {}
    for k in range(model.getNumReactions()):
        reaction = model.getReaction(k)
        found = _unsafe_divisions(reaction.getKineticLaw().getMath(), species)
        if found:
            problems[reaction.getId()] = found
    assert not problems


def test_dead_cells_release_their_enzymes(tmp_path: Path) -> None:
    """Without oxygen and glycolytic substrate ATP falls to 0, the cells die,
    their enzymes and potassium move into the vascular space, the uptake of
    glucose and oxygen stops."""
    params = core.CoreParameters.default()
    starts = dict(params.starts)
    starts["o2"] = dataclasses.replace(starts["o2"], value=0.0)
    starts["glyc"] = dataclasses.replace(starts["glyc"], value=0.0)
    params = dataclasses.replace(
        params,
        GLC2PYR_Vmax=dataclasses.replace(params.GLC2PYR_Vmax, value=0.0),
        starts=FrozenMap(starts),
    )
    composition = dataclasses.replace(
        core.Composition.default(),
        o2=dataclasses.replace(core.Composition.default().o2, value=0.0),
    )
    zone_model(tmp_path / "zone.xml", params, composition)
    r = load_zone(tmp_path / "zone.xml")
    cells, vascular = r[CELLS], r[VASCULAR]
    contents = {q: r[f"[{ids.species(q, CELL)}]"] * cells for q in ids.RELEASED}
    before = {q: r[f"[{ids.species(q, ZONE)}]"] * vascular for q in ids.RELEASED}
    glcim = ids.reaction("GLCIM", ZONE)
    o2im = ids.reaction("O2IM", ZONE)
    r.timeCourseSelections = [
        "time",
        ids.alive(ZONE),
        f"[{ids.species('atp', CELL)}]",
        glcim,
        o2im,
        *[f"[{ids.species(q, ZONE)}]" for q in ids.RELEASED],
    ]
    s = r.simulate(0, 10000, 1001)
    assert s[f"[{ids.species('atp', CELL)}]"][-1] < 0.01
    assert s[ids.alive(ZONE)][-1] < 1e-6
    for q in ids.RELEASED:
        after = s[f"[{ids.species(q, ZONE)}]"][-1] * vascular
        assert after == pytest.approx(before[q] + contents[q], rel=1e-6), q
    assert abs(s[glcim][-1]) < 1e-9
    assert abs(s[o2im][-1]) < 1e-9


CARBON = {"glc": 6, "glyc": 6, "pyr": 3, "lac": 3, "co2": 1, "hco3": 1}


def _totals(
    r: "roadrunner.RoadRunner", cells: float, vascular: float
) -> dict[str, float]:
    def c(q: str, site: str) -> float:
        sid = ids.species(q, site)
        try:
            return r[f"[{sid}]"]
        except RuntimeError:
            return 0.0

    def amount(q: str) -> float:
        return c(q, CELL) * cells + c(q, DEAD) * cells + c(q, ZONE) * vascular

    totals = {
        "carbon": sum(n * amount(q) for q, n in CARBON.items()),
        "adenine": amount("atp") + amount("adp"),
        "nad": amount("nad") + amount("nadh"),
        "nadp": amount("nadp") + amount("nadph"),
        "glutathione": amount("gsh") + 2 * amount("gssg"),
    }
    for q in ids.RELEASED:
        totals[q] = amount(q)
    return totals


def test_atomic_balances(tmp_path: Path) -> None:
    """Carbon, the nucleotides, glutathione, potassium and the enzymes are
    conserved in a closed zone without bile, while the cells consume the
    oxygen, acidify the vascular space and die."""
    zone_model(tmp_path / "zone.xml", bile=False)
    r = load_zone(tmp_path / "zone.xml")
    cells, vascular = r[CELLS], r[VASCULAR]
    start = _totals(r, cells, vascular)
    r.simulate(0, 3000, 301)
    end = _totals(r, cells, vascular)
    assert r[ids.alive(ZONE)] < 0.99, "the zone should die in the closed space"
    for key, value in start.items():
        assert end[key] == pytest.approx(value, rel=1e-9), key


def test_the_check_of_the_denominators() -> None:
    """The check flags the singular ATP inhibition of the ported model and
    passes a Michaelis-Menten term and a product of sums."""
    species = {"atp", "s", "p"}

    def unsafe(formula: str) -> list[str]:
        return _unsafe_divisions(libsbml.parseL3Formula(formula), species)

    assert unsafe("Vmax * power(Ki/atp, n)")
    assert unsafe("s/(s - Km)")
    assert not unsafe("Vmax * s/(Km + s)")
    assert not unsafe("V/(Ks*Kp) * (s*p)/((1 + s/Ks + p/Kp) * (1 + p/Kp))")
    assert not unsafe("1/(1 + power(atp/Ki, n))")


def test_a_zone_id_is_checked() -> None:
    with pytest.raises(ValueError, match="no zone id"):
        core.zone_objects("li", CELLS, VASCULAR, core.CoreParameters.default())


def test_importing_the_library_is_light() -> None:
    """`import pkpdsim.perfusion` imports neither sbmlutils nor the PBPK library."""
    code = (
        "import sys, pkpdsim.perfusion as p; "
        "assert 'sbmlutils' not in sys.modules; "
        "assert 'pkpdsim.pbpk' not in sys.modules; "
        "assert callable(p.zone_objects)"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
