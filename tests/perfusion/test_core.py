"""Tests of the cellular core of a zone of a perfused organ."""

import dataclasses
import math
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


def _replace(params: core.CoreParameters, **values: float) -> core.CoreParameters:
    """The constants with other numbers of some of them."""
    changes = {
        k: dataclasses.replace(getattr(params, k), value=x) for k, x in values.items()
    }
    return dataclasses.replace(params, **changes)


def _starts(params: core.CoreParameters, **values: float) -> core.CoreParameters:
    """The constants with other starts of the cells."""
    starts = dict(params.starts)
    for q, x in values.items():
        starts[q] = dataclasses.replace(starts[q], value=x)
    return dataclasses.replace(params, starts=FrozenMap(starts))


def _composition(**values: float) -> core.Composition:
    """The default perfusate with other numbers of some quantities."""
    c = core.Composition.default()
    changes = {
        k: dataclasses.replace(getattr(c, k), value=x) for k, x in values.items()
    }
    return dataclasses.replace(c, **changes)


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


def test_the_zones_share_the_constants(tmp_path: Path) -> None:
    """Two zones read one parameter per constant, the profile weighs a
    process per zone."""
    result = zone_model(
        tmp_path / "zones.xml",
        zones=("z1", "z2"),
        profiles={"z1": {"OXPHOS": 1.5}, "z2": {"OXPHOS": 0.5}},
    )
    assert result.validation is not None
    assert result.validation.error_count == 0
    assert result.validation.warning_count == 0
    r = load_zone(tmp_path / "zones.xml")
    assert r["GLC2PYR_Vmax"] == 0.2
    assert r[ids.profile("OXPHOS", "z1")] == 1.5
    assert r["OXPHOS_z1"] == pytest.approx(3 * r["OXPHOS_z2"], rel=1e-12)
    assert r["GLC2PYR_z1"] == pytest.approx(r["GLC2PYR_z2"], rel=1e-12)


def test_a_profile_names_processes_of_the_core() -> None:
    params = core.CoreParameters.default()
    with pytest.raises(ValueError, match="No processes of the core"):
        core.zone_objects(ZONE, CELLS, VASCULAR, params, profile={"NEC": 1.0})
    with pytest.raises(ValueError, match="not >= 0"):
        core.zone_objects(ZONE, CELLS, VASCULAR, params, profile={"TCA": -1.0})


def test_every_element_has_a_name_and_an_sbo_term() -> None:
    params = core.CoreParameters.default()
    objects = core.zone_objects(ZONE, CELLS, VASCULAR, params, profile={"OXPHOS": 1.0})
    elements = [
        *objects.species,
        *objects.reactions,
        *objects.parameters,
        *core.core_parameters(params),
        *core.buffer_parameters(),
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
    for q in ids.RELEASABLE:
        assert ids.species(q, DEAD) in sids
    assert not objects.parameters[0].sid.startswith("GLCIM")
    assert ids.zone_id(3) == "z3"
    assert ids.species("o2", "res") == "o2_res"
    with pytest.raises(ValueError, match="numbered from 1"):
        ids.zone_id(0)


def test_every_constant_is_a_value_of_the_database() -> None:
    """Every constant is a value of the perfusion table with a source; the
    module holds no numbers of its own."""
    params = core.CoreParameters.default()
    for field in dataclasses.fields(params):
        v = getattr(params, field.name)
        values = v.values() if isinstance(v, FrozenMap) else [v]
        for x in values:
            assert isinstance(x, Value), field.name
            assert x.source is not None, field.name
            if x.source == "ModelChoice":
                assert x.notes, field.name
    assert params.GLC2PYR_Vmax.notes is not None
    assert "basis changed" in params.GLC2PYR_Vmax.notes
    assert params.LDH_Keq.notes is not None
    assert "1e4" in params.LDH_Keq.notes


def test_no_rate_divides_by_a_species(tmp_path: Path) -> None:
    """No kinetic law or rule has a species in a denominator that can reach 0."""
    zone_model(tmp_path / "zone.xml")
    doc = libsbml.readSBMLFromFile(str(tmp_path / "zone.xml"))
    model = doc.getModel()
    species = _species_ids(model)
    maths = {
        model.getReaction(k).getId(): model.getReaction(k).getKineticLaw().getMath()
        for k in range(model.getNumReactions())
    }
    maths |= {
        f"rule {model.getRule(k).getVariable()}": model.getRule(k).getMath()
        for k in range(model.getNumRules())
    }
    assert any(key.startswith("rule ") for key in maths)
    problems = {
        key: found
        for key, math_ in maths.items()
        if (found := _unsafe_divisions(math_, species))
    }
    assert not problems


def test_hill_terms_are_defined_below_zero(tmp_path: Path) -> None:
    """A Hill coefficient of 2.5 and a species slightly below 0 give finite
    rates (`power(abs(S)/K, n)`)."""
    params = _replace(
        core.CoreParameters.default(),
        GLC2PYR_n_atp=2.5,
        ROSP_n_o2=2.5,
        NEC_n_atp=2.5,
        NEC_n_h=2.5,
    )
    zone_model(tmp_path / "zone.xml", params)
    r = load_zone(tmp_path / "zone.xml")
    for q in ("atp", "o2"):
        r[f"[{ids.species(q, CELL)}]"] = -1e-12
    r[f"[{ids.species('h', ZONE)}]"] = -1e-12
    rates = r.model.getReactionRates()
    assert all(math.isfinite(x) for x in rates)
    assert math.isfinite(r[ids.necrosis_rate(ZONE)])
    assert r[ids.reaction("GLC2PYR", ZONE)] > 0


def test_dead_cells_release_their_enzymes(tmp_path: Path) -> None:
    """Without oxygen ATP falls to 0 (the anaerobic glycolysis of the
    default constants is limited by the NAD and cannot cover the ATP use),
    the cells die, their enzymes and potassium move into the vascular space,
    the uptake of glucose and oxygen stops."""
    params = _starts(core.CoreParameters.default(), o2=0.0)
    zone_model(tmp_path / "zone.xml", params, _composition(o2=0.0))
    r = load_zone(tmp_path / "zone.xml")
    cells, vascular = r[CELLS], r[VASCULAR]
    contents = {q: r[f"[{ids.species(q, CELL)}]"] * cells for q in ids.RELEASABLE}
    before = {q: r[f"[{ids.species(q, ZONE)}]"] * vascular for q in ids.RELEASABLE}
    glcim = ids.reaction("GLCIM", ZONE)
    o2im = ids.reaction("O2IM", ZONE)
    r.timeCourseSelections = [
        "time",
        ids.alive(ZONE),
        f"[{ids.species('atp', CELL)}]",
        glcim,
        o2im,
        *[f"[{ids.species(q, ZONE)}]" for q in ids.RELEASABLE],
    ]
    s = r.simulate(0, 3000, 301)
    assert s[ids.alive(ZONE)][-1] < 1e-6
    assert s[f"[{ids.species('atp', CELL)}]"][-1] < 0.01
    assert s[glcim][0] > 0
    for q in ids.RELEASABLE:
        after = s[f"[{ids.species(q, ZONE)}]"][-1] * vascular
        assert after == pytest.approx(before[q] + contents[q], rel=1e-6), q
    assert abs(s[glcim][-1]) < 1e-9
    assert abs(s[o2im][-1]) < 1e-9


def test_the_buffer_follows_henderson_hasselbalch(tmp_path: Path) -> None:
    """At the buffer equilibrium pH = pK + log10(HCO3/CO2), at the start
    (the initial assignment of the CO2) and while the cells acidify."""
    zone_model(tmp_path / "zone.xml")
    r = load_zone(tmp_path / "zone.xml")
    pk = 3 - math.log10(r[ids.BUFFER_KA])
    assert pk == pytest.approx(6.1, rel=1e-12)

    def henderson_hasselbalch() -> float:
        hco3 = r[f"[{ids.species('hco3', ZONE)}]"]
        co2 = r[f"[{ids.species('co2', ZONE)}]"]
        return pk + math.log10(hco3 / co2)

    assert r[ids.ph(ZONE)] == pytest.approx(7.4, abs=1e-12)
    assert henderson_hasselbalch() == pytest.approx(7.4, abs=1e-12)
    r.simulate(0, 300, 31)
    assert r[ids.ph(ZONE)] < 7.35
    # the CO2 the cells produce keeps the buffer slightly behind its equilibrium
    assert r[ids.ph(ZONE)] == pytest.approx(henderson_hasselbalch(), abs=1e-4)


def test_lactate_export_lowers_the_bicarbonate(tmp_path: Path) -> None:
    """A lactate exported with its proton consumes one bicarbonate:
    HCO3 + lac - h of the perfusate is conserved."""
    params = _starts(core.CoreParameters.default(), lac=20.0)
    zone_model(tmp_path / "zone.xml", params)
    r = load_zone(tmp_path / "zone.xml")

    def c(q: str) -> float:
        return r[f"[{ids.species(q, ZONE)}]"]

    lac0, hco30, h0 = c("lac"), c("hco3"), c("h")
    r.simulate(0, 60, 61)
    assert c("lac") - lac0 > 1.0
    assert c("hco3") - hco30 == pytest.approx(
        -(c("lac") - lac0) + (c("h") - h0), rel=1e-9
    )
    assert c("hco3") - hco30 == pytest.approx(-(c("lac") - lac0), rel=1e-3)


def test_acidosis_drives_necrosis(tmp_path: Path) -> None:
    """At the ATP of the start the cells die only in an acidic perfusate."""
    rates = {}
    for ph in (7.4, 6.8):
        path = tmp_path / f"zone_{ph}.xml"
        zone_model(path, composition=_composition(ph=ph))
        r = load_zone(path)
        assert r[f"[{ids.species('atp', CELL)}]"] == 3.0
        rates[ph] = r[ids.necrosis_rate(ZONE)] / r["NEC_k"]
    assert rates[7.4] < 1e-3
    assert rates[6.8] > 0.9


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
    for q in ids.RELEASABLE:
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
