"""Tests of the kidney and the generic organ of a machine perfusion."""

import dataclasses
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest
from sbmlutils.validation import ValidationOptions, validate_doc

from pkpdsim.database import tables
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion import core, ids
from pkpdsim.perfusion.builder import build
from pkpdsim.perfusion.definition import (
    GenericOrgan,
    KidneyOrgan,
    LiverOrgan,
    PerfusionDefinition,
)
from pkpdsim.perfusion.units import PU
from pkpdsim.species import Species

from .conftest import (
    generic_definition,
    kidney_definition,
    liver_definition,
    load_circuit,
)

if TYPE_CHECKING:
    import roadrunner

TRACER = "inu"
KIDNEY = KidneyOrgan.key


def _ids(path: Path) -> set[str]:
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    elements = model.getListOfAllElements()
    return {e.getId() for e in (elements.get(k) for k in range(elements.getSize()))} - {
        ""
    }


def _zones(n: int) -> list[str]:
    return [ids.zone_id(k) for k in range(1, n + 1)]


def _perfusate_volume(r: "roadrunner.RoadRunner", n: int) -> float:
    """The perfusate in the reservoir, the tubing and the zones (l)."""
    sites = [ids.RESERVOIR, ids.INFLOW, *_zones(n), ids.OUTFLOW]
    return sum(r[ids.volume(s)] for s in sites)


@pytest.mark.parametrize(
    "make",
    [
        kidney_definition,
        lambda: kidney_definition(3, sid="kidney3"),
        generic_definition,
    ],
    ids=["kidney", "kidney-3-zones", "muscle"],
)
def test_an_organ_builds_and_simulates_6_h(
    tmp_path: Path, make: Callable[[], PerfusionDefinition]
) -> None:
    """The rat kidney and a generic organ validate with units and run 6 h."""
    definition = make()
    path = build(definition, tmp_path)
    doc = libsbml.readSBMLFromFile(str(path))
    validation = validate_doc(doc, options=ValidationOptions())
    assert validation.error_count == 0
    assert validation.warning_count == 0
    r = load_circuit(path)
    r.timeCourseSelections = [
        "time",
        *r.model.getFloatingSpeciesIds(),
        *r.model.getGlobalParameterIds(),
        *r.model.getReactionIds(),
    ]
    s = r.simulate(0, 360, 37)
    assert np.isfinite(np.asarray(s)).all()


def test_the_kidney_makes_urine_and_no_bile(tmp_path: Path) -> None:
    known = _ids(build(kidney_definition(), tmp_path))
    assert {ids.URINE, ids.URINE_RATE, ids.Q_URINE, ids.gfr(KIDNEY)} <= known
    assert ids.volume(ids.URINE_SITE) in known
    assert {ids.species(q, ids.URINE_SITE) for q in ids.FILTERED} <= known
    assert ids.BILE not in known
    assert ids.BILE_RATE not in known
    assert ids.bile("z1") not in known


def test_the_generic_organ_has_no_bile_and_no_urine(tmp_path: Path) -> None:
    known = _ids(build(generic_definition(), tmp_path))
    for sid in (ids.BILE, ids.BILE_RATE, ids.URINE, ids.URINE_RATE, ids.Q_URINE):
        assert sid not in known
    assert ids.mass("mu") in known
    assert ids.species("atp", ids.cell_site("z1")) in known


@pytest.mark.parametrize("zones", [1, 3])
def test_the_kidney_conserves_perfusate_water_plus_urine(
    tmp_path: Path, zones: int
) -> None:
    """The urine leaves the reservoir: perfusate and urine keep their volume."""
    r = load_circuit(build(kidney_definition(zones), tmp_path))
    start = _perfusate_volume(r, zones)
    assert r[ids.volume(ids.URINE_SITE)] == 0.0
    rate = r[ids.Q_URINE]
    assert rate > 0
    r.simulate(0, 360, 37)
    urine = r[ids.volume(ids.URINE_SITE)]
    assert urine == pytest.approx(360 * rate, rel=1e-9)
    assert _perfusate_volume(r, zones) + urine == pytest.approx(start, rel=1e-12)
    assert r[ids.URINE] == pytest.approx(urine * 1e6, rel=1e-12)
    assert r[ids.URINE_RATE] == pytest.approx(rate * 1e6, rel=1e-12)
    # the net urine flow is the filtrate which the tubules do not reabsorb
    f_reab = r[ids.reabsorption(KIDNEY)]
    assert rate == pytest.approx((1 - f_reab) * r[ids.gfr(KIDNEY)], rel=1e-12)
    assert 0 < f_reab < 1


def _with_loads(path: Path, quantities: tuple[str, ...], zones: int) -> Path:
    """Add the filtered load of quantities, the integral of GFR*C of the zones."""
    doc = libsbml.readSBMLFromFile(str(path))
    model = doc.getModel()
    for q in quantities:
        sid = f"load_{q}"
        p = model.createParameter()
        p.setId(sid)
        p.setValue(0.0)
        p.setConstant(False)
        p.setUnits("mmole")
        rule = model.createRateRule()
        rule.setVariable(sid)
        terms = " + ".join(
            f"{ids.gfr(KIDNEY)}/{zones} * {ids.species(q, z)}" for z in _zones(zones)
        )
        rule.setMath(libsbml.parseL3Formula(terms))
    out = path.with_name(f"{path.stem}_loads.xml")
    libsbml.writeSBMLToFile(doc, str(out))
    return out


@pytest.mark.parametrize("zones", [1, 2])
def test_the_urine_of_a_filtered_marker_is_its_filtered_load(
    tmp_path: Path, zones: int
) -> None:
    """The inert marker is filtered and not reabsorbed; a small solute leaves
    with the urine water, its load times (1 - reabsorption)."""
    path = build(kidney_definition(zones, tracers=(TRACER,)), tmp_path)
    r = load_circuit(_with_loads(path, (TRACER, "k"), zones))
    r[f"[{ids.species(TRACER, ids.RESERVOIR)}]"] = 1.0
    r.simulate(0, 360, 37)
    f_reab = r[ids.reabsorption(KIDNEY)]
    marker = r[ids.species(TRACER, ids.URINE_SITE)]
    assert marker > 1e-3
    assert marker == pytest.approx(r[f"load_{TRACER}"], rel=1e-6)
    potassium = r[ids.species("k", ids.URINE_SITE)]
    assert potassium == pytest.approx((1 - f_reab) * r["load_k"], rel=1e-6)
    # the marker in the urine and the perfusate is the marker of the start
    total = sum(
        r[f"[{ids.species(TRACER, s)}]"] * r[ids.volume(s)]
        for s in [ids.RESERVOIR, ids.INFLOW, *_zones(zones), ids.OUTFLOW]
    )
    assert total + marker == pytest.approx(1.0 * 1.0, rel=1e-9)


def _finite_and_positive(r: "roadrunner.RoadRunner", end: float) -> None:
    """Simulate to `end` and check every value finite, no species below -1e-9."""
    species = r.model.getFloatingSpeciesIds()
    r.timeCourseSelections = [
        "time",
        *(f"[{s}]" if not s.endswith(f"_{ids.URINE_SITE}") else s for s in species),
        *r.model.getGlobalParameterIds(),
        *r.model.getReactionIds(),
    ]
    s = np.asarray(r.simulate(r.model.getTime(), end, 61))
    assert np.isfinite(s).all()
    lowest = s[:, 1 : 1 + len(species)].min(axis=0)
    assert not [sp for sp, x in zip(species, lowest, strict=True) if x < -1e-9]


@pytest.mark.parametrize("zones", [1, 3])
def test_zero_flow_stays_defined(tmp_path: Path, zones: int) -> None:
    """The kidney stops filtering when the pump stops: an hour at Q_pump = 0
    after an hour at the default flow."""
    r = load_circuit(build(kidney_definition(zones, tracers=(TRACER,)), tmp_path))
    r[f"[{ids.species(TRACER, ids.RESERVOIR)}]"] = 1.0
    start = _perfusate_volume(r, zones)
    r.simulate(0, 60, 7)
    assert r[ids.gfr(KIDNEY)] == pytest.approx(r[ids.gfr_max(KIDNEY)], rel=1e-12)
    r[ids.Q_PUMP] = 0.0
    _finite_and_positive(r, 120)
    assert r[ids.gfr(KIDNEY)] == 0.0
    assert r[ids.Q_URINE] == 0.0
    for site in (*_zones(zones), ids.OUTFLOW):
        assert r[ids.outflow(site)] == 0.0
    urine = r[ids.volume(ids.URINE_SITE)]
    assert _perfusate_volume(r, zones) + urine == pytest.approx(start, rel=1e-12)


@pytest.mark.parametrize("zones", [1, 3])
def test_a_low_flow_caps_the_filtration(tmp_path: Path, zones: int) -> None:
    """Below the flow of the uncapped urine the glomeruli filter the maximal
    filtration fraction of the flow, and the flow after the kidney stays
    positive."""
    r = load_circuit(build(kidney_definition(zones), tmp_path))
    start = _perfusate_volume(r, zones)
    uncapped = r[ids.Q_URINE]
    q = uncapped / 2
    r.simulate(0, 60, 7)
    r[ids.Q_PUMP] = q
    _finite_and_positive(r, 120)
    ff = r[ids.ff_max(KIDNEY)]
    f_reab = r[ids.reabsorption(KIDNEY)]
    assert r[ids.gfr(KIDNEY)] == pytest.approx(ff * q, rel=1e-12)
    assert r[ids.Q_URINE] == pytest.approx((1 - f_reab) * ff * q, rel=1e-12)
    assert r[ids.outflow(ids.OUTFLOW)] > 0
    urine = r[ids.volume(ids.URINE_SITE)]
    assert _perfusate_volume(r, zones) + urine == pytest.approx(start, rel=1e-12)


def test_the_kidney_zones_add_up_to_the_organ(tmp_path: Path) -> None:
    one = load_circuit(build(kidney_definition(1, sid="one"), tmp_path))
    many = load_circuit(build(kidney_definition(4, sid="many"), tmp_path))
    for sid in (ids.gfr(KIDNEY), ids.Q_URINE):
        assert many[sid] == pytest.approx(one[sid], rel=1e-12)
    cells = sum(many[ids.cells(z)] for z in _zones(4))
    assert cells == pytest.approx(one[ids.cells("z1")], rel=1e-12)


def test_the_released_set_is_that_of_the_organ(tmp_path: Path) -> None:
    """The kidney and the generic organ release LDH and potassium, the liver
    ALT, AST, LDH and potassium; the perfusate, the cells and the activities
    follow."""
    for definition, released in (
        (kidney_definition(), ("ldh", "k")),
        (generic_definition(), ("ldh", "k")),
        (liver_definition(), ("alt", "ast", "ldh", "k")),
    ):
        known = _ids(build(definition, tmp_path))
        assert definition.released == released
        for q in ids.RELEASABLE:
            present = q in released
            assert (ids.species(q, ids.dead_site("z1")) in known) is present, q
        for e in ids.ENZYMES:
            present = e in released
            for sid in (
                ids.species(e, ids.RESERVOIR),
                ids.species(e, ids.cell_site("z1")),
                ids.activity(e, ids.OUTFLOW),
                ids.specific_activity(e),
                ids.fresh(e),
            ):
                assert (sid in known) is present, sid
        # potassium stays in the perfusate and the cells
        assert ids.species("k", ids.RESERVOIR) in known
        assert ids.species("k", ids.cell_site("z1")) in known


def test_the_released_quantities_need_their_content() -> None:
    d = kidney_definition()
    organ = dataclasses.replace(d.organ, released=("alt", "ldh", "nad"))
    problems = "\n".join(dataclasses.replace(d, organ=organ).validate())
    assert "'nad'" in problems
    assert "ALT_content" in problems
    assert core.CoreParameters.default("kidney_rat").ALT_content is None
    assert core.CoreParameters.default("kidney_rat").BILE_Vmax is None
    # the liver needs its bile constants
    liver = liver_definition()
    no_bile = dataclasses.replace(liver.core, BILE_Vmax=None)
    problems = "\n".join(dataclasses.replace(liver, core=no_bile).validate())
    assert "BILE_Vmax" in problems


def test_the_validation_of_the_organs() -> None:
    d = kidney_definition()
    assert d.validate() == []
    assert generic_definition().validate() == []
    kidney = d.organ
    assert isinstance(kidney, KidneyOrgan)
    bad = dataclasses.replace(
        kidney,
        gfr=Value(1.0, PU.ml_per_min, "ModelChoice", "wrong unit"),
        reabsorption=Value(1.0, PU.dimensionless, "ModelChoice", "all"),
    )
    problems = "\n".join(dataclasses.replace(d, organ=bad).validate())
    assert "ml_per_min_g" in problems
    assert "reabsorption" in problems
    # a slow pump flow is no problem: the filtration follows it
    slow = dataclasses.replace(
        d.circuit, flow=Value(0.1, PU.ml_per_min, "ModelChoice", "slow")
    )
    assert dataclasses.replace(d, circuit=slow).validate() == []
    for ff in (0.0, 1.5):
        bad = dataclasses.replace(
            kidney, ff_max=Value(ff, PU.dimensionless, "ModelChoice", "test")
        )
        problems = "\n".join(dataclasses.replace(d, organ=bad).validate())
        assert "ff_max" in problems, ff
    # a generic organ: a token of the ids, the fractions of an organ of the
    # physiology table or given
    for key in ("res", "z1", "Mu", "1a"):
        organ = GenericOrgan(key, Value(1.0, PU.g, "ModelChoice", "test"))
        problems = "\n".join(
            dataclasses.replace(generic_definition(), organ=organ).validate()
        )
        assert f"'{key}'" in problems, key
    unknown = GenericOrgan("xx", Value(1.0, PU.g, "ModelChoice", "test"))
    problems = "\n".join(
        dataclasses.replace(generic_definition(), organ=unknown).validate()
    )
    assert "vascular_fraction" in problems
    given = GenericOrgan(
        "xx",
        Value(1.0, PU.g, "ModelChoice", "test"),
        vascular_fraction=Value(0.1, PU.dimensionless, "ModelChoice", "test"),
        cell_fraction=Value(0.7, PU.dimensionless, "ModelChoice", "test"),
        density=Value(1.0, PU.kg_per_l, "ModelChoice", "test"),
    )
    assert dataclasses.replace(generic_definition(), organ=given).validate() == []
    assert isinstance(liver_definition().organ, LiverOrgan)


def test_zones_and_profiles_follow_the_mass() -> None:
    """`zones` is positional after the mass (the token and the mass of a
    generic organ)."""
    mass = Value(1.0, PU.g, "ModelChoice", "test")
    assert LiverOrgan(mass, 3).zones == 3
    assert KidneyOrgan(mass, 2).zones == 2
    assert GenericOrgan("mu", mass, 4, {"OXPHOS": [1, 1, 1, 2]}).zones == 4


def test_the_generic_organ_reads_the_physiology_table() -> None:
    organ = GenericOrgan.default("mu")
    physiology = "physiology_rat"
    section = "male.organs.mu"
    mass = tables.entry(physiology, section, "mass").value
    f_vas = tables.entry(physiology, section, "f_vas").value
    f_int = tables.entry(physiology, section, "f_int").value
    assert isinstance(mass, float)
    assert isinstance(f_vas, float)
    assert isinstance(f_int, float)
    assert organ.mass.value == pytest.approx(mass * 1e3, rel=1e-12)
    resolved = organ.resolved(Species.RAT)
    assert resolved.vascular_fraction is not None
    assert resolved.cell_fraction is not None
    assert resolved.vascular_fraction.value == f_vas
    assert resolved.cell_fraction.value == pytest.approx(1 - f_vas - f_int, rel=1e-12)


@pytest.mark.parametrize("species", ["rat", "human", "pig"])
def test_the_kidney_of_every_species(species: str) -> None:
    """The kidney sections of the table: one kidney, the GFR per g of the
    physiology table; the default rat kidney filters below its maximal
    filtration fraction of the default flow."""
    organ = KidneyOrgan.default(species=Species(species))
    physiology = f"physiology_{species}"
    mass = tables.entry(physiology, "male.organs.ki", "mass").value
    gfr = tables.entry(physiology, "male", "gfr").value
    assert isinstance(mass, float)
    assert isinstance(gfr, float)
    assert organ.mass.value == pytest.approx(mass * 1e3 / 2, rel=1e-3)
    assert organ.gfr is not None
    assert organ.gfr.value == pytest.approx(gfr / (mass * 1e3), rel=1e-3)
    assert organ.problems() == []
    assert organ.ff_max is not None
    assert 0 < organ.ff_max.value <= 1
    if species == "rat":
        gfr_ml = organ.gfr.value * organ.mass.value
        flow = tables.entry("perfusion", "circuit", "flow").value
        assert isinstance(flow, float)
        assert gfr_ml < organ.ff_max.value * flow


def test_the_copied_core_constants_say_so() -> None:
    """The kidney and the generic organ have the constants of the liver,
    marked as such."""
    data = tables.load("perfusion")
    for section in ("core_kidney_rat", "core_generic"):
        for start in ("", "_start"):
            for key in data[section + start]:
                e = tables.entry("perfusion", section + start, key)
                if key == "glyc":
                    # the glycogen store is the liver's
                    assert e.value == 0.0
                    assert e.source == "ModelChoice"
                    continue
                liver = tables.entry("perfusion", "core_liver_rat" + start, key)
                assert e.value == liver.value, (section, key)
                assert e.source == "ModelChoice"
                assert "liver values, to be replaced" in (e.notes or "")
    for section in ("core_kidney_rat", "core_generic"):
        assert "ALT_content" not in data[section]
        assert "AST_content" not in data[section]
        assert "BILE_Km_atp" not in data[section]
        assert "LDH_content" in data[section]
