"""Tests of the definition and the build of a machine perfusion circuit."""

import dataclasses
import math
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
import numpy as np
import pytest
from pbpk.conftest import individual
from sbmlutils.validation import ValidationOptions, validate_doc

from pkpdsim.database.tables import entry
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion import core, ids
from pkpdsim.perfusion.builder import build
from pkpdsim.perfusion.definition import LiverOrgan, Perfusate
from pkpdsim.perfusion.units import PU

from .conftest import liver_definition, load_circuit

if TYPE_CHECKING:
    import roadrunner

TRACER = "tr"


def _sites(n: int) -> list[str]:
    """The sites of the perfusate of a circuit with `n` zones."""
    zones = [ids.zone_id(k) for k in range(1, n + 1)]
    return [ids.RESERVOIR, ids.INFLOW, ids.OUTFLOW, *zones]


def _model(path: Path) -> libsbml.Model:
    return libsbml.readSBMLFromFile(str(path)).getModel()


def _ids(model: libsbml.Model) -> set[str]:
    elements = model.getListOfAllElements()
    return {e.getId() for e in (elements.get(k) for k in range(elements.getSize()))} - {
        ""
    }


def _label(e: libsbml.SBase) -> str:
    if isinstance(e, libsbml.Rule):
        return e.getVariable()
    if isinstance(e, libsbml.InitialAssignment):
        return e.getSymbol()
    return e.getId()


@pytest.mark.parametrize("zones", [1, 8])
def test_build_writes_a_valid_model(tmp_path: Path, zones: int) -> None:
    """A rat liver in the circuit validates with units and loads at O0."""
    path = build(liver_definition(zones), tmp_path)
    assert path == tmp_path / "perfusion_liver.xml"
    doc = libsbml.readSBMLFromFile(str(path))
    assert (doc.getLevel(), doc.getVersion()) == (3, 2)
    assert doc.getModel().getPlugin("comp") is None
    validation = validate_doc(doc, options=ValidationOptions())
    assert validation.error_count == 0
    assert validation.warning_count == 0
    r = load_circuit(path)
    r.simulate(0, 10, 11)


def test_every_element_has_a_name_and_an_sbo_term(tmp_path: Path) -> None:
    model = _model(build(liver_definition(3), tmp_path))
    groups = {
        "compartment": model.getListOfCompartments(),
        "species": model.getListOfSpecies(),
        "parameter": model.getListOfParameters(),
        "reaction": model.getListOfReactions(),
        "rule": model.getListOfRules(),
        "initial assignment": model.getListOfInitialAssignments(),
    }
    missing = [
        (kind, _label(e))
        for kind, elements in groups.items()
        for e in elements
        if not e.getName() or not e.isSetSBOTerm()
    ]
    assert not missing


@pytest.mark.parametrize("n", [1, 4, 10])
def test_zones_add_up_to_the_organ(tmp_path: Path, n: int) -> None:
    """The cell and vascular volume, the enzyme content and the capacity of
    every process of `n` uniform zones are those of one zone."""
    one = load_circuit(build(liver_definition(1, sid="one"), tmp_path))
    many = load_circuit(build(liver_definition(n, sid="many"), tmp_path))
    zones = [ids.zone_id(k) for k in range(1, n + 1)]

    def total(
        r: "roadrunner.RoadRunner", make: Callable[[str], str], zs: list[str]
    ) -> float:
        return sum(r[make(z)] for z in zs)

    for volume in (ids.cells, ids.vascular):
        assert total(many, volume, zones) == pytest.approx(one[volume("z1")], rel=1e-12)
    for q in ("alt", "ast", "ldh", "k"):

        def amount(z: str, q: str = q) -> str:
            return f"{ids.species(q, ids.cell_site(z))}"

        assert total(many, amount, zones) == pytest.approx(
            one[amount("z1")], rel=1e-12
        ), q
    for process in core.PROCESSES:
        if process == "BILE":
            continue

        def rate(z: str, p: str = process) -> str:
            return ids.reaction(p, z)

        assert total(many, rate, zones) == pytest.approx(one[rate("z1")], rel=1e-12), (
            process
        )
    assert many[ids.BILE_RATE] == pytest.approx(one[ids.BILE_RATE], rel=1e-12)
    assert one[ids.BILE_RATE] > 0


def test_a_profile_keeps_the_capacity_of_the_organ(tmp_path: Path) -> None:
    """A periportal gradient moves the oxygen consumption between the zones
    and keeps the total."""
    one = load_circuit(build(liver_definition(1, sid="one"), tmp_path))
    graded = load_circuit(
        build(
            liver_definition(4, sid="graded", profiles={"OXPHOS": [4, 3, 2, 1]}),
            tmp_path,
        )
    )
    rates = [graded[ids.reaction("OXPHOS", ids.zone_id(k))] for k in range(1, 5)]
    assert sum(rates) == pytest.approx(one[ids.reaction("OXPHOS", "z1")], rel=1e-12)
    assert rates[0] == pytest.approx(4 * rates[3], rel=1e-12)
    assert graded[ids.profile("OXPHOS", "z1")] == pytest.approx(1.6, rel=1e-12)


@pytest.mark.parametrize("n", [1, 3, 10])
def test_mass_of_the_perfusate_is_conserved(tmp_path: Path, n: int) -> None:
    """An inert tracer in the perfusate keeps its amount over 6 h."""
    path = build(liver_definition(n, tracers=(TRACER,)), tmp_path)
    r = load_circuit(path)
    r[f"[{ids.species(TRACER, ids.RESERVOIR)}]"] = 1.0

    def amount() -> float:
        return sum(
            r[f"[{ids.species(TRACER, site)}]"] * r[ids.volume(site)]
            for site in _sites(n)
        )

    start = amount()
    assert start == pytest.approx(r[ids.volume(ids.RESERVOIR)], rel=1e-12)
    r.simulate(0, 360, 37)
    assert amount() == pytest.approx(start, rel=1e-9)
    # the tracer has reached the organ
    assert r[f"[{ids.species(TRACER, ids.zone_id(n))}]"] > 0.1


def test_zero_flow_stays_defined(tmp_path: Path) -> None:
    """`Q_pump = 0` simulates an hour without NaN."""
    path = build(liver_definition(3), tmp_path)
    r = individual(path, **{ids.Q_PUMP: 0.0})
    r.timeCourseSelections = [
        "time",
        *r.model.getFloatingSpeciesIds(),
        *r.model.getGlobalParameterIds(),
        *r.model.getReactionIds(),
    ]
    s = r.simulate(0, 60, 61)
    assert np.isfinite(np.asarray(s)).all()
    assert r[ids.Q_PUMP] == 0.0


def test_every_quantity_at_every_site(tmp_path: Path) -> None:
    """The perfusate quantities, the pO2, the pH and the enzyme activities at
    the reservoir, the inflow, the outflow, every zone and the periportal and
    pericentral means; the cells per zone and their means; the bile."""
    n = 3
    path = build(liver_definition(n, tracers=(TRACER,)), tmp_path)
    known = _ids(_model(path))
    sites = [*_sites(n), ids.PERIPORTAL, ids.PERICENTRAL]
    expected = {
        ids.Q_PUMP,
        ids.TEMPERATURE,
        ids.F_T,
        ids.BILE,
        ids.BILE_RATE,
        ids.VO2,
    }
    for site in sites:
        expected |= {ids.species(q, site) for q in (*ids.PERFUSATE, TRACER)}
        expected |= {ids.po2(site), ids.ph(site)}
        expected |= {ids.activity(e, site) for e in ids.ENZYMES}
    for site in [ids.zone_id(k) for k in range(1, n + 1)] + [
        ids.PERIPORTAL,
        ids.PERICENTRAL,
    ]:
        expected |= {ids.species(q, ids.cell_site(site)) for q in ids.CELL}
        expected.add(ids.alive(site))
    for site in _sites(n):
        expected.add(ids.volume(site))
    for k in range(1, n + 1):
        expected.add(ids.cells(ids.zone_id(k)))
    assert not sorted(expected - known)


def test_the_observers_read_the_sites(tmp_path: Path) -> None:
    """The pO2 is the free oxygen over the solubility, the activity the
    enzyme times its specific activity, the means weigh the zones by their
    overlap with the first and last third of the sinusoid."""
    n = 4
    r = load_circuit(build(liver_definition(n), tmp_path))
    target = entry("perfusion", "perfusate", "pO2").value
    assert r[ids.po2(ids.RESERVOIR)] == pytest.approx(target, rel=1e-12)
    r.simulate(0, 120, 13)
    res = ids.RESERVOIR
    assert r[ids.po2(res)] == pytest.approx(
        r[f"[{ids.species('o2', res)}]"] / r[ids.SOLUBILITY], rel=1e-12
    )
    alt = r[f"[{ids.species('alt', res)}]"]
    assert alt > 0
    assert r[ids.activity("alt", res)] == pytest.approx(
        alt * r[ids.specific_activity("alt")], rel=1e-12
    )

    def zone(q: str, k: int) -> float:
        return r[f"[{ids.species(q, ids.zone_id(k))}]"]

    # n = 4: the first third is zone 1 and a third of zone 2
    pp = (zone("glc", 1) + zone("glc", 2) / 3) / (4 / 3)
    pc = (zone("glc", 4) + zone("glc", 3) / 3) / (4 / 3)
    assert r[ids.species("glc", ids.PERIPORTAL)] == pytest.approx(pp, rel=1e-12)
    assert r[ids.species("glc", ids.PERICENTRAL)] == pytest.approx(pc, rel=1e-12)
    # the oxygen falls along the sinusoid
    assert r[ids.po2(ids.PERIPORTAL)] > r[ids.po2(ids.PERICENTRAL)]
    assert r[ids.BILE] == pytest.approx(
        sum(r[ids.bile(ids.zone_id(k))] for k in range(1, n + 1)) * 1e6, rel=1e-8
    )


def test_the_oxygenator_sets_the_inflow(tmp_path: Path) -> None:
    """Without flow the oxygen of the inflow tubing approaches the target
    first order."""
    path = build(liver_definition(1), tmp_path)
    r = individual(path, **{ids.Q_PUMP: 0.0})
    o2 = f"[{ids.species('o2', ids.INFLOW)}]"
    target = r[o2]
    r[o2] = 0.0
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    r.simulate(0, 1.0, 11)
    k = r[ids.OXYGENATOR_K]
    assert r[o2] == pytest.approx(target * (1 - math.exp(-k)), rel=1e-8)


def test_the_validation_of_a_definition() -> None:
    d = liver_definition(2)
    assert d.validate() == []
    bad = dataclasses.replace(
        d,
        organ=LiverOrgan(
            mass=d.organ.mass,
            zones=2,
            profiles={"OXPHOS": [1.0], "NEC": [1.0, 1.0]},
        ),
        perfusate=dataclasses.replace(
            d.perfusate, carrier="rbc", tracers=("glc", "1x")
        ),
        temperature=Value(37.0, PU.kelvin, "ModelChoice", "wrong unit"),
    )
    problems = "\n".join(bad.validate())
    for text in ("OXPHOS", "NEC", "rbc", "'glc'", "'1x'", "degC"):
        assert text in problems, text
    with pytest.raises(ValueError, match="Invalid perfusion definition"):
        build(bad, Path("."))
    with pytest.raises(ValueError, match="at least one zone"):
        LiverOrgan(mass=d.organ.mass, zones=0)


def test_the_parameter_value() -> None:
    d = liver_definition(1)
    assert d.parameter_value(ids.Q_PUMP) == pytest.approx(3e-3, rel=1e-12)
    assert d.parameter_value(ids.fresh("glc")) == Perfusate.default().glucose.value
    with pytest.raises(KeyError, match="no constant"):
        d.parameter_value("unknown")
