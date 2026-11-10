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
from sbmlutils.factory import Creator, Model, Parameter, Reaction, Species
from sbmlutils.factory.core_elements import KineticLaw
from sbmlutils.metadata import SBO
from sbmlutils.validation import ValidationOptions, validate_doc

from pkpdsim.database.tables import entry
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion import circuit, core, ids
from pkpdsim.perfusion.builder import PerfusionExtension, build, zone_cell_model
from pkpdsim.perfusion.definition import LiverOrgan, Perfusate
from pkpdsim.perfusion.units import PU
from pkpdsim.templates import terms_of_use

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
    assert path == tmp_path / "perfusion_liver_flat.xml"
    doc = libsbml.readSBMLFromFile(str(path))
    assert (doc.getLevel(), doc.getVersion()) == (3, 2)
    assert doc.getModel().getPlugin("comp") is None
    validation = validate_doc(doc, options=ValidationOptions())
    assert validation.error_count == 0
    assert validation.warning_count == 0
    r = load_circuit(path)
    r.simulate(0, 10, 11)


def test_a_liver_is_written_with_its_cell_model_and_the_flat_model(
    tmp_path: Path,
) -> None:
    """The composed model, the cell model of its zones and the flat model which is simulated."""
    build(liver_definition(3), tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "perfusion_liver.xml",
        "perfusion_liver_core.xml",
        "perfusion_liver_flat.xml",
    ]
    composed = libsbml.readSBMLFromFile(str(tmp_path / "perfusion_liver.xml"))
    plugin = composed.getModel().getPlugin("comp")
    assert [plugin.getSubmodel(k).getId() for k in range(plugin.getNumSubmodels())] == [
        "M01",
        "M02",
        "M03",
    ]


def _problems(path: Path) -> list[tuple[str, int]]:
    """The libsbml errors and warnings of a file."""
    document = libsbml.readSBMLFromFile(str(path))
    document.checkConsistency()
    return [
        (e.getSeverityAsString(), e.getErrorId())
        for e in map(document.getError, range(document.getNumErrors()))
    ]


def test_the_composed_model_of_a_liver_with_cells_is_valid(tmp_path: Path) -> None:
    """The composed model, which carries the structure and which a reader flattens
    for himself, has no libsbml error and no warning: the observers of the cells
    name elements of the composed model, never the ids a flattening produces."""
    build(liver_definition(3), tmp_path)
    assert _problems(tmp_path / "perfusion_liver.xml") == []


def test_the_flat_ids_of_the_cells_are_the_ids_of_the_perfusion(
    tmp_path: Path,
) -> None:
    """The ids of `pkpdsim.perfusion.ids` name elements of the flat model, whether
    the cells of a zone are a submodel or not."""
    definition = liver_definition(3)
    flat = _ids(_model(build(definition, tmp_path)))
    cell = _model(tmp_path / "perfusion_liver_core.xml")
    reactions = [r.getId() for r in cell.getListOfReactions()]
    for k in (1, 2, 3):
        z = ids.zone_id(k)
        expected = {
            *(
                ids.species(q, ids.cell_site(z))
                for q in ids.cell_quantities(definition.released)
            ),
            *(ids.species(q, ids.dead_site(z)) for q in definition.released),
            *(ids.reaction(r, z) for r in reactions if r != "BUF_ext"),
            ids.buffer(z),
            ids.alive(z),
            ids.necrosis_rate(z),
            ids.bile(z),
            ids.ph(z),
        }
        assert expected <= flat, sorted(expected - flat)
        # the proxies are merged with the elements of the submodel, not added: no
        # exposed state is in the flat model under its submodel prefix as well
        cell_ids = {
            *(f"{q}_cell" for q in ids.cell_quantities(definition.released)),
            *(f"{q}_dead" for q in definition.released),
            *reactions,
            "f_alive",
            "v_nec",
            "Vbile",
            "pH_ext",
        }
        left = {f"M{k:02}__{x}" for x in cell_ids} & flat
        assert not left, sorted(left)


def _mechanism(model: Model) -> Model:
    """Add a mechanism of the cells to a cell model: a species, a reaction, a constant."""
    model.parameters.append(
        Parameter(
            "XP_k",
            1.0,
            U.per_min,
            name="rate constant of the mechanism",
            sboTerm=SBO.KINETIC_CONSTANT,
            port=True,
        )
    )
    model.species.append(
        Species(
            "x_cell",
            "Vli",
            initialConcentration=0.0,
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=False,
            name="product of the mechanism (cells)",
            sboTerm=SBO.SIMPLE_CHEMICAL,
        )
    )
    model.reactions.append(
        Reaction(
            "XP",
            "-> x_cell [atp_cell]",
            formula=("XP_k * Vli * atp_cell", U.mmole_per_min),
            compartment="Vli",
            name="production of the product of the mechanism",
            sboTerm=SBO.BIOCHEMICAL_REACTION,
        )
    )
    return model


def test_an_extension_of_the_cells_is_one_mechanism_in_every_zone(
    tmp_path: Path,
) -> None:
    """A mechanism of a package is written into the cell model once and reaches
    every zone: the engine exposes its species and its reaction under the ids of
    the perfusion, the liver holds its constant (one in the flat model, none in a
    submodel) and what the mechanism adds to the composed model is there too."""

    def composed(model: Model) -> Model:
        model.parameters += [
            Parameter(
                "x_reservoir",
                0.0,
                U.mmole,
                name="product of the mechanism in the reservoir",
                sboTerm=SBO.SIMPLE_CHEMICAL,
            ),
            # the shared constant is the liver's, which replaces the copy of
            # every submodel: it is in both models
            Parameter(
                "XP_k",
                1.0,
                U.per_min,
                name="rate constant of the mechanism",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ]
        return model

    zones = (1, 2, 3)
    extension = PerfusionExtension(cell=_mechanism, shared=("XP_k",), composed=composed)
    path = build(liver_definition(len(zones)), tmp_path, extension)
    flat = _ids(_model(path))
    expected = {
        *(ids.species("x", ids.cell_site(ids.zone_id(k))) for k in zones),
        *(ids.reaction("XP", ids.zone_id(k)) for k in zones),
        "XP_k",
        "x_reservoir",
    }
    assert expected <= flat, sorted(expected - flat)
    # the constant of the mechanism is the liver's, not one per zone
    left = {f"M{k:02}__XP_k" for k in zones} & flat
    assert not left, sorted(left)


def test_an_extension_which_replaces_a_law_of_the_library_drops_its_constants(
    tmp_path: Path,
) -> None:
    """A constant of a law which the extension replaces is gone from both models and
    is shared no more: the flat model holds it neither once nor per zone."""
    gone = ("O2IM_Vmax", "O2IM_Km_o2")

    def cell(model: Model) -> Model:
        model = _mechanism(model)
        model.parameters = [p for p in model.parameters if p.sid not in gone]
        (uptake,) = [r for r in model.reactions if r.sid == "O2IM"]
        uptake.formula = KineticLaw("XP_k * Vli * (o2_ext - o2_cell)", U.mmole_per_min)
        return model

    def composed(model: Model) -> Model:
        model.parameters = [
            *(p for p in model.parameters if p.sid not in gone),
            Parameter(
                "XP_k",
                1.0,
                U.per_min,
                name="rate constant of the mechanism",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ]
        return model

    zones = (1, 2)
    extension = PerfusionExtension(cell=cell, shared=("XP_k",), composed=composed)
    flat = _ids(_model(build(liver_definition(len(zones)), tmp_path, extension)))
    left = {
        *gone,
        *(f"M{k:02}__{sid}" for k in zones for sid in gone),
    } & flat
    assert not left, sorted(left)
    assert {ids.reaction("O2IM", ids.zone_id(k)) for k in zones} <= flat


def test_an_extension_releases_into_a_quantity_of_the_perfusate(
    tmp_path: Path,
) -> None:
    """A quantity of the perfusate which the cells of the library do not exchange is
    wired like the others when the extension gives the cell model its port: the
    reaction of the zone releases into the vascular space of the zone, and no copy
    of the port is left in the submodel."""

    def cell(model: Model) -> Model:
        model.species.append(
            Species(
                "tracer_ext",
                "Vext",
                initialConcentration=0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,
                name="tracer of the perfusate of the zone",
                sboTerm=SBO.SIMPLE_CHEMICAL,
                port=True,
            )
        )
        model.reactions.append(
            Reaction(
                "TRACERR",
                "-> tracer_ext",
                formula=("f_alive * Vli * 1 mM * 1 per_min", U.mmole_per_min),
                compartment="Vext",
                name="release of the tracer",
                sboTerm=SBO.TRANSPORT_REACTION,
            )
        )
        return model

    definition = liver_definition(2, tracers=("tracer",))
    flat = _model(build(definition, tmp_path, PerfusionExtension(cell=cell)))
    sids = _ids(flat)
    for k in (1, 2):
        zone = ids.zone_id(k)
        release = flat.getReaction(ids.reaction("TRACERR", zone))
        assert release is not None, ids.reaction("TRACERR", zone)
        products = [r.getSpecies() for r in release.getListOfProducts()]
        assert products == [ids.species("tracer", zone)], products
    left = {f"M{k:02}__tracer_ext" for k in (1, 2)} & sids
    assert not left, sorted(left)


def test_a_shared_element_an_extension_leaves_out_of_the_liver_is_refused(
    tmp_path: Path,
) -> None:
    """A shared constant which the composed model does not hold is an error of the
    extension, named before anything is written: the liver replaces the copy of
    every submodel by its own element."""
    extension = PerfusionExtension(cell=_mechanism, shared=("XP_k",))
    with pytest.raises(ValueError, match="shared elements of its extension"):
        build(liver_definition(2), tmp_path, extension)


def test_the_start_values_of_the_replaced_vascular_co2_and_the_cell_ph(
    tmp_path: Path,
) -> None:
    """The composed model carries the start of `co2_vas_<z>` (the equilibrium of the
    buffer, `h*hco3/Ka`) and of `pH_<z>` (`3 - log10(h/1 mM)`) of the fresh
    perfusate: the expectation is computed from the composition of the definition
    and its buffer, not copied from the builder. No simulation reads either: the
    flat model starts the first by the initial assignment of the cell's `co2_ext`
    and computes the second by the assignment rule of the cell's `pH_ext`."""
    definition = liver_definition(3)
    build(definition, tmp_path)
    model = _model(tmp_path / "perfusion_liver.xml")
    start = circuit.composition(definition.perfusate).concentrations(
        definition.released
    )
    h, hco3 = start["h"], start["hco3"]
    assert h is not None and hco3 is not None
    co2 = h * hco3 / core.Buffer.default().ka.value
    ph = 3 - math.log10(h)
    assert co2 > 0 and 6.5 < ph < 8
    for k in (1, 2, 3):
        z = ids.zone_id(k)
        assert model.getSpecies(ids.species("co2", z)).getInitialConcentration() == (
            pytest.approx(co2, rel=1e-12)
        )
        assert model.getParameter(ids.ph(z)).getValue() == pytest.approx(ph, rel=1e-12)


def test_an_element_the_liver_does_not_expose_is_addressed_through_its_submodel(
    tmp_path: Path,
) -> None:
    """The spec's addressing, `M<kk>__<id>`, stays for the cell elements the liver
    does not observe: the temperature factor of the cell is one per zone, the
    liver has its own."""
    flat = _ids(_model(build(liver_definition(3), tmp_path)))
    assert {"M01__f_T", "M02__f_T", "M03__f_T"} <= flat
    assert "f_T" in flat
    assert not {"f_T_z1", "f_T_z2", "f_T_z3"} & flat


def test_two_definitions_in_one_folder_have_a_cell_model_each(tmp_path: Path) -> None:
    build(liver_definition(2, sid="one"), tmp_path)
    build(liver_definition(2, sid="two", profiles={"OXPHOS": [2, 1]}), tmp_path)
    assert (tmp_path / "one_core.xml").is_file()
    assert (tmp_path / "two_core.xml").is_file()
    assert "p_OXPHOS" not in (tmp_path / "one_core.xml").read_text()
    assert "p_OXPHOS" in (tmp_path / "two_core.xml").read_text()


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
        ids.BILE_PER_MASS,
        ids.VO2,
    }
    for site in sites:
        expected |= {ids.species(q, site) for q in (*ids.PERFUSATE, TRACER)}
        expected |= {ids.po2(site), ids.pco2(site), ids.ph(site)}
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
        sum(r[ids.bile(ids.zone_id(k))] for k in range(1, n + 1)) * 1e6,
        rel=1e-8,
    )
    # µl per g: the bile (µl) over the mass of the organ (kg) in g
    assert r[ids.BILE] > 0
    assert r[ids.BILE_PER_MASS] == pytest.approx(
        r[ids.BILE] / (1000 * r[ids.mass("li")]), rel=1e-12
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


def test_the_validation_of_a_definition(tmp_path: Path) -> None:
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
        build(bad, tmp_path)
    with pytest.raises(ValueError, match="at least one zone"):
        LiverOrgan(mass=d.organ.mass, zones=0)


@pytest.mark.parametrize(
    ("ph", "refused"),
    [(5.99, True), (6.0, False), (7.4, False), (8.0, False), (74.0, True)],
)
def test_the_ph_of_the_perfusate_is_in_its_range(ph: float, refused: bool) -> None:
    """A pH outside of 6 to 8 is refused (a pH of 74 is a typo)."""
    d = liver_definition(1)
    value = dataclasses.replace(d.perfusate.ph, value=ph)
    bad = dataclasses.replace(d, perfusate=dataclasses.replace(d.perfusate, ph=value))
    problems = [p for p in bad.validate() if "pH" in p]
    assert bool(problems) == refused, problems
    if refused:
        assert "outside of 6.0 to 8.0" in problems[0]


def test_a_tracer_has_no_id_of_the_model() -> None:
    """A tracer whose ids are observers, compartments or constants is refused."""
    d = liver_definition(2)
    for tracer, clash in (("pO2", "pO2_res"), ("V", "V_res"), ("V", "V_vas_z1")):
        bad = dataclasses.replace(
            d, perfusate=dataclasses.replace(d.perfusate, tracers=(tracer,))
        )
        problems = "\n".join(bad.validate())
        assert clash in problems, tracer
    assert liver_definition(2, tracers=(TRACER,)).validate() == []


def test_the_fractions_of_the_table_are_checked() -> None:
    """A cell fraction which exceeds 1 with the vascular fraction of the
    table is refused."""
    d = liver_definition(1)
    organ = LiverOrgan(
        mass=d.organ.mass,
        cell_fraction=Value(0.9, PU.dimensionless, "ModelChoice", "too large"),
    )
    problems = dataclasses.replace(d, organ=organ).validate()
    assert any("exceed 1" in p for p in problems)


def test_the_fresh_perfusate_stays_at_its_ph(tmp_path: Path) -> None:
    """Without flow the inflow tubing holds the fresh perfusate: the
    oxygenator neither adds nor vents CO2 at the pCO2 of the fresh perfusate."""
    path = build(liver_definition(1), tmp_path)
    r = individual(path, **{ids.Q_PUMP: 0.0})
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    ph = r[ids.ph(ids.INFLOW)]
    assert ph == pytest.approx(Perfusate.default().ph.value, abs=1e-12)
    assert r[ids.OXYGENATOR_CO2] == pytest.approx(0.0, abs=1e-15)
    r.simulate(0, 360, 37)
    assert r[ids.ph(ids.INFLOW)] == pytest.approx(ph, abs=1e-9)


def test_the_gas_follows_a_change_of_the_fresh_perfusate(tmp_path: Path) -> None:
    """Without a pCO2 of the gas in the definition, the gas is the pCO2 of the
    fresh perfusate (an initial assignment), so a fresh bicarbonate changed
    before the initialization keeps the fresh perfusate at its pH."""
    path = build(liver_definition(1), tmp_path)
    r = individual(path, **{ids.Q_PUMP: 0.0, ids.fresh("hco3"): 12.0})
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    gas = r[ids.CO2_SOLUBILITY] * r[ids.CO2_GAS]
    assert gas == pytest.approx(r[ids.fresh("co2")], rel=1e-12)
    assert r[ids.OXYGENATOR_CO2] == pytest.approx(0.0, abs=1e-15)
    ph = r[ids.ph(ids.INFLOW)]
    assert ph == pytest.approx(Perfusate.default().ph.value, abs=1e-12)
    r.simulate(0, 360, 37)
    assert r[ids.ph(ids.INFLOW)] == pytest.approx(ph, abs=1e-9)


def test_a_given_gas_is_a_constant(tmp_path: Path) -> None:
    """A pCO2 of the gas in the definition is a constant without an assignment."""
    d = liver_definition(1)
    gas = Value(5.0, PU.kPa, "ModelChoice", "a test")
    d = dataclasses.replace(d, circuit=dataclasses.replace(d.circuit, pco2_gas=gas))
    path = build(d, tmp_path)
    r = individual(path, **{ids.fresh("hco3"): 12.0})
    assert r[ids.CO2_GAS] == pytest.approx(5.0, rel=1e-12)


def test_the_oxygenator_vents_co2(tmp_path: Path) -> None:
    """CO2 added at the inflow decays first order towards the gas (the
    hydration of the buffer switched off, so the decay is analytic)."""
    path = build(liver_definition(1), tmp_path)
    r = individual(path, **{ids.Q_PUMP: 0.0, ids.BUFFER_K_HYDRATION: 0.0})
    r.integrator.absolute_tolerance = 1e-14
    r.integrator.relative_tolerance = 1e-10
    co2 = f"[{ids.species('co2', ids.INFLOW)}]"
    gas = r[ids.CO2_SOLUBILITY] * r[ids.CO2_GAS]
    assert r[co2] == pytest.approx(gas, rel=1e-12)
    r[co2] = gas + 2.0
    r.simulate(0, 1.0, 11)
    k = r[ids.OXYGENATOR_K_CO2]
    assert r[co2] == pytest.approx(gas + 2.0 * math.exp(-k), rel=1e-8)


def test_the_parameter_value() -> None:
    d = liver_definition(1)
    assert d.parameter_value(ids.Q_PUMP) == pytest.approx(3e-3, rel=1e-12)
    assert d.parameter_value(ids.fresh("glc")) == Perfusate.default().glucose.value
    with pytest.raises(KeyError, match="no constant"):
        d.parameter_value("unknown")


def test_the_notes_of_a_definition_follow_the_description(tmp_path: Path) -> None:
    """The notes of a definition (the terms of use of a model) follow its description."""
    plain = _model(build(liver_definition(1), tmp_path / "plain")).getNotesString()
    assert "Machine perfusion of the" in plain
    terms = "This work is licensed under the terms of use of the test."
    # an indented block, as the terms of use of a model package
    d = dataclasses.replace(
        liver_definition(1),
        notes=f"""
    The liver of the test.

    ## Terms of use
    {terms}
    """,
    )
    notes = _model(build(d, tmp_path / "notes")).getNotesString()
    assert "Machine perfusion of the" in notes
    assert "<h2>Terms of use</h2>" in notes
    assert f"<p>{terms}</p>" in notes
    assert notes.index("Machine perfusion of the") < notes.index(terms)


def test_the_cell_model_carries_the_notes_and_the_creators_of_its_definition(
    tmp_path: Path,
) -> None:
    """The cell model is a model of its own, so it carries the terms of use and the
    creators of its package as the composed model does: a reader who has only the
    cell model has the licence it is under."""
    terms = "This work is licensed under the terms of use of the test."
    creator = Creator(
        familyName="Tester", givenName="A", email="a@test.org", organization="Test"
    )
    d = dataclasses.replace(
        liver_definition(3),
        creators=(creator,),
        notes=f"""
    The liver of the test.

    ## Terms of use
    {terms}
    """,
    )
    flat = build(d, tmp_path)
    cell = _model(flat.parent / f"{d.sid}_core.xml")
    notes = cell.getNotesString()
    assert "cellular core" in notes
    assert f"<p>{terms}</p>" in notes
    # the creators of a model are its model history
    history = cell.getModelHistory()
    assert history is not None
    assert history.getNumCreators() == 1
    assert history.getCreator(0).getFamilyName() == "Tester"


def test_a_zone_has_no_empty_bile_compartment(tmp_path: Path) -> None:
    """The compartments of the liver are the circuit and, per zone, its vascular
    space and its cells: no bile compartment is left by an unreplaced port."""
    model = _model(build(liver_definition(3), tmp_path))
    compartments = {
        model.getCompartment(k).getId() for k in range(model.getNumCompartments())
    }
    zones = {f"V_{kind}_z{k}" for kind in ("vas", "cell") for k in (1, 2, 3)}
    assert compartments == {"V_res", "V_in", "V_out"} | zones


def test_the_zones_share_the_arrhenius_constants_of_the_circuit(tmp_path: Path) -> None:
    """One `Ea`, `T_ref` and `R_gas` in the liver: the zones have none of their own."""
    model = _model(build(liver_definition(3), tmp_path))
    parameters = {
        model.getParameter(k).getId() for k in range(model.getNumParameters())
    }
    for sid in ("Ea", "T_ref", "R_gas"):
        assert sid in parameters
        assert not {p for p in parameters if p.endswith(f"__{sid}")}


@pytest.mark.parametrize("flow", [0.5e-3, 3.0e-3, 0.0])
def test_the_flow_through_the_liver_is_the_flow_of_the_pump(
    tmp_path: Path, flow: float
) -> None:
    """The study sets its flow rates through `Q_pump`: the engine's `Q_li` follows
    it exactly, whatever the pump is set to (also 0: no flow, nothing divides by it)."""
    path = build(liver_definition(3), tmp_path)
    r = individual(path, **{ids.Q_PUMP: flow})
    assert r["Q_li"] == pytest.approx(flow, abs=0, rel=0)
    assert r[ids.Q_PUMP] == flow


def test_the_cell_model_carries_the_terms_of_use_without_the_circuit(
    tmp_path: Path,
) -> None:
    """The notes of a package are the description of its circuit and its terms of
    use: the cell model of the zones is one cell, so it takes the terms alone."""
    description = "A reservoir of 500 ml, the pump and the oxygenator."
    definition = dataclasses.replace(
        liver_definition(2, sid="noted"),
        notes=f"    {description}\n{terms_of_use}",
    )
    cell = zone_cell_model(definition)
    assert cell is not None and cell.notes is not None
    assert "The cellular core of the liver" in cell.notes
    # the notes are markdown, which sbmlutils converts: `## Terms of use`
    assert "<h2>Terms of use</h2>" in cell.notes
    assert "creativecommons.org" in cell.notes
    assert description not in cell.notes
    circuit_model = _model(build(definition, tmp_path))
    notes = circuit_model.getNotesString()
    assert description in notes and "Terms of use" in notes


def test_the_notes_of_a_package_without_the_terms_of_use_stay_whole(
    tmp_path: Path,
) -> None:
    """Notes the terms of use cannot be told from are not cut: a cell model which
    says too much loses nothing, one without its licence is unusable."""
    definition = dataclasses.replace(
        liver_definition(1, sid="unlicensed"), notes="    Something else entirely."
    )
    cell = zone_cell_model(definition)
    assert cell is not None and cell.notes is not None
    assert "Something else entirely." in cell.notes


def test_the_liver_keeps_the_provenance_of_its_density_and_fractions(
    tmp_path: Path,
) -> None:
    """`rho_li`, `f_vas` and `f_cell` of a perfused liver name their source.

    The organ of the definition holds these three numbers as `Value`s of the
    physiology table with their citation (the vascular fraction of Brown 1997
    Table 30, the cell fraction derived from it), the engine carries them
    through (`SinusoidGeometry.density`, `ZonationSettings.fractions`) and
    `M_li` shows what the notes of the three looked like before the liver
    moved onto the engine: a reader of the flat model or of the ODE markdown
    has to be able to see that the vascular fraction is a measurement, which
    is the stated purpose of `value_parameter`.
    """
    model = _model(build(liver_definition(), tmp_path / "provenance"))
    notes = {
        pid: model.getParameter(pid).getNotesString()
        for pid in ("M_li", "rho_li", "f_vas", "f_cell")
    }
    for pid, text in notes.items():
        assert "Source:" in text, pid
        assert "doi:" in text or "Choice of a library of pkpdsim" in text, pid
    assert "Brown" in notes["f_vas"], "the vascular fraction is Brown 1997 Table 30"
    assert "Table 30" in notes["f_vas"]
    assert "derived" in notes["f_cell"]
