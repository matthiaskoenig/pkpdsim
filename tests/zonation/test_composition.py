"""The cell model per zone through comp."""

from pathlib import Path

import libsbml
import numpy as np
import pytest
import roadrunner
from sbmlutils.factory import PORT_SUFFIX, Model

from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.zonation import ids
from pkpdsim.zonation.builder import ZonationBuilder
from pkpdsim.zonation.composition import CellPorts, compose, write
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.settings import Geometry, ZonationSettings

from .conftest import mass_g
from .data.toy_cell import TOY_ID, UPTAKE, create_toy_cell

DRUG = Substance("drg", "drug", Value(180.0, U.g_per_mole))


def _builder(n: int = 10) -> ZonationBuilder:
    """A builder of `n` zones with the drug."""
    return ZonationBuilder(
        ZonationSettings(n_zones=n, liver_mass=mass_g(1000.0)), [DRUG], None
    )


def _written(tmp_path: Path, n: int = 10) -> Path:
    """The flat model of the toy cell composed into `n` zones."""
    create_toy_cell(tmp_path)
    return write(compose(_builder(n), TOY_ID, CellPorts()), tmp_path)


def _replaced_value(model: Model, submodel: str, port: str) -> float:
    """The value of the parameter which replaces a port of a submodel."""
    (element,) = (
        r
        for r in model.replaced_elements
        if r.submodelRef == submodel and r.portRef == f"{port}{PORT_SUFFIX}"
    )
    (parameter,) = (p for p in model.parameters if p.sid == element.elementRef)
    return float(parameter.value)  # ty: ignore[invalid-argument-type]


def test_one_submodel_per_zone() -> None:
    """M01 to M10, each of the one external definition."""
    m = compose(_builder(10), TOY_ID, CellPorts())
    assert [s.sid for s in m.submodels] == [f"M{k:02}" for k in range(1, 11)]
    assert len(m.external_model_definitions) == 1
    assert {s.modelRef for s in m.submodels} == {m.external_model_definitions[0].sid}
    assert m.external_model_definitions[0].source == f"{TOY_ID}.xml"


def test_the_disse_compartment_replaces_the_external_port() -> None:
    """The external port is the space of Disse, the cell the hepatocyte, bile the canaliculus."""
    m = compose(_builder(10), TOY_ID, CellPorts())
    found = {
        (r.submodelRef, r.portRef): r.elementRef
        for r in m.replaced_elements
        if r.submodelRef == "M03"
    }
    assert found[("M03", f"Vext{PORT_SUFFIX}")] == ids.volume(Layer.DISSE, "z3")
    assert found[("M03", f"Vli{PORT_SUFFIX}")] == ids.volume(Layer.CELL, "z3")
    assert found[("M03", f"Vbi{PORT_SUFFIX}")] == ids.volume(Layer.BILE, "z3")
    assert found[("M03", f"drg_ext{PORT_SUFFIX}")] == ids.species(
        "drg", Layer.DISSE, "z3"
    )


def test_every_compartment_carries_a_meta_id() -> None:
    """The flattening needs it, of the builder alone as well."""
    for m in (_builder(10).build(), compose(_builder(10), TOY_ID, CellPorts())):
        assert m.compartments
        assert all(c.metaId == f"meta_{c.sid}" for c in m.compartments)


def test_position_is_the_zone_midpoint() -> None:
    """(k - 0.5) / 10 of the zone."""
    m = compose(_builder(10), TOY_ID, CellPorts())
    assert _replaced_value(m, "M01", "position") == pytest.approx(0.05)
    assert _replaced_value(m, "M10", "position") == pytest.approx(0.95)


def test_temperature_is_one_global_parameter() -> None:
    """Every submodel takes the one parameter."""
    m = compose(_builder(10), TOY_ID, CellPorts())
    elements = {
        r.elementRef
        for r in m.replaced_elements
        if r.portRef == f"temperature{PORT_SUFFIX}"
    }
    assert elements == {"temperature"}
    assert len([p for p in m.parameters if p.sid == "temperature"]) == 1


def test_deletions_name_the_port_of_every_submodel() -> None:
    """A deletion per zone."""
    m = compose(_builder(3), TOY_ID, CellPorts(), deletions=("k_up",))
    assert {(d.submodelRef, d.portRef) for d in m.deletions} == {
        (f"M{k:02}", f"k_up{PORT_SUFFIX}") for k in (1, 2, 3)
    }


def test_write_requires_the_cell_model_next_to_it(tmp_path: Path) -> None:
    """A clear error names the expected path."""
    with pytest.raises(FileNotFoundError, match=f"{TOY_ID}.xml"):
        write(compose(_builder(2), TOY_ID, CellPorts()), tmp_path)


def test_the_flat_model_has_no_submodels_left(tmp_path: Path) -> None:
    """Flattened, and the unflat model is kept next to it."""
    path = _written(tmp_path)
    document = libsbml.readSBMLFromFile(str(path))
    plugin = document.getModel().getPlugin("comp")
    assert plugin is None or plugin.getNumSubmodels() == 0
    assert (tmp_path / "zonated_liver_rat.xml").is_file()


ACCEPTED_WARNINGS: dict[int, str] = {}
"""The libsbml warnings (by number, with the reason) the flat model may emit: none."""


def test_the_flat_model_validates(tmp_path: Path) -> None:
    """No libsbml error and no warning, the units checks included."""
    document = libsbml.readSBMLFromFile(str(_written(tmp_path)))
    document.checkConsistency()
    problems = [
        (e.getSeverityAsString(), e.getErrorId(), e.getMessage()[:300])
        for e in map(document.getError, range(document.getNumErrors()))
        if e.getErrorId() not in ACCEPTED_WARNINGS
    ]
    assert problems == []


def test_the_flat_model_simulates(tmp_path: Path) -> None:
    """The toy takes up the drug of the Disse space of zone 1 into its cell."""
    r = roadrunner.RoadRunner(str(_written(tmp_path)))
    assert np.all(np.isfinite(r.simulate(0, 60, 61)))
    r.reset()
    r[f"[{ids.species('drg', Layer.DISSE, 'z1')}]"] = 1.0
    # the rate of the toy at t = 0 is k * C * V of its external compartment,
    # which is the Disse space of zone 1: 0.5 /min * 1 mmol/l * V_disse_z1
    assert r["M01__UP"] == pytest.approx(
        UPTAKE * 1.0 * r[ids.volume(Layer.DISSE, "z1")]
    )
    # each zone has its own cell: zone 10 has no drug yet and takes up none
    assert r["M10__UP"] == 0.0
    # and the hepatocytes of zone 1 hold what the cell model took up
    r.simulate(0, 1.0, 11)
    assert r["M01__drg"] > 0.0
    # the temperature and the positions of the zones are in the flat model
    assert r[ids.TEMPERATURE] == 37.0
    assert r[ids.position("z10")] == pytest.approx(0.95)


def test_without_a_bile_no_port_is_replaced_by_it() -> None:
    """The bile port is left alone, the other four kinds of port stay."""
    b = ZonationBuilder(
        ZonationSettings(n_zones=2, bile=None, liver_mass=mass_g(1000.0)),
        [DRUG],
        None,
    )
    m = compose(b, TOY_ID, CellPorts())
    ports = {r.portRef for r in m.replaced_elements if r.submodelRef == "M01"}
    assert ports == {
        f"{p}{PORT_SUFFIX}"
        for p in ("Vext", "Vli", "position", "temperature", "drg_ext")
    }


def test_the_well_stirred_vascular_site_replaces_the_external_port() -> None:
    """The Disse space is merged into the vascular site there."""
    b = ZonationBuilder(
        ZonationSettings(
            n_zones=2, geometry=Geometry.WELL_STIRRED, liver_mass=mass_g(1000.0)
        ),
        [DRUG],
        None,
    )
    m = compose(b, TOY_ID, CellPorts())
    found = {
        r.portRef: r.elementRef for r in m.replaced_elements if r.submodelRef == "M02"
    }
    assert found[f"Vext{PORT_SUFFIX}"] == ids.volume(ids.VASCULAR, "z2")
    assert found[f"drg_ext{PORT_SUFFIX}"] == ids.species("drg", ids.VASCULAR, "z2")


def test_every_port_is_replaced_exactly_once() -> None:
    """The full set of replacements of two zones, and no pair twice."""
    m = compose(_builder(2), TOY_ID, CellPorts())
    triples = [(r.submodelRef, r.portRef, r.elementRef) for r in m.replaced_elements]
    expected = {
        (
            f"M0{k}",
            f"{port}{PORT_SUFFIX}",
            element,
        )
        for k in (1, 2)
        for port, element in (
            ("Vext", f"V_disse_z{k}"),
            ("Vli", f"V_cell_z{k}"),
            ("Vbi", f"V_bile_z{k}"),
            ("drg_ext", f"drg_disse_z{k}"),
            ("position", f"position_z{k}"),
            ("temperature", "temperature"),
        )
    }
    assert set(triples) == expected
    pairs = [(t[0], t[1]) for t in triples]
    assert len(pairs) == len(set(pairs)) == len(triples)


def test_the_flat_cell_lives_in_the_hepatocyte_of_its_zone(tmp_path: Path) -> None:
    """Every zone has its own cell species in its own hepatocyte, fed from its Disse space."""
    n = 4
    doc = libsbml.readSBMLFromFile(str(_written(tmp_path, n)))
    model = doc.getModel()
    for k in range(1, n + 1):
        z = f"z{k}"
        assert model.getSpecies(f"M{k:02}__drg").getCompartment() == ids.volume(
            Layer.CELL, z
        )
        reaction = model.getReaction(f"M{k:02}__UP")
        assert [
            reaction.getReactant(i).getSpecies()
            for i in range(reaction.getNumReactants())
        ] == [ids.species("drg", Layer.DISSE, z)]
        assert [
            reaction.getProduct(i).getSpecies()
            for i in range(reaction.getNumProducts())
        ] == [f"M{k:02}__drg"]


def _problems(path: Path) -> list[tuple[str, int]]:
    """The libsbml errors and warnings of a file."""
    document = libsbml.readSBMLFromFile(str(path))
    document.checkConsistency()
    return [
        (e.getSeverityAsString(), e.getErrorId())
        for e in map(document.getError, range(document.getNumErrors()))
    ]


def _counts(path: Path) -> tuple[int, int]:
    """The numbers of compartments and species of a file."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    return model.getNumCompartments(), model.getNumSpecies()


@pytest.mark.parametrize(
    ("kwargs", "compartments", "species"),
    [
        ({"n_zones": 10}, 61, 71),
        ({"n_zones": 1}, 7, 8),
        ({"n_zones": 10, "geometry": Geometry.WELL_STIRRED}, 31, 41),
        # the bile port is not replaced: ten local bile compartments remain
        ({"n_zones": 10, "bile": None}, 60, 60),
        ({"n_zones": 10, "haematocrit": 0.0}, 51, 61),
    ],
    ids=["ten", "one", "well_stirred", "no_bile", "no_haematocrit"],
)
def test_configurations_validate_flat_and_unflat(
    tmp_path: Path, kwargs: dict[str, object], compartments: int, species: int
) -> None:
    """Both files validate on their own, and the flat model has these sizes."""
    create_toy_cell(tmp_path)
    builder = ZonationBuilder(
        ZonationSettings(liver_mass=mass_g(1000.0), **kwargs),  # ty: ignore[invalid-argument-type]
        [DRUG],
        None,
    )
    model = compose(builder, TOY_ID, CellPorts())
    flat = write(model, tmp_path)
    assert _problems(flat) == []
    assert _problems(tmp_path / f"{model.sid}.xml") == []
    assert _counts(flat) == (compartments, species)


def test_the_deleted_element_is_absent_from_the_flat_model(tmp_path: Path) -> None:
    """A deletion which does nothing would leave the reaction in."""
    create_toy_cell(tmp_path)
    flat = write(compose(_builder(3), TOY_ID, CellPorts(), deletions=("UP",)), tmp_path)
    model = libsbml.readSBMLFromFile(str(flat)).getModel()
    assert [model.getReaction(i).getId() for i in range(model.getNumReactions())].count(
        "M01__UP"
    ) == 0
    assert not [
        i
        for i in range(model.getNumReactions())
        if model.getReaction(i).getId().endswith("__UP")
    ]
    assert _problems(flat) == []


def test_the_submodel_ids_stay_unique_above_99_zones() -> None:
    """M100 widens the padding and the ids stay distinct."""
    m = compose(_builder(100), TOY_ID, CellPorts())
    sids = [s.sid for s in m.submodels]
    assert len(sids) == len(set(sids)) == 100
    assert sids[-1] == "M100"


def test_a_cell_model_with_annotated_species_flattens_and_validates(
    tmp_path: Path,
) -> None:
    """An annotated species of the cell model has a metaid, which comp requires
    of the species of the liver which replaces it: `compose` gives it one."""
    create_toy_cell(tmp_path, annotated=True)
    model = compose(_builder(3), TOY_ID, CellPorts())
    path = write(model, tmp_path)
    document = libsbml.readSBMLFromFile(str(path))
    assert document.getNumErrors() == 0
    assert document.checkConsistency() == 0
    flat = document.getModel()
    assert flat.getSpecies("M01__drg_ext") is None, "replaced by the Disse species"
    assert flat.getNumSpecies() > 0
