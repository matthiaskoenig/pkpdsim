"""The layer compartments of the zones."""

from dataclasses import replace
from pathlib import Path

import pytest
import roadrunner
from sbmlutils.factory import Compartment, create_model

from pkpdsim.database.tables import entry
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion.units import PU
from pkpdsim.species import Species
from pkpdsim.zonation import ids
from pkpdsim.zonation.builder import ZonationBuilder
from pkpdsim.zonation.geometry import Layer, fractions, geometry
from pkpdsim.zonation.profiles import Constant
from pkpdsim.zonation.settings import BileTransport, Geometry, ZonationSettings

from .conftest import liver_flow


def _total(**kwargs: object) -> float:
    """The total volume (l) of the compartments of the settings."""
    return sum(float(c.value) for c in _compartments(**kwargs))  # ty: ignore[invalid-argument-type]


def _table_mass() -> float:
    """The reference liver mass (kg) of the physiology table of the rat."""
    return entry("physiology_rat", "male.organs.li", "mass").value  # ty: ignore[invalid-return-type]


def _compartments(**kwargs: object) -> list[Compartment]:
    """The compartments of a builder of the settings."""
    kwargs.setdefault("reference_flow", liver_flow())
    return ZonationBuilder(ZonationSettings(**kwargs), (), None).compartments()  # ty: ignore[invalid-argument-type]


def test_sinusoid_has_six_layer_compartments_per_zone() -> None:
    """Six layers in each of ten zones."""
    m = _compartments(n_zones=10)
    assert len(m) == 60
    assert ids.volume(Layer.DISSE, ids.zone(10)) in {c.sid for c in m}


def test_well_stirred_has_three_compartments_per_zone() -> None:
    """A vascular, a cell and a bile compartment per zone."""
    # the spec collapses the erythrocyte, endothelial and Disse compartments only
    m = _compartments(n_zones=10, geometry=Geometry.WELL_STIRRED)
    assert len(m) == 30


def test_without_bile_the_compartments_are_the_old_counts() -> None:
    """`bile=None` omits the bile compartment in both geometries."""
    assert len(_compartments(n_zones=10, bile=None)) == 50
    assert (
        len(_compartments(n_zones=10, bile=None, geometry=Geometry.WELL_STIRRED)) == 20
    )


@pytest.mark.parametrize("bile", [BileTransport.ADVECTIVE, None])
@pytest.mark.parametrize("haematocrit", [None, 0.0])
def test_well_stirred_keeps_the_total_volume_of_the_sinusoid(
    bile: BileTransport | None, haematocrit: float | None
) -> None:
    """Collapsing the layers preserves the volume."""
    totals = [
        _total(n_zones=4, bile=bile, haematocrit=haematocrit, geometry=g)
        for g in (Geometry.SINUSOID, Geometry.WELL_STIRRED)
    ]
    assert totals[0] == pytest.approx(totals[1])


def test_well_stirred_without_erythrocytes_has_three_compartments_per_zone() -> None:
    """The vascular compartment exists at haematocrit 0."""
    m = _compartments(n_zones=10, geometry=Geometry.WELL_STIRRED, haematocrit=0.0)
    assert len(m) == 30


def test_the_geometry_constants_replace_the_table() -> None:
    """A given geometry changes the parenchymal fraction and so the volumes."""
    geo = replace(geometry(Species.RAT), f_par=0.5)
    given = _compartments(n_zones=2, geometry_constants=geo)
    table = _compartments(n_zones=2)
    assert float(given[0].value) != float(table[0].value)  # ty: ignore[invalid-argument-type]


def test_without_haematocrit_there_is_no_erythrocyte_compartment() -> None:
    """Five layers per zone at haematocrit 0."""
    m = _compartments(n_zones=4, haematocrit=0.0)
    assert len(m) == 20
    assert not any("rbc" in c.sid for c in m)


def test_one_zone_is_legal() -> None:
    """A single zone is the unzonated organ."""
    assert len(_compartments(n_zones=1)) == 6


def test_the_default_is_ten_zones() -> None:
    """The default number of zones."""
    assert ZonationSettings().n_zones == 10


def test_less_than_one_zone_raises() -> None:
    """A liver has at least one zone."""
    with pytest.raises(ValueError, match="at least one zone"):
        ZonationSettings(n_zones=0)


def test_an_unused_profile_key_raises() -> None:
    """A profile of a process nothing uses is refused before anything is built."""
    with pytest.raises(ValueError, match="no process"):
        ZonationBuilder(
            ZonationSettings(
                profiles={"axial_diffusion": Constant(), "nothing": Constant()}
            ),
            (),
            None,
        ).build()
    ZonationBuilder(
        ZonationSettings(profiles={"axial_diffusion": Constant()}), (), None
    )


def test_every_compartment_has_unit_name_and_sbo_term() -> None:
    """Each compartment is annotated."""
    settings = ZonationSettings(n_zones=3, reference_flow=liver_flow())
    model = ZonationBuilder(settings, (), None).build()
    assert all(c.unit is not None and c.name and c.sboTerm for c in model.compartments)


@pytest.mark.parametrize("geom", list(Geometry))
def test_the_initial_assignments_evaluate_to_the_layer_volumes(
    geom: Geometry, tmp_path: Path
) -> None:
    """The written model gives every compartment `f_<layer> * f_par * M_li / rho_li / n`."""
    n = 3
    settings = ZonationSettings(n_zones=n, geometry=geom, reference_flow=liver_flow())
    model = ZonationBuilder(settings, (), None).build()
    path = tmp_path / "liver.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    rr = roadrunner.RoadRunner(str(path))
    par = rr[ids.PARENCHYMA] * rr[ids.MASS] / rr[ids.DENSITY]
    # the reference liver of the table, in kg and kg/l
    assert rr[ids.MASS] == pytest.approx(_table_mass())
    # the compartments which are no layer of a zone are exactly the duct
    zone_layers = {
        ids.volume(layer, ids.zone(k))
        for layer in (*Layer, ids.VASCULAR)
        for k in range(1, n + 1)
    }
    assert {c.sid for c in model.compartments} - zone_layers == {
        ids.volume(Layer.BILE, ids.DUCT)
    }
    for c in model.compartments:
        if c.sid == ids.volume(Layer.BILE, ids.DUCT):
            continue  # a collector of amounts, asserted below as the one extra
        layer = c.sid.removeprefix("V_").rsplit("_", 1)[0]
        assert rr[c.sid] == pytest.approx(rr[ids.fraction(layer)] * par / n)
        assert rr[c.sid] > 0.0


def test_the_total_volume_is_the_parenchyma() -> None:
    """The layers of all zones fill the parenchyma, whatever the number of zones."""
    geo = geometry(Species.RAT)
    expected = (
        sum(fractions(Species.RAT).values()) * geo.f_par * _table_mass() / geo.density
    )
    for n in (1, 3, 10):
        assert _total(n_zones=n) == pytest.approx(expected)


def test_a_liver_mass_in_another_mass_unit_is_converted() -> None:
    """A mass in g gives the volumes of the same mass in kg."""
    in_g = _total(liver_mass=Value(7.229, PU.g))
    in_kg = _total(liver_mass=Value(0.007229, U.kg))
    assert in_g == pytest.approx(in_kg)


def test_a_liver_mass_without_a_mass_unit_raises() -> None:
    """The unit is named in the error."""
    with pytest.raises(ValueError, match="needs a mass"):
        _compartments(liver_mass=Value(1.0, U.litre))
