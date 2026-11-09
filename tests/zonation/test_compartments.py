"""The layer compartments of the zones."""

from dataclasses import replace

import pytest
from sbmlutils.factory import Compartment

from pkpdsim.species import Species
from pkpdsim.zonation import ids
from pkpdsim.zonation.builder import ZonationBuilder
from pkpdsim.zonation.geometry import Layer, geometry
from pkpdsim.zonation.profiles import Constant
from pkpdsim.zonation.settings import BileTransport, Geometry, ZonationSettings


def _compartments(**kwargs: object) -> list[Compartment]:
    """The compartments of a builder of the settings."""
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


def test_well_stirred_keeps_the_total_volume_of_the_sinusoid() -> None:
    """Collapsing the layers preserves the volume."""
    for bile in (BileTransport.ADVECTIVE, None):
        totals = [
            sum(float(c.value) for c in _compartments(n_zones=4, bile=bile, geometry=g))  # ty: ignore[invalid-argument-type]
            for g in (Geometry.SINUSOID, Geometry.WELL_STIRRED)
        ]
        assert totals[0] == pytest.approx(totals[1])


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
            ZonationSettings(profiles={"nothing": Constant()}), (), None
        ).build()


def test_every_compartment_follows_the_liver_mass() -> None:
    """Each compartment has an initial assignment, a unit, a name and an SBO term."""
    model = ZonationBuilder(ZonationSettings(n_zones=3), (), None).build()
    assert {a.symbol for a in model.assignments} >= {c.sid for c in model.compartments}
    assert all(c.unit is not None and c.name and c.sboTerm for c in model.compartments)


def test_the_total_volume_is_the_parenchyma() -> None:
    """The layers of all zones fill the parenchyma, whatever the number of zones."""
    from pkpdsim.zonation.geometry import fractions, geometry

    sums = []
    for n in (1, 3, 10):
        model = ZonationBuilder(ZonationSettings(n_zones=n), (), None).build()
        sums.append(sum(float(c.value) for c in model.compartments))  # ty: ignore[invalid-argument-type]
    geo = geometry(Species.RAT)
    assert sums[0] == pytest.approx(sums[1]) == pytest.approx(sums[2])
    assert sums[0] == pytest.approx(
        sum(fractions(Species.RAT).values()) * geo.f_par * 0.007229 / geo.density
    )
