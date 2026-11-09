"""Tests of the geometry, the volume fractions and the haematocrit split."""

import pytest

from pkpdsim.species import Species
from pkpdsim.zonation.geometry import (
    Layer,
    axial_area,
    fractions,
    geometry,
    layer_volumes,
    sinusoid_count,
    table_value,
)

from .conftest import mass_g


def test_rat_fractions_match_blouin() -> None:
    fr = fractions(Species.RAT)
    assert fr[Layer.PLASMA] + fr[Layer.RBC] == pytest.approx(0.106)
    assert fr[Layer.DISSE] == pytest.approx(0.049)
    assert fr[Layer.CELL] == pytest.approx(0.780)
    assert fr[Layer.ENDOTHELIUM] == pytest.approx(0.028)
    assert fr[Layer.BILE] == pytest.approx(0.004)
    # Kupffer 2.1 % and stellate 1.4 % are unassigned
    assert sum(fr.values()) == pytest.approx(0.967)


def test_haematocrit_splits_the_lumen_and_conserves_it() -> None:
    for hct in (0.0, 0.45):
        fr = fractions(Species.RAT, haematocrit=hct)
        assert fr.get(Layer.RBC, 0.0) == pytest.approx(hct * 0.106)
        assert fr[Layer.PLASMA] + fr.get(Layer.RBC, 0.0) == pytest.approx(0.106)


def test_no_erythrocyte_layer_without_erythrocytes() -> None:
    assert Layer.RBC not in fractions(Species.RAT, haematocrit=0.0)


def test_layer_order_is_lumen_to_bile() -> None:
    assert list(Layer) == [
        Layer.RBC,
        Layer.PLASMA,
        Layer.ENDOTHELIUM,
        Layer.DISSE,
        Layer.CELL,
        Layer.BILE,
    ]


def test_human_carries_declared_fallback() -> None:
    assert fractions(Species.HUMAN, haematocrit=0.45) == fractions(
        Species.RAT, haematocrit=0.45
    )


def test_volumes_scale_with_mass_and_split_over_zones() -> None:
    geo, fr = geometry(Species.RAT), fractions(Species.RAT)
    v10 = layer_volumes(mass_g(10.0), 10, geo, fr)
    v1 = layer_volumes(mass_g(10.0), 1, geo, fr)
    assert v10[(Layer.CELL, 1)].value * 10 == pytest.approx(v1[(Layer.CELL, 1)].value)
    v20 = layer_volumes(mass_g(20.0), 1, geo, fr)
    assert v20[(Layer.CELL, 1)].value == pytest.approx(2 * v1[(Layer.CELL, 1)].value)


def test_hepatocyte_volume_is_fraction_of_liver_volume() -> None:
    geo, fr = geometry(Species.RAT), fractions(Species.RAT)
    v = layer_volumes(mass_g(10.0), 1, geo, fr)[(Layer.CELL, 1)]
    # 10 g at 1.07 kg/l is 9.346 ml
    assert v.value == pytest.approx(0.780 * 10.0 / 1070.0)


def test_axial_area_is_volume_over_length() -> None:
    geo, fr = geometry(Species.RAT), fractions(Species.RAT)
    a = axial_area(Layer.PLASMA, mass_g(10.0), geo, fr)
    assert a.value > 0
    v_plasma = 0.55 * 0.106 * 10.0 / 1070.0  # l, the plasma at haematocrit 0.45
    assert a.value == pytest.approx(v_plasma * 1e-3 / geo.L_sin)


def test_sinusoid_count_fills_the_lumen() -> None:
    geo, fr = geometry(Species.RAT), fractions(Species.RAT)
    n = sinusoid_count(mass_g(10.0), geo, fr)
    assert n > 1e5


def test_human_haematocrit_and_density_are_icrp89_not_fallbacks() -> None:
    geo = geometry(Species.HUMAN)
    assert geo.haematocrit == pytest.approx(0.43)
    assert geo.density == pytest.approx(1.0)
    for key in ("haematocrit", "density"):
        e = table_value("liver_human", key)
        assert e.source == "ICRP89"
        assert "fallback" not in (e.notes or "").replace("not a rat fallback", "")
