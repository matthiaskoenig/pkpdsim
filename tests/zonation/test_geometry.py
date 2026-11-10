"""Tests of the geometry, the volume fractions and the haematocrit split."""

import math

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
    fr = {layer: f.value for layer, f in fractions(Species.RAT).items()}
    assert fr[Layer.PLASMA] + fr[Layer.RBC] == pytest.approx(0.106)
    assert fr[Layer.DISSE] == pytest.approx(0.049)
    assert fr[Layer.CELL] == pytest.approx(0.780)
    assert fr[Layer.ENDOTHELIUM] == pytest.approx(0.028)
    assert fr[Layer.BILE] == pytest.approx(0.004)
    # Kupffer 2.1 % and stellate 1.4 % are unassigned
    assert sum(fr.values()) == pytest.approx(0.967)


def test_haematocrit_splits_the_lumen_and_conserves_it() -> None:
    for hct in (0.0, 0.45):
        fr = {layer: f.value for layer, f in fractions(Species.RAT, hct).items()}
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
    """The human fractions are the rat ones, and every value says so.

    The design requires the human entries to fall back to the rat fractions
    "with that stated in the table and in the model notes", so the note of
    every human fraction carries the `fallback: rat` of its row - the derived
    shares of the lumen (the plasma and the erythrocytes) in the notes of the
    row they rest on - and reaches the notes of the parameter the builder
    writes from it (`tests/zonation/test_compartments.py`).
    """
    human = fractions(Species.HUMAN, haematocrit=0.45)
    rat = fractions(Species.RAT, haematocrit=0.45)
    assert {layer: f.value for layer, f in human.items()} == {
        layer: f.value for layer, f in rat.items()
    }
    for layer, f in human.items():
        assert "fallback: rat" in (f.notes or ""), layer


def test_volumes_scale_with_mass_and_split_over_zones() -> None:
    geo, fr = geometry(Species.RAT), fractions(Species.RAT)
    v10 = layer_volumes(mass_g(10.0), 10, geo, fr)
    v1 = layer_volumes(mass_g(10.0), 1, geo, fr)
    assert v10[(Layer.CELL, 1)].value * 10 == pytest.approx(v1[(Layer.CELL, 1)].value)
    v20 = layer_volumes(mass_g(20.0), 1, geo, fr)
    assert v20[(Layer.CELL, 1)].value == pytest.approx(2 * v1[(Layer.CELL, 1)].value)


def test_hepatocyte_volume_is_a_fraction_of_the_parenchyma() -> None:
    geo, fr = geometry(Species.RAT), fractions(Species.RAT)
    v = layer_volumes(mass_g(10.0), 1, geo, fr)[(Layer.CELL, 1)]
    # 10 g at the 1.05 kg/l of the library is 9.524 ml of parenchyma (f_par 1)
    assert v.value == pytest.approx(0.780 * 10.0 / 1050.0)


def test_axial_area_is_volume_over_length() -> None:
    geo, fr = geometry(Species.RAT), fractions(Species.RAT)
    a = axial_area(Layer.PLASMA, mass_g(10.0), geo, fr)
    assert a.value > 0
    v_plasma = 0.54 * 0.106 * 10.0 / 1050.0  # l, the plasma at haematocrit 0.46
    assert a.value == pytest.approx(v_plasma * 1e-3 / geo.L_sin.value)


def test_sinusoid_count_fills_the_lumen() -> None:
    """The count is the lumen over the volume of one sinusoid, by hand.

    10 g of rat liver at 1.05 kg/l is 9.523810 ml of parenchyma, the lumen
    10.6 % of it: 1.009524e-6 m3. One sinusoid is a cylinder of the radius
    4.4 um and the length 300 um, pi * (4.4e-6)**2 * 300e-6 = 1.8246370e-14
    m3, so 55 327 378 of them hold the lumen. The arithmetic is written out
    here and not taken from the geometry, so a radius which is not squared
    (12568 sinusoids) or a missing pi (1.7378e7) fails.
    """
    geo, fr = geometry(Species.RAT), fractions(Species.RAT)
    n = sinusoid_count(mass_g(10.0), geo, fr)
    one_sinusoid = math.pi * (4.4e-6) ** 2 * 300e-6
    assert one_sinusoid == pytest.approx(1.8246370e-14, rel=1e-6)
    assert n == pytest.approx(1.009524e-6 / 1.8246370e-14, rel=1e-5)
    assert n == pytest.approx(5.5327378e7, rel=1e-6)


def test_human_haematocrit_and_density_are_human_values_not_fallbacks() -> None:
    """The two human rows which are not the rat's, each with the source it has.

    The haematocrit is the male reference of ICRP 89; the density is the
    convention of `physiology_human.toml` that an organ to which ICRP 89
    gives no specific gravity has 1.0 kg/l, which its own note states, so it
    is a `ModelChoice` as the same number is in `perfusion.toml` and not a
    value cited to a publication which does not give it.
    """
    geo = geometry(Species.HUMAN)
    assert geo.haematocrit.value == pytest.approx(0.43)
    assert geo.density.value == pytest.approx(1.0)
    assert table_value("liver_human", "haematocrit").source == "ICRP89"
    assert table_value("liver_human", "density").source == "ModelChoice"
    for key in ("haematocrit", "density"):
        notes = table_value("liver_human", key).notes or ""
        assert "fallback: rat" not in notes
