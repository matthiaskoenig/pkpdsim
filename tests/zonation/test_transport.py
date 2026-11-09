"""Tests of the derived diffusion, the override and the sieving of the fenestrae.

The anchors are aqueous diffusion coefficients at 37 C in 1e-9 m**2/s. The
sources of the measured values (insulin is not sourced) are in the docstring of
`pkpdsim.zonation.transport`; the value at 37 C is the measured one scaled
with `T / eta` (Stokes-Einstein) where it was measured at 20 or 25 C.
"""

import pytest

from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion.units import PU
from pkpdsim.species import Species
from pkpdsim.zonation.geometry import geometry
from pkpdsim.zonation.transport import (
    CM2_PER_S,
    M2_PER_S,
    aqueous_diffusion,
    cytosolic_diffusion,
    diffusion_of,
    diffusion_parameter,
    fenestral_permeability,
    molecular_radius,
    renkin,
)

ANCHORS: dict[str, tuple[float, float]] = {
    # name: (Mr in g/mole, D in 1e-9 m**2/s at 37 C)
    "glucose": (180.16, 0.90),  # 6.7e-10 at 25 C, biofilmbook Table 4-1
    "sucrose": (342.30, 0.70),  # 5.2e-10 at 25 C, biofilmbook Table 4-1
    "gd_chelate": (
        558.6,
        0.50,
    ),  # Gd-DOTA, 2.92e-10 in plasma at 37 C, viscosity scaled
    "insulin": (5808.0, 0.20),  # monomer, 1.5e-10 at 25 C
    "albumin": (66500.0, 0.093),  # 6.1e-11 at 20 C, Sober 1970
}


def _mr(x: float) -> Value:
    return Value(x, PU.g_per_mole)


def _drug(**kwargs: object) -> Substance:
    return Substance(sid="drg", name="drug", mr=_mr(300.0), **kwargs)  # ty: ignore[invalid-argument-type]


@pytest.mark.parametrize("name", list(ANCHORS))
def test_derived_diffusion_matches_the_anchor_within_a_factor(name: str) -> None:
    mr, expected = ANCHORS[name]
    got = aqueous_diffusion(_mr(mr)).value * 1e9
    assert 1 / 1.5 <= got / expected <= 1.5


def test_insulin_is_within_a_factor_of_two() -> None:
    # Not sourced: insulin self-associates, its D depends on concentration, pH and
    # zinc and monomer and hexamer differ by a factor of about 1.6, so the gate
    # is a factor of 2 around the indicative 0.20e-9 m2/s (1.5e-10 at 25 C).
    got = aqueous_diffusion(_mr(5808.0)).value * 1e9
    assert 0.5 <= got / 0.20 <= 2.0


def test_diffusion_decreases_with_the_molar_mass() -> None:
    values = [aqueous_diffusion(_mr(x)).value for x in (100, 1e3, 1e4, 1e5)]
    assert values == sorted(values, reverse=True)


def test_a_plain_molar_mass_is_g_per_mole() -> None:
    assert aqueous_diffusion(342.3) == aqueous_diffusion(_mr(342.3))


def test_the_radius_of_albumin_is_a_few_nm() -> None:
    assert 3e-9 <= molecular_radius(_mr(66500.0)).value <= 4.2e-9


def test_cytosolic_diffusion_is_a_fraction_of_the_aqueous_one() -> None:
    d = aqueous_diffusion(_mr(180.16))
    assert cytosolic_diffusion(d).value == pytest.approx(0.27 * d.value)
    assert cytosolic_diffusion(d, 0.5).value == pytest.approx(0.5 * d.value)


def test_renkin_is_monotone_and_closes_at_one() -> None:
    values = [renkin(x) for x in (0.0, 0.25, 0.5, 0.75)]
    assert values == sorted(values, reverse=True)
    assert renkin(0.0) == 1.0
    assert renkin(1.0) == pytest.approx(0.0)
    assert renkin(1.5) == 0.0


def test_renkin_is_never_negative() -> None:
    assert all(renkin(i / 1000) >= 0.0 for i in range(1001))


def test_fenestral_permeability_follows_the_formula() -> None:
    geo = geometry(Species.RAT)
    d = aqueous_diffusion(_mr(180.16))
    r = Value(0.365e-9, U.meter)
    p = fenestral_permeability(d, r, geo, 0.07)
    expected = 0.07 * d.value / geo.t_endo * renkin(r.value / geo.r_fen)
    assert p.value == pytest.approx(expected)


def test_an_explicit_diffusion_wins_over_the_correlation() -> None:
    s = _drug(diffusion=Value(1e-10, M2_PER_S))
    assert diffusion_of(s).value == pytest.approx(1e-10)


def test_diffusion_in_other_units_is_converted_not_reread() -> None:
    s = _drug(diffusion=Value(1e-6, CM2_PER_S))
    d = diffusion_of(s)
    assert d.value == pytest.approx(1e-10)
    assert d.unit == M2_PER_S


def test_without_an_override_the_diffusion_is_derived() -> None:
    assert diffusion_of(_drug()) == aqueous_diffusion(_mr(300.0))


@pytest.mark.parametrize("bad", [Value(1.0, PU.ml), 1e-10])
def test_a_diffusion_that_is_not_an_area_per_time_raises(bad: object) -> None:
    with pytest.raises(ValueError, match="area per time"):
        _drug(diffusion=bad)


def test_a_non_positive_diffusion_raises() -> None:
    with pytest.raises(ValueError, match="positive"):
        _drug(diffusion=Value(0.0, M2_PER_S))


def test_a_derived_parameter_carries_its_provenance() -> None:
    p = diffusion_parameter(_drug())
    assert p.sid == "D_drg"
    assert "derived" in (p.notes or "")
    assert "300" in (p.notes or "")
    assert p.value == pytest.approx(aqueous_diffusion(_mr(300.0)).value)


def test_a_given_parameter_says_so() -> None:
    p = diffusion_parameter(_drug(diffusion=Value(1e-6, CM2_PER_S)))
    assert "given" in (p.notes or "")
    assert "derived" not in (p.notes or "").replace("not derived", "")
    assert p.value == pytest.approx(1e-10)
