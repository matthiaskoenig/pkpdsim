"""Tests of the lobe sets of the liver and their tables."""

import pytest

from pkpdsim.pbpk import U
from pkpdsim.pbpk.lobes import (
    HABR,
    LOBE_SPECIES,
    Lobe,
    check_lobes,
    lobe_table,
)
from pkpdsim.sbml.values import Value
from pkpdsim.species import Species

RAT_TOKENS = ("cl", "lll", "lml", "rl", "rml")
MOUSE_TOKENS = ("cl", "ll", "ml", "ri", "rs")
HUMAN_TOKENS = tuple(f"sg{i}" for i in range(1, 9))


def _share(value: float) -> Value:
    return Value(value, U.dimensionless, "ModelChoice")


def test_lobe_table_rat_tokens() -> None:
    assert tuple(lobe.token for lobe in lobe_table(Species.RAT)) == RAT_TOKENS


def test_lobe_table_mouse_tokens() -> None:
    assert tuple(lobe.token for lobe in lobe_table(Species.MOUSE)) == MOUSE_TOKENS


def test_lobe_table_human_is_couinaud() -> None:
    assert tuple(lobe.token for lobe in lobe_table(Species.HUMAN)) == HUMAN_TOKENS


@pytest.mark.parametrize("species", sorted(LOBE_SPECIES))
def test_lobe_table_mass_fractions_sum_to_one(species: Species) -> None:
    total = sum(lobe.mass_fraction.value for lobe in lobe_table(species))
    assert total == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("species", sorted(LOBE_SPECIES))
def test_every_lobe_has_a_name_and_a_source(species: Species) -> None:
    for lobe in lobe_table(species):
        assert lobe.name
        assert lobe.mass_fraction.source


def test_lobe_table_pig_raises() -> None:
    with pytest.raises(ValueError, match="pig"):
        lobe_table(Species.PIG)


def test_human_segment_order_follows_mise2014() -> None:
    """Mise 2014: the median volumes fall in the order VIII, VII, IV, V, III, VI, II, I."""
    by_token = {
        lobe.token: lobe.mass_fraction.value for lobe in lobe_table(Species.HUMAN)
    }
    order = ["sg8", "sg7", "sg4", "sg5", "sg3", "sg6", "sg2", "sg1"]
    shares = [by_token[t] for t in order]
    assert shares == sorted(shares, reverse=True)
    assert by_token["sg8"] == pytest.approx(0.261, abs=5e-4)


def test_rat_lobes_without_the_caudate_are_ninety_percent() -> None:
    """Madrahimov 2006: the left lateral, median and right lobes are the 90 % hepatectomy."""
    by_token = {
        lobe.token: lobe.mass_fraction.value for lobe in lobe_table(Species.RAT)
    }
    resected = by_token["lll"] + by_token["lml"] + by_token["rml"] + by_token["rl"]
    assert resected == pytest.approx(0.90, abs=1e-9)


def test_check_lobes_refuses_duplicate_token() -> None:
    with pytest.raises(ValueError, match="'rl'"):
        check_lobes([Lobe("rl", "a", _share(0.5)), Lobe("rl", "b", _share(0.5))])


def test_check_lobes_refuses_sum_not_one() -> None:
    with pytest.raises(ValueError, match=r"0\.8"):
        check_lobes([Lobe("rl", "a", _share(0.8))])


def test_check_lobes_refuses_non_identifier_token() -> None:
    with pytest.raises(ValueError, match="'1a'"):
        check_lobes([Lobe("1a", "a", _share(1.0))])


def test_check_lobes_refuses_an_empty_set() -> None:
    with pytest.raises(ValueError, match="at least one"):
        check_lobes([])


def test_check_lobes_refuses_an_upper_case_token() -> None:
    with pytest.raises(ValueError, match="'LLL'"):
        check_lobes([Lobe("LLL", "a", _share(1.0))])


def test_check_lobes_returns_the_tuple() -> None:
    lobes = check_lobes(lobe_table(Species.RAT))
    assert isinstance(lobes, tuple)
    assert len(lobes) == 5


def test_supply_fraction_defaults_to_mass_fraction() -> None:
    lobe = lobe_table(Species.RAT)[0]
    assert lobe.portal_share() == lobe.mass_fraction.value
    assert lobe.arterial_share() == lobe.mass_fraction.value


def test_a_given_supply_fraction_is_used() -> None:
    lobe = Lobe("rl", "right", _share(1.0), portal_fraction=_share(0.4))
    assert lobe.portal_share() == 0.4
    assert lobe.arterial_share() == 1.0


def test_check_lobes_refuses_supply_fractions_which_do_not_sum_to_one() -> None:
    with pytest.raises(ValueError, match="portal"):
        check_lobes(
            [
                Lobe("a", "a", _share(0.5), portal_fraction=_share(0.4)),
                Lobe("b", "b", _share(0.5), portal_fraction=_share(0.4)),
            ]
        )


def test_habr_has_the_measured_default_gain() -> None:
    """Lautt 1985: the buffer response compensates 25.5 +- 2.7 % of the lost portal flow."""
    assert HABR().gain.value == pytest.approx(0.255, abs=1e-9)
    assert HABR().gain.source == "Lautt1985"
