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


def test_a_negative_supply_fraction_is_refused() -> None:
    """A share of an inflow is no more allowed to be negative than a share of the mass.

    A set whose supply fractions sum to 1 with a negative one among them
    gives a lobe a negative blood flow, which is a model of nothing: the
    check refuses it where it refuses a negative mass fraction.
    """
    lobes = [
        Lobe("a", "a", _share(0.5), portal_fraction=_share(-0.2)),
        Lobe("b", "b", _share(0.5), portal_fraction=_share(1.2)),
    ]
    with pytest.raises(ValueError, match=r"portal fraction .* is negative"):
        check_lobes(lobes)


def test_a_token_of_a_lobe_is_no_prefix_of_another() -> None:
    """The tokens join the zone tokens into the site `<lobe>_<zone>`.

    A site is resolved back to its lobe by its prefix, so a set with `a` and
    `a_b` in it would give the sites of one lobe to the other: the tokens
    carry no underscore.
    """
    with pytest.raises(ValueError, match="underscore"):
        check_lobes([Lobe("a", "a", _share(0.5)), Lobe("a_b", "a b", _share(0.5))])


def test_the_growth_and_the_tracking_time_say_what_they_are() -> None:
    """The kind of an element is the term requested for it, so it has to fit.

    `dVli_<lobe>` is the rate of change of a volume and `tau_track` a time;
    neither is a blood flow to a lobe nor a flow per volume, and the
    requested terms of the kinds go upstream to the ontology
    (`docs/pbpk/annotation.md`).
    """
    from pkpdsim.pbpk.annotation import KINDS

    requested = {kind: KINDS[kind].requested or "" for kind in KINDS}
    assert "volume" in requested["lobe_growth"]
    assert "time" in requested["remodelling_time"]
    for kind in ("lobe_growth", "remodelling_time"):
        assert "blood flow" not in requested[kind]
