"""Tests of the profiles over the lobular position and of the lobular weight."""

import math
from collections.abc import Callable
from itertools import pairwise

import pytest

from pkpdsim.zonation.profiles import (
    Constant,
    Exponential,
    Linear,
    Profile,
    Random,
    Sigmoid,
    Weights,
    lobular_weight,
    normalized,
    zone_positions,
)

# The normalization constants of the spt cell model (pkdb_models
# `tpm/spt/models/spt_liver.py`), where `protein = 0.5 * <constant> * <shape>`:
# the constant makes the protein amount of every pattern the same, so
# `<constant> * mean(<shape>) = 0.5` over the weight of the lobule.
SPT_LINEAR_INCREASE = 1.1406136
SPT_LINEAR_DECREASE = 0.8902508
SPT_SHARP_PERICENTRAL = 4.552802
SPT_SHARP_PERIPORTAL = 2.3939743


def _integral(
    f: Callable[[float], float], lo: float = 0.0, hi: float = 1.0, n: int = 2001
) -> float:
    h = (hi - lo) / n
    return sum(f(lo + (i + 0.5) * h) for i in range(n)) * h


def _mean_of(f: Callable[[float], float]) -> float:
    """The mean of `f` over the lobular weight, which integrates to 1."""
    return _integral(lambda x: f(x) * lobular_weight(x))


def _hill_pericentral(x: float) -> float:
    return 1.0 - (1.0 - x) ** 10 / ((1.0 - x) ** 10 + 0.25**10)


def _hill_periportal(x: float) -> float:
    return 1.0 - x**10 / (x**10 + 0.25**10)


def _sampled_mean(f: Callable[[float], float], n: int = 240) -> float:
    """The area mean of `f` of the position, sampled on the hexagon itself.

    An independent check of `lobular_weight`: a grid over one twelfth of the
    lobule (the wedge of the nearest portal triad at `(1, 0)`), the position
    the relative distance of the point.
    """
    apothem = math.sqrt(3.0) / 2.0
    num = den = 0.0
    for i in range(n):
        theta = (i + 0.5) * (math.pi / 6.0) / n
        r_max = apothem / math.cos(math.pi / 6.0 - theta)
        for j in range(n):
            r = r_max * (j + 0.5) / n
            d_pp = math.hypot(r * math.cos(theta) - 1.0, r * math.sin(theta))
            area = r * r_max
            num += f(d_pp / (d_pp + r)) * area
            den += area
    return num / den


def _weighted_mean(values: tuple[float, ...], n: int) -> float:
    weights = [lobular_weight(x) for x in zone_positions(n)]
    return sum(v * w for v, w in zip(values, weights, strict=True)) / sum(weights)


# ------------------------------------------------------------------------------------
# the lobular weight
# ------------------------------------------------------------------------------------


def test_the_lobular_weight_integrates_to_one() -> None:
    assert _integral(lobular_weight) == pytest.approx(1.0, rel=1e-6)


def test_the_lobular_weight_is_positive_inside_and_vanishes_at_the_vessels() -> None:
    assert lobular_weight(0.0) == 0.0
    assert lobular_weight(1.0) == 0.0
    assert all(lobular_weight(k / 100.0) > 0.0 for k in range(1, 100))


@pytest.mark.parametrize("x", [-0.001, 1.001, 2.0])
def test_a_position_outside_the_lobule_raises(x: float) -> None:
    with pytest.raises(ValueError, match="position"):
        lobular_weight(x)


def test_the_bisectors_cut_a_third_of_the_lobule() -> None:
    # the points closer to the centre than to the nearest triad are the inner
    # hexagon of the perpendicular bisectors, a third of the lobule, and the
    # weight at the bisector is 40/27
    assert _integral(lobular_weight, 0.5, 1.0) == pytest.approx(1.0 / 3.0, rel=1e-6)
    assert lobular_weight(0.5) == pytest.approx(40.0 / 27.0, rel=1e-9)


def test_the_lobular_weight_matches_the_hexagon_it_is_derived_from() -> None:
    assert _mean_of(lambda x: x) == pytest.approx(_sampled_mean(lambda x: x), rel=1e-3)
    assert _mean_of(_hill_periportal) == pytest.approx(
        _sampled_mean(_hill_periportal), rel=1e-2
    )
    assert _mean_of(_hill_pericentral) == pytest.approx(
        _sampled_mean(_hill_pericentral), rel=1e-2
    )


def test_the_lobular_weight_lies_between_the_uniform_and_the_disc_weight() -> None:
    # a uniform weight has the mean 0.5, the disc weight `2 (1 - x)` the mean 1/3
    assert 1.0 / 3.0 < _mean_of(lambda x: x) < 0.5


def test_the_lobular_weight_reproduces_the_spt_linear_normalization() -> None:
    # the hexagonal lobule gives 0.434041 and 0.565959 against the 0.438361
    # and 0.561640 of spt, 1.0 % and 0.8 % off; a uniform weight (0.5) is 14 %
    # off and the disc weight (1/3) 24 %
    assert _mean_of(lambda x: x) == pytest.approx(0.5 / SPT_LINEAR_INCREASE, rel=1.5e-2)
    assert _mean_of(lambda x: 1.0 - x) == pytest.approx(
        0.5 / SPT_LINEAR_DECREASE, rel=1.5e-2
    )


@pytest.mark.xfail(
    strict=True,
    reason="the sharp constants of spt are not the moments of one lobular weight: "
    "the hexagonal lobule is 24 % and 11 % off and no vessel radius, Hill exponent "
    "or threshold reconciles both; the spt TODO.md asks for the sharp patterns to "
    "be fixed",
)
def test_the_lobular_weight_reproduces_the_spt_sharp_normalization() -> None:
    assert _mean_of(_hill_pericentral) == pytest.approx(
        0.5 / SPT_SHARP_PERICENTRAL, rel=1e-3
    )
    assert _mean_of(_hill_periportal) == pytest.approx(
        0.5 / SPT_SHARP_PERIPORTAL, rel=1e-3
    )


# ------------------------------------------------------------------------------------
# the zones and the normalization
# ------------------------------------------------------------------------------------


def test_zone_positions_are_the_centres_of_the_zones() -> None:
    assert zone_positions(1) == (0.5,)
    assert zone_positions(4) == (0.125, 0.375, 0.625, 0.875)
    assert len(zone_positions(10)) == 10


@pytest.mark.parametrize("n", [0, -1])
def test_a_lobule_without_zones_raises(n: int) -> None:
    with pytest.raises(ValueError, match="at least one zone"):
        zone_positions(n)


def test_normalized_checks_the_zone_count() -> None:
    with pytest.raises(ValueError, match="10 zones"):
        normalized([1.0, 1.0], 10)


def test_a_profile_without_a_positive_mean_cannot_be_normalized() -> None:
    with pytest.raises(ValueError, match="mean"):
        normalized([0.0] * 10, 10)


# ------------------------------------------------------------------------------------
# the profiles
# ------------------------------------------------------------------------------------


@pytest.mark.parametrize("n", [1, 5, 10, 50])
@pytest.mark.parametrize(
    "profile",
    [
        Constant(),
        Linear(0.0, 1.0),
        Linear(4.0, 1.0),
        Sigmoid(0.25, 10.0),
        Sigmoid(0.25, 10.0, reverse=True),
        Exponential(1.0, 2.0),
    ],
)
def test_every_profile_has_weighted_mean_one(profile: Profile, n: int) -> None:
    values = profile.values(n, lobular_weight)
    assert len(values) == n
    assert _weighted_mean(values, n) == pytest.approx(1.0, rel=1e-9)


@pytest.mark.parametrize(
    "profile", [Weights([0.5, 1.0, 1.5, 2.0] * 2 + [1.0, 1.0]), Random(42, 0.3)]
)
def test_the_explicit_and_the_random_profile_have_weighted_mean_one(
    profile: Profile,
) -> None:
    assert _weighted_mean(profile.values(10, lobular_weight), 10) == pytest.approx(1.0)


def test_every_profile_satisfies_the_protocol() -> None:
    profiles = [
        Constant(),
        Linear(1.0, 2.0),
        Exponential(1.0, 1.0),
        Sigmoid(0.25, 10.0),
        Weights([1.0] * 10),
        Random(1, 0.1),
    ]
    assert all(isinstance(p, Profile) for p in profiles)


def test_gradient_shape_is_independent_of_the_zone_count() -> None:
    at10 = Sigmoid(0.25, 10.0).values(10, lobular_weight)
    at50 = Sigmoid(0.25, 10.0).values(50, lobular_weight)
    assert at50[0] == pytest.approx(at10[0], rel=0.1)
    assert at50[-1] == pytest.approx(at10[-1], abs=1e-3)


def test_the_declared_scale_does_not_matter() -> None:
    assert Constant(7.0).values(10, lobular_weight) == pytest.approx(
        Constant().values(10, lobular_weight)
    )
    assert Linear(2.0, 8.0).values(10, lobular_weight) == pytest.approx(
        Linear(1.0, 4.0).values(10, lobular_weight)
    )


def test_a_linear_profile_keeps_the_declared_gradient() -> None:
    values = Linear(1.0, 4.0).values(10, lobular_weight)
    x = zone_positions(10)
    assert values[-1] / values[0] == pytest.approx(
        (1.0 + 3.0 * x[-1]) / (1.0 + 3.0 * x[0])
    )
    assert all(b > a for a, b in pairwise(values))


def test_an_exponential_profile_decays_from_periportal_to_pericentral() -> None:
    values = Exponential(1.0, 2.0).values(10, lobular_weight)
    x = zone_positions(10)
    assert all(b < a for a, b in pairwise(values))
    assert values[-1] / values[0] == pytest.approx(math.exp(-2.0 * (x[-1] - x[0])))


def test_the_sigmoid_is_periportal_and_its_reverse_pericentral() -> None:
    # the sharp patterns of spt: `Sigmoid(0.25, 10)` is its periportal pattern
    # and the reverse its pericentral one, the shape mirrored in the position
    pp = Sigmoid(0.25, 10.0).values(10, lobular_weight)
    pc = Sigmoid(0.25, 10.0, reverse=True).values(10, lobular_weight)
    assert all(b <= a for a, b in pairwise(pp))
    assert all(b >= a for a, b in pairwise(pc))
    mirrored = pc[::-1]
    scale = mirrored[0] / pp[0]
    assert mirrored == pytest.approx([scale * v for v in pp])


def test_the_sigmoid_is_the_sharp_shape_of_spt() -> None:
    values = Sigmoid(0.25, 10.0).values(10, lobular_weight)
    raw = [_hill_periportal(x) for x in zone_positions(10)]
    scale = values[0] / raw[0]
    assert values == pytest.approx([scale * v for v in raw])


def test_weights_of_the_wrong_length_raise() -> None:
    with pytest.raises(ValueError, match="10 zones"):
        Weights([1.0, 1.0]).values(10, lobular_weight)


def test_equal_weights_equal_a_constant_profile() -> None:
    assert Weights([1.0] * 10).values(10, lobular_weight) == Constant().values(
        10, lobular_weight
    )


def test_the_random_profile_is_reproducible_and_spreads_with_its_cv() -> None:
    assert Random(42, 0.3).values(10, lobular_weight) == Random(42, 0.3).values(
        10, lobular_weight
    )
    assert Random(7, 0.3).values(10, lobular_weight) != Random(42, 0.3).values(
        10, lobular_weight
    )
    assert Random(42, 0.0).values(10, lobular_weight) == Constant().values(
        10, lobular_weight
    )
    spread = Random(42, 0.5).values(200, lobular_weight)
    narrow = Random(42, 0.1).values(200, lobular_weight)
    assert max(spread) - min(spread) > max(narrow) - min(narrow)


@pytest.mark.parametrize(
    "call",
    [
        lambda: Constant(0.0),
        lambda: Linear(0.0, 0.0),
        lambda: Linear(-1.0, 1.0),
        lambda: Exponential(0.0, 1.0),
        lambda: Sigmoid(0.0, 10.0),
        lambda: Sigmoid(1.0, 10.0),
        lambda: Sigmoid(0.25, 0.0),
        lambda: Weights([]),
        lambda: Weights([0.0, 0.0]),
        lambda: Weights([-1.0, 1.0]),
        lambda: Random(1, -0.1),
    ],
)
def test_a_declaration_which_is_no_profile_raises(call: Callable[[], Profile]) -> None:
    with pytest.raises(ValueError):
        call()
