"""Tests of the profiles over the lobular position and of the lobular weight."""

import math
from collections.abc import Callable, Sequence
from itertools import pairwise

import pytest

from pkpdsim.zonation.profiles import (
    Constant,
    Exponential,
    Linear,
    Profile,
    Random,
    Sigmoid,
    Weight,
    Weights,
    lobular_area,
    lobular_weight,
    normalized,
    zone_areas,
    zone_positions,
)

# The normalization constants of the spt cell model (pkdb_models
# `tpm/spt/models/spt_liver.py`), where `protein = 0.5 * <constant> * <shape>`:
# the constant makes the protein amount of every pattern the same, so
# `<constant> * mean(<shape>) = 0.5` over the area of the lobule.
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


def _weighted_mean(values: Sequence[float], weights: Sequence[float]) -> float:
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


@pytest.mark.parametrize("f", [lobular_weight, lobular_area])
@pytest.mark.parametrize("x", [-0.001, 1.001, 2.0])
def test_a_position_outside_the_lobule_raises(
    f: Callable[[float], float], x: float
) -> None:
    with pytest.raises(ValueError, match="position"):
        f(x)


def test_the_bisectors_cut_a_third_of_the_lobule() -> None:
    # the points closer to the centre than to the nearest triad are the inner
    # hexagon of the perpendicular bisectors, a third of the lobule, and the
    # weight at the bisector is 40/27
    assert _integral(lobular_weight, 0.5, 1.0) == pytest.approx(1.0 / 3.0, rel=1e-6)
    assert lobular_area(0.5) == pytest.approx(2.0 / 3.0, rel=1e-9)
    assert lobular_weight(0.5) == pytest.approx(40.0 / 27.0, rel=1e-9)


def test_the_lobule_area_is_the_integral_of_the_lobular_weight() -> None:
    assert lobular_area(0.0) == 0.0
    assert lobular_area(1.0) == 1.0
    for x in (0.1, 1.0 / (1.0 + math.sqrt(3.0)), 0.5, 0.8):
        assert lobular_area(x) == pytest.approx(
            _integral(lobular_weight, 0.0, x), rel=1e-4
        )
    assert all(b > a for a, b in pairwise([lobular_area(k / 20.0) for k in range(21)]))


def test_the_areas_of_the_zones_are_exact_and_not_the_midpoint_density() -> None:
    areas = zone_areas(10)
    assert sum(areas) == pytest.approx(1.0, rel=1e-12)
    assert areas == pytest.approx(
        [lobular_area(k / 10.0) - lobular_area((k - 1) / 10.0) for k in range(1, 11)]
    )
    assert areas[3] == pytest.approx(_integral(lobular_weight, 0.3, 0.4), rel=1e-4)
    # the fourth zone holds the peak of the weight at `1 / (1 + sqrt(3))`, and
    # the density at the centre of the zone overstates its area by 6.7 %
    assert lobular_weight(0.35) / (10.0 * areas[3]) == pytest.approx(1.067, rel=1e-2)


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
    # the area of the hexagonal lobule gives 0.434041 and 0.565959 against the
    # 0.438361 and 0.561640 of spt, -0.99 % and +0.77 % off, so the constants
    # of the linear patterns of spt were computed from a hexagonal lobule; a
    # uniform weight (0.5) is 14 % off and the disc weight (1/3) 24 %
    assert _mean_of(lambda x: x) == pytest.approx(0.5 / SPT_LINEAR_INCREASE, rel=1.5e-2)
    assert _mean_of(lambda x: 1.0 - x) == pytest.approx(
        0.5 / SPT_LINEAR_DECREASE, rel=1.5e-2
    )


@pytest.mark.xfail(
    strict=True,
    reason="no lobular weight reproduces the sharp constants of spt: the area of "
    "the hexagonal lobule, which carries its linear constants to 1 %, is -23.78 % "
    "off the sharp pericentral and -11.08 % off the sharp periportal constant, and "
    "no distance power, vessel radius, Hill exponent or threshold of the geometry "
    "meets both; `pkdb_models/tpm/spt/TODO.md` records those two patterns as too "
    "sharp and not symmetrical",
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


@pytest.mark.parametrize("f", [zone_positions, zone_areas])
@pytest.mark.parametrize("n", [0, -1])
def test_a_lobule_without_zones_raises(
    f: Callable[[int], tuple[float, ...]], n: int
) -> None:
    with pytest.raises(ValueError, match="at least one zone"):
        f(n)


def test_normalized_checks_the_zone_count() -> None:
    with pytest.raises(ValueError, match="10 zones"):
        normalized([1.0, 1.0], 10)


def test_a_weight_of_the_wrong_length_raises() -> None:
    with pytest.raises(ValueError, match="10 zones"):
        normalized([1.0] * 10, 10, [1.0, 2.0])


def test_a_profile_without_a_positive_mean_cannot_be_normalized() -> None:
    with pytest.raises(ValueError, match="mean"):
        normalized([0.0] * 10, 10)


# ------------------------------------------------------------------------------------
# the profiles
# ------------------------------------------------------------------------------------


def _weight_and_weights(kind: str, n: int) -> tuple[Weight, tuple[float, ...]]:
    """A weight of every kind the module takes, and the weights it stands for."""
    if kind == "even":
        return None, (1.0,) * n
    if kind == "areas":
        return zone_areas(n), zone_areas(n)
    return lobular_weight, tuple(lobular_weight(x) for x in zone_positions(n))


@pytest.mark.parametrize("n", [1, 5, 10, 50])
@pytest.mark.parametrize("kind", ["even", "areas", "density"])
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
def test_every_profile_has_weighted_mean_one(
    profile: Profile, kind: str, n: int
) -> None:
    weight, weights = _weight_and_weights(kind, n)
    values = profile.values(n, weight)
    assert len(values) == n
    assert _weighted_mean(values, weights) == pytest.approx(1.0, rel=1e-9)


def test_the_default_weight_is_the_even_split_of_the_zones() -> None:
    profile = Linear(0.0, 1.0)
    values = profile.values(10)
    assert sum(values) / 10.0 == pytest.approx(1.0, rel=1e-12)
    assert values == pytest.approx(profile.values(10, [1.0] * 10))
    assert values == pytest.approx(profile.values(10, lambda _: 1.0))
    assert values != pytest.approx(profile.values(10, zone_areas(10)))


def test_the_exact_areas_weight_the_zones_and_a_density_is_only_sampled() -> None:
    raw = [_hill_periportal(x) for x in zone_positions(10)]
    areas = zone_areas(10)
    midpoints = tuple(lobular_weight(x) for x in zone_positions(10))
    exact = Sigmoid(0.25, 10.0).values(10, areas)
    sampled = Sigmoid(0.25, 10.0).values(10, lobular_weight)
    assert exact == pytest.approx([v / _weighted_mean(raw, areas) for v in raw])
    assert sampled == pytest.approx([v / _weighted_mean(raw, midpoints) for v in raw])
    # the lobular weighting is the sequence of the exact areas; a callable is
    # sampled at the centres of the zones, whatever it is, which for the
    # density of the lobule shifts this profile by 1.48 %
    assert max(abs(a / b - 1.0) for a, b in zip(exact, sampled, strict=True)) == (
        pytest.approx(0.0148, rel=1e-2)
    )
    assert sampled == pytest.approx(
        Sigmoid(0.25, 10.0).values(10, lambda x: lobular_weight(x))
    )


@pytest.mark.parametrize(
    "profile", [Weights([0.5, 1.0, 1.5, 2.0] * 2 + [1.0, 1.0]), Random(42, 0.3)]
)
def test_the_explicit_and_the_random_profile_have_weighted_mean_one(
    profile: Profile,
) -> None:
    areas = zone_areas(10)
    assert _weighted_mean(profile.values(10, areas), areas) == pytest.approx(1.0)
    assert sum(profile.values(10)) / 10.0 == pytest.approx(1.0)


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


def test_the_midpoint_shares_are_not_the_areas() -> None:
    """`lobular_weight` at the centre is a density: it misses the area by a few percent.

    The peak zone (the fourth of ten) is the worst, 5.6 % (measured, relative to the area); the map
    must not be drawn from the midpoint rule.
    """
    weights = [lobular_weight(x) for x in zone_positions(10)]
    midpoint = [w / sum(weights) for w in weights]
    deviation = [
        abs(m / a - 1.0) for m, a in zip(midpoint, zone_areas(10), strict=True)
    ]
    assert max(deviation) == pytest.approx(0.0557, abs=1e-3)
    assert deviation.index(max(deviation)) == 3
    assert max(deviation) > 0.05
