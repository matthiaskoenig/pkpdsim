"""Tests of the physiology over age: the curves, the model and the reference individuals."""

import itertools
import math
from pathlib import Path

import pytest

from pkpdsim.impairment import AGE_CLASSES
from pkpdsim.pbpk import Organ, Sex, build, ids, reference_individual
from pkpdsim.pbpk.age import (
    ADULT_AGE,
    AGE_REFERENCE,
    age_factor,
    factor_curves,
    flow_factors,
    flow_fraction_curves,
    protein_factor,
    protein_factors,
    reference_curves,
)
from pkpdsim.pbpk.curves import (
    Affine,
    Constant,
    Hill,
    Interpolation,
    Logarithmic,
    Product,
    Quotient,
    Relative,
    Spliced,
    pma,
)
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.physiology import _references, normal_gfr
from pkpdsim.pbpk.sex_differences import protein_sex_factor, sex_factor

from .conftest import definition, individuals
from .test_physiology import _assert_balance, _assert_flows

AGES = (0.0, 0.01, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 15.0, 18.0, 20.0, 27.5, 35.0, 42.0)
AGES += (50.0, 57.0, 65.0, 72.0, 80.0, 85.0, 90.0, 100.0)
"""The ages of the comparisons: every growth anchor, splice and the elderly."""


@pytest.fixture(scope="module")
def builds(tmp_path_factory: pytest.TempPathFactory) -> dict[Sex, Path]:
    """The simple definition built for each sex."""
    return {
        sex: build(
            definition(sex=sex, age=True), tmp_path_factory.mktemp(sex.value)
        ).flat_path
        for sex in Sex
    }


def _fractions(sex: Sex) -> dict[Organ, float]:
    """The flow fractions of the adults of a sex."""
    refs = _references(sex, load("physiology_human"))
    return {o: refs[ids.flow_fraction(o)].value for o in Organ if o is not Organ.LUNG}


def test_the_curves_are_their_formulas_at_their_points() -> None:
    """The forms of the curves evaluate as defined."""
    i = Interpolation((0.0, 1.0, 5.0), (1.0, 3.0, 2.0))
    assert [i.value(a) for a in (-1.0, 0.0, 0.5, 1.0, 3.0, 5.0, 9.0)] == [
        1.0,
        1.0,
        2.0,
        3.0,
        2.5,
        2.0,
        2.0,
    ]
    h = Hill(x50=47.7, n=3.4, variable="pma")
    assert h.value(0.0) == pytest.approx(40**3.4 / (40**3.4 + 47.7**3.4))
    assert pma(1.0) == pytest.approx(40 + 365.25 / 7)
    assert Hill(x50=2.0, n=1.0, base=0.5).value(2.0) == 0.75
    assert Hill(x50=1.0, n=1.0, x_min=1.0).value(0.0) == 0.5
    log = Logarithmic(33.746, 1.1287)
    assert log.value(0.0) == 33.746
    assert log.value(1.0) == pytest.approx(33.746 + 1.1287 * math.log(365.25))
    assert Relative(Constant(4.0), 35.0).value(3.0) == 1.0
    assert Affine(1.0, ((2.0, Constant(3.0)),)).value(0.0) == 7.0
    assert Quotient(Constant(3.0), Constant(4.0)).value(0.0) == 0.75
    assert Product((Constant(3.0), Constant(4.0))).value(0.0) == 12.0
    spliced = Spliced((1.0,), (Constant(2.0), Interpolation((1.0, 2.0), (2.0, 4.0))))
    assert spliced.value(0.5) == 2.0
    assert spliced.value(1.5) == 3.0


def test_spliced_curves_must_meet() -> None:
    """Curves which do not meet at a splice are refused."""
    with pytest.raises(ValueError, match=r"do not meet at the age 1\.0"):
        Spliced((1.0,), (Constant(2.0), Constant(3.0)))


@pytest.mark.parametrize("sex", list(Sex))
def test_every_factor_is_one_at_the_reference_age(sex: Sex) -> None:
    """The adult of ICRP 89 keeps its values: every factor of the age is 1 at 35 years."""
    for c in factor_curves(sex).values():
        assert c.curve.value(AGE_REFERENCE) == pytest.approx(1.0, rel=1e-12)
    f = flow_factors(sex)
    for curve in (f.brain, f.liver, f.kidney):
        assert curve.value(AGE_REFERENCE) == pytest.approx(1.0, rel=1e-12)
    adult = load("physiology_human")[sex.value]
    refs = reference_curves(sex)
    assert (
        refs[ids.ref(ids.BW)].curve.value(AGE_REFERENCE) == adult["bodyweight"]["value"]
    )
    assert refs[ids.ref(ids.CO)].curve.value(AGE_REFERENCE) == pytest.approx(
        adult["cardiac_output"]["value"], rel=1e-12
    )
    for organ, protein in protein_factors():
        curve = protein_factor(organ, protein)
        if curve is not None:
            assert curve.curve.value(AGE_REFERENCE) == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("sex", list(Sex))
def test_the_flow_fractions_sum_to_one_at_every_age(sex: Sex) -> None:
    """The fractions of the cardiac output close at every age."""
    curves = flow_fraction_curves(sex, _fractions(sex))
    for age in AGES:
        assert sum(c.value(age) for c in curves.values()) == pytest.approx(
            1.0, abs=1e-12
        )


@pytest.mark.parametrize("sex", list(Sex))
def test_the_reference_children_are_those_of_icrp89(sex: Sex) -> None:
    """The curves pass through the reference children of ICRP 89 Tables 2.8, 2.9 and 2.39."""
    children = load("children_human")
    refs = reference_curves(sex)
    for i, age in enumerate(children["ages"]):
        table = children[sex.value]
        assert refs[ids.ref(ids.BW)].curve.value(age) == table["bodyweight"]["value"][i]
        assert (
            refs[ids.ref(ids.CO)].curve.value(age)
            == table["cardiac_output"]["value"][i]
        )
        for organ, entry in table["organs"].items():
            if organ == "ad":
                continue
            pid = ids.ref(ids.organ_mass(Organ(organ)))
            assert refs[pid].curve.value(age) == pytest.approx(entry["value"][i])


def test_the_maturation_and_the_decline_of_the_gfr() -> None:
    """GFR per body surface area: 36 % at birth, 90 % at 1 year, the decline of the adults."""
    gfr = factor_curves(Sex.MALE)[ids.GFR].curve
    assert gfr.value(0.0) == pytest.approx(40**3.4 / (40**3.4 + 47.7**3.4), rel=1e-4)
    assert gfr.value(1.0) == pytest.approx(0.904, abs=0.002)
    assert gfr.value(40.0) == pytest.approx(1.0, abs=1e-4)
    assert gfr.value(80.0) == pytest.approx(0.9169 - 0.00678 * 30, rel=1e-3)
    women = factor_curves(Sex.FEMALE)[ids.GFR].curve
    assert women.value(80.0) == pytest.approx(0.9169 - 0.00866 * 30, rel=1e-3)
    assert normal_gfr(Sex.MALE, 80.0) == pytest.approx(
        normal_gfr(Sex.MALE) * gfr.value(80.0), rel=1e-12
    )


def test_the_brain_takes_half_of_the_cardiac_output_of_an_infant() -> None:
    """The brain share of the cardiac output, Chiron 1992 with the ICRP 89 masses and flows."""
    curves = flow_fraction_curves(Sex.MALE, _fractions(Sex.MALE))
    brain = curves[Organ.BRAIN]
    assert 0.45 < brain.value(1.0) < 0.6
    assert 0.25 < brain.value(5.0) < 0.35
    assert brain.value(AGE_REFERENCE) == pytest.approx(
        _fractions(Sex.MALE)[Organ.BRAIN]
    )


def test_the_proteins_of_the_age_and_the_sex() -> None:
    """The ontogeny of CYP3A4 and its female factor; a stable protein has none.

    The hepatic CYP3A4 and CYP1A2 are the in vivo ontogeny of Upreti 2016
    (digitized from their Figure 1), which is not monotone: the activity per
    amount of liver is low in the neonate, passes the adult value in the
    first year, peaks above it in the toddler and falls back to exactly 1 at
    20 years. The in vitro ontogeny of the other proteins rises to 1.
    """
    cyp3a4 = protein_factor(Organ.LIVER, "CYP3A4")
    assert cyp3a4 is not None
    curve = cyp3a4.curve
    assert curve.value(0.0) < 0.3 < curve.value(0.1)
    assert curve.value(1.0) == pytest.approx(1.687, abs=5e-3)
    assert curve.value(5.0) < curve.value(1.0)
    assert max(curve.value(a) for a in (0.5, 1.0, 2.0, 3.0)) > 1.6
    for age in (20.0, 35.0, 70.0):
        assert curve.value(age) == pytest.approx(1.0, rel=1e-12)
    cyp1a2 = protein_factor(Organ.LIVER, "CYP1A2")
    assert cyp1a2 is not None
    assert cyp1a2.curve.value(2.0) == pytest.approx(1.516, abs=5e-3)
    assert cyp1a2.curve.value(35.0) == pytest.approx(1.0, rel=1e-12)
    small = protein_factor(Organ.SMALL_INTESTINE, "CYP3A4")
    assert small is not None
    assert small.curve.value(0.0) < small.curve.value(5.0) < small.curve.value(35.0)
    assert protein_factor(Organ.LIVER, "SLCO1B1") is None
    assert protein_sex_factor(Sex.FEMALE, Organ.LIVER, "CYP3A4") is not None
    assert protein_sex_factor(Sex.MALE, Organ.LIVER, "CYP3A4") is None
    assert (Organ.KIDNEY, "SLC22A6") in protein_factors()


@pytest.mark.parametrize("sex", list(Sex))
def test_the_model_computes_the_curves(builds: dict[Sex, Path], sex: Sex) -> None:
    """Every curve of the model is its curve of python at every age, to 1e-12."""
    at = individuals(builds[sex])
    refs = reference_curves(sex)
    factors = factor_curves(sex)
    fractions = flow_fraction_curves(sex, _fractions(sex))
    cyp3a4 = protein_factor(Organ.LIVER, "CYP3A4")
    assert cyp3a4 is not None
    for age in AGES:
        r = at(AGE=age)
        for pid, c in refs.items():
            assert r[pid] == pytest.approx(c.curve.value(age), rel=1e-12), (pid, age)
        for key in (ids.GFR, ids.HCT, ids.FFM, "albumin", "agp"):
            assert r[age_factor(key)] == pytest.approx(
                factors[key].curve.value(age), rel=1e-12
            ), (key, age)
        for o, c in fractions.items():
            assert r[ids.flow_fraction(o)] == pytest.approx(c.value(age), rel=1e-12)
        assert r[age_factor("li_cyp3a4")] == pytest.approx(
            cyp3a4.curve.value(age), rel=1e-12
        )
        female = protein_sex_factor(sex, Organ.LIVER, "CYP3A4")
        factor = 1.0 if female is None else female.curve.value(age)
        assert r[ids.impairment_activity("CYP3A4", Organ.LIVER)] == pytest.approx(
            cyp3a4.curve.value(age) * factor, rel=1e-12
        )


@pytest.mark.parametrize("sex", list(Sex))
def test_the_reference_individual_of_every_age_balances(
    builds: dict[Sex, Path], sex: Sex
) -> None:
    """The reference individual of an age: its masses are its body weight, its flows its output."""
    at = individuals(builds[sex])
    for age in AGES:
        ref = reference_individual(sex, age)
        r = at(AGE=age, BW=ref.bw, HEIGHT=ref.height)
        assert r[ids.SV_TIS] == pytest.approx(1.0, rel=1e-12)
        assert r[ids.organ_mass(Organ.ADIPOSE)] > 0
        _assert_balance(r, sex)
        _assert_flows(r)
        for o in Organ:
            assert r[ids.organ_volume(o)] == pytest.approx(ref.volumes[o], rel=1e-12)
        assert r[ids.GFR_BASE] == pytest.approx(ref.gfr, rel=1e-12)
        assert r[ids.HCT] < 1 - r[ids.GFR] / r[ids.organ_flow(Organ.KIDNEY)]


def test_the_sex_factor_is_in_the_female_build_only(builds: dict[Sex, Path]) -> None:
    """The female build has the sex factor of CYP3A4, the male build none."""
    key = "li_cyp3a4"
    male = individuals(builds[Sex.MALE])()
    female = individuals(builds[Sex.FEMALE])()
    assert female[sex_factor(key)] == 1.2
    with pytest.raises(RuntimeError):
        male[sex_factor(key)]


def test_the_growth_meets_the_adult(builds: dict[Sex, Path]) -> None:
    """The children reach the adult with its aging at 20 years."""
    refs = reference_curves(Sex.MALE)
    for c in refs.values():
        below = c.curve.value(ADULT_AGE - 1e-9)
        assert c.curve.value(ADULT_AGE) == pytest.approx(below, rel=1e-6)


def test_the_age_classes() -> None:
    """The ICH age classes cover the ages from birth without gaps."""
    classes = list(AGE_CLASSES)
    assert classes[0].limits == (0.0, 28 / 365.25)
    for a, b in itertools.pairwise(classes):
        assert a.limits is not None and b.limits is not None
        assert a.limits[1] == b.limits[0]
        assert a.limits[0] <= a.value < a.limits[1]  # ty: ignore[unsupported-operator]
    assert AGE_CLASSES["adult"].value == AGE_REFERENCE
