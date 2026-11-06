"""The differences of the sexes: the factors of the female build over age."""

from pathlib import Path

import libsbml
import pytest

from pkpdsim.pbpk import Organ, Sex, build, ids, reference_individual
from pkpdsim.pbpk.age import AGE_REFERENCE, factor_curves, flow_fraction_curves
from pkpdsim.pbpk.curves import Interpolation
from pkpdsim.pbpk.data import load
from pkpdsim.pbpk.physiology import _references
from pkpdsim.pbpk.sex_differences import (
    ALBUMIN,
    BRAIN_FLOW,
    CARDIAC_OUTPUT,
    COLON_TRANSIT,
    GASTRIC_EMPTYING_FED,
    protein_sex_factor,
    protein_sex_factors,
    sex_curves,
    sex_factor,
)

from .conftest import definition, individuals

AGES = (0.0, 1.0, 5.0, 10.0, 12.0, 14.0, 15.0, 16.0, 18.0, 20.0, 35.0, 45.0)
AGES += (50.0, 55.0, 60.0, 65.0, 80.0, 100.0)
"""The ages of the comparisons: the children, puberty, the adults, the menopause."""


@pytest.fixture(scope="module")
def builds(tmp_path_factory: pytest.TempPathFactory) -> dict[Sex, Path]:
    """The simple definition built for each sex."""
    return {
        sex: build(definition(sex=sex), tmp_path_factory.mktemp(sex.value)).flat_path
        for sex in Sex
    }


def test_the_male_build_is_the_reference() -> None:
    """The male build has no factor of the sex; the female build has every one."""
    assert sex_curves(Sex.MALE) == {}
    assert set(sex_curves(Sex.FEMALE)) == {
        CARDIAC_OUTPUT,
        BRAIN_FLOW,
        ALBUMIN,
        GASTRIC_EMPTYING_FED,
        COLON_TRANSIT,
    }
    for organ, protein in protein_sex_factors():
        assert protein_sex_factor(Sex.MALE, organ, protein) is None
        assert protein_sex_factor(Sex.FEMALE, organ, protein) is not None


def test_every_factor_is_one_in_the_children() -> None:
    """The differences arise at puberty: every factor of a curve is 1 until 9 years."""
    curves = [c for k, c in sex_curves(Sex.FEMALE).items() if k != ALBUMIN]
    curves += [
        c
        for organ, protein in protein_sex_factors()
        if (c := protein_sex_factor(Sex.FEMALE, organ, protein)) is not None
    ]
    for c in curves:
        assert c.sources
        for age in (0.0, 1.0, 5.0, 9.0):
            assert c.curve.value(age) == pytest.approx(1.0, abs=1e-12)


def _nhanes(sex: Sex) -> Interpolation:
    """The albumin of NHANES of a sex over age."""
    rows = load("blood_nhanes")["albumin"][sex.value]
    return Interpolation(tuple(rows["age"]), tuple(rows["value"]))


def test_the_albumin_of_the_women_is_that_of_nhanes_relative_to_the_men() -> None:
    """The female albumin over the male is the ratio of NHANES at every age from 12 years."""
    female = factor_curves(Sex.FEMALE)["albumin"].curve
    male = factor_curves(Sex.MALE)["albumin"].curve
    sex = sex_curves(Sex.FEMALE)[ALBUMIN].curve
    women, men = _nhanes(Sex.FEMALE), _nhanes(Sex.MALE)
    for age in (12.0, 15.0, 18.0, 25.0, 35.0, 42.0, 55.0, 60.0, 70.0, 75.0):
        model = female.value(age) * sex.value(age) / male.value(age)
        assert model == pytest.approx(women.value(age) / men.value(age), rel=1e-12)
    assert sex.value(AGE_REFERENCE) == pytest.approx(41.65 / 44.32, rel=1e-3)
    # the children before puberty have the albumin of the boys
    for age in (0.0, 1.0, 5.0, 10.0):
        model = female.value(age) * sex.value(age) / male.value(age)
        assert model == pytest.approx(1.0, abs=0.005)


def _fractions(sex: Sex) -> dict[Organ, float]:
    """The flow fractions of the adults of a sex."""
    refs = _references(sex, load("physiology_human"))
    return {o: refs[ids.flow_fraction(o)].value for o in Organ if o is not Organ.LUNG}


def test_the_brain_of_the_women_takes_more_of_the_cardiac_output() -> None:
    """The female brain fraction is the adult one times the factors; the sum stays 1."""
    fractions = _fractions(Sex.FEMALE)
    curves = flow_fraction_curves(Sex.FEMALE, fractions)
    sex = sex_curves(Sex.FEMALE)[BRAIN_FLOW].curve
    assert curves[Organ.BRAIN].value(AGE_REFERENCE) == pytest.approx(
        fractions[Organ.BRAIN] * sex.value(AGE_REFERENCE), rel=1e-12
    )
    assert sex.value(AGE_REFERENCE) > 1.2
    assert sex.value(80.0) < sex.value(AGE_REFERENCE)
    for age in AGES:
        assert sum(c.value(age) for c in curves.values()) == pytest.approx(
            1.0, abs=1e-12
        )


def test_the_model_computes_the_factors(builds: dict[Sex, Path]) -> None:
    """Every factor of the female build is its curve of python at every age, to 1e-12."""
    at = individuals(builds[Sex.FEMALE])
    curves = {sex_factor(k): c.curve for k, c in sex_curves(Sex.FEMALE).items()}
    cyp3a4 = protein_sex_factor(Sex.FEMALE, Organ.LIVER, "CYP3A4")
    assert cyp3a4 is not None
    curves[sex_factor("li_cyp3a4")] = cyp3a4.curve
    for age in AGES:
        r = at(AGE=age)
        for pid, curve in curves.items():
            assert r[pid] == pytest.approx(curve.value(age), rel=1e-12), (pid, age)


def test_the_male_build_has_no_factor_of_the_sex(builds: dict[Sex, Path]) -> None:
    """No parameter `fsex_` is in the male build."""
    model = libsbml.readSBMLFromFile(str(builds[Sex.MALE])).getModel()
    assert not [
        p.getId() for p in model.getListOfParameters() if p.getId().startswith("fsex_")
    ]


def _sizes(age: float) -> list[tuple[float, float]]:
    """Body weights and heights of an age: between the reference individuals, of adults more."""
    refs = [reference_individual(sex, age) for sex in Sex]
    sizes = [
        (sum(r.bw for r in refs) / 2, sum(r.height for r in refs) / 2),
        *((r.bw, r.height) for r in refs),
    ]
    if age >= 18.0:
        sizes += [(55.0, 160.0), (90.0, 170.0), (100.0, 190.0)]
    return sizes


def test_women_and_men_of_the_same_size_have_the_same_cardiac_index(
    builds: dict[Sex, Path],
) -> None:
    """At the same age, body weight and height the cardiac output of the sexes is the same.

    Within 1.5 % from 12 years up to a body mass index of 31: the measured
    cardiac index of women and men of the same size is the same (Williams 1994,
    Carlsson 2012); the factor removes the 9 % more of the scaling with the
    lean tissue (the difference of the fat-free mass of the sexes falls with
    the body mass index, 2.3 % at 120 kg and 175 cm).
    """
    women, men = individuals(builds[Sex.FEMALE]), individuals(builds[Sex.MALE])
    for age in (12.0, 15.0, 18.0, 20.0, 35.0, 50.0, 80.0):
        for bw, height in _sizes(age):
            female = women(AGE=age, BW=bw, HEIGHT=height)[ids.CO]
            male = men(AGE=age, BW=bw, HEIGHT=height)[ids.CO]
            assert female / male == pytest.approx(1.0, abs=0.015), (age, bw, height)


def test_the_women_empty_a_meal_slower(builds: dict[Sex, Path]) -> None:
    """The fed transit of the stomach of a woman of 35 years is 1.2 times; fasted alike."""
    women, men = individuals(builds[Sex.FEMALE]), individuals(builds[Sex.MALE])
    for age, fed, ratio in (
        (35.0, 1.0, 1.2),
        (35.0, 0.0, 1.0),
        (8.0, 1.0, 1.0),
        (60.0, 1.0, 1.0),
    ):
        female = women(AGE=age, FED=fed)[ids.T_STOMACH]
        male = men(AGE=age, FED=fed)[ids.T_STOMACH]
        assert female / male == pytest.approx(ratio, rel=1e-12), (age, fed)


COLON = ("cae", "asc", "tra", "des", "sig", "rec")
"""The segments of the large intestine."""


def _rates(path: Path) -> list[str]:
    """The rate laws of the reactions of a built model."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    return [
        libsbml.formulaToL3String(r.getKineticLaw().getMath())
        for r in model.getListOfReactions()
    ]


def test_the_colonic_transit_of_the_women_is_slower(builds: dict[Sex, Path]) -> None:
    """Every transit of a segment of the large intestine reads the factor of the sex."""
    pid = sex_factor(COLON_TRANSIT)
    female = _rates(builds[Sex.FEMALE])
    colon = [r for r in female if any(f"T_{seg}" in r for seg in COLON)]
    assert colon
    assert all(pid in r for r in colon)
    assert all(any(f"T_{seg}" in r for seg in COLON) for r in female if pid in r)
    assert not [r for r in _rates(builds[Sex.MALE]) if pid in r]
    assert sex_curves(Sex.FEMALE)[COLON_TRANSIT].curve.value(35.0) == pytest.approx(
        1.333
    )
