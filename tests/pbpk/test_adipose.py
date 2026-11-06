"""The adipose tissue, the body composition and the physiology of obesity."""

import logging
from itertools import pairwise
from pathlib import Path

import libsbml
import numpy as np
import pytest

from pkpdsim.impairment import BMI_CLASSES, CLASSIFICATIONS
from pkpdsim.pbpk import Adipose, Organ, PBPKSimulationExperiment, build, ids
from pkpdsim.pbpk.adipose import (
    BINDING_FACTORS,
    BMI_NORMAL,
    CARDIAC_OUTPUT,
    FLOW_FACTORS,
    obesity_factor,
    obesity_factors,
    obesity_protein_factors,
    visceral_share,
)
from pkpdsim.pbpk.age import AGE_REFERENCE
from pkpdsim.pbpk.physiology import body_composition, reference_individual
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Formula
from pkpdsim.sex import Sex
from pkpdsim.species import Species

from .conftest import definition, individuals, over_age, total_amount

FLOW_ORGANS = [o for o in Organ if o is not Organ.LUNG]
"""The organs whose flows sum to the cardiac output."""

WEIGHT_EFFECT = Formula(
    "-0.2 dimensionless * Cve_drg / (Cve_drg + 0.0002 mM)", U.dimensionless
)
"""An Emax effect of the drug on the body weight at the set point."""


@pytest.fixture(scope="module")
def adipose(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The simple definition with the body composition and the obesity factors."""
    return build(
        definition(sid="simple_adipose", adipose=Adipose()),
        tmp_path_factory.mktemp("adipose"),
    ).flat_path


@pytest.fixture(scope="module")
def without_obesity(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The same definition with the body composition but without the factors."""
    return build(
        definition(sid="simple_composition", adipose=Adipose(obesity=False)),
        tmp_path_factory.mktemp("composition"),
    ).flat_path


@pytest.fixture(scope="module")
def losing_weight(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The definition whose body weight follows an effect of the drug."""
    return build(
        definition(sid="simple_loss", adipose=Adipose(effect=WEIGHT_EFFECT)),
        tmp_path_factory.mktemp("loss"),
    ).flat_path


def _parameters(path: Path) -> set[str]:
    """The parameters of a built model."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    return {p.getId() for p in model.getListOfParameters()}


def test_a_definition_without_adipose_has_none_of_the_elements(
    tmp_path: Path,
) -> None:
    """The body composition and the factors are a choice of the build."""
    path = build(definition(sid="plain"), tmp_path).flat_path
    parameters = _parameters(path)
    composition = {ids.LBW, ids.FM, ids.FMI, ids.BFF, ids.F_VISCERAL, ids.VAT, ids.SAT}
    assert not parameters & composition
    assert not {p for p in parameters if p.startswith("fobes_")}


def test_the_composition_of_a_build(adipose: Path) -> None:
    """`LBW`, `FM`, `BFF`, `FMI`, `f_vat`, `VAT` and `SAT` of the individual."""
    r = individuals(adipose)(BW=92.9, HEIGHT=176.0)
    composition = body_composition(Sex.MALE, 92.9, 176.0)
    assert r.getValue(ids.FM) == pytest.approx(composition.fm, rel=1e-12)
    assert r.getValue(ids.LBW) == pytest.approx(92.9 - composition.adipose, rel=1e-12)
    assert r.getValue(ids.BFF) == pytest.approx(composition.fm / 92.9, rel=1e-12)
    assert r.getValue(ids.FMI) == pytest.approx(composition.fm / 1.76**2, rel=1e-12)
    share = visceral_share(Sex.MALE).curve.value(AGE_REFERENCE)
    fat = r.getValue(ids.organ_mass(Organ.ADIPOSE))
    assert r.getValue(ids.F_VISCERAL) == pytest.approx(share, rel=1e-12)
    assert r.getValue(ids.VAT) == pytest.approx(share * fat, rel=1e-12)
    assert r.getValue(ids.VAT) + r.getValue(ids.SAT) == pytest.approx(fat, rel=1e-12)


def test_the_visceral_share_rises_with_age_and_is_smaller_in_women() -> None:
    """Shen 2009: 0.112 of the men and 0.056 of the women at 38 years, rising."""
    for sex in Sex:
        share = visceral_share(sex)
        values = [share.curve.value(a) for a in range(18, 71, 2)]
        assert all(b > a for a, b in pairwise(values))
        assert share.curve.value(5.0) == share.curve.value(18.0)
        assert share.curve.value(95.0) == share.curve.value(70.0)
        assert share.sources[0] == "Shen2009"
    assert visceral_share(Sex.MALE).curve.value(38.0) == pytest.approx(0.112, abs=2e-3)
    assert visceral_share(Sex.FEMALE).curve.value(38.0) == pytest.approx(
        0.056, abs=2e-3
    )
    for a in (20.0, 35.0, 60.0):
        assert visceral_share(Sex.FEMALE).curve.value(a) < visceral_share(
            Sex.MALE
        ).curve.value(a)


def test_every_factor_of_obesity_is_one_up_to_a_body_mass_index_of_25() -> None:
    """A normal weight individual has the physiology of a build without the factors."""
    for sex in Sex:
        factors = dict(obesity_factors(sex)) | {
            f"{o.value}_{p}": f for (o, p), f in obesity_protein_factors(sex).items()
        }
        assert set(obesity_factors(sex)) == {
            CARDIAC_OUTPUT,
            *FLOW_FACTORS,
            *BINDING_FACTORS,
        }
        for key, factor in factors.items():
            for bmi in (10.0, 18.5, 22.0, BMI_NORMAL):
                assert factor.value(bmi) == pytest.approx(1.0), key
            assert factor.sources
            assert factor.notes


def test_the_factors_of_obesity_follow_the_evidence() -> None:
    """The size and the direction of every factor at a body mass index of 40."""
    male = obesity_factors(Sex.MALE)
    at_40 = {k: f.value(40.0) for k, f in male.items()}
    # the cardiac output of the library overshoots, the brain, muscle, skin,
    # gut and pancreas fractions fall, the adipose and heart fractions rise
    assert at_40[CARDIAC_OUTPUT] == pytest.approx(0.915, abs=5e-3)
    assert at_40["brain_flow"] == pytest.approx(0.8375, abs=5e-3)
    assert at_40["adipose_flow"] == pytest.approx(1.25, abs=5e-3)
    assert at_40["heart_flow"] > 1.25
    assert all(at_40[k] < 0.9 for k in ("muscle_flow", "skin_flow", "gut_flow"))
    assert at_40["albumin"] == pytest.approx(0.9449, abs=5e-3)
    assert at_40["agp"] == pytest.approx(1.231, abs=5e-3)
    proteins = {
        p: f.value(40.0) for (_o, p), f in obesity_protein_factors(Sex.MALE).items()
    }
    assert proteins["CYP3A4"] == pytest.approx(0.7765, abs=5e-3)
    assert proteins["CYP2E1"] > 1.6


def test_the_factors_of_the_build_at_the_reference_individual(adipose: Path) -> None:
    """Every `fobes_<id>` is exactly 1 for the reference man, who is not obese."""
    reference = reference_individual(Sex.MALE)
    r = individuals(adipose)(BW=reference.bw, HEIGHT=reference.height)
    factors = {p for p in _parameters(adipose) if p.startswith("fobes_")}
    assert factors == {
        obesity_factor(ids.CO),
        *(
            obesity_factor(ids.organ_flow(o))
            for os in FLOW_FACTORS.values()
            for o in os
        ),
        *(obesity_factor(k) for k in BINDING_FACTORS),
        obesity_factor("li_cyp3a4"),
    }
    for pid in factors:
        assert r.getValue(pid) == pytest.approx(1.0, rel=1e-12), pid


def test_the_obese_physiology_of_the_build(
    adipose: Path, without_obesity: Path
) -> None:
    """The factors lower the cardiac output and keep the brain and muscle flows."""
    height = 176.0
    bw = 40.0 * (height / 100.0) ** 2
    make = individuals(adipose)
    reference = make(BW=reference_individual(Sex.MALE).bw, HEIGHT=height)
    brain_reference = reference.getValue(ids.organ_flow(Organ.BRAIN))
    muscle_reference = reference.getValue(ids.organ_flow(Organ.MUSCLE))
    perfusion_reference = reference.getValue(
        ids.organ_flow(Organ.ADIPOSE)
    ) / reference.getValue(ids.organ_mass(Organ.ADIPOSE))
    obese = make(BW=bw, HEIGHT=height)
    plain = individuals(without_obesity)(BW=bw, HEIGHT=height)
    # the cardiac output stays the sum of the flows and is lower than without
    flows = sum(obese.getValue(ids.organ_flow(o)) for o in FLOW_ORGANS)
    assert obese.getValue(ids.CO) == pytest.approx(flows, rel=1e-12)
    assert obese.getValue(ids.CO) < 0.9 * plain.getValue(ids.CO)
    # the brain and the muscle keep the absolute flow of the reference individual
    assert obese.getValue(ids.organ_flow(Organ.BRAIN)) == pytest.approx(
        brain_reference, rel=0.02
    )
    assert obese.getValue(ids.organ_flow(Organ.MUSCLE)) == pytest.approx(
        muscle_reference, rel=0.02
    )
    assert plain.getValue(ids.organ_flow(Organ.BRAIN)) > 1.3 * brain_reference
    # the adipose flow rises and its perfusion per kg falls
    assert obese.getValue(ids.organ_flow(Organ.ADIPOSE)) > 1.1 * plain.getValue(
        ids.organ_flow(Organ.ADIPOSE)
    )
    fat = obese.getValue(ids.organ_mass(Organ.ADIPOSE))
    assert (
        obese.getValue(ids.organ_flow(Organ.ADIPOSE)) / fat < 0.6 * perfusion_reference
    )
    # the GFR needs no factor: it is that of the build without the factors
    assert obese.getValue(ids.GFR) == pytest.approx(plain.getValue(ids.GFR), rel=1e-12)
    # the binding and the enzyme follow their factors
    assert obese.getValue(ids.impairment_factor("albumin")) == pytest.approx(
        0.9449, abs=5e-3
    )
    assert obese.getValue(
        ids.impairment_activity("CYP3A4", Organ.LIVER)
    ) == pytest.approx(0.7765, abs=5e-3)


def test_the_body_weight_of_a_model_with_an_effect_is_a_state(
    losing_weight: Path,
) -> None:
    """`BW` approaches `BW0*(1 + f_bw)` and the physiology follows it."""
    r = individuals(losing_weight)(BW=105.0, HEIGHT=176.0, Ri_drg=5e-3)
    r.timeCourseSelections = [
        "time",
        ids.BW,
        ids.BW_BASE,
        ids.BW_SET,
        ids.BW_EFFECT,
        "BW_change",
        "BW_ratio",
        ids.organ_mass(Organ.ADIPOSE),
        ids.FM,
        ids.FFM,
        ids.CO,
        ids.GFR,
        obesity_factor(ids.CO),
    ]
    result = r.simulate(0.0, 52 * 7 * 24 * 60.0, 400)
    frame = {name: result[:, i] for i, name in enumerate(result.colnames)}
    assert np.allclose(frame[ids.BW_BASE], 105.0)
    assert frame[ids.BW][-1] < frame[ids.BW][0]
    assert np.all(np.diff(frame[ids.BW]) <= 0)
    assert frame[ids.BW][-1] > frame[ids.BW_SET][-1]
    assert -0.2 < frame[ids.BW_EFFECT][-1] < -0.15
    assert frame["BW_change"][-1] == pytest.approx(frame[ids.BW][-1] - 105.0, rel=1e-9)
    assert frame["BW_ratio"][-1] == pytest.approx(frame[ids.BW][-1] / 105.0, rel=1e-9)
    # the fat takes 0.68 to 0.73 of the change, the measured fraction of the
    # DXA substudies of the GLP-1 trials (0.709, SD 0.19, over 40 reports)
    change = frame[ids.BW][-1] - frame[ids.BW][0]
    fat = (
        frame[ids.organ_mass(Organ.ADIPOSE)][-1]
        - frame[ids.organ_mass(Organ.ADIPOSE)][0]
    )
    lean = frame[ids.FFM][-1] - frame[ids.FFM][0]
    assert fat / change == pytest.approx(0.70, abs=0.05)
    assert lean / change == pytest.approx(0.30, abs=0.05)
    # the physiology follows the weight
    assert frame[ids.CO][-1] < frame[ids.CO][0]
    assert frame[ids.GFR][-1] < frame[ids.GFR][0]
    assert frame[obesity_factor(ids.CO)][-1] > frame[obesity_factor(ids.CO)][0]


def test_the_mass_balance_of_a_model_whose_weight_changes(losing_weight: Path) -> None:
    """The dose is conserved while the volumes of the individual change."""
    r = individuals(losing_weight)(BW=105.0, HEIGHT=176.0, IVDOSE_drg=10.0)
    r.timeCourseSelections = ["time", *r.model.getFloatingSpeciesIds()]
    r.simulate(0.0, 8 * 7 * 24 * 60.0, 200)
    mr = r.getValue("Mr_drg")
    assert total_amount(r, "drg") + total_amount(r, "drg_m") == pytest.approx(
        10.0 / mr, rel=1e-4
    )


def test_an_effect_which_is_not_dimensionless_is_refused() -> None:
    """The effect is the fractional change of the body weight."""
    with pytest.raises(ValueError, match="dimensionless"):
        Adipose(effect=Formula("1 l", U.litre))


def test_a_rate_without_an_effect_and_a_rate_which_is_not_positive_are_refused() -> (
    None
):
    """The rate belongs to a model whose weight changes."""
    with pytest.raises(ValueError, match="needs an effect"):
        Adipose(k_bw=1e-6)
    with pytest.raises(ValueError, match="not positive"):
        Adipose(effect=WEIGHT_EFFECT, k_bw=0.0)


def test_a_laboratory_animal_has_no_adipose_model(tmp_path: Path) -> None:
    """The body composition and the obesity of the library are those of the human."""
    with pytest.raises(ValueError, match="body composition"):
        build(
            definition(sid="rat_adipose", species=Species.RAT, adipose=Adipose()),
            tmp_path,
        )


def test_the_classes_of_the_body_mass_index() -> None:
    """The WHO cut-offs, the class of normal weight the reference individual."""
    assert CLASSIFICATIONS[BMI_CLASSES.key] is BMI_CLASSES
    assert BMI_CLASSES.parameter == ids.BMI
    assert BMI_CLASSES.keys() == [
        "normal",
        "underweight",
        "overweight",
        "obesity_1",
        "obesity_2",
        "obesity_3",
    ]
    assert BMI_CLASSES["normal"].value is None
    assert [
        BMI_CLASSES[k].limits for k in ("overweight", "obesity_1", "obesity_3")
    ] == [
        (25.0, 30.0),
        (30.0, 35.0),
        (40.0, None),
    ]


def test_the_changes_of_a_class_of_the_body_mass_index() -> None:
    """They set the body weight at the height of the reference individual."""

    class Experiment(PBPKSimulationExperiment):
        """An experiment of the men."""

    reference = reference_individual(Sex.MALE)
    changes = Experiment.class_changes(BMI_CLASSES, "obesity_2")
    assert changes[ids.HEIGHT].magnitude == pytest.approx(reference.height)
    assert changes[ids.BW].magnitude == pytest.approx(
        37.5 * (reference.height / 100.0) ** 2
    )
    normal = Experiment.class_changes(BMI_CLASSES, "normal")
    assert normal[ids.BW].magnitude == pytest.approx(reference.bw)


@pytest.mark.uncached_build
def test_a_build_whose_weight_changes_logs_no_sbmlutils_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """`BW0` stays a constant with its initial assignment (`StartValue`).

    The body weight is a state, so every quantity downstream of it becomes an
    assignment rule; the body weight at the start must not follow it, and a
    variable parameter which nothing changes is a warning of sbmlutils.
    """
    with caplog.at_level(logging.WARNING, logger="sbmlutils"):
        build(
            definition(sid="quiet_loss", adipose=Adipose(effect=WEIGHT_EFFECT)),
            tmp_path,
        )
    assert [
        r.getMessage() for r in caplog.records if r.name.startswith("sbmlutils")
    ] == []


@over_age
def test_the_visceral_share_of_a_build_over_age_is_a_curve_of_the_age(
    tmp_path: Path,
) -> None:
    """A build over age reads `f_vat` from `AGE`, a build without is 35 years."""
    path = build(
        definition(sid="aged_adipose", adipose=Adipose(), age=True), tmp_path
    ).flat_path
    make = individuals(path)
    share = visceral_share(Sex.MALE).curve
    for age in (25.0, 45.0, 65.0):
        reference = reference_individual(Sex.MALE, age)
        r = make(AGE=age, BW=reference.bw, HEIGHT=reference.height)
        assert r.getValue(ids.F_VISCERAL) == pytest.approx(share.value(age), rel=1e-9)


def test_the_visceral_share_of_a_build_over_age_reads_the_age(tmp_path: Path) -> None:
    """The initial assignment of `f_vat` is the curve of `AGE` (no model is loaded)."""
    paths = build(definition(sid="aged_share", adipose=Adipose(), age=True), tmp_path)
    model = libsbml.readSBMLFromFile(str(paths.flat_path)).getModel()
    assignment = next(
        a
        for a in model.getListOfInitialAssignments()
        if a.getSymbol() == ids.F_VISCERAL
    )
    formula = libsbml.formulaToL3String(assignment.getMath())
    assert ids.AGE in formula
    assert str(visceral_share(Sex.MALE).curve.value(18.0)) in formula
