"""Tests of the virtual populations: the draws of the factors and the scan of a sex."""

from pathlib import Path

import numpy as np
import pytest

from pkpdsim.pbpk import (
    Organ,
    PBPKSimulationExperiment,
    Sex,
    VirtualPopulation,
    build,
    ids,
    virtual_population,
)
from pkpdsim.pbpk.populations import population

from .conftest import definition, individuals


def test_the_draws_are_reproducible() -> None:
    """The same seed draws the same individuals, another seed others."""
    a = virtual_population("white_american_nhanes", 50, seed=1)
    b = virtual_population("white_american_nhanes", 50, seed=1)
    c = virtual_population("white_american_nhanes", 50, seed=2)
    assert a == b
    assert a != c
    assert isinstance(a, VirtualPopulation)
    assert len(a) == 50


def test_the_covariates_follow_the_population_of_the_age() -> None:
    """Height and log body weight have the moments and the correlation of their age."""
    vp = virtual_population(
        "white_american_nhanes", 4000, ages=(60.0, 60.0 + 1e-9), female=0.0
    )
    heights = np.array([i.height for i in vp.individuals])
    log_bw = np.log([i.bw for i in vp.individuals])
    d = population("white_american_nhanes").anthropometry[Sex.MALE].at(60.0)
    assert heights.mean() == pytest.approx(d.height_mean, rel=0.005)
    assert heights.std() == pytest.approx(d.height_sd, rel=0.06)
    assert log_bw.mean() == pytest.approx(d.ln_bw_mean, abs=0.01)
    assert log_bw.std() == pytest.approx(d.ln_bw_sd, rel=0.06)
    assert d.rho is not None
    assert np.corrcoef(heights, log_bw)[0, 1] == pytest.approx(d.rho, abs=0.05)
    # truncated at three standard deviations
    assert np.all(np.abs(heights - d.height_mean) <= 3 * d.height_sd + 1e-9)


def test_the_factors_combine() -> None:
    """Sex by its share, ages of a class, a fixed and a drawn severity."""
    vp = virtual_population(
        "black_american_nhanes",
        400,
        ages="old_old",
        female=0.25,
        conditions={ids.CP_SEVERITY: (1.0, 3.0), ids.HF_SEVERITY: 2.0},
    )
    women = vp.of_sex(Sex.FEMALE)
    assert 0.15 < len(women) / len(vp) < 0.35
    assert all(75.0 <= i.age < 85.0 for i in vp.individuals)
    assert all(1.0 <= i.conditions[ids.CP_SEVERITY] <= 3.0 for i in vp.individuals)
    assert all(i.conditions[ids.HF_SEVERITY] == 2.0 for i in vp.individuals)
    values = vp.values(Sex.FEMALE)
    assert set(values) == {
        ids.AGE,
        ids.BW,
        ids.HEIGHT,
        ids.CP_SEVERITY,
        ids.HF_SEVERITY,
    }
    assert all(len(v) == len(women) for v in values.values())
    assert women[0].changes()[ids.BW] == women[0].bw


def test_the_problems_of_a_population_are_refused() -> None:
    """No individual, a share outside of 0 to 1, a population without distributions or ages."""
    with pytest.raises(ValueError, match="needs individuals"):
        virtual_population("white_american_nhanes", 0)
    with pytest.raises(ValueError, match="share of women"):
        virtual_population("white_american_nhanes", 10, female=1.5)
    with pytest.raises(ValueError, match="no distribution of the covariates"):
        virtual_population("european_icrp", 10)
    with pytest.raises(ValueError, match="outside of the ages"):
        virtual_population("east_asian_tanaka", 10, ages=(5.0, 10.0))
    with pytest.raises(ValueError, match="Unknown condition"):
        virtual_population("white_american_nhanes", 10, conditions={"BW": 3.0})
    with pytest.raises(ValueError, match=r"outside of 0\.0 to 3\.0"):
        virtual_population(
            "white_american_nhanes", 10, conditions={ids.CP_SEVERITY: 4.0}
        )
    vp = virtual_population("white_american_nhanes", 5, female=0.0)
    with pytest.raises(ValueError, match="no individual of the sex female"):
        vp.values(Sex.FEMALE)


def test_the_children_of_a_population_grow() -> None:
    """The children of a class have the size of their age."""
    infants = virtual_population("mexican_american_nhanes", 300, ages="infant")
    children = virtual_population("mexican_american_nhanes", 300, ages="child")
    mean = lambda vp: np.mean([i.height for i in vp.individuals])  # noqa: E731
    assert mean(infants) < 90.0 < mean(children) < 150.0


def test_the_dimension_of_a_sex(tmp_path: Path) -> None:
    """The dimension of a scan has the covariates of the individuals of the sex in the units of the model."""
    vp = virtual_population("white_american_nhanes", 20, ages="child", seed=3)
    dimension = PBPKSimulationExperiment.population_dimension(vp, Sex.MALE)
    men = vp.of_sex(Sex.MALE)
    values = dimension.values
    assert values[ids.AGE].m_as("yr").tolist() == [i.age for i in men]
    assert values[ids.BW].m_as("kg").tolist() == [i.bw for i in men]
    # every child of the draw is a valid individual of the male build
    at = individuals(build(definition(age=True), tmp_path).flat_path)
    for i in men:
        r = at(**i.changes())
        assert r[ids.organ_mass(Organ.ADIPOSE)] > 0
        assert r[ids.AGE] == i.age
