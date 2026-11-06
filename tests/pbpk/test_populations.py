"""Tests of the populations: reference individuals and covariate distributions."""

import dataclasses
from pathlib import Path

import libsbml
import pytest
from sbmlsim.simulation.sampling import parameter_distributions
from sbmlutils.factory import Model, Uncertainty, UncertParameter, create_model

from pkpdsim.pbpk import MODEL_UNITS, KpMethod, U, Value, build, ids, sources
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.physiology import Sex, physiology_objects
from pkpdsim.pbpk.populations import (
    DEFAULT_POPULATION,
    POPULATIONS,
    Anthropometry,
    population,
)

from .conftest import definition, individual, simple_drug
from .test_physiology import _assert_balance, _assert_flows


def _physiology(path: Path, sex: Sex, key: str) -> Path:
    """The physiology of a sex and population as SBML model."""
    model = Model(
        sid=f"physiology_{key}_{sex}",
        units=U,
        model_units=MODEL_UNITS,
        objects=physiology_objects(sex, population=key),
    )
    create_model(model, filepath=path)
    return path


def test_population_keys() -> None:
    """The populations whose primary data were read; ICRP 89 is the default."""
    assert set(POPULATIONS) == {
        "european_icrp",
        "white_american_nhanes",
        "black_american_nhanes",
        "mexican_american_nhanes",
        "east_asian_tanaka",
    }
    assert DEFAULT_POPULATION == "european_icrp"
    assert definition().population == DEFAULT_POPULATION
    for key, p in POPULATIONS.items():
        assert p.key == key
        assert p.name
        assert p.source in sources()
        for sex in Sex:
            assert set(p.reference[sex]) == {ids.BW, ids.HEIGHT, ids.AGE}


def test_european_reproduces_defaults() -> None:
    """The default population is the ICRP 89 reference individual of physiology_human."""
    p = population("european_icrp")
    assert p.reference[Sex.MALE][ids.BW].value == 73.0
    assert p.reference[Sex.MALE][ids.HEIGHT].value == 176.0
    assert p.reference[Sex.FEMALE][ids.BW].value == 60.0
    assert p.reference[Sex.FEMALE][ids.HEIGHT].value == 163.0
    assert p.anthropometry == {}
    for sex in Sex:
        default = {str(o.sid): o for o in physiology_objects(sex)}
        explicit = {
            str(o.sid): o for o in physiology_objects(sex, population="european_icrp")
        }
        assert default.keys() == explicit.keys()
        for pid in (ids.BW, ids.HEIGHT):
            assert default[pid].value == explicit[pid].value  # ty: ignore[unresolved-attribute]
    # the age is a covariate of a build over age only (`PBPKDefinition.age`)
    assert ids.AGE not in default
    over = {str(o.sid): o for o in physiology_objects(Sex.MALE, over_age=True)}
    age = over[ids.AGE]
    assert age.value == 35.0  # ty: ignore[unresolved-attribute]
    assert "the reference individual of the sex" in str(age.notes)


@pytest.mark.parametrize("key", sorted(POPULATIONS))
@pytest.mark.parametrize("sex", list(Sex))
def test_population_reference_sums(tmp_path: Path, key: str, sex: Sex) -> None:
    """The reference individual of every population passes the sums of the physiology.

    Its covariates are the defaults of the model; organ masses and blood are
    the body weight, the flows the cardiac output, adipose stays positive.
    """
    path = _physiology(tmp_path / "physiology.xml", sex, key)
    r = individual(path)
    reference = population(key).reference[sex]
    for pid in (ids.BW, ids.HEIGHT):
        assert r[pid] == reference[pid].value, pid
    assert r[ids.organ_mass(Organ.ADIPOSE)] > 0
    _assert_balance(r, sex)
    _assert_flows(r)


def test_population_distributions() -> None:
    """The NHANES populations describe every age, the East Asian the adults of 20-49 years."""
    for key, p in POPULATIONS.items():
        if key == "european_icrp":
            continue
        for sex in Sex:
            a = p.anthropometry[sex]
            assert isinstance(a, Anthropometry)
            if key.endswith("_nhanes"):
                assert a.age_range == (0.0, 100.0)
                assert a.rho is not None
                adult = a.at(35.0)
                reference = p.reference[sex]
                # the kernel moments of 35 years are close to the mean of 20-49 years
                assert adult.height_mean == pytest.approx(
                    reference[ids.HEIGHT].value, rel=0.01
                )
                assert adult.bw_mean == pytest.approx(reference[ids.BW].value, rel=0.08)
                # the children grow
                assert (
                    a.at(1.0).height_mean < a.at(10.0).height_mean < adult.height_mean
                )
            else:
                assert a.age_range == (20.0, 49.0)
                assert a.rho is None


def test_tanaka_reference_values() -> None:
    """Tanaka 1996 Tables 9 and 10, Asian Reference Man of 20 to 49 years."""
    p = population("east_asian_tanaka")
    male, female = p.reference[Sex.MALE], p.reference[Sex.FEMALE]
    assert (male[ids.BW].value, male[ids.HEIGHT].value) == (60.0, 169.97)
    assert (female[ids.BW].value, female[ids.HEIGHT].value) == (51.28, 155.04)
    bw = p.anthropometry[Sex.MALE].at(35.0)
    assert bw.bw_mean == pytest.approx(63.63, rel=1e-12)
    assert bw.height_sd == 5.68


def test_an_age_outside_of_the_data_is_refused() -> None:
    """The East Asian distribution describes the adults of 20 to 49 years only."""
    with pytest.raises(ValueError, match=r"outside of the ages 20\.0 to 49\.0"):
        population("east_asian_tanaka").anthropometry[Sex.MALE].at(5.0)


def test_population_chosen_at_build(tmp_path: Path) -> None:
    """The population of a definition gives the default individual of the model."""
    d = definition(population="east_asian_tanaka", sex=Sex.FEMALE)
    r = individual(build(d, tmp_path).flat_path)
    assert r[ids.BW] == 51.28
    assert r[ids.HEIGHT] == 155.04


def test_unknown_population_raises() -> None:
    """An unknown population is a problem of the definition."""
    with pytest.raises(ValueError, match=r"population 'japanese' is unknown"):
        definition(population="japanese").validate()
    with pytest.raises(ValueError, match="Unknown population"):
        population("japanese")


def test_no_distrib_from_population(tmp_path: Path) -> None:
    """A built model has no distrib uncertainty but those of its fitted values.

    The covariate distributions of a population are library data:
    `run_simulation --uncertainty` draws every parameter with a distrib
    uncertainty and must draw only the uncertainty of the fits.
    """
    fitted = Value(
        1.2,
        U.dimensionless,
        "fit",
        uncertainties=[
            Uncertainty(
                uncertParameters=[
                    UncertParameter(
                        type=libsbml.DISTRIB_UNCERTTYPE_DISTRIBUTION,
                        math="normal(1.2, 0.1)",
                    )
                ]
            )
        ],
    )
    method = KpMethod(
        method="rodgers_rowland",
        logp=2.0,
        pka=(),
        compound="neutral",
        ftissue=1.0,
        factor=fitted,
    )
    drug = dataclasses.replace(simple_drug(), distribution=method)
    for key in ("white_american_nhanes", "east_asian_tanaka"):
        d = definition(substances=[drug], chemistry=[], population=key)
        result = build(d, tmp_path / key)
        assert set(parameter_distributions(result.flat_path)) == {ids.fkp("drg")}
