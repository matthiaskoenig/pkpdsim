"""Tests of the populations: reference individuals and covariate distributions."""

import dataclasses
from pathlib import Path

import libsbml
import pytest
from sbmlutils.factory import Model, Uncertainty, UncertParameter, create_model

from pkpdsim.pbpk import MODEL_UNITS, KpMethod, U, Value, build, ids, sources
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.physiology import Sex, physiology_objects
from pkpdsim.pbpk.populations import (
    DEFAULT_POPULATION,
    POPULATIONS,
    Distribution,
    population,
)
from pkpdsim.uncertainty import parameter_distributions

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
    assert p.covariates == {Sex.MALE: {}, Sex.FEMALE: {}}
    for sex in Sex:
        default = {str(o.sid): o for o in physiology_objects(sex)}
        explicit = {
            str(o.sid): o for o in physiology_objects(sex, population="european_icrp")
        }
        assert default.keys() == explicit.keys()
        for pid in (ids.BW, ids.HEIGHT, ids.AGE):
            assert default[pid].value == explicit[pid].value  # ty: ignore[unresolved-attribute]
    age = {str(o.sid): o for o in physiology_objects(Sex.MALE)}[ids.AGE]
    notes = str(age.notes)
    assert "the population gives no value" in notes
    assert "reference individual" not in notes


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
    for pid in (ids.BW, ids.HEIGHT, ids.AGE):
        assert r[pid] == reference[pid].value, pid
    assert r[ids.organ_mass(Organ.ADIPOSE)] > 0
    _assert_balance(r, sex)
    _assert_flows(r)


def test_population_distributions() -> None:
    """Body weight is lognormal, height normal, around the reference of NHANES."""
    for key, p in POPULATIONS.items():
        if key == "european_icrp":
            continue
        for sex in Sex:
            covariates = p.covariates[sex]
            assert set(covariates) == {ids.BW, ids.HEIGHT}, key
            bw, height = covariates[ids.BW], covariates[ids.HEIGHT]
            assert isinstance(bw, Distribution)
            assert (bw.kind, bw.unit) == ("lognormal", U.kg)
            assert (height.kind, height.unit) == ("normal", U.cm)
            assert 0 < bw.sd < bw.mean
            assert 0 < height.sd < 0.1 * height.mean
            assert height.mean == pytest.approx(
                p.reference[sex][ids.HEIGHT].value, rel=1e-12
            )
            if key.endswith("_nhanes"):
                assert bw.mean == p.reference[sex][ids.BW].value


def test_tanaka_reference_values() -> None:
    """Tanaka 1996 Tables 9 and 10, Asian Reference Man of 20 to 49 years."""
    p = population("east_asian_tanaka")
    male, female = p.reference[Sex.MALE], p.reference[Sex.FEMALE]
    assert (male[ids.BW].value, male[ids.HEIGHT].value) == (60.0, 169.97)
    assert (female[ids.BW].value, female[ids.HEIGHT].value) == (51.28, 155.04)
    assert p.covariates[Sex.MALE][ids.BW].sd == 8.84


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
