"""Tests of the standard human physiology and its scaling with the covariates."""

import copy
from pathlib import Path
from typing import Any

import pytest
import roadrunner
from sbmlutils.factory import Model, Parameter, create_model

from pkpdsim.pbpk import MODEL_UNITS, U, cite, ids, load, physiology, sources
from pkpdsim.pbpk.organs import PORTAL_ORGANS, Organ, Pool
from pkpdsim.pbpk.physiology import Sex, physiology_objects

LEAN = [o for o in Organ if o is not Organ.ADIPOSE]


def _data(sex: Sex) -> dict[str, Any]:
    return load("physiology_human")[sex.value]


def _v(table: dict[str, Any], key: str) -> float:
    return float(table[key]["value"])


def _density(table: dict[str, Any]) -> float:
    return _v(table, "density") if "density" in table else 1.0


@pytest.fixture(scope="module")
def paths(tmp_path_factory: pytest.TempPathFactory) -> dict[Sex, Path]:
    """The physiology of each sex as SBML model."""
    paths = {}
    for sex in Sex:
        path = tmp_path_factory.mktemp("physiology") / f"physiology_{sex}.xml"
        model = Model(
            sid=f"physiology_{sex}",
            units=U,
            model_units=MODEL_UNITS,
            objects=physiology_objects(sex),
        )
        create_model(model, filepath=path)
        paths[sex] = path
    return paths


def _rr(paths: dict[Sex, Path], sex: Sex) -> roadrunner.RoadRunner:
    return roadrunner.RoadRunner(str(paths[sex]))


def _janmahasatian(bw: float, height: float, sex: Sex) -> float:
    """Fat-free mass of Janmahasatian 2005 (kg)."""
    bmi = bw / (height / 100) ** 2
    if sex is Sex.MALE:
        return 9270 * bw / (6680 + 216 * bmi)
    return 9270 * bw / (8780 + 244 * bmi)


def _assert_balance(r: roadrunner.RoadRunner, sex: Sex) -> None:
    """Tissue masses and blood are the body weight, organ blood and pools the blood.

    The tissue mass of an organ is its volume without its blood times its density
    (the organ volume includes its blood, R11, R12).
    """
    p = _data(sex)
    vblood = r[ids.VBLOOD]
    organ_blood = {o: _v(p["organs"][o], "blood") * vblood for o in Organ}
    tissue = sum(
        (r[ids.organ_volume(o)] - organ_blood[o]) * _density(p["organs"][o])
        for o in Organ
    )
    blood_mass = vblood * _v(p["blood"], "density")
    assert tissue + blood_mass == pytest.approx(r[ids.BW], rel=1e-9)
    pools = sum(r[ids.pool_volume(pool)] for pool in Pool)
    assert sum(organ_blood.values()) + pools == pytest.approx(vblood, rel=1e-9)


def _assert_flows(r: roadrunner.RoadRunner) -> None:
    """The organs in parallel receive the cardiac output, the lung in series."""
    co = r[ids.CO]
    total = sum(r[ids.organ_flow(o)] for o in Organ if o is not Organ.LUNG)
    assert total == pytest.approx(co, rel=1e-9)
    assert r[ids.organ_flow(Organ.LUNG)] == pytest.approx(co, rel=1e-12)
    portal = sum(r[ids.organ_flow(o)] for o in PORTAL_ORGANS)
    assert r[ids.QPO] == pytest.approx(portal, rel=1e-12)
    assert r[ids.QHA] == pytest.approx(r[ids.organ_flow(Organ.LIVER)], rel=1e-12)
    assert r[ids.QH] == pytest.approx(r[ids.QHA] + r[ids.QPO], rel=1e-12)


@pytest.mark.parametrize("sex", list(Sex))
def test_reference_individual_reproduces_reference(
    paths: dict[Sex, Path], sex: Sex
) -> None:
    """The reference individual has the reference volumes, flows and functions.

    Lean organs have the volume of their ICRP 89 tissue and blood, adipose the
    remainder of the body weight (R2).
    """
    r = _rr(paths, sex)
    p = _data(sex)
    organs = p["organs"]
    vblood = _v(p["blood"], "volume")
    assert r[ids.SEX] == (0.0 if sex is Sex.MALE else 1.0)
    assert r[ids.BW] == _v(p, "bodyweight")
    assert r[ids.HEIGHT] == _v(p, "height")
    assert r[ids.VBLOOD] == pytest.approx(vblood, rel=1e-9)
    for o in LEAN:
        expected = (
            _v(organs[o], "mass") / _density(organs[o])
            + _v(organs[o], "blood") * vblood
        )
        assert r[ids.organ_volume(o)] == pytest.approx(expected, rel=1e-9), o
    lean = sum(_v(organs[o], "mass") for o in LEAN)
    adipose = _v(p, "bodyweight") - lean - vblood * _v(p["blood"], "density")
    expected = (
        adipose / _density(organs[Organ.ADIPOSE])
        + _v(organs[Organ.ADIPOSE], "blood") * vblood
    )
    assert r[ids.organ_volume(Organ.ADIPOSE)] == pytest.approx(expected, rel=1e-9)
    assert r[ids.organ_mass(Organ.ADIPOSE)] == pytest.approx(adipose, rel=1e-9)
    for pool, key in [(Pool.ARTERIAL, "arterial"), (Pool.PORTAL, "portal")]:
        assert r[ids.pool_volume(pool)] == pytest.approx(
            _v(p["blood"], key) * vblood, rel=1e-9
        )
    assert r[ids.CO] == pytest.approx(_v(p, "cardiac_output"), rel=1e-9)
    assert r[ids.HCT] == _v(p, "hematocrit")
    assert r[ids.GFR] == pytest.approx(_v(p, "gfr") / 1000, rel=1e-9)
    assert r[ids.FFM] == pytest.approx(r[ids.ref(ids.FFM)], rel=1e-12)


@pytest.mark.parametrize("sex", list(Sex))
def test_volumes_sum_to_bodyweight(paths: dict[Sex, Path], sex: Sex) -> None:
    """Tissue masses and blood mass are the body weight (R12)."""
    _assert_balance(_rr(paths, sex), sex)


@pytest.mark.parametrize("sex", list(Sex))
def test_flows_sum_to_cardiac_output(paths: dict[Sex, Path], sex: Sex) -> None:
    """The flows are the ICRP 89 fractions renormalised to cardiac output (R9)."""
    r = _rr(paths, sex)
    _assert_flows(r)
    organs = _data(sex)["organs"]
    for o in LEAN:
        if o is Organ.LUNG:
            continue
        ratio = r[ids.organ_flow(o)] / r[ids.organ_flow(Organ.MUSCLE)]
        expected = _v(organs[o], "flow") / _v(organs[Organ.MUSCLE], "flow")
        assert ratio == pytest.approx(expected, rel=1e-9), o


def test_scaled_individual_sums(paths: dict[Sex, Path]) -> None:
    """An individual of 100 kg and 180 cm, scaled as in Huisinga 2012.

    Brain is constant, skin follows the body surface area, the other lean
    organs, the blood and the cardiac output follow the fat-free mass less
    brain and skin (eqs. 9, 10, 12, 13); adipose grows.
    """
    reference = _rr(paths, Sex.MALE)
    r = _rr(paths, Sex.MALE)
    r.setValue(ids.BW, 100.0)
    r.setValue(ids.HEIGHT, 180.0)

    ffm = _janmahasatian(100.0, 180.0, Sex.MALE)
    bsa = 0.007184 * 100**0.425 * 180**0.725
    assert r[ids.BMI] == pytest.approx(100.0 / 1.8**2, rel=1e-12)
    assert r[ids.BSA] == pytest.approx(bsa, rel=1e-12)
    assert r[ids.FFM] == pytest.approx(ffm, rel=1e-12)

    br, sk = ids.organ_mass(Organ.BRAIN), ids.organ_mass(Organ.SKIN)
    bsa_ratio = bsa / reference[ids.BSA]
    sv = (ffm - reference[br] - reference[sk] * bsa_ratio) / (
        reference[ids.FFM] - reference[br] - reference[sk]
    )
    assert r[ids.SV_TIS] == pytest.approx(sv, rel=1e-12)
    assert r[br] == pytest.approx(reference[br], rel=1e-12)
    assert r[sk] == pytest.approx(reference[sk] * bsa_ratio, rel=1e-12)
    for o in LEAN:
        if o in (Organ.BRAIN, Organ.SKIN):
            continue
        mid = ids.organ_mass(o)
        assert r[mid] == pytest.approx(reference[mid] * sv, rel=1e-12), o
    assert r[ids.VBLOOD] == pytest.approx(reference[ids.VBLOOD] * sv, rel=1e-12)
    vad = ids.organ_volume(Organ.ADIPOSE)
    assert r[vad] > reference[vad] + 10.0
    assert r[ids.CO] == pytest.approx(reference[ids.CO] * sv, rel=1e-12)
    assert r[ids.GFR] == pytest.approx(reference[ids.GFR] * bsa_ratio, rel=1e-12)
    _assert_balance(r, Sex.MALE)
    _assert_flows(r)


def test_extreme_individual_adipose_positive(paths: dict[Sex, Path]) -> None:
    """A lean, tall man of 45 kg and 190 cm keeps a positive adipose volume."""
    r = _rr(paths, Sex.MALE)
    r.setValue(ids.BW, 45.0)
    r.setValue(ids.HEIGHT, 190.0)
    assert r[ids.organ_mass(Organ.ADIPOSE)] > 0
    assert r[ids.organ_volume(Organ.ADIPOSE)] > 0
    _assert_balance(r, Sex.MALE)


def test_sex_switch(paths: dict[Sex, Path]) -> None:
    """SEX = 1 on the male model gives the values of the female model."""
    female = _rr(paths, Sex.FEMALE)
    r = _rr(paths, Sex.MALE)
    p = _data(Sex.FEMALE)
    r.setValue(ids.SEX, 1.0)
    r.setValue(ids.BW, _v(p, "bodyweight"))
    r.setValue(ids.HEIGHT, _v(p, "height"))
    pids = female.model.getGlobalParameterIds()
    assert set(pids) == set(r.model.getGlobalParameterIds())
    for pid in pids:
        assert r[pid] == pytest.approx(female[pid], rel=1e-12), pid
    _assert_balance(r, Sex.FEMALE)


def test_renal_and_cardiac_factors(paths: dict[Sex, Path]) -> None:
    """The renal factor scales GFR, the cardiac factor CO and every flow."""
    reference = _rr(paths, Sex.MALE)
    r = _rr(paths, Sex.MALE)
    r.setValue(ids.F_RENAL_FUNCTION, 0.5)
    r.setValue(ids.F_CARDIAC_FUNCTION, 0.5)
    assert r[ids.GFR] == pytest.approx(0.5 * reference[ids.GFR], rel=1e-12)
    assert r[ids.CO] == pytest.approx(0.5 * reference[ids.CO], rel=1e-12)
    flows = [ids.organ_flow(o) for o in Organ] + [ids.QHA, ids.QPO, ids.QH]
    for qid in flows:
        assert r[qid] == pytest.approx(0.5 * reference[qid], rel=1e-12), qid
    _assert_flows(r)


def test_reset_keeps_covariates(paths: dict[Sex, Path]) -> None:
    """A covariate changed before the simulation holds through it."""
    r = _rr(paths, Sex.MALE)
    r.setValue(ids.BW, 100.0)
    vad = r[ids.organ_volume(Organ.ADIPOSE)]
    r.simulate(0, 10, 3)
    r.reset()
    assert r[ids.BW] == 100.0
    assert r[ids.organ_volume(Organ.ADIPOSE)] == pytest.approx(vad, rel=1e-12)


@pytest.mark.parametrize("sex", list(Sex))
def test_reference_fractions(paths: dict[Sex, Path], sex: Sex) -> None:
    """The vascular and interstitial fractions of the reference are those of the data.

    Adipose is the remainder of the body weight (R2), larger than its ICRP mass
    of the data, so its vascular fraction is smaller than the `f_vas` of the data.
    """
    r = _rr(paths, sex)
    fractions = {o: physiology.reference_fractions(o)[sex] for o in Organ}
    for o, (f_vas, f_int) in fractions.items():
        assert f_vas == _v(_data(sex)["organs"][o], "f_vas")
        if o is Organ.ADIPOSE:
            continue
        assert r[ids.vascular_fraction(o)] == pytest.approx(f_vas, rel=2e-3), o
        assert r[ids.interstitial_fraction(o)] == pytest.approx(f_int, rel=2e-3), o
    vad = ids.organ_volume(Organ.ADIPOSE)
    blood = _v(_data(sex)["organs"][Organ.ADIPOSE], "blood")
    assert r[ids.vascular_fraction(Organ.ADIPOSE)] == pytest.approx(
        blood * r[ids.VBLOOD] / r[vad], rel=1e-12
    )


def test_vascular_spaces_are_blood_volume(paths: dict[Sex, Path]) -> None:
    """For any individual the vascular spaces of the organs and the pools are the blood."""
    r = _rr(paths, Sex.MALE)
    r.setValue(ids.BW, 100.0)
    r.setValue(ids.HEIGHT, 180.0)
    vascular = sum(r[ids.vascular_fraction(o)] * r[ids.organ_volume(o)] for o in Organ)
    pools = sum(r[ids.pool_volume(pool)] for pool in Pool)
    assert vascular + pools == pytest.approx(r[ids.VBLOOD], rel=1e-9)
    for o in Organ:
        fint = _v(_data(Sex.MALE)["organs"][o], "f_int") / (
            1 - _v(_data(Sex.MALE)["organs"][o], "f_vas")
        )
        assert r[ids.interstitial_fraction(o)] == pytest.approx(
            fint * (1 - r[ids.vascular_fraction(o)]), rel=1e-9
        )


@pytest.mark.parametrize("sex", list(Sex))
def test_parameters_annotated(sex: Sex) -> None:
    """Every parameter has an SBO term, a unit, a name and notes with its source.

    The source is a citation, or the value is a model default.
    """
    citations = [cite(key) for key in sources()]
    parameters = [o for o in physiology_objects(sex) if isinstance(o, Parameter)]
    defaults = {p.sid for p in parameters if physiology.MODEL_DEFAULT in str(p.notes)}
    assert defaults == {ids.AGE, ids.F_RENAL_FUNCTION, ids.F_CARDIAC_FUNCTION}
    for p in parameters:
        assert p.sboTerm, p.sid
        assert p.unit is not None, p.sid
        assert p.name, p.sid
        if p.sid not in defaults:
            assert any(c in str(p.notes) for c in citations), p.sid


def test_negative_reference_adipose_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reference data which leave no room for adipose are rejected."""
    data = copy.deepcopy(load("physiology_human"))
    data["female"]["organs"]["mu"]["mass"]["value"] = 60.0

    def fake_load(name: str) -> dict[str, Any]:
        assert name == "physiology_human"
        return data

    monkeypatch.setattr(physiology, "load", fake_load)
    with pytest.raises(ValueError, match=r"female.*adipose.*kg"):
        physiology_objects(Sex.MALE)


def test_unit_of_data_checked(monkeypatch: pytest.MonkeyPatch) -> None:
    """A unit of the data which the model cannot convert is rejected."""
    data = copy.deepcopy(load("physiology_human"))
    data["male"]["gfr"]["unit"] = "ml/h"
    monkeypatch.setattr(physiology, "load", lambda name: data)
    with pytest.raises(ValueError, match=r"'gfr'.*ml/h.*l/min"):
        physiology_objects(Sex.MALE)
