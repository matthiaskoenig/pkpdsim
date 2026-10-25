"""Tests of the reference physiology data (ICRP 89 male and female)."""

from typing import Any

import pytest

from pkpdsim.pbpk import Organ, load
from pkpdsim.pbpk.organs import PORTAL_ORGANS

SEXES = ["male", "female"]


def _physiology(sex: str) -> dict[str, Any]:
    return load("physiology_human")[sex]


def _v(table: dict[str, Any], key: str) -> float:
    return float(table[key]["value"])


@pytest.mark.parametrize(
    ("sex", "bodyweight", "height"), [("male", 73.0, 176.0), ("female", 60.0, 163.0)]
)
def test_reference_bodyweight(sex: str, bodyweight: float, height: float) -> None:
    """The reference individuals of ICRP 89 Table 2.9."""
    p = _physiology(sex)
    assert _v(p, "bodyweight") == bodyweight
    assert p["bodyweight"]["unit"] == "kg"
    assert _v(p, "height") == height
    assert p["height"]["unit"] == "cm"


@pytest.mark.parametrize("sex", SEXES)
def test_every_organ_present(sex: str) -> None:
    """Every organ has a flow, blood content, f_vas and f_int; all but adipose a mass.

    Adipose is the remainder of the body weight, its ICRP mass is a reference.
    """
    organs = _physiology(sex)["organs"]
    assert set(organs) == {o.value for o in Organ}
    for organ in Organ:
        for key in ("flow", "blood", "f_vas", "f_int"):
            assert key in organs[organ], f"{organ}: {key}"
        if organ is not Organ.ADIPOSE:
            assert organs[organ]["mass"]["unit"] == "kg"


@pytest.mark.parametrize("sex", SEXES)
def test_flow_fractions_sum_to_one(sex: str) -> None:
    """The flows of the organs in parallel and the unmapped tissues sum to one.

    The lung is in series and receives the cardiac output; the liver states
    its arterial flow, the portal organs feed it in addition.
    """
    p = _physiology(sex)
    organs = p["organs"]
    assert _v(organs[Organ.LUNG], "flow") == 1.0
    total = sum(_v(organs[o], "flow") for o in Organ if o is not Organ.LUNG)
    total += sum(_v(t, "flow") for t in p["unmapped"].values() if "flow" in t)
    assert total == pytest.approx(1.0, abs=1e-3)


@pytest.mark.parametrize("sex", SEXES)
def test_blood_fractions_sum_to_one(sex: str) -> None:
    """The blood of the organs, the pools and the unmapped tissues is all blood."""
    p = _physiology(sex)
    total = sum(_v(organ, "blood") for organ in p["organs"].values())
    total += sum(_v(p["blood"], pool) for pool in ("arterial", "venous", "portal"))
    total += sum(_v(t, "blood") for t in p["unmapped"].values() if "blood" in t)
    assert total == pytest.approx(1.0, abs=1e-3)


@pytest.mark.parametrize("sex", SEXES)
def test_organ_masses_below_bodyweight(sex: str) -> None:
    """The organs without adipose and the blood leave room for adipose.

    The remainder, adipose in the physiology, holds at least the ICRP adipose.
    """
    p = _physiology(sex)
    blood = _v(p["blood"], "volume") * _v(p["blood"], "density")
    lean = sum(
        _v(table, "mass")
        for organ, table in p["organs"].items()
        if organ != Organ.ADIPOSE
    )
    remainder = _v(p, "bodyweight") - lean - blood
    assert remainder > 0
    assert remainder >= _v(p["organs"][Organ.ADIPOSE], "mass")


@pytest.mark.parametrize("sex", SEXES)
def test_masses_sum_to_bodyweight(sex: str) -> None:
    """Organs, blood and the unmapped tissues are the ICRP 89 reference body."""
    p = _physiology(sex)
    total = sum(_v(table, "mass") for table in p["organs"].values())
    total += sum(_v(t, "mass") for t in p["unmapped"].values() if "mass" in t)
    total += _v(p["blood"], "volume") * _v(p["blood"], "density")
    assert total == pytest.approx(_v(p, "bodyweight"), rel=2e-3)


@pytest.mark.parametrize("sex", SEXES)
def test_fractions_in_range(sex: str) -> None:
    """Vascular and interstitial space are parts of the organ volume."""
    for organ, table in _physiology(sex)["organs"].items():
        f_vas, f_int = _v(table, "f_vas"), _v(table, "f_int")
        assert 0 < f_vas < 1, organ
        assert 0 < f_int < 1, organ
        assert f_vas + f_int < 1, organ


@pytest.mark.parametrize("sex", SEXES)
def test_vascular_fraction_from_blood_content(sex: str) -> None:
    """f_vas is the organ blood of the organ volume including that blood."""
    p = _physiology(sex)
    volume = _v(p["blood"], "volume")
    for organ, table in p["organs"].items():
        density = _v(table, "density") if "density" in table else 1.0
        blood = _v(table, "blood") * volume
        f_vas = blood / (_v(table, "mass") / density + blood)
        assert _v(table, "f_vas") == pytest.approx(f_vas, rel=1e-3), organ


@pytest.mark.parametrize("sex", SEXES)
def test_hepatic_artery_fraction(sex: str) -> None:
    """The hepatic artery fraction is the arterial share of the liver inflow."""
    p = _physiology(sex)
    organs = p["organs"]
    arterial = _v(organs[Organ.LIVER], "flow")
    portal = sum(_v(organs[o], "flow") for o in PORTAL_ORGANS)
    assert _v(p, "hepatic_artery_fraction") == pytest.approx(
        arterial / (arterial + portal), rel=1e-3
    )


@pytest.mark.parametrize("sex", SEXES)
def test_blood_and_function(sex: str) -> None:
    """Hematocrit, cardiac output and GFR are in the ranges of adults."""
    p = _physiology(sex)
    assert 0.35 < _v(p, "hematocrit") < 0.5
    assert 5.0 < _v(p, "cardiac_output") < 7.0
    assert p["cardiac_output"]["unit"] == "l/min"
    assert 90 < _v(p, "gfr") < 140
    assert p["gfr"]["unit"] == "ml/min"
    assert 3.5 < _v(p["blood"], "volume") < 6.0
