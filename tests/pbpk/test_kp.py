"""Tests of the Kp method of Rodgers & Rowland and its build."""

import copy
import dataclasses
from pathlib import Path

import pytest

from pkpdsim.pbpk import (
    Erythrocytes,
    Excretion,
    KpMethod,
    Mode,
    Organ,
    Route,
    Substance,
    U,
    Value,
    build,
    ids,
    load,
)
from pkpdsim.pbpk import kp as kp_module
from pkpdsim.pbpk.kp import (
    Composition,
    composition,
    kp_sources,
    kp_values,
    rodgers_rowland,
)
from pkpdsim.pbpk.physiology import Sex, reference_hematocrit

from .conftest import definition, dose_amount, individual, simple_drug, total_amount
from .conftest import load as load_model

DAY = 24 * 60.0


def _method(**changes: object) -> KpMethod:
    """A neutral compound of logP 2, merged tissue."""
    fields: dict[str, object] = {
        "method": "rodgers_rowland",
        "logp": 2.0,
        "pka": (),
        "compound": "neutral",
        "ftissue": 1.0,
    }
    fields.update(changes)
    return KpMethod(**fields)  # ty: ignore[invalid-argument-type]


def _drug(method: KpMethod, **changes: object) -> Substance:
    """The simple drug (fu 0.1) with a Kp method."""
    return dataclasses.replace(simple_drug(), distribution=method, **changes)


# Coutinho 2024 (Pharm Res 41:1121, open access), Supplementary Tables S5 and S7:
# Kpu of the Rodgers & Rowland equation (their Eq. 3) for 13 tissues, computed
# with the rat composition of their Table III, which is the composition of
# data/tissue_composition.toml (adipose aside, whose two lipid fractions their
# table exchanges). Drug data of their Table II: itraconazole base pKa 4.57,
# fup 0.002; isavuconazole base pKa 3.28, fup 0.01. Only the columns where the
# plasma protein term `1/fu - 1 - lipids_p` is positive: Coutinho 2024 do not
# set a negative affinity to zero (Yau 2020, line 119). Rodgers & Rowland
# 2005/2006 were not accessible, so their own examples are not reproduced.
_PUBLISHED: list[tuple[str, float, float, float, Organ, float]] = [
    # drug, pKa, fup, logP, organ, Kpu
    ("itraconazole", 4.57, 0.002, 4.0, Organ.MUSCLE, 152.4),  # Table S5
    ("itraconazole", 4.57, 0.002, 4.0, Organ.BONE, 226.2),
    ("itraconazole", 4.57, 0.002, 4.0, Organ.SKIN, 746.7),
    ("itraconazole", 4.57, 0.002, 4.5, Organ.BRAIN, 1269.5),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.SMALL_INTESTINE, 4153.1),  # gut
    ("itraconazole", 4.57, 0.002, 5.0, Organ.LIVER, 2081.4),
    ("itraconazole", 4.57, 0.002, 5.0, Organ.KIDNEY, 1957.8),
    ("isavuconazole", 3.28, 0.01, 3.0, Organ.MUSCLE, 19.1),  # Table S7
    ("isavuconazole", 3.28, 0.01, 4.0, Organ.SKIN, 636.8),
    ("isavuconazole", 3.28, 0.01, 4.5, Organ.BRAIN, 1252.0),
]


@pytest.mark.parametrize(("drug", "pka", "fu", "logp", "organ", "kpu"), _PUBLISHED)
def test_rr_published_values(
    drug: str, pka: float, fu: float, logp: float, organ: Organ, kpu: float
) -> None:
    """The weak bases of Coutinho 2024 within the rounding of their tables.

    The tolerance is the rounding of the publication (one decimal) or 0.1 %,
    the three significant digits of the composition of their Table III.
    """
    method = _method(logp=logp, pka=(pka,), compound="base")
    kp = rodgers_rowland(method, organ, fu, None, 0.45)
    assert kp / fu == pytest.approx(kpu, abs=0.05, rel=1e-3), drug
    assert kp == pytest.approx(kpu * fu, abs=0.05 * fu, rel=1e-3)


def test_rr_neutral_without_lipophilicity() -> None:
    """A hydrophilic neutral compound without binding distributes into the water."""
    method = _method(logp=-10.0)
    c = composition(Organ.MUSCLE)
    kp = rodgers_rowland(method, Organ.MUSCLE, 1.0, None, 0.45)
    assert kp == pytest.approx(c.f_ew + c.f_iw + 0.7 * c.f_np, rel=1e-6)


def test_rr_base_reproduces_erythrocytes(monkeypatch: pytest.MonkeyPatch) -> None:
    """The affinity of a base follows from the erythrocytes (eqs. 3, 5 of Yau 2020).

    The equation of the tissues applied to the composition of the erythrocytes
    (at the pH of the erythrocytes) gives back the erythrocyte to unbound
    plasma coefficient `Kpu_bc = (HCT - 1 + BP)/(fu*HCT)`.
    """
    data = copy.deepcopy(kp_module._data())
    data["intracellular"]["ph"] = data["erythrocytes"]["ph"]
    monkeypatch.setattr(kp_module, "_data", lambda: data)
    method = _method(logp=1.5, pka=(9.2,), compound="base")
    fu, bp, hct = 0.4, 1.3, 0.45
    erythrocytes = Composition(
        "erythrocytes",
        f_ew=0.0,
        f_iw=0.63,
        f_nl=0.0033,
        f_np=0.0012,
        albumin=0.0,
        lipoprotein=0.0,
        ap=0.5,
    )
    kp = rodgers_rowland(method, Organ.MUSCLE, fu, bp, hct, tissue=erythrocytes)
    assert kp / fu == pytest.approx((hct - 1 + bp) / (fu * hct), rel=1e-9)


def test_rr_strong_base_binds_phospholipids() -> None:
    """A strong base is concentrated in the tissues rich in acidic phospholipids."""
    method = _method(logp=1.0, pka=(9.5,), compound="base")
    kp = {o: rodgers_rowland(method, o, 0.5, 2.0, 0.45) for o in Organ}
    assert kp[Organ.KIDNEY] > kp[Organ.MUSCLE] > kp[Organ.BRAIN]
    assert kp[Organ.LUNG] > 10.0
    # a low BP gives a negative affinity, which is zero (Yau 2020)
    low = rodgers_rowland(method, Organ.KIDNEY, 0.5, 0.55, 0.45)
    assert low < kp[Organ.KIDNEY]
    assert low > 0


def test_rr_acid_ionisation() -> None:
    """An acid is excluded from the cells, which are more acidic than plasma."""
    neutral = rodgers_rowland(_method(logp=0.0), Organ.MUSCLE, 1.0, None, 0.45)
    acid = rodgers_rowland(
        _method(logp=0.0, pka=(4.0,), compound="acid"), Organ.MUSCLE, 1.0, None, 0.45
    )
    c = composition(Organ.MUSCLE)
    # X/Y = (1 + 10^(7.0 - 4))/(1 + 10^(7.4 - 4)) for the intracellular water
    x_y = (1 + 10**3.0) / (1 + 10**3.4)
    lipids = (c.f_nl + c.f_np) / (1 + 10**3.4)
    assert acid == pytest.approx(c.f_ew + x_y * c.f_iw + lipids, rel=1e-9)
    assert acid < neutral


def test_rr_strong_base_needs_bp() -> None:
    """Without a blood to plasma ratio a strong base cannot be computed."""
    method = _method(logp=1.0, pka=(9.5,), compound="base")
    with pytest.raises(ValueError, match="blood to plasma"):
        rodgers_rowland(method, Organ.LIVER, 0.5, None, 0.45)


def test_surrogate_tissues() -> None:
    """Organs without a composition of their own take a documented surrogate."""
    assert composition(Organ.STOMACH).tissue == "gut"
    assert composition(Organ.LARGE_INTESTINE).tissue == "gut"
    assert composition(Organ.GONADS).tissue == "muscle"
    assert "ModelChoice" in kp_sources(Organ.GONADS)
    assert "ModelChoice" not in kp_sources(Organ.LIVER)
    assert {composition(o).tissue for o in Organ} >= {"adipose", "liver", "skin"}


def test_kp_values_use_sex_hematocrit() -> None:
    """The coefficients of a base depend on the reference hematocrit of the sex."""
    drug = _drug(
        _method(logp=1.0, pka=(9.5,), compound="base"),
        erythrocytes=Erythrocytes(mode="equilibrium", bp=1.5),
    )
    male, female = kp_values(drug, Sex.MALE), kp_values(drug, Sex.FEMALE)
    assert set(male) == set(Organ)
    assert reference_hematocrit(Sex.MALE) != reference_hematocrit(Sex.FEMALE)
    assert male[Organ.LIVER].value != female[Organ.LIVER].value
    hct = reference_hematocrit(Sex.MALE)
    expected = rodgers_rowland(drug.distribution, Organ.LIVER, 0.1, 1.5, hct)  # ty: ignore[invalid-argument-type]
    assert male[Organ.LIVER].value == pytest.approx(expected, rel=1e-12)


def test_kp_values_of_entered_organs() -> None:
    """Only the organs the substance enters get a coefficient."""
    organs = frozenset({Organ.LIVER, Organ.MUSCLE})
    values = kp_values(_drug(_method(organs=organs)), Sex.MALE)
    assert set(values) == organs


def test_kp_method_factor(tmp_path: Path) -> None:
    """Every `Kp_<o>_<s>` is the factor `fKp_<s>` times the value of the method.

    The factor is a parameter of the top model, a fit value carries its
    uncertainty; a change of the factor before the simulation scales every
    coefficient.
    """
    method = _method(factor=Value(1.5, U.dimensionless, "fit"))
    drug = _drug(method)
    result = build(definition(substances=[drug], chemistry=[]), tmp_path)
    values = kp_values(drug, Sex.MALE)
    r = load_model(result)
    assert r[ids.fkp("drg")] == 1.5
    for organ, v in values.items():
        assert r[ids.kp_organ(organ, "drg")] == pytest.approx(1.5 * v.value, rel=1e-12)
    assert ids.kp("drg") not in r.model.getGlobalParameterIds()

    changed = individual(result.flat_path, **{ids.fkp("drg"): 3.0})
    for organ, v in values.items():
        assert changed[ids.kp_organ(organ, "drg")] == pytest.approx(
            3.0 * v.value, rel=1e-12
        )


@pytest.mark.parametrize("mode", list(Mode))
def test_kp_method_modes_conserve_mass(tmp_path: Path, mode: Mode) -> None:
    """Perfusion, merged and full tissue build with a Kp method and keep the dose."""
    drug = _drug(
        _method(logp=1.0, pka=(9.0,), compound="base", mode=mode),
        erythrocytes=Erythrocytes(mode="equilibrium", bp=1.2),
        excretion=Excretion(),
    )
    d = definition(substances=[drug], chemistry=[], routes=[Route.iv("drg")])
    r = load_model(build(d, tmp_path))
    r[ids.iv_dose("drg")] = 10.0
    start = total_amount(r, "drg") + dose_amount(r, "drg")
    r.simulate(0, DAY, 50)
    assert start > 0
    assert total_amount(r, "drg") + dose_amount(r, "drg") == pytest.approx(
        start, rel=1e-6
    )
    tissue = {
        Mode.TISSUE: "MU__Cti_drg",
        Mode.PERMEABILITY: "MU__Cic_drg",
        Mode.PERFUSION: "MU__Cor_drg",
    }[mode]
    assert r[f"[{tissue}]"] > 0


def test_kp_method_validation() -> None:
    """The definition refuses a strong base without BP and wrong pKa counts."""
    base = _drug(_method(logp=1.0, pka=(9.0,), compound="base"))
    with pytest.raises(ValueError, match=r"drg: Rodgers & Rowland needs the blood"):
        definition(substances=[base], chemistry=[]).validate()
    weak = _drug(_method(logp=1.0, pka=(5.0,), compound="base"))
    definition(substances=[weak], chemistry=[]).validate()
    for pka, compound, match in [
        ((4.0,), "neutral", "neutral compound has no pKa"),
        ((), "acid", "at least one pKa"),
        ((3.0,), "zwitterion", "two pKa"),
    ]:
        drug = _drug(_method(pka=pka, compound=compound))
        with pytest.raises(ValueError, match=match):
            definition(substances=[drug], chemistry=[]).validate()


def test_kp_notes_cite_sources(tmp_path: Path) -> None:
    """The coefficients cite the equations and the composition of their tissue."""
    import libsbml

    drug = _drug(_method())
    result = build(definition(substances=[drug], chemistry=[]), tmp_path)
    doc = libsbml.readSBMLFromFile(str(result.flat_path))
    notes = (
        doc.getModel().getParameter(ids.kp_organ(Organ.LIVER, "drg")).getNotesString()
    )
    for text in ("10.1208/s12248-020-0418-7", "10.1007/s11095-024-03703-4", "fKp_drg"):
        assert text in notes
    assert load("tissue_composition")["organs"]["li"]["tissue"] == "liver"
