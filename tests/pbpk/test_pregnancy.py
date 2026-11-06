"""The pregnancy: the maternal physiology over the gestational age and the fetus."""

from pathlib import Path

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import Organ, Pregnancy, Sex, Species, build, ids
from pkpdsim.pbpk.pregnancy import (
    ORGAN_FLOW,
    PLACENTA_SHARE,
    TERM,
    factors,
    pregnancy_factor,
    pregnancy_objects,
    protein_factors,
)

from .conftest import definition, individuals, total_amount

TRIMESTERS: tuple[float, ...] = (0.0, 13.0, 26.0, TERM)
"""The gestational ages (weeks) of the comparisons."""


@pytest.fixture(scope="module")
def pregnant(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The simple definition as a pregnant female build with the fetoplacental unit."""
    return build(
        definition(sid="simple_pregnant", sex=Sex.FEMALE, pregnancy=Pregnancy()),
        tmp_path_factory.mktemp("pregnancy"),
    ).flat_path


def _parameters(path: Path) -> set[str]:
    """The parameters of a built model."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    return {p.getId() for p in model.getListOfParameters()}


def test_every_factor_is_one_before_the_pregnancy() -> None:
    """Every factor is 1 (a mass 0) at GA 0 and up to the fertilization at GA 2."""
    for f in (*factors().values(), *protein_factors().values()):
        assert f.sources
        expected = 0.0 if f.unit == "kg" or f.key == "upla" else 1.0
        for ga in (0.0, 1.0, 2.0):
            assert f.value(ga) == pytest.approx(expected, abs=1e-12), (f.key, ga)
    assert factors()["upla"].value(0.0) == 0.0


def test_the_factors_at_term() -> None:
    """The published changes at term: the cardiac output, the GFR, the albumin, the fetus."""
    at_term = {k: f.value(TERM) for k, f in factors().items()}
    assert at_term[ids.CO] == pytest.approx(1.40, abs=0.02)
    assert at_term[ids.VBLOOD] == pytest.approx(1.45, abs=0.02)
    assert at_term[ids.HCT] == pytest.approx(0.93, abs=0.02)
    assert at_term[ids.GFR] == pytest.approx(1.33, abs=0.02)
    assert at_term["albumin"] == pytest.approx(0.74, abs=0.02)
    assert at_term["upla"] == pytest.approx(0.12, abs=0.02)
    assert at_term["mass_fetus"] == pytest.approx(3.5, abs=0.1)
    cyp3a4 = protein_factors()[(Organ.LIVER, "CYP3A4")]
    assert cyp3a4.value(TERM) == pytest.approx(1.6, abs=0.01)
    cyp1a2 = protein_factors()[(Organ.LIVER, "CYP1A2")]
    assert cyp1a2.value(TERM) < 0.4


def test_a_pregnancy_is_checked() -> None:
    """The gestational age, the male build, an animal and kinetic erythrocytes are refused."""
    with pytest.raises(ValueError, match="gestational age"):
        Pregnancy(ga=41.0)
    with pytest.raises(ValueError, match="female build"):
        build(definition(sid="male", pregnancy=Pregnancy()), Path("unused"))
    with pytest.raises(ValueError, match="the pregnancy of the library is the human"):
        pregnancy_objects(
            definition(sid="rat", species=Species.RAT, pregnancy=Pregnancy())
        )


def test_the_model_computes_the_factors(pregnant: Path) -> None:
    """Every factor of the model is its value of python at every gestational age."""
    at = individuals(pregnant)
    curves = {pregnancy_factor(k): f for k, f in factors().items()}
    cyp3a4 = protein_factors()[(Organ.LIVER, "CYP3A4")]
    curves[pregnancy_factor(cyp3a4.key)] = cyp3a4
    for ga in (0.0, 5.0, 13.0, 20.0, 26.0, 33.0, TERM):
        r = at(GA=ga)
        for pid, f in curves.items():
            assert r[pid] == pytest.approx(f.value(ga), rel=1e-9), (pid, ga)


def test_the_physiology_follows_the_gestational_age(pregnant: Path) -> None:
    """The cardiac output, the blood, the GFR and fu of the woman change with GA."""
    at = individuals(pregnant)
    before = at(GA=0.0)
    keys = (ids.CO, ids.VBLOOD, ids.HCT, ids.GFR, "fu_drg", ids.BW, ids.BW_TOTAL)
    values = {k: before[k] for k in keys}
    term = at(GA=TERM)
    # the cardiac output includes the placenta, the organs share the rest
    assert term[ids.CO] / values[ids.CO] == pytest.approx(1.37, abs=0.03)
    assert term[ids.VBLOOD] / values[ids.VBLOOD] == pytest.approx(1.45, abs=0.02)
    assert term[ids.HCT] / values[ids.HCT] == pytest.approx(0.93, abs=0.02)
    assert term[ids.GFR] / values[ids.GFR] == pytest.approx(1.33, abs=0.02)
    # the lower albumin raises the unbound fraction
    assert term["fu_drg"] > values["fu_drg"]
    # the body weight of the covariates is the one before the pregnancy
    assert term[ids.BW] == values[ids.BW]
    assert term[ids.BW_TOTAL] - values[ids.BW_TOTAL] == pytest.approx(6.1, abs=0.3)


def test_the_cardiac_output_is_the_sum_of_the_flows(pregnant: Path) -> None:
    """The organs and the placenta receive the cardiac output at every gestational age."""
    at = individuals(pregnant)
    for ga in TRIMESTERS:
        r = at(GA=ga)
        organs = sum(r[ids.organ_flow(o)] for o in Organ if o is not Organ.LUNG)
        assert organs + r[ids.Q_PLACENTA] == pytest.approx(r[ids.CO], rel=1e-9), ga
        assert r[ids.organ_flow(Organ.LUNG)] == pytest.approx(r[ids.CO], rel=1e-12)
    assert at(GA=0.0)[ids.Q_PLACENTA] == 0.0
    assert at(GA=TERM)[ids.Q_PLACENTA] == pytest.approx(0.90, abs=0.05)


def test_the_fetus_tracks_the_maternal_plasma(pregnant: Path) -> None:
    """With the default transfer the fetal concentration follows the maternal one."""
    at = individuals(pregnant)
    selections = ["time", "[Car_drg]", "[Cpla_drg]", "[Cfet_drg]"]
    ratios = []
    for ga in (13.0, 26.0, TERM):
        r = at(GA=ga, **{ids.po_dose("drg"): 10.0})
        s = np.asarray(r.simulate(0, 48 * 60, 1441, selections))
        auc = [float(np.trapezoid(s[:, i], s[:, 0])) for i in range(1, 4)]
        assert auc[1] / auc[0] == pytest.approx(1.0, abs=0.02), ga
        assert auc[2] / auc[0] == pytest.approx(1.0, abs=0.05), ga
        ratios.append(auc[0])
    # the exposure of the mother falls with the gestational age (dilution, induction)
    assert ratios[0] > ratios[1] > ratios[2]


def test_a_dose_is_conserved_with_the_fetus(pregnant: Path) -> None:
    """The drug of an oral dose stays in the mother, the placenta, the fetus and the sinks."""
    at = individuals(pregnant)
    for ga in (0.0, TERM):
        r = at(GA=ga, **{ids.po_dose("drg"): 10.0})
        dose = 10.0 / r[ids.mr("drg")]
        r.simulate(0, 48 * 60, 100)
        total = (
            total_amount(r, "drg")
            + total_amount(r, "drg_m")
            + r[ids.po_dose("drg")] / r[ids.mr("drg")]
        )
        assert total == pytest.approx(dose, rel=1e-6), ga


def test_the_fetus_holds_the_drug_of_its_partition(pregnant: Path) -> None:
    """A fetus with `Kp_fet` 3 holds three times the concentration of the placenta."""
    r = individuals(pregnant)(
        GA=TERM, **{ids.kp_fetus("drg"): 3.0, ids.po_dose("drg"): 10.0}
    )
    s = np.asarray(r.simulate(0, 72 * 60, 721, ["time", "[Cpla_drg]", "[Cfet_drg]"]))
    assert s[-1, 2] / s[-1, 1] == pytest.approx(3.0, rel=0.05)


def test_without_a_pregnancy_the_model_has_none_of_its_elements(tmp_path: Path) -> None:
    """A model without `pregnancy`, and one without the fetus, have no fetal elements."""
    plain = build(definition(sid="plain", sex=Sex.FEMALE), tmp_path / "plain").flat_path
    parameters = _parameters(plain)
    assert ids.GA not in parameters
    assert not [p for p in parameters if p.startswith("fpreg_")]
    mother = build(
        definition(
            sid="mother", sex=Sex.FEMALE, pregnancy=Pregnancy(ga=30.0, fetus=False)
        ),
        tmp_path / "mother",
    ).flat_path
    parameters = _parameters(mother)
    assert pregnancy_factor(ids.CO) in parameters
    assert pregnancy_factor(ORGAN_FLOW) in parameters
    assert pregnancy_factor(PLACENTA_SHARE) in parameters
    assert ids.Q_PLACENTA not in parameters
    assert individuals(mother)()[ids.GA] == 30.0
