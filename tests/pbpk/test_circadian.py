"""The circadian rhythms: the clock time and the factors over the day and the night."""

from pathlib import Path

import libsbml
import numpy as np
import pytest

from pkpdsim.pbpk import Circadian, Organ, PBPKSimulationExperiment, build, ids
from pkpdsim.pbpk.circadian import (
    circadian_factor,
    protein_rhythms,
    rhythms,
)

from .conftest import definition, individuals, total_amount


@pytest.fixture(scope="module")
def circadian(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The simple definition (CYP3A4 of the liver) with the default rhythms."""
    return build(
        definition(sid="simple_circadian", circadian=Circadian()),
        tmp_path_factory.mktemp("circadian"),
    ).flat_path


def _parameters(path: Path) -> set[str]:
    """The parameters of a built model."""
    model = libsbml.readSBMLFromFile(str(path)).getModel()
    return {p.getId() for p in model.getListOfParameters()}


def test_the_rhythms_of_the_table() -> None:
    """GFR, fed gastric emptying and CYP3A4 by default, albumin by choice."""
    assert {k for k, r in rhythms().items() if r.default} == {"GFR", "t_stomach_fed"}
    assert not rhythms()["albumin"].default
    cyp3a4 = next(r for (o, p), r in protein_rhythms().items() if p == "CYP3A4")
    assert cyp3a4.key == "li_cyp3a4"
    for r in [*rhythms().values(), *protein_rhythms().values()]:
        assert r.sources
        mean = np.mean([r.value(h) for h in np.linspace(0, 24, 241)[:-1]])
        assert mean == pytest.approx(1.0, abs=1e-12)
        assert r.value(r.acrophase) == pytest.approx(1.0 + r.amplitude, rel=1e-12)


def test_a_circadian_is_checked() -> None:
    """An unknown rhythm and a clock time outside of the day are refused."""
    with pytest.raises(ValueError, match="unknown rhythms"):
        Circadian(rhythms=("liver_flow",))
    with pytest.raises(ValueError, match="clock time"):
        Circadian(clock=24.0)
    with pytest.raises(ValueError, match="clock time"):
        PBPKSimulationExperiment.clock_changes(-1.0)
    change = PBPKSimulationExperiment.clock_changes(20.0)[ids.CLOCK]
    assert change.to("hr").magnitude == 20.0


def test_the_model_computes_the_factors_over_the_day(circadian: Path) -> None:
    """Every factor is its cosinor of the clock time, `TOD` the clock time modulo 24 h."""
    parameters = _parameters(circadian)
    keys = ["GFR", "t_stomach_fed", "li_cyp3a4"]
    assert {circadian_factor(k) for k in keys} <= parameters
    assert circadian_factor("albumin") not in parameters
    curves = {**rhythms(), **{r.key: r for r in protein_rhythms().values()}}
    for start in (8.0, 20.0):
        r = individuals(circadian)(CLOCK=start)
        selections = ["time", ids.TIME_OF_DAY, *(circadian_factor(k) for k in keys)]
        s = np.asarray(r.simulate(0, 30 * 60, 31, selections))
        clock = start + s[:, 0] / 60
        assert s[:, 1] == pytest.approx(np.mod(clock, 24.0), abs=1e-9)
        for i, key in enumerate(keys, start=2):
            expected = [curves[key].value(c) for c in clock]
            assert s[:, i] == pytest.approx(expected, rel=1e-12)


def test_the_quantities_follow_their_factors(circadian: Path) -> None:
    """The GFR and the activity of CYP3A4 change with their factors during a simulation."""
    r = individuals(circadian)(CLOCK=8.0)
    gfr0 = r[ids.GFR] / r[circadian_factor("GFR")]
    activity = ids.impairment_activity("CYP3A4", Organ.LIVER)
    a0 = r[activity] / r[circadian_factor("li_cyp3a4")]
    s = np.asarray(
        r.simulate(
            0,
            24 * 60,
            25,
            [ids.GFR, circadian_factor("GFR"), activity, circadian_factor("li_cyp3a4")],
        )
    )
    assert s[:, 0] == pytest.approx(gfr0 * s[:, 1], rel=1e-9)
    assert s[:, 2] == pytest.approx(a0 * s[:, 3], rel=1e-9)
    assert s[:, 1].max() > 1.14
    assert s[:, 1].min() < 0.86


def test_a_meal_in_the_evening_leaves_the_stomach_later(circadian: Path) -> None:
    """The fed transit of the stomach follows the clock time; the fasted one does not."""
    at = individuals(circadian)
    rhythm = rhythms()["t_stomach_fed"]
    morning = at(CLOCK=8.0, FED=1.0)[ids.T_STOMACH]
    evening = at(CLOCK=20.0, FED=1.0)[ids.T_STOMACH]
    assert evening / morning == pytest.approx(
        rhythm.value(20.0) / rhythm.value(8.0), rel=1e-9
    )
    fasted = [at(CLOCK=c, FED=0.0)[ids.T_STOMACH] for c in (8.0, 20.0)]
    assert fasted[0] == pytest.approx(fasted[1], rel=1e-12)


def test_a_dose_is_conserved_over_the_day(circadian: Path) -> None:
    """The drug of an oral dose stays in the body and its sinks under the rhythms."""
    r = individuals(circadian)(CLOCK=20.0, FED=1.0, **{ids.po_dose("drg"): 10.0})
    dose = 10.0 / r[ids.mr("drg")]
    r.simulate(0, 24 * 60, 50)
    total = (
        total_amount(r, "drg")
        + total_amount(r, "drg_m")
        + r[ids.po_dose("drg")] / r[ids.mr("drg")]
    )
    assert total == pytest.approx(dose, rel=1e-6)


def test_a_chosen_rhythm_only(tmp_path: Path) -> None:
    """`Circadian(rhythms=...)` builds the chosen factors; without it the model has none."""
    albumin = build(
        definition(sid="simple_albumin", circadian=Circadian(rhythms=("albumin",))),
        tmp_path / "albumin",
    ).flat_path
    parameters = _parameters(albumin)
    assert circadian_factor("albumin") in parameters
    assert circadian_factor("GFR") not in parameters
    plain = build(definition(sid="simple_plain"), tmp_path / "plain").flat_path
    parameters = _parameters(plain)
    assert not {ids.CLOCK, ids.TIME_OF_DAY} & parameters
    assert not [p for p in parameters if p.startswith("fcirc_")]
