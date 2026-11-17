"""The timecourse of a model which is simulated in parts."""

from pathlib import Path

import libsbml
import numpy as np
import pytest
import roadrunner
from sbmlutils.factory import (
    Compartment,
    Model,
    Parameter,
    Reaction,
    Species,
    create_model,
)

from pkpdsim.timecourse import SimPart, Simulation, run_simulation


@pytest.fixture(scope="module")
def model_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A substance `S` which flows in at the rate `k_in` and is degraded (`k_out`)."""
    model = Model(
        sid="inflow",
        compartments=[Compartment("V", 1.0)],
        species=[Species("S", "V", initialConcentration=0.0)],
        parameters=[
            Parameter("k_in", 1.0, constant=False),
            Parameter("k_out", 0.5),
        ],
        reactions=[
            Reaction("IN", equation="-> S", formula="k_in"),
            Reaction("OUT", equation="S ->", formula="k_out * S"),
        ],
    )
    path: Path = tmp_path_factory.mktemp("model") / "inflow.xml"
    create_model(model, filepath=path, validate=False)
    return path


@pytest.fixture
def r(model_path: Path) -> roadrunner.RoadRunner:
    """The model in roadrunner, the time and the concentration of `S` selected."""
    r = roadrunner.RoadRunner(str(model_path))
    r.timeCourseSelections = ["time", "[S]", "k_in"]
    return r


SIM = Simulation(
    sid="inflow_off",
    title="inflow, then no inflow",
    parts=[
        SimPart(tend=10.0, changes={}, steps=10),
        SimPart(tend=5.0, changes={"k_in": 0.0}, steps=20),
    ],
)


def test_the_time_continues_over_the_parts(r: roadrunner.RoadRunner) -> None:
    """A part has `steps + 1` points, its time follows the part before."""
    df = run_simulation(r, SIM)
    assert list(df.columns) == ["time", "[S]", "k_in"]
    assert len(df) == 11 + 21
    assert df["time"].iloc[0] == 0.0
    assert df["time"].iloc[10] == pytest.approx(10.0)
    assert df["time"].iloc[11] == pytest.approx(10.0)
    assert df["time"].iloc[-1] == pytest.approx(15.0)


def test_a_part_continues_from_the_state_of_the_part_before(
    r: roadrunner.RoadRunner,
) -> None:
    """The changes of a part apply from its start, the state continues."""
    df = run_simulation(r, SIM)
    first, second = df.iloc[:11], df.iloc[11:]
    assert (first["k_in"] == 1.0).all()
    assert (second["k_in"] == 0.0).all()
    s_end = 2.0 * (1.0 - np.exp(-0.5 * 10.0))  # k_in / k_out * (1 - exp(-k_out t))
    assert first["[S]"].iloc[-1] == pytest.approx(s_end, rel=1e-5)
    assert second["[S]"].iloc[0] == pytest.approx(s_end, rel=1e-5)
    assert second["[S]"].iloc[-1] == pytest.approx(s_end * np.exp(-0.5 * 5.0), rel=1e-5)


def test_a_simulation_starts_from_the_initial_state(r: roadrunner.RoadRunner) -> None:
    """The model is reset: the changes of a run do not carry over to the next."""
    first = run_simulation(r, SIM)
    second = run_simulation(r, SIM)
    assert second["[S]"].iloc[0] == 0.0
    assert second["k_in"].iloc[0] == 1.0
    np.testing.assert_allclose(first.to_numpy(), second.to_numpy())


def _clock(*, delayed: bool = False) -> roadrunner.RoadRunner:
    document = libsbml.SBMLDocument(3, 2)
    model = document.createModel()
    model.setId("clock")
    for sid in ("x", "flag"):
        parameter = model.createParameter()
        parameter.setId(sid)
        parameter.setValue(0)
        parameter.setConstant(False)
    rule = model.createRateRule()
    rule.setVariable("x")
    rule.setMath(libsbml.parseL3Formula("time"))
    event = model.createEvent()
    event.setUseValuesFromTriggerTime(True)
    trigger = event.createTrigger()
    trigger.setInitialValue(False)
    trigger.setPersistent(True)
    trigger.setMath(libsbml.parseL3Formula("time >= 0.5" if delayed else "time >= 1.5"))
    if delayed:
        event.createDelay().setMath(libsbml.parseL3Formula("1"))
    assignment = event.createEventAssignment()
    assignment.setVariable("flag")
    assignment.setMath(libsbml.parseL3Formula("1"))
    runner = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
    runner.timeCourseSelections = ["time", "x", "flag"]
    runner.integrator.absolute_tolerance = 1e-10
    runner.integrator.relative_tolerance = 1e-10
    return runner


@pytest.mark.parametrize("delayed", [False, True], ids=["event", "delayed-event"])
def test_model_time_and_events_continue_over_parts(delayed: bool) -> None:
    runner = _clock(delayed=delayed)
    simulation = Simulation(
        "clock", "clock", [SimPart(1.0, {}, 10), SimPart(1.0, {}, 10)]
    )
    result = run_simulation(runner, simulation)
    reference = _clock(delayed=delayed).simulate(0, 2, steps=20)
    assert runner.getValue("time") == pytest.approx(2.0)
    np.testing.assert_allclose(result.iloc[-1].to_numpy(), reference[-1], atol=1e-7)
    assert result["x"].iloc[-1] == pytest.approx(2.0, abs=1e-7)
    assert result["flag"].iloc[-1] == 1.0
