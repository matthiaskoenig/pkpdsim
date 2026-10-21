"""The timecourse of a model which is simulated in parts."""

from pathlib import Path

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
