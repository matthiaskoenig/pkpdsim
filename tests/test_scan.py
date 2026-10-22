"""The results of the scans of a model, on the toy model."""

from collections.abc import Sequence

import numpy as np
import pytest
from sbmlsim import Q
from sbmlsim.model import AbstractModel
from sbmlsim.result import ScanResult
from sbmlsim.simulation import Change, Dimension, Scan, Simulation
from sbmlsim.simulator import Simulator

from pkpdsim.scan import (
    end_values,
    point_timecourse,
    single_dimension,
    start_values,
    state_after_changes,
)
from pkpdsim_toy import MODEL_PATH

TOY = AbstractModel(
    source=MODEL_PATH,
    language_type=AbstractModel.LanguageType.SBML,
    selections=["time", "[A]"],
)
DECAY = Simulation(end=10, steps=50)
#: a second dose of A at 5 min: the output has the states before and after it
REDOSE = Simulation(end=10, steps=50, changes=[Change([5.0], {"[A]": Q(1.0, "mM")})])
RATES = Dimension("dim_k", values={"k": Q([0.1, 0.2], "1/min")})


def run(simulation: Simulation, dimensions: Sequence[Dimension] = ()) -> ScanResult:
    return Simulator(n_workers=1).run(TOY, Scan(simulation, list(dimensions)))


def test_the_dimension_of_a_scan() -> None:
    assert single_dimension(run(DECAY, [RATES])) == "dim_k"


def test_a_simulation_has_no_dimension() -> None:
    with pytest.raises(ValueError, match="one dimension"):
        single_dimension(run(DECAY))


def test_the_timecourse_of_a_point_without_changes() -> None:
    res = run(DECAY, [RATES])
    time, values = point_timecourse(res, "[A]", dim_k=1)
    assert np.array_equal(time, res.ds["time"].values)
    assert np.array_equal(values, res["[A]"].isel(dim_k=1).values)


@pytest.mark.parametrize("name", ["res", "key"])
def test_a_dimension_can_have_the_name_of_an_argument(name: str) -> None:
    res = run(DECAY, [Dimension(name, values={"k": Q([0.1, 0.2], "1/min")})])
    time, values = point_timecourse(res, "[A]", **{name: 1})
    assert np.array_equal(time, res.ds["time"].values)
    assert np.array_equal(values, res["[A]"].isel({name: 1}).values)


def test_the_state_before_a_change_is_left_out() -> None:
    res = run(REDOSE, [RATES])
    assert res.ds.sizes["time"] == 52
    time, values = point_timecourse(res, "[A]", dim_k=0)
    assert time.size == 51
    assert np.all(np.diff(time) > 0)
    assert values[time == 5.0] == pytest.approx([1.0])


def test_a_simulation_needs_no_position() -> None:
    time, _ = point_timecourse(run(REDOSE), "[A]")
    assert time.size == 51


def test_a_ragged_result_loses_its_padding() -> None:
    lengths = Dimension(
        "dim_sim",
        simulations={
            "short": Simulation(end=5, steps=10),
            "long": Simulation(end=10, steps=20),
        },
    )
    res = run(DECAY, [lengths])
    assert res.ragged
    time, values = point_timecourse(res, "[A]", dim_sim=0)
    assert time.size == values.size == 11
    assert time[-1] == 5.0


def test_every_dimension_needs_a_position() -> None:
    res = run(DECAY, [RATES])
    with pytest.raises(ValueError, match="dim_k"):
        point_timecourse(res, "[A]")
    with pytest.raises(ValueError, match="dim_x"):
        point_timecourse(res, "[A]", dim_k=0, dim_x=0)


def test_a_variable_which_is_no_timecourse() -> None:
    with pytest.raises(ValueError, match="no timecourse"):
        point_timecourse(run(DECAY, [RATES]), "missing", dim_k=0)


def test_the_state_after_a_change_has_every_time_once() -> None:
    res = run(REDOSE)
    assert res.ds.sizes["time"] == 52
    after = state_after_changes(res)
    time = after.ds["time"].values
    assert time.size == 51
    assert np.all(np.diff(time) > 0)
    assert after["[A]"].values[time == 5.0] == pytest.approx([1.0])


def test_the_state_after_a_change_is_the_state_of_the_point_timecourse() -> None:
    res = run(REDOSE, [RATES])
    after = state_after_changes(res)
    for k in range(2):
        time, values = point_timecourse(res, "[A]", dim_k=k)
        assert np.array_equal(after.ds["time"].values, time)
        assert np.allclose(after["[A]"].isel(dim_k=k).values, values)


def test_a_result_without_a_change_is_unchanged() -> None:
    res = run(DECAY, [RATES])
    after = state_after_changes(res)
    assert np.array_equal(after.ds["time"].values, res.ds["time"].values)
    assert np.allclose(after["[A]"].values, res["[A]"].values)


def test_the_state_after_a_change_of_a_ragged_result() -> None:
    lengths = Dimension(
        "dim_sim",
        simulations={
            "short": Simulation(end=5, steps=10),
            "long": Simulation(
                end=10, steps=20, changes=[Change([3.0], {"[A]": Q(1.0, "mM")})]
            ),
        },
    )
    res = run(DECAY, [lengths])
    assert res.ragged
    after = state_after_changes(res)
    for k in range(2):
        time, values = point_timecourse(res, "[A]", dim_sim=k)
        got = after.isel(dim_sim=k)
        finite = np.isfinite(got.ds["time"].values)
        assert np.allclose(got.ds["time"].values[finite][: time.size], time)
        assert np.allclose(got["[A]"].values[finite][: time.size], values)


def test_the_state_after_a_change_of_a_scan_of_two_dimensions() -> None:
    doses = Dimension("dim_dose", values={"[A]": Q([1.0, 2.0, 3.0], "mM")})
    res = run(REDOSE, [RATES, doses])
    after = state_after_changes(res)
    assert after["[A]"].shape[:2] == (2, 3)
    for i in range(2):
        for j in range(3):
            time, values = point_timecourse(res, "[A]", dim_k=i, dim_dose=j)
            assert np.array_equal(after.ds["time"].values, time)
            assert np.allclose(after["[A]"].isel(dim_k=i, dim_dose=j).values, values)


def test_the_start_and_end_values_on_a_grid() -> None:
    res = run(DECAY, [RATES])
    assert not res.ragged
    start = start_values(res, "[A]")
    end = end_values(res, "[A]")
    assert start.dims == end.dims == ("dim_k",)
    np.testing.assert_array_equal(start.values, res["[A]"].isel(time=0).values)
    np.testing.assert_array_equal(end.values, res["[A]"].isel(time=-1).values)


def test_the_start_and_end_values_of_the_steps_of_the_integrator() -> None:
    """A scan without steps is ragged, the end is the last time of every point."""
    res = run(Simulation(end=10), [RATES])
    assert res.ragged
    start = start_values(res, "[A]")
    end = end_values(res, "[A]")
    assert start.dims == end.dims == ("dim_k",)
    for k in range(2):
        _, values = point_timecourse(res, "[A]", dim_k=k)
        assert start.values[k] == values[0]
        assert end.values[k] == values[-1]
    # the faster decay ends lower
    assert end.values[1] < end.values[0]


def test_the_end_values_of_simulations_of_different_lengths() -> None:
    lengths = Dimension(
        "dim_sim",
        simulations={
            "short": Simulation(end=5, steps=10),
            "long": Simulation(end=10, steps=20),
        },
    )
    res = run(DECAY, [lengths])
    end = end_values(res, "[A]")
    for k in range(2):
        _, values = point_timecourse(res, "[A]", dim_sim=k)
        assert end.values[k] == values[-1]
    assert np.all(np.isfinite(end.values))


def test_the_start_and_end_values_need_a_timecourse() -> None:
    res = run(DECAY, [RATES])
    with pytest.raises(ValueError, match="no timecourse"):
        start_values(res, "missing")
    with pytest.raises(ValueError, match="no timecourse"):
        end_values(res, "missing")
