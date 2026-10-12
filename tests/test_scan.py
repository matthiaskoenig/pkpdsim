"""The results of the scans of a model, on the toy model."""

from collections.abc import Sequence

import numpy as np
import pytest
from sbmlsim import Q
from sbmlsim.model import AbstractModel
from sbmlsim.result import ScanResult
from sbmlsim.simulation import Change, Dimension, Scan, Simulation
from sbmlsim.simulator import Simulator

from pkpdsim.scan import point_timecourse, single_dimension
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
