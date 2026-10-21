"""Repeated doses add to their dose states (`PKPDSimulationExperiment.add_doses`)."""

import math

import pytest
from sbmlsim import Q
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.simulation import Change, Simulation

from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.run import create_simulator
from pkpdsim_toy import DATA_PATHS, MODEL, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyExperiment


class ToyDoses(ToyExperiment):
    """Two doses of A, 0.5 mmole each, at 2 and 4 min."""

    def simulations(self) -> dict[str, Simulation]:
        """Get the decay of A with the two doses."""
        doses = self.add_doses([2.0, 4.0], {"A": Q(500, "µmole")})
        return {"decay": Simulation(end=6, steps=600, changes=[doses])}


def experiment(cls: type[PKPDSimulationExperiment]) -> PKPDSimulationExperiment:
    """Get the experiment with its models loaded and initialized."""
    runner = ExperimentRunner(
        experiment_classes=[cls],
        base_path=TOY_PATH,
        data_path=DATA_PATHS,
        simulator=None,
    )
    exp = runner.experiments[cls.__name__]
    assert isinstance(exp, PKPDSimulationExperiment)
    return exp


def test_a_dose_adds_in_the_unit_of_the_model() -> None:
    """The dose is a formula of the dose state, converted to its unit (mmole)."""
    simulation = experiment(ToyDoses)._simulations["decay"]
    assert isinstance(simulation, Simulation)
    change = simulation.changes[0]
    assert isinstance(change, Change)
    assert change.times == (2.0, 4.0)
    assert dict(change.values) == {"A": "A + 0.5"}


def test_the_doses_keep_the_rest_of_the_earlier_ones() -> None:
    """At a dose the amount jumps by the dose from what is left."""
    exp = experiment(ToyDoses)
    exp._run_tasks(create_simulator(MODEL.simulator))
    result = exp._results["task_decay"]
    times = result["time"].to_numpy().ravel()
    a = result["[A]"].to_numpy().ravel()  # Vc is 1 l: the amount in mmole

    def at(t: float) -> float:
        return float(a[abs(times - t).argmin()])

    # A(t) = A(t0) * exp(-0.2 (t - t0)) between the doses, A(0) = 1 mmole
    decay = math.exp(-0.2 * 2.0)
    a2 = 1.0 * decay + 0.5
    a4 = a2 * decay + 0.5
    assert at(6.0) == pytest.approx(a4 * decay, rel=1e-3)
    # a dose which replaced the state would end at 0.5 * decay
    assert at(6.0) > 0.5 * decay * 1.5


def test_doses_need_a_loaded_model() -> None:
    """The unit of a dose state is the one of a loaded model."""
    exp = experiment(ToyDoses)
    with pytest.raises(ValueError, match="not loaded"):
        exp.add_doses([1.0], {"A": Q(1, "mmole")}, model="other")
