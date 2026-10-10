"""The registry of the experiments of the toy model."""

from sbmlsim.experiment import SimulationExperiment

from pkpdsim_toy.experiments.decay import ToyDecay, ToyDecaySlow

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {
    "studies": [ToyDecay],
    "misc": [ToyDecaySlow],
}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"]

#: a registry without `ToyDecaySlow`, for the tests of the check of the registry
PARTIAL: dict[str, list[type[SimulationExperiment]]] = {"all": [ToyDecay]}
