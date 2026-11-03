"""The registry of the experiments of the credibility test model."""

from sbmlsim.experiment import SimulationExperiment

from pkpdsim_credibility.experiments import CredA, CredB

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {"studies": [CredA, CredB]}
EXPERIMENTS["all"] = EXPERIMENTS["studies"]
