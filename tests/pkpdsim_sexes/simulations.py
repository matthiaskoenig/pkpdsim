"""The registry of the experiments of the toy model per sex."""

from sbmlsim.experiment import SimulationExperiment

from pkpdsim_sexes.experiments.decay import SexesGroups, SexesMen, SexesWomen

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {
    "studies": [SexesMen, SexesWomen, SexesGroups],
}
EXPERIMENTS["all"] = EXPERIMENTS["studies"]
