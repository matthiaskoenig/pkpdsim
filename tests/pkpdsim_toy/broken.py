"""Experiments of the toy model which fail, for the tests of the failures of a run.

They are not below `experiments/`, so the registry of the toy model and the
checks of the registry do not see them; `BROKEN` is a registry with them.
"""

from types import MappingProxyType

from matplotlib.figure import Figure as FigureMPL
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.simulation import Simulation

from pkpdsim_toy.experiments.decay import ToyDecay


class ToyFigureFails(ToyDecay):
    """The decay, whose own figures fail when they are drawn."""

    def figures_mpl(self) -> dict[str, FigureMPL]:
        """Fail when the figures are drawn."""
        raise RuntimeError("the figures of ToyFigureFails fail")


class ToyCannotInitialize(ToyDecay):
    """The decay, whose simulations cannot be created."""

    def simulations(self) -> dict[str, Simulation]:
        """Fail when the experiment is initialized."""
        raise RuntimeError("the simulations of ToyCannotInitialize fail")


class ToyRefusedUnit(ToyDecay):
    """The decay, whose unit of the data cannot be defined (a misspelled unit)."""

    data_units = MappingProxyType({"IU": "IU_toy_refused = 0.6 * microgramm"})


#: a registry with the experiments which fail
BROKEN: dict[str, list[type[SimulationExperiment]]] = {
    "all": [ToyDecay, ToyFigureFails, ToyCannotInitialize, ToyRefusedUnit]
}
