"""The decay of A at the rate of the model and at half of it."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from sbmlsim import Q
from sbmlsim.data import Data, DataSet
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.fit.objects import MappingMetaData
from sbmlsim.model import AbstractModel
from sbmlsim.plot import Axis, Curve, Figure, Plot, SubPlot
from sbmlsim.simulation import Simulation
from sbmlsim.task import Task

from pkpdsim_toy import MODEL_PATH


@dataclass(kw_only=True)
class ToyMappingMetaData(MappingMetaData):
    """The metadata of a fit mapping of the toy model."""

    route: str


class ToyExperiment(SimulationExperiment):
    """The base of the toy experiments: the model, the decay and its figure."""

    #: the rate of the decay [1/min], `None` for the rate of the model
    rate: float | None = None

    def models(self) -> dict[str, AbstractModel | Path]:
        """Get the toy model."""
        return {"model": AbstractModel(source=MODEL_PATH)}

    def simulations(self) -> dict[str, Simulation]:
        """Get the decay over 10 min."""
        changes = None if self.rate is None else {"k": Q(self.rate, "1/min")}
        return {"decay": Simulation(end=10, steps=50, preinit_changes=changes)}

    def tasks(self) -> dict[str, Task]:
        """Get the task of the decay."""
        return {"task_decay": Task(model="model", simulation="decay")}

    def data(self) -> dict[str, Data]:
        """Get the concentration of A."""
        return {"A": Data(index="[A]", task="task_decay")}

    def figures(self) -> dict[str, Figure]:
        """Get the figure of the concentration of A."""
        plot = Plot(
            sid="decay",
            xaxis=Axis("time", unit="min"),
            yaxis=Axis("A", unit="mM"),
        )
        plot.curves.append(
            Curve(
                x=Data(index="time", task="task_decay"),
                y=Data(index="[A]", task="task_decay"),
                name="A",
            )
        )
        return {
            "fig_decay": Figure(
                experiment=self,
                sid="fig_decay",
                num_rows=1,
                num_cols=1,
                subplots=[SubPlot(plot=plot, row=1, col=1)],
            )
        }


class ToyDecay(ToyExperiment):
    """The decay at the rate of the model, with data of A."""

    def datasets(self) -> dict[str, DataSet]:
        """Get the data of the decay."""
        df = pd.DataFrame(
            {
                "time": [0.0, 2.0, 4.0, 6.0, 8.0],
                "time_unit": "min",
                "A": [1.0, 0.67, 0.45, 0.30, 0.20],
                "A_unit": "mM",
                "count": 4,
            }
        )
        return {"data_decay": DataSet.from_df(df, ureg=self.ureg)}

    def fit_mappings(self) -> dict[str, FitMapping]:
        """Get the mapping of the data of A to the simulation."""
        return {
            "fm_decay_A": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset="data_decay",
                    xid="time",
                    yid="A",
                    count="count",
                ),
                observable=FitData(self, task="task_decay", xid="time", yid="[A]"),
                metadata=ToyMappingMetaData(route="iv"),
            )
        }


class ToyDecaySlow(ToyExperiment):
    """The decay at half of the rate of the model."""

    rate = 0.1
