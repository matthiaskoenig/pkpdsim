"""The decay of A at the rate of the model and at half of it."""

from collections.abc import Mapping
from typing import ClassVar

import pandas as pd
from sbmlsim import Q
from sbmlsim.data import Data, DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Curve, Figure, Plot, SubPlot
from sbmlsim.simulation import Simulation

from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.metadata import Dosing, Health, PKPDMappingMetaData, Route, Tissue
from pkpdsim_toy import MODEL


class ToyExperiment(PKPDSimulationExperiment):
    """The base of the toy experiments: the decay and its figure."""

    pkpd_model = MODEL
    selections = ("time", "[A]")
    labels: ClassVar[Mapping[str, str]] = {"time": "time", "[A]": "A"}
    units: ClassVar[Mapping[str, str]] = {"time": "min", "[A]": "mM"}

    #: the rate of the decay [1/min], `None` for the rate of the model
    rate: float | None = None

    def simulations(self) -> dict[str, Simulation]:
        """Get the decay over 10 min."""
        changes = None if self.rate is None else {"k": Q(self.rate, "1/min")}
        return {"decay": Simulation(end=10, steps=50, preinit_changes=changes)}

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
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.IV,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                ),
            )
        }


class ToyDecaySlow(ToyExperiment):
    """The decay at half of the rate of the model."""

    rate = 0.1
