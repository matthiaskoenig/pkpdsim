"""The decay of A in men, in women and in a study of both."""

from collections.abc import Mapping
from typing import ClassVar

import pandas as pd
from sbmlsim.data import Data, DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Curve, Figure, Plot, SubPlot
from sbmlsim.simulation import Simulation

from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.metadata import Dosing, Health, PKPDMappingMetaData, Route, Tissue
from pkpdsim.sex import Sex
from pkpdsim_sexes import FEMALE_A0, MODEL

#: the times of the data [min]
TIMES = [0.0, 2.0, 4.0, 6.0, 8.0]
#: A in men at the times [mM], the decay at the rate 0.2/min from 1 mM
A_MEN = [1.0, 0.67, 0.449, 0.301, 0.202]


class SexesExperiment(PKPDSimulationExperiment):
    """The base of the experiments of the toy model per sex."""

    pkpd_model = MODEL
    selections = ("time", "[A]")
    labels: ClassVar[Mapping[str, str]] = {"time": "time", "[A]": "A"}
    units: ClassVar[Mapping[str, str]] = {"time": "min", "[A]": "mM"}

    def simulations(self) -> dict[str, Simulation]:
        """Get the decay over 10 min."""
        return {"decay": Simulation(end=10)}

    def figures(self) -> dict[str, Figure]:
        """Get the figure of the concentration of A in every task."""
        plot = Plot(sid="decay", xaxis=Axis("time", unit="min"), yaxis=Axis("A"))
        for key in self._simulations:
            plot.curves.append(
                Curve(
                    x=Data(index="time", task=f"task_{key}"),
                    y=Data(index="[A]", task=f"task_{key}"),
                    name=f"A {key}",
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

    def _decay_data(self, factor: float) -> DataSet:
        """Get the data of the decay, A of the men times a factor."""
        df = pd.DataFrame(
            {
                "time": TIMES,
                "time_unit": "min",
                "A": [factor * a for a in A_MEN],
                "A_unit": "mM",
                "count": 4,
            }
        )
        return DataSet.from_df(df, ureg=self.ureg)

    def _mapping(self, dataset: str, task: str) -> FitMapping:
        """Get the mapping of the data of A to a task."""
        return FitMapping(
            self,
            reference=FitData(
                self, dataset=dataset, xid="time", yid="A", count="count"
            ),
            observable=FitData(self, task=task, xid="time", yid="[A]"),
            metadata=PKPDMappingMetaData(
                tissue=Tissue.PLASMA,
                route=Route.IV,
                dosing=Dosing.SINGLE,
                health=Health.HEALTHY,
            ),
        )


class SexesMen(SexesExperiment):
    """The decay in men, with data of A."""

    def datasets(self) -> dict[str, DataSet]:
        """Get the data of the men."""
        return {"data_decay": self._decay_data(1.0)}

    def fit_mappings(self) -> dict[str, FitMapping]:
        """Get the mapping of the data of the men."""
        return {"fm_decay_A": self._mapping("data_decay", "task_decay")}


class SexesWomen(SexesExperiment):
    """The decay in women, on the female build, with data of A."""

    sex = Sex.FEMALE

    def datasets(self) -> dict[str, DataSet]:
        """Get the data of the women."""
        return {"data_decay": self._decay_data(FEMALE_A0)}

    def fit_mappings(self) -> dict[str, FitMapping]:
        """Get the mapping of the data of the women."""
        return {"fm_decay_A": self._mapping("data_decay", "task_decay")}


class SexesGroups(SexesExperiment):
    """A study of a group of men and a group of women, a simulation each."""

    simulation_sexes: ClassVar[Mapping[str, Sex]] = {"women": Sex.FEMALE}

    def simulations(self) -> dict[str, Simulation]:
        """Get the decay of the men and of the women."""
        return {"men": Simulation(end=10), "women": Simulation(end=10)}
