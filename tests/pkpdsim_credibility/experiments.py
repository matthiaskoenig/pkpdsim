"""Two experiments of the decay of A with the mappings of every kind of evidence."""

import pandas as pd
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping

from pkpdsim.metadata import Dosing, Health, PKPDMappingMetaData, Route, Tissue
from pkpdsim_toy.experiments.decay import ToyExperiment

TIMES = [0.0, 2.0, 4.0, 6.0, 8.0]
DECAY = [1.0, 0.67, 0.45, 0.30, 0.20]


def _data(values: list[float], ureg: object) -> DataSet:
    """Get the data of A over the times."""
    df = pd.DataFrame(
        {
            "time": TIMES,
            "time_unit": "min",
            "A": values,
            "A_unit": "mM",
            "count": 4,
            # the group and the intervention of the study tables
            "group": "all",
            "intervention": "A1",
        }
    )
    return DataSet.from_df(df, ureg=ureg)  # ty: ignore[invalid-argument-type]


class CredExperiment(ToyExperiment):
    """The base of the experiments: the decay, the data and the mappings."""

    def datasets(self) -> dict[str, DataSet]:
        """Get the decay, three times the decay and zeros."""
        return {
            "data_decay": _data(DECAY, self.ureg),
            "data_limit": _data([3 * v for v in DECAY], self.ureg),
            "data_zero": _data([0.0] * len(DECAY), self.ureg),
        }

    def mapping(
        self, dataset: str, route: Route, health: Health = Health.HEALTHY
    ) -> FitMapping:
        """Get the mapping of a dataset of A to the simulation."""
        return FitMapping(
            self,
            reference=FitData(
                self, dataset=dataset, xid="time", yid="A", count="count"
            ),
            observable=FitData(self, task="task_decay", xid="time", yid="[A]"),
            metadata=PKPDMappingMetaData(
                tissue=Tissue.PLASMA,
                route=route,
                dosing=Dosing.SINGLE,
                health=health,
            ),
        )


class CredA(CredExperiment):
    """The calibration, validation and outlier mappings."""

    def fit_mappings(self) -> dict[str, FitMapping]:
        """Get the mappings of the study A."""
        return {
            "fm_train": self.mapping("data_decay", Route.IV),
            "fm_dose": self.mapping("data_decay", Route.IV),
            "fm_renal": self.mapping("data_decay", Route.IV, Health.RENAL_IMPAIRMENT),
            "fm_limit": self.mapping("data_limit", Route.IV),
            "fm_subject": self.mapping("data_decay", Route.IV),
            "fm_outlier": self.mapping("data_decay", Route.IV),
        }


class CredB(CredExperiment):
    """A mapping with the key of the study A and a mapping of zeros."""

    def fit_mappings(self) -> dict[str, FitMapping]:
        """Get the mappings of the study B."""
        return {
            "fm_dose": self.mapping("data_decay", Route.IV),
            "fm_zero": self.mapping("data_zero", Route.PO),
        }
