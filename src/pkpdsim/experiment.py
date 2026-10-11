"""The base class of the simulation experiments of a model.

The base of the experiments of a model derives from `PKPDSimulationExperiment`
and binds it to the model object; the experiments of the model derive from that
base:

```python
# not executed
class HctzSimulationExperiment(PKPDSimulationExperiment):
    pkpd_model = MODEL
    selections = ("time", "[Cve_hctz]", "Aurine_hctz")
    Mr: ClassVar[MolecularWeights] = MolecularWeights(hctz=Q(297.74, "g/mole"))
```

The base simulates the model of the model object (`models`), creates a task per
simulation (`tasks`), adds the selections to every task (`data`), loads the
digitized data of the study (`load_dataframe`) and has the classes of impaired
organ function (`pkpdsim.impairment`) and the labels of the pharmacokinetic
parameters (`pkpdsim.pk`) as defaults. A model overrides what it does
differently.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar

from sbmlsim.data import Data, load_pkdb_dataframe, load_pkdb_dataframes_by_substance
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task

from pkpdsim import impairment, pk

if TYPE_CHECKING:
    import pandas as pd
    from sbmlsim.plot import Figure
    from sbmlsim.units import Quantity, UnitRegistry

    from pkpdsim.model import PKPDModel


class PKPDSimulationExperiment(SimulationExperiment):
    """The base of the simulation experiments of a model.

    Attributes:
        pkpd_model: the model object of the model, set once in the base of its
            experiments.
        selections: the selections of every task, see `data`.
        Mr: the molecular weights of the substances, for the conversion of the
            data (a `NamedTuple` of the model).
        labels: the labels of the observables by id.
        units: the units of the observables by id.
        renal_map: the factor of the renal function by class.
        renal_colors: the color of a class of the renal function.
        cirrhosis_map: the degree of cirrhosis by class.
        cirrhosis_colors: the color of a class of the cirrhosis.
        cardiac_map: the factor of the cardiac output by class.
        cardiac_colors: the color of a class of the cardiac function.
        pk_labels: the label of a pharmacokinetic parameter.
        pk_units: the unit of a pharmacokinetic parameter.
    """

    pkpd_model: ClassVar[PKPDModel]
    selections: ClassVar[Sequence[str]] = ()
    Mr: ClassVar[tuple[Quantity, ...]]
    labels: ClassVar[Mapping[str, str]]
    units: ClassVar[Mapping[str, str]]
    renal_map: ClassVar[Mapping[str, float]] = impairment.RENAL_FUNCTION
    renal_colors: ClassVar[Mapping[str, str]] = impairment.RENAL_COLORS
    cirrhosis_map: ClassVar[Mapping[str, float]] = impairment.CIRRHOSIS
    cirrhosis_colors: ClassVar[Mapping[str, str]] = impairment.CIRRHOSIS_COLORS
    cardiac_map: ClassVar[Mapping[str, float]] = impairment.CARDIAC
    cardiac_colors: ClassVar[Mapping[str, str]] = impairment.CARDIAC_COLORS
    pk_labels: ClassVar[Mapping[str, str]] = pk.PK_LABELS
    pk_units: ClassVar[Mapping[str, str]] = pk.PK_UNITS

    def __init__(
        self,
        sid: str | None = None,
        base_path: Path | None = None,
        data_path: Path | Iterable[Path] | None = None,
        ureg: UnitRegistry | None = None,
        **kwargs: object,
    ) -> None:
        """Create the experiment, see `SimulationExperiment`.

        Raises:
            TypeError: if no base of the experiment sets `pkpd_model`.
        """
        if not hasattr(type(self), "pkpd_model"):
            raise TypeError(
                f"{type(self).__name__} has no model object: set `pkpd_model` in "
                "the base of the experiments of the model."
            )
        super().__init__(
            sid=sid, base_path=base_path, data_path=data_path, ureg=ureg, **kwargs
        )

    def models(self) -> dict[str, AbstractModel | Path]:
        """Get the SBML model of the model object, without changes."""
        return {
            "model": AbstractModel(
                source=self.pkpd_model.model_path,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    def tasks(self) -> dict[str, Task]:
        """Get a task of the model for every simulation."""
        return {
            f"task_{key}": Task(model="model", simulation=key)
            for key in self.simulations()
        }

    def default_changes(self) -> Mapping[str, Quantity | float]:
        """Get the changes of every simulation of the experiments, none.

        The parameters are part of the model and are not changed in the
        simulations: a parameter which is fitted must not be set by a
        simulation, it would override the value of the fit.
        """
        return {}

    def data(self) -> dict[str, Data]:
        """Add the selections to every task; returns no further data."""
        self.add_selections_data(selections=self.selections)
        return {}

    def _data_paths(self) -> list[Path]:
        """Get the data paths, which loading the data of a study needs.

        Raises:
            ValueError: if the experiment has no `data_path`.
        """
        if self.data_path is None:
            raise ValueError(f"No 'data_path' set on '{self.__class__.__name__}'.")
        return self.data_path

    def load_dataframe(
        self,
        fig_id: str,
        study: str | None = None,
        *,
        na_values: Sequence[str] | None = None,
    ) -> pd.DataFrame:
        """Load the digitized data of a figure or table of a study.

        Args:
            fig_id: the id of the figure or table in the study, e.g. `Fig1`.
            study: the study, the id of the experiment by default.
            na_values: additional strings of missing values, see
                `pandas.read_csv`.

        Raises:
            ValueError: if the experiment has no `data_path`.
        """
        return load_pkdb_dataframe(
            f"{study or self.sid}_{fig_id}",
            data_path=self._data_paths(),
            na_values=na_values,
        )

    def load_dataframes_by_substance(
        self, fig_id: str, study: str | None = None
    ) -> dict[str, pd.DataFrame]:
        """Load the digitized data of a figure or table by substance.

        Args:
            fig_id: the id of the figure or table in the study, e.g. `Fig1`.
            study: the study, the id of the experiment by default.

        Raises:
            ValueError: if the experiment has no `data_path`.
        """
        return load_pkdb_dataframes_by_substance(
            f"{study or self.sid}_{fig_id}", data_path=self._data_paths()
        )

    @staticmethod
    def figures_by_sid(*figures: Figure) -> dict[str, Figure]:
        """Get figures by their id, which every figure of the experiments has.

        Raises:
            ValueError: for a figure without an id.
        """
        by_sid: dict[str, Figure] = {}
        for figure in figures:
            if figure.sid is None:
                raise ValueError(f"Figure without an id: '{figure}'")
            by_sid[figure.sid] = figure
        return by_sid
