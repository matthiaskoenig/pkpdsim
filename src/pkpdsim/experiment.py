"""The base class of the simulation experiments of a model.

The base of the experiments of a model derives from `PKPDSimulationExperiment`
and binds it to the model object; the experiments of the model derive from that
base:

```python
# not executed
from sbmlsim import Q


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

The keys of an experiment are SIds (sbmlsim checks them): `to_sid` gives the SId
of a text, e.g. of a label of the data, and the method `to_sid` of the base
checks that no two texts of an experiment give the same SId. A unit which the
data use and sbmlsim does not define (an international unit is specific to a
substance) is defined under a name of the model and the substance, `data_units`.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from types import MappingProxyType
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

#: a decimal point between two digits, which `to_sid` writes as `p`
_DECIMAL_POINT = re.compile(r"(?<=[0-9])\.(?=[0-9])")
#: a character which an SId does not allow
_NOT_IN_SID = re.compile(r"[^A-Za-z0-9_]")
#: the metric prefixes a unit of the data can have, the longest first
_PREFIXES: tuple[str, ...] = (
    "da", "Y", "Z", "E", "P", "T", "G", "M", "k", "h", "d", "c", "m",
    "µ", "μ", "u", "n", "p", "f", "a", "z", "y",
)  # fmt: skip


def to_sid(text: str) -> str:
    """Get the SId of a text, e.g. of a label of the data, as a key of an experiment.

    A decimal point between two digits becomes `p` (`0.5mg` is `0p5mg`) and
    every other character which an SId does not allow becomes `_`
    (`po 2.5 mg` is `po_2p5_mg`, `µg/ml` is `_g_ml`). One character becomes
    one character, and an SId stays as it is. The function maps characters, so
    a part of a key gives the same part of the SId of the key
    (`f"task_po_{to_sid(dose)}"`). Whether a whole key is a valid SId, i.e. not
    empty and not starting with a digit, is not checked here (the text of a part
    of a key may start with a digit): sbmlsim checks every key when the
    experiment initializes and names all invalid ones.

    Args:
        text: the text.

    Returns:
        The SId.
    """
    return _NOT_IN_SID.sub("_", _DECIMAL_POINT.sub("p", text))


def define_unit(ureg: UnitRegistry, definition: str) -> str:
    """Define a unit in a registry once, e.g. the international unit of a substance.

    The registry of sbmlsim, one per process, ignores a redefinition, so a
    second model which defined a name with another value would silently get
    the first one. A name which is defined already must therefore have the same
    value (compared in base units); a name of the model and the substance, e.g.
    `IU_ramipril_renin`, keeps the models apart.

    Args:
        ureg: the unit registry, the one of the experiments.
        definition: `"<name> = <expression>"` in the syntax of pint, e.g.
            `"IU_ramipril_renin = 0.6 * microgram"` or `"IU_apixaban_antixa = count"`.

    Returns:
        The name of the unit.

    Raises:
        ValueError: if the definition is not `<name> = <expression>` with an
            identifier as the name, or the name is defined with another value.
    """
    name, separator, expression = (part.strip() for part in definition.partition("="))
    if not separator or not name.isidentifier() or not expression:
        raise ValueError(
            f"A unit is defined as '<name> = <expression>', not '{definition}'."
        )
    value = ureg.Quantity(expression).to_base_units()
    if name in ureg:
        defined = ureg.Quantity(1, name).to_base_units()
        if defined.units != value.units or not math.isclose(
            defined.magnitude, value.magnitude, rel_tol=1e-12
        ):
            raise ValueError(
                f"The unit '{name}' is defined as {defined} already, not as "
                f"{value} ({expression}); give the unit of another substance or "
                "model another name."
            )
        return name
    ureg.define(f"{name} = {expression}")
    return name


def rename_units(df: pd.DataFrame, names: Mapping[str, str]) -> pd.DataFrame:
    """Rename units in the unit columns of a table of the data.

    A name is replaced where it is a unit of a unit expression, alone or after a
    metric prefix: with `{"IU": "IU_x"}`, `IU/ml` is `IU_x/ml` and `µIU/ml` is
    `µIU_x/ml`; `IU_per_ml` and `2IU` stay. The unit columns are `unit` and
    every column whose name ends with `_unit`; a missing value stays missing.

    Args:
        df: the table, changed in place.
        names: the name in the data -> the new name.

    Returns:
        The table.
    """
    if not names:
        return df
    prefixes = "|".join(re.escape(prefix) for prefix in _PREFIXES)
    patterns = [
        (
            re.compile(rf"(?<!\w)({prefixes})?{re.escape(old)}(?!\w)"),
            r"\g<1>" + new.replace("\\", "\\\\"),
        )
        for old, new in names.items()
    ]

    def rename(value: object) -> object:
        if not isinstance(value, str):
            return value
        for pattern, replacement in patterns:
            value = pattern.sub(replacement, value)
        return value

    for column in df.columns:
        if column == "unit" or str(column).endswith("_unit"):
            df[column] = df[column].map(rename)
    return df


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
        data_units: the units of the data which the experiment reads under a
            name of its own, see `define_unit`.
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
    #: the units of the data which the experiment reads under a name of its own:
    #: the name in the data -> the definition of the unit, "<name> = <expression>",
    #: e.g. {"IU": "IU_ramipril_renin = 0.6 * microgram"}, see `define_unit`
    data_units: ClassVar[Mapping[str, str]] = MappingProxyType({})

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
            ValueError: for a definition of `data_units` which `define_unit` refuses.
        """
        if not hasattr(type(self), "pkpd_model"):
            raise TypeError(
                f"{type(self).__name__} has no model object: set `pkpd_model` in "
                "the base of the experiments of the model."
            )
        super().__init__(
            sid=sid, base_path=base_path, data_path=data_path, ureg=ureg, **kwargs
        )
        #: the text of every SId the experiment gave, see `to_sid`
        self._sid_texts: dict[str, str] = {}
        # the units of the data are defined before `initialize` reads the data
        self._unit_names: dict[str, str] = {
            name: define_unit(self.ureg, definition)
            for name, definition in self.data_units.items()
        }

    def to_sid(self, text: str) -> str:
        """Get the SId of a text as a key of the experiment, see `to_sid`.

        The keys of an experiment are unique across its datasets, simulations,
        tasks, figures and fit mappings, so two texts which give one SId in an
        experiment are an error: one entry would replace the other.

        Args:
            text: the text, e.g. a label of the data.

        Returns:
            The SId.

        Raises:
            ValueError: if another text of the experiment gave the same SId.
        """
        sid = to_sid(text)
        first = self._sid_texts.setdefault(sid, text)
        if first != text:
            raise ValueError(
                f"'{first}' and '{text}' are the same SId '{sid}' in the "
                f"experiment '{self.sid}'."
            )
        return sid

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
        """Get a task of the model for every simulation.

        `initialize` creates the simulations before the tasks, so they are not
        created a second time.
        """
        return {
            f"task_{key}": Task(model="model", simulation=key)
            for key in self._simulations
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
        na_values: Sequence[str] | Mapping[str, Sequence[str]] | None = None,
    ) -> pd.DataFrame:
        """Load the digitized data of a figure or table of a study.

        Args:
            fig_id: the id of the figure or table in the study, e.g. `Fig1`.
            study: the study, the id of the experiment by default.
            na_values: additional strings of missing values, for every column or
                by column, e.g. `{"time": ["NR"]}`, see `pandas.read_csv`.

        Raises:
            ValueError: if the experiment has no `data_path`.
        """
        df = load_pkdb_dataframe(
            f"{self.sid if study is None else study}_{fig_id}",
            data_path=self._data_paths(),
            na_values=na_values,
        )
        return rename_units(df, self._unit_names)

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
        frames = load_pkdb_dataframes_by_substance(
            f"{self.sid if study is None else study}_{fig_id}",
            data_path=self._data_paths(),
        )
        return {key: rename_units(df, self._unit_names) for key, df in frames.items()}

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
