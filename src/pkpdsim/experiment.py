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

The base simulates the model of the model object (`models`), the build of the
sex of the study (`sex`, the male build by default, the female build of a study
of women; `simulation_sexes` for a study with a group of each sex), creates a
task per simulation on the build of its sex (`tasks`), adds the selections to
every task (`data`), loads the digitized data of the study (`load_dataframe`) and has the classes of impaired
organ function (`pkpdsim.impairment`) and the labels of the pharmacokinetic
parameters (`pkpdsim.pk`) as defaults. A model overrides what it does
differently.

The keys of an experiment are SIds (sbmlsim checks them): `sbmlsim.experiment.to_sid`
gives the SId of a text, e.g. of a label of the data, and the method `to_sid` of the base
checks that no two texts of an experiment give the same SId. A unit which the
data use and sbmlsim does not define (an international unit is specific to a
substance) is defined under a name of the model and the substance, `data_units`
(`sbmlsim.units.define_unit`).

A run over the draws of the fitted parameters (`run_simulation <package> --uncertainty`)
gives every experiment the setting `uncertainty`, a
`sbmlsim.experiment.draws.UncertaintyRun`: the `initialize` of sbmlsim's
`SimulationExperiment` runs every simulation over the draws and turns the curves
of the simulations into bands (`sbmlsim.experiment.draws.apply_uncertainty`), so
an experiment needs no code of its own for it. A figure which an experiment draws
itself (`figures_mpl`) reads `uncertainty` and draws the bands of the draws, e.g.
with `sbmlsim.result.scan.point_band(result, ..., across=self.across())`; a
function of the results (`add_function`) is a band of the draws as well.

A dose after the first one adds to its dose state (`add_doses`): the rest of the
earlier doses which is not yet dissolved, absorbed or infused stays, which a
change of the dose state to the dose would discard (pkdb_models#105).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, ClassVar, override

from sbmlsim.data import Data, load_pkdb_dataframe, rename_units
from sbmlsim.experiment import SimulationExperiment, to_sid
from sbmlsim.model import AbstractModel
from sbmlsim.plot import Curve
from sbmlsim.simulation import Change
from sbmlsim.task import Task
from sbmlsim.units import define_unit

from pkpdsim import impairment, pk
from pkpdsim.sex import Sex

if TYPE_CHECKING:
    import pandas as pd
    from sbmlsim.plot import Figure, Plot
    from sbmlsim.units import Quantity, UnitRegistry

    from pkpdsim.model import PKPDModel


class PKPDSimulationExperiment(SimulationExperiment):
    """The base of the simulation experiments of a model.

    Attributes:
        pkpd_model: the model object of the model, set once in the base of its
            experiments.
        selections: the selections of every task, see `data`.
        Mr: the molecular weights of the substances, for the conversion of the
            data (a `NamedTuple` of the model). The base does not declare it:
            its type is the `NamedTuple` of the model, which the base of the
            model declares, and the narrowing of a declared class attribute is
            an invalid override of a mutable attribute (ty).
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
        sex: the sex of the subjects of the study, the build of the model
            object its simulations run (`PKPDModel.model_paths`): the male
            build by default, the female build for a study of women. A group
            whose majority is women is simulated on the female build, a mixed
            group without a majority of women on the male one (the study
            module states it). A model without a female build refuses a
            female experiment.
        simulation_sexes: the sex of single simulations which differs from
            `sex`, by the key of the simulation, e.g. the group of women of a
            study of a group of men and a group of women. The models of an
            experiment are loaded before its simulations are created, so the
            sexes are declared and not computed by `simulations`.
    """

    pkpd_model: ClassVar[PKPDModel]
    selections: ClassVar[Sequence[str]] = ()
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
    sex: ClassVar[Sex] = Sex.MALE
    simulation_sexes: ClassVar[Mapping[str, Sex]] = MappingProxyType({})

    def __init__(
        self,
        sid: str | None = None,
        base_path: Path | None = None,
        data_path: Path | Iterable[Path] | None = None,
        ureg: UnitRegistry | None = None,
        **kwargs: object,
    ) -> None:
        """Create the experiment, see `SimulationExperiment`.

        The units of the data (`data_units`) are defined by `initialize`, so a
        definition which is refused fails the initialization of this experiment,
        which a run reports, and not the creation, which stops a run of sbmlsim
        before any experiment ran.

        Raises:
            TypeError: if no base of the experiments sets `pkpd_model`. This is
                an error of the definition of the classes, which every
                experiment of the base has, so it stops a run.
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
        #: the names of the defined units of the data, see `data_unit_names`
        self._unit_names: Mapping[str, str] | None = None

    def initialize(self) -> None:
        """Define the units of the data and initialize the experiment.

        The units of `data_units` are defined before the data are read, see
        `SimulationExperiment.initialize`, which also runs the experiment over
        the draws of its setting `uncertainty`.

        Raises:
            ValueError: for a definition of `data_units` which `define_unit`
                refuses; the runner reports the experiment as failed.
        """
        try:
            self.data_unit_names()
        except ValueError as err:
            err.add_note(f"The experiment '{self.sid}' cannot be initialized.")
            raise
        super().initialize()

    def add_function(
        self,
        plot: Plot,
        task_id: str,
        y: Data,
        label: str | None,
        color: str,
        linestyle: str | None = "-",
    ) -> None:
        """Add a function of the results of a task over the time to a plot.

        A curve of the line width of the curves of the selections (`add_data`,
        2.0); in a run over the draws of the fitted parameters (`uncertainty`)
        the band of the draws instead, with the quantiles of the run: the
        function is no selection of a task, which the run makes a band of
        itself.

        Args:
            plot: the plot.
            task_id: the task whose results the function reads.
            y: the function, a `Data` of a function of the selections of the
                task.
            label: the label of the curve or the band, `None` for none.
            color: the color.
            linestyle: the line style of the curve, `None` for the default of
                sbmlsim (the curve has no line style of its own).
        """
        time = Data("time", task=task_id)
        run = self.uncertainty
        if run is not None:
            plot.band(
                x=time,
                y=y,
                across=run.id,
                quantiles=run.quantiles,
                name=label,
                color=color,
            )
            return
        style: dict[str, Any] = {"color": color, "linewidth": 2.0}
        if linestyle is not None:
            style["linestyle"] = linestyle
        plot.add_curve(Curve(x=time, y=y, label=label, **style))

    def data_unit_names(self) -> Mapping[str, str]:
        """Get the names of the units of the data, defining the units of `data_units` once.

        The loaders rename the units of the data with these names; a table
        which a model reads itself, e.g. the data of a figure of a scan, is
        renamed with them by `sbmlsim.data.rename_units`. `initialize` calls it first; a call
        without `initialize` defines the units on its first call.

        Returns:
            The name of a unit in the data -> the name of its definition,
            read-only.

        Raises:
            ValueError: for a definition which `define_unit` refuses.
        """
        if self._unit_names is None:
            self._unit_names = MappingProxyType(
                {
                    name: define_unit(definition, self.ureg)
                    for name, definition in self.data_units.items()
                }
            )
        return self._unit_names

    def to_sid(self, text: str) -> str:
        """Get the SId of a text as a key of the experiment, see `sbmlsim.experiment.to_sid`.

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

    @classmethod
    def sex_of(cls, simulation: str) -> Sex:
        """Get the sex of the subjects of a simulation: `simulation_sexes`, else `sex`.

        A class method, since the sexes are declared: the class methods of an
        experiment (the renal function of a group, `pkpdsim.pbpk`) resolve the
        sex of a simulation by its key as well.

        Args:
            simulation: the key of the simulation.

        Returns:
            The sex, whose build the simulation runs.

        Raises:
            ValueError: for a sex which is no `Sex`.
        """
        return Sex(cls.simulation_sexes.get(simulation, cls.sex))

    @override
    def variant_of(self, simulation: str) -> str:
        """Get the variant of the model a simulation runs on: the sex of its subjects.

        The variants of a `PKPDModel` are the sexes, so `pkpdsim-verify dump`
        stores the sex per task as `builds`.

        Args:
            simulation: the key of the simulation.

        Returns:
            The value of the sex (`male` or `female`).
        """
        return str(self.sex_of(simulation))

    def model_id(self, sex: Sex) -> str:
        """Get the key of the build of a sex among the models of the experiment.

        The build of `sex` is `model`, the model of an experiment without
        builds per sex; the build of another sex is `model_<sex>`, e.g.
        `model_female` in a study of men with a group of women.

        Args:
            sex: the sex.

        Returns:
            The key of the model.
        """
        sex = Sex(sex)
        return "model" if sex == Sex(self.sex) else f"model_{sex}"

    def models(self) -> dict[str, AbstractModel | Path]:
        """Get the builds of the model object which the experiment simulates.

        The build of `sex` as `model` and the build of every other sex of
        `simulation_sexes` as `model_<sex>` (`model_id`), without changes; a
        model without builds per sex has the one model `model`. Only the builds
        which the experiment simulates are loaded.

        Raises:
            ValueError: if the model object has no build of a sex of the
                experiment, naming the experiment.
        """
        sexes = dict.fromkeys(
            [Sex(self.sex), *(Sex(s) for s in self.simulation_sexes.values())]
        )
        models: dict[str, AbstractModel | Path] = {}
        for sex in sexes:
            try:
                path = self.pkpd_model.model_path_of(sex)
            except ValueError as err:
                raise ValueError(
                    f"The experiment '{self.sid}' simulates subjects of the sex "
                    f"'{sex}': {err}"
                ) from err
            models[self.model_id(sex)] = AbstractModel(
                source=path,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        return models

    def tasks(self) -> dict[str, Task]:
        """Get a task for every simulation, on the build of the sex of the simulation.

        `initialize` creates the simulations before the tasks, so they are not
        created a second time.

        Raises:
            ValueError: if `simulation_sexes` names a simulation which the
                experiment does not have.
        """
        unknown = sorted(set(self.simulation_sexes) - set(self._simulations))
        if unknown:
            raise ValueError(
                f"The simulation_sexes of the experiment '{self.sid}' name "
                f"simulations which it does not have: {unknown}; its "
                f"simulations are {sorted(self._simulations)}."
            )
        return {
            f"task_{key}": Task(model=self.model_id(self.sex_of(key)), simulation=key)
            for key in self._simulations
        }

    def add_doses(
        self,
        times: float | Quantity | Sequence[float | Quantity],
        doses: Mapping[str, Quantity],
        model: str = "model",
    ) -> Change:
        """Get doses at times which add to their dose states.

        A dose is given by its dose state, e.g. `PODOSE_hctz`, from which a
        rate rule dissolves, absorbs or infuses it. A change of the dose state
        to the dose (`Change(times, {"PODOSE_hctz": Q(10, "mg")})`) discards
        what is left of the earlier doses; the change of this method adds the
        dose (`"PODOSE_hctz + 10.0"`, a formula in the unit of the dose state
        in the model), so repeated doses keep the rest of the earlier ones. The
        first dose is usually a change before the initialization
        (`preinit_changes`), which may set the state. The models are loaded
        before `simulations()`, where the doses are created.

        Args:
            times: the times of the doses, in the time unit of the simulation
                or quantities.
            doses: the dose state and the dose, e.g.
                `{"PODOSE_hctz": Q(10, "mg")}`.
            model: the key of the model whose units the dose states have.

        Returns:
            The change of the dose states at the times.

        Raises:
            ValueError: if the model is not loaded.
        """
        loaded = self._models.get(model)
        uinfo = getattr(loaded, "uinfo", None)
        if uinfo is None:
            raise ValueError(
                f"The model '{model}' of the experiment '{self.sid}' is not loaded: "
                "the doses are created in `simulations()` or later."
            )
        values: dict[str, str] = {}
        for sid, dose in doses.items():
            unit = uinfo.quantity(1.0, uinfo[sid]).units
            magnitude = float(dose.to(unit).magnitude)
            values[sid] = f"{sid} + {magnitude!r}"
        return Change(times, values)

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

        The units of the data are renamed (`data_units`, `rename_units`); the
        geometric means get the unit of their row when the table is read as a
        `DataSet` (`DataSet.from_df`).

        Args:
            fig_id: the id of the figure or table in the study, e.g. `Fig1`.
            study: the study, the id of the experiment by default.
            na_values: additional strings of missing values, for every column or
                by column, e.g. `{"time": ["NR"]}`, see `pandas.read_csv`.

        Raises:
            ValueError: if the experiment has no `data_path` or a definition of
                `data_units` is refused.
        """
        names = self.data_unit_names()
        df = load_pkdb_dataframe(
            f"{self.sid if study is None else study}_{fig_id}",
            data_path=self._data_paths(),
            na_values=na_values,
        )
        return rename_units(df, names)

    def load_dataframes_by_substance(
        self, fig_id: str, study: str | None = None
    ) -> dict[str, pd.DataFrame]:
        """Load the digitized data of a figure or table by substance.

        The units of the data are renamed (`data_units`, `rename_units`). The rows without a
        substance (a bile flow, a blood pressure) are in no frame, they are in
        the table of `load_dataframe`.

        Args:
            fig_id: the id of the figure or table in the study, e.g. `Fig1`.
            study: the study, the id of the experiment by default.

        Raises:
            ValueError: if the experiment has no `data_path` or a definition of
                `data_units` is refused.
        """
        names = self.data_unit_names()
        df = load_pkdb_dataframe(
            f"{self.sid if study is None else study}_{fig_id}",
            data_path=self._data_paths(),
            na_values=["na"],
        )
        return {
            str(substance): rename_units(frame.copy(), names)
            for substance, frame in df.groupby("substance", sort=False)
        }

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
