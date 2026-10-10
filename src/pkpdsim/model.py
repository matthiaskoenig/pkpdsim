"""The description of a PBPK/PD model, which the commands of pkpdsim work with.

A model package defines one `PKPDModel` in its `__init__.py`, next to its path
constants:

```python
# not executed
MODEL = PKPDModel(
    name="hctz",
    package=__name__,
    title="HCTZ PBPK/PD model",
    path=HCTZ_PATH,
    model_path=MODEL_PATH,
    data_paths=tuple(DATA_PATHS),
    factory=".models.factory:create_hctz_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
)
```

The factory, the registry of the experiments and the fit definitions are
references `"module:attribute"`, resolved when they are used: the experiments
import the `__init__.py` of their model, so the model object cannot import them,
and importing a model package does not import its experiments. A module which
starts with a dot is relative to `package`.
"""

from __future__ import annotations

import importlib
from collections.abc import Callable, Generator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment
    from sbmlsim.fit.cli import FitDefinition
    from sbmlsim.model.tolerances import AbsoluteTolerance

#: the actions of the `run` command of a model, see `pkpdsim.cli`
RUN_ACTIONS: tuple[str, ...] = ("simulate", "list_experiments", "factory", "all")

#: the class attributes of `sbmlsim.plot.Figure` which `FigureSettings` sets
FIGURE_ATTRIBUTES: tuple[str, ...] = (
    "fig_dpi",
    "fig_titlesize",
    "axes_labelsize",
    "xtick_labelsize",
    "ytick_labelsize",
    "legend_fontsize",
    "legend_position",
    "legend_loc",
    "panel_width",
    "panel_height",
)


def check_reference(reference: str) -> None:
    """Check that a reference has the form `"module:attribute"`.

    The module may start with dots (relative to a package), the attribute may be
    dotted (an attribute of an attribute).

    Args:
        reference: the reference to check.

    Raises:
        ValueError: if the reference does not have the form.
    """
    module, separator, attribute = reference.partition(":")
    names = [*module.lstrip(".").split("."), *attribute.split(".")]
    if not separator or not all(name.isidentifier() for name in names):
        raise ValueError(f"A reference is 'module:attribute', not '{reference}'.")


def resolve(reference: str, package: str | None = None) -> Any:  # noqa: ANN401
    """Resolve a reference `"module:attribute"` to the object it names.

    Args:
        reference: the reference; a module which starts with a dot is relative to
            `package`.
        package: the package of relative references.

    Returns:
        The object.

    Raises:
        ValueError: if the reference is malformed or relative without a package.
        AttributeError: if the module has no such attribute; the message names
            the reference.
    """
    check_reference(reference)
    module_name, _, attribute = reference.partition(":")
    if module_name.startswith(".") and package is None:
        raise ValueError(f"The relative reference '{reference}' needs a package.")
    obj: Any = importlib.import_module(module_name, package)
    for name in attribute.split("."):
        try:
            obj = getattr(obj, name)
        except AttributeError as err:
            raise AttributeError(f"'{reference}' cannot be resolved: {err}") from err
    return obj


@dataclass(frozen=True)
class SimulatorSettings:
    """The integrator settings of the simulation experiments of a model.

    Attributes:
        absolute_tolerance: the absolute tolerance, one value for every state or
            an `AbsoluteTolerance` of sbmlsim with tolerances per kind and state.
        relative_tolerance: the relative tolerance.
    """

    absolute_tolerance: float | AbsoluteTolerance = 1e-10
    relative_tolerance: float = 1e-10


@dataclass(frozen=True)
class FigureSettings:
    """The figures of the simulation experiments of a model.

    The sizes and positions are class attributes of `sbmlsim.plot.Figure` and
    `savefig_dpi` is the `savefig.dpi` of matplotlib: both hold for every figure
    of the process. `applied` sets them for a run and restores them after it, so
    no setting leaks into the run of another model in the same process. `None`
    keeps the value of sbmlsim.

    Attributes:
        formats: the formats of the figures, static images drawn by matplotlib
            (`svg`, `png`) and interactive pages drawn by plotly (`html`).
        savefig_dpi: the resolution of the static images only. `fig_dpi` also sets
            the size of the interactive pages in pixels.
        fig_dpi: the resolution of the figures of sbmlsim.
        fig_titlesize: the font size of the title of a figure.
        axes_labelsize: the font size of the labels of the axes.
        xtick_labelsize: the font size of the ticks of the x axis.
        ytick_labelsize: the font size of the ticks of the y axis.
        legend_fontsize: the font size of the legends.
        legend_position: `inside` or `outside` of the plot.
        legend_loc: the location of the legend in the plot, e.g. `best`.
        panel_width: the width of a panel.
        panel_height: the height of a panel.
    """

    formats: tuple[str, ...] = ("svg", "png")
    savefig_dpi: float | None = None
    fig_dpi: float | None = None
    fig_titlesize: float | None = None
    axes_labelsize: float | None = None
    xtick_labelsize: float | None = None
    ytick_labelsize: float | None = None
    legend_fontsize: float | None = None
    legend_position: str | None = None
    legend_loc: str | None = None
    panel_width: float | None = None
    panel_height: float | None = None

    def __post_init__(self) -> None:
        """Check that there is a format."""
        if not self.formats:
            raise ValueError("A model needs at least one figure format.")

    @contextmanager
    def applied(self) -> Generator[None]:
        """Apply the settings while the context runs, restore them after it."""
        import matplotlib
        from sbmlsim.plot import Figure

        changed = {
            name: getattr(self, name)
            for name in FIGURE_ATTRIBUTES
            if getattr(self, name) is not None
        }
        previous = {name: getattr(Figure, name) for name in changed}
        rc: dict[Any, Any] = {}
        if self.savefig_dpi is not None:
            rc["savefig.dpi"] = self.savefig_dpi
        try:
            for name, value in changed.items():
                setattr(Figure, name, value)
            with matplotlib.rc_context(rc):
                yield
        finally:
            for name, value in previous.items():
                setattr(Figure, name, value)


@dataclass(frozen=True)
class ExtraAction:
    """An action a model adds to its `run` command, e.g. analyses of its own.

    Attributes:
        name: the value of `--action`.
        help: what the action does, shown by `--help`.
        function: reference to the function which runs the action; it is called
            with the model, the parsed arguments and the directory of the results.
        in_all: run the action with `--action all`, after the simulations.
    """

    name: str
    help: str
    function: str
    in_all: bool = False


@dataclass(frozen=True)
class PKPDModel:
    """A PBPK/PD model: what the commands, the verification and the tests need.

    Attributes:
        name: the short name, the suffix of the commands (`run_<name>`).
        package: the dotted name of the model package, `__name__` in its
            `__init__.py`; relative references are relative to it.
        title: the title of the model in the output of the commands.
        path: the directory of the model package, the base path of its
            experiments.
        model_path: the SBML model the experiments simulate.
        factory: reference to the function which creates the models; it is
            called with the output directory.
        experiments: reference to the registry of the experiments, groups of
            experiment classes by name with the group `all`.
        fit_definitions: reference to the fit problems by name
            (`sbmlsim.fit.cli.FitDefinition`).
        data_paths: the directories of the study data.
        results_path: the directory of the results, `<path>/results` by default.
        factory_path: the output directory of the factory, the directory of
            `model_path` by default.
        archive_path: the COMBINE archive the factory creates, linked by the
            `factory` action.
        simulator: the integrator settings of the experiments.
        figures: the formats and settings of the figures.
        reduced_selections: simulate only the selections the experiments need.
        save_results: store the results of the simulations by default.
        fit_identifiability: the fit command analyses the identifiability of the
            fitted parameters by default.
        extra_actions: actions the model adds to its `run` command.
        extra_arguments: reference to a function which adds the options of the
            extra actions to the parser of the `run` command.
    """

    name: str
    package: str
    title: str
    path: Path
    model_path: Path
    factory: str
    experiments: str | None = None
    fit_definitions: str | None = None
    data_paths: tuple[Path, ...] = ()
    results_path: Path | None = None
    factory_path: Path | None = None
    archive_path: Path | None = None
    simulator: SimulatorSettings = field(default_factory=SimulatorSettings)
    figures: FigureSettings = field(default_factory=FigureSettings)
    reduced_selections: bool = True
    save_results: bool = False
    fit_identifiability: bool = False
    extra_actions: tuple[ExtraAction, ...] = ()
    extra_arguments: str | None = None

    def __post_init__(self) -> None:
        """Check the references and the names of the extra actions."""
        references = [
            self.factory,
            self.experiments,
            self.fit_definitions,
            self.extra_arguments,
            *(action.function for action in self.extra_actions),
        ]
        for reference in references:
            if reference is not None:
                check_reference(reference)
        names = [action.name for action in self.extra_actions]
        repeated = sorted({n for n in names if n in RUN_ACTIONS or names.count(n) > 1})
        if repeated:
            raise ValueError(
                f"The extra actions of '{self.name}' repeat an action: {repeated}."
            )

    @property
    def results_dir(self) -> Path:
        """The directory of the results."""
        return self.results_path or self.path / "results"

    @property
    def simulation_dir(self) -> Path:
        """The directory of the results of the simulation experiments."""
        return self.results_dir / "simulation"

    @property
    def fit_dir(self) -> Path:
        """The directory of the results of the fits."""
        return self.results_dir / "fit"

    @property
    def report_dir(self) -> Path:
        """The directory of the reports of stored parameters."""
        return self.results_dir / "report"

    @property
    def identifiability_dir(self) -> Path:
        """The directory of the identifiability analyses of stored parameters."""
        return self.results_dir / "identifiability"

    @property
    def factory_dir(self) -> Path:
        """The output directory of the factory."""
        return self.factory_path or self.model_path.parent

    def load_factory(self) -> Callable[[Path], object]:
        """Get the function which creates the models."""
        return resolve(self.factory, self.package)

    def load_experiments(
        self,
    ) -> Mapping[str, list[type[SimulationExperiment]]]:
        """Get the registry of the experiments."""
        if self.experiments is None:
            raise ValueError(f"The model '{self.name}' has no registry of experiments.")
        return resolve(self.experiments, self.package)

    def load_fit_definitions(self) -> Mapping[str, FitDefinition]:
        """Get the fit problems by name."""
        if self.fit_definitions is None:
            raise ValueError(f"The model '{self.name}' has no fit definitions.")
        return resolve(self.fit_definitions, self.package)

    def create_models(self, output_dir: Path | None = None) -> Path:
        """Create the models with the factory.

        Args:
            output_dir: the output directory, `factory_dir` by default.

        Returns:
            The output directory.
        """
        directory = output_dir or self.factory_dir
        self.load_factory()(directory)
        return directory

    def run_experiments(
        self,
        experiment_classes: type[SimulationExperiment]
        | Sequence[type[SimulationExperiment]],
        output_dir: str | Path,
        *,
        save_results: bool | None = None,
        show_report: bool = False,
        figure_formats: Sequence[str] | None = None,
    ) -> Path:
        """Run simulation experiments and create their report.

        See `pkpdsim.run.run_experiments`.
        """
        from pkpdsim.run import run_experiments

        return run_experiments(
            self,
            experiment_classes,
            output_dir,
            save_results=save_results,
            show_report=show_report,
            figure_formats=figure_formats,
        )

    def run_simulation_experiments(
        self,
        experiment_classes: Sequence[type[SimulationExperiment]],
        output_dir: str | Path,
        *,
        show_report: bool = False,
        figure_formats: Sequence[str] | None = None,
    ) -> Path:
        """Run simulation experiments and collect their figures.

        See `pkpdsim.run.run_simulation_experiments`.
        """
        from pkpdsim.run import run_simulation_experiments

        return run_simulation_experiments(
            self,
            experiment_classes,
            output_dir,
            show_report=show_report,
            figure_formats=figure_formats,
        )

    def run_main(self, args: Sequence[str] | None = None) -> None:
        """The command `run_<name>`, see `pkpdsim.cli.run_command`."""
        from pkpdsim.cli import run_command

        run_command(self, args)

    def fit_main(self, args: Sequence[str] | None = None) -> None:
        """The command `fit_<name>`, see `pkpdsim.cli.fit_command`.

        Returns nothing: a console script exits with the return value of the
        function it calls, and 1 for anything which is not `None` or an int.
        """
        from pkpdsim.cli import fit_command

        fit_command(self, args)

    def report_main(self, args: Sequence[str] | None = None) -> None:
        """The command `report_<name>`, see `pkpdsim.cli.report_command`."""
        from pkpdsim.cli import report_command

        report_command(self, args)

    def identifiability_main(self, args: Sequence[str] | None = None) -> None:
        """The command `identifiability_<name>`, see `pkpdsim.cli`."""
        from pkpdsim.cli import identifiability_command

        identifiability_command(self, args)
