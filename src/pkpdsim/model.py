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

`FigureSettings` is the `sbmlsim.plot.FigureSettings` of sbmlsim with the
formats of the figures and the resolution of the static images; the module
therefore imports sbmlsim's plotting.
"""

from __future__ import annotations

import importlib
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from sbmlsim.plot import FigureSettings as SbmlsimFigureSettings

from pkpdsim.release.settings import ReleaseSettings

if TYPE_CHECKING:
    from sbmlsim.experiment import SimulationExperiment
    from sbmlsim.fit.cli import FitDefinition
    from sbmlsim.model.tolerances import AbsoluteTolerance

    from pkpdsim.uncertainty import FitRecords, UncertaintySettings

#: the actions of the `run` command of a model, see `pkpdsim.cli`
RUN_ACTIONS: tuple[str, ...] = ("simulate", "list_experiments", "factory", "all")


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


def normalize_formats(formats: Sequence[str]) -> tuple[str, ...]:
    """Get figure formats as a lower-case tuple.

    Raises:
        TypeError: if the formats are a single string, e.g. `"svg"` instead of
            `("svg",)`; a string is a sequence of characters.
    """
    if isinstance(formats, str):
        raise TypeError(
            f"The figure formats must be a sequence of strings, not {formats!r}; "
            f"use ({formats!r},)."
        )
    return tuple(f.lower() for f in formats)


@dataclass(frozen=True, kw_only=True)
class SimulatorSettings:
    """The integrator settings of the simulation experiments of a model.

    All fields are keyword-only.

    Attributes:
        absolute_tolerance: the absolute tolerance, one value for every state or
            an `AbsoluteTolerance` of sbmlsim with tolerances per kind and state.
        relative_tolerance: the relative tolerance.
        n_workers: the processes of a scan: `1` runs serially, a number runs a
            pool of that many processes, `None` lets sbmlsim choose (a pool of
            every CPU for a scan of `sbmlsim.parallel.POOL_THRESHOLD` points or
            more, serially otherwise). A pool needs the guard
            `if __name__ == "__main__":` in a script which runs the model.
    """

    absolute_tolerance: float | AbsoluteTolerance = 1e-10
    relative_tolerance: float = 1e-10
    n_workers: int | None = 1

    def __post_init__(self) -> None:
        """Check the number of processes.

        Raises:
            ValueError: if `n_workers` is neither a positive int nor `None`.
        """
        n = self.n_workers
        if n is not None and (isinstance(n, bool) or not isinstance(n, int) or n < 1):
            raise ValueError(
                f"n_workers is a positive number of processes or None, not {n!r}."
            )


@dataclass(frozen=True, kw_only=True)
class FigureSettings(SbmlsimFigureSettings):
    """The figures of the simulation experiments of a model.

    The settings of sbmlsim's figures (`sbmlsim.plot.FigureSettings`: the
    resolution, the sizes of the panels, the fonts and the legends, with the
    defaults of sbmlsim) and what a model adds: the formats of the figures and
    the resolution of the static images. A run hands the settings of sbmlsim to
    every figure (`sbmlsim_settings`); an experiment (`figure_settings`) and a
    figure (`Figure(settings=...)`) override them. `savefig_dpi` is the
    `savefig.dpi` of matplotlib while the figures of a run are drawn
    (`rc_params`).

    The fields of sbmlsim can be given by position, as in sbmlsim; `formats`
    and `savefig_dpi` only by keyword.

    Attributes:
        formats: the formats of the figures, static images drawn by matplotlib
            (`svg`, `png`) and interactive pages drawn by plotly (`html`); SVG
            only by default, a run draws other formats on request
            (`--formats svg,png`).
        savefig_dpi: the resolution of the static images, `None` keeps the one
            of matplotlib. An interactive page has 72 pixels per inch.
    """

    formats: Sequence[str] = ("svg",)
    savefig_dpi: float | None = None

    def __post_init__(self) -> None:
        """Check that there is a format and normalize the formats.

        Raises:
            TypeError: if the formats are a single string.
            ValueError: if there is no format.
        """
        formats = normalize_formats(self.formats)
        if not formats:
            raise ValueError("A model needs at least one figure format.")
        object.__setattr__(self, "formats", formats)

    def sbmlsim_settings(self) -> SbmlsimFigureSettings:
        """Get the settings of sbmlsim, without the formats and the resolution.

        A plain `sbmlsim.plot.FigureSettings`, which sbmlsim compares with its
        defaults (a run with the defaults sets nothing).
        """
        return SbmlsimFigureSettings(
            **{name: getattr(self, name) for name in SbmlsimFigureSettings.names()}
        )

    def rc_params(self) -> dict[str, float]:
        """Get the settings of matplotlib for a run, the resolution of the images."""
        return {} if self.savefig_dpi is None else {"savefig.dpi": self.savefig_dpi}


@dataclass(frozen=True, kw_only=True)
class ExtraAction:
    """An action a model adds to its `run` command, e.g. analyses of its own.

    All fields are keyword-only.

    Attributes:
        name: the value of `--action`.
        help: what the action does, shown by `--help`.
        function: reference to the function which runs the action; it is called
            with the model, the parsed arguments and the directory of the results.
            The arguments have the settings of the simulator of the run as
            `simulator` (`SimulatorSettings`, those of the model with the
            processes of `--workers`), which an action that simulates passes on,
            e.g. `model.run_experiments(..., simulator=options.simulator)` or
            `pkpdsim.run.create_simulator(options.simulator)`.
        in_all: run the action with `--action all`, after the simulations.
    """

    name: str
    help: str
    function: str
    in_all: bool = False


@dataclass(frozen=True, kw_only=True)
class CheckSettings:
    """The known exceptions of the tests of a model, each with its reason.

    All fields are keyword-only. The mappings are copied, so the settings do not
    change with the dictionaries they were created from, and the settings are
    hashable.

    Attributes:
        not_registered: experiment classes which are defined but not in the
            group `all` of the registry on purpose, by class name.
        not_importable: modules below the package which do not import in the
            test environment (an optional dependency), by their name relative to
            the package, e.g. `experiments.zonation.plot_pyvista`.
        one_at_a_time: why the tests run the registered experiments one at a
            time (the memory a run of all of them needs), `None` to run them
            together.
    """

    not_registered: Mapping[str, str] = field(default_factory=dict)
    not_importable: Mapping[str, str] = field(default_factory=dict)
    one_at_a_time: str | None = None

    def __post_init__(self) -> None:
        """Copy the mappings and check the reason to run one at a time.

        Raises:
            ValueError: if one_at_a_time is an empty reason.
        """
        object.__setattr__(self, "not_registered", dict(self.not_registered))
        object.__setattr__(self, "not_importable", dict(self.not_importable))
        if self.one_at_a_time is not None and not self.one_at_a_time.strip():
            raise ValueError(
                "one_at_a_time needs a reason, e.g. the memory a run of all "
                "experiments needs; None runs them together."
            )

    def __hash__(self) -> int:
        """Hash the settings by their content."""
        return hash(
            (
                tuple(sorted(self.not_registered.items())),
                tuple(sorted(self.not_importable.items())),
                self.one_at_a_time,
            )
        )


@dataclass(frozen=True, kw_only=True)
class PKPDModel:
    """A PBPK/PD model: what the commands, the verification and the tests need.

    All fields are keyword-only, so their order is not API. The sequences
    (`data_paths`, `extra_actions` and the formats of the figures) are stored as
    tuples, so a list can be passed and the model is hashable, unless the
    absolute tolerance of the simulator is an `AbsoluteTolerance`, which is not
    hashable.

    Attributes:
        name: the short name of the model.
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
        checks: the known exceptions of the tests of the model, see
            `pkpdsim.testing`.
        fits: the ids of the adopted fits, whose records in `fits_dir` give the
            fitted values and their uncertainty, see `pkpdsim.uncertainty`.
        fits_path: the folder of the adopted fits, `fits` next to the output
            directory of the factory by default.
        release: the release metadata, see `pkpdsim.release.settings`; a model
            without it cannot be released.
    """

    name: str
    package: str
    title: str
    path: Path
    model_path: Path
    factory: str
    experiments: str | None = None
    fit_definitions: str | None = None
    data_paths: Sequence[Path] = ()
    results_path: Path | None = None
    factory_path: Path | None = None
    archive_path: Path | None = None
    simulator: SimulatorSettings = field(default_factory=SimulatorSettings)
    figures: FigureSettings = field(default_factory=FigureSettings)
    reduced_selections: bool = True
    save_results: bool = False
    fit_identifiability: bool = False
    extra_actions: Sequence[ExtraAction] = ()
    extra_arguments: str | None = None
    checks: CheckSettings = field(default_factory=CheckSettings)
    fits: Sequence[str] = ()
    fits_path: Path | None = None
    release: ReleaseSettings | None = None

    def __post_init__(self) -> None:
        """Normalize the sequences, check the references and the extra actions."""
        object.__setattr__(self, "data_paths", tuple(self.data_paths))
        object.__setattr__(self, "fits", tuple(self.fits))
        repeated_fits = sorted({f for f in self.fits if self.fits.count(f) > 1})
        if repeated_fits:
            raise ValueError(
                f"The adopted fits of '{self.name}' repeat a fit: {repeated_fits}."
            )
        object.__setattr__(self, "extra_actions", tuple(self.extra_actions))
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
    def uncertainty_dir(self) -> Path:
        """The directory of the results of the experiments over the draws of the fits."""
        return self.results_dir / "simulation_uncertainty"

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

    @property
    def fits_dir(self) -> Path:
        """The folder of the adopted fits, next to the output of the factory."""
        return self.fits_path or self.factory_dir.parent / "fits"

    def fit_records(self) -> FitRecords:
        """Get the records of the adopted fits, see `pkpdsim.uncertainty.FitRecords`."""
        from pkpdsim.uncertainty import FitRecords

        return FitRecords(self.fits_dir, self.fits)

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
        simulator: SimulatorSettings | None = None,
        uncertainty: UncertaintySettings | None = None,
    ) -> Path:
        """Run simulation experiments and create their report.

        See `pkpdsim.run.run_experiments`; a run with a failed experiment or
        figure raises `sbmlsim.experiment.ExperimentRunError` after its report
        is written.
        """
        from pkpdsim.run import run_experiments

        return run_experiments(
            self,
            experiment_classes,
            output_dir,
            save_results=save_results,
            show_report=show_report,
            figure_formats=figure_formats,
            simulator=simulator,
            uncertainty=uncertainty,
        )

    def run_simulation_experiments(
        self,
        experiment_classes: Sequence[type[SimulationExperiment]],
        output_dir: str | Path,
        *,
        show_report: bool = False,
        figure_formats: Sequence[str] | None = None,
        simulator: SimulatorSettings | None = None,
        uncertainty: UncertaintySettings | None = None,
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
            simulator=simulator,
            uncertainty=uncertainty,
        )
