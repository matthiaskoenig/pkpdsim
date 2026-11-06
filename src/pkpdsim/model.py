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

A `PKPDModel` is a `sbmlsim.project.Project`: the factory, the experiments, the
fit definitions and the credibility declaration are references
`"module:attribute"`, resolved when they are used, the settings of the
simulator are the `sbmlsim.project.SimulatorSettings`, the known exceptions of
the tests the `sbmlsim.project.CheckSettings` and the release metadata the
`sbmlsim.release.settings.ReleaseSettings` of sbmlsim. The keywords of the
model packages (`path`, `model_path`, `female_model_path`, `data_paths`,
`fit_definitions`, `fits`, `fits_path`, `factory_path`, `figures`) are mapped
onto the fields of the project when a model is created, see `PKPDModel`.

A model whose studies include women has a female build next to the male one
(`female_model_path`, the factory builds both from one definition, see
`pkpdsim.pbpk.create_archive`): the builds are the variants of the project by
the value of their `pkpdsim.sex.Sex`, the male build first, the model of every
caller which does not ask for a sex, and an experiment chooses the build of the
sex of its subjects (`pkpdsim.experiment.PKPDSimulationExperiment.sex`).

`FigureSettings` is the `sbmlsim.plot.FigureSettings` of sbmlsim with the
formats of the figures and the resolution of the static images, which a model
splits into the `figure_settings` and the `figures` of its project; the module
therefore imports sbmlsim's plotting.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Self, override

from pydantic import model_validator
from sbmlsim.plot import FigureSettings as SbmlsimFigureSettings
from sbmlsim.project import (
    FigureOutputs,
    Project,
    SimulatorSettings,
    check_reference,
    normalize_formats,
)

from pkpdsim.sex import Sex

if TYPE_CHECKING:
    import roadrunner
    from sbmlsim.experiment import SimulationExperiment

    from pkpdsim.credibility import Credibility
    from pkpdsim.uncertainty import UncertaintySettings

#: the actions of the `run` command of a model, see `pkpdsim.cli`
RUN_ACTIONS: tuple[str, ...] = ("simulate", "list_experiments", "factory", "all")

#: the processes of a run of the experiments of a model by default: serially,
#: so that a script which runs a model needs no guard `if __name__ == "__main__":`
RUN_WORKERS = 1

#: the settings of the simulator of a model which does not give its own
DEFAULT_SIMULATOR = SimulatorSettings(
    absolute_tolerance=1e-10, relative_tolerance=1e-10
)

#: the keywords of a model package which are fields of the project by another name
_RENAMED: dict[str, str] = {
    "path": "root",
    "data_paths": "data_path",
    "fits_path": "fits_dir",
    "factory_path": "factory_dir",
}


@dataclass(frozen=True, kw_only=True)
class FigureSettings(SbmlsimFigureSettings):
    """The figures of the simulation experiments of a model.

    The settings of sbmlsim's figures (`sbmlsim.plot.FigureSettings`: the
    resolution, the sizes of the panels, the fonts and the legends, with the
    defaults of sbmlsim) and what a model adds: the formats of the figures and
    the resolution of the static images. A model splits them into the settings
    of its figures (`sbmlsim_settings`, the `figure_settings` of the project)
    and the files of its figures (`outputs`, the `figures` of the project).

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

    def outputs(self) -> FigureOutputs:
        """Get the files of the figures: the formats and the resolution."""
        return FigureOutputs(formats=tuple(self.formats), savefig_dpi=self.savefig_dpi)


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
            `simulator` (`sbmlsim.project.SimulatorSettings`, those of the
            model) and the processes of `--workers` as `n_workers`, which an
            action that simulates passes on, e.g.
            `model.run_experiments(..., simulator=options.simulator,
            n_workers=options.n_workers)` or
            `pkpdsim.run.create_simulator(options.simulator, options.n_workers)`.
        in_all: run the action with `--action all`, after the simulations.
    """

    name: str
    help: str
    function: str
    in_all: bool = False


def _figures(data: dict[str, Any]) -> None:
    """Split the figure settings of a model into those of its project.

    Args:
        data: the keywords of the model, changed in place.

    Raises:
        ValueError: if the settings of sbmlsim are given twice.
    """
    figures = data.get("figures")
    if not isinstance(figures, SbmlsimFigureSettings):
        return
    if data.get("figure_settings") is not None:
        raise ValueError(
            "The settings of the figures are given twice, as figures and as "
            "figure_settings."
        )
    if isinstance(figures, FigureSettings):
        data["figure_settings"] = figures.sbmlsim_settings()
        data["figures"] = figures.outputs()
    else:
        data["figure_settings"] = figures
        del data["figures"]


def _builds(data: dict[str, Any]) -> None:
    """Map the builds of a model onto the variants of its project.

    Args:
        data: the keywords of the model, changed in place.

    Raises:
        ValueError: if the model has no build, the variants are given as well,
            or the female build is the file of the male build.
    """
    if "model_path" not in data and "female_model_path" not in data:
        if not data.get("variants"):
            raise ValueError("A model needs its build, model_path.")
        return
    if "variants" in data:
        raise ValueError(
            "The builds of a model are given twice, as model_path or "
            "female_model_path and as variants."
        )
    male = data.pop("model_path", None)
    female = data.pop("female_model_path", None)
    if male is None:
        raise ValueError("A model needs its build, model_path.")
    if female is not None and Path(female) == Path(male):
        raise ValueError(
            f"The female build of '{data.get('name')}' is the file of the male "
            f"build {male}: a build per sex is a file per sex."
        )
    variants = {Sex.MALE.value: male}
    if female is not None:
        variants[Sex.FEMALE.value] = female
    data["variants"] = variants


def _fits(data: dict[str, Any]) -> None:
    """Map the fit definitions and the adopted fits of a model onto its project.

    `fit_definitions` are the `fits` of the project; `fits` given as a
    sequence of fit ids are the adopted fits, `adopted_ids`.

    Args:
        data: the keywords of the model, changed in place.

    Raises:
        ValueError: if the adopted fits or the fit definitions are given twice.
    """
    fits = data.get("fits")
    if isinstance(fits, Sequence) and not isinstance(fits, str):
        if "adopted_ids" in data:
            raise ValueError(
                "The adopted fits are given twice, as fits and as adopted_ids."
            )
        data["adopted_ids"] = data.pop("fits")
    if "fit_definitions" in data:
        if "fits" in data:
            raise ValueError(
                "The fit definitions are given twice, as fit_definitions and as fits."
            )
        definitions = data.pop("fit_definitions")
        if definitions is not None:
            data["fits"] = definitions


def _directories(data: dict[str, Any]) -> None:
    """Set the directories of the factory and of the adopted fits of a model.

    The directory of the factory is the one of the male build, the folder of
    the adopted fits `fits` next to it, unless they are given; relative paths
    are relative to the root of the project, as the fields resolve them.

    Args:
        data: the keywords of the model, changed in place.
    """
    variants = data.get("variants")
    if data.get("factory_dir") is None and isinstance(variants, Mapping) and variants:
        data["factory_dir"] = Path(next(iter(variants.values()))).parent
    factory_dir = data.get("factory_dir")
    if data.get("fits_dir") is None and factory_dir is not None:
        data["fits_dir"] = Path(factory_dir).parent / "fits"


class PKPDModel(Project):
    """A PBPK/PD model: what the commands, the verification and the tests need.

    A `sbmlsim.project.Project` with what pkpdsim adds: the COMBINE archive,
    the builds per sex, the adopted fits, the directories of the results below
    the package and the actions of the `run` command. The
    keywords of the model packages are mapped onto the fields of the project
    when a model is created from keywords:

    - `path`: `root`, the directory of the model package, the base path of
      its experiments (the property `path`).
    - `model_path`, `female_model_path`: the `variants` `male` and `female`,
      the male build first: the model of every experiment and every caller
      which does not ask for a sex. The builds have the same ids, so a fit
      applies one parameter set to both, see `model_paths`.
    - `data_paths`: `data_path`, the directories of the study data (the
      property `data_paths`).
    - `fit_definitions`: `fits`, the reference to the fit problems by name
      (`sbmlsim.fit.cli.FitDefinition`).
    - `fits` as a sequence of fit ids: `adopted_ids`.
    - `fits_path`: `fits_dir`, `fits` next to the directory of the factory
      by default.
    - `factory_path`: `factory_dir`, the directory of the male build by
      default.
    - `figures` as a `FigureSettings`: `figure_settings` and `figures`
      (`FigureSettings.sbmlsim_settings`, `FigureSettings.outputs`).

    A keyword and its field are not given together. A model is a frozen
    pydantic model: `replace` gives a changed copy, by the names of the
    fields.

    Attributes:
        archive_path: the COMBINE archive the factory creates, linked by the
            `factory` action.
        results_path: the directory of the results, `<root>/results` by
            default.
        reduced_selections: simulate only the selections the experiments need.
        save_results: store the results of the simulations by default.
        fit_identifiability: the fit command analyses the identifiability of the
            fitted parameters by default.
        extra_actions: actions the model adds to its `run` command.
        extra_arguments: reference to a function which adds the options of the
            extra actions to the parser of the `run` command.
        adopted_ids: the ids of the adopted fits, whose records in `fits_dir`
            give the fitted values and their uncertainty, see
            `sbmlsim.fit.records`; one fit per fit subset.
    """

    # the fields of the project a model always has, narrowed: a model is frozen
    package: str  # ty: ignore[invalid-mutable-override]
    factory_dir: Path  # ty: ignore[invalid-mutable-override]
    fits_dir: Path  # ty: ignore[invalid-mutable-override]
    archive_path: Path | None = None
    results_path: Path | None = None
    reduced_selections: bool = True
    save_results: bool = False
    fit_identifiability: bool = False
    extra_actions: tuple[ExtraAction, ...] = ()
    extra_arguments: str | None = None
    adopted_ids: tuple[str, ...] = ()
    simulator: SimulatorSettings = DEFAULT_SIMULATOR

    def __init__(self, **data: Any) -> None:  # noqa: ANN401
        """Create a model from the keywords of a model package or the fields.

        The keywords of the model packages are mapped onto the fields of the
        project, see the class.

        Raises:
            pydantic.ValidationError: for invalid keywords or fields.
        """
        super().__init__(**data)

    @override
    def replace(self, **changes: Any) -> Self:
        """Get a copy of the model with other fields, validated again.

        See `sbmlsim.project.Project.replace`. A `FigureSettings` as `figures`
        replaces the settings of the figures as well, `figure_settings`. A
        copy with other builds (`variants`) or another `root` derives the
        directory of the factory and the folder of the adopted fits again,
        unless they are given or were not the derived ones: the directory of
        the male build and `fits` next to the directory of the factory.

        Raises:
            pydantic.ValidationError: for invalid fields.
        """
        figures = changes.get("figures")
        if isinstance(figures, SbmlsimFigureSettings):
            changes.setdefault("figure_settings", None)
        if "variants" in changes or "root" in changes:
            derived_factory = self.factory_dir == self.model_path.parent
            if "factory_dir" not in changes and derived_factory:
                changes["factory_dir"] = None
            derived_fits = self.fits_dir == self.factory_dir.parent / "fits"
            if "fits_dir" not in changes and derived_fits:
                changes["fits_dir"] = None
        return super().replace(**changes)

    @model_validator(mode="before")
    @classmethod
    def _keywords_of_a_package(cls, data: Any) -> Any:  # noqa: ANN401
        """Map the keywords of a model package onto the fields of the project.

        Only keywords are mapped, a model given as an object is itself.

        Raises:
            ValueError: if a keyword and its field are given together, see the
                class.
        """
        if not isinstance(data, Mapping):
            return data
        data = dict(data)
        for keyword, name in _RENAMED.items():
            if keyword not in data:
                continue
            if name in data:
                raise ValueError(
                    f"'{keyword}' is the field '{name}' of the project, give one "
                    "of them."
                )
            data[name] = data.pop(keyword)
        _builds(data)
        _fits(data)
        _figures(data)
        _directories(data)
        return data

    @model_validator(mode="after")
    def _check_model(self) -> Self:
        """Check the builds, the adopted fits, the actions and the known failures.

        Raises:
            ValueError: for a model without its build or its factory, a variant
                which is no sex or a female build first, a repeated adopted
                fit, an extra action which repeats an action, a malformed
                reference, or a known failure of the battery which is no check
                of the battery.
        """
        if not self.variants:
            raise ValueError(f"The model '{self.name}' needs its build, model_path.")
        names = list(self.variants)
        if names[0] != Sex.MALE.value or not set(names) <= {s.value for s in Sex}:
            raise ValueError(
                f"The builds of the model '{self.name}' are the variants {names}: "
                "the variants of a model are its builds per sex, the male build "
                "first."
            )
        if self.factory is None:
            raise ValueError(f"The model '{self.name}' needs its factory.")
        repeated_fits = sorted(
            {f for f in self.adopted_ids if self.adopted_ids.count(f) > 1}
        )
        if repeated_fits:
            raise ValueError(
                f"The adopted fits of '{self.name}' repeat a fit: {repeated_fits}."
            )
        for reference in [
            self.extra_arguments,
            *(action.function for action in self.extra_actions),
        ]:
            if reference is not None:
                check_reference(reference)
        actions = [action.name for action in self.extra_actions]
        repeated = sorted(
            {n for n in actions if n in RUN_ACTIONS or actions.count(n) > 1}
        )
        if repeated:
            raise ValueError(
                f"The extra actions of '{self.name}' repeat an action: {repeated}."
            )
        if self.checks.battery:
            from pkpdsim.credibility.verification import BATTERY_CHECKS

            checks = [c.value for c in BATTERY_CHECKS]
            unknown = sorted(set(self.checks.battery) - set(checks))
            if unknown:
                raise ValueError(
                    f"The failures of the battery {unknown} are no checks of the "
                    f"battery, one of {', '.join(checks)}."
                )
        return self

    # ------------------------------------------------------------------------
    # paths
    # ------------------------------------------------------------------------
    @property
    def path(self) -> Path:
        """The directory of the model package, `root`."""
        return self.root

    @property
    def data_paths(self) -> tuple[Path, ...]:
        """The directories of the study data, `data_path` as a tuple."""
        if self.data_path is None:
            return ()
        if isinstance(self.data_path, Path):
            return (self.data_path,)
        return tuple(self.data_path)

    @property
    def model_paths(self) -> dict[Sex, Path]:
        """The builds of the model by sex, the male build first.

        A model is built for one sex (the reference values of the PBPK library
        are those of a sex, `pkpdsim.pbpk.PBPKDefinition.sex`); a model with a
        female build has both, a model without one the male build only. The
        builds have the same ids, so the changes of an experiment and the
        parameters of a fit apply to both.
        """
        return {Sex(name): path for name, path in self.variants.items()}

    @property
    def female_model_path(self) -> Path | None:
        """The female build, `None` for a model without one."""
        return self.variants.get(Sex.FEMALE.value)

    def model_path_of(self, sex: Sex | str) -> Path:
        """Get the build of a sex.

        Args:
            sex: the sex.

        Returns:
            The SBML model of the sex.

        Raises:
            ValueError: if the model has no build of the sex.
        """
        path = self.model_paths.get(Sex(sex))
        if path is None:
            raise ValueError(
                f"The model '{self.name}' has no {Sex(sex)} build: give the model "
                f"object its {Sex(sex)} build (`{Sex(sex)}_model_path`), which the "
                "factory creates, or simulate the study on the male build."
            )
        return path

    @override
    def load_roadrunner(self, variant: str | None = None) -> roadrunner.RoadRunner:
        """Load the build of a sex in roadrunner at the level of the model.

        See `sbmlsim.project.Project.load_roadrunner`: the loader of the tests
        of a model.

        Args:
            variant: the sex of the build (`pkpdsim.sex.Sex`), the male build
                without.

        Returns:
            The roadrunner instance.

        Raises:
            ValueError: if the model has no build of the sex.
        """
        from sbmlsim.model import model_roadrunner

        path = self.model_path_of(Sex.MALE if variant is None else variant)
        return model_roadrunner.load_roadrunner(path, self.simulator.jit_optimization)

    @property
    def results_dir(self) -> Path:
        """The directory of the results."""
        return self.results_path or self.root / "results"

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
    def credibility_dir(self) -> Path:
        """The directory of the credibility assessments."""
        return self.results_dir / "credibility"

    # ------------------------------------------------------------------------
    # fits
    # ------------------------------------------------------------------------
    @override
    def adopted_fits(self) -> dict[str, str]:
        """Get the adopted fits by their subset, those of `adopted_ids`.

        Returns:
            The fit id by subset; the adopted fits of the file of the project
            for a model without `adopted_ids`.

        Raises:
            ValueError: if an adopted fit has no id of a fit, or two adopted
                fits are fits of one subset.
        """
        if not self.adopted_ids:
            return super().adopted_fits()
        from sbmlsim.fit.result import fit_subset

        fits: dict[str, str] = {}
        for fit_id in self.adopted_ids:
            subset = fit_subset(fit_id)
            if subset in fits:
                raise ValueError(
                    f"The adopted fits '{fits[subset]}' and '{fit_id}' of "
                    f"'{self.name}' are fits of one subset '{subset}'."
                )
            fits[subset] = fit_id
        return fits

    # ------------------------------------------------------------------------
    # credibility and experiments
    # ------------------------------------------------------------------------
    @override
    def load_credibility(self) -> Credibility:
        """Get the declaration of the contexts of use.

        Raises:
            ValueError: if the model has no declaration.
            TypeError: if the reference is no `Credibility`.
        """
        from pkpdsim.credibility import Credibility

        declaration = super().load_credibility()
        if not isinstance(declaration, Credibility):
            raise TypeError(
                f"The credibility of '{self.name}' ({self.credibility}) is a "
                f"{type(declaration).__name__}, not a Credibility."
            )
        return declaration

    # the override keeps the signature of the model packages of pkdb_models
    # (`experiment_classes`, `output_dir`, the result a `Path`); the run layer
    # of sbmlsim calls `Project.run_experiments` unbound, so it is not affected
    @override
    def run_experiments(  # ty: ignore[invalid-method-override]
        self,
        experiment_classes: type[SimulationExperiment]
        | Sequence[type[SimulationExperiment]],
        output_dir: str | Path,
        *,
        save_results: bool | None = None,
        show_report: bool = False,
        figure_formats: Sequence[str] | None = None,
        simulator: SimulatorSettings | None = None,
        n_workers: int | None = RUN_WORKERS,
        uncertainty: UncertaintySettings | None = None,
    ) -> Path:
        """Run simulation experiments and create their report.

        The run of the model packages, which returns the directory of the
        results and resolves a relative directory below the results of the
        model (`simulation_dir`), see `pkpdsim.run.run_experiments`; a run with
        a failed experiment or figure raises
        `sbmlsim.experiment.ExperimentRunError` after its report is written.
        It runs serially by default (`RUN_WORKERS`).
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
            n_workers=n_workers,
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
        n_workers: int | None = RUN_WORKERS,
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
            n_workers=n_workers,
            uncertainty=uncertainty,
        )
