"""The tests of a model as functions: the registry, the experiments, the fits.

A repository of models calls them from its own tests, so that a model which is
exported into a repository of its own keeps them:

```python
# not executed
@pytest.mark.parametrize("model", MODELS)
def test_experiments_registered(model: PKPDModel) -> None:
    testing.check_registry(model, not_registered={"Dussol2005"})
```

The functions raise `AssertionError` with a message which names what is wrong.
"""

from __future__ import annotations

import importlib
import inspect
import pkgutil
import traceback
from collections.abc import Collection, Generator, Sequence
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING

from sbmlsim.experiment import ExperimentRunner, SimulationExperiment
from sbmlsim.fit.cli import FitRun, run_fit
from sbmlsim.plot import Figure

if TYPE_CHECKING:
    from pkpdsim.model import PKPDModel


@contextmanager
def restored_figure_settings() -> Generator[None]:
    """Restore the class attributes of sbmlsim's `Figure` after the context.

    A test of a model which sets them directly (not through `FigureSettings`)
    leaks them into the tests of the next model otherwise.
    """
    settings = {
        key: value
        for key, value in vars(Figure).items()
        if not key.startswith("_")
        and not callable(value)
        and not isinstance(value, property | classmethod | staticmethod)
    }
    try:
        yield
    finally:
        for key, value in settings.items():
            setattr(Figure, key, value)


def experiment_modules(
    model: PKPDModel, not_importable: Collection[str] = ()
) -> tuple[list[ModuleType], dict[str, str]]:
    """Import the modules below `experiments` of a model.

    Args:
        model: the model.
        not_importable: modules which are skipped, by their name relative to the
            package of the model, e.g. `experiments.zonation.plot_pyvista`.

    Returns:
        The imported modules and the traceback of every module which does not
        import, by its name relative to the package of the model.
    """
    prefix = f"{model.package}."
    package = importlib.import_module(f"{prefix}experiments")
    modules: list[ModuleType] = [package]
    errors: dict[str, str] = {}

    def on_error(name: str) -> None:
        errors[name.removeprefix(prefix)] = traceback.format_exc()

    for info in pkgutil.walk_packages(
        package.__path__, prefix=f"{package.__name__}.", onerror=on_error
    ):
        if info.name.removeprefix(prefix) in not_importable:
            continue
        try:
            modules.append(importlib.import_module(info.name))
        except Exception:  # every module which does not import is reported
            on_error(info.name)
    return modules, errors


def defined_experiments(
    modules: Sequence[ModuleType],
) -> set[type[SimulationExperiment]]:
    """Get the experiments defined in modules, without the base experiments.

    A base experiment is a class from which another of the found classes derives.
    """
    classes = {
        cls
        for module in modules
        for _, cls in inspect.getmembers(module, inspect.isclass)
        if cls.__module__ == module.__name__ and issubclass(cls, SimulationExperiment)
    }
    return {
        cls
        for cls in classes
        if not any(other is not cls and issubclass(other, cls) for other in classes)
    }


def registered_experiments(model: PKPDModel) -> list[type[SimulationExperiment]]:
    """Get the experiments of the group `all` of the registry of a model."""
    return list(model.load_experiments()["all"])


def check_registry(
    model: PKPDModel,
    not_registered: Collection[str] = (),
    not_importable: Collection[str] = (),
) -> None:
    """Check that every experiment a model defines is in its registry.

    Every module below `experiments` imports, every experiment is registered in
    `all` once, and the known exceptions are still exceptions.

    Args:
        model: the model.
        not_registered: names of experiment classes which are defined but not
            registered on purpose.
        not_importable: modules which are skipped, see `experiment_modules`.

    Raises:
        AssertionError: naming every problem.
    """
    modules, errors = experiment_modules(model, not_importable)
    problems = [f"{name} does not import:\n{tb}" for name, tb in errors.items()]
    names = [cls.__name__ for cls in registered_experiments(model)]
    twice = sorted({name for name in names if names.count(name) > 1})
    if twice:
        problems.append(f"registered twice: {twice}")
    unregistered = sorted(
        cls.__name__
        for cls in defined_experiments(modules)
        if cls.__name__ not in names and cls.__name__ not in not_registered
    )
    if unregistered:
        problems.append(f"defined, but not in the group 'all': {unregistered}")
    stale = sorted(set(not_registered) & set(names))
    if stale:
        problems.append(f"registered, but listed as not registered: {stale}")
    if problems:
        raise AssertionError("\n".join(problems))


def _runner(model: PKPDModel) -> ExperimentRunner:
    """Initialize the registered experiments of a model, as its runs do."""
    return ExperimentRunner(
        experiment_classes=registered_experiments(model),
        base_path=model.path,
        data_path=list(model.data_paths) or None,
        simulator=None,
    )


def initialize_experiments(model: PKPDModel) -> dict[str, SimulationExperiment]:
    """Initialize the registered experiments with their data, figures, mappings.

    Returns:
        The experiments by name.
    """
    return dict(_runner(model).experiments)


def run_all_experiments(
    model: PKPDModel, output_dir: Path, one_at_a_time: bool = False
) -> None:
    """Run the registered experiments, as `run_<model>` does.

    Args:
        model: the model.
        output_dir: the directory of the results.
        one_at_a_time: run every experiment on its own, into a directory named
            after it: the memory of a run of all experiments of a large model
            exceeds that of a runner of a continuous integration (sbmlsim keeps
            the results and figures of a run until its end).

    Raises:
        AssertionError: if a run wrote no report.
    """
    if one_at_a_time:
        for cls in registered_experiments(model):
            path = model.run_experiments(cls, output_dir=output_dir / cls.__name__)
            if not list(path.glob("**/*.html")):
                raise AssertionError(f"{cls.__name__} wrote no report")
        return
    path = output_dir / "all"
    model.run_simulation_experiments(registered_experiments(model), path)
    if not list(path.glob("**/*.html")):
        raise AssertionError("the run of all experiments wrote no report")


def start_fits(
    model: PKPDModel, output_dir: Path, max_nfev: int = 2
) -> dict[str, dict[str, FitRun]]:
    """Start a fit of every fit definition: one start with a few evaluations.

    Args:
        model: the model.
        output_dir: the directory of the results of the fits.
        max_nfev: evaluations of the cost per start.

    Returns:
        The runs of every fit definition.

    Raises:
        AssertionError: if a fit definition did not run.
    """
    runs: dict[str, dict[str, FitRun]] = {}
    for key, definition in model.load_fit_definitions().items():
        runs[key] = run_fit(
            definition,
            opid=f"test_{key}",
            size=1,
            seed=1234,
            output_dir=output_dir,
            max_nfev=max_nfev,
        )
        if not runs[key]:
            raise AssertionError(f"the fit {key} did not run")
    return runs
