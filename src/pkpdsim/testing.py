"""The tests of a model as functions: the registry, the experiments, the fits.

A repository of models calls them from its own tests, so that a model which is
exported into a repository of its own keeps them (`check_model` checks that the
model object belongs to its package):

```python
# not executed
@pytest.mark.parametrize("model", MODELS)
def test_experiments_registered(model: PKPDModel) -> None:
    testing.check_registry(model)
```

The known exceptions of the tests (experiments which are not registered on
purpose, modules which do not import, the reason to run the experiments one at
a time) are the `checks` of the model object, each with its reason.

The functions raise `AssertionError` with a message which names every problem,
also for experiments which fail to initialize or to run.
"""

from __future__ import annotations

import importlib
import importlib.util
import inspect
import pkgutil
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING

from sbmlsim.experiment import (
    ExperimentRunError,
    ExperimentRunner,
    SimulationExperiment,
)
from sbmlsim.fit.cli import FitRun, run_fit

from pkpdsim.model import resolve
from pkpdsim.uncertainty import check_fit_parameters

if TYPE_CHECKING:
    from pkpdsim.model import PKPDModel


def experiment_modules(model: PKPDModel) -> tuple[list[ModuleType], dict[str, str]]:
    """Import the modules below `experiments` of a model.

    The modules in `model.checks.not_importable` are skipped.

    Args:
        model: the model.

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
        if info.name.removeprefix(prefix) in model.checks.not_importable:
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


def _module_exists(name: str) -> bool:
    """Check that a module exists without importing it (its parents are imported)."""
    try:
        return importlib.util.find_spec(name) is not None
    except ModuleNotFoundError:
        return False


def check_registry(model: PKPDModel) -> None:
    """Check that every experiment a model defines is in its registry.

    Every module below `experiments` imports, every experiment is registered in
    `all` once, and the known exceptions of `model.checks` are still exceptions
    (a listed class is defined and not registered, a listed module exists).

    Args:
        model: the model.

    Raises:
        AssertionError: naming every problem.
    """
    modules, errors = experiment_modules(model)
    checks = model.checks
    problems = [f"{name} does not import:\n{tb}" for name, tb in errors.items()]
    names = [cls.__name__ for cls in registered_experiments(model)]
    twice = sorted({name for name in names if names.count(name) > 1})
    if twice:
        problems.append(f"registered twice: {twice}")
    defined = defined_experiments(modules)
    unregistered = sorted(
        cls.__name__
        for cls in defined
        if cls.__name__ not in names and cls.__name__ not in checks.not_registered
    )
    if unregistered:
        problems.append(f"defined, but not in the group 'all': {unregistered}")
    stale = sorted(set(checks.not_registered) & set(names))
    if stale:
        problems.append(f"registered, but listed as not registered: {stale}")
    missing = sorted(
        set(checks.not_registered) - {cls.__name__ for cls in defined} - set(names)
    )
    if missing:
        problems.append(f"listed as not registered, but not defined: {missing}")
    absent = sorted(
        name
        for name in checks.not_importable
        if not _module_exists(f"{model.package}.{name}")
    )
    if absent:
        problems.append(f"listed as not importable, but no module: {absent}")
    if problems:
        raise AssertionError("\n".join(problems))


def check_model(model: PKPDModel) -> None:
    """Check that a model object belongs to its package and its references hold.

    A slip in `package` or `path` of a model resolves the references of another
    model, or names the paths wrongly, without an error otherwise. The name is
    the last component of the package (or the component ends with `_<name>`),
    the builds (`model_path`, `female_model_path`), the factory directory and
    `archive_path` lie below `path`,
    the data paths exist, and every parameter of an adopted fit has its value
    and a distribution in the model (`pkpdsim.uncertainty.check_fit_parameters`).

    Raises:
        AssertionError: naming every problem.
    """
    problems: list[str] = []
    module: ModuleType | None = None
    try:
        module = importlib.import_module(model.package)
    except ImportError as err:
        problems.append(f"the package '{model.package}' does not import: {err}")
    if module is not None:
        if getattr(module, "MODEL", None) is not model:
            problems.append(f"MODEL of '{model.package}' is not this model")
        file = getattr(module, "__file__", None)
        if file is None or model.path != Path(file).parent:
            problems.append(f"path {model.path} is not the directory of {file}")
    component = model.package.rsplit(".", 1)[-1]
    if model.name != component and not component.endswith(f"_{model.name}"):
        problems.append(
            f"the name '{model.name}' is not the name of the package '{component}'"
        )
    for label, path in [
        ("model_path", model.model_path),
        ("female_model_path", model.female_model_path),
        ("the factory directory", model.factory_dir),
        ("archive_path", model.archive_path),
    ]:
        if path is not None and not path.is_relative_to(model.path):
            problems.append(f"{label} {path} is not below the model {model.path}")
    problems += [
        f"the data path {path} does not exist"
        for path in model.data_paths
        if not path.exists()
    ]

    references = {
        "factory": model.factory,
        "experiments": model.experiments,
        "fit_definitions": model.fit_definitions,
        "extra_arguments": model.extra_arguments,
        **{f"extra action '{a.name}'": a.function for a in model.extra_actions},
    }
    resolved: dict[str, object] = {}
    for name, reference in references.items():
        if reference is None:
            continue
        try:
            resolved[name] = resolve(reference, model.package)
        except Exception as err:  # reported, whatever the reference does
            problems.append(f"{name} '{reference}' does not resolve: {err!r}")
    if "experiments" in resolved:
        registry = resolved["experiments"]
        if not isinstance(registry, Mapping) or "all" not in registry:
            problems.append("the registry of experiments is not a mapping with 'all'")
    if "fit_definitions" in resolved:
        definitions = resolved["fit_definitions"]
        if not isinstance(definitions, Mapping) or not definitions:
            problems.append("the fit definitions are not a non-empty mapping")
    # the values of the adopted fits and their uncertainty are in the model
    problems += check_fit_parameters(model)
    if problems:
        raise AssertionError("\n".join(problems))


def _runner(model: PKPDModel) -> ExperimentRunner:
    """Initialize the registered experiments of a model, as its runs do."""
    return ExperimentRunner(
        experiment_classes=registered_experiments(model),
        base_path=model.path,
        data_path=list(model.data_paths) or None,
        simulator=None,
        on_error="log",
        jit_optimization=model.simulator.jit_optimization,
    )


def initialize_experiments(model: PKPDModel) -> dict[str, SimulationExperiment]:
    """Initialize the registered experiments with their data, figures, mappings.

    Every experiment is initialized, also after one which cannot be.

    Returns:
        The experiments by name.

    Raises:
        AssertionError: naming every experiment which cannot be initialized
            with its error; the tracebacks are in the log of sbmlsim.
    """
    runner = _runner(model)
    if runner.failed:
        total = len(runner.failed) + len(runner.experiments)
        lines = [f"{sid}: {result.error}" for sid, result in runner.failed.items()]
        raise AssertionError(
            f"{len(runner.failed)} of {total} experiments cannot be initialized "
            "(the tracebacks are in the log):\n" + "\n".join(lines)
        )
    return dict(runner.experiments)


def run_all_experiments(model: PKPDModel, output_dir: Path) -> None:
    """Run the registered experiments, as `run_simulation` does.

    Every experiment runs, also after one which failed.

    Args:
        model: the model.
        output_dir: the directory of the results. With
            `model.checks.one_at_a_time` every experiment runs on its own, into
            a directory named after it: the memory of a run of all experiments
            of a large model exceeds that of a runner of a continuous
            integration.

    Raises:
        AssertionError: naming every failed experiment and figure, and every
            run which wrote no report.
    """
    failures: list[str] = []
    if model.checks.one_at_a_time is not None:
        for cls in registered_experiments(model):
            path = output_dir / cls.__name__
            try:
                model.run_experiments(cls, output_dir=path)
            except ExperimentRunError as err:
                failures += err.failures
            if not list(path.glob("**/*.html")):
                failures.append(f"{cls.__name__} wrote no report")
    else:
        path = output_dir / "all"
        try:
            model.run_simulation_experiments(registered_experiments(model), path)
        except ExperimentRunError as err:
            failures += err.failures
        if not list(path.glob("**/*.html")):
            failures.append("the run of all experiments wrote no report")
    if failures:
        raise AssertionError("\n".join(failures))


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
        # the problems load their models at the level of the model
        with model.simulator.jit_environment():
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
