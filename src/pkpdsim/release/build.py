"""The build of a release repository.

The build runs the model package inside the release repository (`uv run`, in the
environment of its locked, released dependencies) and copies the generated
models, figures, simulation reports and fit reports to their places in the
repository. The exported package is never imported: the model object of the
source package only gives the paths of the model relative to its package folder
and the adopted fits.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from pkpdsim.model import PKPDModel
from pkpdsim.release.source import ReleaseError, Runner, run, stage

_FACTORY_PATTERNS = ("*.xml", "*.md", "*.ant")
_LICENSES = ("mit.txt", "cc-by-sa-4.0.txt")


class BuildPaths:
    """The folders of the release repository, relative to its root."""

    MODELS = "models"
    FIGURES = "figures"
    SIMULATIONS = "simulations"
    UNCERTAINTY = "simulations_uncertainty"
    FITS = "fit"
    RESULTS = "results"


def _replace(source: Path, target: Path) -> None:
    """Replace the folder `target` by a copy of the folder `source`."""
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(source, target)


def _clear(path: Path) -> None:
    """Remove the result folder of a run, so stale results are never reused."""
    if path.exists():
        shutil.rmtree(path)


def _step(runner: Runner, args: list[str], cwd: Path) -> None:
    """Run one step of the build with its output shown; a failure names the step."""
    try:
        runner(args, cwd, stream=True)
    except ReleaseError as err:
        raise ReleaseError(f"Build step '{' '.join(args)}' failed: {err}") from err


def _require(path: Path, step: str) -> Path:
    """Check that a step created the folder `path`."""
    if not path.is_dir():
        raise ReleaseError(f"The build step '{step}' did not create '{path}'.")
    return path


def _copy_factory_output(
    model: PKPDModel, package_root: Path, release_path: Path
) -> list[Path]:
    """Copy the output of the factory to the places in the release repository."""
    written: list[Path] = []
    factory_dir = package_root / model.factory_dir.relative_to(model.path)
    output_dir = factory_dir.parent

    models_dir = release_path / BuildPaths.MODELS
    if models_dir.exists():
        shutil.rmtree(models_dir)
    models_dir.mkdir(parents=True)
    for pattern in _FACTORY_PATTERNS:
        for source in sorted(factory_dir.glob(pattern)):
            shutil.copy2(source, models_dir / source.name)
            written.append(models_dir.relative_to(release_path) / source.name)
    if not any(models_dir.glob("*.xml")):
        raise ReleaseError(
            f"The factory step created no SBML model in '{factory_dir}'."
        )

    if model.archive_path is not None:
        archive = package_root / model.archive_path.relative_to(model.path)
        if not archive.is_file():
            raise ReleaseError(
                f"The factory step did not create the archive '{archive}'."
            )
        shutil.copy2(archive, release_path / archive.name)
        written.append(Path(archive.name))
    for name in _LICENSES:
        if (output_dir / name).is_file():
            shutil.copy2(output_dir / name, release_path / name)
            written.append(Path(name))

    figures_dir = output_dir / "figures"
    if figures_dir.is_dir():
        _replace(figures_dir, release_path / BuildPaths.FIGURES)
        written.append(Path(BuildPaths.FIGURES))
    return written


def build(
    model: PKPDModel,
    release_path: Path,
    runner: Runner = run,
    workers: str = "auto",
) -> list[Path]:
    """Build the models, simulations and fit reports in the release repository.

    All commands run with `release_path` as working directory. The results of the
    runs are written to the ignored folder `results` and copied to the tracked
    folders `models`, `figures`, `simulations`, `simulations_uncertainty` and
    `fit`; `fit` is replaced, so reports of fits which are no longer adopted
    disappear. `run_report <package> --adopted` writes the report of every
    adopted fit into a folder named by its id, which becomes `fit/<fit id>`.

    A declared output is required, an undeclared one is optional: the factory must
    create at least one SBML model and, if `model.archive_path` is set, the
    archive; the licenses and the `figures` folder are copied if the factory
    creates them.

    Args:
        model: the model object of the source package.
        release_path: the release repository, with the package already exported.
        runner: runs the commands.
        workers: the number of workers of the simulations, `auto` by default.

    Returns:
        The written paths, relative to `release_path`, with the new `uv.lock`.

    Raises:
        ReleaseError: if a step fails; the message names the step.
    """
    package = model.package
    package_root = release_path / "src" / Path(*package.split("."))
    results = BuildPaths.RESULTS
    run_simulation = ["uv", "run", "run_simulation", package]

    stage("Lock and sync of the environment")
    _step(runner, ["uv", "lock"], release_path)
    _step(runner, ["uv", "sync", "--locked"], release_path)

    stage("Models")
    _step(runner, [*run_simulation, "--action", "factory", "-r", results], release_path)
    written = [
        Path("uv.lock"),
        *_copy_factory_output(model, package_root, release_path),
    ]

    simulate = [
        *run_simulation,
        "--action",
        "simulate",
        "--experiments",
        "all",
        "-r",
        results,
        "--workers",
        workers,
    ]
    stage("Simulations")
    _clear(release_path / results / "simulation" / "all")
    _step(runner, simulate, release_path)
    _replace(
        _require(release_path / results / "simulation" / "all", "simulate"),
        release_path / BuildPaths.SIMULATIONS,
    )
    written.append(Path(BuildPaths.SIMULATIONS))

    stage("Simulations with the confidence interval")
    _clear(release_path / results / "simulation_uncertainty" / "all")
    _step(runner, [*simulate, "--uncertainty"], release_path)
    _replace(
        _require(
            release_path / results / "simulation_uncertainty" / "all",
            "simulate --uncertainty",
        ),
        release_path / BuildPaths.UNCERTAINTY,
    )
    written.append(Path(BuildPaths.UNCERTAINTY))

    fits_dir = release_path / BuildPaths.FITS
    if fits_dir.exists():
        shutil.rmtree(fits_dir)
    fits_dir.mkdir()
    reports = Path(results) / "fit_reports"
    _clear(release_path / reports)
    if model.fits:
        # one folder per adopted fit, named by its id, with the HTML report, the
        # methods report (PDF and Typst source) and the figures
        stage("Fit reports")
        _step(
            runner,
            # POSIX form, as every path of a command of the build
            ["uv", "run", "run_report", package, "--adopted", "-o", reports.as_posix()],
            release_path,
        )
    for fit_id in model.fits:
        _replace(
            _require(release_path / reports / fit_id, f"run_report --adopted {fit_id}"),
            fits_dir / fit_id,
        )
        written.append(Path(BuildPaths.FITS) / fit_id)
    return written
