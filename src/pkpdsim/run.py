"""Running the simulation experiments of a model and collecting their figures.

This is the `helpers.run_experiments` and the `run_simulation_experiments` of
`simulations.py` which every model of pkdb_models had as a copy: the
experiments run with the simulator settings of the model, the figures are
drawn with its figure settings, the HTML report is written, and
`run_simulation_experiments` copies the figures into `_figures`.

Every experiment runs, also after another one failed: sbmlsim logs a failure
with its traceback and records it. Once the report (and the collected figures)
are written, a run with a failure raises `ExperimentRunError`, whose
`failures` name every failed experiment and figure.

With `uncertainty` (`pkpdsim.uncertainty.UncertaintySettings`) the experiments
run over draws of the fitted parameters of the model: the parameters are drawn
once per run, every experiment gets the draws as its setting `uncertainty` and
runs every simulation over them, see `pkpdsim.uncertainty.apply_uncertainty`.
A relative output directory is then below `model.uncertainty_dir`.
"""

from __future__ import annotations

import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import matplotlib
from sbmlsim import display, log
from sbmlsim.experiment import (
    ExperimentRunError,
    ExperimentRunner,
    SimulationExperiment,
)
from sbmlsim.report.experiment_report import ExperimentReport, ReportResults
from sbmlsim.simulator import Simulator

from pkpdsim.model import normalize_formats
from pkpdsim.uncertainty import UncertaintyRun

if TYPE_CHECKING:
    from sbmlsim.experiment import ExperimentResult

    from pkpdsim.model import PKPDModel, SimulatorSettings
    from pkpdsim.uncertainty import UncertaintySettings

#: the formats of the figures which are copied into `_figures`
COLLECTED_FORMATS: tuple[str, ...] = ("png", "svg")


def create_simulator(settings: SimulatorSettings) -> Simulator:
    """Get the simulator of sbmlsim with the settings of a model.

    Args:
        settings: the tolerances and the processes of the simulations.

    Returns:
        The simulator.
    """
    return Simulator(
        n_workers=settings.n_workers,
        absolute_tolerance=settings.absolute_tolerance,
        relative_tolerance=settings.relative_tolerance,
    )


def _run(
    model: PKPDModel,
    experiment_classes: type[SimulationExperiment]
    | Sequence[type[SimulationExperiment]],
    output_dir: str | Path,
    *,
    save_results: bool | None,
    show_report: bool,
    figure_formats: Sequence[str] | None,
    simulator: SimulatorSettings | None,
    uncertainty: UncertaintySettings | None,
) -> tuple[Path, list[ExperimentResult]]:
    """Run experiments and write their report, see `run_experiments`.

    Returns:
        The directory of the results and the results of the experiments.
    """
    formats = normalize_formats(figure_formats) if figure_formats is not None else ()
    log.enable_rich_logging()
    settings: dict[str, Any] = {}
    if uncertainty is None:
        output_path = model.simulation_dir / output_dir
    else:
        output_path = model.uncertainty_dir / output_dir
        run = UncertaintyRun.create(model, uncertainty)
        display.key_values(
            {
                "uncertainty": f"{uncertainty.draws} draws of "
                f"{len(run.dimension.values)} parameters, seed {uncertainty.seed}, "
                f"bands {uncertainty.quantiles[0]:g} to {uncertainty.quantiles[1]:g}"
            }
        )
        settings["uncertainty"] = run
    classes: list[type[SimulationExperiment]]
    if isinstance(experiment_classes, Sequence):
        classes = list(experiment_classes)
    else:
        classes = [experiment_classes]
    runner = ExperimentRunner(
        experiment_classes=classes,
        data_path=list(model.data_paths) or None,
        base_path=model.path,
        simulator=create_simulator(simulator or model.simulator),
        # an experiment which cannot be initialized is logged and reported
        on_error="log",
        **settings,
    )
    # the resolution of the static images holds while the figures are drawn
    rc: dict[Any, Any] = model.figures.rc_params()
    with matplotlib.rc_context(rc):
        results = runner.run_experiments(
            output_path=output_path,
            # the figures are written, a window per figure blocks the run
            show_figures=False,
            save_results=model.save_results if save_results is None else save_results,
            figure_formats=list(formats or model.figures.formats),
            reduced_selections=model.reduced_selections,
            figure_settings=model.figures.sbmlsim_settings(),
        )

    report_results = ReportResults()
    for result in results:
        report_results.add_experiment_result(exp_result=result)
    # the HTML report, its section ends the output with the link to the report
    report = ExperimentReport(report_results, metadata=None)
    report.create_report(
        output_path,
        report_type=ExperimentReport.ReportType.HTML,
        show_report=show_report,
    )
    return output_path, results


def _raise_failures(results: list[ExperimentResult]) -> None:
    """Raise the failures of a run, after its outputs are written.

    Raises:
        ExperimentRunError: if an experiment or a figure failed.
    """
    if any(result.failed for result in results):
        raise ExperimentRunError(results)


def run_experiments(
    model: PKPDModel,
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
    """Run simulation experiments of a model and create their report.

    sbmlsim does not configure logging, so the messages of the package are shown
    on the rich console explicitly.

    Args:
        model: the model whose experiments run.
        experiment_classes: the experiment class or classes to run.
        output_dir: directory of the results; a relative directory is below the
            simulation results of the model (`model.simulation_dir`).
        save_results: store the results of the simulations,
            `model.save_results` by default.
        show_report: open the report in a web browser.
        figure_formats: the formats of the figures, those of the model by
            default; an empty sequence falls back to the formats of the model
            as well.
        simulator: the settings of the simulator, those of the model by
            default, e.g. with more processes.
        uncertainty: run the experiments over draws of the fitted parameters
            of the model, see the module; `None` runs them once.

    Returns:
        The directory the results were written to.

    Raises:
        TypeError: if `figure_formats` is a single string.
        ExperimentRunError: after the report is written, if an experiment or a
            figure failed; its `failures` name them, the log of the run has the
            tracebacks.
    """
    output_path, results = _run(
        model,
        experiment_classes,
        output_dir,
        save_results=save_results,
        show_report=show_report,
        figure_formats=figure_formats,
        simulator=simulator,
        uncertainty=uncertainty,
    )
    _raise_failures(results)
    return output_path


def run_simulation_experiments(
    model: PKPDModel,
    experiment_classes: Sequence[type[SimulationExperiment]],
    output_dir: str | Path,
    *,
    show_report: bool = False,
    figure_formats: Sequence[str] | None = None,
    simulator: SimulatorSettings | None = None,
    uncertainty: UncertaintySettings | None = None,
) -> Path:
    """Run simulation experiments and collect their figures.

    The static figures of all experiments (`COLLECTED_FORMATS`) are copied into
    `_figures` of the output directory, also those of a run with failures.

    Args:
        model: the model whose experiments run.
        experiment_classes: the experiments to run.
        output_dir: directory of the results, see `run_experiments`.
        show_report: open the report in a web browser.
        figure_formats: the formats of the figures, those of the model by
            default.
        simulator: the settings of the simulator, those of the model by default.
        uncertainty: run the experiments over draws of the fitted parameters,
            see `run_experiments`.

    Returns:
        The directory of the collected figures.

    Raises:
        ExperimentRunError: after the figures are collected, if an experiment
            or a figure failed.
    """
    output_path, results = _run(
        model,
        list(experiment_classes),
        output_dir,
        save_results=None,
        show_report=show_report,
        figure_formats=figure_formats,
        simulator=simulator,
        uncertainty=uncertainty,
    )
    figures_dir = output_path / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for suffix in COLLECTED_FORMATS:
        for path in output_path.glob(f"**/*.{suffix}"):
            if path.parent != figures_dir:
                shutil.copy2(path, figures_dir / path.name)
    display.link("figures", figures_dir)
    _raise_failures(results)
    return figures_dir
