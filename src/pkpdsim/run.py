"""Running the simulation experiments of a model and collecting their figures.

This is the `helpers.run_experiments` and the `run_simulation_experiments` of
`simulations.py` which every model of pkdb_models had as a copy: the
experiments run with the simulator settings of the model, the figures are
drawn with its figure settings, the HTML report is written, and
`run_simulation_experiments` copies the figures into `_figures`.
"""

from __future__ import annotations

import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from sbmlsim import display, log
from sbmlsim.experiment import ExperimentRunner, SimulationExperiment
from sbmlsim.report.experiment_report import ExperimentReport, ReportResults
from sbmlsim.simulator.simulation_serial import SimulatorSerial

if TYPE_CHECKING:
    from pkpdsim.model import PKPDModel

#: the formats of the figures which are copied into `_figures`
COLLECTED_FORMATS: tuple[str, ...] = ("png", "svg")


def run_experiments(
    model: PKPDModel,
    experiment_classes: type[SimulationExperiment]
    | Sequence[type[SimulationExperiment]],
    output_dir: str | Path,
    *,
    save_results: bool | None = None,
    show_report: bool = False,
    figure_formats: Sequence[str] | None = None,
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
            default.

    Returns:
        The directory the results were written to.
    """
    log.enable_rich_logging()
    output_path = model.simulation_dir / output_dir
    classes: list[type[SimulationExperiment]]
    if isinstance(experiment_classes, Sequence):
        classes = list(experiment_classes)
    else:
        classes = [experiment_classes]
    # the figures of the experiments are created with the runner, the settings
    # hold from there to the end of the drawing
    with model.figures.applied():
        runner = ExperimentRunner(
            experiment_classes=classes,
            data_path=list(model.data_paths) or None,
            base_path=model.path,
            simulator=SimulatorSerial(
                model=model.model_path,
                absolute_tolerance=model.simulator.absolute_tolerance,
                relative_tolerance=model.simulator.relative_tolerance,
            ),
        )
        results = runner.run_experiments(
            output_path=output_path,
            # the figures are written, a window per figure blocks the run
            show_figures=False,
            save_results=model.save_results if save_results is None else save_results,
            figure_formats=list(figure_formats or model.figures.formats),
            reduced_selections=model.reduced_selections,
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
    return output_path


def run_simulation_experiments(
    model: PKPDModel,
    experiment_classes: Sequence[type[SimulationExperiment]],
    output_dir: str | Path,
    *,
    show_report: bool = False,
    figure_formats: Sequence[str] | None = None,
) -> Path:
    """Run simulation experiments and collect their figures.

    The static figures of all experiments (`COLLECTED_FORMATS`) are copied into
    `_figures` of the output directory.

    Args:
        model: the model whose experiments run.
        experiment_classes: the experiments to run.
        output_dir: directory of the results, see `run_experiments`.
        show_report: open the report in a web browser.
        figure_formats: the formats of the figures, those of the model by
            default.

    Returns:
        The directory of the collected figures.
    """
    output_path = run_experiments(
        model,
        list(experiment_classes),
        output_dir,
        show_report=show_report,
        figure_formats=figure_formats,
    )
    figures_dir = output_path / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for suffix in COLLECTED_FORMATS:
        for path in output_path.glob(f"**/*.{suffix}"):
            if path.parent != figures_dir:
                shutil.copy2(path, figures_dir / path.name)
    display.link("figures", figures_dir)
    return figures_dir
