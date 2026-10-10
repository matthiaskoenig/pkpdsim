"""The tolerance and the run of the legacy package, as in pkdb_models before 0.3.0."""

from pathlib import Path

import matplotlib
from sbmlsim.experiment import ExperimentRunner, SimulationExperiment
from sbmlsim.plot import Figure
from sbmlsim.report.experiment_report import ExperimentReport, ReportResults
from sbmlsim.simulator.simulation_serial import SimulatorSerial

from pkpdsim_legacy import DATA_PATHS, MODEL_PATH, PKPDSIM_LEGACY_PATH

ABSOLUTE_TOLERANCE = 1e-10
RELATIVE_TOLERANCE = 1e-10


def run_experiments(
    experiment_classes: type[SimulationExperiment] | list[type[SimulationExperiment]],
    output_dir: Path,
    save_results: bool = False,
) -> Path:
    """Run experiments with the settings of the toy model, the legacy way."""
    Figure.legend_fontsize = 7
    output_path = PKPDSIM_LEGACY_PATH / "results" / "simulation" / output_dir
    if not isinstance(experiment_classes, list):
        experiment_classes = [experiment_classes]
    simulator = SimulatorSerial(
        model=MODEL_PATH,
        absolute_tolerance=ABSOLUTE_TOLERANCE,
        relative_tolerance=RELATIVE_TOLERANCE,
    )
    runner = ExperimentRunner(
        experiment_classes=experiment_classes,
        data_path=DATA_PATHS,
        base_path=PKPDSIM_LEGACY_PATH,
        simulator=simulator,
    )
    with matplotlib.rc_context({"savefig.dpi": 30}):
        results = runner.run_experiments(
            output_path=output_path,
            show_figures=False,
            save_results=save_results,
            figure_formats=["svg", "png"],
            reduced_selections=True,
        )
    report_results = ReportResults()
    for result in results:
        report_results.add_experiment_result(exp_result=result)
    ExperimentReport(report_results, metadata=None).create_report(
        output_path, report_type=ExperimentReport.ReportType.HTML
    )
    return output_path
