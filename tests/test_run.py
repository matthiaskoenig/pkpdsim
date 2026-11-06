"""The run of the simulation experiments of a model.

`PKPDModel.run_experiments` and `run_simulation_experiments` are
`sbmlsim.project.Project.run_experiments` with the settings of the model and its
result directories, tested in sbmlsim; these tests check what the model adds.
"""

from pathlib import Path

import pytest
from sbmlsim.experiment import ExperimentRunError
from sbmlsim.project import SimulatorSettings

from pkpdsim_toy import MODEL
from pkpdsim_toy.broken import ToyCannotInitialize, ToyFigureFails, ToyRefusedUnit
from pkpdsim_toy.experiments.decay import ToyDecay


def test_a_run_writes_the_report_and_the_figures(tmp_path: Path) -> None:
    output = MODEL.run_experiments(ToyDecay, output_dir=tmp_path / "decay")
    assert output == tmp_path / "decay"
    assert list(output.glob("**/*.html"))
    assert list(output.glob("**/*fig_decay*.svg"))
    assert list(output.glob("**/*fig_decay*.png"))


def test_the_results_are_saved_as_the_model_says(tmp_path: Path) -> None:
    MODEL.run_experiments(ToyDecay, output_dir=tmp_path / "default")
    assert not list((tmp_path / "default").glob("**/*.nc"))
    model = MODEL.replace(save_results=True)
    model.run_experiments(ToyDecay, output_dir=tmp_path / "saved")
    assert list((tmp_path / "saved").glob("**/*.nc"))
    MODEL.run_experiments(ToyDecay, output_dir=tmp_path / "given", save_results=True)
    assert list((tmp_path / "given").glob("**/*.nc"))


def test_a_relative_output_dir_is_below_the_simulation_results(
    tmp_path: Path,
) -> None:
    model = MODEL.replace(results_path=tmp_path)
    output = model.run_experiments(ToyDecay, output_dir="relative")
    assert output == tmp_path / "simulation" / "relative"
    assert list(output.glob("**/*.html"))


def test_the_figures_of_a_relative_output_dir_are_collected(tmp_path: Path) -> None:
    model = MODEL.replace(results_path=tmp_path)
    figures = model.run_simulation_experiments([ToyDecay], output_dir="relative")
    assert figures == tmp_path / "simulation" / "relative" / "_figures"
    assert list(figures.glob("*.svg"))


def test_a_run_takes_other_simulator_settings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[tuple[SimulatorSettings, int | None]] = []
    original = SimulatorSettings.simulator

    def create(settings: SimulatorSettings, n_workers: int | None) -> object:
        seen.append((settings, n_workers))
        return original(settings, n_workers)

    monkeypatch.setattr(SimulatorSettings, "simulator", create)
    other = SimulatorSettings(relative_tolerance=1e-9)
    MODEL.run_experiments(
        ToyDecay, output_dir=tmp_path, simulator=other, n_workers=2, scan_workers=None
    )
    assert seen == [(other, 2)]
    # the scans serially by default, the experiments in a pool
    MODEL.run_experiments(ToyDecay, output_dir=tmp_path / "pool", n_workers=2)
    assert seen[-1] == (MODEL.simulator, 1)
    MODEL.run_experiments(ToyDecay, output_dir=tmp_path / "model")
    assert seen[-1] == (MODEL.simulator, 1)


def test_a_failing_experiment_is_reported_after_the_run(tmp_path: Path) -> None:
    with pytest.raises(ExperimentRunError) as info:
        MODEL.run_experiments(
            [ToyCannotInitialize, ToyDecay, ToyFigureFails], output_dir=tmp_path
        )
    failures = "\n".join(info.value.failures)
    assert "ToyCannotInitialize" in failures
    assert "figures_mpl()" in failures
    assert "'ToyDecay'" not in failures
    # the report is written and the experiment between the failures ran
    assert list(tmp_path.glob("*.html"))
    assert list(tmp_path.glob("ToyDecay/*fig_decay*.svg"))


def test_a_refused_unit_of_the_data_fails_its_experiment_only(tmp_path: Path) -> None:
    with pytest.raises(ExperimentRunError) as info:
        MODEL.run_experiments([ToyRefusedUnit, ToyDecay], output_dir=tmp_path)
    failures = "\n".join(info.value.failures)
    assert "ToyRefusedUnit" in failures
    assert "ValueError: The unit 'IU_toy_refused'" in failures
    assert "'microgramm' is not defined" in failures
    # the report is written and the other experiment ran
    assert list(tmp_path.glob("*.html"))
    assert list(tmp_path.glob("ToyDecay/*fig_decay*.svg"))
