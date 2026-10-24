"""The run of the simulation experiments of a model."""

import dataclasses
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.image
import pytest
from sbmlsim.experiment import ExperimentRunError, ExperimentRunner
from sbmlsim.plot import FigureSettings as SbmlsimFigureSettings

import pkpdsim.run
from pkpdsim.model import FigureSettings, SimulatorSettings
from pkpdsim.run import create_simulator
from pkpdsim_toy import MODEL
from pkpdsim_toy.broken import ToyCannotInitialize, ToyFigureFails, ToyRefusedUnit
from pkpdsim_toy.experiments.decay import ToyDecay, ToyDecaySlow


def test_a_run_writes_the_report_and_the_figures(tmp_path: Path) -> None:
    output = MODEL.run_experiments(ToyDecay, output_dir=tmp_path / "decay")
    assert output == tmp_path / "decay"
    assert list(output.glob("**/*.html"))
    assert list(output.glob("**/*fig_decay*.svg"))
    assert list(output.glob("**/*fig_decay*.png"))


def test_the_formats_of_a_run(tmp_path: Path) -> None:
    output = MODEL.run_experiments(
        ToyDecay, output_dir=tmp_path, figure_formats=["png"]
    )
    assert list(output.glob("**/*.png"))
    assert not list(output.glob("**/*.svg"))


def test_a_relative_output_dir_is_below_the_simulation_results(
    tmp_path: Path,
) -> None:
    model = dataclasses.replace(MODEL, results_path=tmp_path)
    output = model.run_experiments(ToyDecay, output_dir="relative")
    assert output == tmp_path / "simulation" / "relative"
    assert list(output.glob("**/*.html"))


def test_the_static_resolution_is_applied(tmp_path: Path) -> None:
    """`savefig_dpi` sets the size of the images in pixels."""
    widths = []
    for dpi in [20, 40]:
        model = dataclasses.replace(
            MODEL, figures=FigureSettings(formats=("png",), savefig_dpi=dpi)
        )
        output = model.run_experiments(ToyDecay, output_dir=tmp_path / str(dpi))
        (png,) = output.glob("**/*fig_decay*.png")
        widths.append(matplotlib.image.imread(png).shape[1])
    assert abs(widths[1] - 2 * widths[0]) <= 2


def test_the_figures_are_collected(tmp_path: Path) -> None:
    figures = MODEL.run_simulation_experiments(
        [ToyDecay, ToyDecaySlow], output_dir=tmp_path / "all"
    )
    assert figures == tmp_path / "all" / "_figures"
    names = sorted(path.name for path in figures.iterdir())
    assert any(name.startswith("ToyDecay_") for name in names)
    assert any(name.startswith("ToyDecaySlow_") for name in names)
    assert {Path(name).suffix for name in names} == {".png", ".svg"}


def test_the_figures_of_a_relative_output_dir_are_collected(tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, results_path=tmp_path)
    figures = model.run_simulation_experiments([ToyDecay], output_dir="relative")
    assert figures == tmp_path / "simulation" / "relative" / "_figures"
    assert list(figures.glob("*.svg"))


def test_a_bare_string_is_not_a_list_of_formats(tmp_path: Path) -> None:
    with pytest.raises(TypeError, match="not 'png'"):
        MODEL.run_experiments(
            ToyDecay,
            output_dir=tmp_path,
            figure_formats="png",
        )


def test_empty_formats_fall_back_to_the_formats_of_the_model(tmp_path: Path) -> None:
    output = MODEL.run_experiments(ToyDecay, output_dir=tmp_path, figure_formats=[])
    assert list(output.glob("**/*.png"))
    assert list(output.glob("**/*.svg"))


def test_a_run_restores_the_matplotlib_settings(tmp_path: Path) -> None:
    dpi = matplotlib.rcParams["savefig.dpi"]
    MODEL.run_experiments(ToyDecay, output_dir=tmp_path)
    assert matplotlib.rcParams["savefig.dpi"] == dpi


def test_the_figure_settings_of_the_model_reach_sbmlsim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}
    original = ExperimentRunner.run_experiments

    def run_experiments(
        self: ExperimentRunner,
        *args: Any,  # noqa: ANN401, the arguments of sbmlsim, passed on
        **kwargs: Any,  # noqa: ANN401
    ) -> object:
        captured.update(kwargs)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(ExperimentRunner, "run_experiments", run_experiments)
    MODEL.run_experiments(ToyDecay, output_dir=tmp_path)
    settings = captured["figure_settings"]
    assert isinstance(settings, SbmlsimFigureSettings)
    assert settings == MODEL.figures.sbmlsim_settings()
    assert settings.legend_fontsize == 7


def test_the_simulator_has_the_settings_of_the_model() -> None:
    simulator = create_simulator(
        SimulatorSettings(absolute_tolerance=1e-9, relative_tolerance=1e-8, n_workers=3)
    )
    assert simulator.n_workers == 3
    assert simulator.integrator_settings["absolute_tolerance"] == 1e-9
    assert simulator.integrator_settings["relative_tolerance"] == 1e-8


def test_a_run_takes_other_simulator_settings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[SimulatorSettings] = []
    original = pkpdsim.run.create_simulator

    def create(settings: SimulatorSettings) -> object:
        seen.append(settings)
        return original(settings)

    monkeypatch.setattr(pkpdsim.run, "create_simulator", create)
    MODEL.run_experiments(
        ToyDecay, output_dir=tmp_path, simulator=SimulatorSettings(n_workers=2)
    )
    assert [s.n_workers for s in seen] == [2]
    MODEL.run_experiments(ToyDecay, output_dir=tmp_path / "model")
    assert seen[-1] == MODEL.simulator


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


def test_the_figures_are_collected_before_the_failures_are_raised(
    tmp_path: Path,
) -> None:
    with pytest.raises(ExperimentRunError):
        MODEL.run_simulation_experiments(
            [ToyDecay, ToyFigureFails], output_dir=tmp_path / "all"
        )
    assert list((tmp_path / "all" / "_figures").glob("ToyDecay_*.svg"))
