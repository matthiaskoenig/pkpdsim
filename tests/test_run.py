"""The run of the simulation experiments of a model."""

import dataclasses
from pathlib import Path

import matplotlib
import matplotlib.image
import pytest
from sbmlsim.plot import Figure

from pkpdsim.model import FIGURE_ATTRIBUTES, FigureSettings
from pkpdsim_toy import MODEL
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


def test_a_run_restores_the_figure_settings(tmp_path: Path) -> None:
    before = {name: getattr(Figure, name) for name in FIGURE_ATTRIBUTES}
    dpi = matplotlib.rcParams["savefig.dpi"]
    MODEL.run_experiments(ToyDecay, output_dir=tmp_path)
    assert {name: getattr(Figure, name) for name in FIGURE_ATTRIBUTES} == before
    assert matplotlib.rcParams["savefig.dpi"] == dpi


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
