"""The model object: its references, settings, paths and validation."""

import dataclasses
from pathlib import Path

import matplotlib
import pytest
from sbmlsim.plot import Figure

from pkpdsim.model import (
    FIGURE_ATTRIBUTES,
    ExtraAction,
    FigureSettings,
    PKPDModel,
    SimulatorSettings,
    check_reference,
    resolve,
)
from pkpdsim_toy import MODEL, MODEL_BASE_PATH, TOY_PATH


def test_resolve_an_absolute_reference() -> None:
    assert resolve("pkpdsim_toy:TOY_PATH") == TOY_PATH


def test_resolve_a_relative_reference() -> None:
    registry = resolve(".simulations:EXPERIMENTS", "pkpdsim_toy")
    assert set(registry) == {"studies", "misc", "all"}


def test_resolve_a_dotted_attribute() -> None:
    assert resolve("pkpdsim_toy:MODEL.name") == "toy"


@pytest.mark.parametrize(
    "reference",
    ["simulations.EXPERIMENTS", ":EXPERIMENTS", ".simulations:", "a:b-c", "a b:c"],
)
def test_a_malformed_reference_is_rejected(reference: str) -> None:
    with pytest.raises(ValueError, match="module:attribute"):
        check_reference(reference)


def test_a_malformed_reference_fails_at_creation() -> None:
    with pytest.raises(ValueError, match=r"simulations\.EXPERIMENTS"):
        dataclasses.replace(MODEL, experiments="simulations.EXPERIMENTS")


def test_a_missing_attribute_names_the_reference() -> None:
    with pytest.raises(AttributeError, match=r"\.simulations:MISSING"):
        resolve(".simulations:MISSING", "pkpdsim_toy")


def test_the_paths_derive_from_the_path_of_the_model() -> None:
    assert MODEL.results_dir == TOY_PATH / "results"
    assert MODEL.simulation_dir == TOY_PATH / "results" / "simulation"
    assert MODEL.fit_dir == TOY_PATH / "results" / "fit"
    assert MODEL.report_dir == TOY_PATH / "results" / "report"
    assert MODEL.identifiability_dir == TOY_PATH / "results" / "identifiability"
    assert MODEL.factory_dir == MODEL_BASE_PATH


def test_given_paths_win(tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, results_path=tmp_path, factory_path=tmp_path)
    assert model.simulation_dir == tmp_path / "simulation"
    assert model.factory_dir == tmp_path


def test_the_references_load() -> None:
    assert set(MODEL.load_experiments()["all"]) == set(
        MODEL.load_experiments()["studies"] + MODEL.load_experiments()["misc"]
    )
    assert callable(MODEL.load_factory())


def test_a_model_without_experiments_says_so() -> None:
    model = dataclasses.replace(MODEL, experiments=None, fit_definitions=None)
    with pytest.raises(ValueError, match="'toy' has no registry of experiments"):
        model.load_experiments()
    with pytest.raises(ValueError, match="'toy' has no fit definitions"):
        model.load_fit_definitions()


def test_create_models(tmp_path: Path) -> None:
    assert MODEL.create_models(tmp_path) == tmp_path
    assert (tmp_path / "toy.xml").exists()


def test_an_extra_action_must_not_repeat_an_action() -> None:
    for name in ["all", "simulate"]:
        with pytest.raises(ValueError, match="repeat an action"):
            dataclasses.replace(
                MODEL,
                extra_actions=(
                    ExtraAction(name=name, help="help", function=".extras:greet"),
                ),
            )
    twice = (ExtraAction(name="greet", help="h", function=".extras:greet"),) * 2
    with pytest.raises(ValueError, match="repeat an action"):
        dataclasses.replace(MODEL, extra_actions=twice)


def test_a_bare_string_is_not_a_list_of_formats() -> None:
    with pytest.raises(TypeError, match="not 'svg'"):
        FigureSettings(formats="svg")


def test_the_formats_are_a_lower_case_tuple() -> None:
    assert FigureSettings(formats=["SVG", "Png"]).formats == ("svg", "png")


def test_the_dataclasses_are_keyword_only() -> None:
    with pytest.raises(TypeError):
        FigureSettings(("svg",))  # ty: ignore[too-many-positional-arguments]
    with pytest.raises(TypeError):
        SimulatorSettings(1e-10)  # ty: ignore[too-many-positional-arguments]
    with pytest.raises(TypeError):
        ExtraAction("a", "b", "c")  # ty: ignore[missing-argument, too-many-positional-arguments]


def test_the_sequences_are_tuples_and_the_model_is_hashable() -> None:
    model = dataclasses.replace(
        MODEL,
        data_paths=[TOY_PATH],
        extra_actions=[],
    )
    assert model.data_paths == (TOY_PATH,)
    assert model.extra_actions == ()
    assert isinstance(hash(model), int)
    assert hash(model) == hash(dataclasses.replace(model))


def test_resolve_a_relative_reference_needs_a_package() -> None:
    with pytest.raises(ValueError, match="needs a package"):
        resolve(".simulations:EXPERIMENTS")


def test_a_model_needs_a_figure_format() -> None:
    with pytest.raises(ValueError, match="figure format"):
        FigureSettings(formats=())


def test_the_settings_are_frozen() -> None:
    with pytest.raises(dataclasses.FrozenInstanceError):
        MODEL.name = "other"  # ty: ignore[invalid-assignment]
    assert SimulatorSettings().absolute_tolerance == 1e-10


def test_the_figure_attributes_exist_in_sbmlsim() -> None:
    for name in FIGURE_ATTRIBUTES:
        assert hasattr(Figure, name), name


def test_settings_are_applied_and_restored() -> None:
    before = {name: getattr(Figure, name) for name in FIGURE_ATTRIBUTES}
    dpi = matplotlib.rcParams["savefig.dpi"]
    settings = FigureSettings(fig_dpi=123, legend_fontsize=4, savefig_dpi=17)
    with settings.applied():
        assert Figure.fig_dpi == 123
        assert Figure.legend_fontsize == 4
        assert matplotlib.rcParams["savefig.dpi"] == 17
        # what is not set keeps its value
        assert Figure.panel_width == before["panel_width"]
    assert {name: getattr(Figure, name) for name in FIGURE_ATTRIBUTES} == before
    assert matplotlib.rcParams["savefig.dpi"] == dpi


def test_settings_are_restored_after_an_error() -> None:
    before = Figure.fig_dpi
    with pytest.raises(RuntimeError), FigureSettings(fig_dpi=7).applied():
        raise RuntimeError("the run failed")
    assert Figure.fig_dpi == before


def test_the_model_of_the_toy_package() -> None:
    assert isinstance(MODEL, PKPDModel)
    assert MODEL.package == "pkpdsim_toy"
    assert MODEL.figures.formats == ("svg", "png")
