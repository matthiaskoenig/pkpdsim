"""The model object: a project of sbmlsim with the keywords of the model packages."""

import dataclasses
import pickle
from pathlib import Path

import pytest
from pydantic import ValidationError
from sbmlsim.plot import FigureSettings as SbmlsimFigureSettings
from sbmlsim.project import (
    CheckSettings,
    FigureOutputs,
    Project,
    SimulatorSettings,
)

from pkpdsim.model import (
    DEFAULT_JIT_OPTIMIZATION,
    DEFAULT_SIMULATOR,
    JIT_OPTIMIZATION_VARIABLE,
    RUN_WORKERS,
    ExtraAction,
    FigureSettings,
    PKPDModel,
)
from pkpdsim.sex import Sex
from pkpdsim_toy import DATA_PATHS, MODEL, MODEL_BASE_PATH, MODEL_PATH, TOY_PATH


def keywords(**changes: object) -> dict[str, object]:
    """Get the keywords of the toy model package, with changes."""
    return {
        "name": "toy",
        "package": "pkpdsim_toy",
        "title": "Toy PK model",
        "path": TOY_PATH,
        "model_path": MODEL_PATH,
        "factory": ".factory:create_toy_models",
        "experiments": ".simulations:EXPERIMENTS",
        "fit_definitions": ".fitting:FIT_DEFINITIONS",
        "data_paths": tuple(DATA_PATHS),
        **changes,
    }


def test_a_model_is_a_project() -> None:
    assert isinstance(MODEL, Project)
    assert MODEL.root == MODEL.path == TOY_PATH
    assert MODEL.package == "pkpdsim_toy"
    assert MODEL.variants == {"male": MODEL_PATH}
    assert MODEL.model_path == MODEL_PATH
    assert MODEL.fits == ".fitting:FIT_DEFINITIONS"
    assert MODEL.data_path == [TOY_PATH]
    assert MODEL.data_paths == (TOY_PATH,)


def test_the_keywords_of_a_package_are_the_fields_of_the_project() -> None:
    model = PKPDModel(
        **keywords(
            female_model_path=MODEL_BASE_PATH / "toy_female.xml",
            fits=("PK_20260101_000000__abcde",),
            figures=FigureSettings(formats=("svg", "png"), savefig_dpi=30, fig_dpi=300),
        )
    )
    assert model.variants == {
        "male": MODEL_PATH,
        "female": MODEL_BASE_PATH / "toy_female.xml",
    }
    assert model.female_model_path == MODEL_BASE_PATH / "toy_female.xml"
    assert model.adopted_ids == ("PK_20260101_000000__abcde",)
    assert model.adopted_fits() == {"PK": "PK_20260101_000000__abcde"}
    assert model.fits == ".fitting:FIT_DEFINITIONS"
    assert model.fits_dir == MODEL_BASE_PATH.parent / "fits"
    assert model.factory_dir == MODEL_BASE_PATH
    assert model.figures == FigureOutputs(formats=("svg", "png"), savefig_dpi=30)
    assert model.figure_settings == SbmlsimFigureSettings(fig_dpi=300)


def test_relative_paths_are_relative_to_the_package() -> None:
    model = PKPDModel(
        **keywords(model_path=Path("models/toy.xml"), data_paths=(Path("data"),))
    )
    assert model.model_path == MODEL_PATH
    assert model.factory_dir == MODEL_BASE_PATH
    assert model.fits_dir == TOY_PATH / "fits"
    assert model.data_paths == (TOY_PATH / "data",)


def test_the_directories_derive_from_the_male_build() -> None:
    assert MODEL.factory_dir == MODEL_BASE_PATH
    assert MODEL.fits_dir == TOY_PATH / "fits"
    assert MODEL.results_dir == TOY_PATH / "results"
    assert MODEL.simulation_dir == TOY_PATH / "results" / "simulation"
    assert MODEL.uncertainty_dir == TOY_PATH / "results" / "simulation_uncertainty"
    assert MODEL.fit_dir == TOY_PATH / "results" / "fit"
    assert MODEL.report_dir == TOY_PATH / "results" / "report"
    assert MODEL.identifiability_dir == TOY_PATH / "results" / "identifiability"
    assert MODEL.credibility_dir == TOY_PATH / "results" / "credibility"


def test_given_paths_win(tmp_path: Path) -> None:
    model = MODEL.replace(
        results_path=tmp_path, factory_dir=TOY_PATH, fits_dir=tmp_path / "fits"
    )
    assert model.simulation_dir == tmp_path / "simulation"
    assert model.factory_dir == TOY_PATH
    assert model.fits_dir == tmp_path / "fits"


def test_the_directories_derive_again_with_other_builds(tmp_path: Path) -> None:
    build = tmp_path / "models" / "toy.xml"
    model = MODEL.replace(root=tmp_path, variants={"male": build})
    assert model.factory_dir == tmp_path / "models"
    assert model.fits_dir == tmp_path / "fits"
    moved = MODEL.replace(root=TOY_PATH.parent)
    assert moved.factory_dir == MODEL.factory_dir
    assert moved.fits_dir == MODEL.fits_dir


def test_given_directories_stay_with_other_builds(tmp_path: Path) -> None:
    build = tmp_path / "models" / "toy.xml"
    given = MODEL.replace(root=tmp_path, factory_dir=tmp_path, fits_dir=TOY_PATH)
    model = given.replace(variants={"male": build})
    assert model.factory_dir == tmp_path
    assert model.fits_dir == TOY_PATH
    derived = MODEL.replace(
        root=tmp_path, variants={"male": build}, fits_dir=tmp_path / "records"
    )
    assert derived.factory_dir == tmp_path / "models"
    assert derived.fits_dir == tmp_path / "records"


@pytest.mark.parametrize(
    ("keyword", "field"),
    [
        ("path", "root"),
        ("data_paths", "data_path"),
    ],
)
def test_a_keyword_and_its_field_are_given_once(keyword: str, field: str) -> None:
    data = keywords(**{keyword: TOY_PATH})
    data[field] = TOY_PATH
    with pytest.raises(ValidationError, match=f"'{keyword}' is the field '{field}'"):
        PKPDModel(**data)


def test_the_builds_and_the_variants_are_given_once() -> None:
    with pytest.raises(ValidationError, match="given twice"):
        PKPDModel(**keywords(variants={"male": MODEL_PATH}))
    with pytest.raises(ValidationError, match="needs its build"):
        PKPDModel(**{k: v for k, v in keywords().items() if k != "model_path"})


def test_the_fits_are_given_once() -> None:
    with pytest.raises(ValidationError, match="adopted fits are given twice"):
        PKPDModel(**keywords(fits=("A_20260101_000000__abcde",), adopted_ids=()))
    with pytest.raises(ValidationError, match="fit definitions are given twice"):
        PKPDModel(**keywords(fits=".fitting:FIT_DEFINITIONS"))


def test_the_figure_settings_are_given_once() -> None:
    with pytest.raises(ValidationError, match="given twice"):
        PKPDModel(
            **keywords(
                figures=FigureSettings(), figure_settings=SbmlsimFigureSettings()
            )
        )
    plain = PKPDModel(**keywords(figures=SbmlsimFigureSettings(legend_fontsize=3)))
    assert plain.figure_settings == SbmlsimFigureSettings(legend_fontsize=3)
    assert plain.figures == FigureOutputs()


def test_a_model_without_figures_has_those_of_sbmlsim() -> None:
    model = PKPDModel(**keywords())
    assert model.figure_settings is None
    assert model.figures == FigureOutputs()


def test_the_defaults_of_a_model() -> None:
    model = PKPDModel(**keywords())
    # the settings of the simulator with the level of pkpdsim, which the
    # creation of a model writes into them
    assert DEFAULT_SIMULATOR.jit_optimization is None
    assert model.simulator == SimulatorSettings(
        absolute_tolerance=1e-10,
        relative_tolerance=1e-10,
        jit_optimization=DEFAULT_JIT_OPTIMIZATION,
    )
    assert RUN_WORKERS == 1
    assert model.reduced_selections is True
    assert model.save_results is False
    assert model.fit_identifiability is False
    assert model.adopted_ids == ()
    assert model.checks == CheckSettings()


def test_a_model_needs_its_factory() -> None:
    with pytest.raises(ValidationError, match="needs its factory"):
        PKPDModel(**keywords(factory=None))


def test_a_malformed_reference_fails_at_creation() -> None:
    with pytest.raises(ValidationError, match=r"simulations\.EXPERIMENTS"):
        MODEL.replace(experiments="simulations.EXPERIMENTS")
    with pytest.raises(ValidationError, match="module:attribute"):
        MODEL.replace(extra_arguments="extras.add_arguments")


def test_the_references_load() -> None:
    groups = MODEL.load_experiments()
    assert set(groups["all"]) == set(groups["studies"] + groups["misc"])
    assert callable(MODEL.load_factory())
    assert set(MODEL.load_fit_definitions()) == {"DECAY"}


def test_a_model_without_experiments_and_fits() -> None:
    model = MODEL.replace(experiments=None, fits={})
    with pytest.raises(ValueError, match="'toy' has no simulation experiments"):
        model.load_experiments()
    assert model.load_fit_definitions() == {}


def test_create_models(tmp_path: Path) -> None:
    """A rebuild into a directory has the layout of the package (sbmlsim)."""
    factory_dir = tmp_path / MODEL.factory_dir.relative_to(MODEL.root)
    assert MODEL.create_models(tmp_path) == factory_dir
    assert (factory_dir / "toy.xml").exists()


def test_an_extra_action_must_not_repeat_an_action() -> None:
    for name in ["all", "simulate"]:
        with pytest.raises(ValidationError, match="repeat an action"):
            MODEL.replace(
                extra_actions=(
                    ExtraAction(name=name, help="help", function=".extras:greet"),
                ),
            )
    twice = (ExtraAction(name="greet", help="h", function=".extras:greet"),) * 2
    with pytest.raises(ValidationError, match="repeat an action"):
        MODEL.replace(extra_actions=twice)


def test_an_adopted_fit_is_given_once() -> None:
    fit = "DECAY_20260101_000000__abcde"
    with pytest.raises(ValidationError, match="repeat a fit"):
        MODEL.replace(adopted_ids=(fit, fit))


def test_the_adopted_fits_are_one_per_subset() -> None:
    model = MODEL.replace(
        adopted_ids=("DECAY_20260101_000000__abcde", "DECAY_20260102_000000__abcde")
    )
    with pytest.raises(ValueError, match="fits of one subset 'DECAY'"):
        model.adopted_fits()


def test_the_sequences_are_tuples() -> None:
    model = MODEL.replace(
        adopted_ids=["DECAY_20260101_000000__abcde"], extra_actions=[]
    )
    assert model.adopted_ids == ("DECAY_20260101_000000__abcde",)
    assert model.extra_actions == ()


def test_the_model_is_frozen() -> None:
    with pytest.raises(ValidationError, match="frozen"):
        MODEL.name = "other"


def test_a_model_pickles() -> None:
    assert pickle.loads(pickle.dumps(MODEL)) == MODEL


def test_a_bare_string_is_not_a_list_of_formats() -> None:
    with pytest.raises(TypeError, match="not 'svg'"):
        FigureSettings(formats="svg")


def test_the_figures_are_svg_by_default() -> None:
    assert FigureSettings().formats == ("svg",)


def test_the_formats_are_a_lower_case_tuple() -> None:
    assert FigureSettings(formats=["SVG", "Png"]).formats == ("svg", "png")


def test_a_model_needs_a_figure_format() -> None:
    with pytest.raises(ValueError, match="figure format"):
        FigureSettings(formats=())


def test_the_dataclasses_are_keyword_only() -> None:
    assert [f.name for f in dataclasses.fields(FigureSettings) if f.kw_only] == [
        "formats",
        "savefig_dpi",
    ]
    with pytest.raises(TypeError):
        ExtraAction("a", "b", "c")  # ty: ignore[missing-argument, too-many-positional-arguments]


def test_the_figure_settings_are_those_of_sbmlsim() -> None:
    settings = FigureSettings(fig_dpi=123, legend_fontsize=4)
    assert isinstance(settings, SbmlsimFigureSettings)
    assert settings.panel_width == SbmlsimFigureSettings().panel_width
    assert settings.fig_dpi == 123


def test_sbmlsim_gets_its_settings_without_the_formats() -> None:
    settings = FigureSettings(formats=("png",), savefig_dpi=17, legend_fontsize=4)
    handed = settings.sbmlsim_settings()
    assert type(handed) is SbmlsimFigureSettings
    assert handed == SbmlsimFigureSettings(legend_fontsize=4)
    assert FigureSettings().sbmlsim_settings() == SbmlsimFigureSettings()


def test_the_outputs_are_the_formats_and_the_resolution() -> None:
    outputs = FigureSettings(formats=("png",), savefig_dpi=17).outputs()
    assert outputs == FigureOutputs(formats=("png",), savefig_dpi=17)
    assert outputs.rc_params() == {"savefig.dpi": 17}
    assert FigureSettings().outputs().rc_params() == {}


def test_the_model_of_the_toy_package() -> None:
    assert isinstance(MODEL, PKPDModel)
    assert MODEL.package == "pkpdsim_toy"
    assert MODEL.figures.formats == ("svg", "png")


def test_a_test_loads_the_model_at_its_level(monkeypatch: pytest.MonkeyPatch) -> None:
    """`load_roadrunner` compiles the build of a sex at the level of the model."""
    import sbmlsim.model.model_roadrunner as model_roadrunner

    calls: list[tuple[object, object]] = []
    load = model_roadrunner.load_roadrunner

    def recording(source: object, optimization: object = None) -> object:
        calls.append((source, optimization))
        return load(source, optimization)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(model_roadrunner, "load_roadrunner", recording)
    simulator = SimulatorSettings(
        absolute_tolerance=1e-10, relative_tolerance=1e-10, jit_optimization="O0"
    )
    r = MODEL.replace(simulator=simulator).load_roadrunner()
    assert calls == [(MODEL.model_path, "O0")]
    assert "k" in r.model.getGlobalParameterIds()
    with pytest.raises(ValueError, match="no female build"):
        MODEL.load_roadrunner(Sex.FEMALE)


def test_the_fit_definitions_carry_the_level_of_the_model() -> None:
    simulator = SimulatorSettings(jit_optimization="O2")
    definitions = MODEL.replace(simulator=simulator).load_fit_definitions()
    assert {d.jit_optimization for d in definitions.values()} == {"O2"}
    # a model which names no level carries the level of pkpdsim
    levels = {d.jit_optimization for d in MODEL.load_fit_definitions().values()}
    assert levels == {DEFAULT_JIT_OPTIMIZATION}


def test_a_model_which_names_no_level_is_compiled_at_the_level_of_pkpdsim(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`DEFAULT_JIT_OPTIMIZATION`, which the environment variable of sbmlsim wins over."""
    monkeypatch.delenv(JIT_OPTIMIZATION_VARIABLE, raising=False)
    assert DEFAULT_JIT_OPTIMIZATION == "O0"
    # the settings of a model which gives its own without a level keep their
    # other values
    tolerances = SimulatorSettings(absolute_tolerance=1e-9)
    model = PKPDModel(**keywords(simulator=tolerances))
    assert model.simulator.jit_optimization == DEFAULT_JIT_OPTIMIZATION
    assert model.simulator.absolute_tolerance == 1e-9
    # the level of a model is its own, also in a copy
    named = PKPDModel(**keywords(simulator=SimulatorSettings(jit_optimization="O3")))
    assert named.simulator.jit_optimization == "O3"
    assert named.replace(save_results=True).simulator.jit_optimization == "O3"
    # the level of the environment of the process is the level of sbmlsim, which
    # resolves it when it loads a model, so the model names none
    monkeypatch.setenv(JIT_OPTIMIZATION_VARIABLE, "O2")
    assert PKPDModel(**keywords()).simulator.jit_optimization is None
    assert (
        PKPDModel(**keywords(simulator=tolerances)).simulator.jit_optimization is None
    )
    assert named.simulator.jit_optimization == "O3"


def test_the_level_variable_is_the_one_of_sbmlsim() -> None:
    """pkpdsim names the variable itself, so that a model does not import roadrunner."""
    from sbmlsim.model.model_roadrunner import JIT_OPTIMIZATION_VARIABLE as variable

    assert variable == JIT_OPTIMIZATION_VARIABLE


def test_a_model_without_credibility() -> None:
    with pytest.raises(ValueError, match="toy"):
        MODEL.load_credibility()


def test_a_credibility_which_is_no_declaration() -> None:
    model = MODEL.replace(credibility="pkpdsim_toy:MODEL")
    with pytest.raises(TypeError, match="not a Credibility"):
        model.load_credibility()
