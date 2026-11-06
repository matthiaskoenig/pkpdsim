"""The tests of a model as functions, run on the toy model."""

from pathlib import Path

import pytest
from pydantic import ValidationError
from sbmlsim.fit.cli import FitDefinition
from sbmlsim.project import CheckSettings, SimulatorSettings

from pkpdsim import testing
from pkpdsim_toy import MODEL
from pkpdsim_toy.broken import ToyCannotInitialize, ToyFigureFails, ToyRefusedUnit
from pkpdsim_toy.experiments.decay import ToyDecay, ToyDecaySlow

PARTIAL = MODEL.replace(experiments=".simulations:PARTIAL")
BROKEN = MODEL.replace(experiments=".broken:BROKEN")
KNOWN = PARTIAL.replace(checks=CheckSettings(not_registered={"ToyDecaySlow": "a test"}))


def test_the_modules_of_the_experiments() -> None:
    modules, errors = testing.experiment_modules(MODEL)
    assert errors == {}
    assert {module.__name__ for module in modules} == {
        "pkpdsim_toy.experiments",
        "pkpdsim_toy.experiments.decay",
    }


def test_the_defined_experiments_leave_out_the_base() -> None:
    modules, _ = testing.experiment_modules(MODEL)
    assert testing.defined_experiments(modules) == {ToyDecay, ToyDecaySlow}


def test_a_complete_registry_passes() -> None:
    testing.check_registry(MODEL)


def test_an_unregistered_experiment_fails() -> None:
    with pytest.raises(AssertionError, match="ToyDecaySlow"):
        testing.check_registry(PARTIAL)


def test_a_known_exception_passes() -> None:
    testing.check_registry(KNOWN)


def test_a_stale_exception_fails() -> None:
    stale = MODEL.replace(checks=CheckSettings(not_registered={"ToyDecay": "a test"}))
    with pytest.raises(AssertionError, match="listed as not registered"):
        testing.check_registry(stale)


def test_a_listed_class_which_is_not_defined_fails() -> None:
    gone = MODEL.replace(checks=CheckSettings(not_registered={"ToyGone": "a test"}))
    with pytest.raises(AssertionError, match=r"not defined.*ToyGone"):
        testing.check_registry(gone)


def test_a_module_which_is_not_importable_is_skipped() -> None:
    skipping = MODEL.replace(
        checks=CheckSettings(not_importable={"experiments.decay": "a test"})
    )
    modules, _ = testing.experiment_modules(skipping)
    assert [module.__name__ for module in modules] == ["pkpdsim_toy.experiments"]


def test_a_listed_module_which_does_not_exist_fails() -> None:
    gone = MODEL.replace(
        checks=CheckSettings(not_importable={"experiments.gone": "a test"})
    )
    with pytest.raises(AssertionError, match=r"no module.*experiments\.gone"):
        testing.check_registry(gone)


def test_the_experiments_initialize() -> None:
    assert set(testing.initialize_experiments(MODEL)) == {"ToyDecay", "ToyDecaySlow"}


@pytest.mark.parametrize("one_at_a_time", [None, "a test"])
def test_all_experiments_run(tmp_path: Path, one_at_a_time: str | None) -> None:
    model = MODEL.replace(checks=CheckSettings(one_at_a_time=one_at_a_time))
    testing.run_all_experiments(model, tmp_path)
    if one_at_a_time:
        assert (tmp_path / "ToyDecay").is_dir()
        assert (tmp_path / "ToyDecaySlow").is_dir()
    else:
        assert (tmp_path / "all" / "_figures").is_dir()


def test_the_fits_start(tmp_path: Path) -> None:
    runs = testing.start_fits(MODEL, tmp_path)
    assert set(runs) == {"DECAY"}


def test_a_model_of_its_package_passes_the_check() -> None:
    testing.check_model(MODEL)


def test_a_broken_model_names_every_problem() -> None:
    broken = MODEL.replace(
        root=MODEL.path.parent,
        factory=".factory:missing",
        experiments=".simulations:not_a_registry",
        fits=".fitting:EMPTY",
        extra_arguments=".extras:missing",
    )
    with pytest.raises(AssertionError) as info:
        testing.check_model(broken)
    message = str(info.value)
    for expected in [
        "MODEL of 'pkpdsim_toy' is not this model",
        "is not the directory",
        "factory '.factory:missing' does not resolve",
        "experiments '.simulations:not_a_registry' does not resolve",
        "fit_definitions '.fitting:EMPTY' does not resolve",
        "extra_arguments '.extras:missing' does not resolve",
    ]:
        assert expected in message, expected


def test_a_registry_without_all_and_empty_fits_fail() -> None:
    model = MODEL.replace(
        experiments=".simulations:WITHOUT_ALL",
        fits=".simulations:NO_FITS",
    )
    with pytest.raises(AssertionError) as info:
        testing.check_model(model)
    assert "'all'" in str(info.value)
    assert "non-empty mapping" in str(info.value)


def test_a_package_which_does_not_import_fails() -> None:
    model = MODEL.replace(package="pkpdsim_missing")
    with pytest.raises(AssertionError, match="does not import"):
        testing.check_model(model)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"name": "other"}, "is not the name of the package"),
        (
            {
                "variants": {"male": Path("/elsewhere/model.xml")},
                "factory_dir": MODEL.factory_dir,
            },
            "model_path",
        ),
        ({"archive_path": Path("/elsewhere/model.omex")}, "archive_path"),
        ({"data_path": (Path("/does/not/exist"),)}, "does not exist"),
    ],
)
def test_check_model_finds_a_slip(changes: dict[str, object], message: str) -> None:
    with pytest.raises(AssertionError, match=message):
        testing.check_model(MODEL.replace(**changes))


def test_the_factory_writes_below_the_model() -> None:
    with pytest.raises(ValidationError, match="outside of the root"):
        MODEL.replace(factory_dir=Path("/elsewhere"))


def test_every_failure_of_the_initialization_is_named() -> None:
    with pytest.raises(AssertionError) as info:
        testing.initialize_experiments(BROKEN)
    message = str(info.value)
    assert "2 of 4 experiments cannot be initialized" in message
    assert (
        "ToyCannotInitialize: RuntimeError: the simulations of ToyCannotInitialize fail"
        in message
    )
    assert "ToyRefusedUnit: ValueError: The unit 'IU_toy_refused'" in message


def test_every_failure_of_a_run_is_named(tmp_path: Path) -> None:
    with pytest.raises(AssertionError) as info:
        testing.run_all_experiments(BROKEN, tmp_path)
    message = str(info.value)
    assert ToyCannotInitialize.__name__ in message
    assert ToyFigureFails.__name__ in message
    assert ToyRefusedUnit.__name__ in message
    assert "figures_mpl()" in message


def test_every_failure_of_a_run_one_at_a_time_is_named(tmp_path: Path) -> None:
    one_by_one = BROKEN.replace(
        checks=CheckSettings(one_at_a_time="a test of the failures")
    )
    with pytest.raises(AssertionError) as info:
        testing.run_all_experiments(one_by_one, tmp_path)
    message = str(info.value)
    assert ToyCannotInitialize.__name__ in message
    assert ToyFigureFails.__name__ in message
    assert ToyRefusedUnit.__name__ in message
    # the experiment after a failing one (ToyFigureFails) ran and wrote its report
    assert list((tmp_path / "ToyCannotInitialize").glob("**/*.html"))


def test_the_initialization_compiles_with_the_level_of_the_model() -> None:
    model = MODEL.replace(
        simulator=SimulatorSettings(
            absolute_tolerance=1e-10, relative_tolerance=1e-10, jit_optimization="O0"
        )
    )
    assert str(testing._runner(model).jit_optimization) == "O0"


def test_the_fits_start_at_the_level_of_the_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`start_fits` builds the fit problems at `jit_optimization` of the model."""
    seen: list[str] = []

    def recording_run_fit(definition: FitDefinition, **kwargs: object) -> dict:
        seen.append(str(definition.problem(opid="level").jit_optimization))
        return {"run": object()}

    monkeypatch.setattr(testing, "run_fit", recording_run_fit)
    simulator = SimulatorSettings(
        absolute_tolerance=1e-10, relative_tolerance=1e-10, jit_optimization="O0"
    )
    testing.start_fits(MODEL.replace(simulator=simulator), tmp_path)
    assert seen == ["O0"]
