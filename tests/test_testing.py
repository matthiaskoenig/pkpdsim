"""The tests of a model as functions, run on the toy model."""

import dataclasses
from pathlib import Path

import pytest
from sbmlsim.plot import Figure

from pkpdsim import testing
from pkpdsim_toy import MODEL
from pkpdsim_toy.experiments.decay import ToyDecay, ToyDecaySlow

PARTIAL = dataclasses.replace(MODEL, experiments=".simulations:PARTIAL")


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
    testing.check_registry(PARTIAL, not_registered={"ToyDecaySlow"})


def test_a_stale_exception_fails() -> None:
    with pytest.raises(AssertionError, match="listed as not registered"):
        testing.check_registry(MODEL, not_registered={"ToyDecay"})


def test_a_module_which_is_not_importable_is_skipped() -> None:
    modules, _ = testing.experiment_modules(MODEL, not_importable={"experiments.decay"})
    assert [module.__name__ for module in modules] == ["pkpdsim_toy.experiments"]


def test_the_experiments_initialize() -> None:
    assert set(testing.initialize_experiments(MODEL)) == {"ToyDecay", "ToyDecaySlow"}


@pytest.mark.parametrize("one_at_a_time", [False, True])
def test_all_experiments_run(tmp_path: Path, one_at_a_time: bool) -> None:
    testing.run_all_experiments(MODEL, tmp_path, one_at_a_time=one_at_a_time)
    if one_at_a_time:
        assert (tmp_path / "ToyDecay").is_dir()
        assert (tmp_path / "ToyDecaySlow").is_dir()
    else:
        assert (tmp_path / "all" / "_figures").is_dir()


def test_the_fits_start(tmp_path: Path) -> None:
    runs = testing.start_fits(MODEL, tmp_path)
    assert set(runs) == {"DECAY"}


def test_the_figure_settings_are_restored() -> None:
    before = Figure.legend_fontsize
    with testing.restored_figure_settings():
        Figure.legend_fontsize = 3
    assert Figure.legend_fontsize == before


def test_a_model_of_its_package_passes_the_check() -> None:
    testing.check_model(MODEL)


def test_a_broken_model_names_every_problem() -> None:
    broken = dataclasses.replace(
        MODEL,
        path=MODEL.path.parent,
        factory=".factory:missing",
        experiments=".simulations:not_a_registry",
        fit_definitions=".fitting:EMPTY",
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
    model = dataclasses.replace(
        MODEL,
        experiments=".simulations:WITHOUT_ALL",
        fit_definitions=".simulations:NO_FITS",
    )
    with pytest.raises(AssertionError) as info:
        testing.check_model(model)
    assert "'all'" in str(info.value)
    assert "non-empty mapping" in str(info.value)


def test_a_package_which_does_not_import_fails() -> None:
    model = dataclasses.replace(MODEL, package="pkpdsim_missing")
    with pytest.raises(AssertionError, match="does not import"):
        testing.check_model(model)
