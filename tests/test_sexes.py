"""The builds per sex of a model: the model object, the experiments and the commands.

The toy model per sex (`tests/pkpdsim_sexes`) has a male and a female build,
which differ in the initial concentration of A (1 and `FEMALE_A0` mM); the
toy model (`tests/pkpdsim_toy`) has the male build only. A task which runs the
wrong build starts at the wrong concentration.
"""

import dataclasses
from collections.abc import Mapping
from pathlib import Path
from typing import ClassVar, cast

import numpy as np
import pytest
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.model import AbstractModel
from sbmlsim.simulator import Simulator

from pkpdsim.cli import fit_command, run_command
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.sex import Sex
from pkpdsim.testing import check_model, initialize_experiments, run_all_experiments
from pkpdsim_sexes import FEMALE_A0, FEMALE_MODEL_PATH, MODEL, MODEL_PATH
from pkpdsim_sexes.experiments.decay import SexesGroups, SexesMen, SexesWomen
from pkpdsim_toy import MODEL as TOY
from pkpdsim_toy.experiments.decay import ToyDecay


def initialized(cls: type[PKPDSimulationExperiment]) -> PKPDSimulationExperiment:
    """Load the models of an experiment and initialize it, as a run does."""
    runner = ExperimentRunner(
        experiment_classes=[cls],
        base_path=cls.pkpd_model.path,
        data_path=list(cls.pkpd_model.data_paths) or None,
        simulator=Simulator(n_workers=1),
        on_error="raise",
    )
    (experiment,) = runner.experiments.values()
    assert isinstance(experiment, PKPDSimulationExperiment)
    experiment.run(runner.simulator, output_path=None)
    return experiment


def loaded_paths(experiment: PKPDSimulationExperiment) -> dict[str, Path]:
    """Get the file of every loaded model of an experiment by its key."""
    paths: dict[str, Path] = {}
    for key, model in experiment._models.items():
        assert model.source.path is not None
        paths[key] = model.source.path
    return paths


def task_models(experiment: PKPDSimulationExperiment) -> dict[str, str]:
    """Get the key of the model of every task of an experiment."""
    return {key: task.model_id for key, task in experiment._tasks.items()}


def start(experiment: PKPDSimulationExperiment, task: str) -> float:
    """Get the concentration of A at the start of a task [mM]."""
    return float(np.asarray(experiment.results[task].ds["[A]"]).ravel()[0])


# ----------------------------------------------------------------------------
# the model object
# ----------------------------------------------------------------------------
def test_a_model_has_a_build_per_sex() -> None:
    assert MODEL.model_paths == {Sex.MALE: MODEL_PATH, Sex.FEMALE: FEMALE_MODEL_PATH}
    assert list(MODEL.model_paths) == [Sex.MALE, Sex.FEMALE]
    assert MODEL.model_path == MODEL_PATH
    assert MODEL.model_path_of(Sex.MALE) == MODEL_PATH
    assert MODEL.model_path_of(Sex.FEMALE) == FEMALE_MODEL_PATH
    assert hash(MODEL) == hash(dataclasses.replace(MODEL))


def test_a_model_without_a_female_build_has_the_male_one() -> None:
    assert TOY.female_model_path is None
    assert TOY.model_paths == {Sex.MALE: TOY.model_path}
    with pytest.raises(ValueError, match="'toy' has no female build") as error:
        TOY.model_path_of(Sex.FEMALE)
    assert "female_model_path" in str(error.value)


def test_a_replaced_male_build_keeps_the_female_one(tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, model_path=tmp_path / "male.xml")
    assert model.model_paths == {
        Sex.MALE: tmp_path / "male.xml",
        Sex.FEMALE: FEMALE_MODEL_PATH,
    }


def test_the_female_build_is_another_file() -> None:
    with pytest.raises(ValueError, match="the file of the male build"):
        dataclasses.replace(MODEL, female_model_path=MODEL_PATH)


# ----------------------------------------------------------------------------
# the choice of the experiments
# ----------------------------------------------------------------------------
def test_an_experiment_runs_the_male_build_by_default() -> None:
    assert SexesMen.sex is Sex.MALE
    men = initialized(SexesMen)
    assert loaded_paths(men) == {"model": MODEL_PATH}
    assert task_models(men) == {"task_decay": "model"}
    assert start(men, "task_decay") == pytest.approx(1.0)


def test_a_female_experiment_runs_the_female_build() -> None:
    women = initialized(SexesWomen)
    # the female build is the model of the experiment, the male one is not loaded
    assert loaded_paths(women) == {"model": FEMALE_MODEL_PATH}
    assert task_models(women) == {"task_decay": "model"}
    assert women.sex_of("decay") is Sex.FEMALE
    assert start(women, "task_decay") == pytest.approx(FEMALE_A0)


def test_a_study_of_both_sexes_runs_a_build_per_simulation() -> None:
    groups = initialized(SexesGroups)
    assert loaded_paths(groups) == {
        "model": MODEL_PATH,
        "model_female": FEMALE_MODEL_PATH,
    }
    assert task_models(groups) == {"task_men": "model", "task_women": "model_female"}
    assert groups.model_id(Sex.MALE) == "model"
    assert groups.model_id(Sex.FEMALE) == "model_female"
    assert start(groups, "task_men") == pytest.approx(1.0)
    assert start(groups, "task_women") == pytest.approx(FEMALE_A0)


def test_a_study_of_women_with_a_group_of_men() -> None:
    class Women(SexesGroups):
        sex = Sex.FEMALE
        simulation_sexes: ClassVar[Mapping[str, Sex]] = {"men": Sex.MALE}

    women = initialized(Women)
    assert loaded_paths(women) == {
        "model": FEMALE_MODEL_PATH,
        "model_male": MODEL_PATH,
    }
    assert task_models(women) == {"task_men": "model_male", "task_women": "model"}
    assert start(women, "task_men") == pytest.approx(1.0)
    assert start(women, "task_women") == pytest.approx(FEMALE_A0)


def test_a_model_without_a_female_build_runs_as_before() -> None:
    toy = initialized(ToyDecay)
    assert loaded_paths(toy) == {"model": TOY.model_path}
    assert task_models(toy) == {"task_decay": "model"}
    # the models of an experiment of a model with one build, as before the builds
    # per sex
    models = cast("dict[str, AbstractModel]", toy.models())
    assert list(models) == ["model"]
    assert models["model"].source.path == TOY.model_path
    assert models["model"].changes == {}


def test_a_female_experiment_of_a_model_without_a_female_build_fails() -> None:
    class Women(ToyDecay):
        sex = Sex.FEMALE

    with pytest.raises(ValueError, match="'toy' has no female build") as error:
        initialized(Women)
    assert "Women" in str(error.value)


def test_the_sex_of_a_simulation_which_the_experiment_has_not_fails() -> None:
    class Typo(SexesGroups):
        simulation_sexes: ClassVar[Mapping[str, Sex]] = {"woman": Sex.FEMALE}

    with pytest.raises(ValueError, match="'woman'") as error:
        initialized(Typo)
    assert "simulation_sexes" in str(error.value)


def test_the_sex_of_an_experiment_is_a_sex() -> None:
    class Text(SexesMen):
        sex = "female"  # ty: ignore[invalid-assignment]

    class Woman(SexesMen):
        sex = "woman"  # ty: ignore[invalid-assignment]

    assert loaded_paths(initialized(Text)) == {"model": FEMALE_MODEL_PATH}
    with pytest.raises(ValueError, match="'woman' is not a valid Sex"):
        initialized(Woman)


# ----------------------------------------------------------------------------
# the commands and the tests of a model
# ----------------------------------------------------------------------------
def test_the_check_of_the_model_sees_both_builds(tmp_path: Path) -> None:
    check_model(MODEL)
    outside = dataclasses.replace(MODEL, female_model_path=tmp_path / "female.xml")
    with pytest.raises(AssertionError, match="female_model_path"):
        check_model(outside)


def test_the_tests_of_a_model_initialize_and_run_both_builds(tmp_path: Path) -> None:
    experiments = initialize_experiments(MODEL)
    assert {
        sid: sorted(loaded_paths(cast("PKPDSimulationExperiment", e)).values())
        for sid, e in experiments.items()
    } == {
        "SexesMen": [MODEL_PATH],
        "SexesWomen": [FEMALE_MODEL_PATH],
        "SexesGroups": sorted([MODEL_PATH, FEMALE_MODEL_PATH]),
    }
    run_all_experiments(MODEL, tmp_path)


def test_run_simulation_runs_the_experiments_of_both_builds(tmp_path: Path) -> None:
    run_command(MODEL, ["-a", "simulate", "-e", "all", "--results-dir", str(tmp_path)])
    figures = tmp_path / "simulation" / "all" / "_figures"
    for name in ("SexesMen", "SexesWomen", "SexesGroups"):
        assert list(figures.glob(f"{name}_*.svg")), name


def test_a_fit_applies_one_parameter_set_to_both_builds(tmp_path: Path) -> None:
    """The data of both sexes decay at 0.2/min from their own start."""
    runs = fit_command(
        MODEL,
        ["--runs", "1", "--no-mapping-figures", "--output_dir", str(tmp_path)],
    )
    (run,) = runs.values()
    (k,) = run.result.xopt
    assert k == pytest.approx(0.2, rel=1e-2)
