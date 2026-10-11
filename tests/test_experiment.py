"""The base of the experiments of a model, on the toy model."""

from pathlib import Path
from typing import ClassVar

import pandas as pd
import pytest
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.model import AbstractModel
from sbmlsim.plot import Figure

from pkpdsim import impairment, pk
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim_toy import DATA_PATHS, MODEL, TOY_PATH
from pkpdsim_toy.experiments.decay import ToyDecay, ToyExperiment


def initialized(
    cls: type[PKPDSimulationExperiment],
    data_path: Path | list[Path] | None = DATA_PATHS,
) -> PKPDSimulationExperiment:
    runner = ExperimentRunner(
        experiment_classes=[cls],
        base_path=TOY_PATH,
        data_path=data_path,
        simulator=None,
    )
    experiment = runner.experiments[cls.__name__]
    assert isinstance(experiment, PKPDSimulationExperiment)
    return experiment


def test_the_base_simulates_the_model_of_the_model_object() -> None:
    model = initialized(ToyDecay).models()["model"]
    assert isinstance(model, AbstractModel)
    assert Path(model.source.source) == MODEL.model_path
    assert model.language_type == AbstractModel.LanguageType.SBML
    assert model.changes == {}


def test_every_simulation_has_a_task_of_the_model() -> None:
    tasks = initialized(ToyDecay).tasks()
    assert list(tasks) == ["task_decay"]
    assert tasks["task_decay"].model_id == "model"
    assert tasks["task_decay"].simulation_id == "decay"


def test_the_selections_are_added_to_every_task() -> None:
    experiment = initialized(ToyDecay)
    indices = {(d.task_id, d.selection) for d in experiment._data.values()}
    assert {("task_decay", "time"), ("task_decay", "[A]")} <= indices


def test_there_are_no_default_changes() -> None:
    assert initialized(ToyDecay).default_changes() == {}


def test_the_defaults_are_the_presets() -> None:
    assert ToyExperiment.renal_map is impairment.RENAL_FUNCTION
    assert ToyExperiment.renal_colors is impairment.RENAL_COLORS
    assert ToyExperiment.cirrhosis_map is impairment.CIRRHOSIS
    assert ToyExperiment.cirrhosis_colors is impairment.CIRRHOSIS_COLORS
    assert ToyExperiment.cardiac_map is impairment.CARDIAC
    assert ToyExperiment.cardiac_colors is impairment.CARDIAC_COLORS
    assert ToyExperiment.pk_labels is pk.PK_LABELS
    assert ToyExperiment.pk_units is pk.PK_UNITS


class ToyStudy(ToyExperiment):
    """An experiment named like a study, whose data is in TSV files."""


def write_study(directory: Path, name: str, content: str) -> None:
    (directory / name.split("_")[0]).mkdir(exist_ok=True)
    (directory / name.split("_")[0] / f".{name}.tsv").write_text(content)


def test_a_figure_of_the_study_is_loaded(tmp_path: Path) -> None:
    write_study(tmp_path, "ToyStudy_Fig1", "time\tA\n0\t1.0\n1\t2.0\n")
    df = initialized(ToyStudy, tmp_path).load_dataframe("Fig1")
    assert list(df.columns) == ["time", "A"]
    assert df["A"].tolist() == [1.0, 2.0]


def test_the_missing_values_are_passed_on(tmp_path: Path) -> None:
    write_study(tmp_path, "ToyStudy_Fig1", "time\tA\n0\t1.0\n1\tna\n")
    df = initialized(ToyStudy, tmp_path).load_dataframe("Fig1", na_values=["na"])
    assert pd.isna(df["A"].tolist()[1])


def test_the_study_can_be_given_by_position(tmp_path: Path) -> None:
    write_study(tmp_path, "Other2000_Tab2", "time\tA\n0\t2.0\n")
    experiment = initialized(ToyStudy, tmp_path)
    assert experiment.load_dataframe("Tab2", "Other2000")["A"].tolist() == [2.0]
    assert experiment.load_dataframe("Tab2", study="Other2000")["A"].tolist() == [2.0]


def test_the_data_is_split_by_substance(tmp_path: Path) -> None:
    write_study(tmp_path, "ToyStudy_Fig2", "substance\ttime\tA\nx\t0\t1.0\ny\t0\tna\n")
    frames = initialized(ToyStudy, tmp_path).load_dataframes_by_substance("Fig2")
    assert set(frames) == {"x", "y"}
    assert pd.isna(frames["y"]["A"].tolist()[0])


def test_an_experiment_without_data_paths_cannot_load() -> None:
    experiment = initialized(ToyStudy, None)
    with pytest.raises(ValueError, match="No 'data_path' set on 'ToyStudy'"):
        experiment.load_dataframe("Fig1")
    with pytest.raises(ValueError, match="No 'data_path' set on 'ToyStudy'"):
        experiment.load_dataframes_by_substance("Fig1")


def test_the_figures_are_keyed_by_their_id() -> None:
    experiment = initialized(ToyDecay)
    figure = experiment.figures()["fig_decay"]
    assert PKPDSimulationExperiment.figures_by_sid(figure) == {"fig_decay": figure}
    anonymous = Figure(experiment=experiment, sid="anonymous")
    anonymous.sid = None
    with pytest.raises(ValueError, match="Figure without an id"):
        PKPDSimulationExperiment.figures_by_sid(anonymous)


class Unbound(PKPDSimulationExperiment):
    """An experiment whose base sets no model object."""


def test_an_experiment_without_a_model_object_is_an_error() -> None:
    with pytest.raises(TypeError, match="Unbound has no model object"):
        Unbound(base_path=TOY_PATH)


class Narrowed(ToyExperiment):
    """A base which types its labels as a dict, as the models do."""

    labels: ClassVar[dict[str, str]] = {"time": "time"}


def test_a_model_types_its_conventions_itself() -> None:
    assert Narrowed.labels == {"time": "time"}
