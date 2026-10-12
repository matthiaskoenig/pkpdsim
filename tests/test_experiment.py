"""The base of the experiments of a model, on the toy model."""

from pathlib import Path
from types import MappingProxyType
from typing import ClassVar

import pandas as pd
import pytest
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.model import AbstractModel
from sbmlsim.plot import Figure
from sbmlsim.simulation import Simulation
from sbmlsim.units import ureg

from pkpdsim import impairment, pk
from pkpdsim.experiment import (
    PKPDSimulationExperiment,
    define_unit,
    rename_units,
    to_sid,
)
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


@pytest.mark.parametrize(
    ("text", "sid"),
    [
        ("task_po", "task_po"),
        ("po 2.5 mg", "po_2p5_mg"),
        ("dset_torasemide-M5_po_2.5", "dset_torasemide_M5_po_2p5"),
        ("tor_iv_chronic renal disease", "tor_iv_chronic_renal_disease"),
        ("Fig1_CIM400, ENA10_B1", "Fig1_CIM400__ENA10_B1"),
        ("µg/ml", "_g_ml"),
        ("1.5.2", "1p5p2"),
        ("a.b", "a_b"),
        ("2.", "2_"),
        (".5", "_5"),
        ("a b", "a_b"),
        ("a‐b", "a_b"),
        ("(a+b)*c&d", "_a_b__c_d"),
        ("", ""),
    ],
)
def test_a_text_becomes_an_sid(text: str, sid: str) -> None:
    assert to_sid(text) == sid
    assert len(to_sid(text)) == len(text)


def test_a_part_of_a_key_becomes_the_same_part_of_its_sid() -> None:
    assert f"task_po_{to_sid('2.5')}" == to_sid("task_po_2.5")


def test_the_experiment_gives_the_sid_of_a_text() -> None:
    experiment = initialized(ToyStudy)
    assert experiment.to_sid("po 2.5") == "po_2p5"
    assert experiment.to_sid("po 2.5") == "po_2p5"
    assert experiment.to_sid("task_decay") == "task_decay"


def test_two_texts_with_one_sid_are_an_error() -> None:
    experiment = initialized(ToyStudy)
    experiment.to_sid("po 2.5")
    with pytest.raises(ValueError, match=r"'po 2.5' and 'po-2.5'.*'po_2p5'.*ToyStudy"):
        experiment.to_sid("po-2.5")


def test_every_experiment_has_its_own_sids() -> None:
    initialized(ToyStudy).to_sid("a b")
    assert initialized(ToyStudy).to_sid("a-b") == "a_b"


def test_a_unit_is_defined_once() -> None:
    assert (
        define_unit(ureg, "IU_pkpdsim_test_once = 0.5 * mg") == "IU_pkpdsim_test_once"
    )
    assert define_unit(ureg, "IU_pkpdsim_test_once = 0.5*mg") == "IU_pkpdsim_test_once"
    assert ureg.Quantity(2, "µIU_pkpdsim_test_once").to(
        "ng"
    ).magnitude == pytest.approx(1.0)


def test_a_unit_with_another_value_is_an_error() -> None:
    define_unit(ureg, "IU_pkpdsim_test_other = 0.5 * mg")
    with pytest.raises(ValueError, match="IU_pkpdsim_test_other"):
        define_unit(ureg, "IU_pkpdsim_test_other = 0.6 * mg")


def test_a_count_is_a_unit() -> None:
    define_unit(ureg, "IU_pkpdsim_test_count = count")
    assert ureg.Quantity(3, "IU_pkpdsim_test_count/ml").to(
        "1/ml"
    ).magnitude == pytest.approx(3.0)


@pytest.mark.parametrize("definition", ["IU_x", "IU_x =", "= 0.5 * mg", "1IU = mg"])
def test_a_definition_needs_a_name_and_an_expression(definition: str) -> None:
    with pytest.raises(ValueError, match="<name> = <expression>"):
        define_unit(ureg, definition)


@pytest.mark.parametrize(
    ("unit", "renamed"),
    [
        ("IU/ml", "IU_x/ml"),
        ("IU/mL", "IU_x/mL"),
        ("µIU/ml", "µIU_x/ml"),
        ("μIU/ml", "μIU_x/ml"),
        ("mIU/l", "mIU_x/l"),
        ("IU", "IU_x"),
        ("IU/(ml*hr)", "IU_x/(ml*hr)"),
        ("IU_per_ml", "IU_per_ml"),
        ("2IU", "2IU"),
        ("mg/IUx", "mg/IUx"),
        ("ng/ml", "ng/ml"),
    ],
)
def test_units_are_renamed_where_they_are_units(unit: str, renamed: str) -> None:
    df = pd.DataFrame({"unit": [unit], "mean_unit": [unit], "label": ["IU"]})
    assert rename_units(df, {"IU": "IU_x"}) is df
    assert df["unit"].tolist() == [renamed]
    assert df["mean_unit"].tolist() == [renamed]
    assert df["label"].tolist() == ["IU"]


def test_a_missing_unit_stays_missing() -> None:
    df = pd.DataFrame({"unit": ["IU/ml", None]})
    rename_units(df, {"IU": "IU_x"})
    assert df["unit"].tolist()[0] == "IU_x/ml"
    assert pd.isna(df["unit"].tolist()[1])


class ToyIUStudy(ToyStudy):
    """A study whose data are in an international unit."""

    data_units = MappingProxyType({"IU": "IU_pkpdsim_test_toy = 0.5 * mg"})


def test_the_units_of_the_data_are_defined_and_renamed(tmp_path: Path) -> None:
    write_study(tmp_path, "ToyIUStudy_Fig1", "time\tA\tunit\n0\t1.0\tµIU/ml\n")
    experiment = initialized(ToyIUStudy, tmp_path)
    assert experiment.load_dataframe("Fig1")["unit"].tolist() == [
        "µIU_pkpdsim_test_toy/ml"
    ]
    quantity = experiment.ureg.Quantity(2, "µIU_pkpdsim_test_toy")
    assert quantity.to("ng").magnitude == pytest.approx(1.0)


def test_the_units_of_the_data_by_substance_are_renamed(tmp_path: Path) -> None:
    write_study(tmp_path, "ToyIUStudy_Fig2", "substance\tA\tunit\nx\t1.0\tIU/ml\n")
    frames = initialized(ToyIUStudy, tmp_path).load_dataframes_by_substance("Fig2")
    assert frames["x"]["unit"].tolist() == ["IU_pkpdsim_test_toy/ml"]


def test_missing_values_by_column(tmp_path: Path) -> None:
    write_study(tmp_path, "ToyStudy_Fig3", "time\ttime_unit\tlabel\nNR\tNR\tNR\n")
    df = initialized(ToyStudy, tmp_path).load_dataframe(
        "Fig3", na_values={"time": ["NR"], "time_unit": ["NR"]}
    )
    assert df["time"].isna().all()
    assert df["time_unit"].isna().all()
    assert df["label"].tolist() == ["NR"]


class CountedStudy(ToyStudy):
    """A study which counts the calls of its simulations."""

    calls: ClassVar[list[str]] = []

    def simulations(self) -> dict[str, Simulation]:
        """Count the call."""
        CountedStudy.calls.append(self.sid)
        return super().simulations()


def test_the_simulations_are_created_once() -> None:
    CountedStudy.calls.clear()
    experiment = initialized(CountedStudy)
    assert CountedStudy.calls == ["CountedStudy"]
    assert set(experiment._tasks) == {"task_decay"}


class DigitKeyStudy(ToyStudy):
    """A study whose simulation is keyed by a text which starts with a digit."""

    def simulations(self) -> dict[str, Simulation]:
        """Key the simulation by the SId of a dose."""
        return {self.to_sid("2.5mg"): Simulation(end=10, steps=50)}


def test_a_whole_key_which_starts_with_a_digit_fails_the_initialization() -> None:
    # the behavior of sbmlsim, which checks every key when it initializes
    experiment = DigitKeyStudy(base_path=TOY_PATH, data_path=DATA_PATHS)
    with pytest.raises(
        ValueError, match=r"must be valid SIds.*_simulations key '2p5mg'"
    ):
        experiment.initialize()
