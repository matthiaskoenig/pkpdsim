"""The base of the experiments of a model, on the toy model."""

from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import ClassVar, NamedTuple

import pandas as pd
import pytest
from sbmlsim.data import DataSet, rename_units
from sbmlsim.experiment import ExperimentRunner
from sbmlsim.model import AbstractModel
from sbmlsim.plot import Figure, Plot
from sbmlsim.simulation import Simulation
from sbmlsim.units import Quantity, define_unit, ureg

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
    (directory / name.split("_")[0] / f".{name}.tsv").write_text(
        content, encoding="utf-8"
    )


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


#: a table of PK-DB with geometric means in two units and a mean
GMEAN_TABLE = (
    "label\ttime\ttime_unit\tmean\tsd\tgmean\tgsd\tunit\n"
    "a\t1\thr\tNA\tNA\t2.0\t1.5\tng/ml\n"
    "a\t2\thr\tNA\tNA\t3.0\t1.4\tµg/ml\n"
    "b\t1\thr\t5.0\t1.0\tNA\tNA\tmg\n"
)


def test_the_geometric_means_get_their_unit(tmp_path: Path) -> None:
    """The table is read as it is; its `DataSet` gives `gmean` the unit of its rows."""
    write_study(tmp_path, "ToyStudy_Fig3", GMEAN_TABLE)
    experiment = initialized(ToyStudy, tmp_path)
    df = experiment.load_dataframe("Fig3")
    assert "gmean_unit" not in df.columns
    dset = DataSet.from_df(df[df["label"] == "a"], experiment.ureg)
    assert dset.uinfo["gmean"] == "ng/ml"
    assert dset["gmean"].tolist() == pytest.approx([2.0, 3000.0])
    # the geometric SD is a ratio
    assert dset.uinfo["gsd"] == "dimensionless"
    assert dset["gsd"].tolist() == pytest.approx([1.5, 1.4])


def test_the_geometric_means_by_substance_get_their_unit(tmp_path: Path) -> None:
    write_study(tmp_path, "ToyStudy_Fig4", "substance\tgmean\tunit\nx\t1.0\tng/ml\n")
    experiment = initialized(ToyStudy, tmp_path)
    frames = experiment.load_dataframes_by_substance("Fig4")
    dset = DataSet.from_df(frames["x"], experiment.ureg)
    assert dset.uinfo["gmean"] == "ng/ml"


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


class MolecularWeights(NamedTuple):
    """The molecular weights of a model, by substance."""

    A: Quantity


class Narrowed(ToyExperiment):
    """A base which types its conventions itself, as the models do."""

    Mr: ClassVar[MolecularWeights] = MolecularWeights(A=ureg.Quantity(100, "g/mole"))
    labels: ClassVar[Mapping[str, str]] = {"time": "time"}


def test_a_model_types_its_conventions_itself() -> None:
    assert Narrowed.Mr.A.m_as("g/mole") == 100
    assert Narrowed.labels == {"time": "time"}


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


def test_the_units_of_the_data_are_defined_when_a_loader_is_called_first(
    tmp_path: Path,
) -> None:
    """A loader called outside of `initialize` defines the units of the data."""

    class ToyIUStudyAlone(ToyStudy):
        data_units = MappingProxyType({"IU": "IU_pkpdsim_test_alone = 0.25 * mg"})

    write_study(tmp_path, "ToyIUStudyAlone_Fig1", "time\tA\tunit\n0\t1.0\tIU/ml\n")
    experiment = ToyIUStudyAlone(base_path=TOY_PATH, data_path=tmp_path)
    for _ in range(2):
        assert experiment.load_dataframe("Fig1")["unit"].tolist() == [
            "IU_pkpdsim_test_alone/ml"
        ]
    quantity = experiment.ureg.Quantity(4, "IU_pkpdsim_test_alone")
    assert quantity.to("mg").magnitude == pytest.approx(1.0)


def test_the_names_of_the_units_of_the_data(tmp_path: Path) -> None:
    """A table which a model reads itself is renamed with the names of its experiment."""
    experiment = initialized(ToyIUStudy, tmp_path)
    names = experiment.data_unit_names()
    assert dict(names) == {"IU": "IU_pkpdsim_test_toy"}
    df = rename_units(pd.DataFrame({"unit": ["IU/ml"]}), names)
    assert df["unit"].tolist() == ["IU_pkpdsim_test_toy/ml"]
    with pytest.raises(TypeError):
        names["IU"] = "IU_other"  # ty: ignore[invalid-assignment]
    assert dict(experiment.data_unit_names()) == {"IU": "IU_pkpdsim_test_toy"}
    assert dict(initialized(ToyStudy, tmp_path).data_unit_names()) == {}


def test_the_names_of_the_units_of_the_data_define_them(tmp_path: Path) -> None:
    """The names, asked for before `initialize`, define the units."""

    class ToyIUStudyNames(ToyStudy):
        data_units = MappingProxyType({"IU": "IU_pkpdsim_test_names = 0.125 * mg"})

    experiment = ToyIUStudyNames(base_path=TOY_PATH, data_path=tmp_path)
    assert dict(experiment.data_unit_names()) == {"IU": "IU_pkpdsim_test_names"}
    quantity = experiment.ureg.Quantity(8, "IU_pkpdsim_test_names")
    assert quantity.to("mg").magnitude == pytest.approx(1.0)


class ToyRefusedStudy(ToyStudy):
    """A study whose unit of the data is defined with another value."""

    data_units = MappingProxyType({"IU": "IU_pkpdsim_test_toy = 0.6 * mg"})


def test_a_refused_unit_of_the_data_fails_the_initialization(tmp_path: Path) -> None:
    define_unit("IU_pkpdsim_test_toy = 0.5 * mg", ureg)
    # the experiment is created, which a run does before its error handling
    experiment = ToyRefusedStudy(base_path=TOY_PATH, data_path=tmp_path)
    with pytest.raises(ValueError, match="IU_pkpdsim_test_toy") as info:
        experiment.initialize()
    assert "The experiment 'ToyRefusedStudy' cannot be initialized." in (
        info.value.__notes__
    )
    with pytest.raises(ValueError, match="IU_pkpdsim_test_toy"):
        experiment.load_dataframe("Fig1")


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


def test_rows_without_a_substance_are_in_no_frame(tmp_path: Path) -> None:
    """The rows without a substance (a bile flow) make no empty frame."""
    write_study(
        tmp_path, "ToyStudy_Fig5", "substance\ttime\tA\nx\t0\t1.0\nNA\t0\t2.0\n"
    )
    frames = initialized(ToyStudy, tmp_path).load_dataframes_by_substance("Fig5")
    assert list(frames) == ["x"]
    assert frames["x"]["A"].tolist() == [1.0]


def _plot() -> Plot:
    """A plot of a new figure."""
    return Figure(experiment=None, sid="fig", num_rows=1, num_cols=1).create_plots()[0]


def test_a_function_of_a_task_is_a_curve() -> None:
    """`add_function` draws a function of the results of a task over the time."""
    from sbmlsim.data import Data

    experiment = initialized(ToyDecay)
    assert experiment.across() is None
    plot = _plot()
    y = Data(
        index="twice_A", function="2 * A", variables={"A": Data("A", task="task_decay")}
    )
    experiment.add_function(plot, "task_decay", y, label="2A", color="black")
    experiment.add_function(plot, "task_decay", y, None, "tab:blue", linestyle="--")
    experiment.add_function(plot, "task_decay", y, "plain", "red", linestyle=None)
    first, second, _ = plot.curves
    assert (first.name, first.x.index, first.y.index) == ("2A", "time", "twice_A")
    assert first.style is not None
    assert second.name is None
    styles = [c.kwargs for c in plot.curves]
    assert styles[0] == {"color": "black", "linestyle": "-", "linewidth": 2.0}
    assert styles[1] == {"color": "tab:blue", "linestyle": "--", "linewidth": 2.0}
    # without a linestyle the curve has the default of sbmlsim
    assert styles[2] == {"color": "red", "linewidth": 2.0}


def test_a_function_of_a_task_is_a_band_over_the_draws() -> None:
    """In a run over the draws of the fitted parameters the function is a band."""
    from sbmlsim import Q
    from sbmlsim.data import Data
    from sbmlsim.simulation import Dimension

    from pkpdsim.uncertainty import UncertaintyRun, UncertaintySettings

    experiment = initialized(ToyDecay)
    draws = Dimension("dim_draws", values={"k": Q([0.1, 0.2], "1/min")})
    experiment.settings["uncertainty"] = UncertaintyRun(UncertaintySettings(), draws)
    assert experiment.across() == "dim_draws"
    plot = _plot()
    y = Data(
        index="twice_A", function="2 * A", variables={"A": Data("A", task="task_decay")}
    )
    experiment.add_function(plot, "task_decay", y, label="2A", color="black")
    assert not plot.curves
    (band,) = plot.bands
    assert band.across == "dim_draws"
    assert band.quantiles == (0.05, 0.95)
