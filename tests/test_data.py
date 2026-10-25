import csv
import math
import re
import shutil
import warnings
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest
from pkdb.schemas.study import CanonicalStudy

from pkpdsim.data import (
    GROUP_COLUMNS,
    INDIVIDUAL_COLUMNS,
    INTERVENTION_COLUMNS,
    MEASUREMENT_COLUMNS,
    PKDB_DATA,
    SCATTER_COLUMNS,
    CollectionWarning,
    _label,
    _statistics,
    collect_tsv_files,
    database_differences,
    mean_of_individuals,
    prepare_study,
    sampling_times,
    studies_dir,
    study_format,
    study_tables,
    write_study_tables,
)


@pytest.fixture
def workdir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A repository as working directory, without the environment variable."""
    repository = tmp_path / "repository"
    repository.mkdir()
    monkeypatch.chdir(repository)
    monkeypatch.delenv(PKDB_DATA, raising=False)
    return repository


def test_studies_dir_from_environment(
    workdir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The variable wins over a checkout next to the working directory."""
    checkout = tmp_path / "checkouts" / "pkdb_data"
    (checkout / "studies").mkdir(parents=True)
    (tmp_path / "pkdb_data" / "studies").mkdir(parents=True)
    monkeypatch.setenv(PKDB_DATA, str(checkout))
    assert studies_dir() == (checkout / "studies").resolve()


def test_studies_dir_next_to_working_directory(workdir: Path) -> None:
    """Without the variable the checkout is `../pkdb_data`."""
    studies = workdir.parent / "pkdb_data" / "studies"
    studies.mkdir(parents=True)
    assert studies_dir() == studies.resolve()
    assert studies_dir().is_absolute()


def test_studies_dir_empty_environment(
    workdir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An empty variable is no variable."""
    studies = workdir.parent / "pkdb_data" / "studies"
    studies.mkdir(parents=True)
    monkeypatch.setenv(PKDB_DATA, "")
    assert studies_dir() == studies.resolve()


def test_studies_dir_missing(workdir: Path) -> None:
    """The error names the variable and the default checkout."""
    with pytest.raises(FileNotFoundError) as error:
        studies_dir()
    message = str(error.value)
    assert PKDB_DATA in message
    assert str((workdir.parent / "pkdb_data").resolve()) in message


def test_studies_dir_environment_missing(
    workdir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A variable which points nowhere is an error, also with a default checkout."""
    (workdir.parent / "pkdb_data" / "studies").mkdir(parents=True)
    missing = tmp_path / "missing"
    monkeypatch.setenv(PKDB_DATA, str(missing))
    with pytest.raises(FileNotFoundError) as error:
        studies_dir()
    message = str(error.value)
    assert f"{PKDB_DATA}='{missing}'" in message
    assert str((workdir.parent / "pkdb_data").resolve()) in message


@pytest.fixture
def source_dir(tmp_path: Path) -> Path:
    """The folder of the studies of a substance in pkdb_data."""
    source = tmp_path / "pkdb_data" / "studies" / "albuterol"
    files = {
        "Boulton1996a": [
            ".Boulton1996a_Fig1.tsv",
            ".Boulton1996a_Tab1.tsv",
            "Boulton1996a_Fig1.png",
            "Boulton1996a.pdf",
            "study.json",
        ],
        "Hindle1992": [".Hindle1992_Fig2.tsv", "Hindle1992_Fig2.tsv"],
        "Ward2000": ["Ward2000.pdf"],
    }
    for study, names in files.items():
        (source / study).mkdir(parents=True)
        for name in names:
            (source / study / name).write_text(f"{study}/{name}", encoding="utf-8")
    # a study of format 1 has no field `format`
    (source / "Boulton1996a" / "study.json").write_text(
        '{"sid": "Boulton1996a"}', encoding="utf-8"
    )
    (source / "README.md").write_text("not a study", encoding="utf-8")
    return source


def tree(directory: Path) -> set[str]:
    return {
        path.relative_to(directory).as_posix()
        for path in directory.rglob("*")
        if path.is_file()
    }


def test_collect_tsv_files(source_dir: Path, tmp_path: Path) -> None:
    """The hidden tsv files are copied, in the folders of their studies."""
    target_dir = tmp_path / "model" / "data"
    target_dir.parent.mkdir()
    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)
    assert tree(target_dir) == {
        "Boulton1996a/.Boulton1996a_Fig1.tsv",
        "Boulton1996a/.Boulton1996a_Tab1.tsv",
        "Hindle1992/.Hindle1992_Fig2.tsv",
    }
    # a study without data gets no folder
    assert not (target_dir / "Ward2000").exists()
    copy = target_dir / "Hindle1992" / ".Hindle1992_Fig2.tsv"
    assert copy.read_text(encoding="utf-8") == "Hindle1992/.Hindle1992_Fig2.tsv"


def test_collect_tsv_files_filter(source_dir: Path, tmp_path: Path) -> None:
    """Only the studies which pass the filter are collected."""
    target_dir = tmp_path / "data"
    collect_tsv_files(
        source_dir=source_dir,
        target_dir=target_dir,
        filter_study=lambda name: name.startswith("Hindle"),
    )
    assert tree(target_dir) == {"Hindle1992/.Hindle1992_Fig2.tsv"}


def test_collect_tsv_files_again(source_dir: Path, tmp_path: Path) -> None:
    """A second collection replaces the copies and keeps what is there."""
    target_dir = tmp_path / "data"
    (target_dir / "Boulton1996a").mkdir(parents=True)
    (target_dir / "Boulton1996a" / ".Boulton1996a_Fig1.tsv").write_text(
        "outdated", encoding="utf-8"
    )
    (target_dir / "collect_data.py").write_text("", encoding="utf-8")
    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)
    copy = target_dir / "Boulton1996a" / ".Boulton1996a_Fig1.tsv"
    assert copy.read_text(encoding="utf-8") == "Boulton1996a/.Boulton1996a_Fig1.tsv"
    assert (target_dir / "collect_data.py").exists()


def test_the_study_format(tmp_path: Path) -> None:
    """A study without a study.json or without a format is in format 1."""
    assert study_format(tmp_path) == 1
    (tmp_path / "study.json").write_text('{"sid": "A"}', encoding="utf-8")
    assert study_format(tmp_path) == 1
    (tmp_path / "study.json").write_text('{"format": 2}', encoding="utf-8")
    assert study_format(tmp_path) == 2


@pytest.mark.parametrize(
    "text", ["no json", "[2]", '{"format": "2"}', '{"format": true}']
)
def test_a_study_format_which_cannot_be_read(tmp_path: Path, text: str) -> None:
    """A study.json which is no JSON object with an integer format is an error."""
    (tmp_path / "study.json").write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=r"study\.json"):
        study_format(tmp_path)


#: a study of format 2, a copy of pkdb_data (see tests/data/pkdb_data/README.md)
HYDROCHLOROTHIAZIDE = (
    Path(__file__).parent / "data" / "pkdb_data" / "studies" / "hydrochlorothiazide"
)
PATEL1984 = HYDROCHLOROTHIAZIDE / "Patel1984"
#: a study of format 2 with scatters
NILSEN1989 = HYDROCHLOROTHIAZIDE / "Nilsen1989"


@pytest.fixture(scope="module")
def patel1984_tables() -> dict[str, pd.DataFrame]:
    """The tables of format 1 of the study of format 2."""
    return study_tables(PATEL1984)


def test_study_tables_sources(patel1984_tables: dict[str, pd.DataFrame]) -> None:
    """A table per source of the measurements, the subjects and the interventions."""
    assert list(patel1984_tables) == [
        "Fig1",
        "Fig2",
        "Tab2",
        "Tab3",
        "Tab4",
        "TabGroups",
        "TabIndividuals",
        "TabInterventions",
    ]
    for source in ("Fig1", "Tab2"):
        assert list(patel1984_tables[source].columns) == list(MEASUREMENT_COLUMNS)
    assert list(patel1984_tables["TabGroups"].columns) == list(GROUP_COLUMNS)
    assert list(patel1984_tables["TabIndividuals"].columns) == list(INDIVIDUAL_COLUMNS)
    assert list(patel1984_tables["TabInterventions"].columns) == list(
        INTERVENTION_COLUMNS
    )


def test_study_tables_measurements(patel1984_tables: dict[str, pd.DataFrame]) -> None:
    """The reported measurements with the columns of format 1, no normalized ones."""
    source = pd.read_csv(PATEL1984 / "timecourses_Fig1.tsv", sep="\t")
    df = patel1984_tables["Fig1"]
    assert len(df) == len(source)
    assert set(df["study"]) == {"Patel1984"}
    assert df["time"].tolist() == source["time"].tolist()
    assert df["mean"].tolist() == pytest.approx(source["mean"].tolist())
    assert df["unit"].tolist() == source["unit"].tolist()
    assert df["label"].tolist() == source["label"].tolist()
    assert df["intervention"].tolist() == source["interventions"].tolist()
    assert df["measurement_type"].tolist() == source["measurement"].tolist()
    # the count of a row without one is the count of its group
    assert set(df["group"]) == {"all"}
    assert set(df["count"]) == {12}


def test_study_tables_outputs(patel1984_tables: dict[str, pd.DataFrame]) -> None:
    """An output gets a label, its missing statistics are calculated."""
    df = patel1984_tables["Tab4"]
    assert df["label"].notna().all()
    assert df["label"].iloc[0] == (
        "cumulative_amount_change_chloride_urine_HCTZ25_sus_all"
    )
    with_sd = df["sd"].notna() & df["count"].notna()
    assert with_sd.any()
    expected = df.loc[with_sd, "sd"] / df.loc[with_sd, "count"] ** 0.5
    assert df.loc[with_sd, "se"].tolist() == pytest.approx(expected.tolist())


@pytest.mark.parametrize(
    ("parts", "expected"),
    [
        (("auc", "hctz", "plasma", "HCTZ25", "all"), "auc_hctz_plasma_HCTZ25_all"),
        (
            ("cumulative amount (change)", None, "urine", "", "g1"),
            "cumulative_amount_change_urine_g1",
        ),
    ],
)
def test_label(parts: tuple[str | None, ...], expected: str) -> None:
    """The label of an output is made of the columns which identify it."""
    assert _label(*parts) == expected


@pytest.mark.parametrize(
    ("statistics", "expected"),
    [
        ({"mean": 10.0, "error_bar": 12.5, "error_type": "sd"}, {"sd": 2.5}),
        ({"mean": 10.0, "error_bar": 9.0, "error_type": "se"}, {"se": 1.0}),
        ({"gmean": 4.0, "error_bar": 8.0, "error_type": "gsd"}, {"gsd": 2.0}),
        # a reported statistic is kept
        ({"mean": 10.0, "sd": 1.0, "error_bar": 12.5, "error_type": "sd"}, {"sd": 1.0}),
    ],
)
def test_error_bars(statistics: dict[str, Any], expected: dict[str, float]) -> None:
    """A digitized error bar gives the statistic it names, the size of the bar."""
    from pkdb.schemas.study import Statistics

    observation = SimpleNamespace(statistics=Statistics(**statistics))
    values = _statistics(observation)  # ty: ignore[invalid-argument-type]
    assert {key: values[key] for key in expected} == pytest.approx(expected)


def test_study_tables_subjects(patel1984_tables: dict[str, pd.DataFrame]) -> None:
    """The groups with their characteristics, the interventions with their doses."""
    groups = patel1984_tables["TabGroups"]
    assert set(groups["name"]) == {"all"}
    assert set(groups["group_count"]) == {12}
    assert {"species", "sex", "healthy"} <= set(groups["measurement_type"])
    interventions = patel1984_tables["TabInterventions"]
    source = pd.read_csv(PATEL1984 / "interventions.tsv", sep="\t")
    assert sorted(interventions["name"]) == sorted(source["name"])
    assert set(interventions["substance"]) == {"hydrochlorothiazide"}
    assert {"interval", "doses"} <= set(interventions.columns)
    assert patel1984_tables["TabIndividuals"].empty


def test_study_tables_scatters() -> None:
    """A point of a scatter is a row with the x and y columns of format 1."""
    tables = study_tables(NILSEN1989)
    source = pd.read_csv(NILSEN1989 / "scatters_Fig5.tsv", sep="\t")
    df = tables["Fig5"]
    assert list(df.columns) == list(SCATTER_COLUMNS)
    assert len(df) == len(source)
    assert df["label"].tolist() == source["name"].tolist()
    assert df["x_value"].tolist() == pytest.approx(source["x_mean"].tolist())
    assert df["y_value"].tolist() == pytest.approx(source["y_mean"].tolist())
    assert df["x_unit"].tolist() == source["x_unit"].tolist()
    assert df["y_measurement_type"].tolist() == source["y_measurement"].tolist()
    # the value also in the columns of sbmlsim, a point of the group counts the group
    assert df["x"].tolist() == df["x_value"].tolist()
    assert set(df["count"]) == {9}


def edit_table(path: Path, edit: Callable[[list[dict[str, str]]], None]) -> None:
    """Edit the rows of a table of format 2, `edit(rows)` on the rows as dicts."""
    with path.open(encoding="utf-8", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        columns = list(reader.fieldnames or [])
        rows = list(reader)
    edit(rows)
    lines = ["\t".join(columns), *("\t".join(row[c] for c in columns) for row in rows)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="")


@pytest.fixture(scope="module")
def edited_tables(tmp_path_factory: pytest.TempPathFactory) -> dict[str, pd.DataFrame]:
    """The tables of a copy of Patel1984 with the cases the ports of the models met.

    The cmax of Tab2 is a geometric mean with a gcv of 30 %, the first point of Fig1
    is below the limit of quantification (only its `max`), and Tab3 has a cmax
    of the same measurement, substance, tissue, interventions and subject as
    Tab2, i.e. the same generated label.
    """
    from pkdb.studyformat.formatter import format_folder

    study = tmp_path_factory.mktemp("studies") / "Patel1984"
    shutil.copytree(PATEL1984, study)

    def geometric_cmax(rows: list[dict[str, str]]) -> None:
        for row in rows:
            if row["measurement"] == "cmax":
                row.update(gmean=row["mean"], mean="", sd="", gcv="30")

    def below_loq(rows: list[dict[str, str]]) -> None:
        rows[0].update(mean="", sd="", max="5", comment="below the LOQ of 5 ng/ml")

    def second_cmax(rows: list[dict[str, str]]) -> None:
        row = dict(rows[0])
        row.update(measurement="cmax", mean="130", sd="40", unit="ng/ml", comment="")
        rows.insert(0, row)

    edit_table(study / "outputs_Tab2.tsv", geometric_cmax)
    edit_table(study / "timecourses_Fig1.tsv", below_loq)
    edit_table(study / "outputs_Tab3.tsv", second_cmax)
    format_folder(study)
    return study_tables(study)


def test_a_geometric_mean(edited_tables: dict[str, pd.DataFrame]) -> None:
    """A row with a geometric mean and no mean is a geometric mean, with its gsd."""
    df = edited_tables["Tab2"]
    cmax = df[df["measurement_type"] == "cmax"].iloc[0]
    assert cmax["calculation_type"] == "geometric mean"
    assert pd.isna(cmax["mean"])
    assert cmax["gmean"] == pytest.approx(134)
    assert cmax["gcv"] == pytest.approx(0.3)
    assert cmax["gsd"] == pytest.approx(math.exp(math.sqrt(math.log(1 + 0.3**2))))
    others = df[df["measurement_type"] != "cmax"]
    assert set(others["calculation_type"]) == {"sample mean"}


def test_a_value_below_the_loq(
    edited_tables: dict[str, pd.DataFrame], tmp_path: Path
) -> None:
    """A row without a central value has none, also in the written table."""
    df = edited_tables["Fig1"]
    first = df.iloc[0]
    assert first["time"] == 0.5
    assert first[["value", "mean", "median", "gmean"]].isna().all()
    assert first["max"] == 5
    assert first["comments"] == "below the LOQ of 5 ng/ml"
    assert df["mean"].iloc[1:].notna().all()
    path = tmp_path / "Fig1.tsv"
    df.to_csv(path, sep="\t", index=False, na_rep="NA")
    written = pd.read_csv(path, sep="\t", comment="#")
    assert written["mean"].isna().tolist() == df["mean"].isna().tolist()
    assert len(written.dropna(subset=["mean"])) == len(df) - 1


def test_a_generated_label_is_unique(edited_tables: dict[str, pd.DataFrame]) -> None:
    """A generated label of several sources gets the source, in every source."""
    label = "cmax_hydrochlorothiazide_plasma_HCTZ25_sus_all"
    tab2, tab3 = edited_tables["Tab2"], edited_tables["Tab3"]
    suspension = tab2[
        (tab2["measurement_type"] == "cmax") & (tab2["intervention"] == "HCTZ25_sus")
    ]
    assert suspension["label"].tolist() == [f"{label}_Tab2"]
    assert tab3.loc[tab3["measurement_type"] == "cmax", "label"].tolist() == [
        f"{label}_Tab3"
    ]
    # the labels which only one source has stay
    assert "cmax_hydrochlorothiazide_plasma_HCTZ25_tab_all" in set(tab2["label"])
    assert "auc_end_hydrochlorothiazide_plasma_HCTZ25_sus_all" in set(tab2["label"])
    assert "kel_hydrochlorothiazide_plasma_HCTZ25_sus_all" in set(tab3["label"])
    labels = [set(df["label"]) for df in (tab2, tab3, edited_tables["Fig1"])]
    assert not labels[0] & labels[1]
    assert not (labels[0] | labels[1]) & labels[2]


def test_the_interventions_in_the_order_of_format_2() -> None:
    """The interventions are the names of the cell in its natural order, comma separated."""
    from pkdb.studyformat.text import natural_key

    source = pd.read_csv(NILSEN1989 / "outputs_Tab1.tsv", sep="\t")
    df = study_tables(NILSEN1989)["Tab1"]
    assert df["intervention"].tolist() == source["interventions"].tolist()
    assert "CIL5,HCTZ25" in set(df["intervention"])
    for interventions in df["intervention"]:
        names = interventions.split(",")
        assert names == sorted(names, key=natural_key)
    # in the generated labels joined by `_`
    combined = df[df["intervention"] == "CIL5,HCTZ25"]
    assert all("_CIL5_HCTZ25_" in label for label in combined["label"])


def test_collect_a_study_of_format_2(source_dir: Path, tmp_path: Path) -> None:
    """A study of format 2 is written as the hidden tables of format 1."""
    shutil.copytree(PATEL1984, source_dir / "Patel1984")
    target_dir = tmp_path / "data"
    (target_dir / "Patel1984").mkdir(parents=True)
    # a source the study no longer has is removed, other files stay
    (target_dir / "Patel1984" / ".Patel1984_Tab9.tsv").write_text(
        "old", encoding="utf-8"
    )
    (target_dir / "Patel1984" / "notes.md").write_text("", encoding="utf-8")
    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)
    assert {path for path in tree(target_dir) if path.startswith("Patel1984/")} == {
        "Patel1984/notes.md",
        *(
            f"Patel1984/.Patel1984_{source}.tsv"
            for source in (
                "Fig1",
                "Fig2",
                "Tab2",
                "Tab3",
                "Tab4",
                "TabGroups",
                "TabIndividuals",
                "TabInterventions",
            )
        ),
    }
    # the format 1 studies are copied as before
    assert "Hindle1992/.Hindle1992_Fig2.tsv" in tree(target_dir)
    df = pd.read_csv(target_dir / "Patel1984" / ".Patel1984_Fig1.tsv", sep="\t")
    assert list(df.columns) == list(MEASUREMENT_COLUMNS)


def test_a_study_of_another_format_is_not_collected(
    source_dir: Path, tmp_path: Path
) -> None:
    """A study which is neither in format 1 nor 2 is an error, nothing is written."""
    study = source_dir / "Mendell2012"
    study.mkdir()
    (study / "study.json").write_text('{"format": 3}', encoding="utf-8")
    target_dir = tmp_path / "data"
    with pytest.raises(ValueError, match=r"Mendell2012 \(format 3\)"):
        collect_tsv_files(source_dir=source_dir, target_dir=target_dir)
    assert not target_dir.exists()


class FakeStudies:
    """The studies of a server, by their study id."""

    def __init__(self, studies: dict[str, CanonicalStudy]) -> None:
        self.studies = studies

    def get(self, sid: str) -> CanonicalStudy:
        from pkdb.errors import ClientError

        if sid not in self.studies:
            raise ClientError("not found", status_code=404)
        return self.studies[sid]


class FakeClient:
    """A client of a server with the given studies."""

    def __init__(self, studies: dict[str, CanonicalStudy]) -> None:
        self.studies = FakeStudies(studies)


def test_database_differences() -> None:
    """A study agrees with the database when the digests of its sources are equal."""
    local = prepare_study(PATEL1984)
    sid = "hydrochlorothiazide/Patel1984"
    client = FakeClient({sid: local})
    assert database_differences([PATEL1984], client=client) == []  # ty: ignore[invalid-argument-type]

    other = local.model_copy(update={"source_digest": "0" * 64})
    differences = database_differences(
        [PATEL1984],
        client=FakeClient({sid: other}),  # ty: ignore[invalid-argument-type]
    )
    assert [(d.sid, d.study) for d in differences] == [(sid, PATEL1984)]
    assert "another version" in differences[0].reason

    missing = database_differences([PATEL1984], client=FakeClient({}))  # ty: ignore[invalid-argument-type]
    assert [d.reason for d in missing] == ["not in the database"]


def collect_warnings(
    function: Callable[[], dict[str, pd.DataFrame]],
) -> tuple[dict[str, pd.DataFrame], list[str]]:
    """Call a function, with the messages of the warnings of the collection it gave."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = function()
    messages = [
        str(warning.message)
        for warning in caught
        if issubclass(warning.category, CollectionWarning)
    ]
    return result, messages


@pytest.fixture(scope="module")
def followups(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[Path, dict[str, pd.DataFrame], list[str]]:
    """The tables and warnings of a copy of Patel1984 with the cases of the model sync.

    Fig2 `amount_cumulative_hctz25_sus` is a curve of geometric means with the
    curator 0 at 0 hr and a point below the LOQ at 2 hr, Tab3 repeats the AUC
    of the suspension of Tab2, a group has no characteristics (no count of a
    characteristic) and a curve of Fig1 has a label which is no SId.
    """
    from pkdb.studyformat.formatter import format_folder

    study = tmp_path_factory.mktemp("followups") / "Patel1984"
    shutil.copytree(PATEL1984, study)

    def geometric_curve(rows: list[dict[str, str]]) -> None:
        for row in rows:
            if row["label"] != "amount_cumulative_hctz25_sus" or row["time"] == "0":
                continue
            if row["time"] == "2":
                row.update(mean="", sd="", max="4", comment="below the LOQ")
            else:
                row.update(gmean=row["mean"], mean="", sd="")

    def repeated_auc(rows: list[dict[str, str]]) -> None:
        with (PATEL1984 / "outputs_Tab2.tsv").open(encoding="utf-8") as file:
            auc = next(csv.DictReader(file, delimiter="\t"))
        assert auc["measurement"] == "auc_end"
        rows.append({**auc, "source": "Tab3"})

    def group_without_characteristics(rows: list[dict[str, str]]) -> None:
        rows.append({**rows[0], "name": "men", "parent": "all", "count": "6"})

    def no_sid(rows: list[dict[str, str]]) -> None:
        for row in rows:
            if row["label"] == "hctz25_sus":
                row["label"] = "hctz-25_sus"

    edit_table(study / "timecourses_Fig2.tsv", geometric_curve)
    edit_table(study / "outputs_Tab3.tsv", repeated_auc)
    edit_table(study / "subjects.tsv", group_without_characteristics)
    edit_table(study / "timecourses_Fig1.tsv", no_sid)
    format_folder(study)
    tables, messages = collect_warnings(lambda: study_tables(study))
    return study, tables, messages


def test_the_counts_are_integers(
    followups: tuple[Path, dict[str, pd.DataFrame], list[str]], tmp_path: Path
) -> None:
    """A count is an integer, also in a column with a missing count (not 12.0)."""
    study, tables, _ = followups
    groups = tables["TabGroups"]
    assert groups["count"].isna().any()
    for source, column in (
        ("TabGroups", "count"),
        ("TabGroups", "group_count"),
        ("TabIndividuals", "count"),
        ("Fig1", "count"),
        ("TabInterventions", "doses"),
    ):
        assert str(tables[source][column].dtype) == "Int64", (source, column)
    with pytest.warns(CollectionWarning):
        paths = write_study_tables(study, tmp_path)
    path = next(path for path in paths if path.name.endswith("_TabGroups.tsv"))
    written = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    assert set(written["count"]) == {"12", "NA"}


def test_a_row_below_the_loq_has_the_calculation_of_its_curve(
    followups: tuple[Path, dict[str, pd.DataFrame], list[str]],
) -> None:
    """A row below the LOQ in a curve of geometric means is a geometric mean."""
    _, tables, _ = followups
    df = tables["Fig2"]
    curve = df[df["label"] == "amount_cumulative_hctz25_sus"]
    loq = curve[curve["time"] == 2].iloc[0]
    assert loq[["value", "mean", "median", "gmean"]].isna().all()
    assert loq["max"] == 4
    assert set(curve["calculation_type"]) == {"geometric mean"}
    # the curves of means stay
    other = df[df["label"] == "amount_cumulative_hctz25_tab"]
    assert set(other["calculation_type"]) == {"sample mean"}


def test_the_curator_zero_of_a_curve_of_geometric_means_is_dropped(
    followups: tuple[Path, dict[str, pd.DataFrame], list[str]],
) -> None:
    """A mean of 0 cannot be a geometric mean: it is left out with a warning."""
    _, tables, messages = followups
    df = tables["Fig2"]
    curve = df[df["label"] == "amount_cumulative_hctz25_sus"]
    assert 0 not in set(curve["time"])
    assert len(curve) == 9
    # the 0 of a curve of means stays
    other = df[df["label"] == "amount_cumulative_hctz25_tab"]
    assert 0 in set(other["time"])
    assert any(
        "Patel1984" in message
        and "Fig2" in message
        and "amount_cumulative_hctz25_sus" in message
        and "geometric" in message
        for message in messages
    )


def test_a_duplicate_across_sources_is_a_warning(
    followups: tuple[Path, dict[str, pd.DataFrame], list[str]],
) -> None:
    """The same output in two sources of a study is reported with both sources."""
    _, tables, messages = followups
    duplicates = [message for message in messages if "duplicate" in message]
    assert len(duplicates) == 1
    assert "Patel1984" in duplicates[0]
    assert "Tab2" in duplicates[0]
    assert "Tab3" in duplicates[0]
    assert "auc_end" in duplicates[0]
    # both rows are kept
    for source in ("Tab2", "Tab3"):
        df = tables[source]
        assert (df["measurement_type"] == "auc_end").any()


def test_a_label_which_is_no_sid_is_a_warning(
    followups: tuple[Path, dict[str, pd.DataFrame], list[str]],
) -> None:
    """A curated label which is no SId stays, with a warning."""
    _, tables, messages = followups
    assert "hctz-25_sus" in set(tables["Fig1"]["label"])
    assert any(
        "Patel1984" in message and "'hctz-25_sus'" in message for message in messages
    )
    # the generated labels are SIds
    for df in tables.values():
        if "label" not in df.columns or "measurement_type" not in df.columns:
            continue
        for label in df["label"].dropna():
            if label != "hctz-25_sus":
                assert re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", label), label


def test_a_study_without_the_cases_gives_no_warning() -> None:
    """The unedited study warns of nothing."""
    _, messages = collect_warnings(lambda: study_tables(PATEL1984))
    assert messages == []


def test_a_generated_label_is_an_sid() -> None:
    """A generated label starting with a digit gets an underscore."""
    assert _label("24h urine", "sodium") == "_24h_urine_sodium"
    assert _label("auc", "hctz") == "auc_hctz"


@pytest.fixture(scope="module")
def nilsen1989_tables() -> dict[str, pd.DataFrame]:
    """The tables of the study with scatters."""
    return study_tables(NILSEN1989)


def test_a_scatter_has_the_times_and_interventions_of_both_axes(
    nilsen1989_tables: dict[str, pd.DataFrame],
) -> None:
    """x and y keep their own time and interventions, the joined ones stay."""
    source = pd.read_csv(NILSEN1989 / "scatters_Fig5.tsv", sep="\t")
    df = nilsen1989_tables["Fig5"]
    assert df["x_time"].tolist() == pytest.approx(source["x_time"].tolist())
    assert df["y_time"].tolist() == pytest.approx(source["y_time"].tolist())
    assert df["time"].tolist() == df["y_time"].tolist()
    assert set(df["x_time_unit"]) == set(df["y_time_unit"]) == {"hr"}
    assert df["x_intervention"].tolist() == source["x_interventions"].tolist()
    assert df["y_intervention"].tolist() == source["y_interventions"].tolist()
    assert df["intervention"].tolist() == source["y_interventions"].tolist()
    assert set(df["x_substance"]) == {"hydrochlorothiazide"}
    assert set(df["y_substance"]) == {"sodium"}


def test_a_scatter_of_group_means_has_the_count_of_the_group(
    nilsen1989_tables: dict[str, pd.DataFrame],
) -> None:
    """A point of a group counts the group, an anonymous subject of a group one."""
    # Fig5: the means of the 9 volunteers at the times of the collections
    assert set(nilsen1989_tables["Fig5"]["count"]) == {9}
    # Fig2: 8 points of the volunteers with the same interventions and times
    assert set(nilsen1989_tables["Fig2"]["count"]) == {1}


def test_collect_further_studies(source_dir: Path, tmp_path: Path) -> None:
    """A study folder of another substance is collected with the substance."""
    other = tmp_path / "pkdb_data" / "studies" / "hydrochlorothiazide"
    shutil.copytree(PATEL1984, other / "Patel1984")
    target_dir = tmp_path / "data"
    collect_tsv_files(
        source_dir=source_dir,
        target_dir=target_dir,
        filter_study=lambda name: name == "Hindle1992",
        studies=[other / "Patel1984"],
    )
    folders = {path.split("/")[0] for path in tree(target_dir)}
    assert folders == {"Hindle1992", "Patel1984"}
    assert (target_dir / "Patel1984" / ".Patel1984_Fig1.tsv").is_file()


def test_collect_keeps_files_on_request(source_dir: Path, tmp_path: Path) -> None:
    """With `keep_files` a table the collection did not write stays."""
    shutil.copytree(PATEL1984, source_dir / "Patel1984")
    target_dir = tmp_path / "data"
    (target_dir / "Patel1984").mkdir(parents=True)
    extra = target_dir / "Patel1984" / ".Patel1984_Fig2_pH.tsv"
    extra.write_text("pH", encoding="utf-8")
    outdated = target_dir / "Patel1984" / ".Patel1984_Fig1.tsv"
    outdated.write_text("old", encoding="utf-8")
    collect_tsv_files(source_dir=source_dir, target_dir=target_dir, keep_files=True)
    assert extra.read_text(encoding="utf-8") == "pH"
    # a table of the study is written
    assert outdated.read_text(encoding="utf-8") != "old"
    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)
    assert not extra.exists()


def test_the_geometric_statistics_of_a_group(tmp_path: Path) -> None:
    """A characteristic of a group with a geometric mean has it in the table."""
    from pkdb.studyformat.formatter import format_folder

    study = tmp_path / "Patel1984"
    shutil.copytree(PATEL1984, study)

    def geometric_age(rows: list[dict[str, str]]) -> None:
        for row in rows:
            if row["measurement"] == "age":
                row.update(mean="", gmean="24", gcv="20")

    edit_table(study / "characteristica.tsv", geometric_age)
    format_folder(study)
    groups = study_tables(study)["TabGroups"]
    assert {"cv", "gmean", "gsd", "gcv"} <= set(groups.columns)
    age = groups[groups["measurement_type"] == "age"].iloc[0]
    assert pd.isna(age["mean"])
    assert age["gmean"] == pytest.approx(24)
    assert age["gcv"] == pytest.approx(0.2)


def test_the_mean_of_individuals() -> None:
    """The mean, SD and count of the individuals of a curve at a time."""
    df = pd.DataFrame(
        {
            "study": "S",
            "group": "all",
            "individual": ["i1", "i2", "i3", "i1", "i2", "i3"],
            "intervention": "D1",
            "substance": "x",
            "tissue": "plasma",
            "measurement_type": "concentration",
            "time": [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
            "time_unit": "hr",
            "value": [1.0, 2.0, 3.0, 4.0, None, 8.0],
            "unit": "mg/l",
        }
    )
    means = mean_of_individuals(df)
    assert means["time"].tolist() == [1.0, 2.0]
    assert means["mean"].tolist() == pytest.approx([2.0, 6.0])
    assert means["sd"].tolist() == pytest.approx([1.0, math.sqrt(8)])
    assert means["count"].tolist() == [3, 2]
    assert set(means["unit"]) == {"mg/l"}
    assert set(means["time_unit"]) == {"hr"}


def test_the_sampling_times() -> None:
    """The points of a curve move in their order to the nearest sampling times."""
    df = pd.DataFrame(
        {
            "label": ["a", "a", "a", "a", "b", "b"],
            "time": [-0.1, 0.48, 0.55, 2.2, 0.9, 0.8],
            "mean": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        }
    )
    moved = sampling_times(df, [0.0, 0.5, 1.0, 2.0])
    a = moved[moved["label"] == "a"]
    # two points near 0.5 go to distinct samples, in their order
    assert a["time"].tolist() == [0.0, 0.5, 1.0, 2.0]
    assert a["mean"].tolist() == [1.0, 2.0, 3.0, 4.0]
    b = moved[moved["label"] == "b"]
    assert b["time"].tolist() == [0.5, 1.0]
    assert b["mean"].tolist() == [6.0, 5.0]
    # the data are not changed
    assert df["time"].tolist()[0] == -0.1
    # the points before `start` keep their times
    kept = sampling_times(df, [0.0, 0.5, 1.0, 2.0], start=0.5)
    assert kept[kept["label"] == "a"]["time"].tolist() == [-0.1, 0.48, 0.5, 2.0]
    # a curve with more points than samples
    with pytest.raises(ValueError, match="'a' has 4 points for 2 sampling times"):
        sampling_times(df, [0.0, 1.0])


def test_a_scatter_of_partly_repeated_points_is_a_warning(tmp_path: Path) -> None:
    """Points of a group of which some have the same times are a warning."""
    from pkdb.studyformat.formatter import format_folder

    study = tmp_path / "Nilsen1989"
    shutil.copytree(NILSEN1989, study)

    def same_times(rows: list[dict[str, str]]) -> None:
        rows[1].update(x_time=rows[0]["x_time"], y_time=rows[0]["y_time"])

    edit_table(study / "scatters_Fig5.tsv", same_times)
    format_folder(study)
    tables, messages = collect_warnings(lambda: study_tables(study))
    df = tables["Fig5"]
    curve = df[df["label"] == "excretion_sodium2_hctz25"]
    assert sorted(curve["count"]) == [1, 1, 9, 9, 9]
    assert any(
        "Nilsen1989" in message
        and "excretion_sodium2_hctz25" in message
        and "volunteers" in message
        for message in messages
    )
