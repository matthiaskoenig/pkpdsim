import shutil
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
    _label,
    _statistics,
    collect_tsv_files,
    database_differences,
    prepare_study,
    studies_dir,
    study_format,
    study_tables,
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
    # the value also in the columns of sbmlsim, a point counts one
    assert df["x"].tolist() == df["x_value"].tolist()
    assert set(df["count"]) == {1}


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
