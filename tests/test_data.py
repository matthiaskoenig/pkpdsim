from pathlib import Path

import pytest

from pkpdsim.data import PKDB_DATA, collect_tsv_files, studies_dir, study_format


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


def test_a_study_of_format_2_is_not_collected(source_dir: Path, tmp_path: Path) -> None:
    """A study of format 2 has no hidden tsv files: nothing is copied (#109)."""
    study = source_dir / "Mendell2012"
    study.mkdir()
    (study / "study.json").write_text('{"format": 2}', encoding="utf-8")
    (study / "timecourses_Fig1.tsv").write_text("study", encoding="utf-8")
    target_dir = tmp_path / "data"
    with pytest.raises(ValueError, match=r"Mendell2012 \(format 2\)"):
        collect_tsv_files(source_dir=source_dir, target_dir=target_dir)
    assert not target_dir.exists()
    # a filter which leaves the study out collects the others
    collect_tsv_files(
        source_dir=source_dir,
        target_dir=target_dir,
        filter_study=lambda name: name != "Mendell2012",
    )
    assert "Hindle1992/.Hindle1992_Fig2.tsv" in tree(target_dir)
