"""The backfill of the study records from the git history (`scripts/backfill_study_records.py`)."""

import importlib.util
import json
import os
import subprocess
import sys
from datetime import date
from pathlib import Path
from types import ModuleType

import pytest

from pkpdsim.data import StudyRecord, read_study_record, study_record_path

ROOT = Path(__file__).parents[1]


@pytest.fixture(scope="module")
def script() -> ModuleType:
    """Import `scripts/backfill_study_records.py` from its path."""
    path = ROOT / "scripts" / "backfill_study_records.py"
    spec = importlib.util.spec_from_file_location("backfill_study_records", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # the dataclasses of the script look their module up
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _commit(repo: Path, day: str, message: str) -> str:
    env = {
        **os.environ,
        "GIT_AUTHOR_NAME": "Test",
        "GIT_AUTHOR_EMAIL": "test@example.org",
        "GIT_COMMITTER_NAME": "Test",
        "GIT_COMMITTER_EMAIL": "test@example.org",
        "GIT_AUTHOR_DATE": f"{day}T12:00:00+00:00",
        # a rebase moves the commit date, never the author date
        "GIT_COMMITTER_DATE": "2026-11-01T12:00:00+00:00",
    }
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True, env=env)
    subprocess.run(
        ["git", "-C", str(repo), "commit", "-q", "-m", message], check=True, env=env
    )
    return subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _table(folder: Path, source: str, text: str) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f".{folder.name}_{source}.tsv").write_text(text, encoding="utf-8")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A repository of models whose model `a` moved to `b`."""
    repo = tmp_path / "models"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    _table(repo / "a" / "data" / "S1", "Fig1", "first\n")
    _table(repo / "a" / "data" / "S3", "Fig1", "first\n")
    (repo / "a" / "data" / "notes.tsv").write_text("", encoding="utf-8")
    _commit(repo, "2020-01-01", "model a")
    subprocess.run(["git", "-C", str(repo), "mv", "a", "b"], check=True)
    _table(repo / "b" / "data" / "S2", "Fig1", "second\n")
    _commit(repo, "2021-01-01", "model a is b, collect S2")
    # a later collection replaces the sources of S1
    (repo / "b" / "data" / "S1" / ".S1_Fig1.tsv").unlink()
    _table(repo / "b" / "data" / "S1", "Tab1", "third\n")
    _commit(repo, "2022-01-01", "collect S1 again")
    # an experiment of S5 before its tables, read from pkdb_data
    (repo / "c" / "experiments").mkdir(parents=True)
    (repo / "c" / "experiments" / "s5.py").write_text("class S5:\n    pass\n")
    _commit(repo, "2019-01-01", "experiment of S5")
    _table(repo / "c" / "data" / "S5", "Fig1", "fifth\n")
    _commit(repo, "2023-01-01", "c keeps its data")
    # a model imported with its data and a fit of before the import
    _table(repo / "d" / "data" / "S6", "Fig1", "sixth\n")
    _fit(repo / "d", "PK_20231201_120000__abcde")
    _commit(repo, "2024-01-01", "import d")
    # a study collected after a fit was adopted
    _fit(repo / "e", "PK_20240115_120000__abcde")
    _commit(repo, "2024-02-01", "adopt the fit of e")
    _table(repo / "e" / "data" / "S7", "Fig1", "seventh\n")
    _commit(repo, "2024-03-01", "collect S7")
    return repo


def _fit(model: Path, fit_id: str) -> None:
    folder = model / "models" / "fits" / fit_id
    folder.mkdir(parents=True)
    (folder / "parameters.json").write_text("{}", encoding="utf-8")


def test_the_days_of_the_studies(script: ModuleType, repo: Path) -> None:
    result = script.backfill(repo, dry_run=True)
    days = {folder: record.collected for folder, record in result.written.items()}
    assert days == {
        "b/data/S1": date(2020, 1, 1),
        "b/data/S2": date(2021, 1, 1),
        "b/data/S3": date(2020, 1, 1),
        # the experiment mentioned the study before its tables
        "c/data/S5": date(2019, 1, 1),
        # the fit was committed with the study, so it knew it
        "d/data/S6": date(2023, 12, 1),
        # the fit was in the repository before any trace of the study
        "e/data/S7": date(2024, 3, 1),
    }
    adopted = script.fits(repo)
    assert {path: day for path, (day, _) in adopted.items()} == {
        "d/models/fits/PK_20231201_120000__abcde/parameters.json": date(2023, 12, 1),
        "e/models/fits/PK_20240115_120000__abcde/parameters.json": date(2024, 1, 15),
    }
    owners = script.fits_of({"b", "d", "e", "e/nested"}, adopted)
    assert [day for day, _ in owners["d"]] == [date(2023, 12, 1)]
    assert owners["b"] == owners["e/nested"] == []


def test_the_datings(script: ModuleType, repo: Path) -> None:
    found = script.datings(repo)
    assert found["b/data/S1"].first == date(2020, 1, 1)
    assert found["b/data/S1"].latest == date(2022, 1, 1)
    assert found["b/data/S2"].first == date(2021, 1, 1)
    assert "a/data" not in found


def test_the_backfill(script: ModuleType, repo: Path, tmp_path: Path) -> None:
    pkdb_data = tmp_path / "pkdb_data" / "studies"
    (pkdb_data / "substance" / "S1").mkdir(parents=True)
    (pkdb_data / "substance" / "S1" / "study.json").write_text(
        json.dumps({"release": {"date": "2024-03-26", "pkdb_id": "PKDB00001"}})
    )
    # a study with a record keeps it, one with a change of a table is left out
    kept = StudyRecord(study="S2", collected=date(2019, 5, 5), dating="collection")
    study_record_path(repo / "b" / "data" / "S2").write_text(json.dumps(kept.to_dict()))
    _table(repo / "b" / "data" / "S3", "Fig1", "edited\n")
    _table(repo / "b" / "data" / "S4", "Fig1", "new\n")
    assert script.main([str(repo), "--pkdb-data", str(pkdb_data), "--dry-run"]) == 0
    assert read_study_record(repo / "b" / "data" / "S1") is None
    result = script.backfill(repo, pkdb_data=pkdb_data)
    assert list(result.written) == ["b/data/S1", "c/data/S5", "d/data/S6", "e/data/S7"]
    assert result.recorded == ["b/data/S2"]
    assert result.changed == ["b/data/S3"]
    assert result.untracked == ["b/data/S4"]
    record = read_study_record(repo / "b" / "data" / "S1")
    assert record is not None
    assert (record.collected, record.dating, record.updated) == (
        date(2020, 1, 1),
        "git",
        date(2022, 1, 1),
    )
    assert record.dating_commit is not None
    assert len(record.dating_commit) == 40
    assert (record.released, record.pkdb_id) == (date(2024, 3, 26), "PKDB00001")
    assert record.changed_tables(repo / "b" / "data" / "S1") == []
    assert read_study_record(repo / "b" / "data" / "S2") == kept
    # forced, the record of S2 is replaced from the history
    script.backfill(repo, force=True)
    replaced = read_study_record(repo / "b" / "data" / "S2")
    assert replaced is not None
    assert (replaced.collected, replaced.dating) == (date(2021, 1, 1), "git")
