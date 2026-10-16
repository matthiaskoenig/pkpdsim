"""Tests of the source repository, the release repository and the export."""

import dataclasses
from pathlib import Path

import pytest
from release_git import git

from pkpdsim.model import PKPDModel
from pkpdsim.release.settings import ReleaseSettings
from pkpdsim.release.source import (
    ReleaseError,
    check_version,
    export_package,
    release_repo,
    run,
    source_repo,
)
from pkpdsim_toy import MODEL

assert MODEL.release is not None
SETTINGS: ReleaseSettings = MODEL.release
pytestmark = pytest.mark.usefixtures("git_env")


def test_source_dirty_refused(source_repo_model: PKPDModel) -> None:
    (source_repo_model.path / "new.txt").write_text("untracked\n")
    with pytest.raises(ReleaseError, match="commit the changes of"):
        source_repo(source_repo_model)


def test_source_commit(source_repo_model: PKPDModel) -> None:
    source = source_repo(source_repo_model)
    assert source.commit == git(source.root, "rev-parse", "HEAD").strip()
    assert source.root == source_repo_model.path.parent.parent.resolve()


def test_release_repo_cloned_when_missing(
    tmp_path: Path, source_repo_model: PKPDModel, release_settings: ReleaseSettings
) -> None:
    source = source_repo(source_repo_model)
    path = tmp_path / "elsewhere" / "clone"
    path.parent.mkdir()
    repo = release_repo(release_settings, path, source)
    assert (repo.path / "README.md").exists()
    assert repo.last_tag == "0.1.0"


def test_release_repo_default_path(
    source_repo_model: PKPDModel, release_settings: ReleaseSettings
) -> None:
    source = source_repo(source_repo_model)
    repo = release_repo(release_settings, None, source)
    assert repo.path == source.root.parent / "toy-model"


def test_release_repo_other_origin_refused(
    source_repo_model: PKPDModel,
    remote_and_clone: tuple[Path, Path],
    release_settings: ReleaseSettings,
) -> None:
    _, clone = remote_and_clone
    git(clone, "remote", "set-url", "origin", "https://github.com/other/thing.git")
    source = source_repo(source_repo_model)
    with pytest.raises(ReleaseError) as err:
        release_repo(release_settings, clone, source)
    assert "other/thing" in str(err.value)
    assert "pkpdsim/toy-model" in str(err.value)


def test_release_repo_wrong_branch_refused(
    source_repo_model: PKPDModel,
    remote_and_clone: tuple[Path, Path],
    release_settings: ReleaseSettings,
) -> None:
    _, clone = remote_and_clone
    git(clone, "checkout", "-b", "feature")
    with pytest.raises(ReleaseError, match="not 'main'"):
        release_repo(release_settings, clone, source_repo(source_repo_model))


def test_release_repo_behind_refused(
    tmp_path: Path,
    source_repo_model: PKPDModel,
    remote_and_clone: tuple[Path, Path],
    release_settings: ReleaseSettings,
) -> None:
    remote, clone = remote_and_clone
    other = tmp_path / "other"
    git(tmp_path, "clone", str(remote), str(other))
    (other / "new.txt").write_text("new\n")
    git(other, "add", "-A")
    git(other, "commit", "-m", "more")
    git(other, "push", "origin", "main")
    with pytest.raises(ReleaseError, match="behind"):
        release_repo(release_settings, clone, source_repo(source_repo_model))


def test_check_version() -> None:
    check_version("0.2.0", "0.1.0")
    check_version("0.10.0", "0.9.0")
    check_version("0.1.0", None)
    with pytest.raises(ReleaseError):
        check_version("0.1.0", "0.1.0")
    with pytest.raises(ReleaseError):
        check_version("0.0.9", "0.1.0")
    with pytest.raises(ReleaseError, match="not a version"):
        check_version("1.0", None)


def test_run_failure_message(tmp_path: Path) -> None:
    with pytest.raises(ReleaseError, match="git") as err:
        run(["git", "rev-parse", "--verify", "nothing"], tmp_path)
    assert "git rev-parse --verify nothing" in str(err.value)


def export(tmp_path: Path, model: PKPDModel) -> tuple[Path, list[Path]]:
    release = tmp_path / "release"
    release.mkdir()
    written = export_package(model, SETTINGS, source_repo(model), release)
    return release, written


def test_export_excludes(tmp_path: Path, source_repo_model: PKPDModel) -> None:
    release, written = export(tmp_path, source_repo_model)
    target = release / "src" / "pkpdsim_toy"
    assert (target / "__init__.py").exists()
    assert not (target / "TODO.md").exists()
    assert not (target / "literature").exists()
    assert Path("src/pkpdsim_toy/__init__.py") in written
    assert all(not p.is_absolute() for p in written)
    assert not (release / "src" / "__init__.py").exists()


def test_export_namespace(tmp_path: Path, source_repo_model: PKPDModel) -> None:
    model = dataclasses.replace(source_repo_model, package="pkg.pkpdsim_toy")
    release = tmp_path / "release"
    (release / "src" / "pkg").mkdir(parents=True)
    (release / "src" / "pkg" / "__init__.py").write_text("")
    (release / "src" / "__init__.py").write_text("")
    export_package(model, SETTINGS, source_repo(model), release)
    assert not (release / "src" / "pkg" / "__init__.py").exists()
    assert not (release / "src" / "__init__.py").exists()
    assert (release / "src" / "pkg" / "pkpdsim_toy" / "__init__.py").exists()


def test_export_replaces(tmp_path: Path, source_repo_model: PKPDModel) -> None:
    release = tmp_path / "release"
    stale = release / "src" / "pkpdsim_toy" / "stale.py"
    stale.parent.mkdir(parents=True)
    stale.write_text("old\n")
    export_package(source_repo_model, SETTINGS, source_repo(source_repo_model), release)
    assert not stale.exists()
