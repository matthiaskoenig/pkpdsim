"""Tests of the release manifest and the pruning of the release repository."""

import json
from pathlib import Path

import pytest

from pkpdsim.release.manifest import (
    Manifest,
    changed_by_build,
    check_file_sizes,
    check_not_ignored,
    finish_build,
    matches,
    prune,
    start_build,
    stop_build,
    write_manifest,
)
from pkpdsim.release.source import ReleaseError, SourceRepo, run

LOCK = """version = 1

[[package]]
name = "pkpdsim"
version = "0.8.0"

[[package]]
name = "sbmlsim"
version = "0.3.1"

[[package]]
name = "sbmlutils"
version = "0.9.2"
"""
SOURCE = SourceRepo(root=Path("/unused"), commit="abc123")


pytestmark = pytest.mark.usefixtures("git_env")


def write(root: Path, relative: str, text: str = "x\n") -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def release(tmp_path: Path) -> tuple[Path, list[Path]]:
    """A release repository with an older commit and an uncommitted build."""
    root = tmp_path / "release"
    root.mkdir()
    run(["git", "init", "-b", "main"], root)
    write(root, ".gitignore", "results/\n.venv/\n")
    write(root, "release-notes/0.7.1.md", "old\n")
    write(root, "docs/RELEASE.md", "old\n")
    write(root, "src/pkg/run_hctz.py", "old\n")
    run(["git", "add", "-A"], root)
    run(["git", "commit", "-m", "old"], root)
    write(root, "results/x", "mine\n")
    write(root, ".venv/lib/y", "env\n")
    write(root, "simulations/a/index.html", "a\n")
    written = [
        write(root, "README.md", "readme\n"),
        write(root, "uv.lock", LOCK),
        write(root, ".gitignore", "results/\n.venv/\n"),
        write(root, "release-notes/0.8.0.md", "new\n"),
        root / "simulations",  # a stage may return a folder
    ]
    return root, written


def test_manifest_roundtrip(release: tuple[Path, list[Path]]) -> None:
    root, written = release
    manifest = write_manifest(root, "0.8.0", SOURCE, written)
    assert manifest.pkpdsim == "0.8.0"
    assert manifest.sbmlsim == "0.3.1"
    assert manifest.sbmlutils == "0.9.2"
    assert "simulations/a/index.html" in manifest.files
    assert "release.json" not in manifest.files
    assert Manifest.read(root / "release.json") == manifest
    assert Manifest.read(root / "missing.json") is None
    assert b"\r" not in (root / "release.json").read_bytes()


def test_manifest_local_pkpdsim(release: tuple[Path, list[Path]]) -> None:
    """A local pkpdsim (`--pkpdsim-path`) has no version in the lock."""
    root, written = release
    lock = LOCK.replace(
        'version = "0.8.0"', 'source = { editable = "/home/user/pkpdsim" }'
    )
    write(root, "uv.lock", lock)
    manifest = write_manifest(root, "0.8.0", SOURCE, written)
    assert manifest.pkpdsim == "editable /home/user/pkpdsim"
    assert not matches(root, "0.8.0", "abc123", "0.8.0")


def test_matches_after_build(release: tuple[Path, list[Path]]) -> None:
    root, written = release
    write_manifest(root, "0.8.0", SOURCE, [p.relative_to(root) for p in written])
    assert matches(root, "0.8.0", "abc123", "0.8.0")


def test_hand_edit_breaks_match(release: tuple[Path, list[Path]]) -> None:
    root, written = release
    write_manifest(root, "0.8.0", SOURCE, written)
    with (root / "README.md").open("a") as f:
        f.write("x")
    assert not matches(root, "0.8.0", "abc123", "0.8.0")


def test_extra_file_breaks_match(release: tuple[Path, list[Path]]) -> None:
    root, written = release
    write_manifest(root, "0.8.0", SOURCE, written)
    write(root, "extra.txt")
    assert not matches(root, "0.8.0", "abc123", "0.8.0")


def test_other_commit_breaks_match(release: tuple[Path, list[Path]]) -> None:
    root, written = release
    write_manifest(root, "0.8.0", SOURCE, written)
    assert not matches(root, "0.8.0", "other", "0.8.0")
    assert not matches(root, "0.9.0", "abc123", "0.8.0")
    assert not matches(root, "0.8.0", "abc123", "0.8.1")
    assert not matches(root, "0.8.0", "abc123", None)


def test_prune_keeps_release_notes_and_ignored(
    release: tuple[Path, list[Path]],
) -> None:
    root, written = release
    removed = prune(root, written)
    assert (root / "release-notes" / "0.7.1.md").exists()
    assert (root / "results" / "x").exists()
    assert (root / ".venv" / "lib" / "y").exists()
    assert (root / ".git").is_dir()
    assert (root / "README.md").exists()
    assert not (root / "docs").exists()
    assert not (root / "src").exists()
    assert root / "docs" / "RELEASE.md" in removed
    assert root / "src" / "pkg" / "run_hctz.py" in removed


def test_empty_folder_of_ignored_directory_survives(
    release: tuple[Path, list[Path]],
) -> None:
    root, written = release
    (root / ".venv" / "include").mkdir()
    (root / "results" / "empty").mkdir()
    prune(root, written)
    assert (root / ".venv" / "include").is_dir()
    assert (root / "results" / "empty").is_dir()


def test_matches_after_prune_of_tracked_files(
    release: tuple[Path, list[Path]],
) -> None:
    root, written = release
    prune(root, written)
    write_manifest(root, "0.8.0", SOURCE, written)
    assert matches(root, "0.8.0", "abc123", "0.8.0")
    with (root / "README.md").open("a") as f:
        f.write("x")
    assert not matches(root, "0.8.0", "abc123", "0.8.0")


def test_missing_written_path_is_a_release_error(
    release: tuple[Path, list[Path]],
) -> None:
    root, written = release
    with pytest.raises(ReleaseError, match=r"nope\.txt"):
        write_manifest(root, "0.8.0", SOURCE, [*written, Path("nope.txt")])


def test_changed_by_build(release: tuple[Path, list[Path]]) -> None:
    root, written = release
    assert not changed_by_build(root, "0.8.0")
    write_manifest(root, "0.8.0", SOURCE, written)
    assert changed_by_build(root, "0.8.0")
    with (root / "README.md").open("a") as f:
        f.write("edited by hand")
    assert changed_by_build(root, "0.8.0")
    assert not changed_by_build(root, "0.9.0")
    write(root, "extra.txt")
    assert not changed_by_build(root, "0.8.0")


def test_changed_by_build_clean(release: tuple[Path, list[Path]]) -> None:
    root, _ = release
    run(["git", "add", "-A"], root)
    run(["git", "commit", "-m", "all"], root)
    assert changed_by_build(root, "0.8.0")
    (root / "docs" / "RELEASE.md").unlink()
    assert changed_by_build(root, "0.8.0")


def test_unfinished_build_accepts_its_changes(
    release: tuple[Path, list[Path]],
) -> None:
    root, written = release
    write_manifest(root, "0.8.0", SOURCE, written)
    write(root, "extra.txt")
    assert not changed_by_build(root, "0.8.0")
    marker = root / ".git" / "pkpd_release_build"
    start_build(root, "0.9.0")
    assert json.loads(marker.read_text()) == {"version": "0.9.0"}
    # a killed build recorded no paths, every change may be its change
    assert changed_by_build(root, "0.8.0")
    assert changed_by_build(root, "1.0.0")
    stop_build(root)
    assert "extra.txt" in json.loads(marker.read_text())["changed"]
    assert changed_by_build(root, "1.0.0")
    write(root, "added_by_hand.txt")
    assert not changed_by_build(root, "1.0.0")
    (root / "added_by_hand.txt").unlink()
    finish_build(root)
    assert not marker.exists()
    assert not changed_by_build(root, "0.8.0")


def test_stop_build_without_marker(release: tuple[Path, list[Path]]) -> None:
    root, _ = release
    stop_build(root)
    assert not (root / ".git" / "pkpd_release_build").exists()


def test_marker_of_a_worktree(release: tuple[Path, list[Path]], tmp_path: Path) -> None:
    """In a worktree `.git` is a file, the marker lies in its git folder."""
    root, _ = release
    run(["git", "add", "-A"], root)
    run(["git", "commit", "-m", "all"], root)
    worktree = tmp_path / "worktree"
    run(["git", "worktree", "add", "-b", "other", str(worktree)], root)
    assert (worktree / ".git").is_file()
    start_build(worktree, "0.9.0")
    git_dir = Path(run(["git", "rev-parse", "--absolute-git-dir"], worktree).strip())
    assert (git_dir / "pkpd_release_build").is_file()
    write(worktree, "extra.txt")
    removed = prune(worktree, [])
    assert (worktree / ".git").is_file()
    assert worktree / "extra.txt" in removed
    finish_build(worktree)
    assert not (git_dir / "pkpd_release_build").exists()


def test_check_not_ignored(release: tuple[Path, list[Path]]) -> None:
    root, written = release
    check_not_ignored(root, written)
    tracked = write(root, "src/pkg/fitting/results/data.xlsx")
    write(root, ".gitignore", "/results/\n.venv/\n")
    check_not_ignored(root, [*written, tracked])
    write(root, ".gitignore", "results/\n.venv/\n*.log\n")
    log = write(root, "simulations/a/run.log")
    with pytest.raises(ReleaseError) as err:
        check_not_ignored(root, [*written, tracked])
    assert "src/pkg/fitting/results/data.xlsx" in str(err.value)
    assert "simulations/a/run.log" in str(err.value)
    assert "README.md" not in str(err.value)
    assert log.is_file()


def test_check_file_sizes(release: tuple[Path, list[Path]]) -> None:
    root, written = release
    check_file_sizes(root, written, max_size=1000)
    large = write(root, "simulations/a/a_task.nc", "x" * 2000)
    with pytest.raises(ReleaseError) as err:
        check_file_sizes(root, [*written, large], max_size=1000)
    assert "simulations/a/a_task.nc (0 MB)" in str(err.value)
    assert "README.md" not in str(err.value)
