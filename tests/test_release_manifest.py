"""Tests of the release manifest and the pruning of the release repository."""

from pathlib import Path

import pytest

from pkpdsim.release.manifest import Manifest, matches, prune, write_manifest
from pkpdsim.release.source import SourceRepo, run

GIT_ENV = {
    "GIT_AUTHOR_NAME": "Test",
    "GIT_AUTHOR_EMAIL": "test@example.org",
    "GIT_COMMITTER_NAME": "Test",
    "GIT_COMMITTER_EMAIL": "test@example.org",
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_SYSTEM": "/dev/null",
}
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


@pytest.fixture(autouse=True)
def git_env(monkeypatch: pytest.MonkeyPatch) -> None:
    for key, value in GIT_ENV.items():
        monkeypatch.setenv(key, value)


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
