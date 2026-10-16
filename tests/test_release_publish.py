"""Tests of the publishing of a release."""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import cast

import pytest

from pkpdsim.release.publish import check_publishable, publish
from pkpdsim.release.settings import ReleaseSettings
from pkpdsim.release.source import ReleaseError
from pkpdsim_toy import MODEL

SETTINGS = cast(ReleaseSettings, MODEL.release)
VERSION = "1.0.0"
IMAGE = SETTINGS.docker_image
REPO = SETTINGS.repository


class FakeRunner:
    """Keeps the state of git, GitHub and the registry between calls."""

    def __init__(self) -> None:
        self.commands: list[list[str]] = []
        self.dirty = True
        self.head = "c0"
        self.head_version: str | None = None
        self.tags: dict[str, str] = {}
        self.remote: dict[str, str] = {}
        self.pages = False
        self.images: set[str] = set()
        self.releases: set[str] = set()
        self.fail_on: str | None = None
        self.missing: set[str] = set()

    def __call__(self, args: Sequence[str], cwd: Path) -> str:
        cmd = [str(a) for a in args]
        self.commands.append(cmd)
        joined = " ".join(cmd)
        if self.fail_on and self.fail_on in joined:
            self.fail_on = None
            raise ReleaseError(f"'{joined}' failed")
        if joined in self.missing:
            raise ReleaseError(f"'{joined}' failed")
        return self._run(cmd, joined)

    def _run(self, cmd: list[str], joined: str) -> str:
        if cmd[:3] == ["git", "status", "--porcelain"]:
            return " M x\n" if self.dirty else ""
        if cmd[:2] == ["git", "show"]:
            if self.head_version is None:
                raise ReleaseError("no manifest")
            return json.dumps({"version": self.head_version})
        if cmd[:2] == ["git", "add"]:
            return ""
        if cmd[:2] == ["git", "commit"]:
            self.dirty = False
            self.head = "c1"
            self.head_version = VERSION
            return ""
        if cmd[:3] == ["git", "tag", "--list"]:
            return cmd[3] if cmd[3] in self.tags else ""
        if cmd[:2] == ["git", "rev-list"]:
            return self.tags[cmd[-1]]
        if cmd[:2] == ["git", "rev-parse"]:
            return self.head
        if cmd[:2] == ["git", "tag"]:
            self.tags[cmd[2]] = self.head
            return ""
        if cmd[:2] == ["git", "ls-remote"]:
            return "".join(f"{sha}\t{ref}\n" for ref, sha in self.remote.items())
        if cmd[:2] == ["git", "push"]:
            self.remote["refs/heads/main"] = self.head
            self.remote[f"refs/tags/{cmd[-1]}"] = self.tags[cmd[-1]]
            return ""
        if cmd[:2] == ["gh", "api"]:
            return self._pages(cmd)
        if cmd[:3] == ["docker", "manifest", "inspect"]:
            if cmd[3] not in self.images:
                raise ReleaseError("no such manifest")
            return ""
        if cmd[:2] == ["docker", "build"]:
            return ""
        if cmd[:2] == ["docker", "push"]:
            self.images.add(cmd[2])
            return ""
        if cmd[:3] == ["gh", "release", "view"]:
            if cmd[3] not in self.releases:
                raise ReleaseError("release not found")
            return ""
        if cmd[:3] == ["gh", "release", "create"]:
            self.releases.add(cmd[3])
            return ""
        return ""

    def _pages(self, cmd: list[str]) -> str:
        if cmd[2] == "-X":
            self.pages = True
            return ""
        if not self.pages:
            raise ReleaseError("gh: Not Found (HTTP 404)")
        return "{}"


ARCHIVE = Path("toy.zip")


def _publish(runner: FakeRunner) -> list[str]:
    return publish(SETTINGS, Path("/release"), VERSION, ARCHIVE, runner)


def test_publish_all_steps() -> None:
    runner = FakeRunner()
    steps = _publish(runner)
    assert steps == ["commit", "tag", "push", "pages", "docker", "github_release"]
    mutating = [
        c
        for c in runner.commands
        if c[:2] in (["git", "add"], ["git", "commit"], ["git", "push"])
        or c[:2] in (["docker", "build"], ["docker", "push"])
        or c[:3] == ["gh", "release", "create"]
        or c == ["git", "tag", VERSION]
        or c[:4] == ["gh", "api", "-X", "POST"]
    ]
    assert mutating == [
        ["git", "add", "-A"],
        ["git", "commit", "-m", f"Release {VERSION}"],
        ["git", "tag", VERSION],
        ["git", "push", "origin", "main", VERSION],
        ["gh", "api", "-X", "POST", f"repos/{REPO}/pages", "-f", "build_type=workflow"],
        ["docker", "build", "-t", f"{IMAGE}:{VERSION}", "-t", f"{IMAGE}:latest", "."],
        ["docker", "push", f"{IMAGE}:{VERSION}"],
        ["docker", "push", f"{IMAGE}:latest"],
        [
            "gh", "release", "create", VERSION,
            "--repo", REPO,
            "--title", f"{SETTINGS.repo_name} {VERSION}",
            "--notes-file", f"release-notes/{VERSION}.md",
            str(ARCHIVE),
        ],
    ]  # fmt: skip


def test_publish_resume_after_docker_failure() -> None:
    runner = FakeRunner()
    runner.fail_on = f"docker push {IMAGE}:{VERSION}"
    with pytest.raises(ReleaseError):
        _publish(runner)
    assert _publish(runner) == ["docker", "github_release"]
    assert runner.releases == {VERSION}


def test_publish_idempotent() -> None:
    runner = FakeRunner()
    _publish(runner)
    assert _publish(runner) == []
    assert _publish(runner) == []


def test_tag_elsewhere_refused() -> None:
    runner = FakeRunner()
    runner.tags[VERSION] = "old"
    with pytest.raises(ReleaseError, match="not at HEAD"):
        _publish(runner)


def test_pages_other_failure_raised() -> None:
    runner = FakeRunner()
    runner.fail_on = "gh api"
    with pytest.raises(ReleaseError, match="gh api"):
        _publish(runner)


def _lock(release: Path, source: str) -> None:
    (release / "uv.lock").write_text(
        'version = 1\n\n[[package]]\nname = "toy-release"\nversion = "1.0.0"\n'
        'source = { editable = "." }\n\n'
        f'[[package]]\nname = "pkpdsim"\nversion = "0.7.1"\nsource = {source}\n'
    )


def test_check_publishable_ok(tmp_path: Path) -> None:
    _lock(tmp_path, '{ registry = "https://pypi.org/simple" }')
    check_publishable(tmp_path, FakeRunner())


def test_check_publishable_local_source(tmp_path: Path) -> None:
    _lock(tmp_path, '{ directory = "../pkpdsim" }')
    with pytest.raises(ReleaseError, match="pkpdsim"):
        check_publishable(tmp_path, FakeRunner())


def test_check_publishable_uv_sources(tmp_path: Path) -> None:
    _lock(tmp_path, '{ registry = "https://pypi.org/simple" }')
    (tmp_path / "pyproject.toml").write_text(
        '[tool.uv.sources]\npkpdsim = { path = "../pkpdsim" }\n'
    )
    with pytest.raises(ReleaseError, match=r"tool\.uv\.sources"):
        check_publishable(tmp_path, FakeRunner())


def test_check_publishable_lists_missing_tools(tmp_path: Path) -> None:
    _lock(tmp_path, '{ registry = "https://pypi.org/simple" }')
    runner = FakeRunner()
    runner.missing = {"docker info", "gh auth status"}
    with pytest.raises(ReleaseError) as err:
        check_publishable(tmp_path, runner)
    assert "docker" in str(err.value)
    assert "gh auth status" in str(err.value)
