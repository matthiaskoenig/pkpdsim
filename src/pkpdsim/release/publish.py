"""Publish a built release: repository, tag, Pages, Docker image, GitHub release.

Every step checks whether it is already done and is skipped then, so a publish
which failed half way is finished by running the command again. The tools are
the real `git`, `gh` and `docker` command line programs, called through a
`Runner`.
"""

from __future__ import annotations

import json
import tomllib
from pathlib import Path

from pkpdsim.release.manifest import MANIFEST_FILE
from pkpdsim.release.settings import ReleaseSettings
from pkpdsim.release.source import ReleaseError, Runner, run

_PYPI = "https://pypi.org/simple"


def _local_packages(release_path: Path) -> list[str]:
    lock = release_path / "uv.lock"
    if not lock.is_file():
        return []
    data = tomllib.loads(lock.read_text(encoding="utf-8"))
    names = []
    for package in data.get("package", []):
        source = package.get("source", {})
        if source.get("registry") == _PYPI or source in (
            {"editable": "."},
            {"virtual": "."},
        ):
            continue
        names.append(package.get("name", "?"))
    return names


def _tool_missing(runner: Runner, args: list[str], cwd: Path) -> bool:
    try:
        runner(args, cwd)
    except ReleaseError:
        return True
    return False


def check_publishable(release_path: Path, runner: Runner = run) -> None:
    """Check that the release repository can be published.

    Args:
        release_path: the checkout of the release repository.
        runner: runs the tool checks.

    Raises:
        ReleaseError: if the lock has a package which is not from PyPI, the
            project uses `tool.uv.sources`, Docker is not running or `gh` is not
            logged in; one message lists everything.
    """
    problems = []
    local = _local_packages(release_path)
    if local:
        problems.append(
            f"uv.lock has packages which are not from PyPI: {', '.join(local)}"
        )
    pyproject = release_path / "pyproject.toml"
    if pyproject.is_file():
        data = tomllib.loads(pyproject.read_text(encoding="utf-8"))
        if data.get("tool", {}).get("uv", {}).get("sources"):
            problems.append("pyproject.toml has tool.uv.sources")
    if _tool_missing(runner, ["docker", "info"], release_path):
        problems.append("docker is not available (docker info failed)")
    if _tool_missing(runner, ["gh", "auth", "status"], release_path):
        problems.append("gh is not logged in (gh auth status failed)")
    if problems:
        raise ReleaseError("cannot publish:\n- " + "\n- ".join(problems))


def _succeeds(runner: Runner, args: list[str], cwd: Path) -> bool:
    return not _tool_missing(runner, args, cwd)


def _head_has_manifest(path: Path, version: str, runner: Runner) -> bool:
    try:
        text = runner(["git", "show", f"HEAD:{MANIFEST_FILE}"], path)
        return json.loads(text).get("version") == version
    except (ReleaseError, ValueError, AttributeError):
        return False


def _commit(path: Path, version: str, runner: Runner) -> bool:
    clean = not runner(["git", "status", "--porcelain"], path).strip()
    if clean and _head_has_manifest(path, version, runner):
        return False
    runner(["git", "add", "-A"], path)
    runner(["git", "commit", "-m", f"Release {version}"], path)
    return True


def _tag(path: Path, version: str, runner: Runner) -> bool:
    if runner(["git", "tag", "--list", version], path).strip():
        tagged = runner(["git", "rev-list", "-n", "1", version], path).strip()
        if tagged != runner(["git", "rev-parse", "HEAD"], path).strip():
            raise ReleaseError(f"tag {version} exists and is not at HEAD")
        return False
    runner(["git", "tag", version], path)
    return True


def _push(path: Path, version: str, runner: Runner) -> bool:
    head = runner(["git", "rev-parse", "HEAD"], path).strip()
    listing = runner(
        ["git", "ls-remote", "origin", "refs/heads/main", f"refs/tags/{version}"],
        path,
    )
    remote = {
        ref: sha
        for sha, ref in (line.split() for line in listing.splitlines() if line.strip())
    }
    if (
        remote.get("refs/heads/main") == head
        and remote.get(f"refs/tags/{version}") == head
    ):
        return False
    runner(["git", "push", "origin", "main", version], path)
    return True


def _pages(path: Path, settings: ReleaseSettings, runner: Runner) -> bool:
    endpoint = f"repos/{settings.repository}/pages"
    try:
        runner(["gh", "api", endpoint], path)
    except ReleaseError as err:
        if "404" not in str(err) and "Not Found" not in str(err):
            raise
        runner(["gh", "api", "-X", "POST", endpoint, "-f", "build_type=workflow"], path)
        return True
    return False


def _docker(
    path: Path, settings: ReleaseSettings, version: str, runner: Runner
) -> bool:
    image = settings.docker_image
    if _succeeds(runner, ["docker", "manifest", "inspect", f"{image}:{version}"], path):
        return False
    runner(
        ["docker", "build", "-t", f"{image}:{version}", "-t", f"{image}:latest", "."],
        path,
    )
    runner(["docker", "push", f"{image}:{version}"], path)
    runner(["docker", "push", f"{image}:latest"], path)
    return True


def _github_release(
    path: Path, settings: ReleaseSettings, version: str, archive: Path, runner: Runner
) -> bool:
    repository = settings.repository
    if _succeeds(
        runner, ["gh", "release", "view", version, "--repo", repository], path
    ):
        return False
    runner(
        [
            "gh", "release", "create", version,
            "--repo", repository,
            "--title", f"{settings.repo_name} {version}",
            "--notes-file", f"release-notes/{version}.md",
            str(archive),
        ],
        path,
    )  # fmt: skip
    return True


def publish(
    settings: ReleaseSettings,
    release_path: Path,
    version: str,
    archive: Path,
    runner: Runner = run,
) -> list[str]:
    """Publish a built release; every step already done is skipped.

    Args:
        settings: the release settings of the model.
        release_path: the checkout of the release repository.
        version: the released version.
        archive: the source archive attached to the GitHub release.
        runner: runs git, gh and docker.

    Returns:
        The names of the steps which ran: `commit`, `tag`, `push`, `pages`,
        `docker`, `github_release`.

    Raises:
        ReleaseError: if a command fails or the tag exists at another commit.
    """
    steps = {
        "commit": lambda: _commit(release_path, version, runner),
        "tag": lambda: _tag(release_path, version, runner),
        "push": lambda: _push(release_path, version, runner),
        "pages": lambda: _pages(release_path, settings, runner),
        "docker": lambda: _docker(release_path, settings, version, runner),
        "github_release": lambda: _github_release(
            release_path, settings, version, archive, runner
        ),
    }
    return [name for name, step in steps.items() if step()]
