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

from sbmlutils.console import console

from pkpdsim.release.manifest import MANIFEST_FILE
from pkpdsim.release.settings import ReleaseSettings
from pkpdsim.release.source import ReleaseError, Runner, run, stage

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


def _ok(runner: Runner, args: list[str], cwd: Path) -> bool:
    try:
        runner(args, cwd)
    except ReleaseError:
        return False
    return True


def _lock_problems(release_path: Path) -> list[str]:
    problems = []
    if not (release_path / "uv.lock").is_file():
        problems.append("no uv.lock")
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
    return problems


def _tool_problems(release_path: Path, runner: Runner) -> list[str]:
    problems = []
    if not _ok(runner, ["docker", "info"], release_path):
        problems.append("docker is not available (docker info failed)")
    if not _ok(runner, ["gh", "auth", "status"], release_path):
        problems.append("gh is not logged in (gh auth status failed)")
    return problems


def _refuse(problems: list[str]) -> None:
    if problems:
        raise ReleaseError("cannot publish:\n- " + "\n- ".join(problems))


def check_tools(release_path: Path, runner: Runner = run) -> None:
    """Check that the tools of the publishing steps are ready.

    The check needs no build, so a release runs it before the build.

    Args:
        release_path: the checkout of the release repository.
        runner: runs the tool checks.

    Raises:
        ReleaseError: if Docker is not running or `gh` is not logged in; one
            message lists both.
    """
    _refuse(_tool_problems(release_path, runner))


def check_lock(release_path: Path) -> None:
    """Check that the built project depends on released packages only.

    Args:
        release_path: the checkout of the release repository, after the build.

    Raises:
        ReleaseError: if there is no `uv.lock`, the lock has a package which is
            not from PyPI or the project uses `tool.uv.sources`; one message
            lists everything.
    """
    _refuse(_lock_problems(release_path))


def head_is_release(path: Path, version: str, runner: Runner = run) -> bool:
    """Tell whether the release of a version is committed at `HEAD`.

    Args:
        path: the checkout of the release repository.
        version: the released version.
        runner: runs git.

    Returns:
        `True` if the working tree is clean and the `release.json` of `HEAD` is
        the manifest of `version`.
    """
    if runner(["git", "status", "--porcelain"], path).strip():
        return False
    try:
        text = runner(["git", "show", f"HEAD:{MANIFEST_FILE}"], path)
        return json.loads(text).get("version") == version
    except (ReleaseError, ValueError, AttributeError):
        return False


def tag_commit(path: Path, version: str, runner: Runner = run) -> str | None:
    """Get the commit of the tag of a version.

    Args:
        path: the checkout of the release repository.
        version: the version, the name of the tag.
        runner: runs git.

    Returns:
        The commit of `refs/tags/<version>`, `None` without the tag.
    """
    if not runner(["git", "tag", "--list", version], path).strip():
        return None
    return runner(["git", "rev-list", "-n", "1", f"refs/tags/{version}"], path).strip()


def _commit(path: Path, version: str, runner: Runner) -> bool:
    if head_is_release(path, version, runner):
        return False
    runner(["git", "add", "-A"], path)
    runner(["git", "commit", "-m", f"Release {version}"], path)
    return True


def _tag(path: Path, version: str, runner: Runner) -> bool:
    tagged = tag_commit(path, version, runner)
    if tagged is not None:
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
    runner(["git", "push", "origin", "main", version], path, stream=True)
    return True


def _pages(path: Path, settings: ReleaseSettings, runner: Runner) -> bool:
    endpoint = f"repos/{settings.repository}/pages"
    try:
        runner(["gh", "api", endpoint], path)
    except ReleaseError as err:
        if "HTTP 404" not in str(err) and "Not Found" not in str(err):
            raise
        runner(["gh", "api", "-X", "POST", endpoint, "-f", "build_type=workflow"], path)
        return True
    return False


def _manifest(runner: Runner, reference: str, path: Path) -> str | None:
    try:
        return runner(["docker", "manifest", "inspect", reference], path)
    except ReleaseError:
        return None


def _docker(
    path: Path, settings: ReleaseSettings, version: str, runner: Runner
) -> bool:
    image = settings.docker_image
    versioned, latest = f"{image}:{version}", f"{image}:latest"
    remote_version = _manifest(runner, versioned, path)
    remote_latest = _manifest(runner, latest, path)
    if remote_version is not None and remote_version == remote_latest:
        return False
    if not _ok(runner, ["docker", "image", "inspect", versioned], path):
        runner(
            ["docker", "build", "-t", versioned, "-t", latest, "."], path, stream=True
        )
    if remote_version is None:
        runner(["docker", "push", versioned], path, stream=True)
    runner(["docker", "push", latest], path, stream=True)
    return True


def _github_release(
    path: Path, settings: ReleaseSettings, version: str, archive: Path, runner: Runner
) -> bool:
    repository = settings.repository
    try:
        view = runner(
            [
                "gh",
                "release",
                "view",
                version,
                "--repo",
                repository,
                "--json",
                "assets",
            ],
            path,
        )
    except ReleaseError:
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
    assets = {a.get("name") for a in json.loads(view or "{}").get("assets", [])}
    if archive.name in assets:
        return False
    runner(
        ["gh", "release", "upload", version, str(archive), "--repo", repository,
         "--clobber"],
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
        "commit": ("Commit", lambda: _commit(release_path, version, runner)),
        "tag": ("Tag", lambda: _tag(release_path, version, runner)),
        "push": ("Push", lambda: _push(release_path, version, runner)),
        "pages": ("GitHub Pages", lambda: _pages(release_path, settings, runner)),
        "docker": (
            "Docker image",
            lambda: _docker(release_path, settings, version, runner),
        ),
        "github_release": (
            "GitHub release",
            lambda: _github_release(release_path, settings, version, archive, runner),
        ),
    }
    done = []
    for name, (label, step) in steps.items():
        stage(label)
        if step():
            done.append(name)
        else:
            console.print(f"{label}: done already")
    return done
