"""The source repository and the release repository of a model.

A release exports the tracked files of a model package from its source
repository (the monorepo) into a clean checkout of the release repository. The
functions here check both repositories and write the package; every refusal is a
`ReleaseError`. Git is called through a `Runner`, so later steps and the tests
can replace it.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from pkpdsim.model import PKPDModel
from pkpdsim.release.settings import ReleaseSettings

_VERSION = re.compile(r"^(\d+)\.(\d+)\.(\d+)$")
_GITHUB_URL = re.compile(
    r"^(?:https://github\.com/|git@github\.com:|ssh://git@github\.com/)"
    r"(?P<repository>[^/]+/[^/]+?)(?:\.git)?/?$"
)
_STDERR_TAIL = 20


class ReleaseError(Exception):
    """A refusal of the release command, the message is for the user."""


Runner = Callable[[Sequence[str], Path], str]


def run(args: Sequence[str], cwd: Path) -> str:
    """Run a command and return its stdout.

    Args:
        args: the command and its arguments.
        cwd: the working directory.

    Returns:
        The standard output of the command.

    Raises:
        ReleaseError: if the command fails or does not exist; the message has the
            command and the tail of its standard error.
    """
    command = " ".join(str(a) for a in args)
    try:
        result = subprocess.run(
            [str(a) for a in args],
            cwd=cwd,
            check=True,
            text=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as err:
        tail = "\n".join((err.stderr or "").strip().splitlines()[-_STDERR_TAIL:])
        raise ReleaseError(f"'{command}' failed in {cwd}:\n{tail}") from err
    except OSError as err:
        raise ReleaseError(f"'{command}' could not run in {cwd}: {err}") from err
    return result.stdout


@dataclass(frozen=True)
class SourceRepo:
    """The clean source repository of a model.

    Attributes:
        root: the top level folder of the git repository.
        commit: the commit the release is built from.
    """

    root: Path
    commit: str


@dataclass(frozen=True)
class ReleaseRepo:
    """The checkout of the release repository.

    Attributes:
        path: the folder of the checkout.
        last_tag: the highest released version, `None` without a release.
    """

    path: Path
    last_tag: str | None


def source_repo(model: PKPDModel, runner: Runner = run) -> SourceRepo:
    """Find the source repository of a model and check it is clean.

    Args:
        model: the model, its `path` lies in the repository.
        runner: runs the git commands.

    Returns:
        The root and the current commit of the repository.

    Raises:
        ReleaseError: if there are uncommitted or untracked files.
    """
    root = Path(runner(["git", "rev-parse", "--show-toplevel"], model.path).strip())
    if runner(["git", "status", "--porcelain"], root).strip():
        raise ReleaseError(f"commit the changes of {root} first")
    commit = runner(["git", "rev-parse", "HEAD"], root).strip()
    return SourceRepo(root=root, commit=commit)


def clone_url(repository: str) -> str:
    """The URL the release repository is cloned from.

    Args:
        repository: the GitHub repository as `owner/name`.

    Returns:
        The https URL of the repository.
    """
    return f"https://github.com/{repository}.git"


def _is_repository(origin: str, repository: str) -> bool:
    if origin == clone_url(repository):
        return True
    match = _GITHUB_URL.match(origin)
    return bool(match) and match["repository"].lower() == repository.lower()


def _parse_version(version: str) -> tuple[int, int, int]:
    match = _VERSION.match(version)
    if match is None:
        raise ReleaseError(f"'{version}' is not a version 'X.Y.Z'")
    return int(match[1]), int(match[2]), int(match[3])


def _last_tag(tags: Sequence[str]) -> str | None:
    versions = [t for t in tags if _VERSION.match(t)]
    return max(versions, key=_parse_version, default=None)


def release_repo(
    settings: ReleaseSettings,
    path: Path | None,
    source: SourceRepo,
    runner: Runner = run,
) -> ReleaseRepo:
    """Find or clone the release repository and check it is ready.

    Args:
        settings: the release settings of the model.
        path: the checkout, default is a sibling of the source repository named
            like the release repository.
        source: the source repository.
        runner: runs the git commands.

    Returns:
        The checkout and its highest released version.

    Raises:
        ReleaseError: if `origin` is another repository, the branch is not `main`
            or the branch is behind `origin/main`.
    """
    path = path or source.root.parent / settings.repo_name
    if not path.exists():
        runner(
            ["git", "clone", clone_url(settings.repository), str(path)],
            path.parent,
        )
    origin = runner(["git", "remote", "get-url", "origin"], path).strip()
    if not _is_repository(origin, settings.repository):
        raise ReleaseError(
            f"{path} has origin {origin}, not the release repository "
            f"{settings.repository}"
        )
    branch = runner(["git", "rev-parse", "--abbrev-ref", "HEAD"], path).strip()
    if branch != "main":
        raise ReleaseError(f"{path} is on branch '{branch}', not 'main'")
    runner(["git", "fetch", "--tags", "origin"], path)
    behind = runner(["git", "rev-list", "--count", "HEAD..origin/main"], path)
    if int(behind.strip() or "0") > 0:
        raise ReleaseError(f"{path} is behind origin/main, pull it first")
    tags = runner(["git", "tag", "--list"], path).split()
    return ReleaseRepo(path=path, last_tag=_last_tag(tags))


def check_version(version: str, last_tag: str | None) -> None:
    """Check that a version is newer than the last release.

    Args:
        version: the new version `X.Y.Z`.
        last_tag: the last released version, `None` without a release.

    Raises:
        ReleaseError: if a version is malformed or `version` is not greater.
    """
    new = _parse_version(version)
    if last_tag is not None and new <= _parse_version(last_tag):
        raise ReleaseError(
            f"version {version} must be greater than the last release {last_tag}"
        )


def _excluded(relative: PurePosixPath, exclude: Sequence[str]) -> bool:
    for entry in exclude:
        entry_path = PurePosixPath(entry.strip("/"))
        if relative == entry_path or entry_path in relative.parents:
            return True
    return False


def export_package(
    model: PKPDModel,
    settings: ReleaseSettings,
    source: SourceRepo,
    release_path: Path,
    runner: Runner = run,
) -> list[Path]:
    """Write the tracked files of the model package into the release repository.

    The files below `src/<package as folders>` are replaced. The parent folders
    of the package (and `src`) are namespace folders, so an `__init__.py` there is removed
    and none is created.

    Args:
        model: the model, `path` is the folder and `package` the dotted name.
        settings: the release settings, `exclude` lists paths of the model folder
            which are not released (a folder excludes everything below it).
        source: the source repository.
        release_path: the checkout of the release repository.
        runner: runs the git commands.

    Returns:
        The written files, relative to `release_path`.
    """
    listing = runner(["git", "ls-files", "-z", "--", "."], model.path)
    files = [
        PurePosixPath(name)
        for name in listing.split("\0")
        if name and not _excluded(PurePosixPath(name), settings.exclude)
    ]
    package_dirs = Path(*model.package.split("."))
    src = release_path / "src"
    target = src / package_dirs
    if target.exists():
        shutil.rmtree(target)
    for parent in package_dirs.parents:
        stale = src / parent / "__init__.py"
        if stale.exists():
            stale.unlink()
    written: list[Path] = []
    for relative in files:
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(model.path / relative, destination)
        written.append(destination.relative_to(release_path))
    return written
