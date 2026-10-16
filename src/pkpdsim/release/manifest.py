"""The manifest of a release and the removal of what the release did not build.

The build writes `release.json` (the version, the commit of the source, the
library versions and a hash of every generated file). A later run finds a build
which is not committed yet and publishes it as it is when the manifest still
matches, so a reviewed `--no-publish` build needs no new simulations. `prune`
makes the release repository hold only what the command generates.
"""

from __future__ import annotations

import hashlib
import json
import os
import tomllib
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path, PurePosixPath

from pkpdsim.release.source import ReleaseError, Runner, SourceRepo, run

MANIFEST_FILE = "release.json"
_RELEASE_NOTES = "release-notes"


@dataclass(frozen=True)
class Manifest:
    """What a build of the release repository produced.

    Attributes:
        version: the released version.
        source_commit: the commit of the source repository.
        pkpdsim: the version of pkpdsim in the lock of the release.
        sbmlsim: the version of sbmlsim in the lock of the release.
        sbmlutils: the version of sbmlutils in the lock of the release.
        files: the sha256 of every generated file by its POSIX path relative to
            the release repository.
    """

    version: str
    source_commit: str
    pkpdsim: str
    sbmlsim: str
    sbmlutils: str
    files: dict[str, str]

    def to_json(self) -> str:
        """Serialise the manifest.

        Returns:
            The JSON text with sorted keys and a final newline.
        """
        return json.dumps(asdict(self), indent=2, sort_keys=True) + "\n"

    @staticmethod
    def read(path: Path) -> Manifest | None:
        """Read a manifest.

        Args:
            path: the `release.json` file.

        Returns:
            The manifest, `None` if the file is missing.

        Raises:
            ReleaseError: if the file is not a manifest.
        """
        if not path.is_file():
            return None
        try:
            return Manifest(**json.loads(path.read_text(encoding="utf-8")))
        except (ValueError, TypeError) as err:
            raise ReleaseError(f"'{path}' is not a release manifest: {err}") from err


def _relative(release_path: Path, path: Path) -> str:
    """Normalise a path of a stage to a POSIX path relative to the release."""
    if path.is_absolute():
        try:
            path = path.relative_to(release_path)
        except ValueError as err:
            raise ReleaseError(
                f"'{path}' is not inside the release repository {release_path}."
            ) from err
    return PurePosixPath(*path.parts).as_posix()


def _generated(release_path: Path, written: Iterable[Path]) -> set[str]:
    """The relative POSIX paths of all files written; folders are expanded."""
    files: set[str] = set()
    for path in written:
        relative = _relative(release_path, path)
        absolute = release_path / relative
        if absolute.is_dir():
            for root, _, names in os.walk(absolute):
                for name in names:
                    files.add(_relative(release_path, Path(root) / name))
        else:
            files.add(relative)
    return files


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _library_versions(release_path: Path) -> dict[str, str]:
    lock_path = release_path / "uv.lock"
    if not lock_path.is_file():
        raise ReleaseError(f"The release has no 'uv.lock' at {lock_path}.")
    lock = tomllib.loads(lock_path.read_text(encoding="utf-8"))
    return {p["name"]: p["version"] for p in lock.get("package", []) if "version" in p}


def write_manifest(
    release_path: Path,
    version: str,
    source: SourceRepo,
    written: Iterable[Path],
) -> Manifest:
    """Write `release.json` with a hash of every generated file.

    Args:
        release_path: the release repository.
        version: the released version.
        source: the source repository the release is built from.
        written: the files (or folders) the stages wrote, absolute or relative to
            the release repository.

    Returns:
        The manifest which was written.

    Raises:
        ReleaseError: if `uv.lock` lacks pkpdsim, sbmlsim or sbmlutils.
    """
    libraries = _library_versions(release_path)
    missing = [n for n in ("pkpdsim", "sbmlsim", "sbmlutils") if n not in libraries]
    if missing:
        raise ReleaseError(f"'uv.lock' of the release lacks {', '.join(missing)}.")
    files = {
        relative: _hash(release_path / relative)
        for relative in sorted(_generated(release_path, written))
        if relative != MANIFEST_FILE
    }
    manifest = Manifest(
        version=version,
        source_commit=source.commit,
        pkpdsim=libraries["pkpdsim"],
        sbmlsim=libraries["sbmlsim"],
        sbmlutils=libraries["sbmlutils"],
        files=files,
    )
    (release_path / MANIFEST_FILE).write_text(manifest.to_json(), encoding="utf-8")
    return manifest


def _status_paths(release_path: Path, runner: Runner) -> set[str]:
    out = runner(
        ["git", "status", "--porcelain", "-z", "--untracked-files=all"], release_path
    )
    paths: set[str] = set()
    entries = iter(out.split("\0"))
    for entry in entries:
        if len(entry) < 4:
            continue
        paths.add(entry[3:])
        if entry[0] in "RC":  # the source of a rename or copy follows
            paths.add(next(entries, ""))
    return paths


def matches(
    release_path: Path,
    version: str,
    source_commit: str,
    pkpdsim_version: str,
    runner: Runner = run,
) -> bool:
    """Tell whether the release repository holds the build of this release.

    Args:
        release_path: the release repository.
        version: the version to release.
        source_commit: the commit of the source repository.
        pkpdsim_version: the version of pkpdsim which runs the release.
        runner: runs git.

    Returns:
        `True` if the manifest exists, version, commit and pkpdsim version are
        equal, every listed file hashes the same and git status lists no file
        outside the manifest (besides `release.json`).
    """
    manifest = Manifest.read(release_path / MANIFEST_FILE)
    if manifest is None:
        return False
    if (manifest.version, manifest.source_commit, manifest.pkpdsim) != (
        version,
        source_commit,
        pkpdsim_version,
    ):
        return False
    for relative, digest in manifest.files.items():
        path = release_path / relative
        if not path.is_file() or _hash(path) != digest:
            return False
    allowed = {*manifest.files, MANIFEST_FILE}
    return _status_paths(release_path, runner) <= allowed


def _ignored(release_path: Path, runner: Runner) -> tuple[set[str], set[str]]:
    """The untracked ignored files and folders (the latter with a trailing /)."""
    out = runner(
        [
            "git",
            "ls-files",
            "--others",
            "--ignored",
            "--exclude-standard",
            "--directory",
            "-z",
        ],
        release_path,
    )
    entries = {e for e in out.split("\0") if e}
    files = {e for e in entries if not e.endswith("/")}
    return files, entries - files


def prune(
    release_path: Path, written: Iterable[Path], runner: Runner = run
) -> list[Path]:
    """Remove every file of the release repository the command did not generate.

    Kept are `.git/`, ignored files (`.gitignore`), `release-notes/*.md`,
    `release.json` and `written`. Folders left empty are removed.

    Args:
        release_path: the release repository.
        written: the files (or folders) the stages wrote.
        runner: runs git.

    Returns:
        The removed files and folders (absolute).
    """
    keep = _generated(release_path, written) | {MANIFEST_FILE}
    ignored_files, ignored_dirs = _ignored(release_path, runner)
    removed: list[Path] = []
    for root, dirs, names in os.walk(release_path):
        here = Path(root)
        relative_root = here.relative_to(release_path)
        dirs[:] = [
            d
            for d in dirs
            if not (relative_root == Path() and d == ".git")
            and (relative_root / d).as_posix() + "/" not in ignored_dirs
        ]
        for name in names:
            relative = (relative_root / name).as_posix()
            if (
                relative in keep
                or relative in ignored_files
                or (relative_root.as_posix() == _RELEASE_NOTES and name.endswith(".md"))
            ):
                continue
            (here / name).unlink()
            removed.append(here / name)
    for root, _, _ in os.walk(release_path, topdown=False):
        here = Path(root)
        relative = here.relative_to(release_path)
        if relative == Path() or relative.parts[0] == ".git":
            continue
        prefix = relative.as_posix() + "/"
        if prefix in ignored_dirs or any(i.startswith(prefix) for i in ignored_files):
            continue
        if not any(here.iterdir()):
            here.rmdir()
            removed.append(here)
    return removed
