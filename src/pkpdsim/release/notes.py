"""The release notes of a release repository.

The notes of a release are generated from the difference between the previous
release, the last tag of the release repository, and the new build in its working
tree: the core libraries with a summary of their own release notes, the locked
libraries, the study data, the models, the experiments and the adopted fits. A
previous release of the old layout lacks parts (e.g. the fit records); what cannot
be compared is stated in the notes, not omitted.

Git and `gh` are called through a `Runner`; the release notes of the core
libraries come from a `Fetch`, so the tests run offline.
"""

from __future__ import annotations

import functools
import json
import logging
import math
import re
import tomllib
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any, Literal

import libsbml

from pkpdsim import __version__
from pkpdsim.fit_report import fit_subset
from pkpdsim.model import PKPDModel
from pkpdsim.release.build import BuildPaths
from pkpdsim.release.project import ProjectContext, render
from pkpdsim.release.source import ReleaseError, ReleaseRepo, Runner, SourceRepo, run

logger = logging.getLogger(__name__)

CORE_LIBRARIES = ("pkpdsim", "sbmlsim", "sbmlutils")
"""The libraries whose release notes are summarized."""

MAX_BULLETS = 5
"""The maximal number of bullets of the summary of a core library."""

NOTES_DIR = "release-notes"
"""The folder of the release notes in the release repository."""

_HASH_CHUNK = 200
_BULLET = re.compile(r"^- (?P<text>.+)$")
_SENTENCE_END = re.compile(
    r"(?<!\be\.g)(?<!\bi\.e)(?<!\bcf)(?<!\bvs)[.!?](?:\*\*|\*|\)|`)?(?=\s)"
)
"""The end of a sentence: a full stop, not of an abbreviation, before a space; a
closing emphasis, parenthesis or backtick belongs to the sentence."""
_RELEASE_VERSION = re.compile(r"^v?(?P<version>\d+(?:\.\d+)*)$")
_ELEMENTS = ("species", "reactions", "parameters")
_REQUIREMENT_NAME = re.compile(r"^\s*(?P<name>[A-Za-z0-9][A-Za-z0-9._-]*)")

Fetch = Callable[[str], dict[str, str]]
"""Get the release notes of a library, the body by version, from its repository."""


@dataclass(frozen=True)
class LibraryChange:
    """The changed version of a locked library.

    Attributes:
        name: the name of the package.
        old: the version in the previous release, `None` if it was added.
        new: the version in the new release, `None` if it was removed.
    """

    name: str
    old: str | None
    new: str | None


@dataclass(frozen=True)
class CoreSummary:
    """The summary of the release notes of a core library between two versions.

    Attributes:
        name: the name of the library.
        old: the version in the previous release, `None` if it was not used.
        new: the version in the new release, `None` if it is no longer used.
        bullets: the first sentences of the top level bullets of the releases.
        fetched: whether the release notes of the library could be fetched.
        total: the number of bullets of the releases, more than `bullets` if
            they were cut to `MAX_BULLETS`.
    """

    name: str
    old: str | None
    new: str | None
    bullets: list[str]
    fetched: bool = True
    total: int = 0


@dataclass(frozen=True)
class DataChanges:
    """The studies of the data, `<substance>/<study>`, which changed.

    Attributes:
        added: the studies which are new.
        removed: the studies which are gone.
        changed: the studies with added, removed or changed files.
    """

    added: list[str] = field(default_factory=list)
    removed: list[str] = field(default_factory=list)
    changed: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class ModelChange:
    """The changes of an SBML model of the release.

    Attributes:
        file: the name of the SBML file in `models/`.
        added: the number of species, reactions and parameters which are new.
        removed: the number of species, reactions and parameters which are gone.
        changed_parameters: the parameters with another value, `(id, old, new)`.
        status: whether the model is new, gone or changed.
    """

    file: str
    added: dict[str, int]
    removed: dict[str, int]
    changed_parameters: list[tuple[str, float | None, float | None]]
    status: Literal["added", "removed", "changed"] = "changed"


@dataclass(frozen=True)
class FitChange:
    """The change of the adopted fit of a fit problem.

    Attributes:
        old_id: the adopted fit of the previous release, `None` without a record.
        new_id: the adopted fit of the new release, `None` if the problem has none.
        parameters: the fitted parameters which moved, `(id, old, new)`; `None`
            for a parameter which is not fitted in one of the fits.
    """

    old_id: str | None
    new_id: str | None
    parameters: list[tuple[str, float | None, float | None]]

    @property
    def subset(self) -> str:
        """The fit problem of the fits."""
        fit_id = self.new_id or self.old_id
        assert fit_id is not None
        return fit_subset(fit_id)


def _tree(release: ReleaseRepo, runner: Runner) -> dict[str, str]:
    """Get the git blob hash of every file of the previous release by its path."""
    if release.last_tag is None:
        return {}
    listing = runner(["git", "ls-tree", "-r", "-z", release.last_tag], release.path)
    tree: dict[str, str] = {}
    for entry in listing.split("\0"):
        if not entry:
            continue
        meta, path = entry.split("\t", 1)
        _, kind, blob = meta.split()
        if kind == "blob":
            tree[path] = blob
    return tree


def _show(release: ReleaseRepo, path: str, runner: Runner) -> str:
    """Get the content of a file of the previous release."""
    return runner(["git", "show", f"{release.last_tag}:{path}"], release.path)


def _hash_objects(
    release: ReleaseRepo, paths: Sequence[str], runner: Runner
) -> list[str]:
    """Get the git blob hashes of files of the working tree."""
    hashes: list[str] = []
    for start in range(0, len(paths), _HASH_CHUNK):
        chunk = list(paths[start : start + _HASH_CHUNK])
        output = runner(["git", "hash-object", "--", *chunk], release.path)
        hashes.extend(output.split())
    return hashes


def _worktree_files(release: ReleaseRepo, prefix: str, runner: Runner) -> list[str]:
    """Get the files below a folder of the working tree which are not ignored."""
    listing = runner(
        [
            "git",
            "ls-files",
            "-z",
            "--cached",
            "--others",
            "--exclude-standard",
            "--",
            prefix,
        ],
        release.path,
    )
    names = {name for name in listing.split("\0") if name}
    return sorted(name for name in names if (release.path / name).is_file())


def _locked(lock: str) -> dict[str, str]:
    """Get the versions of the packages of a lock file, without the project."""
    versions: dict[str, set[str]] = {}
    for package in tomllib.loads(lock).get("package", []):
        source = package.get("source", {})
        if source.get("editable") == "." or source.get("virtual") == ".":
            continue
        versions.setdefault(package["name"], set()).add(package["version"])
    return {name: ", ".join(sorted(found)) for name, found in versions.items()}


def library_changes(old_lock: str | None, new_lock: str) -> list[LibraryChange]:
    """Compare the locked versions of the libraries of two releases.

    Args:
        old_lock: the `uv.lock` of the previous release, `None` without one.
        new_lock: the `uv.lock` of the new release.

    Returns:
        The libraries which were added, removed or changed their version, sorted
        by name; empty without a previous lock file.
    """
    if old_lock is None:
        return []
    old, new = _locked(old_lock), _locked(new_lock)
    return [
        LibraryChange(name, old.get(name), new.get(name))
        for name in sorted(old.keys() | new.keys())
        if old.get(name) != new.get(name)
    ]


def _version_key(version: str) -> tuple[int, ...] | None:
    match = _RELEASE_VERSION.match(version.strip())
    if match is None:
        return None
    return tuple(int(part) for part in match["version"].split("."))


def _in_range(version: str, old: str, new: str) -> bool:
    """Check that a version is after `old` and up to `new`."""
    key, low, high = _version_key(version), _version_key(old), _version_key(new)
    return (
        key is not None and low is not None and high is not None and low < key <= high
    )


def _first_sentence(text: str) -> str:
    text = text.strip()
    end = _SENTENCE_END.search(text)
    return text[: end.end()] if end else text


def _bullets(body: str) -> list[str]:
    """Get the first sentences of the top level bullets of release notes."""
    bullets = []
    for line in body.splitlines():
        match = _BULLET.match(line.rstrip())
        if match:
            bullets.append(_first_sentence(match["text"]))
    return bullets


def github_releases(repository: str, runner: Runner = run) -> dict[str, str]:
    """Get the release notes of a library of matthiaskoenig from GitHub.

    One call of `gh api` lists all releases, so the releases between two
    versions need no further call. A leading `v` of a tag is removed.

    Args:
        repository: the name of the repository, e.g. `sbmlsim`.
        runner: runs `gh`.

    Returns:
        The body of the release notes by version.

    Raises:
        ReleaseError: if `gh` is missing or fails, e.g. without network.
    """
    output = runner(
        [
            "gh",
            "api",
            f"repos/matthiaskoenig/{repository}/releases",
            "--paginate",
            "--jq",
            '.[] | {tag: .tag_name, body: (.body // "")}',
        ],
        Path.cwd(),
    )
    releases: dict[str, str] = {}
    for line in output.splitlines():
        if line.strip():
            release = json.loads(line)
            releases[release["tag"].removeprefix("v")] = release["body"] or ""
    return releases


def core_summaries(
    changes: Iterable[LibraryChange], fetch: Fetch = github_releases
) -> list[CoreSummary]:
    """Summarize the release notes of the core libraries which changed.

    For every release after the old up to the new version, newest first, the
    first sentences of the top level bullets (`- `) are taken, at most
    `MAX_BULLETS` per library. A library whose release notes cannot be fetched
    gets a summary without bullets and `fetched=False`, it never fails the
    release. A library which was added or removed gets a summary without
    bullets.

    Args:
        changes: the changed libraries, see `library_changes`.
        fetch: gets the release notes of a library by version.

    Returns:
        The summaries of the core libraries which changed.
    """
    by_name = {change.name: change for change in changes}
    summaries = []
    for name in CORE_LIBRARIES:
        change = by_name.get(name)
        if change is None:
            continue
        if change.old is None or change.new is None:
            summaries.append(CoreSummary(name, change.old, change.new, []))
            continue
        try:
            releases = fetch(name)
        except (ReleaseError, ValueError, KeyError) as err:
            logger.warning(
                "The release notes of %s could not be fetched: %s", name, err
            )
            summaries.append(CoreSummary(name, change.old, change.new, [], False))
            continue
        versions = sorted(
            (v for v in releases if _in_range(v, change.old, change.new)),
            key=lambda v: _version_key(v) or (),
            reverse=True,
        )
        bullets = [
            bullet for version in versions for bullet in _bullets(releases[version])
        ]
        summaries.append(
            CoreSummary(
                name,
                change.old,
                change.new,
                bullets[:MAX_BULLETS],
                total=len(bullets),
            )
        )
    return summaries


def _studies(files: dict[str, str]) -> dict[str, set[tuple[str, str]]]:
    """Group files of the data folder, `path: hash`, by their study."""
    studies: dict[str, set[tuple[str, str]]] = {}
    for path, blob in files.items():
        parts = PurePosixPath(path).parts
        if len(parts) >= 3:
            study = f"{parts[0]}/{parts[1]}"
            studies.setdefault(study, set()).add(("/".join(parts[2:]), blob))
    return studies


def data_changes(
    release: ReleaseRepo, data_dir: str, runner: Runner = run
) -> DataChanges:
    """Compare the studies of the data of the previous release and the new build.

    A study is a folder `<substance>/<study>` of the data folder; it changed if a
    file was added, removed or changed. The files of the previous release are
    compared by the git blob hashes at the last tag with the hashes git gives the
    files of the working tree (`git hash-object`), ignored files are skipped.

    Args:
        release: the release repository.
        data_dir: the data folder, relative to the release repository.
        runner: runs git.

    Returns:
        The added, removed and changed studies.
    """
    prefix = data_dir.strip("/") + "/"
    old = _studies(
        {
            path.removeprefix(prefix): blob
            for path, blob in _tree(release, runner).items()
            if path.startswith(prefix)
        }
    )
    paths = _worktree_files(release, prefix, runner)
    hashes = _hash_objects(release, paths, runner)
    new = _studies(
        {
            path.removeprefix(prefix): blob
            for path, blob in zip(paths, hashes, strict=True)
        }
    )
    return DataChanges(
        added=sorted(new.keys() - old.keys()),
        removed=sorted(old.keys() - new.keys()),
        changed=sorted(s for s in old.keys() & new.keys() if old[s] != new[s]),
    )


def _is_experiment(name: str) -> bool:
    """Check a folder of the simulation report is an experiment, not `_figures`."""
    return not name.startswith(("_", "."))


def _experiments(release: ReleaseRepo, runner: Runner) -> tuple[set[str], set[str]]:
    """Get the experiments of the previous release and of the new build."""
    old = set()
    for path in _tree(release, runner):
        parts = PurePosixPath(path).parts
        if len(parts) >= 3 and parts[0] == BuildPaths.SIMULATIONS:
            old.add(parts[1])
    folder = release.path / BuildPaths.SIMULATIONS
    new = {p.name for p in folder.iterdir() if p.is_dir()} if folder.is_dir() else set()
    return {e for e in old if _is_experiment(e)}, {e for e in new if _is_experiment(e)}


def experiment_changes(
    release: ReleaseRepo, runner: Runner = run
) -> tuple[list[str], list[str]]:
    """Compare the experiments of the previous release and the new build.

    The experiments are the folders of `simulations/`, without the folders whose
    name starts with `_` (e.g. `_figures`).

    Args:
        release: the release repository.
        runner: runs git.

    Returns:
        The added and the removed experiments.
    """
    old, new = _experiments(release, runner)
    return sorted(new - old), sorted(old - new)


def _model_elements(text: str | None, name: str) -> dict[str, dict[str, Any]]:
    """Get the species, reactions and parameters of an SBML model by id."""
    elements: dict[str, dict[str, Any]] = {kind: {} for kind in _ELEMENTS}
    if text is None:
        return elements
    model = libsbml.readSBMLFromString(text).getModel()
    if model is None:
        raise ReleaseError(f"The model '{name}' of the release is no SBML model.")
    elements["species"] = {s.getId(): None for s in model.getListOfSpecies()}
    elements["reactions"] = {r.getId(): None for r in model.getListOfReactions()}
    elements["parameters"] = {p.getId(): _value(p) for p in model.getListOfParameters()}
    return elements


def _value(parameter: libsbml.Parameter) -> float | None:
    """Get the value of a parameter, `None` if it is not set or not a number."""
    if not parameter.isSetValue() or math.isnan(parameter.getValue()):
        return None
    return parameter.getValue()


def model_changes(release: ReleaseRepo, runner: Runner = run) -> list[ModelChange]:
    """Compare the SBML models of the previous release and the new build.

    The models are the files `models/*.xml`; a model whose file is unchanged is
    not listed.

    Args:
        release: the release repository.
        runner: runs git.

    Returns:
        The models which are new, gone or changed, sorted by file: the number of
        species, reactions and parameters added and removed and the parameters
        with another value.

    Raises:
        ReleaseError: if a file is no SBML model.
    """
    old = {
        PurePosixPath(path).name: blob
        for path, blob in _tree(release, runner).items()
        if PurePosixPath(path).parent == PurePosixPath(BuildPaths.MODELS)
        and path.endswith(".xml")
    }
    folder = release.path / BuildPaths.MODELS
    names = sorted(p.name for p in folder.glob("*.xml")) if folder.is_dir() else []
    paths = [f"{BuildPaths.MODELS}/{name}" for name in names]
    new = dict(zip(names, _hash_objects(release, paths, runner), strict=True))

    changes = []
    for name in sorted(old.keys() | new.keys()):
        if old.get(name) == new.get(name):
            continue
        old_text = (
            _show(release, f"{BuildPaths.MODELS}/{name}", runner)
            if name in old
            else None
        )
        new_text = (folder / name).read_text("utf-8") if name in new else None
        before = _model_elements(old_text, name)
        after = _model_elements(new_text, name)
        old_parameters, new_parameters = before["parameters"], after["parameters"]
        changes.append(
            ModelChange(
                file=name,
                added={k: len(after[k].keys() - before[k].keys()) for k in _ELEMENTS},
                removed={k: len(before[k].keys() - after[k].keys()) for k in _ELEMENTS},
                changed_parameters=[
                    (pid, old_parameters[pid], new_parameters[pid])
                    for pid in sorted(old_parameters.keys() & new_parameters.keys())
                    if old_parameters[pid] != new_parameters[pid]
                ],
                status=(
                    "added"
                    if old_text is None
                    else "removed"
                    if new_text is None
                    else "changed"
                ),
            )
        )
    return changes


def _best_values(text: str, name: str) -> dict[str, float]:
    """Get the values of the best parameter set of a `parameters.json`.

    Raises:
        ReleaseError: if the record is no JSON or has no parameter set with values.
    """
    try:
        return dict(json.loads(text)["sets"][0]["values"])
    except (ValueError, KeyError, IndexError, TypeError) as err:
        raise ReleaseError(
            f"The fit record {name} is malformed, it needs the values of a "
            f"parameter set in 'sets': {err!r}"
        ) from err


def _fits_prefix(model: PKPDModel) -> PurePosixPath:
    """Get the folder of the adopted fits in the release repository."""
    try:
        relative = model.fits_dir.relative_to(model.path)
    except ValueError as err:
        raise ReleaseError(
            f"The fits of '{model.name}' at {model.fits_dir} lie outside of the "
            f"model package {model.path} and are not released."
        ) from err
    return PurePosixPath("src", *model.package.split("."), *relative.parts)


def fit_changes(
    release: ReleaseRepo, model: PKPDModel, runner: Runner = run
) -> list[FitChange]:
    """Compare the adopted fits of the previous release and the new build.

    The fits are matched by their fit problem (`fit_subset`); a fit problem
    without a record in the previous release gives a change with `old_id=None`.
    The values are the best parameter set of the `parameters.json` of a record.

    Args:
        release: the release repository.
        model: the model, `fits` are the adopted fits of the new release.
        runner: runs git.

    Returns:
        The changed fits, sorted by fit problem.

    Raises:
        ReleaseError: if the record of an adopted fit lacks `parameters.json`.
    """
    prefix = _fits_prefix(model)
    old: dict[str, str] = {}
    for path in sorted(_tree(release, runner)):
        relative = PurePosixPath(path)
        if relative.parent.parent != prefix or relative.name != "parameters.json":
            continue
        try:
            old[fit_subset(relative.parent.name)] = relative.parent.name
        except ValueError:
            continue
    new = {fit_subset(fit_id): fit_id for fit_id in model.fits}

    changes = []
    for subset in sorted(old.keys() | new.keys()):
        old_id, new_id = old.get(subset), new.get(subset)
        if old_id == new_id:
            continue
        before: dict[str, float] = {}
        after: dict[str, float] = {}
        if old_id is not None:
            name = f"{prefix}/{old_id}/parameters.json"
            before = _best_values(
                _show(release, name, runner), f"{release.last_tag}:{name}"
            )
        if new_id is not None:
            record = release.path / prefix / new_id / "parameters.json"
            if not record.is_file():
                raise ReleaseError(
                    f"The adopted fit '{new_id}' has no record {record}."
                )
            after = _best_values(record.read_text("utf-8"), str(record))
        moved = (
            [
                (pid, before.get(pid), after.get(pid))
                for pid in sorted(before.keys() | after.keys())
                if before.get(pid) != after.get(pid)
            ]
            if old_id is not None and new_id is not None
            else []
        )
        changes.append(FitChange(old_id, new_id, moved))
    return changes


def _lock_text(
    release: ReleaseRepo, tree: dict[str, str], runner: Runner
) -> str | None:
    """Get the lock file of the previous release, `None` without one."""
    return _show(release, "uv.lock", runner) if "uv.lock" in tree else None


def dependency_names(dependencies: Iterable[str]) -> set[str]:
    """Get the normalized names of the packages of requirements.

    Args:
        dependencies: requirements as in `[project].dependencies`, e.g.
            `"pkpdsim[docs]>=0.7"`.

    Returns:
        The names, lower case with runs of `-`, `_` and `.` as `-` (as in `uv.lock`).
    """
    names = set()
    for dependency in dependencies:
        match = _REQUIREMENT_NAME.match(dependency)
        if match:
            names.add(re.sub(r"[-_.]+", "-", match["name"]).lower())
    return names


def _number(value: float | None) -> str:
    """Format a value of a parameter for the release notes."""
    return "none" if value is None else f"{value:.6g}"


def write_release_notes(
    context: ProjectContext,
    release: ReleaseRepo,
    model: PKPDModel,
    source: SourceRepo,
    runner: Runner = run,
    fetch: Fetch | None = None,
) -> Path:
    """Write the release notes of the new release into the release repository.

    Args:
        context: the project context, with the version of the release.
        release: the release repository with the new build in its working tree
            and the previous release as `last_tag`.
        model: the model which is released.
        source: the source repository the release is made from.
        runner: runs git.
        fetch: gets the release notes of a core library by version, by default
            from GitHub with `gh` (see `github_releases`).

    Returns:
        The written file `release-notes/<version>.md`.

    Raises:
        ReleaseError: if the new build has no `uv.lock`.
    """
    fetch = fetch or functools.partial(github_releases, runner=runner)
    tree = _tree(release, runner)
    lock_path = release.path / "uv.lock"
    if not lock_path.is_file():
        raise ReleaseError(f"The release repository has no 'uv.lock' at {lock_path}.")
    new_lock = lock_path.read_text("utf-8")
    old_lock = _lock_text(release, tree, runner)
    libraries = library_changes(old_lock, new_lock)
    locked = _locked(new_lock)
    direct = dependency_names(context.dependencies)
    data_dir = PurePosixPath("src", *model.package.split("."), "data")
    added_experiments, removed_experiments = experiment_changes(release, runner)
    extra = {
        "commit": source.commit[:10],
        "pkpdsim_version": __version__,
        "previous": release.last_tag,
        "old_lock": old_lock is not None,
        "core_versions": [(n, locked[n]) for n in CORE_LIBRARIES if n in locked],
        "core": core_summaries(libraries, fetch),
        "libraries": libraries,
        "direct_libraries": [c for c in libraries if c.name in direct],
        "data": data_changes(release, str(data_dir), runner),
        "models": model_changes(release, runner),
        "experiments_added": added_experiments,
        "experiments_removed": removed_experiments,
        "fits": fit_changes(release, model, runner),
        "number": _number,
    }
    path = release.path / NOTES_DIR / f"{context.version}.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render("release-notes.md.j2", context, extra), encoding="utf-8")
    return path
