"""Backfill the records of the studies a repository of models collected.

The collection of a study writes its record `.<study>_study.json`
(`pkpdsim.data.StudyRecord`) since pkpdsim 0.19: the day of the first
collection, which the temporal validation compares with the day of an
adopted fit, and the provenance of the tables. The studies collected before
have tables and no record, so their day of collection is unknown. This
one-off script dates them from the git history of the repository of the
models (pkdb_models).

The day is an upper bound of the day the model first had the study, the
earliest of

- the first commit which added a table of the study to its folder, following
  the renames of the folders (a model which moved keeps its days) and keeping
  the earliest day of a folder whose tables were replaced;
- the first commit which mentioned the study in a python file of the model
  (a folder of the name of the model in its path), since a model read the
  tables of many studies from a checkout of pkdb_data before it kept them in
  its own folder, and an experiment of a study is written after its data;
- the day of every adopted fit of the model whose record was committed with
  the study present (the commit of the day above is an ancestor of the commit
  which added the fit, or the same): such a fit is taken to have known the
  study, since a model imported with its data and its fits in one commit, or
  a fit adopted after the study arrived, gives no day before the fit.

The days are author dates, which a rebase keeps. A study is thereby new for a
fit only if the history shows the fit in the repository before any trace of
the study. A model whose history lies in another repository has the day of
its import, so a study can be dated later than it was collected but never new
for a fit imported with it.

A record written here has the dating `"git"` and the commit of its day, the
day of the latest commit which changed a table as the day of the latest
collection, the release of the study in a checkout of pkdb_data
(`--pkdb-data`) and the sha256 of the committed tables; the source digest and
the commit of pkdb_data of the collection are not known. Only studies whose
tables are tracked and have no uncommitted change are recorded; a study with
a record keeps it unless `--force`.

    uv run python scripts/backfill_study_records.py ../pkdb_models \\
        --pkdb-data ../pkdb_data/studies
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path, PurePosixPath

from pkpdsim.data import (
    StudyRecord,
    read_study_record,
    study_record_path,
    table_digests,
)

#: the tables of the studies in the history, by the pathspec of git
TABLES = ":(glob)**/.*.tsv"
#: the code of the models in the history, by the pathspec of git
CODE = ":(glob)**/*.py"
#: the record of an adopted fit, `.../fits/<fit id>/parameters.json`
FIT_RECORD = re.compile(
    r"/fits/(?P<fit>[^/]+_(?P<day>\d{8})_\d{6}__[0-9a-f]{5})/parameters\.json$"
)
#: a word which can be the name of a study
WORD = re.compile(r"[A-Za-z][A-Za-z0-9_-]*")


@dataclass
class Dating:
    """The first and the latest commit of the tables of a study folder."""

    first: date
    commit: str
    latest: date

    def earlier(self, day: date, commit: str) -> None:
        """Keep the earlier of the first day and another day."""
        if day < self.first:
            self.first, self.commit = day, commit


@dataclass
class Backfill:
    """What the backfill did, by the folder of the study relative to the repository."""

    written: dict[str, StudyRecord] = field(default_factory=dict)
    recorded: list[str] = field(default_factory=list)
    changed: list[str] = field(default_factory=list)
    untracked: list[str] = field(default_factory=list)


def _git(repo: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-c", "core.quotepath=off", "-C", str(repo), *args],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def _study_folder(path: str) -> str | None:
    """Get the folder of a table of a study, `None` for another file."""
    pure = PurePosixPath(path)
    if pure.name.startswith(f".{pure.parent.name}_") and pure.suffix == ".tsv":
        return pure.parent.as_posix()
    return None


def datings(repo: Path) -> dict[str, Dating]:
    """Date every study folder which the history of a repository has tables of.

    Args:
        repo: the git repository of the models.

    Returns:
        The dating of every study folder, relative to the repository, also of
        folders which no longer exist.
    """
    log = _git(
        repo,
        "log",
        "--reverse",
        "--date-order",
        "--format=%x00%H %as",
        "--name-status",
        "-M",
        "--",
        TABLES,
    )
    result: dict[str, Dating] = {}

    def touch(folder: str | None, day: date, commit: str, first: Dating | None) -> None:
        if folder is None:
            return
        dating = result.get(folder)
        if dating is None:
            dating = result[folder] = Dating(first=day, commit=commit, latest=day)
        dating.earlier(day, commit)
        if first is not None:
            dating.earlier(first.first, first.commit)
        dating.latest = max(dating.latest, day)

    for block in log.split("\0")[1:]:
        header, _, changes = block.partition("\n")
        commit, text = header.split()
        day = date.fromisoformat(text)
        for line in changes.splitlines():
            if not line.strip():
                continue
            status, *paths = line.split("\t")
            if status[0] in "RC":
                old, new = paths
                source = result.get(_study_folder(old) or "")
                # a renamed table carries the days of its folder, a copy is a
                # new collection of its own
                touch(
                    _study_folder(new),
                    day,
                    commit,
                    source if status[0] == "R" else None,
                )
                touch(_study_folder(old), day, commit, None)
            else:
                touch(_study_folder(paths[0]), day, commit, None)
    return result


def mentions(
    repo: Path, studies: dict[str, set[str]]
) -> dict[tuple[str, str], tuple[date, str]]:
    """Get the first commit which mentioned a study in the code of a model.

    Args:
        repo: the git repository of the models.
        studies: the names of the studies of every model, by the name of the
            folder of the model.

    Returns:
        The day and the commit of the first added line of a python file with
        the name of the study, in any case, below a folder of the name of the
        model, by model and study.
    """
    lowered = {
        model: {name.lower(): name for name in names}
        for model, names in studies.items()
    }
    log = _git(
        repo,
        "log",
        "--reverse",
        "--date-order",
        "--format=%x00%H %as",
        "-p",
        "--unified=0",
        "--no-renames",
        "--",
        CODE,
    )
    found: dict[tuple[str, str], tuple[date, str]] = {}
    for block in log.split("\0")[1:]:
        header, _, patch = block.partition("\n")
        commit, text = header.split()
        day = date.fromisoformat(text)
        models: list[str] = []
        for line in patch.splitlines():
            if line.startswith("+++ "):
                parts = line[6:].split("/") if line.startswith("+++ b/") else []
                models = [part for part in parts[:-1] if part in lowered]
                continue
            if not models or not line.startswith("+"):
                continue
            words = {word.lower() for word in WORD.findall(line)}
            for model in models:
                for word in words & lowered[model].keys():
                    key = (model, lowered[model][word])
                    # the order of the log is by commit date, the day the author's
                    if key not in found or day < found[key][0]:
                        found[key] = (day, commit)
    return found


def fits(repo: Path) -> dict[str, tuple[date, str]]:
    """Get the day and the commit of every adopted fit of the repository.

    Returns:
        The day of the fit (of its id) and the first commit which added its
        record at its path, by the path of the record.
    """
    result: dict[str, tuple[date, str]] = {}
    for path in _git(
        repo, "ls-files", "--", ":(glob)**/fits/*/parameters.json"
    ).splitlines():
        match = FIT_RECORD.search(path)
        if match is None:
            continue
        added = _git(repo, "log", "--diff-filter=A", "--format=%H", "--", path).split()
        if not added:
            continue
        text = match["day"]
        result[path] = (date(int(text[:4]), int(text[4:6]), int(text[6:])), added[-1])
    return result


def fits_of(
    models: set[str], adopted: dict[str, tuple[date, str]]
) -> dict[str, list[tuple[date, str]]]:
    """Get the adopted fits of every model folder.

    Args:
        models: the folders of the models, relative to the repository.
        adopted: the fits by the path of their record, see `fits`.

    Returns:
        The fits below every model folder, a fit of a nested model folder only
        of the nested model.
    """
    result: dict[str, list[tuple[date, str]]] = {model: [] for model in models}
    for path, fit in adopted.items():
        owners = [model for model in models if path.startswith(f"{model}/")]
        if owners:
            result[max(owners, key=len)].append(fit)
    return result


def _is_ancestor(repo: Path, commit: str, of: str) -> bool:
    """Get whether a commit is an ancestor of another one or the same."""
    result = subprocess.run(
        ["git", "-C", str(repo), "merge-base", "--is-ancestor", commit, of],
        capture_output=True,
        check=False,
    )
    return result.returncode == 0


def _model_folder(folder: str) -> str:
    """Get the folder of the model of a study folder, the parent of its `data`."""
    parts = PurePosixPath(folder).parts
    index = max((i for i, part in enumerate(parts) if part == "data"), default=0)
    return PurePosixPath(*parts[:index]).as_posix() if index else ""


def _day(
    repo: Path,
    dating: Dating,
    mentioned: tuple[date, str] | None,
    adopted: list[tuple[date, str]],
) -> tuple[date, str]:
    """Get the day of the first collection of a study folder and its commit.

    Args:
        repo: the git repository of the models.
        dating: the commits of the tables of the folder.
        mentioned: the first mention of the study in the code of its model.
        adopted: the day and the commit of every adopted fit of the model.

    Returns:
        The earliest trace of the study, no later than a fit committed with a
        trace present.
    """
    traces = [(dating.first, dating.commit)]
    if mentioned is not None:
        traces.append(mentioned)
    day, commit = min(traces)
    for fit_day, fit_commit in sorted(adopted):
        if fit_day < day and any(_is_ancestor(repo, c, fit_commit) for _, c in traces):
            return fit_day, fit_commit
    return day, commit


def _releases(pkdb_data: Path | None) -> dict[str, tuple[date | None, str | None]]:
    """Get the release of every study of a checkout of pkdb_data, by its name."""
    releases: dict[str, tuple[date | None, str | None]] = {}
    if pkdb_data is None:
        return releases
    for path in sorted(pkdb_data.glob("*/*/study.json")):
        release = json.loads(path.read_text(encoding="utf-8")).get("release") or {}
        day, pkdb_id = release.get("date"), release.get("pkdb_id")
        releases.setdefault(
            path.parent.name,
            (
                None if day is None else date.fromisoformat(str(day)),
                None if pkdb_id is None else str(pkdb_id),
            ),
        )
    return releases


def backfill(
    repo: Path,
    *,
    pkdb_data: Path | None = None,
    force: bool = False,
    dry_run: bool = False,
) -> Backfill:
    """Write the records of the tracked study folders of a repository from its history.

    Args:
        repo: the git repository of the models.
        pkdb_data: the folder `studies` of a checkout of pkdb_data, for the
            releases of the studies.
        force: replace an existing record.
        dry_run: write nothing.

    Returns:
        What was written and what was left out, by study folder.
    """
    repo = repo.resolve()
    tracked = sorted(
        {
            folder
            for path in _git(repo, "ls-files", "--", TABLES).splitlines()
            if (folder := _study_folder(path)) is not None
        }
    )
    changed = {
        folder
        for line in _git(repo, "status", "--porcelain", "--", TABLES).splitlines()
        if (folder := _study_folder(line[3:].split(" -> ")[-1])) is not None
    }
    found = datings(repo)
    models: dict[str, set[str]] = {}
    for folder in tracked:
        model = PurePosixPath(_model_folder(folder)).name
        models.setdefault(model, set()).add(PurePosixPath(folder).name)
    mentioned = mentions(repo, models)
    adopted = fits_of({_model_folder(f) for f in tracked}, fits(repo))
    releases = _releases(pkdb_data)
    result = Backfill()
    for folder in tracked:
        path = repo / folder
        if folder in changed:
            result.changed.append(folder)
            continue
        if not force and read_study_record(path) is not None:
            result.recorded.append(folder)
            continue
        dating = found[folder]
        model = _model_folder(folder)
        day, commit = _day(
            repo,
            dating,
            mentioned.get((PurePosixPath(model).name, path.name)),
            adopted.get(model, []),
        )
        released, pkdb_id = releases.get(path.name, (None, None))
        record = StudyRecord(
            study=path.name,
            collected=day,
            dating="git",
            dating_commit=commit,
            updated=dating.latest,
            released=released,
            pkdb_id=pkdb_id,
            tables=table_digests(path),
        )
        if not dry_run:
            study_record_path(path).write_text(
                json.dumps(record.to_dict(), indent=2) + "\n",
                encoding="utf-8",
                newline="\n",
            )
        result.written[folder] = record
    result.untracked = sorted(
        {
            folder
            for line in _git(
                repo, "ls-files", "--others", "--exclude-standard", "--", TABLES
            ).splitlines()
            if (folder := _study_folder(line)) is not None
        }
        - set(tracked)
    )
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("repo", type=Path, help="the git repository of the models")
    parser.add_argument(
        "--pkdb-data", type=Path, help="the folder studies of a checkout of pkdb_data"
    )
    parser.add_argument("--force", action="store_true", help="replace existing records")
    parser.add_argument("--dry-run", action="store_true", help="write nothing")
    args = parser.parse_args(argv)
    result = backfill(
        args.repo, pkdb_data=args.pkdb_data, force=args.force, dry_run=args.dry_run
    )
    verb = "would write" if args.dry_run else "wrote"
    print(f"{verb} {len(result.written)} records")
    released = sum(r.released is not None for r in result.written.values())
    print(f"  {released} of them with a release in PK-DB")
    for label, folders in (
        ("kept the record of", result.recorded),
        ("left out, tables with uncommitted changes:", result.changed),
        ("left out, tables not tracked:", result.untracked),
    ):
        if folders:
            print(f"{label} {len(folders)} studies")
            for folder in folders:
                print(f"  {folder}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
