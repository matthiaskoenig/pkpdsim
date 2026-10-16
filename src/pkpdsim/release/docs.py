"""The documentation site and the README of a release repository.

The pages of the Zensical site (`docs/`) and the thin `README.md` are rendered
from the templates in `resources/release/docs` and `resources/release` with the
content of the built release repository: the models, the experiment reports and
the records of the adopted fits. The reports themselves are copied into the site
by the docs workflow and are linked by the URL of the site.
"""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from pkpdsim.release.build import BuildPaths
from pkpdsim.release.notes import NOTES_DIR, fits_prefix
from pkpdsim.release.project import ProjectContext, render
from pkpdsim.release.source import ReleaseError
from pkpdsim.uncertainty import FitRecord

PAGES = ("index", "running", "models", "simulations", "fit", "release-notes")
"""The pages of the site, the entries of its navigation."""

_IMAGES = "images"
_HEADING = re.compile(r"^(#+ )", re.MULTILINE)
_VERSION = re.compile(r"^\d+(\.\d+)*$")


@dataclass(frozen=True)
class ModelEntry:
    """A model of the release repository as the models page shows it.

    Attributes:
        stem: the name of the model file without extension.
        file: the file name of the SBML model.
        equations: if the model has an equations file.
        figure: the name of the figure in the images of the site, if any.
    """

    stem: str
    file: str
    equations: bool
    figure: str | None


@dataclass(frozen=True)
class FitRow:
    """A fitted parameter as the fit page shows it."""

    pid: str
    value: str
    unit: str
    interval: str
    identifiability: str


@dataclass(frozen=True)
class FitEntry:
    """An adopted fit as the fit page shows it."""

    id: str
    rows: list[FitRow]


def _given_family(name: str) -> str:
    """Turn `"Family, Given"` into `"Given Family"`."""
    family, _, given = name.partition(",")
    return f"{given.strip()} {family.strip()}".strip()


def _join(names: list[str], last: str) -> str:
    if len(names) <= 1:
        return "".join(names)
    return ", ".join(names[:-1]) + last + names[-1]


def citation(context: ProjectContext) -> str:
    """Get the citation of the release.

    Args:
        context: the project context.

    Returns:
        The authors, the year, the title and the DOI of the Zenodo record.
    """
    settings = context.settings
    authors = _join([a.citation_name for a in settings.authors], " & ")
    return (
        f"{authors} ({context.year}). *{context.model.title}*. Zenodo. "
        f"https://doi.org/{settings.doi}"
    )


def _number(value: float | None) -> str:
    return "open" if value is None else f"{value:.4g}"


def _fit_entries(context: ProjectContext, release_path: Path) -> list[FitEntry]:
    """Read the records of the adopted fits from the exported package."""
    directory = release_path / fits_prefix(context.model)
    entries = []
    for fit_id in context.model.fits:
        try:
            record = FitRecord.load(directory / fit_id)
        except (FileNotFoundError, ValueError) as err:
            raise ReleaseError(
                f"The adopted fit '{fit_id}' has no readable record: {err}"
            ) from err
        rows = [
            FitRow(
                pid=pid,
                value=_number(p.value),
                unit=p.unit or "-",
                interval=f"{_number(p.ci_lower)} - {_number(p.ci_upper)}",
                identifiability=str(p.identifiability or "-"),
            )
            for pid, p in sorted(record.parameters.items())
        ]
        entries.append(FitEntry(fit_id, rows))
    return entries


def _copy_image(release_path: Path, name: str, images: Path) -> bool:
    """Copy a figure of the release into the images of the site, if it exists."""
    source = release_path / BuildPaths.FIGURES / name
    if not source.is_file():
        return False
    images.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, images / name)
    return True


def _models(release_path: Path, images: Path) -> list[ModelEntry]:
    entries = []
    for path in sorted((release_path / BuildPaths.MODELS).glob("*.xml")):
        figure = f"{path.stem}.png"
        entries.append(
            ModelEntry(
                stem=path.stem,
                file=path.name,
                equations=path.with_suffix(".md").is_file(),
                figure=figure if _copy_image(release_path, figure, images) else None,
            )
        )
    return entries


def _experiments(release_path: Path) -> list[str]:
    folder = release_path / BuildPaths.SIMULATIONS
    if not folder.is_dir():
        return []
    return sorted(
        path.name
        for path in folder.iterdir()
        if not path.name.startswith("_") and (path / "index.html").is_file()
    )


def _version_key(path: Path) -> tuple[int, ...]:
    return tuple(int(part) for part in path.stem.split("."))


def _release_notes(release_path: Path) -> list[str]:
    """Get the release notes, newest version first, with the headings demoted."""
    files = [
        path
        for path in (release_path / NOTES_DIR).glob("*.md")
        if _VERSION.match(path.stem)
    ]
    return [
        _HEADING.sub(r"#\1", path.read_text("utf-8").strip())
        for path in sorted(files, key=_version_key, reverse=True)
    ]


def write_docs(context: ProjectContext, release_path: Path) -> list[Path]:
    """Write the documentation site and the README into the release repository.

    The folder `docs` is replaced. The release repository must be built, the
    pages are made from its models, figures, reports, fit records and release
    notes.

    Args:
        context: the project context.
        release_path: the folder of the release repository.

    Returns:
        The written files.

    Raises:
        ReleaseError: if the record of an adopted fit cannot be read.
    """
    docs = release_path / "docs"
    if docs.exists():
        shutil.rmtree(docs)
    docs.mkdir()
    images = docs / _IMAGES

    figure = f"{context.model.name}_model.png"
    omex = next(iter(sorted(release_path.glob("*.omex"))), None)
    extra = {
        "figure": figure if _copy_image(release_path, figure, images) else None,
        "citation": citation(context),
        "authors": _join(
            [_given_family(a.name) for a in context.settings.authors], " and "
        ),
        "models": _models(release_path, images),
        "omex": omex.name if omex else None,
        "experiments": _experiments(release_path),
        "fits": _fit_entries(context, release_path),
        "fit_ids": list(context.model.fits),
        "notes": _release_notes(release_path),
    }

    written = []
    for page in PAGES:
        path = docs / f"{page}.md"
        path.write_text(render(f"docs/{page}.md.j2", context, extra), "utf-8")
        written.append(path)
    readme = release_path / "README.md"
    readme.write_text(render("README.md.j2", context, extra), "utf-8")
    written.append(readme)
    return written
