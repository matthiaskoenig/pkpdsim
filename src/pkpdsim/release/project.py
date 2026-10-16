"""The project files of a release repository, written from templates.

The files which are the same in structure for every model (`pyproject.toml`,
`Dockerfile`, the workflows, the zensical configuration, `.zenodo.json`) are
rendered from the templates in `resources/release` with the settings of the
model and the dependencies of the source project.
"""

from __future__ import annotations

import datetime
import html
import json
import re
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Any

from jinja2 import Environment, StrictUndefined, Template

from pkpdsim.model import PKPDModel
from pkpdsim.release.settings import ReleaseSettings
from pkpdsim.release.source import ReleaseError, SourceRepo

_PYTHON_FLOOR = re.compile(r">=\s*(\d+\.\d+)")
_PYTHON_VERSION = re.compile(r"^\d+\.\d+")

_TEMPLATES = {
    "pyproject.toml": "pyproject.toml.j2",
    "Dockerfile": "Dockerfile.j2",
    ".gitignore": "gitignore.j2",
    ".python-version": "python-version.j2",
    "zensical.toml": "zensical.toml.j2",
    ".github/workflows/ci.yml": "workflows/ci.yml.j2",
    ".github/workflows/docs.yml": "workflows/docs.yml.j2",
}
"""The files of the release repository by their template."""


@dataclass(frozen=True)
class ProjectContext:
    """The values the templates of the project files are rendered with.

    Attributes:
        model: the model which is released.
        settings: the release settings of the model.
        version: the version of the release.
        python: the Python version of the source, e.g. `"3.14"`.
        dependencies: the dependencies of the source project.
        constraints: `name==version` for every registry package of the source lock.
        pkpdsim_path: local checkout of pkpdsim to use instead of the registry.
        year: the year of the release.
        requires_python: the `requires-python` of the source project.
    """

    model: PKPDModel
    settings: ReleaseSettings
    version: str
    python: str
    dependencies: tuple[str, ...]
    constraints: tuple[str, ...]
    pkpdsim_path: Path | None
    year: int
    requires_python: str

    @property
    def repo_name(self) -> str:
        """The name of the release repository."""
        return self.settings.repo_name

    @property
    def top_package(self) -> str:
        """The first component of the model package, the folder below `src`."""
        return self.model.package.split(".")[0]


def _python_version(source: SourceRepo, requires_python: str) -> str:
    version_file = source.root / ".python-version"
    if version_file.exists():
        match = _PYTHON_VERSION.match(version_file.read_text().strip())
        if match:
            return match.group(0)
    floor = _PYTHON_FLOOR.search(requires_python)
    if floor is None:
        raise ReleaseError(
            "The Python version of the source is unknown, add a '.python-version' "
            f"or a '>=X.Y' requires-python (is '{requires_python}')."
        )
    return floor.group(1)


def _constraints(source: SourceRepo, pkpdsim_path: Path | None) -> tuple[str, ...]:
    """Get `name==version` of the registry packages of the source lock.

    Packages locked at several versions (resolution forks) are skipped, as
    pins to two versions of one package cannot be satisfied.
    """
    lock_path = source.root / "uv.lock"
    if not lock_path.exists():
        raise ReleaseError(f"The source has no 'uv.lock' at {lock_path}.")
    lock = tomllib.loads(lock_path.read_text())
    versions: dict[str, set[str]] = {}
    for package in lock.get("package", []):
        if "registry" not in package.get("source", {}):
            continue
        name = package["name"]
        if pkpdsim_path is not None and name == "pkpdsim":
            continue
        versions.setdefault(name, set()).add(package["version"])
    # a package locked at several versions (resolution forks) gets no constraint
    constraints = [
        f"{name}=={next(iter(found))}"
        for name, found in versions.items()
        if len(found) == 1
    ]
    return tuple(constraints)


def project_context(
    model: PKPDModel,
    settings: ReleaseSettings,
    source: SourceRepo,
    version: str,
    pkpdsim_path: Path | None = None,
) -> ProjectContext:
    """Read the dependencies of the source project into a project context.

    Args:
        model: the model which is released.
        settings: the release settings of the model.
        source: the source repository with `pyproject.toml` and `uv.lock`.
        version: the version of the release.
        pkpdsim_path: local checkout of pkpdsim, which then replaces the
            registry pkpdsim in the release project.

    Returns:
        The context of the templates.

    Raises:
        ReleaseError: if the source lacks `uv.lock` or its Python version is unknown.
    """
    pyproject = tomllib.loads((source.root / "pyproject.toml").read_text())
    project = pyproject["project"]
    requires_python = project["requires-python"]
    return ProjectContext(
        model=model,
        settings=settings,
        version=version,
        python=_python_version(source, requires_python),
        dependencies=tuple(project.get("dependencies", ())),
        constraints=_constraints(source, pkpdsim_path),
        pkpdsim_path=pkpdsim_path,
        year=datetime.date.today().year,
        requires_python=requires_python,
    )


def _toml_string(value: object) -> str:
    """Quote a value as TOML basic string (JSON strings are valid TOML)."""
    return json.dumps(str(value), ensure_ascii=False)


def _environment() -> Environment:
    environment = Environment(
        undefined=StrictUndefined,
        keep_trailing_newline=True,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    environment.filters["toml"] = _toml_string
    return environment


def _template(name: str) -> Template:
    source = files("pkpdsim").joinpath("resources/release", name).read_text("utf-8")
    return _environment().from_string(source)


def render(
    name: str, context: ProjectContext, extra: Mapping[str, Any] | None = None
) -> str:
    """Render a template of `resources/release`.

    Args:
        name: the path of the template below `resources/release`.
        context: the values of the template, available by field name.
        extra: further values of the template, e.g. the changes of the release
            notes.

    Returns:
        The rendered text.

    Raises:
        jinja2.UndefinedError: if the template uses a value the context lacks.
    """
    values: dict[str, Any] = {
        "model": context.model,
        "settings": context.settings,
        "version": context.version,
        "python": context.python,
        "dependencies": context.dependencies,
        "constraints": context.constraints,
        "pkpdsim_path": str(context.pkpdsim_path) if context.pkpdsim_path else None,
        "year": context.year,
        "requires_python": context.requires_python,
        "repo_name": context.repo_name,
        "top_package": context.top_package,
        **(extra or {}),
    }
    return _template(name).render(**values)


def zenodo_json(context: ProjectContext) -> dict[str, Any]:
    """Get the metadata of the Zenodo record of the release.

    Args:
        context: the project context.

    Returns:
        The content of `.zenodo.json`.
    """
    settings = context.settings
    creators = []
    for author in settings.authors:
        creator: dict[str, str] = {}
        if author.orcid is not None:
            creator["orcid"] = author.orcid
        if author.affiliation is not None:
            creator["affiliation"] = author.affiliation
        creator["name"] = author.name
        creators.append(creator)
    issues = f"https://github.com/{settings.repository}/issues"
    description = (
        f"<p>{html.escape(settings.description, quote=False)}</p>\n"
        f'<p>If you have any questions or problems please <a href="{issues}">'
        "open an issue</a></p>\n"
        f"<h2>Funding</h2><p>{html.escape(settings.funding, quote=False)}</p>"
    )
    return {
        "upload_type": "software",
        "title": context.model.title,
        "creators": creators,
        "description": description,
        "access_right": "open",
        "license": "MIT",
        "keywords": list(settings.keywords),
        "contributors": [],
        "communities": [],
        "grants": [],
    }


def write_project_files(context: ProjectContext, release_path: Path) -> list[Path]:
    """Write the project files into the release repository.

    Workflows in `.github/workflows` which are not written are removed.

    Args:
        context: the project context.
        release_path: the folder of the release repository.

    Returns:
        The written files.
    """
    written: list[Path] = []
    for target, template in _TEMPLATES.items():
        path = release_path / target
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render(template, context), encoding="utf-8")
        written.append(path)
    zenodo = release_path / ".zenodo.json"
    zenodo.write_text(
        json.dumps(zenodo_json(context), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    written.append(zenodo)
    for workflow in (release_path / ".github" / "workflows").glob("*"):
        if workflow.is_file() and workflow not in written:
            workflow.unlink()
    return written
