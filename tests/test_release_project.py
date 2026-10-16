"""Tests of the project files of a release repository."""

import dataclasses
import json
import tomllib
from pathlib import Path

import pytest
from jinja2 import UndefinedError

from pkpdsim.release import project as project_module
from pkpdsim.release.project import (
    ProjectContext,
    project_context,
    render,
    write_project_files,
    zenodo_json,
)
from pkpdsim.release.settings import Author
from pkpdsim.release.source import SourceRepo
from pkpdsim_toy import MODEL

assert MODEL.release is not None
SETTINGS = MODEL.release

PYPROJECT = """\
[project]
name = "source"
requires-python = ">=3.13"
dependencies = ["sbmlsim>=0.9.7", "sbmlutils>=0.16.0"]
"""

LOCK = """\
version = 1

[[package]]
name = "source"
version = "0.1.0"
source = { editable = "." }

[[package]]
name = "sbmlsim"
version = "0.9.7"
source = { registry = "https://pypi.org/simple" }

[[package]]
name = "pkpdsim"
version = "0.2.0"
source = { registry = "https://pypi.org/simple" }

[[package]]
name = "gitdep"
version = "1.0.0"
source = { git = "https://github.com/x/gitdep" }
"""


@pytest.fixture
def source(tmp_path: Path) -> SourceRepo:
    root = tmp_path / "source"
    root.mkdir()
    (root / "pyproject.toml").write_text(PYPROJECT)
    (root / "uv.lock").write_text(LOCK)
    (root / ".python-version").write_text("3.14\n")
    return SourceRepo(root=root, commit="abc")


@pytest.fixture
def context(source: SourceRepo) -> ProjectContext:
    return project_context(MODEL, SETTINGS, source, "0.1.0")


def test_constraints_from_lock(context: ProjectContext) -> None:
    assert context.constraints == ("sbmlsim==0.9.7", "pkpdsim==0.2.0")
    assert context.dependencies == ("sbmlsim>=0.9.7", "sbmlutils>=0.16.0")
    assert context.python == "3.14"


def test_python_from_requires_python(source: SourceRepo) -> None:
    (source.root / ".python-version").unlink()
    assert project_context(MODEL, SETTINGS, source, "0.1.0").python == "3.13"


def test_pyproject_parses(context: ProjectContext) -> None:
    data = tomllib.loads(render("pyproject.toml.j2", context))
    assert data["project"]["version"] == "0.1.0"
    assert data["project"]["name"] == "toy-model"
    assert data["project"]["requires-python"] == ">=3.13"
    assert "scripts" not in data["project"]
    assert "sources" not in data["tool"]["uv"]
    assert data["tool"]["uv"]["constraint-dependencies"] == [
        "sbmlsim==0.9.7",
        "pkpdsim==0.2.0",
    ]
    assert data["tool"]["hatch"]["build"]["targets"]["wheel"]["packages"] == [
        "src/pkpdsim_toy"
    ]
    assert data["dependency-groups"]["docs"] == ["zensical"]


def test_pyproject_local_pkpdsim(source: SourceRepo) -> None:
    context = project_context(MODEL, SETTINGS, source, "0.1.0", Path("/x/pkpdsim"))
    data = tomllib.loads(render("pyproject.toml.j2", context))
    assert data["tool"]["uv"]["sources"]["pkpdsim"] == {
        "path": "/x/pkpdsim",
        "editable": True,
    }
    assert "pkpdsim==0.2.0" not in data["tool"]["uv"]["constraint-dependencies"]


def test_zenodo_json(context: ProjectContext) -> None:
    data = zenodo_json(context)
    assert data["creators"] == [
        {"orcid": "0000-0000-0000-0000", "affiliation": "Lab", "name": "Doe, Jane"}
    ]
    assert "No funding." in data["description"]
    assert "https://github.com/pkpdsim/toy-model/issues" in data["description"]
    assert data["license"] == "MIT"
    settings = dataclasses.replace(SETTINGS, authors=(Author("Roe, Rick"),))
    anonymous = dataclasses.replace(context, settings=settings)
    assert zenodo_json(anonymous)["creators"] == [{"name": "Roe, Rick"}]


def test_workflows_replaced(context: ProjectContext, tmp_path: Path) -> None:
    workflows = tmp_path / ".github" / "workflows"
    workflows.mkdir(parents=True)
    (workflows / "python.yml").write_text("old")
    written = write_project_files(context, tmp_path)
    assert not (workflows / "python.yml").exists()
    assert {p.name for p in workflows.iterdir()} == {"ci.yml", "docs.yml"}
    assert tmp_path / ".zenodo.json" in written
    assert (tmp_path / ".python-version").read_text() == "3.14\n"
    assert "toy:0.1.0" in (tmp_path / "Dockerfile").read_text()
    assert (
        "${{ steps.deployment.outputs.page_url }}"
        in (workflows / "docs.yml").read_text()
    )
    assert tomllib.loads((tmp_path / "zensical.toml").read_text())["project"]["nav"]


def test_no_em_dash(context: ProjectContext, tmp_path: Path) -> None:
    for path in write_project_files(context, tmp_path):
        assert "—" not in path.read_text(), path


def test_strict_undefined(
    context: ProjectContext, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        project_module,
        "_template",
        lambda name: project_module._environment().from_string("{{ unknown }}"),
    )
    with pytest.raises(UndefinedError):
        render("x.j2", context)


def test_zenodo_file_is_json(context: ProjectContext, tmp_path: Path) -> None:
    write_project_files(context, tmp_path)
    assert json.loads((tmp_path / ".zenodo.json").read_text())["title"] == MODEL.title
