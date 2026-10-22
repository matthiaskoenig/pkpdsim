"""Tests of the project files of a release repository."""

import dataclasses
import json
import os
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path, PureWindowsPath

import pytest
import yaml
from jinja2 import UndefinedError

from pkpdsim.release import project as project_module
from pkpdsim.release.project import (
    UV_VERSION,
    ProjectContext,
    project_context,
    render,
    write_project_files,
    zenodo_json,
)
from pkpdsim.release.settings import Author
from pkpdsim.release.source import ReleaseError, SourceRepo
from pkpdsim_toy import MODEL

assert MODEL.release is not None
SETTINGS = MODEL.release

PYPROJECT = """\
[project]
name = "source"
description = "Model of Łódź"
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
    # TOML files are UTF-8, whatever the encoding of the locale
    (root / "pyproject.toml").write_text(PYPROJECT, encoding="utf-8")
    (root / "uv.lock").write_text(LOCK, encoding="utf-8")
    (root / ".python-version").write_text("3.14\n", encoding="utf-8")
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
    assert "allow-direct-references" not in render("pyproject.toml.j2", context)
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


def test_pyproject_local_pkpdsim_on_windows(context: ProjectContext) -> None:
    """A Windows path is a TOML string with forward slashes, which uv accepts."""
    windows = dataclasses.replace(
        context, pkpdsim_path=PureWindowsPath(r"C:\Users\me\pkpdsim")
    )
    data = tomllib.loads(render("pyproject.toml.j2", windows))
    assert data["tool"]["uv"]["sources"]["pkpdsim"]["path"] == "C:/Users/me/pkpdsim"


def test_pyproject_one_entry_per_line(source: SourceRepo) -> None:
    context = project_context(MODEL, SETTINGS, source, "0.1.0", Path("/x/pkpdsim"))
    text = render("pyproject.toml.j2", context)
    assert (
        '\ndependencies = [\n    "sbmlsim>=0.9.7",\n    "sbmlutils>=0.16.0",\n]\n'
        in text
    )
    assert '\nconstraint-dependencies = [\n    "sbmlsim==0.9.7",\n]\n' in text
    assert '\nauthors = [\n    { name = "Doe, Jane" },\n]\n' in text
    assert "]\n\n[tool.uv.sources]\n" in text


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
    dockerfile = (tmp_path / "Dockerfile").read_text()
    assert dockerfile.startswith(
        "# The image is built and pushed by pkpd_release.\n\nFROM python:3.14-slim\n"
    )
    assert "docker build" not in dockerfile
    assert (
        "${{ steps.deployment.outputs.page_url }}"
        in (workflows / "docs.yml").read_text()
    )
    zensical = tomllib.loads((tmp_path / "zensical.toml").read_text())["project"]
    assert zensical["nav"]
    assert zensical["site_author"] == "Jane Doe"
    assert zensical["copyright"].endswith(" Jane Doe")


def test_no_em_dash(context: ProjectContext, tmp_path: Path) -> None:
    for path in write_project_files(context, tmp_path):
        assert "\u2014" not in path.read_text(encoding="utf-8"), path


def test_project_files_have_unix_line_endings(
    context: ProjectContext, tmp_path: Path
) -> None:
    """The files of the release are the same on every platform."""
    for path in write_project_files(context, tmp_path):
        assert b"\r" not in path.read_bytes(), path


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
    zenodo = json.loads((tmp_path / ".zenodo.json").read_text())
    assert zenodo["title"] == "Toy model of the decay of a substance"
    site_name = tomllib.loads((tmp_path / "zensical.toml").read_text())["project"]
    assert site_name["site_name"] == "Toy model of the decay of a substance"


def test_forked_package_not_constrained(source: SourceRepo) -> None:
    fork = (
        '\n[[package]]\nname = "numpy"\nversion = "1.0.0"\n'
        'source = { registry = "https://pypi.org/simple" }\n'
        '\n[[package]]\nname = "numpy"\nversion = "2.0.0"\n'
        'source = { registry = "https://pypi.org/simple" }\n'
    )
    (source.root / "uv.lock").write_text(LOCK + fork)
    constraints = project_context(MODEL, SETTINGS, source, "0.1.0").constraints
    assert constraints == ("sbmlsim==0.9.7", "pkpdsim==0.2.0")


def test_workflows_are_yaml(context: ProjectContext) -> None:
    ci = yaml.safe_load(render("workflows/ci.yml.j2", context))
    docs = yaml.safe_load(render("workflows/docs.yml.j2", context))
    runs = [s.get("run", "") for s in ci["jobs"]["simulate"]["steps"]]
    assert f"uv run run_simulation {MODEL.package} --action factory" in runs
    assert "deploy" in docs["jobs"]
    # a tag push does not run the tests again (YAML 1.1 reads the key `on` as True)
    assert ci[True]["push"] == {"branches": ["main"]}


def test_uv_version(context: ProjectContext) -> None:
    """The Docker image and the workflows use one uv, at least the fix of the
    CA store (GHSA-hq85-5cpj-mr62)."""
    assert UV_VERSION == "0.13.4"
    dockerfile = render("Dockerfile.j2", context)
    assert "COPY --from=ghcr.io/astral-sh/uv:0.13.4 /uv /bin/uv\n" in dockerfile
    for name in ("ci.yml.j2", "docs.yml.j2"):
        workflow = yaml.safe_load(render(f"workflows/{name}", context))
        steps = [
            step
            for job in workflow["jobs"].values()
            for step in job["steps"]
            if step.get("uses", "").startswith("astral-sh/setup-uv@")
        ]
        assert steps, name
        assert {step["with"]["version"] for step in steps} == {"0.13.4"}, name


@pytest.mark.skipif(sys.platform == "win32", reason="the workflow runs on linux")
@pytest.mark.skipif(shutil.which("uv") is None, reason="the step runs uv")
def test_ci_workflow_puts_libpython_on_the_library_path(
    context: ProjectContext, tmp_path: Path
) -> None:
    """The extension of libroadrunner links against libpython, which the
    interpreters of uv keep in their `lib`: the step puts it on the library
    search path before the first command which imports roadrunner."""
    ci = yaml.safe_load(render("workflows/ci.yml.j2", context))
    runs = [s.get("run", "") for s in ci["jobs"]["simulate"]["steps"]]
    (index,) = [k for k, run in enumerate(runs) if "LD_LIBRARY_PATH" in run]
    first = next(k for k, run in enumerate(runs) if run.startswith("uv run "))
    assert index < first
    env = tmp_path / "github_env"
    result = subprocess.run(
        ["bash", "-e", "-c", runs[index]],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        env={**os.environ, "GITHUB_ENV": str(env)},
    )
    assert result.returncode == 0, result.stderr
    (line,) = env.read_text().splitlines()
    assert line.startswith("LD_LIBRARY_PATH=/")
    assert line.split("=", 1)[1].split(":")[0].endswith("/lib")


@pytest.mark.skipif(sys.platform == "win32", reason="the workflow runs on linux")
@pytest.mark.parametrize(
    "folders",
    [
        ("simulations", "simulations_uncertainty", "fit"),
        ("simulations", "simulations_uncertainty"),
    ],
)
def test_docs_workflow_copies_existing_reports(
    context: ProjectContext, tmp_path: Path, folders: tuple[str, ...]
) -> None:
    """A model without adopted fits has no `fit` folder in the checkout."""
    docs = yaml.safe_load(render("workflows/docs.yml.j2", context))
    (step,) = [s for s in docs["jobs"]["build"]["steps"] if "site/" in s.get("run", "")]
    for folder in folders:
        (tmp_path / folder).mkdir()
        (tmp_path / folder / "index.html").write_text("<html/>")
    (tmp_path / "site").mkdir()
    result = subprocess.run(
        ["bash", "-e", "-c", step["run"]], cwd=tmp_path, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    copied = {p.parent.name for p in (tmp_path / "site").rglob("index.html")}
    assert copied == set(folders)


def test_gitignore(context: ProjectContext, tmp_path: Path) -> None:
    write_project_files(context, tmp_path)
    ignored = (tmp_path / ".gitignore").read_text().splitlines()
    assert {"/results/", "/build/", "/dist/", "/site/", ".venv/", ".idea/"} <= set(
        ignored
    )
    assert not {"results/", "build/", "dist/", "site/"} & set(ignored)


def test_readme_for_the_build(context: ProjectContext, tmp_path: Path) -> None:
    """The build installs the project, whose `readme` must exist before the docs."""
    written = write_project_files(context, tmp_path)
    readme = tmp_path / "README.md"
    assert readme in written
    assert readme.read_text() == f"# {SETTINGS.title}\n"
    readme.write_text("the README of the last release\n")
    assert readme not in write_project_files(context, tmp_path)
    assert readme.read_text() == "the README of the last release\n"


@pytest.mark.parametrize(
    ("first_year", "years"), [(2026, "2026"), (2024, "2024-2026"), (2027, "2027")]
)
def test_copyright_years(
    context: ProjectContext, tmp_path: Path, first_year: int, years: str
) -> None:
    settings = dataclasses.replace(SETTINGS, first_year=first_year)
    context = dataclasses.replace(context, settings=settings, year=2026)
    assert context.copyright_years == years
    write_project_files(context, tmp_path)
    zensical = tomllib.loads((tmp_path / "zensical.toml").read_text())["project"]
    assert zensical["copyright"] == f"Copyright &copy; {years} Jane Doe"


def test_source_without_requires_python(source: SourceRepo) -> None:
    pyproject = source.root / "pyproject.toml"
    pyproject.write_text('[project]\nname = "source"\n', encoding="utf-8")
    with pytest.raises(ReleaseError, match="requires-python"):
        project_context(MODEL, SETTINGS, source, "0.1.0")
    pyproject.unlink()
    with pytest.raises(ReleaseError, match=r"no 'pyproject\.toml'"):
        project_context(MODEL, SETTINGS, source, "0.1.0")
