"""Tests of the documentation site and the README of a release repository."""

import dataclasses
import json
import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

from pkpdsim.release.build import BuildPaths
from pkpdsim.release.docs import PAGES, write_docs
from pkpdsim.release.notes import fits_prefix
from pkpdsim.release.project import ProjectContext, write_project_files
from pkpdsim.uncertainty import FitRecord, adopt_fit
from pkpdsim_toy import MODEL, TOY_PATH

assert MODEL.release is not None
SETTINGS = MODEL.release


def write(path: Path, text: str | bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(text, bytes):
        path.write_bytes(text)
    else:
        path.write_text(text)


@pytest.fixture
def release(toy_fit: Path, tmp_path: Path) -> tuple[ProjectContext, Path]:
    """A built release repository of the toy model with an adopted fit."""
    model = dataclasses.replace(MODEL, fits=[toy_fit.name])
    path = tmp_path / "toy-model"
    adopt_fit(toy_fit, path / fits_prefix(model))
    (path / BuildPaths.MODELS).mkdir(parents=True)
    shutil.copy(TOY_PATH / "models" / "toy.xml", path / BuildPaths.MODELS / "toy.xml")
    write(path / BuildPaths.MODELS / "toy.md", "# toy")
    write(path / BuildPaths.FIGURES / "toy_model.png", b"png")
    write(path / BuildPaths.FIGURES / "toy.png", b"png")
    write(path / "toy_model.omex", b"omex")
    # the simulation report of sbmlsim: `index.html` and `<experiment>/<experiment>.html`
    write(path / BuildPaths.SIMULATIONS / "index.html", "<html/>")
    for experiment in ("Decay", "Slow"):
        write(
            path / BuildPaths.SIMULATIONS / experiment / f"{experiment}.html", "<html/>"
        )
    write(path / BuildPaths.SIMULATIONS / "_figures" / "Decay_fig.svg", "<svg/>")
    for version in ("0.9.0", "0.10.0", "0.2.0"):
        write(path / "release-notes" / f"{version}.md", f"# Notes {version}\n\n## A\n")
    context = ProjectContext(
        model=model,
        settings=SETTINGS,
        version="0.10.0",
        python="3.13",
        dependencies=(),
        constraints=(),
        pkpdsim_path=None,
        year=2026,
        requires_python=">=3.13",
    )
    return context, path


def test_docs_pages(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    write_project_files(context, path)
    written = write_docs(context, path)
    assert {p.name for p in (path / "docs").glob("*.md")} == {f"{p}.md" for p in PAGES}
    assert path / "README.md" in written
    nav = tomllib.loads((path / "zensical.toml").read_text())["project"]["nav"]
    in_nav = {file for entry in nav for file in entry.values()}
    assert in_nav == {f"{page}.md" for page in PAGES}
    assert (path / "docs" / "images" / "toy_model.png").is_file()
    assert path / "docs" / "images" / "toy_model.png" in written
    index = (path / "docs" / "index.md").read_text()
    assert "images/toy_model.png" in index
    assert (
        "Doe, J. (2026). *Toy model of the decay of a substance*. Zenodo. "
        "https://doi.org/10.5281/zenodo.0000000" in index
    )
    assert index.startswith("# Toy model of the decay of a substance\n")
    readme = (path / "README.md").read_text()
    assert "# Toy model of the decay of a substance\n" in readme
    assert "*Toy model of the decay of a substance*" in readme
    assert "Toy PK model" not in index + readme


def test_docs_replaced(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    write(path / "docs" / "stale.md", "stale")
    write_docs(context, path)
    assert not (path / "docs" / "stale.md").exists()


def test_models_links(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    write_docs(context, path)
    models = (path / "docs" / "models.md").read_text()
    raw = SETTINGS.raw_url
    assert f"https://sbml4humans.de/model_url?url={raw}models/toy.xml" in models
    assert f"({raw}models/toy.xml)" in models
    assert "https://github.com/pkpdsim/toy-model/blob/main/models/toy.md" in models
    assert (
        '![toy model diagram](images/toy.png){ style="max-height: 420px; '
        'width: auto; max-width: 100%" }' in models
    )
    assert f"({raw}toy_model.omex)" in models


def test_simulations_links(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    write_docs(context, path)
    page = (path / "docs" / "simulations.md").read_text()
    site = SETTINGS.site_url
    assert f"{site}simulations/index.html" in page
    assert f"{site}simulations_uncertainty/index.html" in page
    assert f"{site}simulations/Decay/Decay.html" in page
    assert f"{site}simulations/Slow/Slow.html" in page
    assert "_figures" not in page
    assert "no experiment reports" not in page


def test_fit_table(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    write_docs(context, path)
    fit_id = context.model.fits[0]
    record = FitRecord.load(path / fits_prefix(context.model) / fit_id)
    page = (path / "docs" / "fit.md").read_text()
    assert f"{SETTINGS.site_url}fit/{fit_id}/index.html" in page
    for pid, parameter in record.parameters.items():
        (row,) = [line for line in page.splitlines() if line.startswith(f"| `{pid}`")]
        if parameter.ci_lower is None or parameter.ci_upper is None:
            assert "open" in row
        else:
            assert f"{parameter.ci_lower:.4g}" in row


def test_fit_table_open_interval(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    fit_id = context.model.fits[0]
    directory = path / fits_prefix(context.model) / fit_id
    record = FitRecord.load(directory)
    pid = next(iter(record.parameters))
    profiles = directory / "identifiability.json"
    text = profiles.read_text()
    assert pid in text
    # open the interval of the parameter at both ends
    data = json.loads(text)
    entry = data["profiles"][pid]
    for key in ("ci_lower", "ci_upper"):
        entry[key] = None
    profiles.write_text(json.dumps(data))
    reloaded = FitRecord.load(directory).parameters[pid]
    if reloaded.ci_lower is not None:
        pytest.skip("the record format keeps the interval of the profile")
    write_docs(context, path)
    page = (path / "docs" / "fit.md").read_text()
    (row,) = [x for x in page.splitlines() if x.startswith(f"| `{pid}`")]
    assert "open - open" in row


def test_release_notes_sorted(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    write_docs(context, path)
    page = (path / "docs" / "release-notes.md").read_text()
    order = [page.index(f"Notes {v}") for v in ("0.10.0", "0.9.0", "0.2.0")]
    assert order == sorted(order)
    assert "\n## Notes 0.10.0" in page


def test_readme_thin(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    write_docs(context, path)
    readme = (path / "README.md").read_text()
    assert SETTINGS.site_url in readme
    assert len(readme.splitlines()) < 60
    assert "uv run run_simulation pkpdsim_toy --action all -r results" in readme
    assert "10.5281/zenodo.0000000" in readme


@pytest.mark.slow
def test_zensical_strict(release: tuple[ProjectContext, Path]) -> None:
    context, path = release
    write_project_files(context, path)
    write_docs(context, path)
    result = subprocess.run(
        ["zensical", "build", "--clean", "--strict"],
        cwd=path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    site = path / "site"
    assert (site / "index.html").is_file()

    # the docs workflow copies the reports into the site, they must not
    # replace the pages of the docs
    write(path / BuildPaths.FITS / "fit1" / "index.html", "<html/>")
    for folder in (BuildPaths.SIMULATIONS, BuildPaths.FITS):
        shutil.copytree(path / folder, site / folder, dirs_exist_ok=True)
    # the figure of a model is bounded (attr_list), the overview is not
    models = (site / "models.html").read_text()
    assert 'alt="toy model diagram"' in models
    assert "max-height: 420px" in models
    assert "max-height" not in (site / "index.html").read_text()

    texts = [page.read_text() for page in site.rglob("*.html")]
    assert any("Simulation report with the confidence interval" in t for t in texts)
    assert any("Fit report" in t for t in texts)
