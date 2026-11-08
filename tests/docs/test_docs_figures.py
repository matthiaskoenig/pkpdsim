"""The generated figures of the docs are up to date.

The figures of the page `docs/pbpk/figures.md` and the figures A and B of every
model a page of `docs/pbpk/` or `docs/perfusion.md` creates, a PBPK model or a
machine perfusion circuit.
"""

import importlib.util
import re
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from scipy.spatial import KDTree

from pkpdsim.perfusion import PerfusionDefinition

ROOT = Path(__file__).parents[2]
NUMBER = re.compile(r"-?\d+\.\d+")
# the path data of an element, compared as geometry
GEOMETRY = re.compile(r'(?<= d=")[^"]*(?=")')
# the traced outlines of the bodies, closed polylines `M x,y x,y ...Z`
POLYLINES = re.compile(r"(M[-\d., ]+Z?)+")
# a cell of the raster of a body (0.25 cm) can flip with the last bit of a
# number on another platform (macOS), which moves its outline by about a cell
# and the numbers which depend on it
TOLERANCE = 0.5
# any version, so that a release does not outdate the stored figures
VERSION = re.compile(r"pkpdsim</a> \d+\.\d+\.\d+[^\s<]*")


def _script() -> ModuleType:
    """Import `scripts/visual/docs_figures.py` from its path."""
    path = ROOT / "scripts" / "visual" / "docs_figures.py"
    spec = importlib.util.spec_from_file_location("docs_figures", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _normal(text: str) -> str:
    """The text without the version of pkpdsim (the credit of the interactive pages)."""
    return VERSION.sub("pkpdsim</a> VERSION", text)


def _points(d: str) -> np.ndarray:
    """The vertices of path data."""
    return np.array(NUMBER.findall(d), dtype=float).reshape(-1, 2)


def _hausdorff(a: np.ndarray, b: np.ndarray) -> float:
    """The symmetric Hausdorff distance of two sets of vertices."""
    if len(a) == 0 or len(b) == 0:
        return 0.0 if len(a) == len(b) else np.inf
    return max(KDTree(b).query(a)[0].max(), KDTree(a).query(b)[0].max())


def _same(fresh: str, stored: str) -> bool:
    """Whether two figures agree within the rounding of another platform.

    The text without the numbers is the same, the numbers agree within
    `TOLERANCE`, and the outlines of a body, which a flipped cell of its raster
    changes in their vertices, within `TOLERANCE` as sets of points.
    """
    if fresh == stored:
        return True
    fresh_d, stored_d = GEOMETRY.findall(fresh), GEOMETRY.findall(stored)
    if len(fresh_d) != len(stored_d):
        return False
    fresh_rest, stored_rest = GEOMETRY.sub("", fresh), GEOMETRY.sub("", stored)
    if NUMBER.sub("#", fresh_rest) != NUMBER.sub("#", stored_rest):
        return False
    numbers = zip(NUMBER.findall(fresh_rest), NUMBER.findall(stored_rest), strict=True)
    if any(abs(float(a) - float(b)) > TOLERANCE for a, b in numbers):
        return False
    for a, b in zip(fresh_d, stored_d, strict=True):
        if a == b:
            continue
        if not (POLYLINES.fullmatch(a) and POLYLINES.fullmatch(b)):
            if NUMBER.sub("#", a) != NUMBER.sub("#", b):
                return False
            pairs = zip(NUMBER.findall(a), NUMBER.findall(b), strict=True)
            if any(abs(float(x) - float(y)) > TOLERANCE for x, y in pairs):
                return False
        elif _hausdorff(_points(a), _points(b)) > TOLERANCE:
            return False
    return True


def test_a_moved_cell_of_an_outline_is_the_same_figure() -> None:
    stored = '<path d="M0.0,0.0 1.0,0.0 1.0,1.0 0.0,1.0Z" transform="scale(1.20)"/>'
    moved = (
        '<path d="M0.0,0.0 0.2,0.1 1.0,0.0 1.2,1.0 0.0,1.0Z" transform="scale(1.30)"/>'
    )
    assert _same(moved, stored)
    assert not _same(moved.replace("1.2,1.0", "2.0,1.0"), stored)
    assert not _same(moved.replace("1.30", "1.80"), stored)
    assert not _same(moved.replace("scale", "rotate"), stored)
    assert not _same(stored + "<g/>", stored)


FIGURES = [
    "anatomy.html",
    "schematic.html",
    "publication.svg",
    "publication_panels.svg",
    "individuals.svg",
    "mouse.html",
    "rat.html",
    "pig.html",
]


@pytest.fixture(scope="module")
def fresh(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, dict[Path, str]]:
    """The figures and the pages of the models, written once into a folder."""
    out = tmp_path_factory.mktemp("docs_figures")
    script = _script()
    script.write(out)
    return out, script.write_pages(out)


@pytest.mark.parametrize("name", FIGURES)
def test_the_figures_of_the_docs_are_up_to_date(
    name: str, fresh: tuple[Path, dict[Path, str]]
) -> None:
    out, _ = fresh
    stored = _normal((_script().OUT / name).read_text(encoding="utf-8"))
    assert _same(_normal((out / name).read_text(encoding="utf-8")), stored), (
        f"{name} is out of date: run `uv run python scripts/visual/docs_figures.py`"
    )


def test_the_pages_of_the_models_are_up_to_date(
    fresh: tuple[Path, dict[Path, str]],
) -> None:
    out, pages = fresh
    script = _script()
    assert pages, "no page of the PBPK library creates a model"
    command = "run `uv run python scripts/visual/docs_figures.py`"
    for page, text in pages.items():
        assert _same(text, page.read_text(encoding="utf-8")), f"{page.name}: {command}"
    written = sorted(p.name for p in (out / script.MODELS).glob("*.html"))
    stored = sorted(p.name for p in (script.OUT / script.MODELS).glob("*.html"))
    assert written == stored, command
    for name in written:
        fresh_text = _normal((out / script.MODELS / name).read_text(encoding="utf-8"))
        stored_text = _normal(
            (script.OUT / script.MODELS / name).read_text(encoding="utf-8")
        )
        assert _same(fresh_text, stored_text), f"{name}: {command}"


def test_every_model_of_a_page_shows_its_figures() -> None:
    script = _script()
    for page in script.model_pages():
        text = page.read_text(encoding="utf-8")
        for number in script.model_fences(page):
            stem = f"{page.stem}-{number}"
            assert f"figures-{stem} " in text, page.name
            for figure in ("anatomy", "schematic"):
                assert f"img/figures/models/{stem}-{figure}.html" in text, page.name


def test_the_perfusion_page_shows_the_figures_of_its_circuit() -> None:
    """A fence which creates a machine perfusion gets the figures of its circuit."""
    script = _script()
    page = ROOT / "docs" / "perfusion.md"
    assert page in script.model_pages()
    definitions = script.model_fences(page)
    assert definitions, "the page creates no circuit"
    text = page.read_text(encoding="utf-8")
    for number, definition in definitions.items():
        assert isinstance(definition, PerfusionDefinition)
        for figure in ("anatomy", "schematic"):
            # relative to the page, the folder of the figures of the models
            src = f"pbpk/img/figures/models/perfusion-{number}-{figure}.html"
            assert f'src="{src}"' in text
            assert (ROOT / "docs" / src).exists()
            html = (ROOT / "docs" / src).read_text(encoding="utf-8")
            assert definition.name in html or definition.sid in html


def test_the_docs_show_the_interactive_figures_only() -> None:
    for page in (ROOT / "docs").rglob("*.md"):
        text = page.read_text(encoding="utf-8")
        for static in (
            "anatomy.svg",
            "schematic.svg",
            "mouse.svg",
            "rat.svg",
            "pig.svg",
        ):
            assert f"img/figures/{static}" not in text, page.name


def test_the_page_shows_every_figure() -> None:
    page = (ROOT / "docs" / "pbpk" / "figures.md").read_text(encoding="utf-8")
    for name in (
        "anatomy.html",
        "schematic.html",
        "publication.svg",
        "individuals.svg",
    ):
        assert f"img/figures/{name}" in page, name
    species = (ROOT / "docs" / "pbpk" / "species.md").read_text(encoding="utf-8")
    for name in ("mouse.html", "rat.html", "pig.html"):
        assert f"img/figures/{name}" in species, name
