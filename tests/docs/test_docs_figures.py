"""The generated figures of the docs page `docs/pbpk/figures.md` are up to date."""

import importlib.util
import re
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).parents[2]
NUMBER = re.compile(r"-?\d+\.\d+")
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


@pytest.mark.parametrize(
    "name",
    [
        "anatomy.svg",
        "anatomy.html",
        "schematic.svg",
        "schematic.html",
        "publication.svg",
        "publication_panels.svg",
        "individuals.svg",
        "mouse.svg",
        "rat.svg",
        "pig.svg",
    ],
)
def test_the_figures_of_the_docs_are_up_to_date(
    name: str, tmp_path_factory: pytest.TempPathFactory
) -> None:
    script = _script()
    out = tmp_path_factory.getbasetemp() / "docs_figures"
    if not (out / name).exists():
        script.write(out)
    fresh = _normal((out / name).read_text(encoding="utf-8"))
    stored = _normal((script.OUT / name).read_text(encoding="utf-8"))
    assert fresh == stored, (
        f"{name} is out of date: run `uv run python scripts/visual/docs_figures.py`"
    )


def test_the_page_shows_every_figure() -> None:
    page = (Path(__file__).parents[2] / "docs" / "pbpk" / "figures.md").read_text(
        encoding="utf-8"
    )
    for name in (
        "anatomy.svg",
        "anatomy.html",
        "schematic.svg",
        "schematic.html",
        "publication.svg",
        "individuals.svg",
    ):
        assert f"img/figures/{name}" in page, name
