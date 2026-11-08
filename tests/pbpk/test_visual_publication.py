"""Tests of the publication figure of a model."""

import builtins
import re
from pathlib import Path
from xml.etree import ElementTree

import pytest

from pkpdsim.pbpk.visual.publication import PAGE, publication_figure

from .test_figure import seven_substances

SVG = "{http://www.w3.org/2000/svg}"


@pytest.fixture(scope="module")
def figure(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return publication_figure(
        seven_substances(), tmp_path_factory.mktemp("pub") / "figure.svg"
    )


def test_the_panels_fit_on_one_page(figure: Path) -> None:
    root = ElementTree.parse(figure.with_name("figure_panels.svg")).getroot()
    assert root.attrib["width"] == "210.0mm" and root.attrib["height"] == "297.0mm"
    panels = [e for e in root if e.tag == f"{SVG}svg"]
    assert len(panels) == 2
    for p in panels:
        assert (
            float(p.attrib["y"]) + float(p.attrib["height"]) <= PAGE[1] - PAGE[2] + 1e-6
        )
        assert (
            float(p.attrib["x"]) + float(p.attrib["width"]) <= PAGE[0] - PAGE[2] + 1e-6
        )
    letters = [e.text for e in root if e.tag == f"{SVG}text"]
    assert letters == ["A", "B"]


def test_the_legend_is_below_the_panels(figure: Path) -> None:
    root = ElementTree.parse(figure).getroot()
    legend = [e for e in root if e.tag == f"{SVG}text" and e.text not in ("A", "B")]
    assert legend and all(float(e.attrib["y"]) > PAGE[1] for e in legend)


def test_the_legend_names_the_credit_once(figure: Path) -> None:
    md = figure.with_name("figure_legend.md").read_text()
    assert md.count("[pkpdsim](https://github.com/matthiaskoenig/pkpdsim)") == 1
    assert md.count("https://matthiaskoenig.github.io/pkpdsim/pbpk/") == 1
    assert "Routes: iv api" in md
    svg = figure.read_text()
    assert len(re.findall(r"Created with", svg)) == 1


def test_pdf_needs_the_extra(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    real = builtins.__import__

    def refuse(name: str, *args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        if name == "cairosvg":
            raise ImportError(name)
        return real(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", refuse)
    with pytest.raises(ImportError, match=r"pkpdsim\[pdf\]"):
        publication_figure(seven_substances(), tmp_path / "f.svg", pdf=True)
