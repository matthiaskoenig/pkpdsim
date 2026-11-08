"""Tests of the anatomical figure of a model and of the interactive pages."""

import json
import re
from pathlib import Path
from xml.etree import ElementTree

from pkpdsim.pbpk import Route
from pkpdsim.pbpk.visual.anatomy import TITLES, AnatomyFigure, anatomy_figure
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.interactive import interactive_figure

from .conftest import definition
from .test_figure import seven_substances

SVG = "{http://www.w3.org/2000/svg}"
ORGANS = (
    "br",
    "lu",
    "he",
    "li",
    "sk",
    "ki",
    "la",
    "st",
    "sp",
    "pa",
    "si",
    "ad",
    "bo",
    "mu",
    "go",
)


def _texts(path: Path) -> list[str]:
    root = ElementTree.parse(path).getroot()
    return ["".join(t.itertext()) for t in root.iter(f"{SVG}text")]


def test_fig1_shows_both_builds(tmp_path: Path) -> None:
    texts = _texts(anatomy_figure(definition(), tmp_path / "a.svg"))
    assert "male build" in texts and "female build" in texts


def test_every_organ_has_a_label_and_an_icon(tmp_path: Path) -> None:
    path = anatomy_figure(definition(), tmp_path / "a.svg")
    texts = " ".join(_texts(path))
    root = ElementTree.parse(path).getroot()
    labels = {
        g.attrib["data-k"]: g
        for g in root.iter(f"{SVG}g")
        if g.attrib.get("data-k") in ORGANS and "transform" in g.attrib
    }
    for key in ORGANS:
        assert TITLES[key] in texts, key
        assert key in labels, key
        # the circled icon: a ring and the drawing of the organ
        assert labels[key].find(f".//{SVG}circle") is not None, key


def test_organs_with_processes_show_their_submodel() -> None:
    figure = AnatomyFigure(ModelContent.of(seven_substances()))
    for key in ("li", "ki", "si"):
        assert figure.has_submodel(key), key
        assert "interstitium" in figure.label(key, 0, 0, 310, 1)[0]
    assert not figure.has_submodel("br")
    assert "interstitium" not in figure.label("br", 0, 0, 310, 1)[0]


def test_the_figure_has_no_credit(tmp_path: Path) -> None:
    svg = anatomy_figure(seven_substances(), tmp_path / "a.svg").read_text()
    assert "Created with" not in svg and "CC BY" not in svg


def test_fig1_is_deterministic(tmp_path: Path) -> None:
    a = anatomy_figure(definition(), tmp_path / "a.svg").read_bytes()
    b = anatomy_figure(definition(), tmp_path / "b.svg").read_bytes()
    assert a == b


def test_an_iv_only_definition_has_no_oral_route(tmp_path: Path) -> None:
    d = definition(sid="iv_only", routes=[Route.iv("drg")])
    svg = anatomy_figure(d, tmp_path / "a.svg").read_text()
    assert 'data-k="route:oral"' not in svg
    assert 'data-k="route:iv"' in svg
    assert ">oral<" not in svg


def test_a_female_definition_shows_both_builds(tmp_path: Path) -> None:
    texts = _texts(anatomy_figure(definition(sex="female"), tmp_path / "a.svg"))
    assert "male build" in texts and "female build" in texts


def test_the_interactive_page_has_cards_and_the_credit(tmp_path: Path) -> None:
    html = interactive_figure(seven_substances(), tmp_path / "a.html").read_text()
    match = re.search(r"const CARDS = (\{.*?\});\n", html, re.S)
    assert match is not None
    cards = json.loads(match.group(1))
    assert any("UBERON:0002107" in row for row in cards["li"])
    assert any("reference man" in row for row in cards["li"])
    assert "sub:api" in cards
    figure, legend = html.split('<p class="legend">', 1)
    assert "Created with" not in figure
    assert legend.count("https://github.com/matthiaskoenig/pkpdsim") == 1
    assert legend.count("https://matthiaskoenig.github.io/pkpdsim/pbpk/") == 1


def test_the_interactive_schematic_has_process_cards(tmp_path: Path) -> None:
    html = interactive_figure(
        seven_substances(), tmp_path / "s.html", figure="schematic"
    ).read_text()
    assert 'data-k="txt:MET: API3A4: api -&gt; api_m1, api_m2 by CYP3A4' in html
    assert 'data-k="li"' in html
