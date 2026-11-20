"""Tests of the anatomical figure of a model and of its interactive figures."""

import json
import re
from pathlib import Path
from xml.etree import ElementTree

import pytest

from pkpdsim.pbpk import Route
from pkpdsim.pbpk.visual import palette
from pkpdsim.pbpk.visual.anatomy import TITLES, AnatomyFigure, anatomy_figure
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.interactive import (
    CREDIT,
    ICON_CREDIT,
    html_page,
    html_snippet,
    interactive_figure,
)

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
    svg = anatomy_figure(seven_substances(), tmp_path / "a.svg").read_text(
        encoding="utf-8"
    )
    assert "Created with" not in svg and "CC BY" not in svg


def test_fig1_is_deterministic(tmp_path: Path) -> None:
    a = anatomy_figure(definition(), tmp_path / "a.svg").read_bytes()
    b = anatomy_figure(definition(), tmp_path / "b.svg").read_bytes()
    assert a == b


def test_an_iv_only_definition_has_no_oral_route(tmp_path: Path) -> None:
    d = definition(sid="iv_only", routes=[Route.iv("drg")])
    svg = anatomy_figure(d, tmp_path / "a.svg").read_text(encoding="utf-8")
    assert 'data-k="route:oral"' not in svg
    assert 'data-k="route:iv"' in svg
    assert ">oral<" not in svg


def test_a_female_definition_shows_both_builds(tmp_path: Path) -> None:
    texts = _texts(anatomy_figure(definition(sex="female"), tmp_path / "a.svg"))
    assert "male build" in texts and "female build" in texts


def test_the_interactive_page_has_cards_and_the_credit(tmp_path: Path) -> None:
    html = interactive_figure(seven_substances(), tmp_path / "a.html").read_text(
        encoding="utf-8"
    )
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
    assert CREDIT in legend


def test_a_card_never_ends_the_script_of_the_page() -> None:
    """The HTML and the line separators of a card are escapes in the script block."""
    rows = ["</script><img src=x onerror=alert(1)>", "a\u2028b", "x & <b>y</b>"]
    html = html_page("<svg/>", {"li": rows}, title="t", width=600, legend="l")
    # the only end of a script element in the page is the end of its block
    assert html.count("</script") == 1
    match = re.search(r"const CARDS = (\{.*?\});\n", html, re.S)
    assert match is not None
    assert json.loads(match.group(1)) == {"li": rows}


def test_the_snippets_of_a_model_are_scoped_to_their_div(tmp_path: Path) -> None:
    """A snippet is one `<div>`: its styles, its ids and its script stay inside it."""
    d = seven_substances()
    snippets = {
        kind: interactive_figure(
            d, tmp_path / f"{kind}.html", figure=kind, snippet=True
        ).read_text(encoding="utf-8")
        for kind in ("anatomy", "schematic")
    }
    ids: dict[str, set[str]] = {}
    for kind, html in snippets.items():
        scope = f"pkf-{d.sid}-{kind}"
        assert html.startswith(f'<div id="{scope}" class="pkf-figure">')
        assert html.rstrip().endswith("</div>")
        assert "<html" not in html and "<body" not in html
        # a blank line would end the snippet in a page of markdown
        assert all(line.strip() for line in html.splitlines())
        for style in re.findall(r"<style>(.*?)</style>", html, re.S):
            selectors = [r.split("{", 1)[0] for r in style.split("}") if "{" in r]
            assert selectors
            assert all(
                s.strip() == f"#{scope}" or s.strip().startswith(f"#{scope} ")
                for selector in selectors
                for s in selector.split(",")
            )
        ids[kind] = set(re.findall(r'\bid="([^"]+)"', html))
        assert ids[kind] and all(i.startswith(scope) for i in ids[kind])
        # the script is a closure over the div, not one of the page
        assert f"document.getElementById('{scope}')" in html
    # the two snippets of a model stand on one page
    assert not ids["anatomy"] & ids["schematic"]
    match = re.search(r"const CARDS = (\{.*?\});\n", snippets["anatomy"], re.S)
    assert match is not None
    assert any("UBERON:0002107" in row for row in json.loads(match.group(1))["li"])
    assert CREDIT in snippets["anatomy"]
    assert ICON_CREDIT in snippets["schematic"]


def test_a_style_sheet_which_cannot_be_scoped_is_refused() -> None:
    """A snippet never leaves a rule of its figure unscoped."""
    with pytest.raises(ValueError, match="at-rule or a comment"):
        html_snippet(
            "<svg><style>@media print{.a{fill:#000}}</style></svg>", {}, "s", 10, "l"
        )


def test_the_interactive_schematic_has_process_cards(tmp_path: Path) -> None:
    html = interactive_figure(
        seven_substances(), tmp_path / "s.html", figure="schematic"
    ).read_text(encoding="utf-8")
    assert 'data-k="txt:MET: API3A4: api -&gt; api_m1, api_m2 by CYP3A4' in html
    assert 'data-k="li"' in html
    # the schematic draws no anatomy, its organ icons are the anatomograms
    assert ICON_CREDIT in html
    assert CREDIT not in html


def test_the_bodies_alone_have_no_annotations(tmp_path: Path) -> None:
    """The bodies with their organs, for a slide: no text, transparent."""
    path = anatomy_figure(definition(), tmp_path / "a.svg", annotations=False)
    assert _texts(path) == []
    svg = path.read_text(encoding="utf-8")
    assert 'fill="#fff"/>' not in svg.split("</style>", 1)[1][:200]
    for organ in ("br", "lu", "he", "li", "ki"):
        assert palette.ORGAN_COLORS[organ] in svg, organ
    annotated = anatomy_figure(definition(), tmp_path / "b.svg")
    width = float(ElementTree.parse(path).getroot().get("width") or 0)
    assert 0 < width < float(ElementTree.parse(annotated).getroot().get("width") or 0)
    # deterministic
    again = anatomy_figure(definition(), tmp_path / "c.svg", annotations=False)
    assert again.read_bytes() == path.read_bytes()
