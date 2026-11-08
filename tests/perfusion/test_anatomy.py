"""Tests of the anatomical figure of a machine perfusion, the organ in its circuit (panel A).

The figure is checked as a drawing (`svg_labels`): valid and deterministic
SVG, every label inside the canvas, no two labels overlapping; the organ in
its drawn shape and tissue colour, the apparatus with grey icons, the zones of
the inset and a hover card for every key of the interactive page.
"""

import re
from pathlib import Path
from xml.etree import ElementTree

import pytest
from svg_labels import SVG, assert_inside_and_apart, parse_labels

from pkpdsim.pbpk.visual import palette
from pkpdsim.perfusion import ids
from pkpdsim.perfusion.anatomy import (
    NEUTRAL,
    anatomy_cards,
    circuit_anatomy_figure,
    circuit_anatomy_svg,
)
from pkpdsim.perfusion.definition import KidneyOrgan, LiverOrgan, PerfusionDefinition
from pkpdsim.perfusion.figure import interactive_circuit

from .conftest import generic_definition, kidney_definition, liver_definition
from .test_figure import PROFILES, _unknown_organ


def _definitions() -> dict[str, PerfusionDefinition]:
    """The definitions of the anatomy tests by name."""
    return {
        "liver1": liver_definition(),
        "liver8": liver_definition(8, profiles=PROFILES, tracers=["inulin"]),
        "liver30": liver_definition(30),
        "kidney": kidney_definition(tracers=["inulin", "fitc"]),
        "kidney3": kidney_definition(3),
        "muscle": generic_definition(),
        "unknown": _unknown_organ(),
    }


@pytest.fixture(params=list(_definitions()))
def figure(
    request: pytest.FixtureRequest, tmp_path: Path
) -> tuple[PerfusionDefinition, Path]:
    """The anatomical figure of a definition of `_definitions`."""
    d = _definitions()[request.param]
    return d, circuit_anatomy_figure(d, tmp_path / f"{d.sid}_anatomy.svg")


def _keys(root: ElementTree.Element) -> list[str]:
    """The hover keys of the elements of a figure, in their order."""
    return [e.attrib["data-k"] for e in root.iter() if "data-k" in e.attrib]


def _groups(root: ElementTree.Element, key: str) -> list[ElementTree.Element]:
    """The groups of a key."""
    return [g for g in root.iter(f"{SVG}g") if g.attrib.get("data-k") == key]


def _grey(color: str) -> bool:
    """Whether a colour is a grey (or white)."""
    c = color.lstrip("#")
    c = c if len(c) == 6 else "".join(x * 2 for x in c)
    return c[0:2] == c[2:4] == c[4:6]


def test_the_figure_is_valid_and_deterministic(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """An SVG of the size of its view box with a title, the same text for the same definition."""
    d, path = figure
    root, _ = parse_labels(path)
    assert root.tag == f"{SVG}svg"
    width, height = root.attrib["width"], root.attrib["height"]
    assert root.attrib["viewBox"] == f"0 0 {width} {height}"
    title = root.find(f"{SVG}title")
    assert title is not None and title.text == d.name
    assert circuit_anatomy_svg(d) == path.read_text(encoding="utf-8")


def test_the_labels_are_inside_the_canvas_and_do_not_overlap(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """Every label box lies in the canvas and no two label boxes overlap."""
    _, path = figure
    root, labels = parse_labels(path)
    assert_inside_and_apart(root, labels)


def test_the_organ_is_drawn_in_its_shape_and_colour(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """The organ in the shape of the anatomy and its tissue colour, or a neutral shape."""
    d, path = figure
    root, _ = parse_labels(path)
    shapes = [g for g in _groups(root, d.organ.key) if g.attrib["class"] == "pka-organ"]
    assert len(shapes) == 1
    if d.organ.key == "xo":
        rect = shapes[0].find(f"{SVG}rect")
        assert rect is not None and f"fill:{NEUTRAL}" in rect.attrib["style"]
        return
    shape = shapes[0].find(f"{SVG}path")
    assert shape is not None
    assert shape.attrib["fill"] == palette.ORGAN_COLORS[d.organ.key]
    # a drawn outline of many points, not a box
    assert shape.attrib["d"].count(",") > 10


def test_the_apparatus_and_its_labels(figure: tuple[PerfusionDefinition, Path]) -> None:
    """The reservoir, the pump, the oxygenator and its gas, the tubing, the collector and the labels."""
    d, path = figure
    root, labels = parse_labels(path)
    keys = _keys(root)
    for key in (
        "reservoir",
        "pump",
        "inflow",
        "gas",
        "outflow",
        "medium",
        "temperature",
    ):
        assert key in keys, key
    for key in ("reservoir", "pump", "inflow", d.organ.key):
        # the drawn element and its label
        assert len(_groups(root, key)) == 2, key
    texts = [label.text for label in labels]
    vein = {"li": "portal vein", "ki": "renal artery"}.get(d.organ.key, "artery")
    assert vein in texts
    assert "1000 ml" in texts
    assert "temperature 37 °C" in texts
    sink = (
        "bile"
        if isinstance(d.organ, LiverOrgan)
        else "urine"
        if isinstance(d.organ, KidneyOrgan)
        else None
    )
    for name in ("bile", "urine"):
        assert (name in keys) == (name == sink), name


def test_the_icons_are_grey(figure: tuple[PerfusionDefinition, Path]) -> None:
    """The icons of the labels and of the apparatus have no colour."""
    _, path = figure
    root, _ = parse_labels(path)
    for g in root.iter(f"{SVG}g"):
        if g.attrib.get("class") not in ("pka-label", "pka-pump", "pka-inflow"):
            continue
        # the nested drawings of the organ icons and the drawn icons on their rings
        for e in g.iter():
            if (
                e.tag in (f"{SVG}path", f"{SVG}circle", f"{SVG}rect")
                and "fill" in e.attrib
            ) and e.attrib["fill"] != "none":
                assert _grey(e.attrib["fill"]), (g.attrib["data-k"], e.attrib)


def test_the_zones_of_the_inset(figure: tuple[PerfusionDefinition, Path]) -> None:
    """The inset shows the zones in the order of the flow: always for the liver, else from two zones."""
    d, path = figure
    root, labels = parse_labels(path)
    n = d.organ.zones
    shown = isinstance(d.organ, LiverOrgan) or n > 1
    zones = [k for k in _keys(root) if k.startswith("zone:")]
    assert zones == (
        [f"zone:{ids.zone_id(k)}" for k in range(1, n + 1)] if shown else []
    )
    if shown:
        texts = [label.text for label in labels]
        first = "periportal" if isinstance(d.organ, LiverOrgan) else "arterial"
        assert any(t.startswith(first) for t in texts)


def test_the_chips_of_the_organ(figure: tuple[PerfusionDefinition, Path]) -> None:
    """The chips of the organ label: filled where the cells hold the quantity, dashed else."""
    d, path = figure
    root, _ = parse_labels(path)
    label = next(
        g for g in _groups(root, d.organ.key) if g.attrib["class"] == "pka-label"
    )
    for g in label.iter(f"{SVG}g"):
        key = g.attrib.get("data-k", "")
        if key.startswith("sub:"):
            rect = g.find(f"{SVG}rect")
            assert rect is not None
            q = key.removeprefix("sub:")
            assert ("stroke-dasharray" in rect.attrib["style"]) == (q not in ids.CELL)


def test_every_key_has_a_card(figure: tuple[PerfusionDefinition, Path]) -> None:
    """Every key of the figure has a hover card."""
    d, path = figure
    root, _ = parse_labels(path)
    cards = anatomy_cards(d)
    for key in _keys(root):
        assert key in cards, key


def test_the_interactive_page_of_the_anatomy(tmp_path: Path) -> None:
    """The anatomical figure is the default of the interactive page; another figure is refused."""
    d = liver_definition(8)
    page = interactive_circuit(d, tmp_path / "a.html").read_text(encoding="utf-8")
    assert circuit_anatomy_svg(d) in page
    assert '"medium"' in page
    width = re.search(r'<svg[^>]* width="([0-9.]+)"', page)
    assert width is not None
    assert f"max-width:{int(float(width.group(1)))}px" in page
    with pytest.raises(ValueError, match="No interactive figure"):
        interactive_circuit(d, tmp_path / "b.html", figure="body")
