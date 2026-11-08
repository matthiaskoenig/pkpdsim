"""Tests of the figure of a machine perfusion circuit, the SVG schematic of a definition.

The figure is checked as a drawing, as the figure of the PBPK library
(`svg_labels`): valid and deterministic SVG, every label inside the canvas
and inside its box, no two labels overlapping; the circuit (reservoir, pump,
oxygenator, tubing), the organ with its grey icon and its zones, the
perfusate quantities as chips, bile or urine, and a hover card for every
key of an interactive page.
"""

import dataclasses
import re
from pathlib import Path
from xml.etree import ElementTree

import pytest
from svg_labels import SVG, assert_inside_and_apart, assert_inside_box, parse_labels

from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.pbpk.visual import palette
from pkpdsim.pbpk.visual.interactive import CREDIT, ICON_CREDIT
from pkpdsim.perfusion import ids
from pkpdsim.perfusion.builder import compose
from pkpdsim.perfusion.circuit import quantities
from pkpdsim.perfusion.definition import GenericOrgan, PerfusionDefinition
from pkpdsim.perfusion.figure import (
    circuit_cards,
    circuit_figure,
    circuit_svg,
    interactive_circuit,
    process_groups,
)
from pkpdsim.perfusion.units import PU

from .conftest import (
    _definition,
    generic_definition,
    kidney_definition,
    liver_definition,
)

PROFILES = {"OXPHOS": (1.5, 1.5, 1.25, 1.0, 1.0, 0.75, 0.5, 0.5), "GLC2PYR": (1.0,) * 8}


def _unknown_organ() -> PerfusionDefinition:
    """An organ without an icon and without a table, its fractions given."""
    fraction = Value(0.1, U.dimensionless, "ModelChoice", "test")
    organ = GenericOrgan(
        "xo",
        Value(2.0, PU.g, "ModelChoice", "test"),
        3,
        vascular_fraction=fraction,
        cell_fraction=Value(0.7, U.dimensionless, "ModelChoice", "test"),
        density=Value(1.05, U.kg_per_l, "ModelChoice", "test"),
    )
    return _definition(
        organ, "generic", tracers=(), sid="perfusion_xo", name="perfused xo"
    )


def _definitions() -> dict[str, PerfusionDefinition]:
    """The definitions of the figure tests by name."""
    return {
        "liver1": liver_definition(),
        "liver8": liver_definition(8, profiles=PROFILES, tracers=["inulin"]),
        "liver30": liver_definition(30),
        "kidney": kidney_definition(tracers=["inulin", "fitc"]),
        "kidney3": kidney_definition(3),
        "muscle": generic_definition(),
        "unknown": _unknown_organ(),
        # a name longer than a line of the legend wraps
        "long_name": dataclasses.replace(
            liver_definition(4),
            name="rat liver in a normothermic machine perfusion with two medium changes",
        ),
    }


@pytest.fixture(params=list(_definitions()))
def figure(
    request: pytest.FixtureRequest, tmp_path: Path
) -> tuple[PerfusionDefinition, Path]:
    """The figure of a definition of `_definitions`."""
    d = _definitions()[request.param]
    return d, circuit_figure(d, tmp_path / f"{d.sid}.svg")


def _keys(root: ElementTree.Element) -> list[str]:
    """The hover keys of the elements of a figure, in their order."""
    return [e.attrib["data-k"] for e in root.iter() if "data-k" in e.attrib]


def _group(root: ElementTree.Element, key: str) -> ElementTree.Element:
    """The group of a key, the first element which carries it."""
    for g in root.iter(f"{SVG}g"):
        if g.attrib.get("data-k") == key:
            return g
    raise KeyError(key)


def _texts(g: ElementTree.Element) -> list[str]:
    """The texts of a group."""
    return [t.text or "" for t in g.iter(f"{SVG}text")]


def test_the_figure_is_valid_and_deterministic(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """An SVG of the size of its view box with a title, the same text for the same definition."""
    d, path = figure
    root, _ = parse_labels(path)
    assert root.tag == f"{SVG}svg"
    width, height = root.attrib["width"], root.attrib["height"]
    assert root.attrib["viewBox"] == f"0 0 {width} {height}"
    assert float(width) > 0 and float(height) > 0
    title = root.find(f"{SVG}title")
    assert title is not None and title.text == d.name
    assert circuit_svg(d) == path.read_text(encoding="utf-8")
    assert circuit_svg(d) == circuit_svg(d)


def test_the_labels_are_inside_the_canvas_and_do_not_overlap(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """Every label box lies in the canvas and no two label boxes overlap."""
    _, path = figure
    root, labels = parse_labels(path)
    assert_inside_and_apart(root, labels)


def test_the_labels_of_a_box_are_inside_it(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """The labels of a row, the reservoir and a compartment lie inside its rectangle."""
    _, path = figure
    root, _ = parse_labels(path)
    boxes = [
        g
        for g in root.iter(f"{SVG}g")
        if g.attrib.get("class") in ("pkf-box", "pkf-sink")
    ]
    assert len(boxes) >= 5
    for g in boxes:
        assert_inside_box(g)


def test_the_circuit_is_drawn(figure: tuple[PerfusionDefinition, Path]) -> None:
    """The reservoir, the pump with the flow, the oxygenator, the tubing and the temperature."""
    d, path = figure
    root, _ = parse_labels(path)
    keys = _keys(root)
    organ = d.organ.key
    for key in ("reservoir", "pump", "inflow", "gas", organ, "outflow", "temperature"):
        assert key in keys, key
    # in the order of the flow: reservoir, pump, inflow tubing, organ, outflow tubing
    order = [k for k in keys if k in ("pump", "inflow", organ, "outflow")]
    assert order == ["pump", "inflow", organ, "outflow"]
    assert "1000 ml" in _texts(_group(root, "reservoir"))
    assert any("3 ml/min" in t for t in _texts(_group(root, "pump")))
    gas = _texts(_group(root, "gas"))
    assert gas == ["sweep gas", "o2", "co2"]
    assert any("37 °C" in t for t in _texts(_group(root, "temperature")))


def test_the_organ_icon_is_grey(figure: tuple[PerfusionDefinition, Path]) -> None:
    """The icon of the organ has no colour: its fills are greys, its outline the ink."""
    d, path = figure
    root, _ = parse_labels(path)
    icon = _group(root, d.organ.key).find(f"{SVG}svg")
    if d.organ.key == "xo":
        # an organ without an icon of the anatomy: grey cells on the circle
        assert icon is None
        assert "organ xo" in _texts(_group(root, "xo"))
        cells = [
            e
            for e in _group(root, "xo").iter(f"{SVG}circle")
            if e.attrib.get("stroke") == palette.ICON_INK
        ]
        assert len(cells) == 3
        return
    assert icon is not None
    paths = list(icon.iter(f"{SVG}path"))
    assert paths
    for p in paths:
        fill = p.attrib["fill"].lstrip("#")
        fill = fill if len(fill) == 6 else "".join(c * 2 for c in fill)
        assert fill[0:2] == fill[2:4] == fill[4:6], p.attrib["fill"]
        assert p.attrib["stroke"] == palette.ICON_INK


def test_the_zones_are_drawn_once_each(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """The organ shows its zones as a strip, the first zone first, each zone once."""
    d, path = figure
    root, _ = parse_labels(path)
    zones = [k for k in _keys(root) if k.startswith("zone:")]
    assert zones == [f"zone:{ids.zone_id(k)}" for k in range(1, d.organ.zones + 1)]
    texts = _texts(_group(root, d.organ.key))
    first, last = (
        ("periportal", "pericentral")
        if d.organ.key == "li"
        else (
            "arterial",
            "venous",
        )
    )
    assert any(t.startswith(first) for t in texts)
    assert any(t.endswith(last) for t in texts)


def test_the_perfusate_quantities_are_chips(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """Every quantity of the perfusate is a chip of the reservoir and of the organ.

    A quantity of the perfusate which the cells hold as well is a filled chip
    in the organ, one which stays in the perfusate (CO2, bicarbonate, the
    protons, a tracer) a dashed chip.
    """
    d, path = figure
    root, _ = parse_labels(path)
    expected = list(quantities(d.perfusate, d.released))
    assert [t for t in _texts(_group(root, "reservoir")) if t in expected] == expected
    organ = _group(root, d.organ.key)
    # the chips of the row, before those of the urine
    chips: dict[str, ElementTree.Element] = {}
    for g in organ.iter(f"{SVG}g"):
        if g.attrib.get("data-k", "").startswith("sub:"):
            chips.setdefault(g.attrib["data-k"].removeprefix("sub:"), g)
    assert list(chips) == expected
    for q, g in chips.items():
        rect = g.find(f"{SVG}rect")
        assert rect is not None
        dashed = "stroke-dasharray" in rect.attrib["style"]
        assert dashed == (q not in ids.CELL), q


def test_bile_or_urine_is_a_compartment(
    figure: tuple[PerfusionDefinition, Path],
) -> None:
    """The liver secretes into the bile, the kidney into the urine with the filtered quantities."""
    d, path = figure
    root, _ = parse_labels(path)
    sinks = {
        g.attrib["data-k"]: _texts(g)
        for g in root.iter(f"{SVG}g")
        if g.attrib.get("class") == "pkf-sink"
    }
    if d.organ.key == "li":
        assert list(sinks) == ["gas", "bile"]
    elif d.organ.key == "ki":
        assert list(sinks) == ["gas", "urine"]
        filtered = [*ids.FILTERED, *d.perfusate.tracers]
        assert sinks["urine"] == ["urine", *filtered]
    else:
        assert list(sinks) == ["gas"]


def test_every_process_of_the_core_is_drawn() -> None:
    """Every reaction of the cells of a zone belongs to a group of a badge of the figure."""
    for d in (liver_definition(), kidney_definition(), generic_definition()):
        sids = [r.sid or "" for r in compose(d).reactions or []]
        reactions = [
            sid.removesuffix("_z1")
            for sid in sids
            if sid.endswith("_z1") and not sid.startswith(("FLOW_", "BUF_"))
        ]
        grouped = {sid for rows in process_groups(d).values() for sid, _ in rows}
        missing = [r for r in reactions if r not in grouped]
        assert not missing, missing


def test_the_profiles_are_drawn(tmp_path: Path) -> None:
    """A zonal profile is a line of the organ and part of the card of each zone."""
    d = liver_definition(8, profiles=PROFILES)
    _, labels = parse_labels(circuit_figure(d, tmp_path / "p.svg"))
    text = " ".join(label.text for label in labels)
    assert "OXPHOS 1.5 (z1) to 0.5 (z8)" in text
    cards = circuit_cards(d)
    assert any("p_OXPHOS_z1</code> 1.5" in row for row in cards["zone:z1"])


def test_every_key_has_a_card(figure: tuple[PerfusionDefinition, Path]) -> None:
    """Every key of the figure but a process line (`txt:`) has a hover card."""
    d, path = figure
    root, _ = parse_labels(path)
    cards = circuit_cards(d)
    for key in _keys(root):
        if not key.startswith("txt:"):
            assert key in cards, key


def test_the_cards_are_escaped_once() -> None:
    """A placeholder of a card is escaped once, `alt_&lt;site&gt;`."""
    d = liver_definition(2, tracers=["inulin"])
    cards = circuit_cards(d)
    rows = [row for card in cards.values() for row in card]
    assert not [row for row in rows if "&amp;lt;" in row or "&amp;gt;" in row]
    assert any("<code>alt_&lt;site&gt;</code>" in row for row in cards["sub:alt"])


def test_the_interactive_page(tmp_path: Path) -> None:
    """The page carries the figure, its cards and the credit."""
    d = liver_definition(8, tracers=["inulin"])
    page = interactive_circuit(
        d, tmp_path / "circuit.html", figure="schematic"
    ).read_text(encoding="utf-8")
    assert page.startswith("<!doctype html>")
    assert circuit_svg(d) in page
    assert '"zone:z8"' in page
    assert "pkpdFigure" in page
    width = re.search(r'<svg[^>]* width="([0-9.]+)"', page)
    assert width is not None
    assert f"max-width:{int(float(width.group(1)))}px" in page
    # the schematic draws no anatomy, its organ icons are the anatomograms
    assert ICON_CREDIT in page
    assert CREDIT not in page
    anatomy = interactive_circuit(d, tmp_path / "anatomy.html").read_text(
        encoding="utf-8"
    )
    assert CREDIT in anatomy


def test_an_invalid_definition_has_no_figure() -> None:
    """The figure and the cards refuse a definition which does not build."""
    d = liver_definition(sid="no id")
    with pytest.raises(ValueError, match="is no SId"):
        circuit_svg(d)
    with pytest.raises(ValueError, match="is no SId"):
        circuit_cards(d)
