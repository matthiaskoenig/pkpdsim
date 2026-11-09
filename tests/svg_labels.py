"""The labels of a figure of the library, checked as a drawing.

A label is a box from its position, its font size (`data-size`) and the
width of a character of the monospace fonts of the figures (`CHAR_WIDTH`);
the figures of the PBPK library (`tests/pbpk/test_figure.py`) and of the
machine perfusion (`tests/perfusion/test_figure.py`) check that every label
lies inside the canvas and inside its box and that no two labels overlap.
"""

from pathlib import Path
from typing import NamedTuple
from xml.etree import ElementTree

from pkpdsim.pbpk.figure import CHAR_WIDTH

__all__ = [
    "SVG",
    "Label",
    "assert_inside_and_apart",
    "assert_inside_box",
    "label_box",
    "overlap",
    "parse_labels",
]

SVG = "{http://www.w3.org/2000/svg}"


class Label(NamedTuple):
    """A text of the figure: its content and box."""

    text: str
    x0: float
    y0: float
    x1: float
    y1: float


def label_box(element: ElementTree.Element) -> Label:
    """The box of a text element from its baseline, font size and length.

    The box spans the ascent (0.8 of the size) above and the descent (0.25)
    below the baseline, the width `CHAR_WIDTH` per character.
    """
    text = element.text or ""
    size = float(element.attrib["data-size"])
    x, y = float(element.attrib["x"]), float(element.attrib["y"])
    assert element.attrib.get("text-anchor", "start") == "start"
    return Label(
        text, x, y - 0.8 * size, x + len(text) * size * CHAR_WIDTH, y + 0.25 * size
    )


def parse_labels(svg: Path) -> tuple[ElementTree.Element, list[Label]]:
    """The root and the labels of a figure."""
    root = ElementTree.parse(svg).getroot()
    return root, [label_box(t) for t in root.iter(f"{SVG}text")]


def overlap(a: Label, b: Label) -> bool:
    """Whether two label boxes overlap."""
    return a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1


def assert_inside_and_apart(root: ElementTree.Element, labels: list[Label]) -> None:
    """Every label box lies in the canvas and no two label boxes overlap."""
    width, height = float(root.attrib["width"]), float(root.attrib["height"])
    for label in labels:
        assert label.x0 >= 0 and label.y0 >= 0, label
        assert label.x1 <= width and label.y1 <= height, label
    for i, a in enumerate(labels):
        for b in labels[i + 1 :]:
            assert not overlap(a, b), (a, b)


def assert_inside_box(g: ElementTree.Element) -> None:
    """The labels of a group lie inside its first rectangle, the box of the group."""
    rect = g.find(f"{SVG}rect")
    assert rect is not None
    x0, y0 = float(rect.attrib["x"]), float(rect.attrib["y"])
    x1, y1 = x0 + float(rect.attrib["width"]), y0 + float(rect.attrib["height"])
    for t in g.iter(f"{SVG}text"):
        label = label_box(t)
        assert x0 < label.x0 and label.x1 < x1, label
        assert y0 < label.y0 and label.y1 < y1, label
