"""The organ icons of the figures are grey."""

import re

import numpy as np
import pytest

from pkpdsim.pbpk.visual import icons, palette

FILL = re.compile(r'fill="#([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})"')


def test_every_organ_has_a_grey() -> None:
    assert set(palette.ORGAN_GREYS) == set(palette.ORGAN_COLORS)
    assert len(set(palette.ORGAN_GREYS.values())) == len(palette.ORGAN_GREYS)


@pytest.mark.parametrize("sex", ["male", "female"])
def test_the_icons_have_no_colour(sex: str) -> None:
    for key in icons.organ_keys(sex):
        svg = icons.organ(key, 0.0, 0.0, sex=sex)
        fills = FILL.findall(svg)
        assert fills, key
        assert all(r == g == b for r, g, b in fills), (key, fills)


def test_the_parts_of_the_bowel_differ() -> None:
    assert palette.ORGAN_GREYS["si"] != palette.ORGAN_GREYS["la"]


PATH = re.compile(r'<path d="([^"]+)" fill="([^"]+)"( fill-opacity="0.15")?/>')


def _parts(svg: str) -> dict[str, float]:
    """The area of the kept and the lost part of a pictogram."""
    out = {"kept": 0.0, "lost": 0.0}
    for d, _, pale in PATH.findall(svg):
        points = np.array([[float(v) for v in p.split(",")] for p in d[1:-1].split()])
        out["lost" if pale else "kept"] += icons._area(points)
    return out


@pytest.mark.parametrize("key", ["li", "ki", "he"])
@pytest.mark.parametrize("loss", [0.0, 0.4, 0.81, 1.0])
def test_the_pale_part_of_a_pictogram_is_the_lost_fraction(
    key: str, loss: float
) -> None:
    svg = icons.impaired_organ(key, 50.0, 30.0, 60.0, 40.0, loss, "#2b8cbe", title="B")
    parts = _parts(svg)
    assert parts["lost"] / (parts["kept"] + parts["lost"]) == pytest.approx(
        loss, abs=1e-3
    )
    assert "<title>B</title>" in svg
    # the outline is the color of the theme
    assert 'stroke="currentColor"' in svg


def test_the_pictogram_fills_its_box() -> None:
    svg = icons.impaired_organ("li", 50.0, 30.0, 60.0, 40.0, 0.5, "#2b8cbe")
    outline = re.findall(r'<path d="([^"]+)" fill="none"', svg)[0]
    points = np.array([[float(v) for v in p.split(",")] for p in outline[1:-1].split()])
    (x0, y0), (x1, y1) = points.min(axis=0), points.max(axis=0)
    assert (x0 + x1) / 2 == pytest.approx(50.0, abs=0.01)
    assert (y0 + y1) / 2 == pytest.approx(30.0, abs=0.01)
    assert max(x1 - x0 - 60.0, y1 - y0 - 40.0) == pytest.approx(0.0, abs=0.01)


def test_the_liver_loses_from_the_right_the_kidney_from_the_top() -> None:
    def lost_centre(key: str) -> np.ndarray:
        svg = icons.impaired_organ(key, 0.0, 0.0, 60.0, 60.0, 0.3, "#000000")
        d = PATH.findall(svg)[1][0]
        return np.array(
            [[float(v) for v in p.split(",")] for p in d[1:-1].split()]
        ).mean(axis=0)

    assert lost_centre("li")[0] > 0.0
    assert lost_centre("ki")[1] < 0.0


@pytest.mark.parametrize(("key", "loss"), [("lu", 0.5), ("li", -0.1), ("li", 1.1)])
def test_a_pictogram_refuses_what_it_cannot_draw(key: str, loss: float) -> None:
    with pytest.raises(ValueError):
        icons.impaired_organ(key, 0.0, 0.0, 10.0, 10.0, loss, "#000000")
