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


PART = re.compile(r'<path d="([^"]+)"[^>]*data-part="(organ|lost|outline)"')


def _points(d: str) -> np.ndarray:
    return np.array([[float(v) for v in p.split(",")] for p in d[1:-1].split()])


def _parts(svg: str) -> dict[str, float]:
    """The area of the organ and of its lost part of a pictogram."""
    out = {"organ": 0.0, "lost": 0.0, "outline": 0.0}
    for d, part in PART.findall(svg):
        out[part] += icons._area(_points(d))
    return out


@pytest.mark.parametrize("key", ["li", "ki", "he"])
@pytest.mark.parametrize("loss", [0.0, 0.4, 0.81, 1.0])
def test_the_faded_part_of_a_pictogram_is_the_lost_fraction(
    key: str, loss: float
) -> None:
    svg = icons.impaired_organ(key, 50.0, 30.0, 60.0, 40.0, loss, "#2b8cbe", title="B")
    parts = _parts(svg)
    assert parts["lost"] / parts["organ"] == pytest.approx(loss, abs=1e-3)
    assert parts["outline"] == pytest.approx(parts["organ"])
    assert "<title>B</title>" in svg
    # the organ in its color, the card framed in the color of the class
    assert f'fill="{palette.ORGAN_COLORS[key]}"' in svg
    assert 'stroke="#2b8cbe"' in svg


def test_the_healthy_card_is_grey() -> None:
    svg = icons.impaired_organ("li", 0.0, 0.0, 60.0, 40.0, 0.0, "black")
    assert 'stroke="#9e9e9e"' in svg
    assert 'data-part="lost"' not in svg


def test_the_pictogram_sits_in_its_card() -> None:
    for key in ("li", "ki", "he"):
        svg = icons.impaired_organ(key, 50.0, 30.0, 60.0, 40.0, 0.5, "#2b8cbe")
        outline = _points(
            re.findall(r'<path d="([^"]+)"[^>]*data-part="outline"', svg)[0]
        )
        (x0, y0), (x1, y1) = outline.min(axis=0), outline.max(axis=0)
        assert 20.0 < x0 < x1 < 80.0, key
        assert 10.0 < y0 < y1 < 50.0, key


def test_the_ids_of_two_pictograms_differ() -> None:
    a = icons.impaired_organ("li", 0.0, 0.0, 60.0, 40.0, 0.4, "#2b8cbe")
    b = icons.impaired_organ("li", 100.0, 0.0, 60.0, 40.0, 0.7, "#045a8d")
    ids = re.compile(r'id="([^"]+)"')
    assert ids.findall(a) != ids.findall(b)
    assert a == icons.impaired_organ("li", 0.0, 0.0, 60.0, 40.0, 0.4, "#2b8cbe")


def test_the_kidney_has_its_vessels_at_the_hilum() -> None:
    svg = icons.impaired_organ("ki", 0.0, 0.0, 60.0, 60.0, 0.3, "#41ae76")
    for color in (palette.BLOOD, palette.URINE):
        assert f'stroke="{color}"' in svg


def test_the_liver_loses_from_the_right_the_kidney_from_the_top() -> None:
    def lost_centre(key: str) -> np.ndarray:
        svg = icons.impaired_organ(key, 0.0, 0.0, 60.0, 60.0, 0.3, "#000000")
        d = re.findall(r'<path d="([^"]+)"[^>]*data-part="lost"', svg)[0]
        return _points(d).mean(axis=0)

    assert lost_centre("li")[0] > 0.0
    assert lost_centre("ki")[1] < 0.0


@pytest.mark.parametrize(("key", "loss"), [("lu", 0.5), ("li", -0.1), ("li", 1.1)])
def test_a_pictogram_refuses_what_it_cannot_draw(key: str, loss: float) -> None:
    with pytest.raises(ValueError):
        icons.impaired_organ(key, 0.0, 0.0, 10.0, 10.0, loss, "#000000")
