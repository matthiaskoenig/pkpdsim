"""The organ icons of the figures are grey."""

import re

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
