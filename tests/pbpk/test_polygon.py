"""The point in polygon test of the outlines, against the answers of matplotlib.

The expected values are those of `matplotlib.path.Path(polygon).contains_points`
(matplotlib 3.11), which the tests of the outlines used before.
"""

import math

import numpy as np
import pytest

from .polygon import contains_points

SQUARE = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]

#: an L of the unit square without its upper right quarter
ELL = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.5), (0.5, 0.5), (0.5, 1.0), (0.0, 1.0)]

#: a pentagram, whose edges cross: its inner pentagon is crossed twice
PENTAGRAM = [
    (
        math.cos(math.pi / 2 + k * 4 * math.pi / 5),
        math.sin(math.pi / 2 + k * 4 * math.pi / 5),
    )
    for k in range(5)
]


@pytest.mark.parametrize(
    ("polygon", "points", "inside"),
    [
        (
            SQUARE,
            [(0.5, 0.5), (0.1, 0.9), (1.5, 0.5), (-0.1, 0.5), (0.5, 2.0)],
            [True, True, False, False, False],
        ),
        # the first vertex repeated at the end, as a closed path writes it
        ([*SQUARE, SQUARE[0]], [(0.5, 0.5), (2.0, 0.5)], [True, False]),
        # the notch of the L is outside, a ray through its corner crosses once
        (
            ELL,
            [(0.25, 0.75), (0.75, 0.25), (0.75, 0.75), (0.25, 0.5)],
            [True, True, False, True],
        ),
        # even-odd: the centre of a pentagram is outside, a point of it inside
        (PENTAGRAM, [(0.0, 0.0), (0.0, 0.8), (0.0, -0.9)], [False, True, False]),
    ],
)
def test_the_answers_of_matplotlib(
    polygon: list[tuple[float, float]],
    points: list[tuple[float, float]],
    inside: list[bool],
) -> None:
    assert contains_points(polygon, points).tolist() == inside


def test_one_point_and_no_points() -> None:
    assert contains_points(SQUARE, (0.5, 0.5)).tolist() == [True]
    assert contains_points(SQUARE, np.empty((0, 2))).tolist() == []
