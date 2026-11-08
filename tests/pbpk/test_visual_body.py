"""Tests of the body of an individual drawn from the anatomy resources."""

import numpy as np
import pytest
from matplotlib.path import Path as MplPath

from pkpdsim.pbpk.visual.body import Body, Individual, body_strip, reference_body
from pkpdsim.sex import Sex


def _area(subs: list[np.ndarray]) -> float:
    return sum(
        0.5
        * abs(
            float(
                np.dot(s[:, 0], np.roll(s[:, 1], 1))
                - np.dot(s[:, 1], np.roll(s[:, 0], 1))
            )
        )
        for s in subs
    )


def test_the_reference_body_is_deterministic() -> None:
    assert Body(Individual(Sex.MALE)).svg() == reference_body("male").svg()


@pytest.mark.parametrize(
    ("individual", "height"),
    [
        (Individual(Sex.MALE), 176.0),
        (Individual(Sex.FEMALE), 163.0),
        (Individual(Sex.MALE, age=1.0), 76.0),
    ],
)
def test_the_height_of_the_body_is_the_height_of_the_individual(
    individual: Individual, height: float
) -> None:
    body = Body(individual)
    assert body.bounds()[2] == pytest.approx(height, abs=1.0)
    assert body.quantities.height == pytest.approx(height, abs=1.0)


def test_a_smaller_liver_is_drawn_smaller() -> None:
    ref = reference_body("male")
    cirrhosis = Body(Individual(Sex.MALE, liver=0.6655))
    ratio = _area(cirrhosis.organs["li"]) / _area(ref.organs["li"])
    assert ratio == pytest.approx(0.6655 ** (2 / 3), rel=0.01)


@pytest.mark.parametrize(
    "individual",
    [
        Individual(Sex.MALE),
        Individual(Sex.FEMALE),
        Individual(Sex.MALE, age=1.0),
        Individual(Sex.MALE, bw=40 * 1.76**2),
    ],
)
def test_the_organs_stay_in_the_body(individual: Individual) -> None:
    body = Body(individual)
    outline = MplPath(body.outline)
    for key, pieces in body.organs.items():
        for piece in pieces:
            assert outline.contains_points(piece).all(), key


def test_obesity_widens_the_trunk_not_the_organs() -> None:
    ref = reference_body("male")
    obese = Body(Individual(Sex.MALE, bw=40 * 1.76**2))
    y = ref.landmarks["waist"]

    def trunk_width(b: Body) -> float:
        left, right = b.edge_point(y, -1, inset=0.0), b.edge_point(y, 1, inset=0.0)
        return float(right[0] - left[0])

    assert trunk_width(obese) > 1.15 * trunk_width(ref)
    v, v_ref = obese.quantities.volumes["li"], ref.quantities.volumes["li"]
    ratio = _area(obese.organs["li"]) / _area(ref.organs["li"])
    assert ratio == pytest.approx((v / v_ref) ** (2 / 3), rel=0.05)
    assert obese.free["abdomen"] > 1.5 * ref.free["abdomen"]


def test_the_quantities_are_those_of_the_library() -> None:
    from pkpdsim.pbpk.physiology import reference_individual

    q = Individual(Sex.FEMALE).quantities()
    ref = reference_individual(Sex.FEMALE)
    assert q.bw == ref.bw and q.height == ref.height
    assert q.volumes["li"] == pytest.approx(
        ref.volumes[next(o for o in ref.volumes if o.value == "li")]
    )
    assert q.vat > 0 and q.sat > q.vat


def test_the_strip_has_a_caption_per_individual() -> None:
    svg = body_strip(
        [
            Individual(Sex.MALE, label="reference man"),
            Individual(Sex.MALE, age=85.0, label="85 years"),
        ]
    )
    assert svg.startswith("<svg") and "reference man" in svg and "85 years" in svg
