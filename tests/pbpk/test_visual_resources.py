"""Tests of the anatomy resources of the figures."""

import re

import pytest
from matplotlib.path import Path as MplPath

from pkpdsim.pbpk.visual.resources import anatomy, parse_d

ORGANS = (
    "br",
    "lu",
    "he",
    "li",
    "ki",
    "st",
    "si",
    "la",
    "sp",
    "pa",
    "go",
    "mu",
    "bo",
    "ad",
)


@pytest.mark.parametrize("sex", ["male", "female"])
def test_every_human_resource_has_the_organs_of_the_library(sex: str) -> None:
    a = anatomy("human", sex)
    assert set(ORGANS) <= set(a.organs)
    assert {"bladder", "gallbladder", "aorta"} <= set(a.anatomy)


@pytest.mark.parametrize(
    ("species", "sex"),
    [
        ("human", "male"),
        ("human", "female"),
        ("mouse", "male"),
        ("mouse", "female"),
        ("rat", "male"),
    ],
)
def test_the_uberon_ids_are_terms_of_uberon(species: str, sex: str) -> None:
    a = anatomy(species, sex)
    for shape in [*a.organs.values(), *a.anatomy.values()]:
        assert shape.uberon, shape.key
        for u in shape.uberon:
            assert re.fullmatch(r"UBERON_\d{7,8}", u), u


@pytest.mark.parametrize("sex", ["male", "female"])
def test_the_cavities_lie_inside_the_outline(sex: str) -> None:
    a = anatomy("human", sex)
    body = MplPath(a.outline)
    for cavity in ("thorax", "abdomen"):
        for poly in a.cavities[cavity]:
            assert body.contains_points(poly).all(), cavity


def test_the_attribution_names_the_source_and_the_license() -> None:
    a = anatomy()
    assert "CC BY 4.0" in a.attribution
    assert "ebi-gene-expression-group/anatomogram" in a.attribution


@pytest.mark.parametrize("sex", ["male", "female"])
def test_the_landmarks_are_ordered(sex: str) -> None:
    m = anatomy("human", sex).landmarks
    order = ["top", "neck", "shoulder", "armpit", "waist", "crotch", "bottom"]
    assert [m[k] for k in order] == sorted(m[k] for k in order)


def test_a_species_without_a_resource_is_refused() -> None:
    with pytest.raises(ValueError, match="No anatomy for dog male"):
        anatomy("dog", "male")


def test_parse_d_reads_subpaths() -> None:
    subs = parse_d("M0,0 1,0 1,1Z M2,2 3,2 3,3Z")
    assert len(subs) == 2
    assert subs[1].tolist() == [[2.0, 2.0], [3.0, 2.0], [3.0, 3.0]]
