"""Tests of the anatomical figure of a model built on a laboratory animal."""

import dataclasses
from pathlib import Path
from xml.etree import ElementTree

import pytest

from pkpdsim.pbpk import Species
from pkpdsim.pbpk.visual.anatomy import anatomy_figure
from pkpdsim.pbpk.visual.animal import AnimalBody, reference_volumes
from pkpdsim.pbpk.visual.interactive import interactive_figure
from pkpdsim.pbpk.visual.publication import publication_figure

from .test_figure import seven_substances

SVG = "{http://www.w3.org/2000/svg}"


def _texts(path: Path) -> list[str]:
    root = ElementTree.parse(path).getroot()
    return ["".join(t.itertext()) for t in root.iter(f"{SVG}text")]


@pytest.mark.parametrize(
    ("species", "grams"), [(Species.MOUSE, "25"), (Species.RAT, "250")]
)
def test_an_animal_is_drawn_with_its_reference_weight(
    species: Species, grams: str, tmp_path: Path
) -> None:
    d = dataclasses.replace(seven_substances(), species=species)
    texts = _texts(anatomy_figure(d, tmp_path / "a.svg"))
    assert f"{species.value} male, {grams} g" in texts
    assert "male build" not in texts and "female build" not in texts
    assert "tail vein: the observer Carm of every substance" in " ".join(texts)


def test_the_mouse_drawing_has_every_organ_of_the_library() -> None:
    body = AnimalBody("rat")
    for key in (
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
        "mu",
        "bo",
        "ad",
    ):
        assert body.drawn(key), key


def test_a_species_without_a_drawing_is_refused() -> None:
    with pytest.raises(ValueError, match="No drawing of the pig"):
        AnimalBody("pig")


def test_the_cards_of_an_animal_carry_its_organ_masses(tmp_path: Path) -> None:
    d = dataclasses.replace(seven_substances(), species=Species.RAT)
    html = interactive_figure(d, tmp_path / "a.html").read_text()
    liver = reference_volumes("rat")["li"]
    assert f"{liver:.3g} (reference rat)" in html


def test_the_publication_figure_of_an_animal(tmp_path: Path) -> None:
    d = dataclasses.replace(seven_substances(), species=Species.MOUSE)
    path = publication_figure(d, tmp_path / "f.svg")
    assert "mouse male, 25 g" in path.with_name("f_panels.svg").read_text()
