"""Tests of the figure of a model, the SVG schematic of a definition.

The figure is checked as a drawing: valid SVG, every organ and substance
labelled, every label inside the canvas and inside its box, no two labels
overlapping. A label is a box from its position, its font size and the
width of a character of the monospace fonts of the figure (`CHAR_WIDTH`).
"""

from pathlib import Path
from typing import NamedTuple
from xml.etree import ElementTree

import pytest
from pymetadata.omex import EntryFormat

from pkpdsim.pbpk import (
    MM,
    Binding,
    CustomReaction,
    Direction,
    Excretion,
    FirstOrderDissolution,
    Inhibition,
    Kp,
    LumenConversion,
    Membrane,
    Metabolism,
    Mode,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    Route,
    SegmentedIntestine,
    SlowPool,
    Substance,
    Transport,
    U,
    Variants,
    create_archive,
    model_figure,
)
from pkpdsim.pbpk.figure import CHAR_WIDTH

from .conftest import cardiac_pd, definition

SVG = "{http://www.w3.org/2000/svg}"

#: the organ names of the rows of the figure
ORGANS = (
    "lung",
    "heart",
    "brain",
    "muscle",
    "skin",
    "adipose",
    "bone",
    "gonads",
    "kidney",
    "liver",
    "spleen",
    "pancreas",
    "stomach",
)


def seven_substances() -> PBPKDefinition:
    """A definition of seven substances with every kind of process of the figure.

    A drug with three metabolites (CYP3A4 into two, CYP1A2 into one, which
    the liver excretes into the bile; P-gp efflux into the lumen), an
    inhibitor of CYP3A4, an endogenous substance with a production in the
    muscle and a slow pool, given orally and converted in the stomach into a
    seventh substance; the zonated liver, the nephron kidney and a PD
    component.
    """
    substances = [
        Substance(
            "api",
            name="apixaban",
            mr=459.5,
            binding=Binding(0.13, protein="albumin", ckd=False),
            distribution=Kp(2.0, ftissue=1.0),
            excretion=Excretion(urine=True),
        ),
        Substance(
            "api_m1",
            name="O-desmethyl apixaban",
            mr=445.5,
            excretion=Excretion(urine=True),
        ),
        Substance(
            "api_m2",
            name="hydroxy apixaban",
            mr=475.5,
            excretion=Excretion(urine=True, fraction_reabsorbed=0.2),
        ),
        Substance("api_m7", name="O-desmethyl apixaban sulfate", mr=525.6),
        Substance(
            "keto",
            name="ketoconazole",
            mr=531.4,
            binding=Binding(0.01, protein="albumin"),
            distribution=Kp(1.0, mode=Mode.PERFUSION),
            inhibits=(Inhibition("CYP3A4", "competitive", ki=1e-5),),
        ),
        Substance(
            "cr",
            name="creatine",
            mr=131.1,
            distribution=Kp(1.0, ftissue=1.0, overrides={Organ.MUSCLE: 50.0}),
            baseline=0.03,
        ),
        Substance(
            "crtn",
            name="creatinine",
            mr=113.1,
            distribution=Kp(0.8, ftissue=1.0),
            excretion=Excretion(urine=True),
        ),
    ]
    chemistry = [
        Metabolism(
            "API3A4",
            Organ.LIVER,
            "api",
            {"api_m1": 0.6, "api_m2": 0.4},
            MM(vmax=1e-3, km=0.01),
            protein="CYP3A4",
        ),
        Metabolism(
            "API1A2",
            Organ.LIVER,
            "api",
            {"api_m7": 1.0},
            MM(vmax=1e-4, km=0.01),
            protein="CYP1A2",
        ),
        Transport(
            "APIM7BILE",
            Organ.LIVER,
            "api_m7",
            Direction.EFFLUX,
            Membrane.APICAL,
            MM(vmax=1e-3, km=0.1),
            protein="ABCC2",
        ),
        Transport(
            "APIPGP",
            Organ.SMALL_INTESTINE,
            "api",
            Direction.EFFLUX,
            Membrane.APICAL,
            MM(vmax=1e-4, km=0.05),
            protein="ABCB1",
        ),
        CustomReaction(
            "CR2CRTN",
            Organ.MUSCLE,
            Part.TISSUE,
            "Cti_cr -> Cti_crtn",
            "k_cr2crtn * Vti * Cti_cr",
        ),
    ]
    return PBPKDefinition(
        sid="seven",
        name="seven substances",
        substances=substances,
        chemistry=chemistry,
        variants=Variants(
            {
                Organ.LIVER: "zonated",
                Organ.KIDNEY: "nephron",
                Organ.SMALL_INTESTINE: SegmentedIntestine(
                    conversions=[LumenConversion("cr", "crtn", 0.5, ["st"])]
                ),
            }
        ),
        routes=[
            Route.iv("api"),
            Route.oral("api", FirstOrderDissolution(k=2.0, fed=0.5)),
            Route.oral_solution("keto"),
            Route.oral_solution("cr"),
        ],
        pd=[cardiac_pd(reads="api")],
        parameters=[ModelParameter("k_cr2crtn", 1e-5, U.per_min)],
        slow_pools=[SlowPool("cr", frozenset({Organ.MUSCLE}), k=0.001, kp=2.0)],
    )


class Label(NamedTuple):
    """A text of the figure: its content and box."""

    text: str
    x0: float
    y0: float
    x1: float
    y1: float


def _label(element: ElementTree.Element) -> Label:
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


def _parse(svg: Path) -> tuple[ElementTree.Element, list[Label]]:
    """The root and the labels of a figure."""
    root = ElementTree.parse(svg).getroot()
    return root, [_label(t) for t in root.iter(f"{SVG}text")]


def _overlap(a: Label, b: Label) -> bool:
    """Whether two label boxes overlap."""
    return a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1


def _definitions() -> dict[str, PBPKDefinition]:
    """The definitions of the figure tests by id.

    The simple definition (segmented intestine, iv and oral), the seven
    substances, the single lumen with an oral route only, and the standard
    intestine without a lumen (small and large intestine as organs).
    """
    return {
        "simple": definition(),
        "seven": seven_substances(),
        "single": definition(
            sid="single",
            variants=Variants({Organ.SMALL_INTESTINE: "single_lumen"}),
            routes=[Route.oral("drg")],
        ),
        "standard": definition(
            sid="standard",
            variants=Variants({Organ.SMALL_INTESTINE: "standard"}),
            routes=[Route.iv("drg")],
        ),
    }


@pytest.fixture(params=["simple", "seven", "single", "standard"])
def figure(
    request: pytest.FixtureRequest, tmp_path: Path
) -> tuple[PBPKDefinition, Path]:
    """The figure of a definition of `_definitions`."""
    d = _definitions()[request.param]
    return d, model_figure(d, tmp_path / f"{d.sid}.svg")


def test_the_seven_substances_are_a_valid_definition() -> None:
    """The definition of the figure test validates, so it is a model the library builds."""
    seven_substances().validate()


def test_the_figure_is_valid_svg(
    figure: tuple[PBPKDefinition, Path], tmp_path: Path
) -> None:
    """An SVG of the size of its view box with a title, the same bytes for the same definition."""
    d, path = figure
    root, _ = _parse(path)
    assert root.tag == f"{SVG}svg"
    width, height = float(root.attrib["width"]), float(root.attrib["height"])
    assert (
        root.attrib["viewBox"] == f"0 0 {root.attrib['width']} {root.attrib['height']}"
    )
    assert width > 0 and height > 0
    assert root.find(f"{SVG}title") is not None
    again = model_figure(d, tmp_path / "again.svg")
    assert again.read_bytes() == path.read_bytes()


def test_every_organ_and_substance_is_labelled(
    figure: tuple[PBPKDefinition, Path],
) -> None:
    """Every organ has a row with its name, every substance a chip in every organ row."""
    d, path = figure
    _, labels = _parse(path)
    texts = [label.text for label in labels]
    for organ in ORGANS:
        assert any(t.startswith(f"{organ} ") for t in texts), organ
    if d.sid == "standard":
        intestine = ("small intestine si", "large intestine la")
    else:
        intestine = ("intestine si la",)
    for title in intestine:
        assert title in texts, title
    rows = 10 + 2 + len(intestine)  # lung, 8 systemic, liver, spleen, pancreas, ...
    for s in d.substances:
        # the chip of the legend and one chip per organ row
        assert texts.count(s.sid) >= 1 + rows, s.sid


def test_the_labels_are_inside_the_canvas_and_do_not_overlap(
    figure: tuple[PBPKDefinition, Path],
) -> None:
    """Every label box lies in the canvas and no two label boxes overlap."""
    _, path = figure
    root, labels = _parse(path)
    width, height = float(root.attrib["width"]), float(root.attrib["height"])
    for label in labels:
        assert label.x0 >= 0 and label.y0 >= 0, label
        assert label.x1 <= width and label.y1 <= height, label
    for i, a in enumerate(labels):
        for b in labels[i + 1 :]:
            assert not _overlap(a, b), (a, b)


def test_the_labels_of_a_box_are_inside_it(figure: tuple[PBPKDefinition, Path]) -> None:
    """The labels of an organ, a blood column or the lumen lie inside its rectangle."""
    _, path = figure
    root, _ = _parse(path)
    boxes = [g for g in root.iter(f"{SVG}g") if g.attrib.get("class") == "pkf-box"]
    assert len(boxes) > 15
    for g in boxes:
        rect = g.find(f"{SVG}rect")
        assert rect is not None
        x0, y0 = float(rect.attrib["x"]), float(rect.attrib["y"])
        x1, y1 = x0 + float(rect.attrib["width"]), y0 + float(rect.attrib["height"])
        for t in g.iter(f"{SVG}text"):
            label = _label(t)
            assert x0 < label.x0 and label.x1 < x1, label
            assert y0 < label.y0 and label.y1 < y1, label


def test_the_figure_shows_the_processes_and_routes(tmp_path: Path) -> None:
    """The processes, the variants, the routes and the lumen of the seven substances."""
    path = model_figure(seven_substances(), tmp_path / "seven.svg")
    text = " ".join(label.text for label in _parse(path)[1])
    for expected in (
        "API3A4: api -> api_m1, api_m2 by CYP3A4",
        "APIM7BILE: api_m7 cell -> bile by ABCC2",
        "APIPGP: api cell -> lumen by ABCB1",
        "CR2CRTN: Cti_cr -> Cti_crtn",
        "slow pool of cr",
        "api_m2 (reabsorbed)",
        "zonated, 8 zones",
        "nephron",
        "12 segments",
        "cr -> crtn in st",
        "bile into du: api_m7",
        "iv: api",
        "inhibits CYP3A4 (competitive)",
        "PD cardio",
    ):
        assert expected in text, expected


def test_the_archive_carries_the_figure(tmp_path: Path) -> None:
    """`create_archive` writes the figure next to the models and adds it; `figure=False` leaves it out."""
    readme = tmp_path / "README.md"
    readme.write_text("# simple\n", encoding="utf-8")
    results = create_archive(
        definition(),
        tmp_path / "models",
        model="simple",
        copyright="2024-2026 Matthias König",
        readme=readme,
    )
    entry = results["simple_svg"]["entry"]
    assert entry.location == "./models/simple.svg"
    assert entry.format == EntryFormat.SVG
    assert results["simple_svg"]["path"] == tmp_path / "models" / "simple.svg"
    plain = tmp_path / "plain"
    results = create_archive(
        definition(),
        plain / "models",
        model="simple",
        copyright="2024-2026 Matthias König",
        figure=False,
    )
    assert "simple_svg" not in results
    assert not (plain / "models" / "simple.svg").exists()


def test_without_an_oral_route_the_lumen_absorbs_nothing(tmp_path: Path) -> None:
    """The lumen of a model without an oral route has no absorption and no oral arrow.

    A substance without an oral route has `ka_<s> = 0`; the lumen only
    carries the bile to the feces.
    """
    path = model_figure(definition(routes=[Route.iv("drg")]), tmp_path / "iv.svg")
    texts = [label.text for label in _parse(path)[1]]
    assert "lumen" in texts
    assert "absorption" not in texts
    assert "oral" not in texts
    path = model_figure(definition(), tmp_path / "oral.svg")
    texts = [label.text for label in _parse(path)[1]]
    assert {"absorption", "oral"} <= set(texts)
