"""The panels of the body figure of a model, restricted to its build."""

from pathlib import Path

import pytest
from sbmlsim.panel import Knob, PanelSimulation

from pkpdsim.pbpk import build
from pkpdsim.pbpk.visual.anatomy import anatomy_figure
from pkpdsim.pbpk.visual.content import ModelContent
from pkpdsim.pbpk.visual.interactive import write_readouts
from pkpdsim.webapp import body_figure

from .conftest import definition

SIMULATION = PanelSimulation(end=24.0, points=200, x_unit="hr")
KNOBS = (
    Knob(
        id="BW",
        label="body weight",
        target="BW",
        unit="kg",
        kind="range",
        min=40.0,
        max=170.0,
        default=75.0,
    ),
)


def _figure(tmp_path: Path) -> tuple[Path, Path]:
    """The drawing with its sidecar and the flat build of the fixture.

    The drug and its one metabolite, since a panel of sbmlsim draws at most
    16 curves of one unit and the figure of a model of seven substances has
    18 in the liver.
    """
    model = definition()
    svg = anatomy_figure(model, tmp_path / "m_body.svg")
    write_readouts(ModelContent.of(model), svg)
    return svg, build(model, tmp_path).flat_path


def test_every_panel_selects_an_entity_of_the_build(tmp_path: Path) -> None:
    """A key whose ids are in no build is dropped instead of being declared."""
    svg, sbml = _figure(tmp_path)
    panels, figure, dropped = body_figure(
        svg, sbml, simulation=SIMULATION, title="fixture"
    )
    assert panels
    assert figure.title == "fixture"
    entities = set()
    doc = sbml.read_text(encoding="utf-8")
    for panel in panels.values():
        for curve in panel.curves:
            sid = curve.id[1:-1] if curve.id.startswith("[") else curve.id
            assert f'id="{sid}"' in doc, f"{panel.id} selects {curve.id}"
            entities.add(sid)
    assert entities
    assert all(key not in panels for key in dropped)


def test_the_knobs_of_the_app_are_on_every_panel(tmp_path: Path) -> None:
    """A derived target of the simulation reads a knob of its own panel."""
    svg, sbml = _figure(tmp_path)
    panels, _, _ = body_figure(
        svg, sbml, simulation=SIMULATION, title="fixture", knobs=KNOBS
    )
    assert panels
    assert {tuple(p.knobs) for p in panels.values()} == {KNOBS}


def test_a_curve_is_labelled_as_a_panel_of_an_app_labels_it(tmp_path: Path) -> None:
    """The compartment alone, the substance beside it where a panel draws several."""
    svg, sbml = _figure(tmp_path)
    panels, _, _ = body_figure(
        svg,
        sbml,
        simulation=SIMULATION,
        title="fixture",
        substances={"metabolite": "M"},
    )
    labels = {c.label for p in panels.values() for c in p.curves}
    assert labels
    assert not any("," in label for label in labels), labels
    assert not any(label.endswith(" unbound") for label in labels), labels
    assert max(len(label) for label in labels) <= 30, labels


def test_a_path_which_does_not_exist_is_named(tmp_path: Path) -> None:
    """A missing drawing, sidecar or build raises instead of being skipped."""
    svg, sbml = _figure(tmp_path)
    with pytest.raises(FileNotFoundError, match="drawing"):
        body_figure(tmp_path / "no.svg", sbml, simulation=SIMULATION, title="x")
    with pytest.raises(FileNotFoundError, match="build"):
        body_figure(svg, tmp_path / "no.xml", simulation=SIMULATION, title="x")
    with pytest.raises(FileNotFoundError, match="readouts"):
        body_figure(
            svg,
            sbml,
            simulation=SIMULATION,
            title="x",
            readouts=tmp_path / "no.json",
        )
