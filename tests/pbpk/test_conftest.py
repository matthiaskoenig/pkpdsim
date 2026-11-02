"""Tests of the shared helpers on a tiny hand-built model."""

from pathlib import Path
from typing import TYPE_CHECKING

import pytest
import roadrunner
from sbmlutils.factory import Creator

from pkpdsim.pbpk.builder import build

from .conftest import dateless, definition, total_amount

if TYPE_CHECKING:
    from conftest import _FlatModelCache

_SBML = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="m">
    <listOfCompartments>
      <compartment id="Vli" size="2" constant="true" spatialDimensions="3"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="Cpl_drg" compartment="Vli" initialConcentration="3"
        hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="Aurine_drg" compartment="Vli" initialAmount="5"
        hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="Cpl_other" compartment="Vli" initialConcentration="100"
        hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
  </model>
</sbml>
"""


def test_total_amount() -> None:
    """A concentration is multiplied with the volume, an amount is taken as is."""
    r = roadrunner.RoadRunner(_SBML)
    assert total_amount(r, "drg") == pytest.approx(3 * 2 + 5)


def test_a_cached_flat_model_is_the_flattened_one(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    flat_model_cache: "_FlatModelCache",
) -> None:
    """A build of a comp model which was flattened before takes its flat model.

    The definition has creators, so the comp model has the dates of the model
    history, which the key of the cache leaves out; an id of its own keeps it
    apart from the builds of the other tests.
    """
    d = definition(
        sid="cached",
        creators=[Creator("Doe", "Jane", "jane@example.org", "Example")],
    )
    flattened = []
    flat_model = flat_model_cache.flat_model

    def counted(comp_path: Path, flat_path: Path) -> None:
        flattened.append(flat_path)
        flat_model(comp_path, flat_path)

    monkeypatch.setattr(flat_model_cache, "flat_model", counted)
    results = []
    for name in ("first", "second"):
        folder = tmp_path / name
        folder.mkdir()
        results.append(build(d, folder))
    first, second = (r.flat_path.read_text() for r in results)
    assert "<dcterms:W3CDTF>" in first
    assert len(flattened) == 1
    assert dateless(second) == dateless(first)
