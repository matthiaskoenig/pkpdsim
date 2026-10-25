"""Tests of the shared helpers on a tiny hand-built model."""

import pytest
import roadrunner

from .conftest import total_amount

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
