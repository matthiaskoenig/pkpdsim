"""The toy model package of the tests: its factory creates a valid model."""

import math
from pathlib import Path

import libsbml
import pytest
import roadrunner

from pkpdsim_toy import MODEL_PATH
from pkpdsim_toy.factory import create_toy_models


def test_the_factory_creates_a_valid_model(tmp_path: Path) -> None:
    path = create_toy_models(tmp_path)
    assert path == tmp_path / "toy.xml"
    document = libsbml.readSBMLFromFile(str(path))
    assert document.getNumErrors(libsbml.LIBSBML_SEV_ERROR) == 0
    assert document.getModel().getId() == "toy"


def test_the_model_of_the_tests_exists() -> None:
    """The controller of the session creates the model before the tests start."""
    assert MODEL_PATH.exists()


def test_a_decays_with_the_rate_k() -> None:
    r = roadrunner.RoadRunner(str(MODEL_PATH))
    result = r.simulate(0, 5, 6, ["time", "[A]"])
    assert result[-1, 1] == pytest.approx(math.exp(-0.2 * 5), rel=1e-5)
