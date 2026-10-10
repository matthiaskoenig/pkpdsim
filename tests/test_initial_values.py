from pathlib import Path

import pytest
import roadrunner
from sbmlutils.factory import Compartment, Model, Species, create_model

from pkpdsim.initial_values import (
    _set_initial_values,
    get_species_keys,
    init_amounts_changes,
    init_concentrations_changes,
)


def test_get_species_keys() -> None:
    """Species of a substance are found by the naming pattern of the ids."""
    species_ids = ["Cve_x", "Cli_plasma_x", "Aurine_x", "Cve_xy", "Cve_x_bound", "x"]
    assert get_species_keys("x", species_ids) == ["Cve_x", "Cli_plasma_x", "Aurine_x"]


def test_get_species_keys_no_match() -> None:
    """An unknown substance has no species."""
    assert get_species_keys("y", ["Cve_x"]) == []


@pytest.fixture(scope="module")
def r(tmp_path_factory: pytest.TempPathFactory) -> roadrunner.RoadRunner:
    """A model with a substance `x` in the venous blood, the liver and the urine."""
    model = Model(
        sid="distribution",
        compartments=[
            Compartment("Vve", 5.0),
            Compartment("Vli", 1.5),
            Compartment("Vurine", 1.0),
        ],
        species=[
            Species("Cve_x", "Vve", initialConcentration=0.0),
            Species("Cli_plasma_x", "Vli", initialConcentration=0.0),
            Species(
                "Aurine_x", "Vurine", initialAmount=0.0, hasOnlySubstanceUnits=True
            ),
            Species("Cve_y", "Vve", initialConcentration=0.0),
            Species("Cli_x", "Vli", initialConcentration=0.0, boundaryCondition=True),
        ],
    )
    path: Path = tmp_path_factory.mktemp("model") / "distribution.xml"
    create_model(model, filepath=path, validate=False)
    return roadrunner.RoadRunner(str(path))


def test_init_concentrations_changes(r: roadrunner.RoadRunner) -> None:
    """The concentration of the substance is set in every tissue but the urine."""
    changes = init_concentrations_changes(r, "x", 2.0)
    assert changes == {"[Cve_x]": 2.0, "[Cli_plasma_x]": 2.0, "[Cli_x]": 2.0}


def test_init_amounts_changes(r: roadrunner.RoadRunner) -> None:
    """The amount of the substance is set in every tissue but the urine."""
    changes = init_amounts_changes(r, "x", 3.0)
    assert changes == {"Cve_x": 3.0, "Cli_plasma_x": 3.0, "Cli_x": 3.0}


def test_changes_apply(r: roadrunner.RoadRunner) -> None:
    """The changes are selections of roadrunner, which sets them on the model."""
    r.resetAll()
    for key, value in init_concentrations_changes(r, "x", 2.0).items():
        r.setValue(f"init({key})", value)
    r.reset()
    assert r.getValue("[Cve_x]") == pytest.approx(2.0)
    assert r.getValue("[Cli_plasma_x]") == pytest.approx(2.0)
    assert r.getValue("Aurine_x") == pytest.approx(0.0)
    assert r.getValue("[Cve_y]") == pytest.approx(0.0)


def test_unknown_method(r: roadrunner.RoadRunner) -> None:
    """Only concentrations and amounts are initial values."""
    with pytest.raises(ValueError, match="method must be"):
        _set_initial_values(r, "x", 1.0, method="volume")
