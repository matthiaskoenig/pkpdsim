"""The templates and the annotations build a valid SBML model with sbmlutils."""

from datetime import datetime

import libsbml
import pytest
from sbmlutils.factory import (
    Compartment,
    Model,
    ModelUnits,
    UnitDefinition,
    Units,
    create_model,
)
from sbmlutils.metadata import BQB
from sbmlutils.validation import ValidationOptions, validate_doc

from pkpdsim.annotations import compartments
from pkpdsim.templates import creators, terms_of_use


def test_creators() -> None:
    assert [creator.familyName for creator in creators] == ["König"]
    assert creators[0].site == "https://livermetabolism.com"


def test_terms_of_use() -> None:
    assert "## Terms of use" in terms_of_use
    assert f"Copyright © {datetime.now().year}" in terms_of_use
    assert "Matthias König" in terms_of_use
    assert "{" not in terms_of_use


def test_annotations() -> None:
    """Every annotation is a list of qualifier and resource pairs."""
    assert {"li", "ki", "gu", "ve", "ar", "po", "plasma", "urine"} <= set(compartments)
    for annotations in compartments.values():
        for qualifier, resource in annotations:
            assert isinstance(qualifier, BQB)
            assert isinstance(resource, str) and resource


class U(Units):
    """The units of the model."""

    mmole = UnitDefinition("mmole")
    min = UnitDefinition("min")
    m2 = UnitDefinition("m2", "meter^2")


@pytest.fixture(scope="module")
def document(tmp_path_factory: pytest.TempPathFactory) -> libsbml.SBMLDocument:
    """A model with a compartment for every annotation, the creators and the terms."""
    model = Model(
        sid="annotated",
        notes="A compartment for every annotation of pkpdsim.\n" + terms_of_use,
        creators=creators,
        units=U,
        model_units=ModelUnits(
            time=U.min,
            extent=U.mmole,
            substance=U.mmole,
            length=U.meter,
            area=U.m2,
            volume=U.liter,
        ),
        compartments=[
            Compartment(f"V{k}", 1.0, name=key, unit=U.liter, annotations=annotations)
            for k, (key, annotations) in enumerate(compartments.items())
        ],
    )
    path = tmp_path_factory.mktemp("sbml") / "annotated.xml"
    result = create_model(model, filepath=path, sbml_level=3, sbml_version=2)
    assert result.sbml_path == path
    return libsbml.readSBMLFromFile(str(path))


def test_model_is_valid(document: libsbml.SBMLDocument) -> None:
    validation = validate_doc(document, options=ValidationOptions())
    assert validation.error_count == 0


def test_model_has_creators_and_terms(document: libsbml.SBMLDocument) -> None:
    model = document.getModel()
    history = model.getModelHistory()
    assert history.getNumCreators() == 1
    assert history.getCreator(0).getFamilyName() == "König"
    assert "Terms of use" in model.getNotesString()


def test_compartments_are_annotated(document: libsbml.SBMLDocument) -> None:
    model = document.getModel()
    for k, (key, annotations) in enumerate(compartments.items()):
        compartment = model.getCompartment(f"V{k}")
        assert compartment.getName() == key
        resources = [
            compartment.getCVTerm(n).getResourceURI(m)
            for n in range(compartment.getNumCVTerms())
            for m in range(compartment.getCVTerm(n).getNumResources())
        ]
        assert len(resources) == len(annotations), key
