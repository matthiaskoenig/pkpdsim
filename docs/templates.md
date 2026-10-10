# Templates and annotations

Every model of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models) names its creators, carries the terms of use under which it is published and annotates its compartments with the terms of ontologies, so that a reader of the SBML, a person or a program, knows who made it, how it may be used and what a compartment is. `pkpdsim.templates` holds the creators and the terms of use, `pkpdsim.annotations` the annotations of the compartments which recur in the models: the liver `li`, the kidneys `ki`, the gut `gu`, the venous and arterial blood `ve` and `ar`, the urine, the feces and more.

## Usage

A model of the liver with its creators, the terms of use behind its description and an annotated compartment, created with the factory of [sbmlutils](https://github.com/matthiaskoenig/sbmlutils):

```python
from pathlib import Path

import libsbml
from sbmlutils.factory import Compartment, Model, ModelUnits, Units, create_model

from pkpdsim.annotations import compartments
from pkpdsim.templates import creators, terms_of_use


class U(Units):
    """The units of the model, `liter` and `second` are predefined."""


model = Model(
    sid="liver",
    name="Model of the liver",
    notes="A model of the liver with a single compartment.\n" + terms_of_use,
    creators=creators,
    units=U,
    model_units=ModelUnits(time=U.second, volume=U.liter),
    compartments=[
        Compartment(
            "Vli", 1.5, unit=U.liter, name="liver", annotations=compartments["li"]
        ),
    ],
)
result = create_model(model, filepath=Path("liver.xml"), sbml_level=3, sbml_version=2)

liver = libsbml.readSBMLFromFile(str(result.sbml_path)).getModel()
history = liver.getModelHistory()
print(history.getCreator(0).getGivenName(), history.getCreator(0).getFamilyName())
compartment = liver.getCompartment("Vli")
for k in range(compartment.getNumCVTerms()):
    term = compartment.getCVTerm(k)
    for m in range(term.getNumResources()):
        print(term.getResourceURI(m))
```

```text
Matthias König
https://identifiers.org/FMA:7197
https://identifiers.org/BTO:0000759
https://identifiers.org/ncit:C12392
```

## Creators

`creators` is a list of `sbmlutils.factory.Creator`; sbmlutils writes it as the model history of the SBML. A model with further creators, e.g. the students who built it, extends the list:

```python
# not executed: the creators of a model with a further creator
from sbmlutils.factory import Creator

from pkpdsim.templates import creators

model_creators = [
    *creators,
    Creator(familyName="Doe", givenName="Jane", email="jane.doe@example.org"),
]
```

## Terms of use

`terms_of_use` is the markdown of the terms of use, which the models append to the notes behind their description: the copyright of the current year, the [Creative Commons Attribution 4.0](http://creativecommons.org/licenses/by/4.0/) license of the model and the disclaimer. The year is the one of the import, so a model created in a new year carries the new year.

## Annotations

`compartments` maps the abbreviation of a compartment to its annotation, a list of `(qualifier, resource)` pairs: the qualifier is a biological qualifier of `sbmlutils.metadata.BQB` (`IS`, `IS_PART_OF`, `IS_VERSION_OF`) and the resource the term of an ontology as `collection/term` of [identifiers.org](https://identifiers.org), e.g. `fma/FMA:7197` for the liver in the Foundational Model of Anatomy. A model passes the list as the `annotations` of its compartment, or of a species which lives in the compartment, and sbmlutils resolves the terms and writes the RDF of the element.
