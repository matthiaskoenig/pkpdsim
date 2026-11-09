"""A toy cell model, the fixture of the composition tests.

It carries the ports a cell model of the zonation has: the compartments `Vext`
(space of Disse), `Vli` (hepatocyte), `Vbi` (bile), the concentration `drg_ext` in `Vext`,
the lobular `position` and the `temperature`, all with `port=True`. One
reaction takes up `drg_ext` into the hepatocyte with a first-order rate, so
that the flat model has something to integrate. It imports nothing of
`pkpdsim.zonation`.
"""

from pathlib import Path

from sbmlutils.factory import (
    Compartment,
    Model,
    Parameter,
    Reaction,
    Species,
    create_model,
)
from sbmlutils.metadata import BQB

from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.perfusion.units import PU

TOY_ID = "toy"
"""The id of the toy cell model, and the stem of its file."""

UPTAKE = 0.5
"""The first-order uptake rate (1/min) of the toy."""


def create_toy_cell(output_dir: Path, annotated: bool = False) -> Path:
    """Write the toy cell model as `toy.xml` into `output_dir`.

    Args:
        output_dir: the folder.
        annotated: annotate the external drug with a term, which gives the
            species a metaid as the species of every real cell model have.

    Returns:
        The path of the file.
    """
    model = Model(
        sid=TOY_ID,
        name="toy cell model",
        units=PU,
        model_units=MODEL_UNITS,
        compartments=[
            Compartment("Vext", 1.0, unit=U.litre, name="external", port=True),
            Compartment("Vli", 1.0, unit=U.litre, name="cell", port=True),
            Compartment("Vbi", 1.0, unit=U.litre, name="bile", port=True),
        ],
        species=[
            Species(
                "drg_ext",
                "Vext",
                initialConcentration=0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,
                name="drug external",
                annotations=[(BQB.IS, "chebi/CHEBI:17234")] if annotated else [],
                port=True,
            ),
            Species(
                "drg",
                "Vli",
                initialConcentration=0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,
                name="drug in the cell",
            ),
        ],
        parameters=[
            Parameter(
                "position",
                0.5,
                unit=U.dimensionless,
                name="lobular position",
                constant=True,
                port=True,
            ),
            Parameter(
                "temperature",
                37.0,
                unit=PU.degC,
                name="temperature",
                constant=True,
                port=True,
            ),
            Parameter(
                "k_up", UPTAKE, unit=U.per_min, name="uptake rate", constant=True
            ),
        ],
        reactions=[
            Reaction(
                "UP",
                "drg_ext -> drg",
                formula=("k_up * drg_ext * Vext", U.mmole_per_min),
                name="uptake",
                port=True,
            )
        ],
    )
    path = output_dir / f"{TOY_ID}.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    return path
