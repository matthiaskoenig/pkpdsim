"""The factory of the toy model: A decays with the rate k."""

from pathlib import Path

from sbmlutils.factory import (
    Compartment,
    Model,
    ModelUnits,
    Parameter,
    Reaction,
    Species,
    UnitDefinition,
    Units,
    create_model,
)


class U(Units):
    """The units of the toy model."""

    mmole = UnitDefinition("mmole")
    min = UnitDefinition("min")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    per_min = UnitDefinition("per_min", "1/min")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")


def toy_model() -> Model:
    """Get the toy model: A in the compartment Vc decays with the rate k."""
    return Model(
        sid="toy",
        name="toy decay model",
        units=U,
        model_units=ModelUnits(
            time=U.min,
            extent=U.mmole,
            substance=U.mmole,
            length=U.meter,
            area=U.m2,
            volume=U.liter,
        ),
        compartments=[Compartment("Vc", 1.0, name="central", unit=U.liter)],
        species=[
            Species(
                "A",
                initialConcentration=1.0,
                compartment="Vc",
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,
                name="A",
            )
        ],
        parameters=[Parameter("k", 0.2, unit=U.per_min, name="rate of the decay")],
        reactions=[
            Reaction(
                "DECAY",
                equation="A ->",
                formula=("k * A * Vc", U.mmole_per_min),
                compartment="Vc",
            )
        ],
    )


def create_toy_models(output_dir: Path) -> Path:
    """Create the SBML of the toy model in a directory.

    Args:
        output_dir: directory of the model, created if it does not exist.

    Returns:
        The path of the SBML file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    result = create_model(
        model=toy_model(),
        filepath=output_dir / "toy.xml",
        sbml_level=3,
        sbml_version=2,
    )
    return result.sbml_path
