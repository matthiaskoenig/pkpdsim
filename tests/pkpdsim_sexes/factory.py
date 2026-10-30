"""The factory of the toy model per sex: the male and the female build."""

from pathlib import Path

from sbmlutils.factory import Model, create_model

from pkpdsim.sex import Sex
from pkpdsim_sexes import FEMALE_A0
from pkpdsim_toy.factory import toy_model


def sexes_model(sex: Sex) -> Model:
    """Get the build of a sex: the toy model, A starts lower in women."""
    model = toy_model()
    model.sid = "sexes" if sex is Sex.MALE else f"sexes_{sex}"
    model.name = f"toy decay model ({sex})"
    if sex is Sex.FEMALE:
        (species,) = model.species
        species.initialConcentration = FEMALE_A0
    return model


def create_sexes_models(output_dir: Path) -> dict[Sex, Path]:
    """Create the SBML of both builds in a directory.

    Args:
        output_dir: directory of the models, created if it does not exist.

    Returns:
        The SBML file of every sex.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[Sex, Path] = {}
    for sex in Sex:
        model = sexes_model(sex)
        result = create_model(
            model=model,
            filepath=output_dir / f"{model.sid}.xml",
            sbml_level=3,
            sbml_version=2,
        )
        paths[sex] = result.sbml_path
    return paths
