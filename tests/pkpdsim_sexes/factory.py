"""The factory of the toy model per sex: the male and the female build."""

from pathlib import Path

from sbmlutils.factory import InitialAssignment, Model, Parameter, create_model

from pkpdsim.sex import Sex
from pkpdsim_sexes import FEMALE_A0
from pkpdsim_toy.factory import U, toy_model


def sexes_model(sex: Sex) -> Model:
    """Get the build of a sex: the toy model, A starts at `A0`, lower in women.

    `A0` is a parameter of the build whose value differs by sex, as the
    reference values of a library model (the body weight, the hematocrit) do.
    """
    model = toy_model()
    model.sid = "sexes" if sex is Sex.MALE else f"sexes_{sex}"
    model.name = f"toy decay model ({sex})"
    a0 = FEMALE_A0 if sex is Sex.FEMALE else 1.0
    model.parameters = [
        *model.parameters,
        Parameter("A0", a0, unit=U.mM, name="initial concentration of A"),
    ]
    model.assignments = [InitialAssignment("A", "A0", unit=U.mM)]
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
