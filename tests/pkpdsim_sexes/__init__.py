"""A toy model package with a build per sex: the decay of A in men and in women.

The male build is the toy model (`pkpdsim_toy`), the female build starts at a
lower concentration of A (`FEMALE_A0`); both decay with the rate `k`, the
parameter of the fit, so one value of `k` fits the data of both sexes. The
tests of pkpdsim run the per-sex builds through the commands, the fits, the
draws of the fitted parameters, the sensitivity analyses, the verification
and the test functions with it. The SBML models are created by
`factory.create_sexes_models` before the tests start (`tests/conftest.py`) and
are not tracked.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings

SEXES_PATH = Path(__file__).parent
MODEL_BASE_PATH = SEXES_PATH / "models"
MODEL_PATH = MODEL_BASE_PATH / "sexes.xml"
FEMALE_MODEL_PATH = MODEL_BASE_PATH / "sexes_female.xml"

#: the initial concentration of A in the female build [mM], 1 in the male one
FEMALE_A0 = 0.8

MODEL = PKPDModel(
    name="sexes",
    package=__name__,
    title="Toy PK model per sex",
    path=SEXES_PATH,
    model_path=MODEL_PATH,
    female_model_path=FEMALE_MODEL_PATH,
    factory=".factory:create_sexes_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting:FIT_DEFINITIONS",
    simulator=SimulatorSettings(absolute_tolerance=1e-10, relative_tolerance=1e-10),
    figures=FigureSettings(formats=("svg",)),
)
