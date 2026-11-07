"""A test model package with every kind of evidence of the credibility declaration.

It has the SBML model of the toy (`pkpdsim_toy`), built by its factory into the
folder `models` of this package by `tests/conftest.py`, and has two experiments with mappings which are calibration, dose, population and
limitation evidence, unclaimed validation data and an outlier. The renal
impairment is a label of the metadata only, the model does not change. The
decay does not start at its steady state, so the battery does not check it
(`battery_settings`).
"""

from pathlib import Path

from sbmlsim.project import SimulatorSettings

from pkpdsim.model import CREDIBILITY_SETTINGS, FigureSettings, PKPDModel
from pkpdsim_toy import DATA_PATHS

CREDIBILITY_PATH = Path(__file__).parent
MODEL_BASE_PATH = CREDIBILITY_PATH / "models"
MODEL_PATH = MODEL_BASE_PATH / "toy.xml"

MODEL = PKPDModel(
    name="credibility",
    package=__name__,
    title="Toy model with every kind of evidence",
    path=CREDIBILITY_PATH,
    model_path=MODEL_PATH,
    factory="pkpdsim_toy.factory:create_toy_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting:FIT_DEFINITIONS",
    credibility=".credibility:CREDIBILITY",
    credibility_settings=CREDIBILITY_SETTINGS.model_copy(
        update={"battery_settings": {"steady_state": False}}
    ),
    data_paths=tuple(DATA_PATHS),
    simulator=SimulatorSettings(absolute_tolerance=1e-10, relative_tolerance=1e-10),
    figures=FigureSettings(formats=("svg",), savefig_dpi=30),
)
