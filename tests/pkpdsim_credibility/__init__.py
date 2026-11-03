"""A test model package with every kind of evidence of the credibility declaration.

It reuses the SBML model of the toy (`pkpdsim_toy`, built by `tests/conftest.py`)
and has two experiments with mappings which are calibration, dose, population and
limitation evidence, unclaimed validation data and an outlier. The renal
impairment is a label of the metadata only, the model does not change.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings
from pkpdsim_toy import DATA_PATHS, MODEL_PATH

CREDIBILITY_PATH = Path(__file__).parent

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
    data_paths=tuple(DATA_PATHS),
    simulator=SimulatorSettings(absolute_tolerance=1e-10, relative_tolerance=1e-10),
    figures=FigureSettings(formats=("svg",), savefig_dpi=30),
)
