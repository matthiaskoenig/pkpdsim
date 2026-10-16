"""A toy model package: the decay of a substance A in one compartment.

The tests of pkpdsim run the model object, the commands, the verification and the
test functions on it. The SBML model is created by `factory.create_toy_models`
before the tests start (`tests/conftest.py`) and is not tracked.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings
from pkpdsim.release.settings import Author, ReleaseSettings

TOY_PATH = Path(__file__).parent
MODEL_BASE_PATH = TOY_PATH / "models"
MODEL_PATH = MODEL_BASE_PATH / "toy.xml"
DATA_PATHS = [TOY_PATH]

MODEL = PKPDModel(
    name="toy",
    package=__name__,
    title="Toy PK model",
    path=TOY_PATH,
    model_path=MODEL_PATH,
    factory=".factory:create_toy_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting:FIT_DEFINITIONS",
    data_paths=tuple(DATA_PATHS),
    simulator=SimulatorSettings(absolute_tolerance=1e-10, relative_tolerance=1e-10),
    figures=FigureSettings(formats=("svg", "png"), savefig_dpi=30, legend_fontsize=7),
    release=ReleaseSettings(
        repository="pkpdsim/toy-model",
        doi="10.5281/zenodo.0000000",
        authors=[Author("Doe, Jane", orcid="0000-0000-0000-0000", affiliation="Lab")],
        keywords=["toy", "pharmacokinetics"],
        docker_image="pkpdsim/toy",
        funding="No funding.",
        description="The decay of a substance A in one compartment.",
        first_year=2026,
    ),
)
