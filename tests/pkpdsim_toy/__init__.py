"""A toy model package: the decay of a substance A in one compartment.

The tests of pkpdsim run the model object, the commands, the verification and the
test functions on it. The SBML model is created by `factory.create_toy_models`
before the tests start (`tests/conftest.py`) and is not tracked.
"""

from pathlib import Path

TOY_PATH = Path(__file__).parent
MODEL_BASE_PATH = TOY_PATH / "models"
MODEL_PATH = MODEL_BASE_PATH / "toy.xml"
DATA_PATHS = [TOY_PATH]
