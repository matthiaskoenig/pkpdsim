"""A model package without a model object: the conventions before pkpdsim 0.3.0.

The paths are constants of the package, the registry is `simulations.EXPERIMENTS`,
the tolerance `helpers.ABSOLUTE_TOLERANCE` and the fit definitions
`fitting.fitting.FIT_DEFINITIONS`. It reuses the toy model.
"""

from pkpdsim_toy import DATA_PATHS, MODEL_PATH, TOY_PATH

PKPDSIM_LEGACY_PATH = TOY_PATH

__all__ = ["DATA_PATHS", "MODEL_PATH", "PKPDSIM_LEGACY_PATH"]
