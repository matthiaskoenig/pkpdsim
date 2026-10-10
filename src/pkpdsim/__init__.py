"""pkpdsim: the building blocks of physiologically based PK/PD models in SBML.

The models of a substance, e.g. the models of
[pkdb_models](https://github.com/matthiaskoenig/pkdb_models), are written in
python with the factory of sbmlutils and simulated with sbmlsim. What these
models share, and what is specific to the simulation of pharmacokinetics and
pharmacodynamics, lives in pkpdsim, so that the folder of a model imports
nothing but itself and released packages.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
