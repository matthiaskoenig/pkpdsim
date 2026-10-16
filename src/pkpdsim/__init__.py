"""pkpdsim: the building blocks of physiologically based PK/PD models in SBML.

The models of a substance, e.g. the models of
[pkdb_models](https://github.com/matthiaskoenig/pkdb_models), are written in
python with the factory of sbmlutils and simulated with sbmlsim. What these
models share, and what is specific to the simulation of pharmacokinetics and
pharmacodynamics, lives in pkpdsim, so that the folder of a model imports
nothing but itself and released packages.

Every topic is a module of its own, which a model imports directly:
`pkpdsim.templates` (the creators and the terms of use of a model),
`pkpdsim.annotations` (the annotations of the recurring compartments),
`pkpdsim.data` (the study data of a model from the study folders of PK-DB),
`pkpdsim.initial_values` (the initial values of a substance in every tissue)
and `pkpdsim.licenses` (the license texts of the COMBINE archive of a model).
`import pkpdsim` imports none of them.
"""

__version__ = "0.7.0.dev0"

__all__ = ["__version__"]
