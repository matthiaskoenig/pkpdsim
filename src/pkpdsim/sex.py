"""The sex of a model build.

A model is built for one sex: the PBPK library takes the reference values of
the sex when it builds a model (`pkpdsim.pbpk.PBPKDefinition.sex`, the ICRP 89
anatomy, the hematocrit, the partition coefficients, the GFR), so a study of
women runs on a female build and a study of men on a male one. The model
object holds a build per sex (`pkpdsim.model.PKPDModel.model_paths`, the male
build the default) and an experiment chooses one
(`pkpdsim.experiment.PKPDSimulationExperiment.sex`).

The module imports nothing else of pkpdsim, so that the model object and the
library (which imports no simulator) share the enum.
"""

from enum import StrEnum

__all__ = ["Sex"]


class Sex(StrEnum):
    """The sex of a model build, the sex of its reference individual."""

    MALE = "male"
    FEMALE = "female"
