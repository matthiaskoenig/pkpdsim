"""The species of a model build.

A model is built for one species: the PBPK library takes the reference
physiology of the species when it builds a model
(`pkpdsim.pbpk.PBPKDefinition.species`, `database/tables/physiology_<species>.toml`),
the human of ICRP 89 by default and the laboratory animals of the preclinical
studies (rat, mouse, pig) with their reference animal.

The module imports nothing else of pkpdsim, so that the model object and the
library (which imports no simulator) share the enum.
"""

from enum import StrEnum

__all__ = ["Species"]


class Species(StrEnum):
    """The species of a model build, the species of its reference individual."""

    HUMAN = "human"
    RAT = "rat"
    MOUSE = "mouse"
    PIG = "pig"
