"""The container of the SBML objects of a part of a model.

A generator of a library creates the species, reactions, parameters, rules
and assignments of one part of a model and hands them over as one object
(`SbmlObjects`), which adds to the objects of the other parts; the builder
puts them into the model. The container is the SBML writing kit's, so that
the generators of the libraries hand over the same object: the machine
perfusion library knows it as `pkpdsim.perfusion.core.CoreObjects`.
"""

import dataclasses
from dataclasses import dataclass, field

from sbmlutils.factory import (
    AssignmentRule,
    InitialAssignment,
    Parameter,
    RateRule,
    Reaction,
    Species,
)

__all__ = ["SbmlObjects"]


@dataclass
class SbmlObjects:
    """The SBML objects of a part of a model.

    Attributes:
        species: the species.
        reactions: the reactions.
        parameters: the parameters.
        rules: the assignment rules.
        rate_rules: the rate rules.
        assignments: the initial assignments.
    """

    species: list[Species] = field(default_factory=list)
    reactions: list[Reaction] = field(default_factory=list)
    parameters: list[Parameter] = field(default_factory=list)
    rules: list[AssignmentRule] = field(default_factory=list)
    rate_rules: list[RateRule] = field(default_factory=list)
    assignments: list[InitialAssignment] = field(default_factory=list)

    def __add__(self, other: SbmlObjects) -> SbmlObjects:
        """Get the objects of both."""
        return SbmlObjects(
            *(
                getattr(self, f.name) + getattr(other, f.name)
                for f in dataclasses.fields(self)
            )
        )
