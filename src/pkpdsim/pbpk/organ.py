"""The interface of an organ component (data only, the SBML follows with the builder)."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Protocol

from sbmlutils.factory import Model

from pkpdsim.pbpk.chemistry import CustomReaction, Metabolism, Transport
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.proteins import Induction, Inhibition
from pkpdsim.pbpk.substance import Substance

__all__ = ["Boundary", "ModelDefinition", "OrganComponent", "OrganContext"]

type ModelDefinition = Model
"""The sbmlutils model an organ component defines."""


@dataclass(frozen=True)
class OrganContext:
    """What an organ component is built from.

    Attributes:
        organ: the organ.
        substances: the substances of the definition.
        parts: the parts of the organ by substance id (`Any` until the parts
            exist, tightened with the standard organ).
        chemistry: the processes of the definition.
        inhibitions: the inhibitions as (inhibitor id, inhibition).
        inductions: the inductions as (inducer id, induction).
    """

    organ: Organ
    substances: tuple[Substance, ...]
    parts: Mapping[str, Any]
    chemistry: tuple[Metabolism | Transport | CustomReaction, ...]
    inhibitions: tuple[tuple[str, Inhibition], ...]
    inductions: tuple[tuple[str, Induction], ...]


@dataclass(frozen=True)
class Boundary:
    """What an organ component exposes to the body.

    Attributes:
        replacements: the ids of the submodel replaced by the body.
        outputs: the ids of the outputs of the organ by name.
    """

    replacements: Mapping[str, str]
    outputs: Mapping[str, str]


class OrganComponent(Protocol):
    """An organ of the body, the standard one or a variant."""

    name: str

    def build(self, ctx: OrganContext) -> tuple[ModelDefinition, Boundary]:
        """Build the model of the organ and its boundary."""
        ...
