"""Pharmacodynamic components (data only, the coupling follows with the builder)."""

from collections.abc import Mapping
from dataclasses import dataclass

from sbmlutils.factory import Model

__all__ = ["PDComponent"]


@dataclass(frozen=True)
class PDComponent:
    """A model with declared ports, connected by name.

    Attributes:
        model: the model of the component.
        inputs: the ids of the input ports of the model by what they read,
            e.g. the free plasma concentration of a substance.
        outputs: the ids of the output ports of the model by what they change.
    """

    model: Model
    inputs: Mapping[str, str]
    outputs: Mapping[str, str]
