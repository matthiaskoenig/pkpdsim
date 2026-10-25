"""The routes of administration (data only, the SBML follows with the builder)."""

from dataclasses import dataclass
from typing import Literal

__all__ = ["Route"]


@dataclass(frozen=True)
class Route:
    """The administration of a substance.

    Attributes:
        kind: intravenous, oral tablet (with dissolution) or oral solution.
        substance: the id of the administered substance.
        formulation: the formulation of the tablet, if any.
    """

    kind: Literal["iv", "oral", "oral_solution"]
    substance: str
    formulation: object | None = None

    @classmethod
    def iv(cls, s: str) -> "Route":
        """Get the intravenous route of a substance."""
        return cls("iv", s)

    @classmethod
    def oral(cls, s: str, formulation: object | None = None) -> "Route":
        """Get the oral route of a substance, with a formulation."""
        return cls("oral", s, formulation)

    @classmethod
    def oral_solution(cls, s: str) -> "Route":
        """Get the oral route of a substance given as solution."""
        return cls("oral_solution", s)
