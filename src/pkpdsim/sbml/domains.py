"""The domain of a quantity of a model.

A domain is the interval of values the model accepts for a quantity: a fraction
lies in [0, 1], a volume or a rate in [0, inf). It is written into the SBML as a
`distrib` uncertainty named `DOMAIN` whose single span is of type `range`, on the
parameter, compartment or species it constrains, and read back by `domain_of` and
`domains`.

A domain is not the uncertainty of a source (`pkpdsim.pbpk.data.sbml_uncertainty`,
a confidence interval around a literature value) and not the distribution of a
fit: those are statements in their own uncertainty elements, which coexist with
the domain on the same element. The reader finds the domain by its name and its
span type and ignores every other uncertainty.

A reader of a model which states no domain learns nothing: `domain_of` gives
`None` and `domains` leaves the element out, there is no default. This module
imports no simulator and no library of pkpdsim, it belongs to the layer below
them (see `pkpdsim.sbml`).
"""

import math
from dataclasses import dataclass
from pathlib import Path

import libsbml
from sbmlutils.factory import Uncertainty, UncertSpan

__all__ = [
    "DOMAIN",
    "FRACTION",
    "NONNEGATIVE",
    "REAL",
    "Domain",
    "domain_of",
    "domain_uncertainty",
    "domains",
]

DOMAIN = "domain"
"""Name of the uncertainty which carries a domain."""


def _bound(x: float) -> str:
    """Render a bound: an integral value without the trailing `.0`, else its text."""
    return str(int(x)) if math.isfinite(x) and x == int(x) else str(x)


@dataclass(frozen=True)
class Domain:
    """Closed interval of accepted values, a bound may be infinite."""

    lower: float
    upper: float

    def clip(self, lower: float, upper: float) -> tuple[float, float]:
        """Move a bound of an interval which lies outside the domain onto it.

        Args:
            lower: Lower bound of the interval.
            upper: Upper bound of the interval.

        Returns:
            The interval with each bound limited to the domain.
        """
        return max(lower, self.lower), min(upper, self.upper)

    def contains(self, value: float) -> bool:
        """Whether a value lies in the domain, the bounds included."""
        return self.lower <= value <= self.upper

    def __str__(self) -> str:
        """Interval notation: a finite bound is inclusive, an infinite one open."""
        left = "[" if math.isfinite(self.lower) else "("
        right = "]" if math.isfinite(self.upper) else ")"
        return f"{left}{_bound(self.lower)}, {_bound(self.upper)}{right}"


FRACTION = Domain(0.0, 1.0)
NONNEGATIVE = Domain(0.0, math.inf)
REAL = Domain(-math.inf, math.inf)


class MetadataUncertainty(Uncertainty):
    """An uncertainty which asks for no SBO term.

    An SBO term classifies a component of a model, an uncertainty is metadata
    about a value; sbmlutils says so for `UncertParameter` and `UncertSpan`
    and this is the local stand-in until `Uncertainty` says it too.
    """

    _hint_sbo_term = False


def domain_uncertainty(d: Domain) -> Uncertainty:
    """Build the uncertainty which writes a domain.

    Args:
        d: The domain.

    Returns:
        An uncertainty named `DOMAIN` with one `range` span.
    """
    return MetadataUncertainty(
        uncertParameters=[
            UncertSpan(
                type=libsbml.DISTRIB_UNCERTTYPE_RANGE,
                valueLower=d.lower,
                valueUpper=d.upper,
            )
        ],
        name=DOMAIN,
    )


def domain_of(sbase: libsbml.SBase) -> Domain | None:
    """Read the domain of an element.

    Args:
        sbase: A parameter, compartment or species.

    Returns:
        The domain, or None if the element states none. A bound which is not
        set reads as the infinity on its side.
    """
    plugin = sbase.getPlugin("distrib")
    if plugin is None:
        return None
    for i in range(plugin.getNumUncertainties()):
        uncertainty = plugin.getUncertainty(i)
        if uncertainty.getName() != DOMAIN:
            continue
        for j in range(uncertainty.getNumUncertParameters()):
            span = uncertainty.getUncertParameter(j)
            if not isinstance(span, libsbml.UncertSpan):
                continue
            if span.getType() != libsbml.DISTRIB_UNCERTTYPE_RANGE:
                continue
            lower = span.getValueLower() if span.isSetValueLower() else -math.inf
            upper = span.getValueUpper() if span.isSetValueUpper() else math.inf
            return Domain(lower, upper)
    return None


def domains(model_path: Path) -> dict[str, Domain]:
    """Read the domains of the parameters, compartments and species of a model.

    Args:
        model_path: Path of an SBML file.

    Returns:
        The domain by id, for the elements which state one.
    """
    doc = libsbml.readSBMLFromFile(str(model_path))
    model = doc.getModel()
    result: dict[str, Domain] = {}
    if model is None:
        return result
    for elements in (
        model.getListOfParameters(),
        model.getListOfCompartments(),
        model.getListOfSpecies(),
    ):
        for element in elements:
            domain = domain_of(element)
            if domain is not None:
                result[element.getId()] = domain
    return result
