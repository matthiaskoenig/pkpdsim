"""The domain of a quantity of a model.

A domain is the interval of values the model accepts for a quantity: a fraction
lies in [0, 1], a volume or a rate in [0, inf). It is written into the SBML as a
`distrib` uncertainty named `DOMAIN` whose single span is of type `range`, on the
parameter, compartment or species it constrains, and read back by `domain_of` and
`domains`.

A domain also states which values of its interval the model accepts (`Values`):
every real number, the integers (an index, a count, a class, `integers`) or the
two truth values (a flag, `BOOLEAN`). An integer or boolean domain writes one
further child of the uncertainty, an `externalParameter` whose `definitionURL`
is the XSD datatype of the value set; a continuous domain writes none, so a
document written before the value sets keeps its meaning. Nothing enforces the
value set at run time, it is the contract for the tools which read it.

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
from enum import StrEnum
from pathlib import Path
from typing import overload

import libsbml
from sbmlutils.factory import Uncertainty, UncertParameter, UncertSpan

__all__ = [
    "BOOLEAN",
    "DOMAIN",
    "FRACTION",
    "NONNEGATIVE",
    "REAL",
    "Domain",
    "Values",
    "domain_of",
    "domain_uncertainty",
    "domains",
    "integers",
]

DOMAIN = "domain"
"""Name of the uncertainty which carries a domain."""


def _bound(x: float) -> str:
    """Render a bound: an integral value without the trailing `.0`, else its text."""
    return str(int(x)) if math.isfinite(x) and x == int(x) else str(x)


class Values(StrEnum):
    """The values inside the interval of a domain which the model accepts.

    A flag is kept apart from an integer domain of the same interval on
    purpose: a tool which offers a flag as a checkbox and a class as a spinner
    has to know which it has.
    """

    CONTINUOUS = "continuous"
    """Every real number of the interval, the default."""
    INTEGER = "integer"
    """The integers of the interval: an index, a count, a stage or a class."""
    BOOLEAN = "boolean"
    """The integers 0 and 1 as the truth values: a flag which is off or on."""


VALUE_SET_URL: dict[Values, str] = {
    Values.INTEGER: "http://www.w3.org/2001/XMLSchema#integer",
    Values.BOOLEAN: "http://www.w3.org/2001/XMLSchema#boolean",
}
"""The `definitionURL` which writes a value set, the XSD datatype; `CONTINUOUS`
has none, a continuous domain writes no such child."""


@dataclass(frozen=True)
class Domain:
    """Closed interval of accepted values, a bound may be infinite.

    Attributes:
        lower: Lower bound.
        upper: Upper bound.
        values: Which values of the interval are accepted.
    """

    lower: float
    upper: float
    values: Values = Values.CONTINUOUS

    @overload
    def clip(self, lower: float, upper: float) -> tuple[float, float]: ...

    @overload
    def clip(self, lower: float, upper: None = None) -> float: ...

    def clip(
        self, lower: float, upper: float | None = None
    ) -> tuple[float, float] | float:
        """Move a value or an interval which lies outside the domain onto it.

        An interval has each bound limited to the interval of the domain. A
        single value is limited to it and, for an integer or boolean domain,
        rounded to the nearest integer the domain contains.

        Args:
            lower: The value, or the lower bound of an interval.
            upper: Upper bound of the interval, `None` for a single value.

        Returns:
            The interval with each bound limited to the domain, or the value
            which `contains` accepts.
        """
        if upper is not None:
            return max(lower, self.lower), min(upper, self.upper)
        x = min(max(lower, self.lower), self.upper)
        if self.values is Values.CONTINUOUS:
            return x
        x = float(round(x))
        if x < self.lower:
            x = float(math.ceil(self.lower))
        if x > self.upper:
            x = float(math.floor(self.upper))
        return x

    def contains(self, value: float) -> bool:
        """Whether a value is one the domain accepts, the bounds included.

        An integer or boolean domain accepts only an integral value.
        """
        if not self.lower <= value <= self.upper:
            return False
        if self.values is Values.CONTINUOUS:
            return True
        return math.isfinite(value) and value == int(value)

    def __str__(self) -> str:
        """Interval notation: a finite bound is inclusive, an infinite one open.

        The value set follows when it is not continuous: `[0, 3] integer`.
        """
        left = "[" if math.isfinite(self.lower) else "("
        right = "]" if math.isfinite(self.upper) else ")"
        text = f"{left}{_bound(self.lower)}, {_bound(self.upper)}{right}"
        if self.values is Values.CONTINUOUS:
            return text
        return f"{text} {self.values.value}"


FRACTION = Domain(0.0, 1.0)
NONNEGATIVE = Domain(0.0, math.inf)
REAL = Domain(-math.inf, math.inf)
BOOLEAN = Domain(0, 1, Values.BOOLEAN)


def integers(lower: float, upper: float) -> Domain:
    """The integer domain of an interval, e.g. a class of 0 to 3.

    Args:
        lower: Lower bound.
        upper: Upper bound.

    Returns:
        The domain which accepts the integers of the interval.
    """
    return Domain(lower, upper, Values.INTEGER)


class _MetadataUncertainty(Uncertainty):
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
        An uncertainty named `DOMAIN` with one `range` span and, for an
        integer or boolean domain, one `externalParameter` whose
        `definitionURL` is the XSD datatype of the value set. The child
        carries the value 1 although only its URL means something: libsbml
        validates a bare one, but libAntimony, which the archive of a model
        loads it with, refuses it.
    """
    parameters: list[UncertParameter | UncertSpan] = [
        UncertSpan(
            type=libsbml.DISTRIB_UNCERTTYPE_RANGE,
            valueLower=float(d.lower),
            valueUpper=float(d.upper),
        )
    ]
    if d.values in VALUE_SET_URL:
        parameters.append(
            UncertParameter(
                type=libsbml.DISTRIB_UNCERTTYPE_EXTERNALPARAMETER,
                definitionURL=VALUE_SET_URL[d.values],
                value=1.0,
            )
        )
    return _MetadataUncertainty(uncertParameters=parameters, name=DOMAIN)


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
        values = _values_of(uncertainty)
        for j in range(uncertainty.getNumUncertParameters()):
            span = uncertainty.getUncertParameter(j)
            if not isinstance(span, libsbml.UncertSpan):
                continue
            if span.getType() != libsbml.DISTRIB_UNCERTTYPE_RANGE:
                continue
            lower = span.getValueLower() if span.isSetValueLower() else -math.inf
            upper = span.getValueUpper() if span.isSetValueUpper() else math.inf
            return Domain(lower, upper, values)
    return None


def _values_of(uncertainty: libsbml.Uncertainty) -> Values:
    """The value set a domain uncertainty states, by the URL of its child.

    A domain with no such child is continuous.
    """
    by_url = {url: values for values, url in VALUE_SET_URL.items()}
    for j in range(uncertainty.getNumUncertParameters()):
        child = uncertainty.getUncertParameter(j)
        if child.getType() != libsbml.DISTRIB_UNCERTTYPE_EXTERNALPARAMETER:
            continue
        values = by_url.get(child.getDefinitionURL())
        if values is not None:
            return values
    return Values.CONTINUOUS


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
