"""The domain of a quantity of a model: the writer of the contract.

A domain is the interval of values the model accepts for a quantity, and the
value set inside that interval: a fraction lies in [0, 1], a volume or a rate
in [0, inf), a severity in the integers of [0, 3], a switch is 0 or 1. It is
written into the SBML as a `distrib` uncertainty named `DOMAIN` whose single
span is of type `range`, on the parameter, compartment or species it
constrains, with an `externalParameter` whose `definitionURL` is the XSD
datatype of the value set for an integer or boolean domain; a continuous domain
writes none, so a document written before the value sets keeps its meaning.
Nothing enforces the value set at run time, it is the contract for the tools
which read it.

pkpdsim writes the models, so `domain_uncertainty` lives here, and the reader of
the contract is the one of sbmlsim (`sbmlsim.domains`), which every tool that
samples, scans, fits or checks a model reads it with: `sampling.local` of
sbmlsim clamps the step of a parameter to the domain of its target, and a design
must not import pkpdsim. `Domain`, `Values`, `domain_of` and `domains` are that
reader, re-exported here next to the writer, so one contract has one
implementation of its reader (#245). The named domains of the library and
`integers` are the vocabulary of the writer and stay here.

A domain is not the uncertainty of a source (`pkpdsim.pbpk.data.sbml_uncertainty`,
a confidence interval around a literature value) and not the distribution of a
fit: those are statements in their own uncertainty elements, which coexist with
the domain on the same element. The reader finds the domain by its name and its
span type and ignores every other uncertainty.

A reader of a model which states no domain learns nothing: `domain_of` gives
`None` and `domains` leaves the element out, there is no default. This module
imports libsbml, sbmlutils and the reader of sbmlsim, which imports libsbml and
the standard library only, no simulator and no library of pkpdsim: it belongs to
the layer below them (see `pkpdsim.sbml`).
"""

import math

import libsbml
from sbmlsim.domains import (
    DOMAIN,
    VALUE_SET_URL,
    Domain,
    Values,
    domain_of,
    domains,
)
from sbmlutils.factory import Uncertainty, UncertParameter, UncertSpan

__all__ = [
    "BOOLEAN",
    "DOMAIN",
    "FRACTION",
    "NONNEGATIVE",
    "REAL",
    "VALUE_SET_URL",
    "Domain",
    "Values",
    "domain_of",
    "domain_uncertainty",
    "domains",
    "integers",
]

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

    Raises:
        ValueError: If no integer lies in `[lower, upper]`.
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
        does not need it, but libAntimony, which the archive of a model loads
        it with, fails without it (sys-bio/antimony#190).
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
