"""The numbers of the PBPK library from the parameter database.

The library holds no numbers of its own: they are the tables and the studies
of `pkpdsim.database`. This module reads them for the library, `value` as a
`Value` with the unit of the model, its source and its uncertainty.
"""

import libsbml
from sbmlutils.factory import Uncertainty as SbmlUncertainty
from sbmlutils.factory import UncertParameter, UncertSpan, UnitDefinition

from pkpdsim.database.tables import Source, Uncertainty, cite, entry, load, sources
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value

__all__ = ["Source", "cite", "load", "sbml_uncertainty", "sources", "value"]

_UNITS: dict[str, UnitDefinition] = {
    "dimensionless": U.dimensionless,
    "l": U.litre,
    "s": U.second,
    "min": U.min,
    "1/min": U.per_min,
    "1/hr": U.per_hr,
}
"""The units of the model for the units of the tables read by `value`."""


def sbml_uncertainty(
    uncertainty: Uncertainty, unit: UnitDefinition, source: str
) -> SbmlUncertainty:
    """Get the SBML uncertainty of a number of the database.

    Args:
        uncertainty: the uncertainty the source states.
        unit: the unit of the number in the model.
        source: the key of the source.

    Returns:
        The distrib uncertainty: the standard deviation, the standard error,
        the coefficient of variation, the sample size and the 95 % confidence
        interval which the source states.
    """
    children: list[UncertParameter | UncertSpan] = []
    if uncertainty.sd is not None:
        children.append(
            UncertParameter(
                type=libsbml.DISTRIB_UNCERTTYPE_STANDARDDEVIATION,
                value=uncertainty.sd,
                unit=unit,
            )
        )
    if uncertainty.se is not None:
        children.append(
            UncertParameter(
                type=libsbml.DISTRIB_UNCERTTYPE_STANDARDERROR,
                value=uncertainty.se,
                unit=unit,
            )
        )
    if uncertainty.cv is not None:
        children.append(
            UncertParameter(
                type=libsbml.DISTRIB_UNCERTTYPE_COEFFIENTOFVARIATION,
                value=uncertainty.cv,
                unit=U.dimensionless,
            )
        )
    if uncertainty.n is not None:
        children.append(
            UncertParameter(
                type=libsbml.DISTRIB_UNCERTTYPE_SAMPLESIZE,
                value=float(uncertainty.n),
                unit=U.dimensionless,
            )
        )
    if uncertainty.ci is not None:
        children.append(
            UncertSpan(
                type=libsbml.DISTRIB_UNCERTTYPE_CONFIDENCEINTERVAL,
                valueLower=uncertainty.ci[0],
                valueUpper=uncertainty.ci[1],
                unit=unit,
                name="95 % confidence interval",
            )
        )
    return SbmlUncertainty(uncertParameters=children, name=f"source {source}")


def value(name: str, section: str, key: str) -> Value:
    """Get a number of a table with its unit, source, notes and uncertainty.

    Args:
        name: the name of the table without extension, e.g. `routes`.
        section: the table of the number.
        key: the key of the number in the table.

    Returns:
        The value, with the uncertainty its source states.

    Raises:
        ValueError: if the unit of the number is not a unit of the model or
            the number is a vector.
    """
    e = entry(name, section, key)
    try:
        unit = _UNITS[e.unit]
    except KeyError:
        raise ValueError(
            f"tables/{name}.toml [{section}] {key}: unit '{e.unit}' is not a unit "
            "of the model."
        ) from None
    if isinstance(e.value, tuple):
        raise ValueError(f"tables/{name}.toml [{section}] {key}: a vector, no number.")
    uncertainties = (
        [sbml_uncertainty(e.uncertainty, unit, e.source)]
        if e.uncertainty is not None
        else None
    )
    return Value(e.value, unit, e.source, e.notes, uncertainties)
