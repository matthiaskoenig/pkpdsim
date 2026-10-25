"""Values with provenance.

A number of a PBPK model is a `Value`: the number, its unit and where it
comes from (a key of `data/sources.toml`, see `pkpdsim.pbpk.data`). A fitted
value is a `Value` as well, from the keyword arguments of a fit record.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from sbmlutils.factory import Uncertainty
from sbmlutils.factory.units import UnitType

__all__ = ["Value", "as_value", "fitted"]


@dataclass(frozen=True)
class Value:
    """A number with unit and provenance.

    Attributes:
        value: the number.
        unit: the unit of the number.
        source: the key of the source in `data/sources.toml`, `fit` for a
            fitted value.
        notes: a remark on the value.
        uncertainties: the uncertainty of a fitted value.
    """

    value: float
    unit: UnitType
    source: str | None = None
    notes: str | None = None
    uncertainties: list[Uncertainty] | None = None


def fitted(kwargs: Mapping[str, Any], source: str = "fit") -> Value:
    """Get the value of a fitted parameter.

    Args:
        kwargs: the keyword arguments of `FitRecords.parameter(pid, unit=...)`:
            `value`, `uncertainties` and `unit`.
        source: the provenance of the value.

    Returns:
        The value with the uncertainty of the fit.
    """
    return Value(
        value=kwargs["value"],
        unit=kwargs["unit"],
        source=source,
        uncertainties=kwargs.get("uncertainties"),
    )


def as_value(x: float | Value, unit: UnitType) -> Value:
    """Get a value of a number or a value.

    Args:
        x: a number, which gets the unit, or a value, which is returned as it is.
        unit: the unit of a number.

    Returns:
        The value.
    """
    if isinstance(x, Value):
        return x
    return Value(x, unit)
