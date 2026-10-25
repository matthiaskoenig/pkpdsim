"""Values with provenance.

A number of a PBPK model is a `Value`: the number, its unit and where it
comes from (a key of `data/sources.toml`, see `pkpdsim.pbpk.data`). A fitted
value is a `Value` as well, from the keyword arguments of a fit record.
"""

import functools
import typing
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Literal

from sbmlutils.factory import Uncertainty
from sbmlutils.factory.units import UnitType

__all__ = ["Q", "Value", "as_value", "check_literals", "fitted"]


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


type Q = float | Value
"""A number, plain or with unit and provenance."""


def fitted(
    kwargs: Mapping[str, Any], source: str = "fit", notes: str | None = None
) -> Value:
    """Get the value of a fitted parameter.

    Args:
        kwargs: the keyword arguments of `FitRecords.parameter(pid, unit=...)`:
            `value`, `uncertainties` and `unit`.
        source: the provenance of the value.
        notes: a remark on the value, e.g. the literature of the process; the
            provenance of the fit is in the uncertainties.

    Returns:
        The value with the uncertainty of the fit.
    """
    return Value(
        value=kwargs["value"],
        unit=kwargs["unit"],
        source=source,
        notes=notes,
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


def check_literals(obj: object) -> None:
    """Check the fields of a dataclass which are typed with a `Literal`.

    A definition written in a script or a notebook may give any text: a text
    which is none of the choices would select no branch of the builder (a
    route without doses, an erythrocyte mode treated as another), so it
    raises.

    Args:
        obj: the dataclass instance.

    Raises:
        ValueError: if a field holds none of its choices, naming the field and
            the allowed values.
    """
    cls = type(obj)
    for name, allowed in _literal_fields(cls):
        value = getattr(obj, name)
        if value not in allowed:
            raise ValueError(
                f"{cls.__name__}.{name}: {value!r} is none of the allowed values "
                f"{list(allowed)}"
            )


@functools.cache
def _literal_fields(cls: type) -> tuple[tuple[str, tuple[object, ...]], ...]:
    """The fields of a class typed with a `Literal`, with their choices."""
    return tuple(
        (name, typing.get_args(hint))
        for name, hint in typing.get_type_hints(cls).items()
        if typing.get_origin(hint) is Literal
    )
