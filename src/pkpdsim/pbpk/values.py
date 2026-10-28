"""Values with provenance.

A number of a PBPK model is a `Value`: the number, its unit and where it
comes from (a key of `data/sources.toml`, see `pkpdsim.pbpk.data`). A fitted
value is a `Value` as well, from the keyword arguments of a fit record. A
value which is tied to other parameters (a metabolite with the Kp of its
parent, a fixed ratio) is a `Formula` over ids of the top model, written as
the initial assignment of its parameter.
"""

import functools
import typing
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import Any, Literal

from sbmlutils.factory import Uncertainty
from sbmlutils.factory.units import UnitType

__all__ = [
    "Formula",
    "FrozenMap",
    "Q",
    "Value",
    "as_value",
    "check_literals",
    "fitted",
]


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


@dataclass(frozen=True)
class Formula:
    """A value tied to other parameters: a formula over ids of the top model.

    The parameter of the value gets the formula as its initial assignment, as
    the covariates do (R21): a change of a parameter it reads before the
    simulation is followed, and an id of `PBPKDefinition.dynamic` keeps it an
    assignment rule. The ids of the formula are constant parameters and
    compartments of the top model (the parameters of the substances, the
    covariates, the physiology), also for a value of an organ, which reads
    them through input ports under the same ids; the build refuses an id the
    top model does not have and the species and observers (they change during
    the simulation). The unit is the unit of the parameter, it is not inferred
    from the formula.

    Attributes:
        expression: the formula (SBML L3 infix), e.g. `Kp_gli` or
            `ftissue_gli * fu_ref_gli / fu_ref_gli_m1`.
        unit: the unit of the parameter.
        notes: a remark on the value, e.g. why it is tied.
    """

    expression: str
    unit: UnitType
    notes: str | None = None


type Q = float | Value | Formula
"""A number, plain or with unit and provenance, or a formula."""


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


def as_value(x: Q, unit: UnitType) -> Value | Formula:
    """Get a value of a number, a value or a formula.

    Args:
        x: a number, which gets the unit, or a value or a formula, which is
            returned as it is.
        unit: the unit of a number.

    Returns:
        The value or the formula.
    """
    if isinstance(x, Value | Formula):
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


class FrozenMap(Mapping[str, Q]):
    """A read-only, hashable mapping from an id to a value.

    The mapping fields of the frozen dataclasses of a definition (the values
    per substance of the intestine, `PerSubstance`; the starts of the
    species, `PBPKDefinition.starts`) are kept as this copy: it reads back as
    a `Mapping` (`ka["a"]`, `.items()`), cannot be changed after the
    validation and is hashed by its items (a dict is unhashable).
    """

    __slots__ = ("_items",)

    def __init__(self, items: Mapping[str, Q]) -> None:
        """Copy the items."""
        self._items: dict[str, Q] = dict(items)

    def __getitem__(self, key: str) -> Q:
        """Get the value of a substance."""
        return self._items[key]

    def __iter__(self) -> Iterator[str]:
        """Iterate over the substance ids."""
        return iter(self._items)

    def __len__(self) -> int:
        """Get the number of substances."""
        return len(self._items)

    def __hash__(self) -> int:
        """Hash the items, independent of their order as the equality is."""
        return hash(frozenset(self._items.items()))

    def __eq__(self, other: object) -> bool:
        """Compare to a mapping by its items."""
        if isinstance(other, Mapping):
            return self._items == dict(other)
        return NotImplemented

    def __repr__(self) -> str:
        """Get the items as a dict."""
        return repr(self._items)
