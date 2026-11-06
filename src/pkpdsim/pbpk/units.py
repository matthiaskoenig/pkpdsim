"""The units of the PBPK models.

Time is in min, substance in mmole, volume in litre and area in m2, the
`model_units` of the models. Every unit of `U` is a unit definition of every
model; a unit outside `U` (`EXTRA_UNITS`, `ML_PER_MIN_1P73M2`) is one only of
a model which uses it.
"""

from collections.abc import Mapping
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

from sbmlutils.factory import ModelUnits, UnitDefinition, Units

if TYPE_CHECKING:
    from sbmlsim.plot import Axis

__all__ = [
    "EXTRA_UNITS",
    "HOUR",
    "KG_PER_MIN",
    "ML_PER_MIN_1P73M2",
    "MODEL_UNITS",
    "UNIT_LABELS",
    "WEEK",
    "U",
    "axis",
    "unit_label",
]


class U(Units):
    """UnitDefinitions."""

    min = UnitDefinition("min")
    mmole = UnitDefinition("mmole")
    mg = UnitDefinition("mg")
    kg = UnitDefinition("kg")
    cm = UnitDefinition("cm")
    yr = UnitDefinition("yr", "year")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    per_min = UnitDefinition("per_min", "1/min")
    per_hr = UnitDefinition("per_hr", "1/hr")
    s_per_min = UnitDefinition("s_per_min", "s/min")
    min_per_hr = UnitDefinition("min_per_hr", "min/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")
    l_per_min = UnitDefinition("l_per_min", "l/min")
    l_per_kg = UnitDefinition("l_per_kg", "l/kg")
    kg_per_l = UnitDefinition("kg_per_l", "kg/l")
    kg_per_m2 = UnitDefinition("kg_per_m2", "kg/meter^2")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    mmole_per_min_l = UnitDefinition("mmole_per_min_l", "mmole/min/l")
    mmole_min = UnitDefinition("mmole_min", "mmole*min")
    g_per_mole = UnitDefinition("g_per_mole", "g/mole")


MODEL_UNITS = ModelUnits(
    time=U.min,
    extent=U.mmole,
    substance=U.mmole,
    length=U.meter,
    area=U.m2,
    volume=U.litre,
)


ML_PER_MIN_1P73M2 = UnitDefinition(
    "ml_per_min_1p73m2", "ml/min/(1.73*meter^2)", name="ml/min/1.73 m2"
)
"""Millilitre per minute per 1.73 m2 body surface area, the unit of a clearance or a
GFR per body surface area (eGFR, `CLrenal/BSA`).

`1 ml/min/1.73 m2 = 1e-3/1.73 l/min/m2`: a quantity in l/min divided by the
body surface area `BSA` (m2) is `1730` times the number in this unit, in
consistent units `GFR / BSA * 1730 ml_per_min_1p73m2 * 1 m2 / 1 l_per_min`.
Not a unit of `U`, so the models without it keep their unit definitions; a
model defines it wherever it uses it (an observer, a parameter of the
definition, a custom reaction, a number of a formula).
"""

HOUR = UnitDefinition("hr", "hr", name="hour")
"""The hour, the unit of the clock time of the circadian rhythms (`CLOCK`, the
acrophases); a model defines it only when it has the rhythms."""

WEEK = UnitDefinition("week", "week", name="week")
"""The week, the unit of the gestational age `GA` of a pregnancy
(`pkpdsim.pbpk.pregnancy`); a model defines it only when it is pregnant."""

KG_PER_MIN = UnitDefinition("kg_per_min", "kg/min", name="kg/min")
"""Kilogram per minute, the unit of the rate rule of the body weight of a model
with an effect on the weight (`pkpdsim.pbpk.adipose`); a model defines it only
when its weight changes."""

EXTRA_UNITS: tuple[UnitDefinition, ...] = (ML_PER_MIN_1P73M2, HOUR, WEEK, KG_PER_MIN)
"""The units of the library outside `U`: a model defines one when it uses it, as
the unit of an element or by its id in a formula."""

UNIT_LABELS: Mapping[str, str] = MappingProxyType(
    {
        "ml_per_min_1p73m2": "ml/min/1.73\u00a0m²",
        "hr": "h",
        "week": "week",
        "kg_per_min": "kg/min",
    }
)
"""The text of a unit of `EXTRA_UNITS` on a figure, by the unit id.

sbmlsim writes a unit which pint does not read (a unit id of the model,
`ml_per_min_1p73m2`) on an axis as it is given; the figures show this text
instead (`unit_label`, `axis`). The words of a text are joined by a
non-breaking space (U+00A0, `1.73 m²`): sbmlsim wraps a label which is longer
than its axis at its spaces, and a unit is not broken. The unit definition
keeps its name in the SBML."""


def unit_label(unit: str | None) -> str:
    """The text of a unit on a figure.

    Args:
        unit: the unit of an axis, a unit id of the model or a unit which
            pint reads.

    Returns:
        The text of `UNIT_LABELS` for a unit of `EXTRA_UNITS` by id, else the
        short form of sbmlsim (`h` for `hr`, empty for `dimensionless`).
    """
    if unit is not None and unit in UNIT_LABELS:
        return UNIT_LABELS[unit]
    from sbmlsim.units import short_unit

    return short_unit(unit)


def axis(label: str, unit: str | None, **kwargs: Any) -> "Axis":  # noqa: ANN401
    """An axis of sbmlsim in a unit, with the text of the unit on the figure.

    The unit stays what the data and the model are converted to
    (`ml_per_min_1p73m2`, a unit of the model), the label of the axis is
    `<label> [<unit_label(unit)>]`: a clearance per body surface area reads
    `[ml/min/1.73 m²]`, not the unit id which sbmlsim writes.

    Args:
        label: the label of the axis.
        unit: the unit of the axis.
        **kwargs: further arguments of `sbmlsim.plot.Axis` (`scale`, `min`, ...).

    Returns:
        The axis, its `name` the label with the text of the unit.
    """
    from sbmlsim.plot import Axis

    text = unit_label(unit)
    name = f"{label} [{text}]" if text else label
    return Axis(label=label, unit=unit, name=name, **kwargs)
