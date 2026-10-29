"""The units of the PBPK models.

Time is in min, substance in mmole, volume in litre and area in m2, the
`model_units` of the models. Every unit of `U` is a unit definition of every
model; a unit outside `U` (`ML_PER_MIN_1P73M2`) is one only of a model whose
observers use it.
"""

from sbmlutils.factory import ModelUnits, UnitDefinition, Units

__all__ = ["ML_PER_MIN_1P73M2", "MODEL_UNITS", "U"]


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
model defines it when an observer of its definition uses it.
"""
