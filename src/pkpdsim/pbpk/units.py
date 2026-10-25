"""The units of the PBPK models.

Time is in min, substance in mmole, volume in litre and area in m2, the
`model_units` of the models.
"""

from sbmlutils.factory import ModelUnits, UnitDefinition, Units

__all__ = ["MODEL_UNITS", "U"]


class U(Units):
    """UnitDefinitions."""

    min = UnitDefinition("min")
    mmole = UnitDefinition("mmole")
    mg = UnitDefinition("mg")
    kg = UnitDefinition("kg")
    cm = UnitDefinition("cm")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    per_min = UnitDefinition("per_min", "1/min")
    l_per_min = UnitDefinition("l_per_min", "l/min")
    l_per_kg = UnitDefinition("l_per_kg", "l/kg")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    mmole_per_min_l = UnitDefinition("mmole_per_min_l", "mmole/min/l")
    g_per_mole = UnitDefinition("g_per_mole", "g/mole")


MODEL_UNITS = ModelUnits(
    time=U.min,
    extent=U.mmole,
    substance=U.mmole,
    length=U.meter,
    area=U.m2,
    volume=U.litre,
)
