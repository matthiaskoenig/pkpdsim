"""The units of the models of pkpdsim.

Time is in min, substance in mmole, volume in litre and area in m2, the
`model_units` of the models (`MODEL_UNITS`). Every unit of `U` is a unit
definition of every model of the PBPK library; a unit outside `U`
(`pkpdsim.pbpk.units.EXTRA_UNITS`) is one only of a model which uses it.

`PU` adds the units of a perfused organ which `U` lacks: the volumes and the
flow of a circuit as an experimenter gives them (ml, ml/min), the mass of an
organ (g), the temperature, the oxygen partial pressure (kPa), the bile and
the urine (µl, the bile per mass µl/g), the GFR per mass (ml/min/g,
l/min/kg) and the catalytic activity of the enzymes (µkat/l = µmol/(l*s)).
Every unit of `PU` is a unit definition of every machine perfusion model
(`pkpdsim.perfusion.units`) and of every model of the zonation engine
(`pkpdsim.zonation.builder.ZonationUnits`).

SBML Level 3 has no degree Celsius: the temperature `T` is in `degC`, a unit
of the size of a kelvin, so a temperature difference is right and the
absolute temperature is `T + 273.15 degC`, which the formulas write out.
"""

from sbmlutils.factory import ModelUnits, UnitDefinition, Units

__all__ = ["MODEL_UNITS", "PU", "U"]


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


class PU(U):
    """The unit definitions of a perfusion model: those of the PBPK library and more."""

    g = UnitDefinition("g", "g")
    ml = UnitDefinition("ml", "ml")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")
    ml_per_min_g = UnitDefinition("ml_per_min_g", "ml/min/g")
    l_per_min_kg = UnitDefinition("l_per_min_kg", "liter/min/kg")
    ul = UnitDefinition("ul", "microliter", name="µl")
    ul_per_min = UnitDefinition("ul_per_min", "microliter/min", name="µl/min")
    ul_per_l = UnitDefinition("ul_per_l", "microliter/liter", name="µl/l")
    ul_per_g = UnitDefinition("ul_per_g", "microliter/g", name="µl/g")
    g_per_kg = UnitDefinition("g_per_kg", "g/kg")
    degC = UnitDefinition("degC", "kelvin", name="degree Celsius (a kelvin step)")
    per_K = UnitDefinition("per_K", "1/kelvin")
    J_per_mole = UnitDefinition("J_per_mole", "joule/mole")
    J_per_mole_K = UnitDefinition("J_per_mole_K", "joule/mole/kelvin")
    kPa = UnitDefinition("kPa", "kPa")
    mM_per_kPa = UnitDefinition("mM_per_kPa", "mmole/liter/kPa")
    ukat_per_l = UnitDefinition("ukat_per_l", "micromole/second/liter", name="µkat/l")
    ukat_per_mmole = UnitDefinition(
        "ukat_per_mmole", "micromole/second/mmole", name="µkat/mmole"
    )
