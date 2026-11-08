"""The units of the machine perfusion models.

The models have the units of the PBPK library (`pkpdsim.pbpk.units`: time in
min, substance in mmole, volume in litre) and the units of a perfusion which
the library lacks: the volumes and the flow of a circuit as an experimenter
gives them (ml, ml/min), the mass of an organ (g), the temperature, the
oxygen partial pressure (kPa), the bile and the urine (µl, the bile per
mass µl/g), the GFR per
mass (ml/min/g, l/min/kg) and the catalytic activity of
the enzymes (µkat/l = µmol/(l*s)). Every unit of `PU` is a unit definition of
every perfusion model.

SBML Level 3 has no degree Celsius: the temperature `T` is in `degC`, a unit
of the size of a kelvin, so a temperature difference is right and the
absolute temperature is `T + 273.15 degC`, which the formulas write out.
"""

from sbmlutils.factory import UnitDefinition

from pkpdsim.pbpk.units import U

__all__ = ["PU"]


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
