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

`PU` is the SBML writing kit's (`pkpdsim.sbml.units.OrganUnits`, of which `PU`
is the alias this library has always called it by): the zonation engine writes
the organs of a perfusion and needs the same units, and the engine is the
lower layer. This module is the import path of this library and re-exports the
class unchanged, so every model of this library is the model it was.
"""

from pkpdsim.sbml.units import PU

__all__ = ["PU"]
