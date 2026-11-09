"""The geometry of the sinusoid, the volume fractions of the layers and the haematocrit.

The fractions are those of the liver volume (stereology, Blouin1977), the
geometry is that of one sinusoid. The numbers are in
`pkpdsim/database/tables/zonation.toml`, a section per species. The volume of
a layer is its fraction of the volume of the parenchyma, `f_par * mass /
density`, split evenly over the zones.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum

from sbmlutils.factory.units import UnitDefinition

from pkpdsim.database.tables import entry
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, Value, number
from pkpdsim.species import Species

__all__ = [
    "Layer",
    "SinusoidGeometry",
    "axial_area",
    "fractions",
    "geometry",
    "layer_volumes",
    "sinusoid_count",
    "table_value",
]

_TABLE = "zonation"

_TABLE_UNITS: dict[str, UnitDefinition] = {
    "m": U.meter,
    "1/m**2": UnitDefinition("per_m2", "1/meter^2"),
    "kg/l": U.kg_per_l,
    "dimensionless": U.dimensionless,
}
"""The unit definitions of the units of the zonation table."""

_MASS_TO_KG: dict[str, float] = {"g": 1e-3, "kg": 1.0, "mg": 1e-6}
"""The mass units a liver mass may have, in kg."""


class Layer(StrEnum):
    """The radial layers of a zone, from the lumen to the bile."""

    RBC = "rbc"
    PLASMA = "plasma"
    ENDOTHELIUM = "endothelium"
    DISSE = "disse"
    HEPATOCYTE = "hepatocyte"
    BILE = "bile"


@dataclass(frozen=True)
class SinusoidGeometry:
    """The geometry of the sinusoid of a species.

    Attributes:
        r_sin: the sinusoid radius (m).
        L_sin: the sinusoid length (m).
        t_endo: the thickness of the endothelium (m).
        t_disse: the width of the space of Disse (m).
        r_fen: the radius of a fenestra (m).
        n_fen: the number of fenestrae per endothelial area (1/m**2).
        porosity_pp: the fenestral porosity periportal (fraction).
        porosity_pc: the fenestral porosity pericentral (fraction).
        f_par: the parenchymal fraction of the liver mass (fraction).
        haematocrit: the haematocrit of the lumen (fraction).
        density: the density of the liver (kg/l).
    """

    r_sin: float
    L_sin: float
    t_endo: float
    t_disse: float
    r_fen: float
    n_fen: float
    porosity_pp: float
    porosity_pc: float
    f_par: float
    haematocrit: float
    density: float


def table_value(section: str, key: str) -> Value:
    """Get a number of the zonation table with its unit, source and notes.

    Args:
        section: the section of `tables/zonation.toml`.
        key: the key of the number.

    Returns:
        The value in the unit definition of its unit.

    Raises:
        TypeError: if the entry is no number.
        ValueError: if its unit is no unit of `_TABLE_UNITS`.
    """
    e = entry(_TABLE, section, key)
    if not isinstance(e.value, float):
        raise TypeError(f"tables/{_TABLE}.toml [{section}] {key} is no number.")
    try:
        unit = _TABLE_UNITS[e.unit]
    except KeyError:
        raise ValueError(
            f"tables/{_TABLE}.toml [{section}] {key}: unit '{e.unit}' is no unit "
            "of the zonation table."
        ) from None
    return Value(e.value, unit, e.source, e.notes)


def _section(species: Species) -> str:
    """The section of the table of a species."""
    match species:
        case Species.RAT:
            return "liver_rat"
        case Species.HUMAN:
            return "liver_human"
        case _:
            raise ValueError(
                f"The zonation table has no liver of the species {species}."
            )


def _number(section: str, key: str) -> float:
    """The number of an entry of the table."""
    return table_value(section, key).value


def fractions(
    species: Species, haematocrit: float | None = None
) -> Mapping[Layer, float]:
    """Get the volume fractions of the layers of a species.

    The lumen is split into plasma and erythrocytes by the haematocrit, which
    keeps the lumen; without erythrocytes (haematocrit 0) there is no `RBC`
    layer. Kupffer and stellate cells are no layer and left out.

    Args:
        species: the species.
        haematocrit: the haematocrit of the lumen, the one of the table if
            `None`.

    Returns:
        The fraction of the liver volume of every layer.

    Raises:
        ValueError: for a haematocrit outside [0, 1] or a species without a liver.
    """
    section = _section(species)
    hct = _number(section, "haematocrit") if haematocrit is None else haematocrit
    if not 0.0 <= hct <= 1.0:
        raise ValueError(f"The haematocrit {hct} is outside [0, 1].")
    lumen = _number(section, "lumen")
    fr: dict[Layer, float] = {}
    if hct > 0.0:
        fr[Layer.RBC] = hct * lumen
    fr[Layer.PLASMA] = (1.0 - hct) * lumen
    fr[Layer.ENDOTHELIUM] = _number(section, "endothelium")
    fr[Layer.DISSE] = _number(section, "disse")
    fr[Layer.HEPATOCYTE] = _number(section, "hepatocyte")
    fr[Layer.BILE] = _number(section, "bile")
    return fr


def geometry(species: Species) -> SinusoidGeometry:
    """Get the sinusoid geometry of a species from the table."""
    section = _section(species)
    return SinusoidGeometry(
        **{
            key: _number(section, key)
            for key in (
                "r_sin",
                "L_sin",
                "t_endo",
                "t_disse",
                "r_fen",
                "n_fen",
                "porosity_pp",
                "porosity_pc",
                "f_par",
                "haematocrit",
                "density",
            )
        }
    )


def _parenchyma_litre(mass: Q, geo: SinusoidGeometry) -> float:
    """The volume of the parenchyma (l) of a liver mass."""
    if not isinstance(mass, Value):
        raise TypeError("The liver mass needs a unit: a Value in g, kg or mg.")
    unit = getattr(mass.unit, "sid", mass.unit)
    try:
        to_kg = _MASS_TO_KG[str(unit)]
    except KeyError:
        raise ValueError(
            f"The liver mass has the unit {unit}, not g, kg or mg."
        ) from None
    m = number(mass)
    assert m is not None
    return geo.f_par * m * to_kg / geo.density


def layer_volumes(
    mass: Q, n_zones: int, geo: SinusoidGeometry, fr: Mapping[Layer, float]
) -> Mapping[tuple[Layer, int], Value]:
    """Get the volume of every layer of every zone.

    Args:
        mass: the liver mass (g, kg or mg).
        n_zones: the number of zones, which share every layer evenly.
        geo: the geometry.
        fr: the fractions of the layers.

    Returns:
        The volume (l) of the layer of a zone, the zones numbered from 1 (periportal).

    Raises:
        ValueError: for less than one zone.
    """
    if n_zones < 1:
        raise ValueError("A liver has at least one zone.")
    v_par = _parenchyma_litre(mass, geo)
    return {
        (layer, z): Value(f * v_par / n_zones, U.litre)
        for layer, f in fr.items()
        for z in range(1, n_zones + 1)
    }


def sinusoid_count(mass: Q, geo: SinusoidGeometry, fr: Mapping[Layer, float]) -> float:
    """Get the number of sinusoids which carry the lumen of a liver mass."""
    f_lumen = fr[Layer.PLASMA] + fr.get(Layer.RBC, 0.0)
    v_lumen_m3 = f_lumen * _parenchyma_litre(mass, geo) * 1e-3
    return v_lumen_m3 / (math.pi * geo.r_sin**2 * geo.L_sin)


def axial_area(
    layer: Layer, mass: Q, geo: SinusoidGeometry, fr: Mapping[Layer, float]
) -> Value:
    """Get the axial cross-section of a layer of the whole liver (m**2).

    The volume of the layer in all zones over the sinusoid length.
    """
    v_m3 = fr[layer] * _parenchyma_litre(mass, geo) * 1e-3
    return Value(v_m3 / geo.L_sin, U.m2)
