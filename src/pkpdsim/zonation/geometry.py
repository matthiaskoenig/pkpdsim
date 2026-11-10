"""The geometry of the sinusoid, the volume fractions of the layers and the haematocrit.

The fractions are those of the liver parenchyma (stereology, Blouin1977), the
geometry is that of one sinusoid. The two are two sources and not one
cylinder: a shell of the endothelial thickness `t_endo` (0.165 um, the
representative sinusoid of the predecessor spt model) around the lumen of
radius 4.4 um is 0.81 % of the parenchyma where the stereology assigns the
endothelium 2.8 %, a factor of 3.5 (the thickness consistent with the
fraction would be 0.547 um; the 1.0 um shell of the space of Disse gives
5.55 % against the 4.9 % of the table, which is close). So the fenestral
conductance, whose path length is the thin wall
(`transport.fenestral_permeability`), is about 3.3 times the one consistent
with the volume it feeds. Nothing rests on it today: the transcellular
permeabilities are 0 until a model sets them, and for a small solute neither
conductance is limiting. The numbers are in
`pkpdsim/database/tables/zonation.toml`, a section per species. The volume of
a layer is its fraction of the volume of the parenchyma, `f_par * mass /
density`, split evenly over the zones.

**Every number keeps its provenance.** A number of the table is a `Value`
with the source and the notes of its row, which travels into the notes of the
parameter the builder writes from it (`pkpdsim.sbml.parameters.value_parameter`,
the convention of the whole library), so a released model says where its
stereology, its sinusoid geometry and its fenestral porosity come from and a
reader can tell a measurement from a choice. A number a caller hands over
instead of the table (`SinusoidGeometry`, `ZonationSettings.fractions`) is
wrapped as the model choice it is (`model_choice`) and never looks sourced,
and a number computed from several rows (the plasma and the erythrocytes of
the lumen) carries the sources of the rows it rests on (`derived`).
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, fields
from enum import StrEnum

from sbmlutils.factory.units import UnitDefinition

from pkpdsim.database.tables import entry
from pkpdsim.sbml.units import U
from pkpdsim.sbml.values import (
    Q,
    Value,
    derived,
    model_choice,
    number,
    unit_text,
)
from pkpdsim.species import Species
from pkpdsim.zonation.profiles import zone_positions

__all__ = [
    "Layer",
    "SinusoidGeometry",
    "axial_area",
    "bile_flow",
    "field_unit",
    "fractions",
    "geometry",
    "layer_volumes",
    "lumen_fraction",
    "lumen_haematocrit",
    "mass_kg",
    "porosities",
    "reference_liver_mass",
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
"""The mass units a liver mass may have, in kg (`mass_kg`)."""


class Layer(StrEnum):
    """The radial layers of a zone, from the lumen to the bile."""

    RBC = "rbc"
    PLASMA = "plasma"
    ENDOTHELIUM = "endothelium"
    DISSE = "disse"
    CELL = "cell"
    """The cells of the organ, the parenchymal cells (for the liver the hepatocytes).

    The engine serves any zoned organ, so the layer, its compartment `V_cell_<zone>` and
    its fraction `f_cell` carry no organ name."""
    BILE = "bile"


_GEOMETRY_UNITS: Mapping[str, UnitDefinition] = {
    "r_sin": U.meter,
    "L_sin": U.meter,
    "t_endo": U.meter,
    "t_disse": U.meter,
    "r_fen": U.meter,
    "n_fen": _TABLE_UNITS["1/m**2"],
    "porosity_pp": U.dimensionless,
    "porosity_pc": U.dimensionless,
    "f_par": U.dimensionless,
    "haematocrit": U.dimensionless,
    "density": U.kg_per_l,
}
"""The unit of every field of `SinusoidGeometry`, the unit of its row of the table."""


def field_unit(name: str) -> UnitDefinition:
    """Get the unit of a field of `SinusoidGeometry` by its name.

    Raises:
        KeyError: for a name which is no field of the geometry.
    """
    try:
        return _GEOMETRY_UNITS[name]
    except KeyError:
        raise KeyError(f"'{name}' is no field of the sinusoid geometry.") from None


@dataclass(frozen=True)
class SinusoidGeometry:
    """The geometry of the sinusoid of a species, every number with its provenance.

    A field is the `Value` of its row of `tables/zonation.toml` (`geometry`),
    with the source and the notes of the row, which reach the notes of the
    parameter the builder writes from it. A caller which replaces a field
    (`dataclasses.replace`) may give a plain number, which is wrapped as the
    model choice it is (`model_choice`, the source `ModelChoice` and a note
    saying the caller supplied it) and never looks sourced.

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

    r_sin: Value
    L_sin: Value
    t_endo: Value
    t_disse: Value
    r_fen: Value
    n_fen: Value
    porosity_pp: Value
    porosity_pc: Value
    f_par: Value
    haematocrit: Value
    density: Value

    def __post_init__(self) -> None:
        """Wrap a plain number a caller gave as the model choice it is."""
        for f in fields(self):
            given = getattr(self, f.name)
            if not isinstance(given, Value):
                object.__setattr__(
                    self,
                    f.name,
                    model_choice(
                        float(given),
                        field_unit(f.name),
                        f"{f.name} of the sinusoid geometry, supplied by the "
                        "caller of the engine and no value of a table",
                    ),
                )


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


def lumen_haematocrit(species: Species, haematocrit: float | None = None) -> Value:
    """Get the haematocrit which splits the lumen into plasma and erythrocytes.

    Args:
        species: the species.
        haematocrit: the haematocrit a caller gives, which is the model choice
            it is; `None` takes the row of the table with its source.

    Returns:
        The haematocrit.

    Raises:
        ValueError: for a haematocrit outside [0, 1] or a species without a liver.
    """
    hct = (
        table_value(_section(species), "haematocrit")
        if haematocrit is None
        else model_choice(
            haematocrit,
            U.dimensionless,
            "the haematocrit of the sinusoidal lumen, supplied by the caller "
            "of the engine and no value of a table",
        )
    )
    if not 0.0 <= hct.value <= 1.0:
        raise ValueError(f"The haematocrit {hct.value} is outside [0, 1].")
    return hct


def fractions(
    species: Species, haematocrit: float | None = None
) -> Mapping[Layer, Value]:
    """Get the volume fractions of the layers of a species, with their provenance.

    The lumen is split into plasma and erythrocytes by the haematocrit, which
    keeps the lumen; without erythrocytes (haematocrit 0) there is no `RBC`
    layer. The two shares of the lumen rest on two rows (the lumen and the
    haematocrit), so they are derived values naming both (`derived`); the
    other layers are the rows of the table unchanged. Kupffer and stellate
    cells are no layer and left out.

    Args:
        species: the species.
        haematocrit: the haematocrit of the lumen, the one of the table if
            `None`.

    Returns:
        The fraction of the parenchyma of every layer.

    Raises:
        ValueError: for a haematocrit outside [0, 1] or a species without a liver.
    """
    section = _section(species)
    hct = lumen_haematocrit(species, haematocrit)
    lumen = table_value(section, "lumen")
    parts = (lumen, hct)
    fr: dict[Layer, Value] = {}
    if hct.value > 0.0:
        fr[Layer.RBC] = derived(
            hct.value * lumen.value,
            U.dimensionless,
            parts,
            f"the erythrocytes of the sinusoidal lumen, the haematocrit "
            f"{hct.value:g} of the lumen. Lumen: {lumen.notes}. Haematocrit: "
            f"{hct.notes}",
        )
    fr[Layer.PLASMA] = derived(
        (1.0 - hct.value) * lumen.value,
        U.dimensionless,
        parts,
        f"the plasma of the sinusoidal lumen, the lumen less its erythrocytes "
        f"(haematocrit {hct.value:g}). Lumen: {lumen.notes}. Haematocrit: "
        f"{hct.notes}",
    )
    fr[Layer.ENDOTHELIUM] = table_value(section, "endothelium")
    fr[Layer.DISSE] = table_value(section, "disse")
    fr[Layer.CELL] = table_value(section, "hepatocyte")
    fr[Layer.BILE] = table_value(section, "bile")
    return fr


def reference_liver_mass(species: Species) -> Value:
    """Get the mass (kg) of the reference liver of a species.

    The liver of the reference male of `physiology_<species>.toml`, the liver
    the sourced flows of the zonation table (`bile_flow`) belong to.
    """
    e = entry(f"physiology_{species}", "male.organs.li", "mass")
    if not isinstance(e.value, float):
        raise TypeError(f"The reference liver mass of the {species} is no number.")
    return Value(e.value, U.kg, e.source, e.notes)


def bile_flow(species: Species) -> Value:
    """Get the canalicular bile flow of the reference liver of a species (l/min).

    A sourced value of each species, never scaled from the other one: the rat
    secretes about 14 times more bile per g of liver than the human canalicular flow, an upper bound on the species difference (the rat value includes a ductular contribution; like for like with the human total it is about 9 times; see the
    entries `bile_flow` of `tables/zonation.toml`). The flow of a liver of
    another mass is this flow times the mass over `reference_liver_mass`.

    Args:
        species: the species.

    Returns:
        The flow (l/min) of the liver of `reference_liver_mass`.

    Raises:
        ValueError: if the table gives the flow in no unit of ml/min and l/min.
    """
    e = entry(_TABLE, _section(species), "bile_flow")
    if not isinstance(e.value, float):
        raise TypeError(f"tables/{_TABLE}.toml bile_flow is no number.")
    try:
        factor = {"ml/min": 1e-3, "l/min": 1.0}[e.unit]
    except KeyError:
        raise ValueError(
            f"The bile flow has the unit {e.unit}, not ml/min or l/min."
        ) from None
    return Value(e.value * factor, U.l_per_min, e.source, e.notes)


def geometry(species: Species) -> SinusoidGeometry:
    """Get the sinusoid geometry of a species from the table.

    Every field is the `Value` of its row, with the source and the notes of
    the row (`table_value`): the provenance travels into the parameters the
    builder and the flow write from the geometry.
    """
    section = _section(species)
    return SinusoidGeometry(
        **{key: table_value(section, key) for key in _GEOMETRY_UNITS}
    )


def mass_kg(mass: Q) -> float:
    """Get the mass of an organ in kg, the one conversion of the engine.

    The engine reads a liver mass in one place, so that a mass in an unexpected
    unit is refused with one message wherever it enters (the volumes of the
    layers, the number of sinusoids, the flow of the chain).

    Args:
        mass: the mass, a `Value` in g, kg or mg.

    Returns:
        The mass in kg.

    Raises:
        TypeError: for a mass without a unit or without a number.
        ValueError: for a mass whose unit is no mass of `_MASS_TO_KG`.
    """
    if not isinstance(mass, Value):
        raise TypeError(
            f"The liver mass needs a unit: a Value in {sorted(_MASS_TO_KG)}, "
            f"not {mass!r}."
        )
    unit = unit_text(mass.unit)
    try:
        to_kg = _MASS_TO_KG[unit]
    except KeyError:
        raise ValueError(
            f"The liver mass has the unit {unit}, not one of {sorted(_MASS_TO_KG)}."
        ) from None
    m = number(mass)
    if m is None:
        raise TypeError(f"The liver mass needs a number, not a formula: {mass!r}.")
    return m * to_kg


def _parenchyma_litre(mass: Q, geo: SinusoidGeometry) -> float:
    """The volume of the parenchyma (l) of a liver mass."""
    return geo.f_par.value * mass_kg(mass) / geo.density.value


def layer_volumes(
    mass: Q, n_zones: int, geo: SinusoidGeometry, fr: Mapping[Layer, Value]
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
        (layer, z): Value(f.value * v_par / n_zones, U.litre)
        for layer, f in fr.items()
        for z in range(1, n_zones + 1)
    }


def porosities(n_zones: int, geo: SinusoidGeometry) -> tuple[Value, ...]:
    """Get the fenestral porosity of every zone, the gradient of the table.

    Args:
        n_zones: the number of zones.
        geo: the geometry with the periportal and the pericentral porosity.

    Returns:
        The porosity of every zone, periportal first: the linear gradient
        from `porosity_pp` to `porosity_pc` over the lobular positions of the
        centres of the zones (`zone_positions`), each a value derived from the
        two ends with their provenance.

    Raises:
        ValueError: for less than one zone.
    """
    pp, pc = geo.porosity_pp, geo.porosity_pc
    return tuple(
        derived(
            pp.value + (pc.value - pp.value) * x,
            U.dimensionless,
            (pp, pc),
            f"the fenestral porosity of the endothelium at the lobular "
            f"position {x:g} of the centre of the zone, the linear gradient "
            f"from {pp.value:g} periportal to {pc.value:g} pericentral. "
            f"Periportal: {pp.notes}. Pericentral: {pc.notes}",
        )
        for x in zone_positions(n_zones)
    )


def sinusoid_count(mass: Q, geo: SinusoidGeometry, fr: Mapping[Layer, Value]) -> float:
    """Get the number of sinusoids which carry the lumen of a liver mass."""
    v_lumen_m3 = lumen_fraction(fr) * _parenchyma_litre(mass, geo) * 1e-3
    return v_lumen_m3 / (math.pi * geo.r_sin.value**2 * geo.L_sin.value)


def lumen_fraction(fr: Mapping[Layer, Value]) -> float:
    """The volume fraction of the lumen, the plasma and the erythrocytes of it."""
    rbc = fr.get(Layer.RBC)
    return fr[Layer.PLASMA].value + (0.0 if rbc is None else rbc.value)


def axial_area(
    layer: Layer, mass: Q, geo: SinusoidGeometry, fr: Mapping[Layer, Value]
) -> Value:
    """Get the axial cross-section of a layer of the whole liver (m**2).

    The volume of the layer in all zones over the sinusoid length.

    Raises:
        ValueError: for a layer without a volume in `fr`, the erythrocytes at
            a haematocrit of 0.
    """
    if layer not in fr:
        raise ValueError(
            f"The layer {layer} has no volume at this haematocrit, so it has no "
            f"axial area (layers with a volume: {sorted(str(x) for x in fr)})."
        )
    v_m3 = fr[layer].value * _parenchyma_litre(mass, geo) * 1e-3
    return Value(v_m3 / geo.L_sin.value, U.m2)
