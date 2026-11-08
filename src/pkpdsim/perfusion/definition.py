"""The definition of a machine perfusion: the circuit, the perfusate and the organ.

A machine perfusion model is built from a `PerfusionDefinition`, as a PBPK
model from a `pkpdsim.pbpk.PBPKDefinition`: the circuit (`Circuit`, the
reservoir, the tubing, the pump flow and the oxygenator), the perfusate
(`Perfusate`, the fresh composition, the oxygenation target, the oxygen
carrier and the inert tracers), the mounted organ (a chain of zones with the
cellular core: `LiverOrgan` with bile, `KidneyOrgan` with the filtration and
the urine, `GenericOrgan` without excretion; each releases its quantities of
`ids.RELEASABLE`), the temperature and the constants of the cellular core
(`pkpdsim.perfusion.core.CoreParameters`). `pkpdsim.perfusion.builder.build`
writes its SBML.

The numbers of a definition are `Value`s in the units an experimenter gives
them (ml, ml/min, g, kPa, degC, `pkpdsim.perfusion.units.PU`), which the
build converts into the units of the model; `validate` refuses a value in
another unit. The defaults (`Circuit.default`, `Perfusate.default`,
`LiverOrgan.default`, `KidneyOrgan.default`) are the numbers of the
parameter database (`tables/perfusion.toml`), `GenericOrgan.default` those of
an organ of the physiology table of the species.
"""

import dataclasses
import math
import re
from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from dataclasses import KW_ONLY, dataclass, field
from typing import TYPE_CHECKING, ClassVar, Literal, Self

from sbmlutils.factory import Creator
from sbmlutils.factory.units import UnitDefinition, UnitType

from pkpdsim.database.tables import entry
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import SID, FrozenMap, Value, check_literals
from pkpdsim.perfusion import ids
from pkpdsim.perfusion.core import PROCESSES, CoreParameters
from pkpdsim.perfusion.units import PU
from pkpdsim.species import Species

__all__ = [
    "PH_RANGE",
    "Circuit",
    "GenericOrgan",
    "KidneyOrgan",
    "LiverOrgan",
    "Organ",
    "Perfusate",
    "PerfusionDefinition",
    "table_value",
]

_TABLE = "perfusion"

PH_RANGE: tuple[float, float] = (6.0, 8.0)
"""The range of the pH of a fresh perfusate which `validate` accepts.

A perfusate is buffered near the physiological pH (7.4); a pH outside of 6 to
8 is a typo (a pH of 74 or of a hydrogen ion concentration) and not a
perfusion, and the protons of the model would start far from any condition of
the constants of the core.
"""

_TABLE_UNITS: dict[str, UnitDefinition] = {
    "ml": PU.ml,
    "ml/min": PU.ml_per_min,
    "ml/min/g": PU.ml_per_min_g,
    "g": PU.g,
    "kg/l": U.kg_per_l,
    "kPa": PU.kPa,
    "mM": U.mM,
    "mM/kPa": PU.mM_per_kPa,
    "1/min": U.per_min,
    "1/K": PU.per_K,
    "degC": PU.degC,
    "J/mol/K": PU.J_per_mole_K,
    "ukat/mmol": PU.ukat_per_mmole,
    "dimensionless": U.dimensionless,
}
"""The unit definitions of the units of the perfusion table."""


def table_value(section: str, key: str) -> Value:
    """Get a number of the perfusion table with its unit, source and notes.

    Args:
        section: the section of `tables/perfusion.toml`.
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
            "of a perfusion definition."
        ) from None
    return Value(e.value, unit, e.source, e.notes)


@dataclass(frozen=True)
class Perfusate:
    """The fresh perfusate of a circuit.

    The perfusate fills the reservoir, the tubing and the vascular space of
    the organ at the start, and a medium change brings it in again; the
    model has its composition as the constants `<q>_fresh` (mM) of every
    quantity of the perfusate (`ids.PERFUSATE`, the tracers), the free
    protons of its pH and the dissolved CO2 of the buffer equilibrium; the
    free oxygen is that of `o2_target_kpa` at the temperature.

    Attributes:
        glucose: glucose (mM).
        lactate: lactate (mM).
        hco3: bicarbonate (mM).
        k: potassium (mM).
        o2_target_kpa: the oxygen partial pressure which the oxygenator sets
            and of the fresh perfusate (kPa).
        ph: the pH of the fresh perfusate, within `PH_RANGE` (6 to 8).
        carrier: the oxygen carrier, `none` (an acellular perfusate, the free
            oxygen only) or `rbc` (red blood cells, declared and not built
            yet, `validate` refuses it).
        hct: the hematocrit of a perfusate with red blood cells.
        tracers: inert quantities of the perfusate (SIds), which flow with it
            and do not enter the cells (a mass balance, a filtration
            marker); fresh at 0 mM.
    """

    glucose: Value
    lactate: Value
    hco3: Value
    k: Value
    o2_target_kpa: Value
    ph: Value
    carrier: Literal["none", "rbc"] = "none"
    hct: Value | None = None
    tracers: Sequence[str] = ()

    def __post_init__(self) -> None:
        """Keep the tracers as a tuple and check the carrier.

        Raises:
            ValueError: for a carrier which is no choice.
        """
        object.__setattr__(self, "tracers", tuple(self.tracers))
        check_literals(self)

    @classmethod
    def default(cls) -> "Perfusate":
        """Get the perfusate of the perfusion table (`[perfusate]`)."""
        return cls(
            glucose=table_value("perfusate", "glucose"),
            lactate=table_value("perfusate", "lactate"),
            hco3=table_value("perfusate", "HCO3"),
            k=table_value("perfusate", "K"),
            o2_target_kpa=table_value("perfusate", "pO2"),
            ph=table_value("perfusate", "pH"),
        )

    def problems(self) -> list[str]:
        """Get what is wrong with the perfusate."""
        where = "perfusate"
        out: list[str] = []
        for name in ("glucose", "lactate", "hco3", "k"):
            out += _unit(where, name, getattr(self, name), U.mM)
            out += _at_least(where, name, getattr(self, name), 0.0)
        out += _unit(where, "o2_target_kpa", self.o2_target_kpa, PU.kPa)
        out += _at_least(where, "o2_target_kpa", self.o2_target_kpa, 0.0)
        out += _unit(where, "ph", self.ph, U.dimensionless)
        if not PH_RANGE[0] <= self.ph.value <= PH_RANGE[1]:
            out.append(
                f"{where}: the pH {self.ph.value} is outside of {PH_RANGE[0]} to "
                f"{PH_RANGE[1]}, the range of a perfusate"
            )
        if self.carrier == "rbc":
            out.append(
                f"{where}: the oxygen carrier 'rbc' is not built yet, only 'none'"
            )
        if self.carrier == "none" and self.hct is not None:
            out.append(f"{where}: a hematocrit needs the carrier 'rbc'")
        taken = set(ids.PERFUSATE) | set(ids.CELL)
        for t in self.tracers:
            if not SID.fullmatch(t):
                out.append(f"{where}: the tracer '{t}' is no SId")
            elif t in taken:
                out.append(f"{where}: the tracer '{t}' is a quantity of the model")
        if len(set(self.tracers)) != len(self.tracers):
            out.append(f"{where}: a tracer is given twice")
        return out


@dataclass(frozen=True)
class Circuit:
    """The circuit of a machine perfusion.

    The pump drives the perfusate from the reservoir through the inflow
    tubing, in which the membrane oxygenator exchanges its oxygen and CO2
    with the sweep gas, through the organ and the outflow tubing back into
    the reservoir; every site is well stirred.

    Attributes:
        reservoir: the volume of the perfusate in the reservoir (ml).
        tubing_in: the volume of the inflow tubing with the oxygenator (ml).
        tubing_out: the volume of the outflow tubing (ml).
        flow: the total flow of the pump (ml/min), the parameter `Q_pump`.
        oxygenator_k: the first order rate constant with which the
            oxygenator sets the oxygen of the inflow tubing towards the
            target of the perfusate (1/min).
        oxygenator_k_co2: the first order rate constant of the exchange of
            CO2 with the gas (1/min), `None` for `oxygenator_k`.
        pco2_gas: the partial pressure of CO2 of the sweep gas (kPa), `None`
            for the pCO2 of the fresh perfusate, at which the fresh
            perfusate stays at its pH.
    """

    reservoir: Value
    tubing_in: Value
    tubing_out: Value
    flow: Value
    oxygenator_k: Value
    oxygenator_k_co2: Value | None = None
    pco2_gas: Value | None = None

    @classmethod
    def default(cls) -> "Circuit":
        """Get the circuit of the perfusion table (`[circuit]`).

        The pCO2 of the gas is that of the fresh perfusate (`None`): the
        sweep gas of the study is not known.
        """
        names = ("reservoir", "tubing_in", "tubing_out", "flow", "oxygenator_k")
        return cls(
            **{name: table_value("circuit", name) for name in names},
            oxygenator_k_co2=table_value("circuit", "oxygenator_k_co2"),
        )

    def problems(self) -> list[str]:
        """Get what is wrong with the circuit."""
        where = "circuit"
        out: list[str] = []
        for name in ("reservoir", "tubing_in", "tubing_out"):
            out += _unit(where, name, getattr(self, name), PU.ml)
            out += _positive(where, name, getattr(self, name))
        out += _unit(where, "flow", self.flow, PU.ml_per_min)
        out += _at_least(where, "flow", self.flow, 0.0)
        for name in ("oxygenator_k", "oxygenator_k_co2"):
            value = getattr(self, name)
            out += _unit(where, name, value, U.per_min)
            if value is not None:
                out += _at_least(where, name, value, 0.0)
        out += _unit(where, "pco2_gas", self.pco2_gas, PU.kPa)
        if self.pco2_gas is not None:
            out += _at_least(where, "pco2_gas", self.pco2_gas, 0.0)
        for name in ("reservoir", "tubing_in", "tubing_out", "flow"):
            if getattr(self, name).uncertainties:
                out.append(
                    f"{where}: {name} has an uncertainty, which the build does "
                    "not convert into the units of the model"
                )
        return out


_SITE_TOKENS = frozenset(
    {ids.RESERVOIR, ids.INFLOW, ids.OUTFLOW, ids.PERIPORTAL, ids.PERICENTRAL}
    | {ids.URINE_SITE}
)
"""The site tokens of the circuit, which no organ token may be."""

_ORGAN_TOKEN = re.compile(r"[a-z][a-z0-9]*")


@dataclass(frozen=True)
class _Zonated(ABC):
    """The zones, the fractions and the released quantities of an organ.

    The keyword-only fields and the methods shared by the organs; each organ
    has the positional fields `mass`, `zones` and `profiles` (after its token,
    a generic organ), so `zones` follows the mass.

    An organ is a chain of zones in series with the flow; each zone is a
    vascular space and its cells with the cellular core, and has `1/n` of
    the vascular and cell volume of the organ. A profile weighs a process of
    the core per zone (`p_<process>_<zone>`), normalised so that the total
    capacity of the organ is that of one zone of the whole volume, the
    pattern of the zonated liver of the PBPK library
    (`pkpdsim.pbpk.variants.liver_zonated`).

    Attributes:
        zones: the number of zones, 1 the single compartment.
        profiles: the relative weight of a process of `core.PROCESSES` in
            every zone, the first zone first; a process without a profile is
            uniform.
        vascular_fraction: the volume fraction of the vascular space,
            `None` for the value of the table of the species.
        cell_fraction: the volume fraction of the cells, `None` for the
            table.
        density: the density of the organ (kg/l), `None` for the table.
        released: the quantities of `ids.RELEASABLE` which its necrotic
            cells release into the perfusate, `None` for those of the organ
            (`RELEASED`); each needs its content in the constants of the core.
    """

    RELEASED: ClassVar[tuple[str, ...]] = ids.RELEASABLE
    """The quantities which the organ releases by default."""

    what: ClassVar[str] = "organ"
    """The organ in the problems of a definition."""

    TABLE_FIELDS: ClassVar[tuple[str, ...]] = (
        "vascular_fraction",
        "cell_fraction",
        "density",
    )
    """The constants which the table of the organ gives where they are `None`."""

    _: KW_ONLY
    vascular_fraction: Value | None = None
    cell_fraction: Value | None = None
    density: Value | None = None
    released: Sequence[str] | None = None

    if TYPE_CHECKING:
        # positional fields of every organ, after its token: mass, zones, profiles

        @property
        def mass(self) -> Value:
            """The mass of the organ (g)."""
            ...

        @property
        def zones(self) -> int:
            """The number of zones, 1 the single compartment."""
            ...

        @property
        def profiles(self) -> Mapping[str, Sequence[float]]:
            """The relative weight of a process of `core.PROCESSES` per zone."""
            ...

    def __post_init__(self) -> None:
        """Keep the profiles as a read-only map of tuples, the released as a tuple.

        Raises:
            ValueError: for less than one zone.
        """
        if self.zones < 1:
            raise ValueError(f"A {self.what} has at least one zone, not {self.zones}.")
        profiles = {k: tuple(float(x) for x in v) for k, v in self.profiles.items()}
        object.__setattr__(self, "profiles", FrozenMap(profiles))
        if self.released is not None:
            object.__setattr__(self, "released", tuple(self.released))

    def releases(self) -> tuple[str, ...]:
        """Get the released quantities, the default of the organ for `None`."""
        return self.RELEASED if self.released is None else tuple(self.released)

    @abstractmethod
    def _table(self, species: Species) -> Mapping[str, Value]:
        """Get the fractions and the density of the table of the species.

        Raises:
            KeyError: if the table has no section of the organ and species.
        """

    def resolved(self, species: Species) -> "Self":
        """Get the organ with the constants of the table where it has none.

        Args:
            species: the species of the table.

        Returns:
            The organ with every fraction and the density given.

        Raises:
            KeyError: if a constant is missing and the table has no section
                of the organ and species.
        """
        names = [n for n in self.TABLE_FIELDS if getattr(self, n) is None]
        if not names:
            return self
        table = self._table(species)
        return dataclasses.replace(self, **{n: table[n] for n in names})

    def problems(self) -> list[str]:
        """Get what is wrong with the organ."""
        where = self.what
        out = _unit(where, "mass", self.mass, PU.g)
        out += _positive(where, "mass", self.mass)
        if self.mass.uncertainties:
            out.append(
                f"{where}: the mass has an uncertainty, which the build does not "
                "convert into the units of the model"
            )
        for name in ("vascular_fraction", "cell_fraction"):
            out += _unit(where, name, getattr(self, name), U.dimensionless)
            out += _positive(where, name, getattr(self, name))
        out += _unit(where, "density", self.density, U.kg_per_l)
        out += _positive(where, "density", self.density)
        f_vas, f_cell = self.vascular_fraction, self.cell_fraction
        if f_vas is not None and f_cell is not None and f_vas.value + f_cell.value > 1:
            out.append(f"{where}: the vascular and cell fractions exceed 1")
        for process, values in self.profiles.items():
            if process not in PROCESSES:
                out.append(
                    f"{where}: the profile '{process}' is no process of the core "
                    f"({', '.join(PROCESSES)})"
                )
            if len(values) != self.zones:
                out.append(
                    f"{where}: the profile '{process}' has {len(values)} values "
                    f"for {self.zones} zones"
                )
            if any(not (x >= 0 and math.isfinite(x)) for x in values):
                out.append(f"{where}: the profile '{process}' has a value below 0")
            elif not sum(values) > 0:
                out.append(f"{where}: the profile '{process}' is 0 in every zone")
        released = self.releases()
        for q in released:
            if q not in ids.RELEASABLE:
                out.append(
                    f"{where}: the released quantity '{q}' is none of "
                    f"{', '.join(ids.RELEASABLE)}"
                )
        if len(set(released)) != len(released):
            out.append(f"{where}: a released quantity is given twice")
        return out


@dataclass(frozen=True)
class LiverOrgan(_Zonated):
    """The liver, a chain of zones from periportal to pericentral, with bile.

    The perfusate enters the first zone (a single portal inflow, the rat
    setup; a split between the hepatic artery and the portal vein is not
    built) and leaves the last one; each zone is a vascular space (the
    sinusoids) in series with the flow and its hepatocytes with the cellular
    core, so oxygen, glucose and the damage form gradients along the
    sinusoid, and secretes bile. The constants without a value are those of
    `[liver_<species>]` of the perfusion table; the liver releases ALT, AST,
    LDH and potassium.

    Attributes:
        mass: the mass of the liver (g).
        zones: the number of zones (`_Zonated`, as the fields below).
        profiles: the zonal profiles of the processes.
        vascular_fraction: the vascular fraction, `None` for the table.
        cell_fraction: the cell fraction, `None` for the table.
        density: the density (kg/l), `None` for the table.
        released: the released quantities, `None` for `RELEASED`.
    """

    key: ClassVar[str] = "li"
    """The organ token of the ids (`M_li`)."""

    what: ClassVar[str] = "liver"

    mass: Value
    zones: int = 1
    profiles: Mapping[str, Sequence[float]] = field(default_factory=dict)

    @classmethod
    def default(
        cls,
        zones: int = 1,
        profiles: Mapping[str, Sequence[float]] | None = None,
        species: Species = Species.RAT,
    ) -> "LiverOrgan":
        """Get the reference liver of a species (`[liver_<species>]`)."""
        return cls(
            mass=table_value(f"liver_{species}", "mass"),
            zones=zones,
            profiles=profiles or {},
        ).resolved(species)

    def _table(self, species: Species) -> Mapping[str, Value]:
        """Get the constants of `[liver_<species>]`."""
        return _section(f"liver_{species}", self.TABLE_FIELDS)


@dataclass(frozen=True)
class KidneyOrgan(_Zonated):
    """The kidney, with the filtration of the perfusate and the urine.

    The perfusate enters the kidney by the renal artery and leaves by the
    renal vein; the zones are in series with the flow (one zone, the cortex,
    by default). The glomeruli filter the water and the small solutes of the
    perfusate (`ids.FILTERED`, not the enzymes) at the GFR `gfr*mass`,
    capped at the fraction `ff_max` of the pump flow (so the filtration
    follows a low flow and stops with the pump), shared by the zones by
    their volume, and the tubules reabsorb the
    fraction `reabsorption` of the filtered water and small solutes, so the
    urine has the composition of the perfusate; the tracers are filtered and
    not reabsorbed (a filtration marker). The urine leaves the perfusate, so
    the flow after the kidney is the pump flow less the urine flow, at least
    `(1 - ff_max)` of it at any flow. The kidney
    releases LDH and potassium and secretes no bile. The constants without
    a value are those of `[kidney_<species>]` of the perfusion table.

    Attributes:
        mass: the mass of the kidney (g).
        zones: the number of zones (`_Zonated`, as the fields below).
        profiles: the zonal profiles of the processes.
        vascular_fraction: the vascular fraction, `None` for the table.
        cell_fraction: the cell fraction, `None` for the table.
        density: the density (kg/l), `None` for the table.
        released: the released quantities, `None` for `RELEASED`.
        gfr: the glomerular filtration rate per g of kidney (ml/min/g),
            `None` for the table.
        reabsorption: the fraction of the filtered water and small solutes
            which the tubules reabsorb, from 0 to below 1, `None` for the
            table.
        ff_max: the largest fraction of the pump flow which the glomeruli
            filter, above 0 and at most 1, `None` for the table.
    """

    key: ClassVar[str] = "ki"
    """The organ token of the ids (`M_ki`)."""

    RELEASED: ClassVar[tuple[str, ...]] = ("ldh", "k")
    what: ClassVar[str] = "kidney"
    TABLE_FIELDS: ClassVar[tuple[str, ...]] = (
        *_Zonated.TABLE_FIELDS,
        "gfr",
        "reabsorption",
        "ff_max",
    )

    mass: Value
    zones: int = 1
    profiles: Mapping[str, Sequence[float]] = field(default_factory=dict)
    _: KW_ONLY
    gfr: Value | None = None
    reabsorption: Value | None = None
    ff_max: Value | None = None

    @classmethod
    def default(
        cls,
        zones: int = 1,
        profiles: Mapping[str, Sequence[float]] | None = None,
        species: Species = Species.RAT,
    ) -> "KidneyOrgan":
        """Get the reference kidney of a species (`[kidney_<species>]`)."""
        return cls(
            mass=table_value(f"kidney_{species}", "mass"),
            zones=zones,
            profiles=profiles or {},
        ).resolved(species)

    def _table(self, species: Species) -> Mapping[str, Value]:
        """Get the constants of `[kidney_<species>]`."""
        return _section(f"kidney_{species}", self.TABLE_FIELDS)

    def problems(self) -> list[str]:
        """Get what is wrong with the kidney."""
        where = self.what
        out = super().problems()
        out += _unit(where, "gfr", self.gfr, PU.ml_per_min_g)
        out += _positive(where, "gfr", self.gfr)
        out += _unit(where, "reabsorption", self.reabsorption, U.dimensionless)
        f = self.reabsorption
        if f is not None and not 0 <= f.value < 1:
            out.append(
                f"{where}: the reabsorption must be at least 0 and below 1, got "
                f"{f.value}"
            )
        out += _unit(where, "ff_max", self.ff_max, U.dimensionless)
        ff = self.ff_max
        if ff is not None and not 0 < ff.value <= 1:
            out.append(f"{where}: ff_max must be above 0 and at most 1, got {ff.value}")
        for name in ("gfr", "reabsorption", "ff_max"):
            value = getattr(self, name)
            if value is not None and value.uncertainties:
                out.append(
                    f"{where}: {name} has an uncertainty, which the build does "
                    "not convert into the units of the model"
                )
        return out


@dataclass(frozen=True)
class GenericOrgan(_Zonated):
    """An organ with the cellular core and no excretion (no bile, no urine).

    Another organ in a circuit is a definition and not a new model: the
    zones, the core and the release of LDH and potassium, with the constants
    of the core of the organ (`core.CoreParameters`, `[core_generic]` copies
    of the liver values by default). The constants without a value are those
    of the organ `key` of the physiology table of the species (`[male.organs.
    <key>]` of `physiology_<species>.toml`: the vascular fraction `f_vas`,
    the cell fraction `1 - f_vas - f_int` and the density), so an organ of
    the PBPK library needs its mass only.

    Attributes:
        key: the organ token of the ids (`M_<key>`), lower case letters and
            digits, an organ of the physiology table (`mu`, `he`) or another.
        mass: the mass of the organ (g).
        zones: the number of zones (`_Zonated`, as the fields below).
        profiles: the zonal profiles of the processes.
        vascular_fraction: the vascular fraction, `None` for the table.
        cell_fraction: the cell fraction, `None` for the table.
        density: the density (kg/l), `None` for the table.
        released: the released quantities, `None` for `RELEASED`.
    """

    RELEASED: ClassVar[tuple[str, ...]] = ("ldh", "k")

    key: str
    mass: Value
    zones: int = 1
    profiles: Mapping[str, Sequence[float]] = field(default_factory=dict)

    @classmethod
    def default(cls, key: str, species: Species = Species.RAT) -> "GenericOrgan":
        """Get an organ of the physiology table with its mass and fractions.

        Args:
            key: the organ of the physiology table (`mu`, `he`, `br`).
            species: the species of the table.

        Returns:
            The organ of the reference individual (male) of the species.

        Raises:
            KeyError: if the table has no organ `key`.
        """
        mass = entry(f"physiology_{species}", f"male.organs.{key}", "mass")
        assert isinstance(mass.value, float)
        return cls(
            key,
            Value(
                mass.value * 1e3,
                PU.g,
                mass.source,
                f"physiology_{species}.toml [male.organs.{key}] mass in g; "
                f"{mass.notes}",
            ),
        ).resolved(species)

    def _table(self, species: Species) -> Mapping[str, Value]:
        """Get the fractions and the density of the physiology table."""
        return _physiology(self.key, species)

    def token_problems(self) -> list[str]:
        """Get what is wrong with the token of the organ."""
        if not _ORGAN_TOKEN.fullmatch(self.key) or self.key in _SITE_TOKENS:
            return [
                f"organ: the token '{self.key}' is no organ token (lower case "
                "letters and digits, no site of the circuit)"
            ]
        if ids.is_zone(self.key):
            return [f"organ: the token '{self.key}' is a zone id"]
        return []

    def problems(self) -> list[str]:
        """Get what is wrong with the organ, its token included."""
        return self.token_problems() + super().problems()


def _section(section: str, names: Sequence[str]) -> dict[str, Value]:
    """Get numbers of a section of the perfusion table.

    Raises:
        KeyError: if the table has no such section or number.
    """
    return {k: table_value(section, k) for k in names}


def _physiology(key: str, species: Species) -> dict[str, Value]:
    """Get the fractions and the density of an organ of the physiology table.

    Raises:
        KeyError: if the table has no organ `key`.
    """
    name, section = f"physiology_{species}", f"male.organs.{key}"
    f_vas = entry(name, section, "f_vas")
    f_int = entry(name, section, "f_int")
    rho = entry(name, section, "density")
    assert isinstance(f_vas.value, float)
    assert isinstance(f_int.value, float)
    assert isinstance(rho.value, float)
    where = f"{name}.toml [{section}]"
    return {
        "vascular_fraction": Value(
            f_vas.value, U.dimensionless, f_vas.source, f"{where} f_vas; {f_vas.notes}"
        ),
        "cell_fraction": Value(
            1 - f_vas.value - f_int.value,
            U.dimensionless,
            "ModelChoice",
            f"derived: 1 - f_vas - f_int of {where}, f_int = {f_int.value}; the "
            "tissue space of the cells",
        ),
        "density": Value(
            rho.value, U.kg_per_l, rho.source, f"{where} density; {rho.notes}"
        ),
    }


type Organ = LiverOrgan | KidneyOrgan | GenericOrgan
"""The organs which a circuit mounts."""


@dataclass(frozen=True)
class PerfusionDefinition:
    """A machine perfusion: an organ mounted into a circuit.

    Attributes:
        sid: the id of the model, the name of its file.
        name: the name of the model.
        species: the species of the organ, the taxon of the model.
        organ: the mounted organ.
        circuit: the circuit.
        perfusate: the fresh perfusate.
        temperature: the temperature of the perfusion (degC), the parameter
            `T`, which an experiment changes.
        core: the constants of the cellular core of the organ.
        creators: the creators of the model (its model history).
    """

    sid: str
    name: str
    species: Species
    organ: Organ
    circuit: Circuit
    perfusate: Perfusate
    temperature: Value
    core: CoreParameters
    creators: Sequence[Creator] = ()

    def __post_init__(self) -> None:
        """Normalise the species from text and keep the creators as a tuple."""
        object.__setattr__(self, "species", Species(self.species))
        object.__setattr__(self, "creators", tuple(self.creators))

    def validate(self) -> list[str]:
        """Get what is wrong with the definition.

        Returns:
            The problems, empty for a definition which builds.
        """
        out: list[str] = []
        if not SID.fullmatch(self.sid):
            out.append(f"the id '{self.sid}' is no SId")
        out += self.circuit.problems()
        out += self.perfusate.problems()
        out += _unit("temperature", "temperature", self.temperature, PU.degC)
        if not self.temperature.value > -273.15:
            out.append("temperature: below the absolute zero")
        out += self._organ_problems()
        if not out and self.perfusate.tracers:
            out += self._tracer_clashes()
        return out

    @property
    def released(self) -> tuple[str, ...]:
        """The quantities which the necrotic cells of the organ release."""
        return self.organ.releases()

    def _organ_problems(self) -> list[str]:
        """Get what is wrong with the organ and its core."""
        organ = self.organ
        if isinstance(organ, GenericOrgan):
            token = organ.token_problems()
            if token:
                return token
        try:
            resolved = organ.resolved(self.species)
        except KeyError:
            return [
                f"{organ.what}: the table of the {self.species} has no "
                f"{organ.what} '{organ.key}' for the vascular_fraction, the "
                "cell_fraction and the density, give them"
            ]
        out = resolved.problems()
        out += self.core.problems(
            tuple(q for q in self.released if q in ids.RELEASABLE),
            bile=isinstance(organ, LiverOrgan),
        )
        return out

    def _tracer_clashes(self) -> list[str]:
        """Get the ids of the tracers which are ids of the model without them.

        A tracer `t` has the ids `t_<site>` (every site and mean), `t_fresh`
        and the flows `FLOW_t_<site>_<next>`; one of them may be an observer,
        a compartment or a constant of the model (`pO2`, `Vvas`, `V`).
        """
        from pkpdsim.perfusion.builder import compose, model_ids, tracer_ids

        plain = dataclasses.replace(
            self, perfusate=dataclasses.replace(self.perfusate, tracers=())
        )
        taken = model_ids(compose(plain))
        out: list[str] = []
        for t in self.perfusate.tracers:
            clashes = sorted(tracer_ids(self, t) & taken)
            if clashes:
                out.append(
                    f"perfusate: the ids of the tracer '{t}' are ids of the "
                    f"model: {clashes}"
                )
        return out

    def parameter_value(self, sid: str) -> float:
        """Get the number of a constant parameter of the model of the definition.

        The experiments read a number of the model here (the flow, the fresh
        perfusate of a medium change), the one place of the number.

        Args:
            sid: the id of the parameter.

        Returns:
            The number in the unit of the parameter of the model.

        Raises:
            KeyError: if the model has no constant parameter of the id with a
                number (a parameter of an initial assignment or a rule).
        """
        from pkpdsim.perfusion.builder import compose

        model = compose(self)
        assigned = {a.symbol for a in model.assignments or []}
        for p in model.parameters or []:
            if (
                p.sid == sid
                and p.constant
                and sid not in assigned
                and isinstance(p.value, int | float)
            ):
                return float(p.value)
        raise KeyError(f"The model '{self.sid}' has no constant number '{sid}'.")


def _unit(where: str, name: str, value: Value | None, unit: UnitType) -> list[str]:
    """Check that a value is given in a unit, compared by the id of the unit.

    The ids tell a temperature in `degC` from one in `kelvin`, which have
    the same definition.
    """
    if value is None or _unit_id(value.unit) == _unit_id(unit):
        return []
    return [f"{where}: {name} must be in {_unit_id(unit)}, got {_unit_id(value.unit)}"]


def _unit_id(unit: UnitType) -> str:
    """Get the id of a unit definition, the text of any other unit."""
    return str(getattr(unit, "sid", unit))


def _positive(where: str, name: str, value: Value | None) -> list[str]:
    """Check that a value is a positive number."""
    if value is None or (value.value > 0 and math.isfinite(value.value)):
        return []
    return [f"{where}: {name} must be positive, got {value.value}"]


def _at_least(where: str, name: str, value: Value, lower: float) -> list[str]:
    """Check that a value is a finite number of at least `lower`."""
    if value.value >= lower and math.isfinite(value.value):
        return []
    return [f"{where}: {name} must be at least {lower}, got {value.value}"]
