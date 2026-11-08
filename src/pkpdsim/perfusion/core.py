"""The cellular core of a zone of a perfused organ.

A zone of an organ in a machine perfusion is a vascular space, through which
the perfusate flows (the circuit, not this module), and its cells. The core
is what happens in the cells and between the cells and the vascular space:
the import of glucose and oxygen, the energy metabolism (glycogenolysis,
glycolysis, lactate dehydrogenase, TCA cycle, oxidative phosphorylation, ATP
use), the reactive oxygen species (ROS) at low oxygen and their
detoxification by glutathione, regenerated with NADPH of the pentose
phosphate pathway, the death of the cells (necrosis) and the release of
their enzymes (ALT, AST, LDH) and potassium into the vascular space, the
secretion of bile, the export of lactate with a proton and the bicarbonate
buffer and pH of the perfusate. The processes and their values are ported
from the liver model of `pkdb_models.misc.machine_perfusion` (from the PVL
liver model) with the defects of its TODO fixed: the rates are per litre of
cells (no `f_scale`), the release moves amounts, the ATP inhibition of the
glycolysis is defined at ATP 0, the glutathione is regenerated with NADPH,
and every element has its own name.

The species of the cells are concentrations in mmole per l of cells, the
species of the vascular space concentrations of the perfusate (mM); every
rate is in mmole/min, a rate per volume of cells times the volume of the
cells (`cells`) and the living fraction `f_alive_<zone>`. The ids are those
of `pkpdsim.perfusion.ids`.

The constants of the core are ids of the whole model, shared by its zones
(`core_parameters`, the numbers of the parameter database, the sections
`core_<organ>` and `core_<organ>_start` of `tables/perfusion.toml`): a fit
has one value per constant. A zonal profile weighs a process per zone (the
pattern of the zonated liver of the PBPK library), the dimensionless factor
`p_<process>_<zone>` in the rate of the process (`zone_objects(profile=)`),
which the organ normalises; a process without a profile has none. The
temperature factor `f_T` (`ids.F_T`) and the constants of the buffer
(`buffer_parameters`) are ids of the whole model as well; the circuit defines
`f_T` as the Arrhenius factor of the temperature, and the start of the
dissolved CO2 of the perfusate is the initial assignment of the buffer
equilibrium, `CO2 = h*HCO3/Ka_hco3`, so it follows the one buffer of the
model.

The rate laws, one choice each, with the reason:

- Biological rates scale with `f_T * f_alive`: the temperature acts on the
  enzymes and transporters, only the living cells metabolise.
- A Hill term reads the absolute value of its species, `(|S|/K)^n`, so it
  stays defined for the small negative values of the integrator at a
  non-integer `n`.
- Glucose, oxygen import: the reversible facilitated transport of the ported
  model, `Vmax/Km*(S_out - S_in)/(1 + S_out/Km + S_in/Km)`, saturable both ways.
- Lactate export: the monocarboxylate transporter co-transports a proton, so
  the reverse term reads the perfusate protons over the cytosolic protons,
  `(lac_in - lac_out*h_out/h_cell)/(...)`; the cytosolic pH is a constant
  (`LACEX_h_cell`, pH 7.2), not a state, and the protons of glycolysis leave
  the cell with the lactate; the reverse rate vanishes with the protons, so
  the protons stay positive.
- Glycogenolysis: Michaelis-Menten in glycogen (glucosyl units), a new store,
  since the perfusate glucose of the study rises from 6 to 11 mM.
- Glycolysis: `glc + 2 ADP + 2 NAD -> 2 pyr + 2 ATP + 2 NADH`, carbon
  balanced (the ported model made one pyruvate), Michaelis-Menten in glucose,
  ADP and NAD, inhibited by ATP with `1/(1 + (|ATP|/Ki)^n)` instead of the
  singular `(Ki/ATP)^n`; both are about 0.06 at the start (ATP 3 mM).
- LDH: `pyr + NADH <-> lac + NAD`, a reversible random bi-bi law
  `Vmax/(Km_pyr*Km_nadh)*(pyr*NADH - lac*NAD/Keq)/((1 + pyr/Km_pyr + lac/Km_lac)
  *(1 + NADH/Km_nadh + NAD/Km_nad))`, which consumes NADH (the ported model
  read it as an activator, so a reactant could become negative).
- TCA: `pyr + 5 NAD -> 5 NADH + 3 CO2`, the pyruvate dehydrogenase and the
  cycle with FADH2 counted as NADH and the GTP left out; the CO2 goes to the
  perfusate directly, it crosses the membrane within seconds.
- OXPHOS: `NADH + 0.5 O2 + 3 ADP -> NAD + 3 ATP`, the ported law and P/O ratio.
- ATP use: Michaelis-Menten in ATP, ported.
- ROS production: `NADH + O2 -> NAD + ROS` (ROS as H2O2), the electron leak of
  the respiratory chain at low oxygen; the ported step function of oxygen
  becomes `O2/(Km + O2) * 1/(1 + (O2/Ki)^n)`, smooth, 0 without oxygen (which
  it consumes) and small at normoxia.
- ROS detoxification: `ROS + 2 GSH -> GSSG` (glutathione peroxidase), the
  ported law `k*GSH*ROS/(Km + ROS)` with the stoichiometry of the enzyme.
- Glutathione reductase: `GSSG + NADPH -> 2 GSH + NADP`, the ported law with
  NADPH, which the ported model replaced by NADH.
- Pentose phosphate pathway: `glc + 12 NADP -> 12 NADPH + 6 CO2`, the net
  oxidative pathway with the recycling of the pentoses, Michaelis-Menten in
  glucose and NADP, the source of NADPH.
- Necrosis: the living fraction dies at `v_nec = NEC_k*f_T*(1 - (1 - s_atp)*
  (1 - s_h))`, either drive alone gives the full rate: the energy failure
  `s_atp = 1/(1 + (ATP/NEC_atp)^n)` (a smooth threshold, instead of the
  instantaneous switch of the ported model) and the acidosis of the perfusate
  `s_h = (h/NEC_h)^n/(1 + (h/NEC_h)^n)`. The acidosis drives the damage, and
  does not inhibit the energy metabolism, since the medium changes of the
  study protect the liver far beyond the dilution of the enzymes, at an
  outflow pO2 of 30 to 70 kPa which points at no energy deficit.
- Necrotic cells: the dying cells move their enzymes and potassium at
  `v_nec` into the dead cells of the zone, which release them first order
  (`PROTEINR_k`, `KR_k`) into the vascular space, so the totals are
  conserved and a dead zone releases its whole content; the release is a
  leak of a ruptured membrane and does not scale with `f_T`.
- Bile: the volume secreted at `BILE_Vmax*ATP/(Km + ATP)` per litre of
  living cells (the ATP-dependent export of the bile salts), a volume only.
- Bicarbonate buffer of the perfusate: `H+ + HCO3- <-> CO2` with the mass
  action `k_co2*(h*hco3/Ka - CO2)` (the dehydration of carbonic acid is fast,
  the hydration of CO2 sets the time), and the pH from the free protons,
  `pH = 3 - log10(h/1 mM)`; the CO2 of the perfusate stands for CO2 and
  H2CO3 together, no further buffer is modelled.
"""

import dataclasses
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import NamedTuple

from sbmlutils.factory import (
    AssignmentRule,
    InitialAssignment,
    Parameter,
    RateRule,
    Reaction,
    Species,
)
from sbmlutils.factory.units import UnitDefinition
from sbmlutils.metadata import BQB, SBO

from pkpdsim.database.tables import entry, load
from pkpdsim.pbpk.annotation import initial_assignment
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import FrozenMap, Value
from pkpdsim.perfusion import ids

__all__ = [
    "PROCESSES",
    "Buffer",
    "Composition",
    "CoreObjects",
    "CoreParameters",
    "buffer_objects",
    "buffer_parameters",
    "core_parameters",
    "perfusate_objects",
    "quantity_name",
    "zone_objects",
]

_TABLE = "perfusion"

_UNITS: dict[str, UnitDefinition] = {
    "mM": U.mM,
    "mmol/min/l": U.mmole_per_min_l,
    "1/min": U.per_min,
    "dimensionless": U.dimensionless,
}
"""The units of the model for the units of the perfusion table."""


def _number(section: str, key: str) -> float:
    """Get a number of the perfusion table."""
    value = entry(_TABLE, section, key).value
    if not isinstance(value, float):
        raise TypeError(f"tables/{_TABLE}.toml [{section}] {key} is no number.")
    return value


def _value(section: str, key: str) -> Value:
    """Get a number of the perfusion table with its unit, source and notes.

    Raises:
        ValueError: if the unit of the number is no unit of `_UNITS`.
    """
    e = entry(_TABLE, section, key)
    try:
        unit = _UNITS[e.unit]
    except KeyError:
        raise ValueError(
            f"tables/{_TABLE}.toml [{section}] {key}: unit '{e.unit}' is not a "
            "unit of the core."
        ) from None
    return Value(_number(section, key), unit, e.source, e.notes)


class _Spec(NamedTuple):
    """The name and SBO term of a constant of the core."""

    name: str
    sbo: str


_V, _KM, _KI = SBO.MAXIMAL_VELOCITY, SBO.MICHAELIS_CONSTANT, SBO.INHIBITORY_CONSTANT
_N, _K, _P = (
    SBO.HILL_CONSTANT,
    SBO.KINETIC_CONSTANT,
    SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
)

_SPECS: dict[str, _Spec] = {
    "GLCIM_Vmax": _Spec("Vmax glucose import", _V),
    "GLCIM_Km_glc": _Spec("Km glucose of the glucose import", _KM),
    "O2IM_Vmax": _Spec("Vmax oxygen import", _V),
    "O2IM_Km_o2": _Spec("Km oxygen of the oxygen import", _KM),
    "LACEX_Vmax": _Spec("Vmax lactate export", _V),
    "LACEX_Km_lac": _Spec("Km lactate of the lactate export", _KM),
    "LACEX_h_cell": _Spec("cytosolic protons of the lactate export", _P),
    "GLGP_Vmax": _Spec("Vmax glycogenolysis", _V),
    "GLGP_Km_glyc": _Spec("Km glycogen of the glycogenolysis", _KM),
    "GLC2PYR_Vmax": _Spec("Vmax glycolysis", _V),
    "GLC2PYR_Km_glc": _Spec("Km glucose of the glycolysis", _KM),
    "GLC2PYR_Km_adp": _Spec("Km ADP of the glycolysis", _KM),
    "GLC2PYR_Km_nad": _Spec("Km NAD of the glycolysis", _KM),
    "GLC2PYR_Ki_atp": _Spec("Ki ATP of the glycolysis", _KI),
    "GLC2PYR_n_atp": _Spec("Hill coefficient ATP of the glycolysis", _N),
    "LDH_Vmax": _Spec("Vmax lactate dehydrogenase", _V),
    "LDH_Km_pyr": _Spec("Km pyruvate of the lactate dehydrogenase", _KM),
    "LDH_Km_lac": _Spec("Km lactate of the lactate dehydrogenase", _KM),
    "LDH_Km_nadh": _Spec("Km NADH of the lactate dehydrogenase", _KM),
    "LDH_Km_nad": _Spec("Km NAD of the lactate dehydrogenase", _KM),
    "LDH_Keq": _Spec("Keq lactate dehydrogenase", SBO.EQUILIBRIUM_CONSTANT),
    "ATPU_Vmax": _Spec("Vmax ATP use", _V),
    "ATPU_Km_atp": _Spec("Km ATP of the ATP use", _KM),
    "TCA_Vmax": _Spec("Vmax TCA cycle", _V),
    "TCA_Km_pyr": _Spec("Km pyruvate of the TCA cycle", _KM),
    "TCA_Km_nad": _Spec("Km NAD of the TCA cycle", _KM),
    "OXPHOS_Vmax": _Spec("Vmax oxidative phosphorylation", _V),
    "OXPHOS_Km_nadh": _Spec("Km NADH of the oxidative phosphorylation", _KM),
    "OXPHOS_Km_o2": _Spec("Km oxygen of the oxidative phosphorylation", _KM),
    "OXPHOS_Km_adp": _Spec("Km ADP of the oxidative phosphorylation", _KM),
    "ROSP_Vmax": _Spec("Vmax ROS production", _V),
    "ROSP_Km_nadh": _Spec("Km NADH of the ROS production", _KM),
    "ROSP_Km_o2": _Spec("Km oxygen of the ROS production", _KM),
    "ROSP_Ki_o2": _Spec("Ki oxygen of the ROS production", _KI),
    "ROSP_n_o2": _Spec("Hill coefficient oxygen of the ROS production", _N),
    "ROSDETOX_k": _Spec("rate constant ROS detoxification", _K),
    "ROSDETOX_Km_ros": _Spec("Km ROS of the ROS detoxification", _KM),
    "GR_k": _Spec("rate constant glutathione reductase", _K),
    "GR_Km_nadph": _Spec("Km NADPH of the glutathione reductase", _KM),
    "PPP_Vmax": _Spec("Vmax pentose phosphate pathway", _V),
    "PPP_Km_glc": _Spec("Km glucose of the pentose phosphate pathway", _KM),
    "PPP_Km_nadp": _Spec("Km NADP of the pentose phosphate pathway", _KM),
    "NEC_k": _Spec("maximal rate of necrosis", _K),
    "NEC_atp": _Spec("ATP threshold of necrosis", _KI),
    "NEC_n_atp": _Spec("Hill coefficient ATP of necrosis", _N),
    "NEC_h": _Spec(
        "protons of the half drive of necrosis by acidosis", SBO.ACTIVATION_CONSTANT
    ),
    "NEC_n_h": _Spec("Hill coefficient protons of necrosis", _N),
    "PROTEINR_k": _Spec("rate constant of the enzyme release of necrotic cells", _K),
    "KR_k": _Spec("rate constant of the potassium release of necrotic cells", _K),
    "ALT_content": _Spec("ALT content of the cells", _P),
    "AST_content": _Spec("AST content of the cells", _P),
    "LDH_content": _Spec("LDH content of the cells", _P),
    "K_content": _Spec("potassium content of the cells", _P),
    "BILE_Vmax": _Spec("maximal bile secretion per volume of cells", _K),
    "BILE_Km_atp": _Spec("Km ATP of the bile secretion", _KM),
    "f_alive_start": _Spec("living fraction of the cells at the start", _P),
}
"""The name and SBO term of every constant of `CoreParameters` but `starts`."""

_CONTENTS = {
    "alt": "ALT_content",
    "ast": "AST_content",
    "ldh": "LDH_content",
    "k": "K_content",
}
"""The constant of the content of the cells of every released quantity."""

_RELEASE = {"alt": "PROTEINR_k", "ast": "PROTEINR_k", "ldh": "PROTEINR_k", "k": "KR_k"}
"""The release constant of every released quantity."""

_NAMES = {
    "glc": "glucose", "glyc": "glycogen", "o2": "oxygen", "pyr": "pyruvate",
    "lac": "lactate", "atp": "ATP", "adp": "ADP", "nad": "NAD", "nadh": "NADH",
    "nadp": "NADP", "nadph": "NADPH", "ros": "ROS", "gsh": "glutathione",
    "gssg": "glutathione disulfide", "co2": "CO2", "hco3": "bicarbonate",
    "h": "protons", "k": "potassium", "alt": "ALT", "ast": "AST", "ldh": "LDH",
}  # fmt: skip
"""The name of every quantity."""

_TERMS: dict[str, list[tuple[BQB, str]]] = {
    "glc": [(BQB.IS_VERSION_OF, "chebi/CHEBI:17634"), (BQB.IS, "chebi/CHEBI:4167")],
    "glyc": [(BQB.IS, "chebi/CHEBI:28087")],
    "o2": [(BQB.IS, "chebi/CHEBI:15379")],
    "pyr": [(BQB.IS, "chebi/CHEBI:15361")],
    "lac": [(BQB.IS, "chebi/CHEBI:16004")],
    "atp": [(BQB.IS, "chebi/CHEBI:15422")],
    "adp": [(BQB.IS, "chebi/CHEBI:16761")],
    "nad": [(BQB.IS, "chebi/CHEBI:44215")],
    "nadh": [(BQB.IS, "chebi/CHEBI:16908")],
    "nadp": [(BQB.IS, "chebi/CHEBI:18009")],
    "nadph": [(BQB.IS, "chebi/CHEBI:16474")],
    "ros": [(BQB.IS_VERSION_OF, "chebi/CHEBI:26523"), (BQB.IS, "chebi/CHEBI:16240")],
    "gsh": [(BQB.IS, "chebi/CHEBI:16856")],
    "gssg": [(BQB.IS, "chebi/CHEBI:17858")],
    "co2": [(BQB.IS, "chebi/CHEBI:16526")],
    "hco3": [(BQB.IS, "chebi/CHEBI:17544")],
    "h": [(BQB.IS, "chebi/CHEBI:15378")],
    "k": [(BQB.IS, "chebi/CHEBI:29103")],
    "alt": [(BQB.IS, "ncit/C25293")],
    "ast": [(BQB.IS, "ncit/C25202")],
    "ldh": [(BQB.IS, "snomedct/259319003")],
}
"""The annotations of every quantity (ROS as hydrogen peroxide equivalents)."""

_ENZYMES = ("alt", "ast", "ldh")


PROCESSES: tuple[str, ...] = (
    "GLCIM",
    "O2IM",
    "LACEX",
    "GLGP",
    "GLC2PYR",
    "LDH",
    "ATPU",
    "TCA",
    "OXPHOS",
    "ROSP",
    "ROSDETOX",
    "GR",
    "PPP",
    "BILE",
)
"""The processes of the living cells which a zonal profile may weigh."""


@dataclass(frozen=True)
class CoreParameters:
    """The constants of the cellular core, shared by the zones of an organ.

    Every field but `starts` is a constant of `_SPECS`, a parameter of the
    model under its name (`core_parameters`), in its unit: the maximal rates
    per litre of cells (mmole/min/l), the first order constants (1/min), the
    Michaelis and inhibition constants (mM), the Hill coefficients, the ATP
    threshold and rate of necrosis, the contents of the released quantities
    in mmole per litre of cells, the bile secretion in litre of bile per
    minute and litre of cells and the living fraction at the start. `starts`
    are the initial concentrations of the other species of the cells (mM) by
    quantity.
    """

    GLCIM_Vmax: Value
    GLCIM_Km_glc: Value
    O2IM_Vmax: Value
    O2IM_Km_o2: Value
    LACEX_Vmax: Value
    LACEX_Km_lac: Value
    LACEX_h_cell: Value
    GLGP_Vmax: Value
    GLGP_Km_glyc: Value
    GLC2PYR_Vmax: Value
    GLC2PYR_Km_glc: Value
    GLC2PYR_Km_adp: Value
    GLC2PYR_Km_nad: Value
    GLC2PYR_Ki_atp: Value
    GLC2PYR_n_atp: Value
    LDH_Vmax: Value
    LDH_Km_pyr: Value
    LDH_Km_lac: Value
    LDH_Km_nadh: Value
    LDH_Km_nad: Value
    LDH_Keq: Value
    ATPU_Vmax: Value
    ATPU_Km_atp: Value
    TCA_Vmax: Value
    TCA_Km_pyr: Value
    TCA_Km_nad: Value
    OXPHOS_Vmax: Value
    OXPHOS_Km_nadh: Value
    OXPHOS_Km_o2: Value
    OXPHOS_Km_adp: Value
    ROSP_Vmax: Value
    ROSP_Km_nadh: Value
    ROSP_Km_o2: Value
    ROSP_Ki_o2: Value
    ROSP_n_o2: Value
    ROSDETOX_k: Value
    ROSDETOX_Km_ros: Value
    GR_k: Value
    GR_Km_nadph: Value
    PPP_Vmax: Value
    PPP_Km_glc: Value
    PPP_Km_nadp: Value
    NEC_k: Value
    NEC_atp: Value
    NEC_n_atp: Value
    NEC_h: Value
    NEC_n_h: Value
    PROTEINR_k: Value
    KR_k: Value
    ALT_content: Value
    AST_content: Value
    LDH_content: Value
    K_content: Value
    BILE_Vmax: Value
    BILE_Km_atp: Value
    f_alive_start: Value
    starts: FrozenMap[Value]

    @classmethod
    def default(cls, organ: str = "liver_rat") -> "CoreParameters":
        """Get the constants of an organ from the parameter database.

        Args:
            organ: the organ of the sections `core_<organ>` (the constants)
                and `core_<organ>_start` (the starts of the cells) of the
                perfusion table.

        Returns:
            The constants of the core; the bile secretion is derived from
            the bile rate of the study (`[<organ>] bile_rate`, per g of
            organ) per litre of cells (`cell_fraction`, 1 kg/l) at the ATP
            of the start.
        """
        section = f"core_{organ}"
        start = f"{section}_start"
        starts = FrozenMap({q: _value(start, q) for q in load(_TABLE)[start]})
        atp0 = starts["atp"].value
        km = _value(section, "BILE_Km_atp")
        # ul/min/g -> l/min per l of organ (1 kg/l) -> per l of cells, at ATP0
        rate = _number(organ, "bile_rate") * 1e-3 / _number(organ, "cell_fraction")
        bile = Value(
            rate * (km.value + atp0) / atp0,
            U.per_min,
            "ModelChoice",
            f"derived: [{organ}] bile_rate of the perfusion table (ul/min/g) at "
            "1 kg/l over cell_fraction, divided by the ATP saturation at the "
            f"start (ATP {atp0} mM, BILE_Km_atp {km.value} mM); to be fitted",
        )
        values = {
            f.name: _value(section, f.name)
            for f in dataclasses.fields(cls)
            if f.name not in ("starts", "BILE_Vmax")
        }
        return cls(**values, BILE_Vmax=bile, starts=starts)


@dataclass(frozen=True)
class Buffer:
    """The constants of the bicarbonate buffer of the perfusate.

    Attributes:
        ka: the apparent dissociation constant `[H+]*[HCO3-]/[CO2]` (mM).
        k_hydration: the rate constant of the hydration of CO2 (1/min).
    """

    ka: Value
    k_hydration: Value

    @classmethod
    def default(cls) -> "Buffer":
        """Get the constants of the perfusion table (`acid_base`)."""
        pk = entry(_TABLE, "acid_base", "pK")
        return cls(
            ka=Value(
                10 ** (3 - _number("acid_base", "pK")),
                U.mM,
                pk.source,
                f"derived: 10^-pK M of [acid_base] pK of the perfusion table; "
                f"{pk.notes}",
            ),
            k_hydration=_value("acid_base", "k_hydration"),
        )


@dataclass(frozen=True)
class Composition:
    """The composition of a perfusate, the start of the perfusate species.

    The free protons follow from the pH; the dissolved CO2 is the initial
    assignment of the buffer equilibrium (`perfusate_objects`); a fresh
    perfusate has no enzymes.

    Attributes:
        glc: glucose (mM).
        o2: free oxygen (mM).
        lac: lactate (mM).
        hco3: bicarbonate (mM).
        k: potassium (mM).
        ph: the pH.
    """

    glc: Value
    o2: Value
    lac: Value
    hco3: Value
    k: Value
    ph: Value

    @classmethod
    def default(cls) -> "Composition":
        """Get the perfusate of the perfusion table (`perfusate`)."""
        return cls(
            glc=_value("perfusate", "glucose"),
            o2=Value(
                _number("perfusate", "pO2") * _number("oxygen", "solubility"),
                U.mM,
                "ModelChoice",
                "derived: [perfusate] pO2 times [oxygen] solubility of the "
                "perfusion table, the oxygenated perfusate at 37 C",
            ),
            lac=_value("perfusate", "lactate"),
            hco3=_value("perfusate", "HCO3"),
            k=_value("perfusate", "K"),
            ph=_value("perfusate", "pH"),
        )

    def concentrations(self) -> dict[str, float | None]:
        """Get the start of every quantity of `ids.PERFUSATE` (mM).

        Returns:
            The concentrations, `None` for the CO2, whose start is the
            buffer equilibrium.
        """
        values: dict[str, float | None] = {
            "glc": self.glc.value,
            "o2": self.o2.value,
            "lac": self.lac.value,
            "co2": None,
            "hco3": self.hco3.value,
            "h": 10 ** (3 - self.ph.value),
            "k": self.k.value,
        }
        return {q: values.get(q, 0.0) for q in ids.PERFUSATE}


@dataclass
class CoreObjects:
    """The SBML objects of a part of a perfusion model.

    Attributes:
        species: the species.
        reactions: the reactions.
        parameters: the parameters.
        rules: the assignment rules.
        rate_rules: the rate rules.
        assignments: the initial assignments.
    """

    species: list[Species] = field(default_factory=list)
    reactions: list[Reaction] = field(default_factory=list)
    parameters: list[Parameter] = field(default_factory=list)
    rules: list[AssignmentRule] = field(default_factory=list)
    rate_rules: list[RateRule] = field(default_factory=list)
    assignments: list[InitialAssignment] = field(default_factory=list)

    def __add__(self, other: "CoreObjects") -> "CoreObjects":
        """Get the objects of both."""
        return CoreObjects(
            *(
                getattr(self, f.name) + getattr(other, f.name)
                for f in dataclasses.fields(self)
            )
        )


def _species(
    sid: str,
    compartment: str,
    quantity: str,
    where: str,
    start: float | None,
    notes: str | None = None,
) -> Species:
    """Get a concentration species of a quantity."""
    return Species(
        sid,
        compartment,
        initialConcentration=start,
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        name=f"{_NAMES[quantity]} ({where})",
        sboTerm=SBO.POLYPEPTIDE_CHAIN if quantity in _ENZYMES else SBO.SIMPLE_CHEMICAL,
        annotations=list(_TERMS[quantity]),
        notes=notes,
    )


def _state(
    sid: str, value: float | None, unit: UnitDefinition, name: str, sbo: str
) -> Parameter:
    """Get a parameter which a rule sets."""
    return Parameter(sid, value, unit, constant=False, name=name, sboTerm=sbo)


def core_parameters(params: CoreParameters) -> list[Parameter]:
    """Get the constants of the core, ids of the whole model.

    Args:
        params: the constants.

    Returns:
        A constant parameter of every field but `starts`, under its name,
        with the citation of its source in the notes.

    Raises:
        TypeError: if the unit of a constant is no unit definition.
    """
    parameters: list[Parameter] = []
    for f in dataclasses.fields(params):
        if f.name == "starts":
            continue
        spec = _SPECS[f.name]
        value: Value = getattr(params, f.name)
        if not isinstance(value.unit, UnitDefinition):
            raise TypeError(f"The unit of '{f.name}' is no UnitDefinition.")
        parameters.append(
            value_parameter(f.name, value, value.unit, spec.name, spec.sbo)
        )
    return parameters


def buffer_parameters(buffer: Buffer | None = None) -> list[Parameter]:
    """Get the constants of the bicarbonate buffer, ids of the whole model.

    Args:
        buffer: the constants, `Buffer.default()` by default.

    Returns:
        The parameters `ids.BUFFER_KA` and `ids.BUFFER_K_HYDRATION`.
    """
    b = buffer or Buffer.default()
    return [
        value_parameter(
            ids.BUFFER_KA,
            b.ka,
            U.mM,
            "Ka bicarbonate buffer",
            SBO.ACID_DISSOCIATION_CONSTANT,
        ),
        value_parameter(
            ids.BUFFER_K_HYDRATION,
            b.k_hydration,
            U.per_min,
            "rate constant CO2 hydration",
            _K,
        ),
    ]


def perfusate_objects(
    site: str, compartment: str, composition: Composition | None = None
) -> CoreObjects:
    """Get the perfusate species of a site, `ids.species(quantity, site)`.

    Args:
        site: the site (a zone for its vascular space, `res`).
        compartment: the compartment of the perfusate at the site.
        composition: the start, `Composition.default()` by default.

    Returns:
        A species of every quantity of `ids.PERFUSATE`, in mM, and the
        initial assignment of the CO2 at the equilibrium of the buffer with
        the protons and the bicarbonate, `CO2 = h*HCO3/Ka_hco3`.
    """
    start = (composition or Composition.default()).concentrations()
    h, hco3, co2 = (ids.species(q, site) for q in ("h", "hco3", "co2"))
    return CoreObjects(
        species=[
            _species(ids.species(q, site), compartment, q, f"perfusate {site}", c)
            for q, c in start.items()
        ],
        assignments=[
            initial_assignment(
                co2, f"{h} * {hco3}/{ids.BUFFER_KA}", U.mM, f"CO2 (perfusate {site})"
            )
        ],
    )


def buffer_objects(site: str, compartment: str) -> CoreObjects:
    """Get the bicarbonate buffer and the pH of the perfusate at a site.

    Args:
        site: the site.
        compartment: the compartment of the perfusate at the site.

    Returns:
        The reaction `BUF_<site>`, `H+ + HCO3- <-> CO2`, and the pH
        `pH_<site>` with its assignment rule.
    """
    h, hco3, co2 = (ids.species(q, site) for q in ("h", "hco3", "co2"))
    ph = ids.ph(site)
    return CoreObjects(
        reactions=[
            Reaction(
                ids.reaction("BUF", site),
                f"{h} + {hco3} <-> {co2}",
                formula=(
                    f"{compartment} * {ids.BUFFER_K_HYDRATION} * "
                    f"({h} * {hco3}/{ids.BUFFER_KA} - {co2})",
                    U.mmole_per_min,
                ),
                compartment=compartment,
                name=f"bicarbonate buffer ({site})",
                sboTerm=SBO.BIOCHEMICAL_REACTION,
                notes="H+ + HCO3- <-> CO2 + H2O, at the rate of the hydration of CO2.",
            )
        ],
        parameters=[_state(ph, None, U.dimensionless, f"pH ({site})", SBO.PH)],
        rules=[
            AssignmentRule(
                ph,
                f"3 dimensionless - log10({h}/(1 mM))",
                U.dimensionless,
                name=f"pH ({site})",
                sboTerm=SBO.MATHEMATICAL_EXPRESSION,
            )
        ],
    )


def zone_objects(
    zone: str,
    cells: str,
    vascular: str,
    params: CoreParameters,
    *,
    profile: Mapping[str, float] | None = None,
    composition: Composition | None = None,
    bile: bool = True,
) -> CoreObjects:
    """Get the cellular core of a zone, its vascular space included.

    Args:
        zone: the zone id (`ids.zone_id`).
        cells: the compartment of the cells of the zone.
        vascular: the compartment of the vascular space of the zone.
        params: the constants of the core, of which the zone reads the starts
            of its cells; the constants themselves are `core_parameters`.
        profile: the factor of a process of `PROCESSES` in this zone, the
            constant `p_<process>_<zone>` in its rate (normalised by the
            organ); a process without a factor has none.
        composition: the start of the perfusate species of the vascular space.
        bile: the zone secretes bile (a liver), `Vbile_<zone>`.

    Returns:
        The species of the cells, the necrotic cells and the vascular space,
        the reactions, the factors of the profile, the rules of the pH and of
        the rate of necrosis, the rate rules of the living fraction and of
        the bile, the initial assignments of the contents of the cells, of
        the living fraction and of the CO2 of the perfusate. The rate laws
        read the constants of the core, `f_T` and the buffer constants, which
        the model defines once.

    Raises:
        ValueError: if `zone` is no zone id, or the profile has a process
            which is not in `PROCESSES` or a factor which is negative.
    """
    if not ids.is_zone(zone):
        raise ValueError(f"'{zone}' is no zone id (z1, z2, ...).")
    profile = dict(profile or {})
    unknown = sorted(set(profile) - set(PROCESSES))
    if unknown:
        raise ValueError(f"No processes of the core: {unknown}, see PROCESSES.")
    negative = sorted(k for k, x in profile.items() if not x >= 0)
    if negative:
        raise ValueError(f"The profile of {negative} in {zone} is not >= 0.")
    cell, dead = ids.cell_site(zone), ids.dead_site(zone)
    c = _ids(ids.species, cell)
    v = _ids(ids.species, zone)
    alive, nec = ids.alive(zone), ids.necrosis_rate(zone)
    base = f"{ids.F_T} * {alive} * {cells}"

    def live(process: str) -> str:
        """The scale of a process: temperature, living cells, profile."""
        if process in profile:
            return f"{base} * {ids.profile(process, zone)}"
        return base

    out = CoreObjects(
        parameters=[
            Parameter(
                ids.profile(process, zone),
                factor,
                U.dimensionless,
                constant=True,
                name=f"zonal factor of {process} ({zone})",
                sboTerm=_P,
            )
            for process, factor in profile.items()
        ]
    )
    for q in ids.CELL:
        if q in _CONTENTS:
            content = _CONTENTS[q]
            value = getattr(params, content).value
            out.species.append(_species(c(q), cells, q, f"cells {zone}", value))
            out.assignments.append(
                initial_assignment(c(q), content, U.mM, f"{_NAMES[q]} (cells {zone})")
            )
        else:
            s = params.starts[q]
            out.species.append(
                _species(c(q), cells, q, f"cells {zone}", s.value, s.notes)
            )
    for q in ids.RELEASED:
        out.species.append(
            _species(ids.species(q, dead), cells, q, f"necrotic cells {zone}", 0.0)
        )
    out += perfusate_objects(zone, vascular, composition)
    out += buffer_objects(zone, vascular)

    def mm(s: str, km: str) -> str:
        return f"{s}/({km} + {s})"

    def inhibition(s: str, ki: str, n: str) -> str:
        return f"1 dimensionless/(1 dimensionless + power(abs({s})/{ki}, {n}))"

    def exchange(process: str, km: str, source: str, sink: str, back: str) -> str:
        """The reversible facilitated transport from `source` to `sink`."""
        return (
            f"{live(process)} * {process}_Vmax/{km} * ({source} - {back})"
            f"/(1 dimensionless + {source}/{km} + {sink}/{km})"
        )

    t = SBO.TRANSPORT_REACTION
    b = SBO.BIOCHEMICAL_REACTION
    reactions: list[tuple[str, str, str, str, str]] = [
        ("GLCIM", "glucose import", f"{v('glc')} <-> {c('glc')}",
         exchange("GLCIM", "GLCIM_Km_glc", v("glc"), c("glc"), c("glc")), t),
        ("O2IM", "oxygen import", f"{v('o2')} <-> {c('o2')}",
         exchange("O2IM", "O2IM_Km_o2", v("o2"), c("o2"), c("o2")), t),
        ("LACEX", "lactate export with a proton", f"{c('lac')} <-> {v('lac')} + {v('h')}",
         exchange("LACEX", "LACEX_Km_lac", c("lac"), v("lac"),
                  f"{v('lac')} * {v('h')}/LACEX_h_cell"), t),
        ("GLGP", "glycogenolysis", f"{c('glyc')} -> {c('glc')}",
         f"{live('GLGP')} * GLGP_Vmax * {mm(c('glyc'), 'GLGP_Km_glyc')}", b),
        ("GLC2PYR", "glycolysis",
         f"{c('glc')} + 2 {c('adp')} + 2 {c('nad')} -> 2 {c('pyr')} + 2 {c('atp')} + 2 {c('nadh')}",
         f"{live('GLC2PYR')} * GLC2PYR_Vmax * {mm(c('glc'), 'GLC2PYR_Km_glc')}"
         f" * {mm(c('adp'), 'GLC2PYR_Km_adp')} * {mm(c('nad'), 'GLC2PYR_Km_nad')}"
         f" * {inhibition(c('atp'), 'GLC2PYR_Ki_atp', 'GLC2PYR_n_atp')}", b),
        ("LDH", "lactate dehydrogenase", f"{c('pyr')} + {c('nadh')} <-> {c('lac')} + {c('nad')}",
         f"{live('LDH')} * LDH_Vmax/(LDH_Km_pyr * LDH_Km_nadh)"
         f" * ({c('pyr')} * {c('nadh')} - {c('lac')} * {c('nad')}/LDH_Keq)"
         f"/((1 dimensionless + {c('pyr')}/LDH_Km_pyr + {c('lac')}/LDH_Km_lac)"
         f" * (1 dimensionless + {c('nadh')}/LDH_Km_nadh + {c('nad')}/LDH_Km_nad))", b),
        ("ATPU", "ATP use", f"{c('atp')} -> {c('adp')}",
         f"{live('ATPU')} * ATPU_Vmax * {mm(c('atp'), 'ATPU_Km_atp')}", b),
        ("TCA", "TCA cycle", f"{c('pyr')} + 5 {c('nad')} -> 5 {c('nadh')} + 3 {v('co2')}",
         f"{live('TCA')} * TCA_Vmax * {mm(c('pyr'), 'TCA_Km_pyr')} * {mm(c('nad'), 'TCA_Km_nad')}", b),
        ("OXPHOS", "oxidative phosphorylation",
         f"{c('nadh')} + 0.5 {c('o2')} + 3 {c('adp')} -> {c('nad')} + 3 {c('atp')}",
         f"{live('OXPHOS')} * OXPHOS_Vmax * {mm(c('nadh'), 'OXPHOS_Km_nadh')}"
         f" * {mm(c('o2'), 'OXPHOS_Km_o2')} * {mm(c('adp'), 'OXPHOS_Km_adp')}", b),
        ("ROSP", "ROS production", f"{c('nadh')} + {c('o2')} -> {c('nad')} + {c('ros')}",
         f"{live('ROSP')} * ROSP_Vmax * {mm(c('nadh'), 'ROSP_Km_nadh')} * {mm(c('o2'), 'ROSP_Km_o2')}"
         f" * {inhibition(c('o2'), 'ROSP_Ki_o2', 'ROSP_n_o2')}", b),
        ("ROSDETOX", "ROS detoxification", f"{c('ros')} + 2 {c('gsh')} -> {c('gssg')}",
         f"{live('ROSDETOX')} * ROSDETOX_k * {c('gsh')} * {mm(c('ros'), 'ROSDETOX_Km_ros')}", b),
        ("GR", "glutathione reductase", f"{c('gssg')} + {c('nadph')} -> 2 {c('gsh')} + {c('nadp')}",
         f"{live('GR')} * GR_k * {c('gssg')} * {mm(c('nadph'), 'GR_Km_nadph')}", b),
        ("PPP", "pentose phosphate pathway",
         f"{c('glc')} + 12 {c('nadp')} -> 12 {c('nadph')} + 6 {v('co2')}",
         f"{live('PPP')} * PPP_Vmax * {mm(c('glc'), 'PPP_Km_glc')} * {mm(c('nadp'), 'PPP_Km_nadp')}", b),
    ]  # fmt: skip
    for q in ids.RELEASED:
        name, d = _NAMES[q], ids.species(q, dead)
        reactions += [
            (f"{q.upper()}D", f"{name} of dying cells", f"{c(q)} -> {d}",
             f"{nec} * {c(q)} * {cells}", SBO.PROCESS),
            (f"{q.upper()}R", f"{name} release of necrotic cells", f"{d} -> {v(q)}",
             f"{_RELEASE[q]} * {d} * {cells}", t),
        ]  # fmt: skip
    for sid, name, equation, rate, sbo in reactions:
        out.reactions.append(
            Reaction(
                ids.reaction(sid, zone),
                equation,
                formula=(rate, U.mmole_per_min),
                compartment=cells,
                name=f"{name} ({zone})",
                sboTerm=sbo,
            )
        )

    atp, h = c("atp"), v("h")
    power_atp = f"power(abs({atp})/NEC_atp, NEC_n_atp)"
    power_h = f"power(abs({h})/NEC_h, NEC_n_h)"
    out.parameters += [
        _state(nec, 0.0, U.per_min, f"rate of necrosis ({zone})", _K),
        _state(
            alive,
            params.f_alive_start.value,
            U.dimensionless,
            f"living fraction of the cells ({zone})",
            _P,
        ),
    ]
    out.assignments.append(
        initial_assignment(
            alive, "f_alive_start", U.dimensionless, f"living fraction ({zone})"
        )
    )
    out.rules.append(
        AssignmentRule(
            nec,
            f"NEC_k * {ids.F_T} * (1 dimensionless - {power_atp}"
            f"/((1 dimensionless + {power_atp}) * (1 dimensionless + {power_h})))",
            U.per_min,
            name=f"rate of necrosis ({zone})",
            sboTerm=SBO.MATHEMATICAL_EXPRESSION,
            notes="Either the energy failure or the acidosis of the perfusate drives it.",
        )
    )
    out.rate_rules.append(
        RateRule(
            alive,
            f"-{nec} * {alive}",
            U.per_min,
            name=f"necrosis ({zone})",
            sboTerm=SBO.MATHEMATICAL_EXPRESSION,
        )
    )
    if bile:
        vbile = ids.bile(zone)
        out.parameters.append(_state(vbile, 0.0, U.litre, f"bile volume ({zone})", _P))
        out.rate_rules.append(
            RateRule(
                vbile,
                f"{live('BILE')} * BILE_Vmax * {mm(atp, 'BILE_Km_atp')}",
                U.l_per_min,
                name=f"bile secretion ({zone})",
                sboTerm=SBO.MATHEMATICAL_EXPRESSION,
            )
        )
    return out


def quantity_name(quantity: str) -> str:
    """Get the name of a quantity of the core (`glucose` of `glc`).

    Raises:
        KeyError: for a token which is no quantity of the core.
    """
    return _NAMES[quantity]


def _ids(make: Callable[[str, str], str], suffix: str) -> Callable[[str], str]:
    """Get the id function of a zone or site."""
    return lambda name: make(name, suffix)
