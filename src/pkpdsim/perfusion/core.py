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
of `pkpdsim.perfusion.ids`. The temperature factor `f_T` (`ids.F_T`) and the
constants of the buffer (`buffer_parameters`) are ids of the whole model,
defined once by the caller; the circuit defines `f_T` as the Arrhenius factor
of the temperature.

The rate laws, one choice each, with the reason:

- Biological rates scale with `f_T * f_alive`: the temperature acts on the
  enzymes and transporters, only the living cells metabolise.
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
  ADP and NAD, inhibited by ATP with `1/(1 + (ATP/Ki)^n)` instead of the
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
  ported law `Vmax*GSH*ROS/(Km + ROS)` with the stoichiometry of the enzyme.
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
from collections.abc import Callable
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

from pkpdsim.database.tables import entry
from pkpdsim.pbpk.annotation import initial_assignment
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import FrozenMap, Value
from pkpdsim.perfusion import ids

__all__ = [
    "Buffer",
    "Composition",
    "CoreObjects",
    "CoreParameters",
    "buffer_objects",
    "buffer_parameters",
    "perfusate_species",
    "zone_objects",
]

PORTED = "PVL liver model, pkdb_models misc/machine_perfusion, to be fitted"
"""The note of a value of the ported liver model."""

_TABLE = "perfusion"


def _number(section: str, key: str) -> float:
    """Get a number of the perfusion table."""
    value = entry(_TABLE, section, key).value
    if not isinstance(value, float):
        raise TypeError(f"tables/{_TABLE}.toml [{section}] {key} is no number.")
    return value


def _table(section: str, key: str, unit: UnitDefinition) -> Value:
    """Get a number of the perfusion table as a value in `unit`."""
    e = entry(_TABLE, section, key)
    return Value(_number(section, key), unit, e.source, e.notes)


def _ported(value: float, unit: UnitDefinition, note: str | None = None) -> Value:
    """Get a value of the ported liver model, with a remark on a change."""
    return Value(value, unit, "ModelChoice", PORTED + (f"; {note}" if note else ""))


def _new(value: float, unit: UnitDefinition, note: str) -> Value:
    """Get an assumed value of a process which the ported model lacks."""
    return Value(value, unit, "ModelChoice", f"{note}; to be fitted")


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
    "ROSDETOX_Vmax": _Spec("rate constant ROS detoxification", _K),
    "ROSDETOX_Km_ros": _Spec("Km ROS of the ROS detoxification", _KM),
    "GR_Vmax": _Spec("rate constant glutathione reductase", _K),
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


@dataclass(frozen=True)
class CoreParameters:
    """The constants of the cellular core of a zone.

    Every field but `starts` is a constant of `_SPECS`, a parameter
    `<field>_<zone>` of the zone (`ids.parameter`), in its unit: the
    maximal rates per litre of cells (mmole/min/l, `Vmax` of a mass action
    law in 1/min), the Michaelis and inhibition constants (mM), the Hill
    coefficients, the ATP threshold and rate of necrosis, the contents of the
    released quantities in mmole per litre of cells, the bile secretion in
    litre of bile per minute and litre of cells and the living fraction at
    the start. `starts` are the initial concentrations of the other species
    of the cells (mM) by quantity. A zonal profile is a `CoreParameters` per
    zone (`dataclasses.replace`).
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
    ROSDETOX_Vmax: Value
    ROSDETOX_Km_ros: Value
    GR_Vmax: Value
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
    starts: FrozenMap[Value] = field(default_factory=lambda: _starts())

    @classmethod
    def default(cls) -> "CoreParameters":
        """Get the ported values of the liver model and the assumed new ones.

        Returns:
            The constants of the core; the bile secretion is the bile rate of
            the study (`liver_rat` of the perfusion table) per litre of cells
            at the ATP of the start.
        """
        mm, d, per_min, vmax = U.mM, U.dimensionless, U.per_min, U.mmole_per_min_l
        starts = _starts()
        atp0 = starts["atp"].value
        km_bile = 0.1
        cell_fraction = _number("liver_rat", "cell_fraction")
        # ul/min/g -> l/min per l of liver (1 kg/l) -> per l of cells, at ATP0
        bile = _number("liver_rat", "bile_rate") * 1e-3 / cell_fraction
        bile_vmax = bile * (km_bile + atp0) / atp0
        return cls(
            GLCIM_Vmax=_ported(100.0, vmax),
            GLCIM_Km_glc=_ported(0.1, mm),
            O2IM_Vmax=_ported(100.0, vmax),
            O2IM_Km_o2=_ported(0.001, mm),
            LACEX_Vmax=_ported(100.0, vmax),
            LACEX_Km_lac=_ported(0.5, mm),
            LACEX_h_cell=_new(
                6.310e-5,
                mm,
                "the cytosolic protons of a hepatocyte at pH 7.2, held constant",
            ),
            GLGP_Vmax=_new(
                0.2, vmax, "glycogenolysis of a fed or stressed donor liver"
            ),
            GLGP_Km_glyc=_new(10.0, mm, "glycogen in glucosyl units"),
            GLC2PYR_Vmax=_ported(0.2, vmax),
            GLC2PYR_Km_glc=_ported(0.5, mm),
            GLC2PYR_Km_adp=_ported(0.01, mm),
            GLC2PYR_Km_nad=_new(0.01, mm, "the glycolysis reduces NAD (GAPDH)"),
            GLC2PYR_Ki_atp=_ported(1.5, mm),
            GLC2PYR_n_atp=_ported(4.0, d),
            LDH_Vmax=_ported(1.0, vmax),
            LDH_Km_pyr=_ported(1.0, mm),
            LDH_Km_lac=_ported(1.0, mm),
            LDH_Km_nadh=_ported(0.1, mm, "the activation constant of NADH, a Km now"),
            LDH_Km_nad=_new(0.1, mm, "the reverse LDH reduces NAD"),
            LDH_Keq=_ported(1.5, d, "of lac*NAD/(pyr*NADH) now, of lac/pyr before"),
            ATPU_Vmax=_ported(0.1, vmax),
            ATPU_Km_atp=_ported(0.01, mm),
            TCA_Vmax=_ported(0.01, vmax),
            TCA_Km_pyr=_ported(0.1, mm),
            TCA_Km_nad=_ported(0.01, mm),
            OXPHOS_Vmax=_ported(0.03, vmax),
            OXPHOS_Km_nadh=_ported(0.01, mm),
            OXPHOS_Km_o2=_ported(0.001, mm),
            OXPHOS_Km_adp=_ported(0.01, mm),
            ROSP_Vmax=_ported(0.005, vmax),
            ROSP_Km_nadh=_ported(0.01, mm, "OXPHOS_Km_nadh before"),
            ROSP_Km_o2=_new(0.001, mm, "the ROS production consumes oxygen"),
            ROSP_Ki_o2=_ported(0.01, mm, "ros_o2_threshold of a step function before"),
            ROSP_n_o2=_new(4.0, d, "the steepness of the smooth oxygen threshold"),
            ROSDETOX_Vmax=_ported(0.005, per_min),
            ROSDETOX_Km_ros=_ported(0.5, mm),
            GR_Vmax=_ported(0.005, per_min, "GSH_Vmax before"),
            GR_Km_nadph=_ported(0.1, mm, "GSH_Km_nadh of NADH before"),
            PPP_Vmax=_new(0.005, vmax, "the source of NADPH, as fast as its use"),
            PPP_Km_glc=_new(0.5, mm, "as the glycolysis"),
            PPP_Km_nadp=_new(0.01, mm, "saturated at the NADP of the start"),
            NEC_k=_new(0.02, per_min, "a dying zone dies within about an hour"),
            NEC_atp=_ported(
                0.1, mm, "atp_necrosis_threshold of a step function before"
            ),
            NEC_n_atp=_new(4.0, d, "the steepness of the smooth ATP threshold"),
            NEC_h=_new(
                1e-4, mm, "the protons of pH 7.0 (the study: 6.9 to 7.1 at 5 h)"
            ),
            NEC_n_h=_new(
                8.0, d, "steep, so the fresh perfusate (pH 7.4) drives no damage"
            ),
            PROTEINR_k=_ported(0.2, per_min, "in 1/min, 1/min/l before"),
            KR_k=_new(
                1.0, per_min, "a small ion leaves a ruptured cell faster than an enzyme"
            ),
            ALT_content=_ported(1.0 / 1.5, mm, "1.0 mmole in the 1.5 l liver before"),
            AST_content=_ported(0.7 / 1.5, mm, "0.7 mmole in the 1.5 l liver before"),
            LDH_content=_ported(0.2 / 1.5, mm, "0.2 mmole in the 1.5 l liver before"),
            K_content=_new(
                140.0, mm, "the intracellular potassium of a mammalian cell"
            ),
            BILE_Vmax=Value(
                bile_vmax,
                per_min,
                "ModelChoice",
                "derived: bile_rate of [liver_rat] of the perfusion table (ul/min/g) "
                "at 1 kg/l over cell_fraction, divided by the ATP saturation at the "
                f"start ({atp0} mM, Km {km_bile} mM); to be fitted",
            ),
            BILE_Km_atp=_new(km_bile, mm, "the ATP-dependent export of the bile salts"),
            f_alive_start=_new(
                1.0, d, "a zone starts alive; the warm ischemia may lower it"
            ),
            starts=starts,
        )


def _starts() -> FrozenMap[Value]:
    """Get the initial concentrations of the cells but the released quantities."""
    mm = U.mM
    return FrozenMap(
        {
            "glc": _ported(5.0, mm),
            "glyc": _new(
                200.0,
                mm,
                "glucosyl units of the glycogen of a fed rat liver (about 300 mM) after the ischemia",
            ),
            "o2": _ported(0.085, mm),
            "pyr": _ported(0.8, mm),
            "lac": _ported(1.0, mm),
            "atp": _ported(3.0, mm),
            "adp": _ported(1.0, mm),
            "nad": _ported(0.1, mm),
            "nadh": _ported(3.9, mm),
            "nadp": _new(0.03, mm, "the NADP of a liver, mostly reduced"),
            "nadph": _new(0.3, mm, "the NADPH of a liver, mostly reduced"),
            "ros": _ported(0.0, mm),
            "gsh": _ported(1.0, mm),
            "gssg": _ported(0.0, mm),
        }
    )


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
                f"derived: 10^-pK M of [acid_base] pK of the perfusion table; {pk.notes}",
            ),
            k_hydration=_table("acid_base", "k_hydration", U.per_min),
        )


@dataclass(frozen=True)
class Composition:
    """The composition of a perfusate, the start of the perfusate species.

    The free protons and the dissolved CO2 follow from the pH and the
    bicarbonate by the buffer; a fresh perfusate has no enzymes.

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
        mm = U.mM
        return cls(
            glc=_table("perfusate", "glucose", mm),
            o2=_table("perfusate", "O2", mm),
            lac=_table("perfusate", "lactate", mm),
            hco3=_table("perfusate", "HCO3", mm),
            k=_table("perfusate", "K", mm),
            ph=_table("perfusate", "pH", U.dimensionless),
        )

    def concentrations(self, buffer: Buffer) -> dict[str, float]:
        """Get the concentration of every quantity of `ids.PERFUSATE` (mM)."""
        h = 10 ** (3 - self.ph.value)
        values = {
            "glc": self.glc.value,
            "o2": self.o2.value,
            "lac": self.lac.value,
            "co2": h * self.hco3.value / buffer.ka.value,
            "hco3": self.hco3.value,
            "h": h,
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
    start: float,
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
    sid: str, value: float, unit: UnitDefinition, name: str, sbo: str
) -> Parameter:
    """Get a parameter which a rule sets."""
    return Parameter(sid, value, unit, constant=False, name=name, sboTerm=sbo)


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


def perfusate_species(
    site: str,
    compartment: str,
    composition: Composition | None = None,
    buffer: Buffer | None = None,
) -> list[Species]:
    """Get the perfusate species of a site, `ids.species(quantity, site)`.

    Args:
        site: the site (a zone for its vascular space, `res`).
        compartment: the compartment of the perfusate at the site.
        composition: the start, `Composition.default()` by default.
        buffer: the buffer which gives the start of the protons and CO2.

    Returns:
        A species of every quantity of `ids.PERFUSATE`, in mM.
    """
    start = (composition or Composition.default()).concentrations(
        buffer or Buffer.default()
    )
    return [
        _species(ids.species(q, site), compartment, q, f"perfusate {site}", start[q])
        for q in ids.PERFUSATE
    ]


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
        parameters=[_state(ph, 7.4, U.dimensionless, f"pH ({site})", SBO.PH)],
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
    composition: Composition | None = None,
    buffer: Buffer | None = None,
    bile: bool = True,
) -> CoreObjects:
    """Get the cellular core of a zone, its vascular space included.

    Args:
        zone: the zone id (`ids.zone_id`).
        cells: the compartment of the cells of the zone.
        vascular: the compartment of the vascular space of the zone.
        params: the constants of the core of the zone.
        composition: the start of the perfusate species of the vascular space.
        buffer: the buffer which gives the start of the protons and CO2; the
            model has its constants (`buffer_parameters`) once.
        bile: the zone secretes bile (a liver), `Vbile_<zone>`.

    Returns:
        The species of the cells, the necrotic cells and the vascular space,
        the reactions, the constants `<name>_<zone>`, the rules of the pH and
        of the rate of necrosis, the rate rules of the living fraction and of
        the bile, the initial assignments of the contents of the cells. The
        rate laws read `f_T` and the buffer constants, which the model defines.

    Raises:
        ValueError: if `zone` is no zone id.
    """
    if not ids.is_zone(zone):
        raise ValueError(f"'{zone}' is no zone id (z1, z2, ...).")
    p = _ids(ids.parameter, zone)
    cell, dead = ids.cell_site(zone), ids.dead_site(zone)
    c = _ids(ids.species, cell)
    v = _ids(ids.species, zone)
    alive, nec = ids.alive(zone), ids.necrosis_rate(zone)
    live = f"{ids.F_T} * {alive} * {cells}"

    out = CoreObjects()
    for f in dataclasses.fields(params):
        if f.name == "starts":
            continue
        spec = _SPECS[f.name]
        value: Value = getattr(params, f.name)
        if not isinstance(value.unit, UnitDefinition):
            raise TypeError(f"The unit of '{f.name}' is no UnitDefinition.")
        out.parameters.append(
            value_parameter(p(f.name), value, value.unit, spec.name, spec.sbo)
        )
    for q in ids.CELL:
        if q in _CONTENTS:
            content = p(_CONTENTS[q])
            out.species.append(
                _species(
                    c(q), cells, q, f"cells {zone}", getattr(params, _CONTENTS[q]).value
                )
            )
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
    out += CoreObjects(species=perfusate_species(zone, vascular, composition, buffer))
    out += buffer_objects(zone, vascular)

    def mm(s: str, km: str) -> str:
        return f"{s}/({p(km)} + {s})"

    def inhibition(s: str, ki: str, n: str) -> str:
        return f"1 dimensionless/(1 dimensionless + power({s}/{p(ki)}, {p(n)}))"

    def exchange(process: str, km: str, source: str, sink: str, back: str) -> str:
        """The reversible facilitated transport from `source` to `sink`."""
        k = p(km)
        return (
            f"{live} * {p(process + '_Vmax')}/{k} * ({source} - {back})"
            f"/(1 dimensionless + {source}/{k} + {sink}/{k})"
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
                  f"{v('lac')} * {v('h')}/{p('LACEX_h_cell')}"), t),
        ("GLGP", "glycogenolysis", f"{c('glyc')} -> {c('glc')}",
         f"{live} * {p('GLGP_Vmax')} * {mm(c('glyc'), 'GLGP_Km_glyc')}", b),
        ("GLC2PYR", "glycolysis",
         f"{c('glc')} + 2 {c('adp')} + 2 {c('nad')} -> 2 {c('pyr')} + 2 {c('atp')} + 2 {c('nadh')}",
         f"{live} * {p('GLC2PYR_Vmax')} * {mm(c('glc'), 'GLC2PYR_Km_glc')}"
         f" * {mm(c('adp'), 'GLC2PYR_Km_adp')} * {mm(c('nad'), 'GLC2PYR_Km_nad')}"
         f" * {inhibition(c('atp'), 'GLC2PYR_Ki_atp', 'GLC2PYR_n_atp')}", b),
        ("LDH", "lactate dehydrogenase", f"{c('pyr')} + {c('nadh')} <-> {c('lac')} + {c('nad')}",
         f"{live} * {p('LDH_Vmax')}/({p('LDH_Km_pyr')} * {p('LDH_Km_nadh')})"
         f" * ({c('pyr')} * {c('nadh')} - {c('lac')} * {c('nad')}/{p('LDH_Keq')})"
         f"/((1 dimensionless + {c('pyr')}/{p('LDH_Km_pyr')} + {c('lac')}/{p('LDH_Km_lac')})"
         f" * (1 dimensionless + {c('nadh')}/{p('LDH_Km_nadh')} + {c('nad')}/{p('LDH_Km_nad')}))", b),
        ("ATPU", "ATP use", f"{c('atp')} -> {c('adp')}",
         f"{live} * {p('ATPU_Vmax')} * {mm(c('atp'), 'ATPU_Km_atp')}", b),
        ("TCA", "TCA cycle", f"{c('pyr')} + 5 {c('nad')} -> 5 {c('nadh')} + 3 {v('co2')}",
         f"{live} * {p('TCA_Vmax')} * {mm(c('pyr'), 'TCA_Km_pyr')} * {mm(c('nad'), 'TCA_Km_nad')}", b),
        ("OXPHOS", "oxidative phosphorylation",
         f"{c('nadh')} + 0.5 {c('o2')} + 3 {c('adp')} -> {c('nad')} + 3 {c('atp')}",
         f"{live} * {p('OXPHOS_Vmax')} * {mm(c('nadh'), 'OXPHOS_Km_nadh')}"
         f" * {mm(c('o2'), 'OXPHOS_Km_o2')} * {mm(c('adp'), 'OXPHOS_Km_adp')}", b),
        ("ROSP", "ROS production", f"{c('nadh')} + {c('o2')} -> {c('nad')} + {c('ros')}",
         f"{live} * {p('ROSP_Vmax')} * {mm(c('nadh'), 'ROSP_Km_nadh')} * {mm(c('o2'), 'ROSP_Km_o2')}"
         f" * {inhibition(c('o2'), 'ROSP_Ki_o2', 'ROSP_n_o2')}", b),
        ("ROSDETOX", "ROS detoxification", f"{c('ros')} + 2 {c('gsh')} -> {c('gssg')}",
         f"{live} * {p('ROSDETOX_Vmax')} * {c('gsh')} * {mm(c('ros'), 'ROSDETOX_Km_ros')}", b),
        ("GR", "glutathione reductase", f"{c('gssg')} + {c('nadph')} -> 2 {c('gsh')} + {c('nadp')}",
         f"{live} * {p('GR_Vmax')} * {c('gssg')} * {mm(c('nadph'), 'GR_Km_nadph')}", b),
        ("PPP", "pentose phosphate pathway",
         f"{c('glc')} + 12 {c('nadp')} -> 12 {c('nadph')} + 6 {v('co2')}",
         f"{live} * {p('PPP_Vmax')} * {mm(c('glc'), 'PPP_Km_glc')} * {mm(c('nadp'), 'PPP_Km_nadp')}", b),
    ]  # fmt: skip
    for q in ids.RELEASED:
        name, d = _NAMES[q], ids.species(q, dead)
        reactions += [
            (f"{q.upper()}D", f"{name} of dying cells", f"{c(q)} -> {d}",
             f"{nec} * {c(q)} * {cells}", SBO.PROCESS),
            (f"{q.upper()}R", f"{name} release of necrotic cells", f"{d} -> {v(q)}",
             f"{p(_RELEASE[q])} * {d} * {cells}", t),
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
    power_atp = f"power({atp}/{p('NEC_atp')}, {p('NEC_n_atp')})"
    power_h = f"power({h}/{p('NEC_h')}, {p('NEC_n_h')})"
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
            alive, p("f_alive_start"), U.dimensionless, f"living fraction ({zone})"
        )
    )
    out.rules.append(
        AssignmentRule(
            nec,
            f"{p('NEC_k')} * {ids.F_T} * (1 dimensionless - {power_atp}"
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
                f"{live} * {p('BILE_Vmax')} * {mm(atp, 'BILE_Km_atp')}",
                U.l_per_min,
                name=f"bile secretion ({zone})",
                sboTerm=SBO.MATHEMATICAL_EXPRESSION,
            )
        )
    return out


def _ids(make: Callable[[str, str], str], suffix: str) -> Callable[[str], str]:
    """Get the id function of a zone or site."""
    return lambda name: make(name, suffix)
