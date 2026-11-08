"""The id scheme of the machine perfusion models.

Every id of a perfusion model is built here, so that the cellular core, the
circuit and the observers agree on them.

A quantity at a site is `<quantity>_<site>` (`species`): the site of the
vascular space of a zone is the zone (`glc_z3`, the glucose of the
perfusate in the sinusoids of zone 3), the cells of a zone are the site
`<zone>_cell` (`cell_site`, `atp_z3_cell`) and the necrotic cells of a zone
the site `<zone>_dead` (`dead_site`, `alt_z3_dead`, the ALT which a dead
cell has not released yet). A quantity which exists in the perfusate and in
the cells (glucose, oxygen, lactate, potassium, the enzymes) has the same
quantity token at both sites. The sites of the circuit are the reservoir
`res`, the inflow tubing `in` (after the oxygenator, before the organ) and
the outflow tubing `out` (after the organ); the periportal and pericentral
means of the zones are the sites `pp` and `pc` (`pp_cell` for the cells), so
every quantity of a zone has its means under the same pattern (`glc_pp`,
`atp_pc_cell`, `pO2_pc`). The volume of a site is `V_<site>`
(`volume`), of the vascular space of a zone `Vvas_<zone>` and of its cells
`Vcell_<zone>`.

The reactions and rules of the core of a zone carry the zone as suffix,
`<name>_<zone>` (`OXPHOS_z3`, `v_nec_z3`). The constants of the core are
shared by the zones and are ids of the whole model under their name
(`GLC2PYR_Vmax`), as the constants of the perfusate (the bicarbonate buffer)
and the temperature factor `F_T`; a zonal profile of a process is the factor
`p_<process>_<zone>` (`profile`), the pattern of the zonated liver of the
PBPK library.
"""

import re

__all__ = [
    "ACTIVATION_ENERGY",
    "BILE",
    "BILE_RATE",
    "BUFFER_KA",
    "BUFFER_K_HYDRATION",
    "CELL",
    "CO2_GAS",
    "CO2_SOLUBILITY",
    "ENZYMES",
    "F_T",
    "GAS_CONSTANT",
    "INFLOW",
    "O2_TARGET",
    "OUTFLOW",
    "OXYGENATOR",
    "OXYGENATOR_CO2",
    "OXYGENATOR_K",
    "OXYGENATOR_K_CO2",
    "PERFUSATE",
    "PERICENTRAL",
    "PERIPORTAL",
    "Q_PUMP",
    "RELEASED",
    "RESERVOIR",
    "SOLUBILITY",
    "SOLUBILITY_COEFFICIENT",
    "SOLUBILITY_REF",
    "SOLUBILITY_TEMPERATURE",
    "TEMPERATURE",
    "TEMPERATURE_REF",
    "VO2",
    "activity",
    "alive",
    "bile",
    "cell_fraction",
    "cell_site",
    "cells",
    "dead_site",
    "density",
    "flow",
    "fresh",
    "is_zone",
    "mass",
    "necrosis_rate",
    "organ_volume",
    "pco2",
    "ph",
    "po2",
    "profile",
    "reaction",
    "species",
    "specific_activity",
    "vascular",
    "vascular_fraction",
    "volume",
    "zone_id",
]

F_T = "f_T"
"""The temperature factor of the rates (dimensionless, 1 at the reference
temperature), defined by the circuit (an Arrhenius initial assignment of the
temperature) and read by the rate laws of the core."""

TEMPERATURE = "T"
"""The temperature of the perfusion (degC), a constant which an experiment
changes."""

TEMPERATURE_REF = "T_ref"
"""The reference temperature of the Arrhenius factor (degC)."""

ACTIVATION_ENERGY = "Ea"
"""The activation energy of the rates (J/mole)."""

GAS_CONSTANT = "R_gas"
"""The molar gas constant (J/mole/K)."""

SOLUBILITY = "alpha_o2"
"""The solubility of oxygen in the perfusate at the temperature (mM/kPa),
the initial assignment of `T`."""

SOLUBILITY_REF = "alpha_o2_ref"
"""The solubility of oxygen at `T_alpha_o2` (mM/kPa)."""

SOLUBILITY_TEMPERATURE = "T_alpha_o2"
"""The temperature of `alpha_o2_ref` (degC)."""

SOLUBILITY_COEFFICIENT = "k_alpha_o2"
"""The relative change of the solubility of oxygen per K (1/K)."""

Q_PUMP = "Q_pump"
"""The total flow of the pump (l/min)."""

O2_TARGET = "pO2_target"
"""The oxygen partial pressure which the oxygenator sets (kPa)."""

OXYGENATOR_K = "k_oxy"
"""The first order rate constant of the oxygenation of the inflow (1/min)."""

OXYGENATOR = "OXY_in"
"""The reaction of the oxygenator, which sets the oxygen of the inflow."""

CO2_SOLUBILITY = "alpha_co2"
"""The solubility of CO2 in the perfusate (mM/kPa)."""

CO2_GAS = "pCO2_gas"
"""The partial pressure of CO2 of the sweep gas of the oxygenator (kPa)."""

OXYGENATOR_K_CO2 = "k_oxy_co2"
"""The first order rate constant of the exchange of CO2 with the gas (1/min)."""

OXYGENATOR_CO2 = "OXY_co2_in"
"""The reaction of the oxygenator, which sets the CO2 of the inflow."""

BILE = "bile"
"""The volume of the bile of the organ (µl), a state."""

BILE_RATE = "bile_rate"
"""The bile secretion of the organ (µl/min)."""

VO2 = "VO2"
"""The oxygen uptake of the organ from the perfusate (mmole/min)."""

RESERVOIR = "res"
"""The site of the reservoir."""

INFLOW = "in"
"""The site of the inflow tubing, after the oxygenator and before the organ."""

OUTFLOW = "out"
"""The site of the outflow tubing, after the organ and before the reservoir."""

PERIPORTAL = "pp"
"""The site of the periportal mean of the zones (the first third of the
sinusoid)."""

PERICENTRAL = "pc"
"""The site of the pericentral mean of the zones (the last third of the
sinusoid)."""

ENZYMES: tuple[str, ...] = ("alt", "ast", "ldh")
"""The enzymes of the perfusate, whose catalytic activity is observed."""

BUFFER_KA = "Ka_hco3"
"""The apparent dissociation constant of the bicarbonate buffer (mM),
`[H+]*[HCO3-]/[CO2]`, a constant of the perfusate."""

BUFFER_K_HYDRATION = "k_co2"
"""The rate constant of the hydration of CO2 (1/min), a constant of the
perfusate."""

PERFUSATE: tuple[str, ...] = (
    "glc",
    "o2",
    "lac",
    "co2",
    "hco3",
    "h",
    "k",
    "alt",
    "ast",
    "ldh",
)
"""The quantities of the perfusate at every site: glucose, oxygen, lactate,
dissolved CO2, bicarbonate, free protons, potassium and the enzymes ALT, AST
and LDH (in mM)."""

CELL: tuple[str, ...] = (
    "glc",
    "glyc",
    "o2",
    "pyr",
    "lac",
    "atp",
    "adp",
    "nad",
    "nadh",
    "nadp",
    "nadph",
    "ros",
    "gsh",
    "gssg",
    "k",
    "alt",
    "ast",
    "ldh",
)
"""The quantities of the cells of a zone (in mmole per l of cells); `glyc`
is glycogen in glucosyl units, `ros` the reactive oxygen species as H2O2
equivalents."""

RELEASED: tuple[str, ...] = ("alt", "ast", "ldh", "k")
"""The quantities which necrotic cells release into the vascular space."""

_ZONE = re.compile(r"z[1-9][0-9]*")


def zone_id(k: int) -> str:
    """Get the id of a zone of an organ.

    Args:
        k: the number of the zone, from 1 (periportal) to `n` (pericentral).

    Returns:
        The zone id, `z<k>`.

    Raises:
        ValueError: if `k` is smaller than 1.
    """
    if k < 1:
        raise ValueError(f"A zone is numbered from 1, not {k}.")
    return f"z{k}"


def species(quantity: str, site: str) -> str:
    """Get the id of a quantity at a site, `<quantity>_<site>`.

    Args:
        quantity: the quantity token (`glc`, `o2`, `atp`).
        site: the site (`res`, the zone `z3` for its vascular space, the
            cells `z3_cell`).

    Returns:
        The id of the species.
    """
    return f"{quantity}_{site}"


def cell_site(zone: str) -> str:
    """Get the site of the living cells of a zone, `<zone>_cell`."""
    return f"{zone}_cell"


def dead_site(zone: str) -> str:
    """Get the site of the necrotic cells of a zone, `<zone>_dead`."""
    return f"{zone}_dead"


def profile(process: str, zone: str) -> str:
    """Get the id of the zonal factor of a process, `p_<process>_<zone>`."""
    return f"p_{process}_{zone}"


def reaction(name: str, zone: str) -> str:
    """Get the id of a reaction of the core of a zone, `<name>_<zone>`."""
    return f"{name}_{zone}"


def alive(zone: str) -> str:
    """Get the id of the living fraction of the cells of a zone."""
    return f"f_alive_{zone}"


def necrosis_rate(zone: str) -> str:
    """Get the id of the rate of necrosis of the cells of a zone (1/min)."""
    return f"v_nec_{zone}"


def bile(zone: str) -> str:
    """Get the id of the volume of bile which a zone has secreted (l)."""
    return f"Vbile_{zone}"


def volume(site: str) -> str:
    """Get the volume of the perfusate at a site, `V_<site>` or `Vvas_<zone>`."""
    return vascular(site) if is_zone(site) else f"V_{site}"


def vascular(zone: str) -> str:
    """Get the compartment of the vascular space of a zone, `Vvas_<zone>`."""
    return f"Vvas_{zone}"


def cells(zone: str) -> str:
    """Get the compartment of the cells of a zone, `Vcell_<zone>`."""
    return f"Vcell_{zone}"


def mass(organ: str) -> str:
    """Get the mass of an organ (kg), `M_<organ>`."""
    return f"M_{organ}"


def density(organ: str) -> str:
    """Get the density of an organ (kg/l), `rho_<organ>`."""
    return f"rho_{organ}"


def organ_volume(organ: str) -> str:
    """Get the volume of an organ (l), `V_<organ>`."""
    return f"V_{organ}"


def vascular_fraction(organ: str) -> str:
    """Get the vascular fraction of the volume of an organ, `f_vas_<organ>`."""
    return f"f_vas_{organ}"


def cell_fraction(organ: str) -> str:
    """Get the cell fraction of the volume of an organ, `f_cell_<organ>`."""
    return f"f_cell_{organ}"


def flow(quantity: str, source: str, target: str) -> str:
    """Get the reaction of the flow of a quantity between two sites."""
    return f"FLOW_{quantity}_{source}_{target}"


def fresh(quantity: str) -> str:
    """Get the concentration of a quantity in the fresh perfusate, `<q>_fresh`."""
    return f"{quantity}_fresh"


def po2(site: str) -> str:
    """Get the oxygen partial pressure of the perfusate at a site (kPa)."""
    return f"pO2_{site}"


def pco2(site: str) -> str:
    """Get the CO2 partial pressure of the perfusate at a site (kPa)."""
    return f"pCO2_{site}"


def activity(enzyme: str, site: str) -> str:
    """Get the catalytic activity of an enzyme at a site (µkat/l)."""
    return f"{enzyme}_act_{site}"


def specific_activity(enzyme: str) -> str:
    """Get the specific activity of an enzyme (µkat/mmole), `a_<enzyme>`."""
    return f"a_{enzyme}"


def ph(site: str) -> str:
    """Get the id of the pH of the perfusate at a site."""
    return f"pH_{site}"


def is_zone(zone: str) -> bool:
    """Check whether a text is a zone id (`z1`, `z12`)."""
    return _ZONE.fullmatch(zone) is not None
