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
quantity token at both sites. The sites of the circuit (the reservoir `res`,
the tubing) are named by the circuit.

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
    "BUFFER_KA",
    "BUFFER_K_HYDRATION",
    "CELL",
    "F_T",
    "PERFUSATE",
    "RELEASED",
    "alive",
    "bile",
    "cell_site",
    "dead_site",
    "is_zone",
    "necrosis_rate",
    "ph",
    "profile",
    "reaction",
    "species",
    "zone_id",
]

F_T = "f_T"
"""The temperature factor of the rates (dimensionless, 1 at the reference
temperature), defined by the circuit (an Arrhenius initial assignment of the
temperature) and read by the rate laws of the core."""

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


def ph(site: str) -> str:
    """Get the id of the pH of the perfusate at a site."""
    return f"pH_{site}"


def is_zone(zone: str) -> bool:
    """Check whether a text is a zone id (`z1`, `z12`)."""
    return _ZONE.fullmatch(zone) is not None
