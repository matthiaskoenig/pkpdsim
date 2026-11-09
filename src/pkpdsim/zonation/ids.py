"""The id scheme of the zonated liver.

Every id of the engine is built here, so that the zones, the layers, the
transport and the observers agree on them. The zone token is `z<k>` from
`z1` (periportal) to `z<n>` (pericentral), the one of the PBPK library
(`pkpdsim.pbpk.ids.zone`) and of the perfusion (`pkpdsim.perfusion.ids.zone_id`).

A layer of a zone is the site `<layer>_<zone>` (`disse_z3`, `site`), its
volume the compartment `V_<layer>_<zone>` (`volume`) and the amount or
concentration of a quantity there `<quantity>_<layer>_<zone>` (`species`,
`glc_cell_z2`). The zonal factor of a process is `p_<process>_<zone>`
(`profile`), the pattern of `pkpdsim.perfusion.ids.profile`, and the submodel
of zone `k` is `M<kk>` (`submodel`, `M03`).

The liver parameters the layer volumes are computed from are the mass `M_li`
(kg), the density `rho_li` (kg/l), the parenchymal fraction `f_par` and the
volume fraction `f_<layer>` of every layer. The well-stirred geometry has the
vascular site `vas` (lumen, endothelium and Disse) instead of those layers.

The zone token and the reaction of a flow between two sites are the engine's
and not the perfusion library's (`is_zone`, `flow`): the engine names the
flows of its chain, and a perfused organ names the flows of its circuit, with
the same function, so the two agree. `pkpdsim.perfusion.ids` re-exports them.
"""

import re

from pkpdsim.zonation.geometry import Layer

__all__ = [
    "BILE_FLOW",
    "BILE_FLOW_REFERENCE",
    "DECAY_LENGTH",
    "DENSITY",
    "DUCT",
    "FENESTRA_COUNT",
    "FENESTRA_RADIUS",
    "FLOW",
    "HAEMATOCRIT",
    "INFLOW",
    "LENGTH",
    "LUMEN_FRACTION",
    "MASS",
    "MASS_REFERENCE",
    "N_SINUSOIDS",
    "OUTFLOW",
    "PARENCHYMA",
    "POROSITY_MEAN",
    "PRESSURE_CENTRAL",
    "PRESSURE_ONCOTIC",
    "PRESSURE_PORTAL",
    "RADIUS",
    "RESISTANCE_LUMEN",
    "RESISTANCE_PORES",
    "TEMPERATURE",
    "T_DISSE",
    "T_ENDO",
    "VASCULAR",
    "VISCOSITY",
    "VISCOSITY_FACTOR",
    "VISCOSITY_LIVER",
    "area",
    "axial",
    "exchange",
    "flat",
    "flow",
    "fraction",
    "is_zone",
    "layer_flow",
    "meta_id",
    "permeability",
    "porosity",
    "position",
    "profile",
    "sinusoid_flow",
    "site",
    "species",
    "submodel",
    "volume",
    "zone",
    "zone_flow",
]

MASS = "M_li"
"""The id of the liver mass (kg)."""

DENSITY = "rho_li"
"""The id of the liver density (kg/l)."""

PARENCHYMA = "f_par"
"""The id of the parenchymal fraction of the liver mass."""

VASCULAR = "vas"
"""The token of the vascular site of the well-stirred geometry."""

FLOW = "Q_li"
"""The id of the blood flow through the whole liver (l/min)."""

HAEMATOCRIT = "Hct_li"
"""The id of the haematocrit of the sinusoidal lumen (fraction)."""

LENGTH = "L_sin"
"""The id of the length of a sinusoid (m), which the zones divide."""

RADIUS = "r_sin"
"""The id of the radius of the lumen of a sinusoid (m)."""

T_ENDO = "t_endo"
"""The id of the thickness of the endothelium (m)."""

T_DISSE = "t_disse"
"""The id of the thickness of the space of Disse (m)."""

TEMPERATURE = "temperature"
"""The id of the temperature of the liver (degC), one for all zones."""

INFLOW = "in"
"""The token of the boundary upstream of the first zone (`plasma_in`)."""

OUTFLOW = "out"
"""The token of the boundary downstream of the last zone (`plasma_out`)."""

DUCT = "duct"
"""The token of the bile duct, downstream of the canaliculus of the first zone
(`bile_duct`): the collector of the bile, a boundary and no zone."""

PRESSURE_PORTAL = "P_pp"
"""The id of the pressure at the periportal end of the sinusoid (Pa)."""

PRESSURE_CENTRAL = "P_pc"
"""The id of the pressure at the pericentral end of the sinusoid (Pa)."""

PRESSURE_ONCOTIC = "P_onc"
"""The id of the oncotic pressure of the leaky tube (Pa), `Poc - Pot` of spt."""

VISCOSITY = "nu_plasma"
"""The id of the viscosity of the plasma (Pa s)."""

VISCOSITY_FACTOR = "f_nu_li"
"""The id of the factor of the hepatic resistance over the plasma viscosity."""

VISCOSITY_LIVER = "nu_li"
"""The id of the viscosity of the blood in the sinusoid (Pa s), the factor times the plasma viscosity."""

FENESTRA_RADIUS = "r_fen"
"""The id of the radius of a fenestra (m)."""

FENESTRA_COUNT = "N_fen"
"""The id of the number of fenestrae per endothelial area (1/m**2), derived from the porosity."""

POROSITY_MEAN = "eps_fen_li"
"""The id of the fenestral porosity the leaky tube is evaluated at (the mean over the sinusoid)."""

LUMEN_FRACTION = "f_lumen"
"""The id of the volume fraction of the lumen (plasma and erythrocytes) in the liver."""

N_SINUSOIDS = "N_sin"
"""The id of the number of sinusoids of the liver (dimensionless)."""

RESISTANCE_LUMEN = "W_sin"
"""The id of the specific hydraulic resistance of the lumen (Pa s/m**4)."""

RESISTANCE_PORES = "w_sin"
"""The id of the specific hydraulic resistance of all pores (Pa s/m**2)."""

DECAY_LENGTH = "lambda_sin"
"""The id of the length (m) over which the pressure of the leaky tube relaxes."""

BILE_FLOW = "Q_bile"
"""The id of the bile flow of the liver (l/min), through every canaliculus."""

BILE_FLOW_REFERENCE = "Q_bile_ref"
"""The id of the bile flow (l/min) of the reference liver of the species."""

MASS_REFERENCE = "M_li_ref"
"""The id of the mass (kg) of the reference liver the bile flow is given for."""

_ZONE = re.compile(r"z[1-9][0-9]*")
"""The pattern of the token of a zone, `z1` to `z<n>`."""


def zone(k: int) -> str:
    """Get the token of zone `k`, `z<k>`.

    Args:
        k: the number of the zone, from 1 (periportal).

    Returns:
        The zone token.

    Raises:
        ValueError: for a zone below 1.
    """
    if k < 1:
        raise ValueError(f"The zones are numbered from 1, not {k}.")
    return f"z{k}"


def is_zone(zone: str) -> bool:
    """Check whether a text is a zone id (`z1`, `z12`)."""
    return _ZONE.fullmatch(zone) is not None


def site(layer: Layer | str, zone: str) -> str:
    """Get the site of a layer in a zone, `<layer>_<zone>`."""
    return f"{layer}_{zone}"


def volume(layer: Layer | str, zone: str) -> str:
    """Get the compartment of a layer in a zone, `V_<layer>_<zone>`."""
    return f"V_{site(layer, zone)}"


def species(quantity: str, layer: Layer | str, zone: str) -> str:
    """Get the id of a quantity in a layer of a zone, `<quantity>_<layer>_<zone>`."""
    return f"{quantity}_{site(layer, zone)}"


def meta_id(sid: str) -> str:
    """Get the metaid of an element, `meta_<sid>`.

    A compartment which a cell model replaces needs one, or the comp model
    cannot be flattened.
    """
    return f"meta_{sid}"


def position(zone: str) -> str:
    """Get the id of the lobular position of a zone, `position_<zone>`."""
    return f"position_{zone}"


def profile(process: str, zone: str) -> str:
    """Get the id of the zonal factor of a process, `p_<process>_<zone>`."""
    return f"p_{process}_{zone}"


def submodel(k: int) -> str:
    """Get the id of the submodel of zone `k`, `M<kk>`.

    The two-digit padding is a minimum width for readability. Above 99 zones
    the ids widen (`M100`) and no longer sort lexically; they stay unique.

    Raises:
        ValueError: for a zone below 1.
    """
    if k < 1:
        raise ValueError(f"The zones are numbered from 1, not {k}.")
    return f"M{k:02}"


def flat(k: int, sid: str) -> str:
    """Get the id of an element of the cell model of zone `k` in the flat model.

    The flattening of a comp model prefixes the elements of a submodel with its
    id, `M<kk>__<id>`.
    """
    return f"{submodel(k)}__{sid}"


def fraction(layer: Layer | str) -> str:
    """Get the id of the volume fraction of a layer of the liver, `f_<layer>`."""
    return f"f_{layer}"


def layer_flow(layer: Layer | str) -> str:
    """Get the flow of a layer through the chain, `Q_<layer>_li` (l/min).

    The plasma and the erythrocytes carry the shares `1 - Hct` and `Hct` of the
    blood flow `Q_li`, so they have one velocity.
    """
    return f"Q_{layer}_li"


def zone_flow(layer: Layer | str | None, zone: str) -> str:
    """Get the flow of a layer in a zone, `Q_<layer>_li_<zone>` (l/min).

    A pressure-driven sinusoid has a flow which depends on the position, so
    every zone has its own flow. `None` (and the vascular site of the
    well-stirred geometry) is the blood flow of the zone, `Q_li_<zone>`.
    """
    if layer is None or layer == VASCULAR:
        return f"Q_li_{zone}"
    return f"Q_{layer}_li_{zone}"


def sinusoid_flow(zone: str) -> str:
    """Get the flow of one sinusoid in a zone, `Q_sin_<zone>` (m**3/s)."""
    return f"Q_sin_{zone}"


def flow(quantity: str, source: str, target: str) -> str:
    """Get the reaction of the flow of a quantity between two sites.

    The vascular space of a zone is the site `vas_<zone>` of the engine.
    """
    return f"FLOW_{quantity}_{_flow_site(source)}_{_flow_site(target)}"


def _flow_site(site: str) -> str:
    """Get the site of a flow, the vascular site of the engine for a zone."""
    return f"{VASCULAR}_{site}" if is_zone(site) else site


def diffusion(quantity: str) -> str:
    """Get the id of the diffusion coefficient of a quantity in water, `D_<q>` (m**2/s)."""
    return f"D_{quantity}"


def axial(quantity: str, source: str, target: str) -> str:
    """Get the axial diffusion of a quantity between two sites, `D_<q>_<from>_<to>`."""
    return f"D_{quantity}_{source}_{target}"


def porosity(zone: str) -> str:
    """Get the id of the fenestral porosity of a zone, `eps_fen_<zone>`."""
    return f"eps_fen_{zone}"


def area(interface: str, zone: str) -> str:
    """Get the id of the lateral area of an interface in a zone, `A_<interface>_<zone>`.

    The interfaces are `pd` (plasma to the space of Disse, the lumen), `pe` (the
    abluminal face of the endothelial cell) and `dh` (the space of Disse to the
    hepatocyte).
    """
    return f"A_{interface}_{zone}"


def permeability(path: str, quantity: str) -> str:
    """Get the id of a permeability of a substance, `P_<path>_<quantity>`.

    The paths are `fen` (the fenestrae, of an open pore), `endo_in` and
    `endo_out` (the luminal and the abluminal membrane of the endothelial cell)
    and `hep` (the basolateral membrane of the hepatocyte).
    """
    return f"P_{path}_{quantity}"


def exchange(kind: str, quantity: str, zone: str, suffix: str = "") -> str:
    """Get the id of a radial exchange reaction, `<kind>_<q>_<zone>[_<suffix>]`.

    The kinds are `FEN` (plasma to Disse through the fenestrae), `END` (plasma to
    endothelial cell to Disse; the suffix `in` or `out`) and `PS` (Disse to
    hepatocyte).
    """
    base = f"{kind}_{quantity}_{zone}"
    return f"{base}_{suffix}" if suffix else base
