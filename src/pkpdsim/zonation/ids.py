"""The id scheme of the zonated liver.

Every id of the engine is built here, so that the zones, the layers, the
transport and the observers agree on them. The zone token is `z<k>` from
`z1` (periportal) to `z<n>` (pericentral), the one of the PBPK library
(`pkpdsim.pbpk.ids.zone`) and of the perfusion (`pkpdsim.perfusion.ids.zone_id`).

A layer of a zone is the site `<layer>_<zone>` (`disse_z3`, `site`), its
volume the compartment `V_<layer>_<zone>` (`volume`) and the amount or
concentration of a quantity there `<quantity>_<layer>_<zone>` (`species`,
`glc_hepatocyte_z2`). The zonal factor of a process is `p_<process>_<zone>`
(`profile`), the pattern of `pkpdsim.perfusion.ids.profile`, and the submodel
of zone `k` is `M<kk>` (`submodel`, `M03`).

The liver parameters the layer volumes are computed from are the mass `M_li`
(kg), the density `rho_li` (kg/l), the parenchymal fraction `f_par` and the
volume fraction `f_<layer>` of every layer. The well-stirred geometry has the
vascular site `vas` (lumen, endothelium and Disse) instead of those layers.
"""

from pkpdsim.zonation.geometry import Layer

__all__ = [
    "DENSITY",
    "FLOW",
    "HAEMATOCRIT",
    "INFLOW",
    "LENGTH",
    "MASS",
    "OUTFLOW",
    "PARENCHYMA",
    "RADIUS",
    "T_DISSE",
    "T_ENDO",
    "VASCULAR",
    "area",
    "axial",
    "exchange",
    "fraction",
    "layer_flow",
    "permeability",
    "porosity",
    "profile",
    "site",
    "species",
    "submodel",
    "volume",
    "zone",
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

INFLOW = "in"
"""The token of the boundary upstream of the first zone (`plasma_in`)."""

OUTFLOW = "out"
"""The token of the boundary downstream of the last zone (`plasma_out`)."""


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


def site(layer: Layer | str, zone: str) -> str:
    """Get the site of a layer in a zone, `<layer>_<zone>`."""
    return f"{layer}_{zone}"


def volume(layer: Layer | str, zone: str) -> str:
    """Get the compartment of a layer in a zone, `V_<layer>_<zone>`."""
    return f"V_{site(layer, zone)}"


def species(quantity: str, layer: Layer | str, zone: str) -> str:
    """Get the id of a quantity in a layer of a zone, `<quantity>_<layer>_<zone>`."""
    return f"{quantity}_{site(layer, zone)}"


def profile(process: str, zone: str) -> str:
    """Get the id of the zonal factor of a process, `p_<process>_<zone>`."""
    return f"p_{process}_{zone}"


def submodel(k: int) -> str:
    """Get the id of the submodel of zone `k`, `M<kk>`.

    Raises:
        ValueError: for a zone below 1.
    """
    if k < 1:
        raise ValueError(f"The zones are numbered from 1, not {k}.")
    return f"M{k:02}"


def fraction(layer: Layer | str) -> str:
    """Get the id of the volume fraction of a layer of the liver, `f_<layer>`."""
    return f"f_{layer}"


def layer_flow(layer: Layer | str) -> str:
    """Get the flow of a layer through the chain, `Q_<layer>_li` (l/min).

    The plasma and the erythrocytes carry the shares `1 - Hct` and `Hct` of the
    blood flow `Q_li`, so they have one velocity.
    """
    return f"Q_{layer}_li"


def axial(quantity: str, source: str, target: str) -> str:
    """Get the axial diffusion of a quantity between two sites, `D_<q>_<from>_<to>`."""
    return f"D_{quantity}_{source}_{target}"


def porosity(zone: str) -> str:
    """Get the id of the fenestral porosity of a zone, `eps_fen_<zone>`."""
    return f"eps_fen_{zone}"


def area(interface: str, zone: str) -> str:
    """Get the id of the lateral area of an interface in a zone, `A_<interface>_<zone>`.

    The interfaces are `pd` (plasma to the space of Disse) and `dh` (the space
    of Disse to the hepatocyte).
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
