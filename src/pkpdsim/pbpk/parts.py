"""The collapsing of the parts of a substance in an organ."""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

from pkpdsim.pbpk.chemistry import CustomReaction, Metabolism, Transport
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.substance import Mode, Substance

__all__ = ["Parts", "parts_for"]


@dataclass(frozen=True)
class Parts:
    """The parts which are built for a substance in an organ.

    Attributes:
        blood_cells: a blood cell species is built (kinetic erythrocytes).
        flow_carrier: what the blood flows carry, the plasma concentration or
            the blood concentration `BP*C_p` (erythrocytes in equilibrium).
        tissue: the tissue species: `none` (vascular only), `organ` (one
            species for the whole organ, flow-limited), `merged`
            (interstitial and intracellular merged, finite exchange) or
            `full` (interstitial and intracellular species).
    """

    blood_cells: bool
    flow_carrier: Literal["plasma", "blood_bp"]
    tissue: Literal["none", "organ", "merged", "full"]


_TISSUE: dict[Mode, Literal["organ", "merged", "full"]] = {
    Mode.PERFUSION: "organ",
    Mode.TISSUE: "merged",
    Mode.PERMEABILITY: "full",
}


def parts_for(
    substance: Substance,
    organ: Organ,
    chemistry: Sequence[Metabolism | Transport | CustomReaction],
) -> Parts:
    """Get the parts of a substance in an organ.

    The collapsing table of the spec: the flags of the substance and the
    chemistry of the organ decide, so only what is needed is built. A
    transport of the substance in the organ (basolateral, or an apical
    process) forces the full form, also when the substance has no
    distribution or does not enter the organ. A product of a metabolism
    gets no tissue by this, it is placed in the plasma part.

    Args:
        substance: the substance.
        organ: the organ.
        chemistry: the processes of the definition (of all organs).

    Returns:
        The parts.
    """
    ery = substance.erythrocytes
    if ery is None:
        blood_cells, carrier = False, "plasma"
    elif ery.mode == "equilibrium":
        blood_cells, carrier = False, "blood_bp"
    else:
        blood_cells, carrier = True, "plasma"

    dist = substance.distribution
    transported = any(
        isinstance(p, Transport) and p.organ is organ and p.substance == substance.sid
        for p in chemistry
    )
    tissue: Literal["none", "organ", "merged", "full"]
    if transported:
        tissue = "full"
    elif dist is None or (dist.organs is not None and organ not in dist.organs):
        tissue = "none"
    else:
        tissue = _TISSUE[dist.mode]
    return Parts(blood_cells=blood_cells, flow_carrier=carrier, tissue=tissue)
