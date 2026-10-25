"""Organs, blood pools and the parts of an organ: the topology of the body."""

from enum import StrEnum

__all__ = [
    "PORTAL_ORGANS",
    "SYSTEMIC_ORGANS",
    "Organ",
    "Part",
    "Pool",
    "inflow",
    "outflow",
]


class Organ(StrEnum):
    """The organs of the whole-body model, by their id prefix."""

    LUNG = "lu"
    HEART = "he"
    BRAIN = "br"
    LIVER = "li"
    KIDNEY = "ki"
    STOMACH = "st"
    SMALL_INTESTINE = "si"
    LARGE_INTESTINE = "la"
    SPLEEN = "sp"
    PANCREAS = "pa"
    MUSCLE = "mu"
    SKIN = "sk"
    ADIPOSE = "ad"
    BONE = "bo"
    GONADS = "go"


class Pool(StrEnum):
    """The blood pools of the whole-body model."""

    ARTERIAL = "ar"
    VENOUS = "ve"
    PORTAL = "po"


class Part(StrEnum):
    """The parts of an organ (compartments and their species)."""

    PLASMA = "pl"
    BLOOD_CELLS = "bc"
    INTERSTITIAL = "is"
    INTRACELLULAR = "ic"
    TISSUE = "ti"
    ORGAN = "or"


PORTAL_ORGANS: frozenset[Organ] = frozenset(
    {
        Organ.STOMACH,
        Organ.SMALL_INTESTINE,
        Organ.LARGE_INTESTINE,
        Organ.SPLEEN,
        Organ.PANCREAS,
    }
)
"""The organs which drain into the portal vein."""

SYSTEMIC_ORGANS: frozenset[Organ] = (
    frozenset(Organ) - PORTAL_ORGANS - {Organ.LUNG, Organ.LIVER}
)
"""The organs in parallel between arterial and venous pool."""


def inflow(organ: Organ) -> tuple[Pool, ...]:
    """The pools which supply an organ.

    Args:
        organ: the organ.

    Returns:
        The pools: arterial and portal for the liver, venous for the lung,
        arterial for all others.
    """
    if organ is Organ.LIVER:
        return (Pool.ARTERIAL, Pool.PORTAL)
    if organ is Organ.LUNG:
        return (Pool.VENOUS,)
    return (Pool.ARTERIAL,)


def outflow(organ: Organ) -> Pool:
    """The pool into which an organ drains.

    Args:
        organ: the organ.

    Returns:
        Portal for the portal organs, arterial for the lung, venous for all others.
    """
    if organ in PORTAL_ORGANS:
        return Pool.PORTAL
    if organ is Organ.LUNG:
        return Pool.ARTERIAL
    return Pool.VENOUS
