"""Organs, blood pools and the parts of an organ: the topology of the body."""

from collections.abc import Iterable, Mapping
from enum import StrEnum
from typing import Any

__all__ = [
    "PORTAL_ORGANS",
    "SYSTEMIC_ORGANS",
    "Organ",
    "Part",
    "Pool",
    "as_member",
    "as_organ_keys",
    "as_organs",
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


def as_member[E: StrEnum](enum: type[E], value: object, name: str) -> E:
    """Get the member of an enum given as member, value or name.

    A definition written in a script or a notebook may give an organ or the
    sex as text: the value (`"ad"`, `"female"`) or the name in any case
    (`"adipose"`, `"FEMALE"`).

    Args:
        enum: the enum.
        value: the member, its value or its name.
        name: the field which holds the value, for the message.

    Returns:
        The member.

    Raises:
        ValueError: if the value is no member, value or name of the enum,
            naming the field and the allowed values.
    """
    if isinstance(value, enum):
        return value
    if isinstance(value, str):
        text = value.lower()
        for member in enum:
            if text in (member.value.lower(), member.name.lower()):
                return member
    allowed = ", ".join(f"'{m.value}' ({m.name.lower()})" for m in enum)
    raise ValueError(
        f"{name}: {value!r} is no {enum.__name__}, the allowed values are "
        f"{allowed}, as {enum.__name__} member, value or name"
    )


def as_organ_keys[V](mapping: Mapping[Any, V], name: str) -> dict[Organ, V]:
    """Get a mapping by organ from a mapping whose keys are organs or their texts.

    Args:
        mapping: the mapping, its keys `Organ` members, values or names.
        name: the field which holds the mapping, for the messages.

    Returns:
        The mapping with `Organ` keys, in the order of the given one.

    Raises:
        ValueError: if the value is no mapping, if a key is no organ or if two
            keys name the same organ.
    """
    if not isinstance(mapping, Mapping):
        raise ValueError(
            f"{name}: {mapping!r} is no mapping by organ, e.g. {{Organ.LIVER: ...}}"
        )
    result: dict[Organ, V] = {}
    for key, value in mapping.items():
        organ = as_member(Organ, key, name)
        if organ in result:
            raise ValueError(f"{name}: the organ '{organ.value}' is given twice")
        result[organ] = value
    return result


def as_organs(organs: Iterable[object], name: str) -> frozenset[Organ]:
    """Get a set of organs from a collection of organs or their texts.

    Args:
        organs: a set, frozenset, list or tuple of `Organ` members, values or
            names.
        name: the field which holds the organs, for the messages.

    Returns:
        The organs.

    Raises:
        ValueError: if the value is no such collection or an element is no
            organ.
    """
    if not isinstance(organs, set | frozenset | list | tuple):
        raise ValueError(
            f"{name}: {organs!r} is no collection of organs, e.g. "
            "frozenset({Organ.LIVER, Organ.KIDNEY})"
        )
    return frozenset(as_member(Organ, organ, name) for organ in organs)
