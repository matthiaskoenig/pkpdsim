"""The registry of the organ variants.

A variant is an `OrganComponent` registered under a name for an organ; a
definition selects it by name (`Variants(organs={Organ.LIVER: "zonated"})`)
or gives a component. Every organ has the standard organ (`standard`); the
intestine has the segmented intestine (`segmented`, the default) and the
single lumen (`single_lumen`), which replace small and large intestine with
one submodel `GU`.
"""

from pkpdsim.pbpk.organ import OrganComponent, StandardOrgan
from pkpdsim.pbpk.organs import Organ
from pkpdsim.pbpk.variants.intestine import SegmentedIntestine, SingleLumen

__all__ = ["DEFAULTS", "lookup", "names", "register"]

_REGISTRY: dict[tuple[Organ, str], OrganComponent] = {}
"""The variants by organ and name."""

DEFAULTS: dict[Organ, str] = {Organ.SMALL_INTESTINE: "segmented"}
"""The variant of an organ without a choice of the definition, `standard` otherwise."""


def register(organ: Organ, name: str, component: OrganComponent) -> None:
    """Register a variant of an organ under a name.

    Args:
        organ: the organ.
        name: the name of the variant.
        component: the component.

    Raises:
        ValueError: if the organ has a variant of the name already.
    """
    if (organ, name) in _REGISTRY:
        raise ValueError(f"Organ '{organ.value}' has a variant '{name}' already.")
    _REGISTRY[(organ, name)] = component


def names(organ: Organ) -> list[str]:
    """Get the names of the variants of an organ.

    Args:
        organ: the organ.

    Returns:
        The sorted names.
    """
    return sorted(name for o, name in _REGISTRY if o is organ)


def lookup(organ: Organ, name: str) -> OrganComponent:
    """Get the variant of an organ by name.

    Args:
        organ: the organ.
        name: the name of the variant.

    Returns:
        The component.

    Raises:
        ValueError: if the organ has no variant of the name, naming the organ,
            the name and the variants of the organ.
    """
    try:
        return _REGISTRY[(organ, name)]
    except KeyError:
        raise ValueError(
            f"Unknown variant '{name}' of organ '{organ.value}', the variants "
            f"are {names(organ)}."
        ) from None


for _organ in Organ:
    register(_organ, "standard", StandardOrgan())
register(Organ.SMALL_INTESTINE, "segmented", SegmentedIntestine())
register(Organ.SMALL_INTESTINE, "single_lumen", SingleLumen())
