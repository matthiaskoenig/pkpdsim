"""The units a flattening copies, merged in the flat model.

A composed model holds the whole unit definitions of its submodels, so that a
submodel is a model which simulates alone, and libsbml's flattening prefixes
the id of every element of a submodel, its unit definitions included: the flat
model of a liver of ten zones carried eleven copies of every unit, 429 unit
definitions against the 39 the model has, and the math of a zone read `1
M01__mM` where the top model reads `1 mM`. Nothing computes differently, and a
reviewer of the ODEs of the model saw ten different `mM` and could not tell
whether the zones share a unit.

`merge_unit_definitions` is the step which drops the copies: it runs on the
flat model of the composition of the zonation engine
(`pkpdsim.zonation.composition.write`), the flattener of the models which have
submodels of their own file, so the flat model of a perfused liver and of a
zonated liver of the PBPK library carries every unit once. The flattener of
the PBPK library (`pkpdsim.pbpk.builder.flat_model`, whose organs are model
definitions of the same file) copies the units as well and does not run it
yet: that is a change of the tracked SBML of every model of the library.
"""

import libsbml

__all__ = ["merge_unit_definitions"]

_UNIT_ATTRIBUTES = (
    "Units",
    "SubstanceUnits",
    "TimeUnits",
    "VolumeUnits",
    "AreaUnits",
    "LengthUnits",
    "ExtentUnits",
)
"""The attributes which name a unit definition: the units of a compartment, a
parameter, a local parameter and an uncertainty of `distrib`, the substance
units of a species and the units of the model itself. They are the unit
SIdRefs of SBML Level 3 together with the `units` of a number of the math,
i.e. what `renameUnitSIdRefs` of libsbml renames; renaming them one by one is
what makes the merge of a few hundred copies a second instead of a minute."""


def merge_unit_definitions(model: libsbml.Model) -> dict[str, str]:
    """Merge the unit definitions of a model which define the same unit.

    Two definitions are the same unit when their name and their units are the
    same, unit for unit (the kind, the exponent, the scale and the multiplier
    of every one), which is what a copy of a flattening is. The definition
    which stays is the first one whose id carries no prefix of the flattening,
    i.e. the id the model wrote itself (`mM` against the copies `M01__mM`,
    `M02__mM`), else the first one; every reference to a dropped id is
    rewritten, in the attributes of the elements, in the uncertainties and in
    the math.

    Args:
        model: the flat model, which is changed in place.

    Returns:
        The id of every dropped definition and the id it was merged into,
        empty for a model whose units are already its own.
    """
    groups: dict[tuple[object, ...], list[str]] = {}
    for k in range(model.getNumUnitDefinitions()):
        definition = model.getUnitDefinition(k)
        groups.setdefault(_key(definition), []).append(definition.getId())
    renames: dict[str, str] = {}
    for sids in groups.values():
        if len(sids) < 2:
            continue
        keep = next((sid for sid in sids if "__" not in sid), sids[0])
        renames.update({sid: keep for sid in sids if sid != keep})
    if not renames:
        return {}
    _rename_attributes(model, renames)
    elements = model.getListOfAllElements()
    for k in range(elements.getSize()):
        element = elements.get(k)
        if isinstance(element, libsbml.UnitDefinition | libsbml.Unit):
            continue
        _rename_attributes(element, renames)
        _rename_math(element, renames)
    for sid in renames:
        model.removeUnitDefinition(sid)
    return renames


def _key(definition: libsbml.UnitDefinition) -> tuple[object, ...]:
    """Get what makes two unit definitions the same unit: the name and the units."""
    units = tuple(
        sorted(
            (
                unit.getKind(),
                unit.getExponentAsDouble(),
                unit.getScale(),
                unit.getMultiplier(),
            )
            for unit in (definition.getUnit(k) for k in range(definition.getNumUnits()))
        )
    )
    return (definition.getName(), units)


def _rename_attributes(element: libsbml.SBase, renames: dict[str, str]) -> None:
    """Rename the unit definitions an element names in its attributes."""
    for attribute in _UNIT_ATTRIBUTES:
        getter = getattr(element, f"get{attribute}", None)
        if getter is None:
            continue
        new = renames.get(getter())
        if new is not None:
            getattr(element, f"set{attribute}")(new)


def _rename_math(element: libsbml.SBase, renames: dict[str, str]) -> None:
    """Rename the unit definitions the numbers of the math of an element name."""
    get_math = getattr(element, "getMath", None)
    set_math = getattr(element, "setMath", None)
    if get_math is None or set_math is None:
        return
    math = get_math()
    if math is None:
        return
    changed = False
    stack = [math]
    while stack:
        node = stack.pop()
        if node.isSetUnits():
            new = renames.get(node.getUnits())
            if new is not None:
                node.setUnits(new)
                changed = True
        stack.extend(node.getChild(k) for k in range(node.getNumChildren()))
    if changed:
        # the math of an element is changed in place (`getMath` of libsbml gives
        # the node of the element, not a copy) and set again, so that the change
        # holds either way; libsbml keeps the node when it is set to its own
        set_math(math)
