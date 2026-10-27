"""Tests of the organ topology."""

from pkpdsim.pbpk.organs import (
    PORTAL_ORGANS,
    SYSTEMIC_ORGANS,
    Organ,
    Part,
    Pool,
    inflow,
    outflow,
)


def test_values() -> None:
    """The enum values are the id prefixes."""
    assert Organ.LIVER == "li"
    assert Pool.PORTAL == "po"
    assert Part.TISSUE == "ti"
    assert len(Organ) == 15


def test_sets() -> None:
    """The portal and systemic organs partition the rest."""
    assert {
        Organ.STOMACH,
        Organ.SMALL_INTESTINE,
        Organ.LARGE_INTESTINE,
        Organ.SPLEEN,
        Organ.PANCREAS,
    } == PORTAL_ORGANS
    assert Organ.LUNG not in SYSTEMIC_ORGANS
    assert Organ.LIVER not in SYSTEMIC_ORGANS
    assert not SYSTEMIC_ORGANS & PORTAL_ORGANS
    assert SYSTEMIC_ORGANS | PORTAL_ORGANS | {Organ.LUNG, Organ.LIVER} == set(Organ)


def test_inflow() -> None:
    """Liver has two inflows, lung the venous pool, others the arterial."""
    assert inflow(Organ.LIVER) == (Pool.ARTERIAL, Pool.PORTAL)
    assert inflow(Organ.LUNG) == (Pool.VENOUS,)
    for organ in SYSTEMIC_ORGANS | PORTAL_ORGANS:
        assert inflow(organ) == (Pool.ARTERIAL,)


def test_outflow() -> None:
    """Every organ drains into exactly one pool."""
    for organ in PORTAL_ORGANS:
        assert outflow(organ) is Pool.PORTAL
    assert outflow(Organ.LUNG) is Pool.ARTERIAL
    for organ in SYSTEMIC_ORGANS | {Organ.LIVER}:
        assert outflow(organ) is Pool.VENOUS
