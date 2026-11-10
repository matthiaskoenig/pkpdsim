"""Shared helpers of the zonation tests, plain functions and not fixtures."""

from pkpdsim.pbpk.physiology import reference_flow
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion.units import PU
from pkpdsim.species import Species


def mass_g(x: float) -> Value:
    """A liver mass in g."""
    return Value(x, PU.g)


def liver_flow(species: Species = Species.RAT) -> Value:
    """The reference blood flow of the whole liver of a species (l/min).

    The engine reads no physiology: `ZonationSettings.reference_flow` is the
    number its caller gives, and the caller of these tests is the PBPK
    library (`pkpdsim.pbpk.physiology.reference_flow`), as the zonated liver
    of the library is. A test which exercises a flow mode at its default (an
    `ImposedFlow` without a value, a `PressureFlow` without a resistance
    factor) states it.
    """
    return reference_flow(species)
