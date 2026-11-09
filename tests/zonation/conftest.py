"""Shared helpers of the zonation tests, plain functions and not fixtures."""

from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion.units import PU


def mass_g(x: float) -> Value:
    """A liver mass in g."""
    return Value(x, PU.g)
