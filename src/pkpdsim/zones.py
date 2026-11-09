"""The thirds of a chain of zones and the weighted mean over them.

A model which resolves an organ into `n` zones along its length reports the
periportal and the pericentral mean of a quantity: the mean over the first and
the last third of the length, every zone weighted by its overlap with the
third. The weights are the same for every such model, so they are written once
here: the machine perfusion library (`pkpdsim.perfusion.observers`), the
zonation engine (`pkpdsim.zonation.observers`) and the zonated liver of the
PBPK library (`pkpdsim.pbpk.variants.liver_zonated`) take them from this
module, and a zonated model of any of them reports the same two numbers.

The module imports nothing else of pkpdsim, so that the three libraries share
the weights without importing each other (`pkpdsim.species` is a leaf for the
same reason).
"""

from collections.abc import Sequence
from fractions import Fraction

__all__ = ["third_weights", "weighted_mean"]


def third_weights(n: int, periportal: bool) -> list[tuple[int, Fraction]]:
    """Get the weights of the zones in the mean of a third of the sinusoid.

    Zone `k` of `n` covers `[(k - 1)/n, k/n]` of the length of the sinusoid;
    its weight is its overlap with the first (periportal) or the last
    (pericentral) third over the length of the third, so the weights sum to 1.

    Args:
        n: the number of zones.
        periportal: the first third, else the last.

    Returns:
        The zones (from 1) with a weight above 0 and their weights.
    """
    third = Fraction(1, 3)
    out = []
    for k in range(1, n + 1):
        lower, upper = Fraction(k - 1, n), Fraction(k, n)
        if not periportal:
            lower, upper = 1 - upper, 1 - lower
        overlap = min(upper, third) - lower
        if overlap > 0:
            out.append((k, overlap / third))
    return out


def weighted_mean(of_zones: Sequence[str], weights: list[tuple[int, Fraction]]) -> str:
    """Get the weighted mean of an id of the zones (periportal first)."""
    terms = [
        of_zones[k - 1] if w == 1 else f"{float(w)!r} dimensionless * {of_zones[k - 1]}"
        for k, w in weights
    ]
    return " + ".join(terms)
