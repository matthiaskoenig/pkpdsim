"""The slow (deep) pool of a tissue.

A substance whose terminal phase is the slow release from a deep
compartment of the tissues (a marker of the extracellular volume, a contrast
agent in the days after a dose) has in each organ of the pool the species
`Cslow_<s>` in the tissue compartment `Vti` and its exchange with the merged
tissue,

```
J_slow_<s>: Cti_<s> <-> Cslow_<s> = kslow_<s> * Vti * (Kpslow_<s> * Cti_<s> - Cslow_<s>)
```

with the exchange `kslow_<s>` per litre of tissue (1/min), so the exchange of
an organ follows the tissue volume of the individual, as the share of
`ftissue` does, and the pool to tissue coefficient `Kpslow_<s>`: at
equilibrium `Cslow = Kpslow*Cti`, the pool holds `Kpslow` times the amount of
the tissue. Both are constants of the top model, one value for every organ
of the pool, which a fit changes.
"""

from collections.abc import Iterable
from dataclasses import dataclass

from pkpdsim.pbpk.organs import Organ, as_organs
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, check_positive, check_unit

__all__ = ["SlowPool", "check_slow_pools"]


@dataclass(frozen=True)
class SlowPool:
    """The slow (deep) pool of the tissue of a substance in some organs.

    Attributes:
        substance: the id of the substance.
        organs: the organs with a pool, standard organs in which the
            substance has the merged tissue (`Mode.TISSUE`); given as members,
            values or names, kept as a frozenset of `Organ` members.
        k: the exchange between the tissue and the pool per litre of tissue
            (1/min), `kslow_<s>`: a number, a `Value` (fitted) or a `Formula`.
        kp: the concentration of the pool to that of the tissue at
            equilibrium (dimensionless), `Kpslow_<s>`.
    """

    substance: str
    organs: frozenset[Organ]
    k: Q
    kp: Q

    def __post_init__(self) -> None:
        """Normalise the organs to `Organ` members."""
        object.__setattr__(self, "organs", as_organs(self.organs, "SlowPool.organs"))


def check_slow_pools(pools: Iterable[SlowPool], known: set[str]) -> list[str]:
    """Check the slow pools of a definition, each on its own.

    A known substance, one pool per substance (its constants are shared by
    its organs), at least one organ, `k` and `kp` positive in their units.
    That every organ is a standard organ in which the substance has the
    merged tissue is checked by `PBPKDefinition.validate`.

    Args:
        pools: the slow pools.
        known: the ids of the substances of the definition.

    Returns:
        The problems, naming the substance.
    """
    pools = list(pools)
    problems: list[str] = []
    substances = [p.substance for p in pools]
    twice = sorted({s for s in substances if substances.count(s) > 1})
    if twice:
        problems.append(
            f"several slow pools of {twice}: one pool per substance, its organs "
            "share kslow and Kpslow"
        )
    for p in pools:
        where = f"slow pool of '{p.substance}'"
        if p.substance not in known:
            problems.append(f"{where}: unknown substance")
        if not p.organs:
            problems.append(f"{where}: no organs")
        problems += check_positive(where, "k", p.k)
        problems += check_unit(where, "k", p.k, U.per_min)
        problems += check_positive(where, "kp", p.kp)
        problems += check_unit(where, "kp", p.kp, U.dimensionless)
    return problems
