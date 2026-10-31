"""The slow (deep) pool of a tissue.

A substance whose terminal phase is the slow release from a deep
compartment of the tissues (a marker of the extracellular volume, a contrast
agent in the days after a dose) has in each organ of the pool the species
`Cslow_<s>` in the tissue compartment `Vti` and its exchange with the merged
tissue (`source=Part.TISSUE`, the default),

```
J_slow_<s>: Cti_<s> <-> Cslow_<s> = kslow_<s> * Vti * (Kpslow_<s> * Cti_<s> - Cslow_<s>)
```

or with the plasma of the organ (`source=Part.PLASMA`), in parallel with the
merged tissue: a deep extracellular space which the substance reaches
through the capillary wall, not through the fast interstitium,

```
J_slow_<s>: Cpl_<s> <-> Cslow_<s> = kslow_<s> * Vti * (Kpslow_<s> * Cpl_<s> - Cslow_<s>)
```

with the exchange `kslow_<s>` per litre of tissue (1/min) in both forms, so
the exchange of an organ follows the tissue volume of the individual, as the
share of `ftissue` does, and the coefficient `Kpslow_<s>` of the pool to its
source: at equilibrium `Cslow = Kpslow*Cti` (`Cslow = Kpslow*Cpl`), the pool
holds `Kpslow` times the concentration of its source in the tissue volume.
Both are constants of the top model, one value for every organ of the pool,
which a fit changes.

The plasma form reads the total plasma concentration, as the Kp of a tissue
does: `Kpslow_<s>` is relative to the total plasma and the rate has no fu,
so a change of fu (the CKD profile of the binding) changes neither the
equilibrium nor the rate of the pool of a bound substance.
"""

from collections.abc import Iterable
from dataclasses import dataclass

from pkpdsim.pbpk.organs import Organ, Part, as_member, as_organs
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Q, check_positive, check_unit

__all__ = ["SOURCES", "OrganCollection", "SlowPool", "check_slow_pools"]

SOURCES: tuple[Part, ...] = (Part.TISSUE, Part.PLASMA)
"""What a slow pool exchanges with: the merged tissue or the plasma of the organ."""

type OrganCollection = (
    frozenset[Organ | str]
    | set[Organ | str]
    | list[Organ | str]
    | tuple[Organ | str, ...]
)
"""The organs as `SlowPool.organs` takes them: a frozenset, set, list or tuple of
`Organ` members, values or names (`as_organs`)."""


@dataclass(frozen=True)
class SlowPool:
    """The slow (deep) pool of the tissue of a substance in some organs.

    Attributes:
        substance: the id of the substance.
        organs: the organs with a pool, standard organs in which the
            substance has the merged tissue (`Mode.TISSUE`); a frozenset,
            set, list or tuple of members, values or names
            (`OrganCollection`), normalised to a frozenset of `Organ` members
            in `__post_init__`, which is what the field holds afterwards.
        k: the exchange between the tissue and the pool per litre of tissue
            (1/min), `kslow_<s>`: a number, a `Value` (fitted) or a `Formula`.
        kp: the concentration of the pool to that of its source at
            equilibrium (dimensionless), `Kpslow_<s>`.
        source: what the pool exchanges with (`SOURCES`): the merged tissue
            `Cti_<s>` (`Part.TISSUE`, the default) or the plasma of the organ
            `Cpl_<s>` (`Part.PLASMA`, a deep space reached through the
            capillary wall in parallel with the merged tissue); a member,
            value or name, normalised to the `Part` member in
            `__post_init__`.
    """

    substance: str
    organs: OrganCollection
    k: Q
    kp: Q
    source: Part | str = Part.TISSUE

    def __post_init__(self) -> None:
        """Normalise the organs to a frozenset of `Organ` members, the source to a `Part`."""
        object.__setattr__(self, "organs", as_organs(self.organs, "SlowPool.organs"))
        object.__setattr__(
            self, "source", as_member(Part, self.source, "SlowPool.source")
        )


def check_slow_pools(pools: Iterable[SlowPool], known: set[str]) -> list[str]:
    """Check the slow pools of a definition, each on its own.

    A known substance, one pool per substance (its constants are shared by
    its organs), at least one organ, the tissue or the plasma as its source
    (`SOURCES`), `k` and `kp` positive in their units.
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
        if p.source not in SOURCES:
            problems.append(
                f"{where}: the source '{p.source}' is neither the tissue nor the "
                "plasma of the organ, the pool exchanges with one of them "
                "(Part.TISSUE, Part.PLASMA)"
            )
        problems += check_positive(where, "k", p.k)
        problems += check_unit(where, "k", p.k, U.per_min)
        problems += check_positive(where, "kp", p.kp)
        problems += check_unit(where, "kp", p.kp, U.dimensionless)
    return problems
