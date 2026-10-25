"""The definition of a PBPK model: the substances, the chemistry and the choices."""

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.chemistry import (
    MM,
    CustomReaction,
    FirstOrder,
    Hill,
    Metabolism,
    Transport,
)
from pkpdsim.pbpk.organ import OrganComponent
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.pd import PDComponent
from pkpdsim.pbpk.physiology import Sex
from pkpdsim.pbpk.proteins import PROTEINS
from pkpdsim.pbpk.routes import Route
from pkpdsim.pbpk.substance import Erythrocytes, Kp, KpMethod, Mode, Substance
from pkpdsim.pbpk.values import Q, Value

__all__ = ["PBPKDefinition", "Variants"]

_SID = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_RESERVED = (
    *(p.value for p in Part),
    *(p.value for p in Pool),
    "free",
    "blood",
)
"""The tokens of the id scheme which a substance id must not start with."""


def _physiology_ids() -> frozenset[str]:
    """Get all ids of the physiology which a substance id must not take."""
    names = {
        value
        for name, value in vars(ids).items()
        if name.isupper() and name != "__all__" and isinstance(value, str)
    }
    for organ in Organ:
        names |= {
            ids.organ_volume(organ),
            ids.organ_flow(organ),
            ids.organ_mass(organ),
            ids.flow_fraction(organ),
        }
    for pool in Pool:
        names |= {ids.pool_volume(pool), ids.pool_fraction(pool)}
    return frozenset(names)


@dataclass(frozen=True)
class Variants:
    """The organ variants of a definition, the standard organ everywhere else.

    Attributes:
        organs: the variant of an organ, by name or as component; the names
            are resolved from the variant registry when the model is built.
    """

    organs: Mapping[Organ, str | OrganComponent] = field(default_factory=dict)


@dataclass(frozen=True)
class PBPKDefinition:
    """The definition of a PBPK model.

    Attributes:
        sid: the id of the model.
        substances: the substances.
        chemistry: the processes, each defined once per organ.
        variants: the organ variants.
        routes: the routes of administration.
        pd: the pharmacodynamic components.
        sex: the sex of the model, chosen when it is built: its reference
            values and default individual; a study of both sexes uses one
            model per sex.
        population: the population of the physiology.
        arm_vein: the arm vein sampling site.
        name: the name of the model.
        notes: the notes of the model.
    """

    sid: str
    substances: Sequence[Substance]
    chemistry: Sequence[Metabolism | Transport | CustomReaction] = ()
    variants: Variants = field(default_factory=Variants)
    routes: Sequence[Route] = ()
    pd: Sequence[PDComponent] = ()
    sex: Sex = Sex.MALE
    population: str = "european_icrp"
    arm_vein: str = "pksim_default"
    name: str | None = None
    notes: str | None = None

    def validate(self) -> None:
        """Check the definition.

        Raises:
            ValueError: with all problems found, naming the offending ids.
        """
        problems: list[str] = []
        problems += self._check_substances()
        known = {s.sid for s in self.substances}
        problems += self._check_chemistry(known)
        problems += self._check_routes(known)
        if problems:
            raise ValueError(
                f"Invalid PBPK definition '{self.sid}':\n- " + "\n- ".join(problems)
            )

    def _check_substances(self) -> list[str]:
        """Check the ids and the mechanisms of the substances."""
        problems: list[str] = []
        sids = [s.sid for s in self.substances]
        duplicates = sorted({s for s in sids if sids.count(s) > 1})
        if duplicates:
            problems.append(f"duplicate substance ids: {duplicates}")
        taken = _physiology_ids()
        for s in self.substances:
            if not _SID.fullmatch(s.sid):
                problems.append(f"substance id '{s.sid}' is not an SBML SId")
            if "__" in s.sid:
                problems.append(f"substance id '{s.sid}' contains '__'")
            for token in _RESERVED:
                if s.sid.startswith(f"{token}_"):
                    problems.append(
                        f"substance id '{s.sid}' starts with the reserved token "
                        f"'{token}_' of the id scheme"
                    )
            if s.sid in taken:
                problems.append(f"substance id '{s.sid}' is an id of the physiology")
            problems += self._check_substance(s)
        return problems

    def _check_substance(self, s: Substance) -> list[str]:
        """Check the mechanisms of a substance."""
        problems: list[str] = []
        for inh in s.inhibits:
            if inh.protein not in PROTEINS:
                problems.append(
                    f"{s.sid}: unknown protein '{inh.protein}' (inhibition)"
                )
        for ind in s.induces:
            if ind.protein not in PROTEINS:
                problems.append(f"{s.sid}: unknown protein '{ind.protein}' (induction)")
        problems += _positive(s.sid, "mr", s.mr)
        if s.binding is not None:
            problems += _in_unit_interval(s.sid, "fu", s.binding.fu, open_low=True)
        problems += _in_unit_interval(
            s.sid, "fraction_reabsorbed", s.excretion.fraction_reabsorbed
        )
        ery = s.erythrocytes
        if ery is not None:
            problems += _check_erythrocytes(s.sid, ery)
        dist = s.distribution
        if dist is not None:
            if isinstance(dist, Kp | KpMethod):
                problems += _check_distribution(s.sid, dist)
            if (
                ery is not None
                and ery.mode == "kinetic"
                and dist.mode is Mode.PERFUSION
            ):
                problems.append(
                    f"{s.sid}: kinetic erythrocytes need a distribution with "
                    "exchange, not mode 'perfusion'"
                )
        return problems

    def _check_chemistry(self, known: set[str]) -> list[str]:
        """Check the processes of the chemistry."""
        problems: list[str] = []
        sids = [p.sid for p in self.chemistry]
        duplicates = sorted({s for s in sids if sids.count(s) > 1})
        if duplicates:
            problems.append(f"duplicate process ids: {duplicates}")
        for p in self.chemistry:
            if not _SID.fullmatch(p.sid):
                problems.append(f"process id '{p.sid}' is not an SBML SId")
            used: list[str] = []
            if isinstance(p, Metabolism):
                used = [p.substrate, *p.products]
                for product, fraction in p.products.items():
                    if fraction <= 0:
                        problems.append(
                            f"{p.sid}: product fraction of '{product}' must be "
                            f"positive, got {fraction}"
                        )
                total = sum(p.products.values())
                if p.products and abs(total - 1.0) > 1e-9:
                    problems.append(f"{p.sid}: product fractions sum to {total}, not 1")
            elif isinstance(p, Transport):
                used = [p.substance]
            if isinstance(p, Metabolism | Transport):
                problems += _check_kinetics(p.sid, p.kinetics)
            unknown = sorted(u for u in used if u not in known)
            if unknown:
                problems.append(f"{p.sid}: unknown substances {unknown}")
            protein = getattr(p, "protein", None)
            if protein is not None and protein not in PROTEINS:
                problems.append(f"{p.sid}: unknown protein '{protein}'")
        return problems

    def _check_routes(self, known: set[str]) -> list[str]:
        """Check the routes."""
        return [
            f"route {r.kind}: unknown substance '{r.substance}'"
            for r in self.routes
            if r.substance not in known
        ]


def _check_erythrocytes(sid: str, ery: Erythrocytes) -> list[str]:
    """Check the parameters of the erythrocytes of a substance."""
    problems = [
        *_positive(sid, "bp", ery.bp),
        *_positive(sid, "kbc", ery.kbc),
        *_positive(sid, "k_exchange", ery.k_exchange),
    ]
    if problems:
        return problems
    if ery.mode == "equilibrium":
        if (ery.bp is None) == (ery.kbc is None):
            return [f"{sid}: equilibrium erythrocytes need exactly one of bp and kbc"]
        if ery.k_exchange is not None:
            return [f"{sid}: equilibrium erythrocytes have no k_exchange"]
    elif ery.kbc is None or ery.k_exchange is None:
        return [f"{sid}: kinetic erythrocytes need kbc and k_exchange"]
    elif ery.bp is not None:
        return [f"{sid}: kinetic erythrocytes have no bp"]
    return []


def _check_distribution(sid: str, dist: Kp | KpMethod) -> list[str]:
    """Check the distribution of a substance."""
    problems: list[str] = []
    if dist.mode in {Mode.TISSUE, Mode.PERMEABILITY} and dist.ftissue is None:
        problems.append(f"{sid}: mode '{dist.mode.value}' requires ftissue")
    entered = dist.organs
    named: dict[str, Mapping[Organ, Q]] = {"ps": dist.ps}
    if isinstance(dist, Kp):
        named["override"] = dist.overrides
        problems += _positive(sid, "Kp", dist.value)
        for organ, value in dist.overrides.items():
            problems += _positive(sid, f"Kp override '{organ.value}'", value)
    else:
        problems += _positive(sid, "factor", dist.factor)
    if dist.ftissue is not None:
        problems += _positive(sid, "ftissue", dist.ftissue)
    for organ, value in dist.ps.items():
        problems += _positive(sid, f"ps '{organ.value}'", value)
    if entered is not None:
        for what, mapping in named.items():
            outside = sorted(o.value for o in mapping if o not in entered)
            if outside:
                problems.append(f"{sid}: {what} for organs not entered: {outside}")
    return problems


def _number(q: Q | None) -> float | None:
    """Get the number of a plain float or a `Value`, `None` if there is none."""
    return q.value if isinstance(q, Value) else q


def _positive(where: str, name: str, q: Q | None) -> list[str]:
    """Check that a number is positive."""
    x = _number(q)
    if x is not None and not x > 0:
        return [f"{where}: {name} must be positive, got {x}"]
    return []


def _in_unit_interval(
    where: str, name: str, q: Q | None, *, open_low: bool = False
) -> list[str]:
    """Check that a number is a fraction, in (0, 1] or [0, 1]."""
    x = _number(q)
    if x is None:
        return []
    low_ok = x > 0 if open_low else x >= 0
    if not (low_ok and x <= 1):
        interval = "(0, 1]" if open_low else "[0, 1]"
        return [f"{where}: {name} must be in {interval}, got {x}"]
    return []


def _check_kinetics(where: str, k: FirstOrder | MM | Hill) -> list[str]:
    """Check that the parameters of a rate law are positive."""
    problems: list[str] = []
    for name in ("k", "vmax", "km", "n"):
        problems += _positive(where, name, getattr(k, name, None))
    return problems
