"""The generated steady state of the substances with a declared steady state.

A custom reaction may declare its rate at the steady state
(`CustomReaction.stationary`): a `Formula` of the rate, `"linear"` for a
rate law which is linear in its substrate (a first order clearance), or
`Saturable(vmax, km)` for a Michaelis-Menten rate of a plasma species. A
substance which a declared reaction moves has a generated steady state: the
builder starts every species of it at the steady state of the equations of
the library with the declared rates, so that a model writes only what the
library cannot know (the start of a pool of a custom species, the rates of
its own mechanisms) and every simulation starts stationary.

The steady state is that of the circulation of the standard topology at a
given arterial plasma `Car`:

- an organ takes its inflow (the arterial pool, the venous pool for the
  lung, the mix of the hepatic artery and the portal vein past the shunt for
  the liver) and its declared rates; with the tissue at its steady state
  (`Cti = Kp*(PS*fu*Cpl + r_ti)/(PS*fu + Kp*cl_ti)`) its plasma obeys
  `(R + CLss)*C + S(C) = R*Cin + Rss`, with the flow `R` of the plasma
  leaving it (the kidney: the plasma flow less the net filtrate), the
  clearance `CLss` and the net source `Rss` of its plasma and tissue
  (`steady_clearance`, `steady_rate`, constants of the organ), and a
  saturable rate `S(C) = Vss*C/(Kss + C)`: a closed form, the root of a
  quadratic with `S`. Every other part of the organ starts at the
  equilibrium with the plasma of the organ (or of its zone), a slow pool at
  `Kpslow` times its tissue, the tubule of the nephron kidney at the flow of
  the filtrate (`fu*Car/f_tub_<seg>`), the mucosa of the intestine at the
  equilibrium with the arterial plasma. The starts read the species of the
  inflow, so each is a short formula of its organ, and every parameter which
  an experiment sets before the simulation (an activity of the organ, an
  individual) is seen by them.
- the portal vein starts at the mix of the organs which drain into it, the
  venous pool at the inflow of the lung which gives the arterial plasma (the
  arterial plasma itself when the lung moves nothing), the arterial pool at
  the baseline.

The arterial plasma is the baseline of the substance (`Substance.baseline`):
a substance with a baseline is *pinned*, its declared rates must balance at
the baseline (the model derives the rate which closes the balance, e.g. a
supply, from the baseline). A substance without a baseline is *solved*: the
library solves the mass balance of the body for the arterial plasma and
writes the baseline (`baseline_<s>`): every organ eliminates
`L = lambda*Cin - sigma` of its inflow (`lambda = G + R*CLss/(R + CLss)`,
`sigma = R*Rss/(R + CLss)`), the portal organs lower the inflow of the
liver, the lung eliminates at the arterial plasma, and the sum is 0; a
saturable process of one systemic organ makes it the root of a quadratic in
the plasma of that organ. The organs give their constants to the top model
as outputs (`CLss_<o>_<s>`, ...), so the baseline follows every parameter
of an organ set before the simulation.

Where no closed form exists the build is refused with the reason: a custom
reaction which moves the substance without a declaration, a process of the
library which moves it (but the uptake from the empty lumen), kinetic
erythrocytes, a declared process in an organ part without a closed form (the
full form, the flow-limited organ, the blood cells, a slow pool), a linear
transfer between two parts of the substance, a threshold reabsorption with
another process of the kidney or without a baseline, and of a solved
substance a saturable process outside the systemic organs, two saturable
processes or no elimination.
"""

import warnings
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING

import libsbml
from sbmlutils.reaction_equation import ReactionEquation

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.chemistry import (
    LINEAR,
    CustomReaction,
    Direction,
    Membrane,
    Metabolism,
    Saturable,
    Transport,
)
from pkpdsim.pbpk.organs import PORTAL_ORGANS, SYSTEMIC_ORGANS, Organ, Part, Pool
from pkpdsim.pbpk.parts import parts_for
from pkpdsim.pbpk.placement import ONE, formula_ids
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Formula

if TYPE_CHECKING:
    from sbmlutils.factory import UnitDefinition

    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "BASELINE_SIZE_LIMIT",
    "IMBALANCE_LIMIT",
    "OrganSteady",
    "SteadyPlan",
    "SteadyStateWarning",
    "SteadyValue",
    "carrier",
    "check_built",
    "expanded_size",
    "linear_clearance",
    "organ_exports",
    "plasma_start",
    "pool_start",
    "portal_moves",
    "reduced",
    "saturable_root",
    "solved_baseline",
    "steady_plan",
    "terms",
    "tissue_start",
]

_PARTS = (Part.PLASMA, Part.TISSUE)
"""The parts whose declared processes have a closed steady state."""


@dataclass(frozen=True)
class OrganSteady:
    """The declared steady state of a substance in an organ: what its custom reactions move.

    Formulas of the organ (they read its elements, the parameters of the
    reactions, the species of other substances and constants of the top
    model).

    Attributes:
        sources: the rates into the plasma and the tissue at the steady state
            (mmole/min, signed: a sink is negative), by part.
        clearances: the first order clearances (l/min) of the plasma and the
            tissue, by part.
        saturable: the saturable rate of the plasma, `vmax` with the
            stoichiometry of the reaction.
    """

    sources: Mapping[Part, tuple[str, ...]] = field(default_factory=dict)
    clearances: Mapping[Part, tuple[str, ...]] = field(default_factory=dict)
    saturable: Saturable | None = None

    @property
    def has_clearance(self) -> bool:
        """Whether the plasma has a clearance at the steady state (`CLss`)."""
        return any(self.clearances.values())

    @property
    def has_rate(self) -> bool:
        """Whether the plasma has a net source at the steady state (`Rss`)."""
        return any(self.sources.values())


@dataclass(frozen=True)
class SteadyPlan:
    """The generated steady state of a definition (`steady_plan`).

    Attributes:
        solved: the substances without a baseline, whose arterial plasma
            the library solves and writes as their baseline.
        pinned: the substances with a baseline, the arterial plasma of their
            steady state.
        organs: the declared steady state by organ and substance id, of the
            organs where a declared reaction moves the substance.
    """

    solved: frozenset[str] = frozenset()
    pinned: frozenset[str] = frozenset()
    organs: Mapping[tuple[Organ, str], OrganSteady] = field(default_factory=dict)

    @property
    def substances(self) -> frozenset[str]:
        """The substances with a generated steady state."""
        return self.solved | self.pinned

    def of(self, organ: Organ, sid: str) -> OrganSteady | None:
        """The declared steady state of a substance in an organ, `None` for none."""
        return self.organs.get((organ, sid))


@dataclass(frozen=True)
class SteadyValue:
    """A constant of the top model of the generated steady state.

    Attributes:
        sid: the id.
        formula: the formula over ids of the top model.
        unit: the unit.
        name: the name.
        kind: the kind of its annotation.
        notes: the notes.
    """

    sid: str
    formula: str
    unit: "UnitDefinition"
    name: str
    kind: str
    notes: str


# --- the algebra: formulas from terms ----------------------------------------


def terms(*parts: str) -> str:
    """The sum of the terms which are not `0`, each in parentheses; `0` for none."""
    kept = [t for t in parts if t != "0"]
    return " + ".join(f"({t})" for t in kept) if kept else "0"


def saturable_root(a: str, vmax: str, km: str, d: str) -> str:
    """The positive root `C` of `a*C + vmax*C/(km + C) = d`.

    `a*C^2 + b*C - d*km = 0` with `b = a*km + vmax - d`, in the form
    `2*d*km/(b + sqrt(b^2 + 4*a*d*km))`, which has no cancellation for
    `b >= 0` and no division by `a` (`a` may be 0, then `C = d*km/(vmax - d)`).
    """
    b = f"({a}) * ({km}) + ({vmax}) - ({d})"
    root = f"sqrt(({b})^2 + 4 dimensionless * ({a}) * ({d}) * ({km}))"
    return f"2 dimensionless * ({d}) * ({km}) / ({b} + {root})"


def reduced(steady: OrganSteady, ps_fu: str | None, kp: str | None) -> tuple[str, str]:
    """Get the clearance and the net source of the plasma with the tissue at its steady state.

    The tissue exchanges `PS*fu*(Cpl - Cti/Kp)` with the plasma and has the
    source `r_ti` and the clearance `cl_ti`, so at its steady state it passes
    `r_ti*PS*fu/(PS*fu + Kp*cl_ti)` to the plasma and clears
    `cl_ti*Kp*PS*fu/(PS*fu + Kp*cl_ti)` of it.

    Args:
        steady: the declared steady state of the substance in the organ.
        ps_fu: the exchange `PS*fu` of the merged tissue, `None` without one.
        kp: the tissue to plasma coefficient, `None` without a tissue.

    Returns:
        The clearance `CLss` (l/min) and the net source `Rss` (mmole/min)
        of the plasma, `0` for none.
    """
    r_pl = terms(*steady.sources.get(Part.PLASMA, ()))
    r_ti = terms(*steady.sources.get(Part.TISSUE, ()))
    cl_pl = terms(*steady.clearances.get(Part.PLASMA, ()))
    cl_ti = terms(*steady.clearances.get(Part.TISSUE, ()))
    if cl_ti == "0":
        return cl_pl, terms(r_pl, r_ti)
    if ps_fu is None or kp is None:
        raise ValueError(
            "A declared clearance of the tissue needs the merged tissue (its "
            "exchange and Kp)."
        )
    # the exchange and the clearance of the tissue in series, each read once
    # (roadrunner inlines every initial assignment a start reads)
    cleared = f"{ONE} / ({ONE} / ({kp} * ({cl_ti})) + {ONE} / ({ps_fu}))"
    passed = (
        "0" if r_ti == "0" else f"({r_ti}) / ({ONE} + {kp} * ({cl_ti}) / ({ps_fu}))"
    )
    return terms(cl_pl, cleared), terms(r_pl, passed)


def plasma_start(
    flow: str,
    inflow: str,
    clearance: str,
    rate: str,
    saturable: tuple[str, str] | None = None,
) -> str:
    """Get the plasma of an organ at its steady state [mM].

    `(R + CLss)*C + S(C) = R*Cin + Rss`: `C = (R*Cin + Rss)/(R + CLss)`, the
    root of a quadratic with a saturable rate `S(C) = vmax*C/(km + C)`, the
    inflow itself without any process.

    Args:
        flow: the flow `R` of the plasma leaving the organ (l/min).
        inflow: the inflow concentration `Cin`.
        clearance: the clearance `CLss` (l/min), `0` for none.
        rate: the net source `Rss` (mmole/min), `0` for none.
        saturable: `vmax` and `km` of a saturable rate of the plasma.

    Returns:
        The formula.
    """
    load = terms(f"{flow} * {inflow}", rate)
    linear = terms(flow, clearance)
    if saturable is not None:
        return saturable_root(linear, saturable[0], saturable[1], load)
    if clearance == "0" and rate == "0":
        return inflow
    return f"({load}) / ({linear})"


def tissue_start(
    kp: str, plasma: str, ps_fu: str | None, rate: str, clearance: str
) -> str:
    """Get the merged tissue of an organ at its steady state [mM].

    `Cti = Kp*(PS*fu*Cpl + r_ti)/(PS*fu + Kp*cl_ti)`, `Kp*Cpl` without a
    process of the tissue.

    Args:
        kp: the tissue to plasma coefficient.
        plasma: the plasma of the organ (its species).
        ps_fu: the exchange `PS*fu`.
        rate: the net source of the tissue (mmole/min), `0` for none.
        clearance: the clearance of the tissue (l/min), `0` for none.

    Returns:
        The formula.
    """
    if rate == "0" and clearance == "0":
        return f"{kp} * {plasma}"
    if ps_fu is None:
        raise ValueError(
            "A declared process of the tissue needs the merged tissue (its exchange)."
        )
    if clearance == "0":
        return f"{kp} * ({plasma} + ({rate}) / ({ps_fu}))"
    passed = terms(f"{ps_fu} * {plasma}", rate)
    return f"{kp} * ({passed}) / ({ps_fu} + {kp} * ({clearance}))"


# --- the linear rate law -------------------------------------------------------


def _path(node: libsbml.ASTNode, name: str) -> list[list[int]]:
    """The paths of child indices from a node to every occurrence of a name."""
    if node.getType() == libsbml.AST_NAME and node.getName() == name:
        return [[]]
    found: list[list[int]] = []
    for k in range(node.getNumChildren()):
        found += [[k, *p] for p in _path(node.getChild(k), name)]
    return found


def _linear_path(node: libsbml.ASTNode, path: Sequence[int]) -> bool:
    """Whether the ancestors of an occurrence keep the formula linear in it.

    A product and the numerator of a quotient keep it linear; a negation
    is refused (a negative clearance is a production, `linear` is a
    clearance of its substrate).
    """
    for k in path:
        kind = node.getType()
        if kind == libsbml.AST_TIMES:
            pass
        elif kind == libsbml.AST_DIVIDE:
            if k != 0:
                return False
        else:
            return False
        node = node.getChild(k)
    return True


def linear_clearance(formula: str, species: str) -> str | None:
    """Get the clearance of a rate law which is linear in a species.

    The rate law `k*C` (a product with the species as one factor, or the
    numerator of a quotient) with the species `C` once: the
    clearance is the rate law without the factor, `k` (l/min for a rate in
    mmole/min of a concentration in mM).

    Args:
        formula: the rate law (SBML L3 infix).
        species: the id of the species.

    Returns:
        The clearance, `None` if the rate law is not linear in the species
        (it reads it more than once, in a sum, a denominator, a function or
        a negation, or not at all).
    """
    ast = libsbml.parseL3Formula(formula)
    if ast is None:
        return None
    paths = _path(ast, species)
    if len(paths) != 1 or not paths[0] or not _linear_path(ast, paths[0]):
        return None
    path = paths[0]
    parent = ast
    for k in path[:-1]:
        parent = parent.getChild(k)
    last = path[-1]
    if parent.getType() == libsbml.AST_TIMES and parent.getNumChildren() > 2:
        parent.removeChild(last)
    elif parent.getType() == libsbml.AST_TIMES:
        other = parent.getChild(1 - last).deepCopy()
        if len(path) == 1:
            ast = other
        else:
            grand = ast
            for k in path[:-2]:
                grand = grand.getChild(k)
            grand.replaceChild(path[-2], other, True)
    else:
        one = libsbml.ASTNode(libsbml.AST_REAL)
        one.setValue(1.0)
        one.setUnits("dimensionless")
        parent.replaceChild(last, one, True)
    return libsbml.formulaToL3String(ast)


# --- the plan --------------------------------------------------------------------


def _species_map(
    substances: Sequence[Substance],
) -> dict[str, tuple[str, Part | Pool | None]]:
    """The species of the substances by local id: of an organ part, a pool, a slow pool."""
    found: dict[str, tuple[str, Part | Pool | None]] = {}
    for s in substances:
        for part in Part:
            found[ids.species(part, s.sid)] = (s.sid, part)
        for pool in Pool:
            found[ids.pool_species(pool, s.sid)] = (s.sid, pool)
            found[ids.pool_species(pool, s.sid, Part.BLOOD_CELLS)] = (s.sid, pool)
        found[ids.slow_species(s.sid)] = (s.sid, None)
    return found


def _moved(
    p: CustomReaction, species: Mapping[str, tuple[str, Part | Pool | None]]
) -> list[tuple[str, float, str, Part | Pool | None]]:
    """The species of substances which a custom reaction moves: (id, signed stoichiometry, substance, place)."""
    eq = ReactionEquation.from_str(p.equation)
    moved: list[tuple[str, float, str, Part | Pool | None]] = []
    for sign, items in ((-1.0, eq.reactants), (1.0, eq.products)):
        for item in items:
            sid = item.species
            if sid in species:
                sub, place = species[sid]
                moved.append((sid, sign * float(item.stoichiometry or 1.0), sub, place))
    return moved


def _times(stoichiometry: float, formula: str) -> str:
    """A formula times a signed stoichiometry."""
    factor = abs(stoichiometry)
    scaled = formula if factor == 1.0 else f"{factor:g} dimensionless * ({formula})"
    return scaled if stoichiometry > 0 else f"-({scaled})"


class _Collector:
    """The contributions of the declared reactions, by organ and substance."""

    def __init__(self) -> None:
        self.sources: dict[tuple[Organ, str, Part], list[str]] = {}
        self.clearances: dict[tuple[Organ, str, Part], list[str]] = {}
        self.saturable: dict[tuple[Organ, str], list[Saturable]] = {}
        self.touched: set[tuple[Organ, str]] = set()

    def source(self, organ: Organ, sid: str, part: Part, formula: str) -> None:
        self.sources.setdefault((organ, sid, part), []).append(formula)
        self.touched.add((organ, sid))

    def clearance(self, organ: Organ, sid: str, part: Part, formula: str) -> None:
        self.clearances.setdefault((organ, sid, part), []).append(formula)
        self.touched.add((organ, sid))

    def organs(self) -> dict[tuple[Organ, str], OrganSteady]:
        result: dict[tuple[Organ, str], OrganSteady] = {}
        for organ, sid in sorted(self.touched):
            saturable = self.saturable.get((organ, sid), [])
            result[(organ, sid)] = OrganSteady(
                sources={
                    part: tuple(self.sources[(organ, sid, part)])
                    for part in _PARTS
                    if (organ, sid, part) in self.sources
                },
                clearances={
                    part: tuple(self.clearances[(organ, sid, part)])
                    for part in _PARTS
                    if (organ, sid, part) in self.clearances
                },
                saturable=saturable[0] if saturable else None,
            )
        return result


def steady_plan(definition: "PBPKDefinition") -> SteadyPlan:
    """Get the generated steady state of a definition, refusing what has no closed form.

    Args:
        definition: the definition.

    Returns:
        The plan, empty without a declared steady state.

    Raises:
        ValueError: naming every reason why the steady state of a substance
            has no closed form (see the module).
    """
    declared = [
        p
        for p in definition.chemistry
        if isinstance(p, CustomReaction) and p.stationary is not None
    ]
    if not declared:
        return SteadyPlan()
    by_sid = {s.sid: s for s in definition.substances}
    species = _species_map(definition.substances)
    custom = [p for p in definition.chemistry if isinstance(p, CustomReaction)]
    moves = {p.sid: _moved(p, species) for p in custom}
    generated = sorted({sub for p in declared for _, _, sub, _ in moves[p.sid]})
    problems: list[str] = []
    collector = _Collector()
    depends: dict[str, set[str]] = {sid: set() for sid in generated}
    for p in custom:
        mine = sorted({sub for _, _, sub, _ in moves[p.sid] if sub in depends})
        if mine and p.stationary is None:
            problems.append(
                f"custom reaction '{p.sid}' moves {mine}, which have a generated "
                "steady state: declare its rate at the steady state "
                "(`stationary=`: a Formula, 'linear' or Saturable)"
            )
        elif p.stationary is not None:
            problems += _declared(p, moves[p.sid], species, collector, depends)
    problems += _library_processes(definition, set(generated))
    for sid in generated:
        problems += _substance_problems(by_sid[sid], collector)
    organs = collector.organs()
    problems += _part_problems(definition, organs, by_sid)
    solved = frozenset(sid for sid in generated if by_sid[sid].baseline is None)
    for sid in sorted(solved):
        problems += _solved_problems(by_sid[sid], organs)
    problems += _cycles(depends)
    if problems:
        raise ValueError(
            f"Invalid PBPK definition '{definition.sid}': no closed steady "
            "state:\n- " + "\n- ".join(problems)
        )
    return SteadyPlan(
        solved=solved,
        pinned=frozenset(generated) - solved,
        organs=organs,
    )


def _declared(
    p: CustomReaction,
    moved: Sequence[tuple[str, float, str, Part | Pool | None]],
    species: Mapping[str, tuple[str, Part | Pool | None]],
    collector: _Collector,
    depends: dict[str, set[str]],
) -> list[str]:
    """Collect the contributions of a declared reaction, with its problems."""
    problems: list[str] = []
    where = f"custom reaction '{p.sid}'"
    subs = {sub for _, _, sub, _ in moved}
    for sid, _, _, place in moved:
        if isinstance(place, Pool):
            problems.append(
                f"{where} moves the species '{sid}' of a blood pool: a declared "
                "steady state acts in the parts of its organ"
            )
        elif place not in _PARTS:
            what = (
                "slow pool" if place is None else place.name.lower().replace("_", " ")
            )
            problems.append(
                f"{where} moves '{sid}' of the {what}: a declared steady state has a "
                "closed form in the plasma and the merged tissue of an organ only"
            )
    if problems:
        return problems
    parts = [
        (sid, nu, sub, place)
        for sid, nu, sub, place in moved
        if isinstance(place, Part)
    ]

    def reads(expression: str, what: str) -> set[str]:
        """The substances whose species an expression reads, refusing the moved ones."""
        read = {species[x][0] for x in formula_ids(expression) if x in species}
        own = sorted(read & subs)
        if own:
            problems.append(
                f"{where}: the {what} reads species of {own}, which the reaction "
                "moves; it is the rate at their steady state, a formula of constants "
                "and of the starts of other species"
            )
        return read - subs

    s = p.stationary
    if isinstance(s, Formula):
        others = reads(s.expression, "stationary formula")
        for _, nu, sub, place in parts:
            collector.source(p.organ, sub, place, _times(nu, s.expression))
            depends[sub] |= others
        return problems
    if s != LINEAR and not isinstance(s, Saturable):
        return [
            f"{where}: stationary is a {type(s).__name__}, not a Formula, 'linear' "
            "or a Saturable"
        ]
    reactants = [m for m in parts if m[1] < 0]
    if len(reactants) != 1:
        problems.append(
            f"{where}: a {'linear' if s == LINEAR else 'saturable'} steady state "
            f"needs one substrate of a substance, the reaction consumes {len(reactants)}"
        )
        return problems
    substrate, nu, sub, place = reactants[0]
    products = [m for m in parts if m[1] > 0]
    for sid, _, other, _ in products:
        if other == sub:
            problems.append(
                f"{where} moves '{sub}' from '{substrate}' to '{sid}': a linear or "
                "saturable transfer between two parts of a substance has no closed "
                "steady state, declare its rate as a Formula"
            )
    if problems:
        return problems
    if s == LINEAR:
        clearance = linear_clearance(p.formula, substrate)
        if clearance is None:
            return [
                f"{where}: the rate law `{p.formula}` is not linear in its substrate "
                f"'{substrate}' (one factor of a product or the numerator of a "
                "quotient), declare its rate as a Formula or a Saturable"
            ]
        others = reads(clearance, "clearance of the linear rate law")
        collector.clearance(p.organ, sub, place, _times(-nu, clearance))
        rate = f"({clearance}) * {substrate}"
    else:
        if place is not Part.PLASMA:
            return [
                f"{where}: a saturable steady state is a rate of the plasma, the "
                f"substrate '{substrate}' is in the {place.name.lower()}"
            ]
        if collector.saturable.get((p.organ, sub)):
            return [
                f"{where}: a second saturable rate of '{sub}' in the organ "
                f"'{p.organ.name.lower()}', its plasma has a closed steady state "
                "with one"
            ]
        others = reads(s.vmax, "vmax of the saturable rate") | reads(
            s.km, "km of the saturable rate"
        )
        collector.saturable.setdefault((p.organ, sub), []).append(
            Saturable(vmax=_times(-nu, s.vmax), km=s.km)
        )
        collector.touched.add((p.organ, sub))
        rate = f"({s.vmax}) * {substrate} / ({s.km} + {substrate})"
    depends[sub] |= others
    for _, nu_p, other, part in products:
        collector.source(p.organ, other, part, _times(nu_p, rate))
        depends[other] |= {sub} | others
    return problems


def _library_processes(definition: "PBPKDefinition", generated: set[str]) -> list[str]:
    """Refuse a metabolism or transport of the library which moves a generated substance.

    The uptake from the lumen of the small intestine (an apical influx) is
    kept: the lumen of an endogenous substance is empty at the start.
    """
    problems: list[str] = []
    for p in definition.chemistry:
        if isinstance(p, Metabolism):
            moved = sorted({p.substrate, *p.products} & generated)
        elif isinstance(p, Transport):
            if (
                p.organ is Organ.SMALL_INTESTINE
                and p.membrane is Membrane.APICAL
                and p.direction is Direction.INFLUX
            ):
                continue
            moved = sorted({p.substance} & generated)
        else:
            continue
        if moved:
            problems.append(
                f"the {type(p).__name__.lower()} '{p.sid}' of the library moves "
                f"{moved}, which have a generated steady state: a process of the "
                "library has no declared steady state (a custom reaction with "
                "`stationary=` has)"
            )
    return problems


def _substance_problems(s: Substance, collector: _Collector) -> list[str]:
    """Refuse the properties of a substance without a closed steady state."""
    problems: list[str] = []
    if s.erythrocytes is not None and s.erythrocytes.mode == "kinetic":
        problems.append(
            f"{s.sid}: kinetic erythrocytes exchange with the plasma of every organ, "
            "no closed steady state (erythrocytes in equilibrium have one)"
        )
    if (
        s.excretion.reabsorption is not None
        and (Organ.KIDNEY, s.sid) in collector.touched
    ):
        problems.append(
            f"{s.sid}: a threshold reabsorption with a declared process of the "
            "kidney: its net filtrate depends on the concentration, the kidney "
            "plasma has no closed steady state"
        )
    return problems


def _part_problems(
    definition: "PBPKDefinition",
    organs: Mapping[tuple[Organ, str], OrganSteady],
    by_sid: Mapping[str, Substance],
) -> list[str]:
    """Refuse a declared process in an organ whose parts have no closed steady state."""
    problems: list[str] = []
    chemistry = tuple(definition.chemistry)
    for (organ, sid), _ in organs.items():
        parts = parts_for(by_sid[sid], organ, chemistry)
        if parts.tissue in ("full", "organ"):
            form = "full form" if parts.tissue == "full" else "flow-limited organ"
            problems.append(
                f"{sid}: a declared process of the organ '{organ.name.lower()}', "
                f"whose tissue is the {form}: a closed steady state has the merged "
                "tissue (Mode.TISSUE) or none"
            )
    return problems


def _solved_problems(
    s: Substance, organs: Mapping[tuple[Organ, str], OrganSteady]
) -> list[str]:
    """Refuse what the solve of the arterial plasma of a substance without baseline cannot do."""
    problems: list[str] = []
    where = f"{s.sid} (no baseline, the library solves its arterial plasma)"
    if s.excretion.reabsorption is not None:
        problems.append(
            f"{where}: a threshold reabsorption makes the mass balance nonlinear "
            "in the arterial plasma; give the substance a baseline"
        )
    mine = {organ: steady for (organ, sid), steady in organs.items() if sid == s.sid}
    saturable = sorted(o.name.lower() for o, st in mine.items() if st.saturable)
    if len(saturable) > 1:
        problems.append(
            f"{where}: saturable processes in {saturable}, the solve has a closed "
            "form with one; give the substance a baseline"
        )
    outside = sorted(
        o.name.lower()
        for o, st in mine.items()
        if st.saturable is not None and o not in SYSTEMIC_ORGANS
    )
    if outside:
        problems.append(
            f"{where}: a saturable process of {outside}, the solve has a closed "
            "form for a saturable process of a systemic organ (not the lung, the "
            "liver or a portal organ); give the substance a baseline"
        )
    eliminated = s.excretion.urine or any(
        st.has_clearance or st.saturable is not None for st in mine.values()
    )
    if not eliminated:
        problems.append(
            f"{where}: nothing eliminates it (no filtration, no linear or saturable "
            "clearance), it has no steady state"
        )
    return problems


def _cycles(depends: Mapping[str, set[str]]) -> list[str]:
    """Refuse substances whose steady states read each other."""
    problems: list[str] = []
    for sid in sorted(depends):
        seen: set[str] = set()
        stack = sorted(depends[sid])
        while stack:
            x = stack.pop()
            if x == sid:
                problems.append(
                    f"{sid}: its steady state reads its own through the declared "
                    "reactions of other substances (a cycle of products)"
                )
                break
            if x in seen:
                continue
            seen.add(x)
            stack += sorted(depends.get(x, ()))
    return problems


# --- the top model: exports, the solved baseline, the pools --------------------


def organ_exports(plan: SteadyPlan, organ: Organ, s: str) -> dict[str, str]:
    """Get the constants of the steady state which an organ gives to the top model.

    The clearance `CLss_<s>`, the net source `Rss_<s>` and the saturable
    `Vss_<s>`, `Kss_<s>` of the plasma of the organ, which it has, as
    `CLss_<o>_<s>`, ...: of every organ of a solved substance (the solve of
    the arterial plasma), of the lung of every substance (the venous pool).

    Args:
        plan: the plan.
        organ: the organ.
        s: the substance id.

    Returns:
        The outputs, local id -> id of the top model.
    """
    steady = plan.of(organ, s)
    if steady is None or not (s in plan.solved or organ is Organ.LUNG):
        return {}
    found: dict[str, str] = {}
    if steady.has_clearance:
        found[ids.steady_clearance(s)] = ids.steady_clearance(s, organ)
    if steady.has_rate:
        found[ids.steady_rate(s)] = ids.steady_rate(s, organ)
    if steady.saturable is not None:
        found[ids.steady_vmax(s)] = ids.steady_vmax(s, organ)
        found[ids.steady_km(s)] = ids.steady_km(s, organ)
    return found


def carrier(s: Substance) -> str:
    """The share of the blood flow which carries the plasma concentration, in the top model."""
    if s.erythrocytes is not None and s.erythrocytes.mode == "equilibrium":
        return ids.bp(s.sid)
    return f"({ONE} - {ids.HCT})"


def _filtered(s: Substance) -> str:
    """The net filtered plasma flow of a constant reabsorbed fraction, `0` without filtration."""
    if not s.excretion.urine:
        return "0"
    return f"{ids.fu(s.sid)} * {ids.GFR} * ({ONE} - {ids.fraction_reabsorbed(s.sid)})"


def _or_zero(formula: str, unit: "UnitDefinition") -> str:
    """A formula, `0` with its unit for none."""
    return f"0 {getattr(unit, 'sid', unit)}" if formula == "0" else formula


@dataclass(frozen=True)
class _Organ:
    """An organ in the solve: inflow flow `F`, net filtrate `G`, its constants."""

    flow: str
    filtrate: str
    clearance: str
    rate: str
    vmax: str | None
    km: str | None

    @property
    def leaving(self) -> str:
        """The flow `R = F - G` of the plasma leaving the organ."""
        return self.flow if self.filtrate == "0" else f"({self.flow} - {self.filtrate})"

    def eliminated(self) -> tuple[str, str]:
        """`lambda` and `sigma` of the elimination `L = lambda*Cin - sigma` of the organ.

        `lambda = G + 1/(1/CLss + 1/R)`, `sigma = Rss/(1 + CLss/R)`: the
        forms of `G + R*CLss/(R + CLss)` and `R*Rss/(R + CLss)` which read
        the clearance and the flow once (roadrunner inlines every initial assignment it reads,
        the load time grows with the copies of the flows).
        """
        if self.clearance == "0":
            return self.filtrate, self.rate
        cleared = f"{ONE} / ({ONE} / {self.clearance} + {ONE} / ({self.leaving}))"
        released = (
            "0"
            if self.rate == "0"
            else f"{self.rate} / ({ONE} + {self.clearance} / ({self.leaving}))"
        )
        return terms(self.filtrate, cleared), released


def _solve_organ(
    plan: SteadyPlan, s: Substance, organ: Organ, flow: str, filtrate: str = "0"
) -> _Organ | None:
    """The organ in the solve of a substance, `None` if it eliminates nothing."""
    steady = plan.of(organ, s.sid)
    if steady is None and filtrate == "0":
        return None
    exports = organ_exports(plan, organ, s.sid) if steady is not None else {}
    return _Organ(
        flow=flow,
        filtrate=filtrate,
        clearance=exports.get(ids.steady_clearance(s.sid), "0"),
        rate=exports.get(ids.steady_rate(s.sid), "0"),
        vmax=exports.get(ids.steady_vmax(s.sid)),
        km=exports.get(ids.steady_km(s.sid)),
    )


def solved_baseline(
    s: Substance, plan: SteadyPlan, organs: Iterable[Organ]
) -> tuple[list[SteadyValue], str]:
    """Get the arterial plasma of the steady state of a solved substance, its baseline.

    The sum of the eliminations of the organs is 0:
    `sum lambda_o*Cin_o - sigma_o = 0`, the inflow of the liver lowered by
    the eliminations of the portal organs (`Cin_li = Car - sum L_p/(Qh*c)`),
    the lung at the arterial plasma; with a saturable process of a systemic
    organ `k` the rest of the body gives `Car = (Rss + R_k*C_k)/(CLss + F_k)`
    and `C_k` is the root of `a*C + q*S(C) = d` with `q = (CLss + F_k)/R_k`,
    `a = CLss + G_k + q*CLss_k` and `d = Rss + q*Rss_k`.

    Args:
        s: the substance, without baseline.
        plan: the plan.
        organs: the built organs (the intestine of a variant has no process
            of the substance).

    Returns:
        The constants of the top model (`CLss_<s>`, `Rss_<s>`, `Css_<o>_<s>`)
        and the formula of the baseline.
    """
    sid = s.sid
    c = carrier(s)
    built = set(organs)
    systemic: dict[Organ, _Organ] = {}
    portal: list[_Organ] = []
    liver: _Organ | None = None
    lung: _Organ | None = None
    for organ in Organ:
        if organ not in built:
            continue
        if organ is Organ.LIVER:
            flow = f"({ONE} - {ids.F_SHUNT}) * {ids.QH} * {c}"
            liver = _solve_organ(plan, s, organ, flow)
            continue
        flow = f"{ids.organ_flow(organ)} * {c}"
        if organ is Organ.LUNG:
            lung = _solve_organ(plan, s, organ, flow)
        elif organ in PORTAL_ORGANS:
            if (o := _solve_organ(plan, s, organ, flow)) is not None:
                portal.append(o)
        elif organ in SYSTEMIC_ORGANS:
            filtrate = _filtered(s) if organ is Organ.KIDNEY else "0"
            if (o := _solve_organ(plan, s, organ, flow, filtrate)) is not None:
                systemic[organ] = o
    saturable = next((o for o, x in systemic.items() if x.vmax is not None), None)
    k = systemic.pop(saturable) if saturable is not None else None
    lam: list[str] = []
    sig: list[str] = []
    for o in systemic.values():
        lam_o, sig_o = o.eliminated()
        lam.append(lam_o)
        sig.append(sig_o)
    lam_p, sig_p = "0", "0"
    if portal:
        pairs = [o.eliminated() for o in portal]
        lam_p = terms(*(x for x, _ in pairs))
        sig_p = terms(*(y for _, y in pairs))
        lam += [lam_p]
        sig += [sig_p]
    if liver is not None:
        lam_li, sig_li = liver.eliminated()
        mixed = f"{ids.QH} * {c}"
        if lam_p != "0":
            lam.append(f"({lam_li}) * ({ONE} - ({lam_p}) / ({mixed}))")
        else:
            lam.append(lam_li)
        sig.append(sig_li)
        if sig_p != "0":
            sig.append(f"-({lam_li}) * ({sig_p}) / ({mixed})")
    if lung is not None:
        # the lung plasma is the arterial plasma: it eliminates CLss*Car - Rss
        lam.append(lung.clearance)
        sig.append(lung.rate)
    clearance = _or_zero(terms(*lam), U.l_per_min)
    rate = _or_zero(terms(*sig), U.mmole_per_min)
    cl_id, rate_id = ids.steady_clearance(sid), ids.steady_rate(sid)
    rest = " besides the organ of the saturable process" if k is not None else ""
    values = [
        SteadyValue(
            cl_id,
            clearance,
            U.l_per_min,
            f"clearance of the arterial plasma at the steady state {sid} [l/min]",
            "steady_clearance",
            "The clearance of the arterial plasma by the body at the steady state"
            f"{rest}: the eliminations of the organs, the portal organs before "
            "the liver (`pkpdsim.pbpk.stationary.solved_baseline`).",
        ),
        SteadyValue(
            rate_id,
            rate,
            U.mmole_per_min,
            f"net source at the steady state {sid} [mmole/min]",
            "steady_rate",
            f"The net source of the body at the steady state{rest}: the declared "
            "rates which leave the organs (`pkpdsim.pbpk.stationary.solved_baseline`).",
        ),
    ]
    if k is None:
        return values, f"{rate_id} / {cl_id}"
    if saturable is None or k.vmax is None or k.km is None:
        raise RuntimeError(f"{sid}: the saturable organ of the solve has no constants.")
    r = k.leaving
    total = f"{cl_id} + {k.flow}"
    # `a*C + S(C) = d` of the organ, times `(CLss + F_k)/R_k`: the form which
    # reads the flows least (without other clearances that of a production
    # eliminated by the organ, `G*C + F/R*S(C) = Rss`)
    q = f"({total}) / ({r})"
    a = terms(cl_id, k.filtrate, "0" if k.clearance == "0" else f"{k.clearance} * {q}")
    d = terms(rate_id, "0" if k.rate == "0" else f"{k.rate} * {q}")
    css = ids.steady_plasma(saturable, sid)
    values.append(
        SteadyValue(
            css,
            saturable_root(a, f"{k.vmax} * {q}", k.km, d),
            U.mM,
            f"plasma of the organ {saturable.name.lower()} at the steady state "
            f"{sid} [mM]",
            "steady_plasma",
            "The plasma of the organ of the saturable process at the steady state: "
            "the positive root of the quadratic of its balance with the rest of the "
            "body (`pkpdsim.pbpk.stationary.solved_baseline`).",
        )
    )
    return values, f"({rate_id} + {r} * {css}) / ({total})"


def portal_moves(plan: SteadyPlan, s: str, gut: Iterable[Organ] = ()) -> bool:
    """Whether an organ which drains into the portal vein moves the substance.

    Without, the portal vein starts at the arterial plasma, and so does the
    inflow of the liver.

    Args:
        plan: the plan.
        s: the substance id.
        gut: the organs which the intestine of a variant replaces (it moves
            nothing).
    """
    replaced = set(gut)
    return any(
        plan.of(o, s) is not None
        for o in Organ
        if o in PORTAL_ORGANS and o not in replaced
    )


def pool_start(
    s: Substance,
    plan: SteadyPlan,
    pool: Pool,
    organs: Iterable[Organ],
    gut: Iterable[Organ] = (),
) -> str | None:
    """Get the start of a blood pool of a substance with a generated steady state.

    The arterial pool at the baseline (`None`, the start of the library);
    the portal vein at the mix of the organs which drain into it (the organs
    of the intestine of a variant at the arterial plasma, it moves nothing),
    the arterial plasma when no portal organ moves the substance; the venous
    pool at the inflow of the lung whose plasma is the arterial plasma,
    `Car + (CLss*Car + S(Car) - Rss)/F` with the constants of the lung, the
    arterial plasma when the lung moves nothing.

    Args:
        s: the substance.
        plan: the plan.
        pool: the pool.
        organs: the built organs.
        gut: the organs which the intestine of a variant replaces, the small
            intestine included.

    Returns:
        The formula, `None` for the start of the library.
    """
    sid = s.sid
    car = ids.pool_species(Pool.ARTERIAL, sid)
    if pool is Pool.ARTERIAL:
        return None
    if pool is Pool.PORTAL:
        replaced = set(gut)
        draining = [o for o in Organ if o in PORTAL_ORGANS]
        if not portal_moves(plan, sid, replaced):
            return car
        built = set(organs)
        moving = [
            o
            for o in draining
            if o not in replaced and o in built and plan.of(o, sid) is not None
        ]
        # the arterial plasma read once, by the flows of the organs which move
        # nothing of the substance
        rest = " + ".join(ids.organ_flow(o) for o in draining if o not in moving)
        mix = [f"({rest}) * {car}"] if rest else []
        mix += [f"{ids.organ_flow(o)} * {ids.organ_outflow(o, sid)}" for o in moving]
        return f"({' + '.join(mix)}) / {ids.QPO}"
    lung = plan.of(Organ.LUNG, sid)
    if lung is None:
        return car
    exports = organ_exports(plan, Organ.LUNG, sid)
    removed = []
    if (cl := exports.get(ids.steady_clearance(sid))) is not None:
        removed.append(f"{cl} * {car}")
    if (vmax := exports.get(ids.steady_vmax(sid))) is not None:
        km = exports[ids.steady_km(sid)]
        removed.append(f"{vmax} * {car} / ({km} + {car})")
    if (rate := exports.get(ids.steady_rate(sid))) is not None:
        removed.append(f"-({rate})")
    flow = f"{ids.organ_flow(Organ.LUNG)} * {carrier(s)}"
    return f"{car} + ({terms(*removed)}) / ({flow})"


# --- the checks of the built model ---------------------------------------------

BASELINE_SIZE_LIMIT = 20_000
"""The expanded size (nodes) of a solved baseline above which the build warns.

roadrunner inlines every initial assignment a start reads, without
memoizing it: every species of the substance compiles a copy of the
expanded baseline (a parameter in between is inlined as well, it is no
cheaper). The baseline of creatinine expands to about 3,800 nodes (the model
about 220,000, 4 s at LLJIT O0); a solved saturable process of a body which
also eliminates in the liver and a portal organ reached 31,000 nodes and
1.2 million in the model (37 s at O0)."""

IMBALANCE_LIMIT = 1e-9
"""The relative rate (1/min) of the venous pool of a pinned substance at the
start above which the build warns that its declared rates do not balance at
its baseline; a balanced start is stationary to the rounding (about 1e-16)."""


class SteadyStateWarning(UserWarning):
    """A generated steady state which loads slowly or does not balance."""


def expanded_size(model: libsbml.Model, sid: str) -> int:
    """Get the size of the expanded initial value of an id (nodes of its math tree).

    Every id the math reads which has an initial assignment or an assignment
    rule is replaced by its own expanded math, as roadrunner compiles it;
    every other id and number counts one.

    Args:
        model: the flat model.
        sid: the id.

    Returns:
        The number of nodes.
    """
    maths = {a.getSymbol(): a.getMath() for a in model.getListOfInitialAssignments()}
    maths |= {
        r.getVariable(): r.getMath()
        for r in model.getListOfRules()
        if r.isAssignment() and r.getVariable() not in maths
    }

    @cache
    def size(x: str) -> int:
        ast = maths.get(x)
        if ast is None:
            return 1
        total, stack = 0, [ast]
        while stack:
            node = stack.pop()
            total += size(node.getName()) if node.getType() == libsbml.AST_NAME else 1
            stack.extend(node.getChild(k) for k in range(node.getNumChildren()))
        return total

    return size(sid)


def check_built(plan: SteadyPlan, flat_path: Path) -> None:
    """Warn about a generated steady state of a built model which loads slowly or does not balance.

    A solved baseline whose expansion exceeds `BASELINE_SIZE_LIMIT`: every
    start reads it, give the substance a baseline (pinned) or fewer declared
    organs. A pinned substance whose venous pool moves at the start of the
    reference individual (`IMBALANCE_LIMIT`): its declared rates do not
    balance at its baseline by the net rate of the venous pool (mmole/min),
    the model derives the rate which closes the balance (e.g. a supply) from
    the baseline (the formulas of `pkpdsim.pbpk.steady_state`).

    Args:
        plan: the plan of the definition.
        flat_path: the flat model.

    Warns:
        SteadyStateWarning: naming the substance and the size or the imbalance.
    """
    if not plan.substances:
        return
    model = libsbml.readSBMLFromFile(str(flat_path)).getModel()
    for sid in sorted(plan.solved):
        n = expanded_size(model, ids.baseline(sid))
        if n > BASELINE_SIZE_LIMIT:
            warnings.warn(
                f"{sid}: the solved baseline expands to {n} nodes (above "
                f"{BASELINE_SIZE_LIMIT}); roadrunner compiles a copy into every "
                "start of the substance, the model loads slowly. Give the "
                "substance a baseline (pinned) or declare fewer organs.",
                SteadyStateWarning,
                stacklevel=3,
            )
    if not plan.pinned:
        return
    # a factory without pinned substances loads no simulator
    import roadrunner

    options = roadrunner.LoadSBMLOptions()
    options.setLLVMBackend(roadrunner.LoadSBMLOptions.LLJIT)
    options.setLLJitOptimizationLevel(roadrunner.LoadSBMLOptions.NONE)
    r = roadrunner.RoadRunner(str(flat_path), options)
    rates = list(r.getRatesOfChange())
    state = {r.model.getStateVectorId(i): rates[i] for i in range(len(rates))}
    for sid in sorted(plan.pinned):
        venous = ids.pool_species(Pool.VENOUS, sid)
        amount = r.model[venous]
        rate = state.get(venous, 0.0)
        relative = abs(rate) / amount if amount > 0 else abs(rate)
        if relative > IMBALANCE_LIMIT:
            warnings.warn(
                f"{sid}: the declared rates do not balance at the baseline: the "
                f"venous pool changes by {rate:.4g} mmole/min ({relative:.3g} of "
                "its amount per minute) at the start of the reference individual. "
                "Derive the rate which closes the balance (e.g. a supply) from "
                "the baseline (`pkpdsim.pbpk.steady_state`).",
                SteadyStateWarning,
                stacklevel=3,
            )
