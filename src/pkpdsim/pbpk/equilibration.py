"""The drug-free equilibration of a condition of a model of the library.

`PBPKSimulationExperiment.equilibrated` (`pkpdsim.pbpk.experiment`) and
`equilibrated_simulation` outside an experiment (the groups of a sensitivity
analysis) share `_equilibrate`: the presimulation of sbmlsim to the steady
state of the condition, the doses moved to the start, the sinks and the
resets set back, the species of the drugs set to exactly 0
(`drug_free_zeros`) and the states which grow without a steady state refused
(`accumulating_states`, `growing_states`). A model file is parsed once per
process (`_model_file`, by its path and modification time); the docs are
the section "Equilibration of a condition" of the impairment page.

The module imports sbmlsim; `pkpdsim.pbpk` exports its names lazily.
"""

from __future__ import annotations

import math
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import overload

import libsbml
from sbmlsim.model.model_roadrunner import RoadrunnerSBMLModel
from sbmlsim.simulation import Change, Scan, Simulation, SteadyState
from sbmlsim.simulation.definition import Value
from sbmlsim.simulation.scan import DimensionKind
from sbmlsim.simulator import Simulator
from sbmlsim.uncertainty import Uncertain
from sbmlsim.units import Quantity

from pkpdsim.pbpk import ids
from pkpdsim.scan import start_values

__all__ = [
    "DOSE_PREFIXES",
    "EQUILIBRATION_TOLERANCES",
    "accumulating_states",
    "drug_free_zeros",
    "equilibrated_simulation",
    "growing_states",
]

#: the absolute and the relative tolerance of the rates of change (per minute)
#: and the longest time (minutes) of the equilibration of a condition, see
#: `PBPKSimulationExperiment.equilibrated`
EQUILIBRATION_TOLERANCES: tuple[float, float, float] = (1e-12, 1e-10, 1e14)


def accumulating_states(sbml: str) -> tuple[str, ...]:
    """Get the states of a model which nothing drains: they only accumulate.

    A state (the variable of a rate rule, a species which reactions change)
    which no rate reads, neither a rate rule nor a kinetic law, directly or
    through the assignment rules, and which no reaction consumes, e.g. the
    urine and the feces of the library or the urine of a PD component: with an
    inflow at the steady state of the rest of the model it grows without a
    steady state. A species which a rate reads is read with its compartment.

    Args:
        sbml: the SBML of the model, flat.

    Returns:
        The ids, sorted.
    """
    return _accumulating(libsbml.readSBMLFromString(sbml).getModel())


def _accumulating(model: libsbml.Model) -> tuple[str, ...]:
    """The accumulating states of a parsed model, see `accumulating_states`."""
    names = _names
    assigned: dict[str, set[str]] = {}
    rated: dict[str, set[str]] = {}
    for rule in model.getListOfRules():
        if rule.isAssignment():
            assigned[rule.getVariable()] = names(rule.getMath())
        elif rule.isRate():
            rated[rule.getVariable()] = names(rule.getMath())
    read = set().union(*rated.values())
    consumed: set[str] = set()
    changed: set[str] = set()
    for reaction in model.getListOfReactions():
        law = reaction.getKineticLaw()
        if law is not None:
            local = {p.getId() for p in law.getListOfLocalParameters()}
            read |= names(law.getMath()) - local
        consumed |= {r.getSpecies() for r in reaction.getListOfReactants()}
        changed |= {r.getSpecies() for r in reaction.getListOfProducts()}
    # what the rates read through the assignment rules, and the compartments of
    # the species they read
    compartments = {s.getId(): s.getCompartment() for s in model.getListOfSpecies()}
    stack = list(read)
    while stack:
        sid = stack.pop()
        for further in (assigned.get(sid, set()) | {compartments.get(sid, "")}) - read:
            if further:
                read.add(further)
                stack.append(further)
    states = set(rated) | (changed - consumed)
    return tuple(sorted(states - read - consumed))


@lru_cache(maxsize=4)
def _accumulation_model(
    model_path: Path, mtime_ns: int
) -> tuple[RoadrunnerSBMLModel, tuple[str, ...]]:
    """Load a model once per process with its accumulating states as selections."""
    states = _model_file(model_path, mtime_ns).accumulating
    model = RoadrunnerSBMLModel(
        source=model_path, selections=["time", *states], jit_optimization="O0"
    )
    return model, states


#: the growing states of a model file by its path, modification time and the
#: tolerances of the equilibration, see `growing_states`
_GROWING: dict[tuple[Path, int, float, float, float], tuple[str, ...]] = {}


def growing_states(
    model_path: Path,
    condition: Mapping[str, Value],
    steady_state: SteadyState,
) -> tuple[str, ...]:
    """Get the states which grow without a steady state during an equilibration.

    The accumulating states of the model (`accumulating_states`) whose value
    after the equilibration of the condition differs from their value at the
    initialization by more than the absolute tolerance (and 1e-9 of the
    value): an inflow at the steady state, which the presimulation reaches
    only relative to their growing amount (after about `1/relative_tolerance`
    minutes), so their value at the start of the simulation is meaningless
    unless they are set back. They are found once per model file (and its
    modification time) and tolerances in a process, in the condition of the
    first equilibration: whether a state has no outflow is a property of the
    model.

    Args:
        model_path: the flat SBML model.
        condition: the changes before the initialization.
        steady_state: the presimulation.

    Returns:
        The ids of the growing states.

    Raises:
        SteadyStateError: if the model reaches no steady state in the
            condition by `max_time` (sbmlsim).
    """
    path = Path(model_path)
    key = (
        path,
        path.stat().st_mtime_ns,
        steady_state.absolute_tolerance,
        steady_state.relative_tolerance,
        steady_state.max_time,
    )
    if key not in _GROWING:
        model, states = _accumulation_model(path, key[1])
        if not states:
            _GROWING[key] = ()
            return ()
        simulator = Simulator(n_workers=1)
        initial = simulator.run(
            model, Simulation(end=1, preinit_changes=dict(condition))
        )
        equilibrated = simulator.run(
            model,
            Simulation(
                end=1, preinit_changes=dict(condition), presimulation=steady_state
            ),
        )
        growing = []
        for sid in states:
            before = float(start_values(initial, sid).item())
            after = float(start_values(equilibrated, sid).item())
            if abs(after - before) > steady_state.absolute_tolerance + 1e-9 * abs(
                before
            ):
                growing.append(sid)
        _GROWING[key] = tuple(growing)
    return _GROWING[key]


#: the prefixes of the ids of the doses of the library, which the
#: equilibration moves to the start
DOSE_PREFIXES: tuple[str, ...] = tuple(
    f("") for f in (ids.po_dose, ids.iv_dose, ids.sol_dose, ids.infusion_rate)
)


@dataclass(frozen=True)
class _ModelFile:
    """What the equilibration reads from a flat model file, parsed once per file.

    Attributes:
        species: the ids of the species.
        parameters: the ids of the parameters.
        accumulating: the states which nothing drains (`accumulating_states`).
        candidates: the species of the substances without a baseline which
            start at 0 (no initial assignment, no rule, not constant), see
            `drug_free_zeros`.
        started: the species with a start: an initial value other than 0,
            an initial assignment or a rule.
        concentrations: the species whose id means their concentration
            (`hasOnlySubstanceUnits` false): a change sets them as `[sid]`.
        reactions: (reactants, products, reversible, a source of a dose) of
            every reaction; a reaction without reactants whose rate reads a
            dose is a dose of the library, 0 without it.
    """

    species: frozenset[str]
    parameters: frozenset[str]
    accumulating: tuple[str, ...]
    candidates: frozenset[str]
    started: frozenset[str]
    concentrations: frozenset[str]
    reactions: tuple[tuple[frozenset[str], frozenset[str], bool, bool], ...]


@lru_cache(maxsize=8)
def _model_file(model_path: Path, mtime_ns: int) -> _ModelFile:
    """Parse a flat model file once, `mtime_ns` a key of the cache."""
    model = libsbml.readSBMLFromFile(str(model_path)).getModel()
    parameters = {p.getId() for p in model.getListOfParameters()}
    substances = {
        pid.removeprefix("Mr_") for pid in parameters if pid.startswith("Mr_")
    }
    exogenous = {s for s in substances if ids.baseline(s) not in parameters}
    assigned = {a.getSymbol() for a in model.getListOfInitialAssignments()}
    ruled = {r.getVariable() for r in model.getListOfRules()}
    species = list(model.getListOfSpecies())
    all_ids = [sp.getId() for sp in species]
    of_exogenous = {sid for s in exogenous for sid in ids.species_of(s, all_ids)}
    candidates, started = set(), set()
    for sp in species:
        sid = sp.getId()
        value = (
            sp.getInitialConcentration()
            if sp.isSetInitialConcentration()
            else sp.getInitialAmount()
        )
        if sid in assigned or sid in ruled or value != 0.0:
            started.add(sid)
        elif (
            sid in of_exogenous
            and not sp.getConstant()
            and not (sp.getBoundaryCondition())
        ):
            candidates.add(sid)
    reactions = []
    for reaction in model.getListOfReactions():
        reactants = frozenset(r.getSpecies() for r in reaction.getListOfReactants())
        products = frozenset(r.getSpecies() for r in reaction.getListOfProducts())
        law = reaction.getKineticLaw()
        reads = set() if law is None else _names(law.getMath())
        dose = not reactants and any(n.startswith(DOSE_PREFIXES) for n in reads)
        reactions.append((reactants, products, reaction.getReversible(), dose))
    return _ModelFile(
        species=frozenset(all_ids),
        parameters=frozenset(parameters),
        accumulating=_accumulating(model),
        candidates=frozenset(candidates),
        started=frozenset(started),
        concentrations=frozenset(
            sp.getId() for sp in species if not sp.getHasOnlySubstanceUnits()
        ),
        reactions=tuple(reactions),
    )


def _names(ast: libsbml.ASTNode | None) -> set[str]:
    """The names a formula reads."""
    found: set[str] = set()
    stack = [] if ast is None else [ast]
    while stack:
        node = stack.pop()
        if node.getType() == libsbml.AST_NAME:
            found.add(node.getName())
        stack.extend(node.getChild(k) for k in range(node.getNumChildren()))
    return found


def drug_free_zeros(model_path: Path, set_ids: Iterable[str] = ()) -> tuple[str, ...]:
    """Get the species which are exactly 0 at the drug-free steady state of a model.

    The species of the substances without a baseline (`Mr_<s>` without
    `baseline_<s>`: the drugs and their metabolites, not an endogenous
    substance), which start at 0 and which no reaction feeds from a species
    with a start: a reaction moves a substance only from a species which
    holds some (a reactant, or a product of a reversible reaction), a
    reaction without reactants (a source) feeds its products unless it is a
    dose of the library (its rate reads `PODOSE_<s>`, `IVDOSE_<s>`, ...,
    which the equilibration moves to the start). The direction of a reaction
    is its SBML `reversible` flag: a `CustomReaction` whose rate law can run
    backwards (an exchange `k*(A - B)`) must be written reversible
    (`A <-> B`), else a species it feeds backwards from a started product is
    taken for 0; a source without reactants whose rate reads a dose id is
    taken for a dose. In exact arithmetic they
    stay 0 during a drug-free equilibration; a presimulation leaves a
    residue of the rounding in them (1.5e-29 mM in hctz), which
    `PBPKSimulationExperiment.equilibrated` sets back to 0.

    Args:
        model_path: the flat SBML model.
        set_ids: the ids which a simulation sets before the equilibration
            (its condition, the dimensions of a scan): a species among them
            holds what it is set to and feeds what it reaches; a
            concentration `[sid]` is the species `sid`.

    Returns:
        The ids, sorted.
    """
    path = Path(model_path)
    graph = _model_file(path, path.stat().st_mtime_ns)
    fed = set(graph.started) | ({_bare(sid) for sid in set_ids} & graph.candidates)
    changed = True
    while changed:
        changed = False
        for reactants, products, reversible, dose in graph.reactions:
            sides = reactants | products if reversible else reactants
            active = bool(sides & fed) or (not reactants and not dose)
            reached = (reactants | products if reversible else products) - fed
            if active and reached:
                fed |= reached
                changed = True
    return tuple(sorted(graph.candidates - fed))


@overload
def equilibrated_simulation(
    simulation: Simulation,
    model_path: Path,
    *,
    resets: Mapping[str, Quantity | float] | None = None,
    absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0],
    relative_tolerance: float = EQUILIBRATION_TOLERANCES[1],
    max_time: float = EQUILIBRATION_TOLERANCES[2],
) -> Simulation: ...


@overload
def equilibrated_simulation(
    simulation: Scan,
    model_path: Path,
    *,
    resets: Mapping[str, Quantity | float] | None = None,
    absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0],
    relative_tolerance: float = EQUILIBRATION_TOLERANCES[1],
    max_time: float = EQUILIBRATION_TOLERANCES[2],
) -> Scan: ...


@overload
def equilibrated_simulation(
    simulation: Simulation | Scan,
    model_path: Path,
    *,
    resets: Mapping[str, Quantity | float] | None = None,
    absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0],
    relative_tolerance: float = EQUILIBRATION_TOLERANCES[1],
    max_time: float = EQUILIBRATION_TOLERANCES[2],
) -> Simulation | Scan: ...


def equilibrated_simulation(
    simulation: Simulation | Scan,
    model_path: Path,
    *,
    resets: Mapping[str, Quantity | float] | None = None,
    absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0],
    relative_tolerance: float = EQUILIBRATION_TOLERANCES[1],
    max_time: float = EQUILIBRATION_TOLERANCES[2],
) -> Simulation | Scan:
    """Get a simulation of a model file which starts from the drug-free steady state of its condition.

    The equilibration of `PBPKSimulationExperiment.equilibrated` outside an
    experiment, e.g. for the groups of a sensitivity analysis
    (`pkpdsim.sensitivity.SensitivityGroup(equilibrate=True)`): the ids of
    the doses, the sinks and the deltas are read from the model file instead
    of a loaded model, everything else is the same (the doses move to the
    start, the sinks and `resets` are set back, the species of the drugs are
    exactly 0, a growing state without a reset is refused, every point of a
    scan equilibrates in its own condition). It replaces a long pre-run from
    a negative start, which leaves a slow state off its steady state.

    Args:
        simulation: the simulation or the scan, in the condition.
        model_path: the flat SBML model the simulation runs.
        resets: the values at the start of the further states which grow
            without a steady state, e.g. the urine of a PD component.
        absolute_tolerance: the absolute tolerance of the rates of the
            equilibration (per minute).
        relative_tolerance: the relative tolerance of the rates.
        max_time: the longest equilibration (minutes).

    Returns:
        The equilibrated simulation or scan, a new object.

    Raises:
        ValueError: as `PBPKSimulationExperiment.equilibrated`.
    """
    path = Path(model_path)
    parsed = _model_file(path, path.stat().st_mtime_ns)
    return _equilibrate(
        simulation,
        species=parsed.species,
        parameters=parsed.parameters,
        path=path,
        resets=resets,
        steady_state=SteadyState(
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
            max_time=max_time,
        ),
        where=f"of the model '{path.name}'",
    )


def _equilibrate(
    simulation: Simulation | Scan,
    *,
    species: Collection[str],
    parameters: Collection[str],
    path: Path | None,
    resets: Mapping[str, Quantity | float] | None,
    steady_state: SteadyState,
    where: str,
) -> Simulation | Scan:
    """Equilibrate a simulation or a scan of a model, see `PBPKSimulationExperiment.equilibrated`.

    Args:
        simulation: the simulation or the scan.
        species: the species of the model.
        parameters: the parameters of the model.
        path: the model file, `None` for SBML text: then no growing state is
            refused and no species is set to 0.
        resets: the further values at the start.
        steady_state: the presimulation.
        where: the experiment or the model, for the messages.
    """
    sinks = sorted(s for s in species if s.startswith((ids.urine(""), ids.feces(""))))
    # a species whose id is its concentration is set as `[sid]`: the same 0,
    # and the quantity of the species for PEtab (a condition of PEtab sets the
    # concentration of such a species, the amount of an amount species)
    concentrations: frozenset[str] = frozenset()
    # the bases of the deltas, evaluated after the equilibration
    deltas = sorted(
        sid
        for sid in parameters
        if ids.delta_base(sid) in parameters and ids.delta(sid) in parameters
    )
    base = simulation.simulation if isinstance(simulation, Scan) else simulation
    zeros: tuple[str, ...] = ()
    if path is not None:
        concentrations = _model_file(path, path.stat().st_mtime_ns).concentrations
        # what a dimension of values or the condition sets is no zero
        scanned = {
            target
            for dimension in (
                simulation.dimensions if isinstance(simulation, Scan) else ()
            )
            if dimension.kind == DimensionKind.VALUES
            for target in dimension.values
        }
        preinit = {**base.preinit_changes, **base.uncertain_preinit}
        set_at_start = {_bare(target) for target in _set_at_start(base)}
        zeros = tuple(
            sid
            for sid in drug_free_zeros(path, set(preinit) | scanned)
            if sid not in set_at_start
        )
    at_start: dict[str, float | Quantity | str] = {
        **dict.fromkeys(_quantities(sinks, concentrations), 0.0),
        **dict.fromkeys(_quantities(zeros, concentrations), 0.0),
        **(resets or {}),
        **{ids.delta_base(sid): sid for sid in deltas},
    }
    equilibrated = _equilibrated(base, DOSE_PREFIXES, at_start, steady_state)
    if path is not None:
        growing = sorted(
            set(growing_states(path, equilibrated.preinit_changes, steady_state))
            - {_bare(target) for target in at_start}
        )
        if growing:
            raise ValueError(
                f"The states {growing} grow without a steady state during the "
                f"equilibration {where} (they have no outflow): give their values "
                "at the start in `resets`."
            )
    if not isinstance(simulation, Scan):
        return equilibrated
    given = {target for change in equilibrated.changes for target in change.values}
    for dimension in simulation.dimensions:
        if dimension.kind != DimensionKind.VALUES:
            raise ValueError(
                f"The dimension '{dimension.id}' varies {dimension.kind}: "
                "equilibrate every simulation of it instead."
            )
        missing = sorted(
            target
            for target in dimension.values
            if target.startswith(DOSE_PREFIXES)
            and dimension.at is None
            and target not in given
        )
        if missing:
            raise ValueError(
                f"The dimension '{dimension.id}' sets the doses {missing}, "
                "which the simulation does not give: they would be doses "
                "before the equilibration; give them in the simulation."
            )
    return Scan(simulation=equilibrated, dimensions=simulation.dimensions)


def _bare(target: str) -> str:
    """The id of a target of a change, `sid` of a concentration `[sid]`."""
    return target[1:-1] if target.startswith("[") and target.endswith("]") else target


def _quantities(sids: Iterable[str], concentrations: Collection[str]) -> list[str]:
    """The targets of changes which set the species, `[sid]` of a concentration."""
    return [f"[{sid}]" if sid in concentrations else sid for sid in sids]


def _set_at_start(simulation: Simulation) -> set[str]:
    """The targets which a change of a simulation sets at its start."""
    start = _time_value(simulation, simulation.start)
    return {
        target
        for change in simulation.changes
        if any(
            math.isclose(_time_value(simulation, t), start, abs_tol=1e-12)
            for t in change.times
        )
        for target in change.values
    }


def _equilibrated(
    simulation: Simulation,
    dose_prefixes: tuple[str, ...],
    at_start: Mapping[str, float | Quantity | str],
    steady_state: SteadyState,
) -> Simulation:
    """Equilibrate a simulation, see `PBPKSimulationExperiment.equilibrated`.

    Args:
        simulation: the simulation.
        dose_prefixes: the prefixes of the ids of the doses.
        at_start: the further values at the start, after the doses.
        steady_state: the presimulation.

    Raises:
        ValueError: for a simulation with a presimulation or a dose before the
            initialization and in a change at the start.
    """
    if simulation.presimulation is not None:
        raise ValueError(
            "The simulation has a presimulation already; it is equilibrated once."
        )
    preinit = {**simulation.preinit_changes, **simulation.uncertain_preinit}
    doses = {k: v for k, v in preinit.items() if k.startswith(dose_prefixes)}
    condition = {k: v for k, v in preinit.items() if k not in doses}
    start = simulation.start
    twice = sorted(set(doses) & _set_at_start(simulation))
    if twice:
        raise ValueError(
            f"The doses {twice} are set before the initialization and by a change at "
            "its start; give each dose once, the equilibration moves the doses "
            "before the initialization to the start."
        )
    values: dict[str, float | Quantity | str | Uncertain] = {**doses, **at_start}
    return Simulation(
        time_unit=simulation.time_unit,
        start=start,
        end=simulation.end,
        preinit_changes=condition,
        changes=[Change(start, values), *simulation.changes],
        presimulation=steady_state,
        times=None if simulation.times is None else list(simulation.times),
        steps=simulation.steps,
        time_shift=simulation.time_shift,
    )


def _time_value(simulation: Simulation, time: object) -> float:
    """Get a time of a simulation as a number in its time unit."""
    if isinstance(time, Quantity):
        unit = simulation.time_unit
        return float(time.m_as(unit) if unit else time.magnitude)
    return float(time)  # ty: ignore[invalid-argument-type]
