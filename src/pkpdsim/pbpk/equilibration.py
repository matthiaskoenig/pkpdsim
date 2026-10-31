"""The drug-free equilibration of a condition of a model of the library.

`equilibrated_simulation` equilibrates a simulation of a model file, the
equilibration of `PBPKSimulationExperiment.equilibrated`
(`pkpdsim.pbpk.experiment`, on the file of its loaded model) and of the groups
of a sensitivity analysis: the presimulation of sbmlsim to the steady
state of the condition, the doses moved to the start, the sinks, the
readouts of the intestine which nothing drains and the resets set back, the
species of the drugs set to exactly 0
(`drug_free_zeros`) and the states which grow without a steady state refused
(`accumulating_states`, `growing_states`). A model file is parsed once per
process (`pkpdsim.pbpk.flat_model`, by its path and modification time, shared
with the validity of the kidney outflow); the docs are
the section "Equilibration of a condition" of the impairment page.

The module imports sbmlsim; `pkpdsim.pbpk` exports its names lazily.
"""

from __future__ import annotations

import math
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import overload

import libsbml
from sbmlsim.simulation import Change, Scan, Simulation, SteadyState
from sbmlsim.simulation.definition import Value
from sbmlsim.simulation.scan import DimensionKind
from sbmlsim.simulator import Simulator
from sbmlsim.uncertainty import Uncertain
from sbmlsim.units import Quantity

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.flat_model import accumulating, flat_model, loaded_model
from pkpdsim.pbpk.organs import Organ
from pkpdsim.scan import start_values

__all__ = [
    "DOSE_PREFIXES",
    "EQUILIBRATION_TOLERANCES",
    "Equilibration",
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
    return accumulating(libsbml.readSBMLFromString(sbml).getModel())


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
        states = flat_model(path).accumulating
        if not states:
            _GROWING[key] = ()
            return ()
        model = loaded_model(path)
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
    graph = flat_model(path)
    fed = set(graph.started) | ({_bare(sid) for sid in set_ids} & graph.candidates)
    changed = True
    while changed:
        changed = False
        for reactants, products, reversible, reads in graph.reactions:
            dose = not reactants and any(n.startswith(DOSE_PREFIXES) for n in reads)
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

    The equilibration of `PBPKSimulationExperiment.equilibrated`, which calls
    it with the file of its loaded model, and of a simulation outside an
    experiment, e.g. the groups of a sensitivity analysis
    (`pkpdsim.sensitivity.SensitivityGroup(prepare=Equilibration(...))`): the ids of
    the doses, the sinks, the readouts of the intestine and the deltas are
    read from the model file (the doses move to the start, the sinks, the
    readouts which nothing drains and `resets` are set back, the species of
    the drugs are exactly 0, a growing state without a reset is refused,
    every point of a scan equilibrates in its own condition), see
    `PBPKSimulationExperiment.equilibrated`. It replaces a long pre-run from
    a negative start, which leaves a slow state off its steady state.

    Args:
        simulation: the simulation or the scan, in the condition.
        model_path: the flat SBML model the simulation runs.
        resets: the values at the start of the further states which grow
            without a steady state, e.g. the urine of a PD component; a
            species by its id (a concentration species is set as `[sid]`).
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
    parsed = flat_model(path)
    steady_state = SteadyState(
        absolute_tolerance=absolute_tolerance,
        relative_tolerance=relative_tolerance,
        max_time=max_time,
    )
    sinks = sorted(
        s for s in parsed.species if s.startswith((ids.urine(""), ids.feces("")))
    )
    # the readouts of the intestine which nothing drains, the amount absorbed
    # and the oral dose applied (`Fa_<s>`, in the intestine or at the top)
    readouts = sorted(_intestine_readouts(parsed.parameters))
    # a species whose id is its concentration is set as `[sid]`: the same 0,
    # and the quantity of the species for PEtab (a condition of PEtab sets the
    # concentration of such a species, the amount of an amount species)
    concentrations = parsed.concentrations
    # the bases of the deltas, evaluated after the equilibration
    deltas = sorted(
        sid
        for sid in parsed.parameters
        if ids.delta_base(sid) in parsed.parameters
        and ids.delta(sid) in parsed.parameters
    )
    base = simulation.simulation if isinstance(simulation, Scan) else simulation
    # what a dimension of values or the condition sets is no zero
    scanned = {
        target
        for dimension in (simulation.dimensions if isinstance(simulation, Scan) else ())
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
    # a reset names a species by its id or as `[sid]`, set as the zeros are
    reset = {
        _quantities([_bare(target)], concentrations)[0]: value
        for target, value in (resets or {}).items()
    }
    reset_ids = {_bare(target) for target in reset}
    at_start: dict[str, float | Quantity | str] = {
        **dict.fromkeys(
            _quantities(
                (s for s in (*sinks, *readouts, *zeros) if s not in reset_ids),
                concentrations,
            ),
            0.0,
        ),
        **reset,
        **{ids.delta_base(sid): sid for sid in deltas},
    }
    equilibrated = _equilibrated(base, DOSE_PREFIXES, at_start, steady_state)
    growing = sorted(
        set(growing_states(path, equilibrated.preinit_changes, steady_state))
        - {_bare(target) for target in at_start}
    )
    if growing:
        raise ValueError(
            f"The states {growing} grow without a steady state during the "
            f"equilibration of the model '{path.name}' (they have no outflow): "
            "give their values at the start in `resets`."
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


def _intestine_readouts(parameters: Collection[str]) -> set[str]:
    """The readouts of the intestine of a model which only accumulate.

    The amount absorbed from the lumen `Aabs_<s>` and the oral dose applied
    `Aapp_<s>` of every substance (`Mr_<s>`), of the intestine
    (`GU__Aabs_<s>`, the readout `Fa_<s>` alone) or of the top model: states
    which nothing drains, set back at the start as the sinks.
    """
    substances = [
        sid.removeprefix(ids.mr("")) for sid in parameters if sid.startswith(ids.mr(""))
    ]
    locals_ = (
        f(s) for s in substances for f in (ids.absorbed_amount, ids.oral_dose_applied)
    )
    return {
        sid
        for local in locals_
        for sid in (local, ids.flat(Organ.SMALL_INTESTINE, local))
        if sid in parameters
    }


@dataclass(frozen=True)
class Equilibration:
    """The drug-free equilibration of the scans of a group of a sensitivity analysis.

    The function of a `pkpdsim.sensitivity.SensitivityGroup` which prepares
    its scans (`prepare`): every point of a design starts from the drug-free
    steady state of its condition on the build of the group
    (`equilibrated_simulation`), which replaces a long pre-run from a
    negative start:

    ```python
    # not executed
    SensitivityGroup(
        "renal_mild", "Mild renal impairment", simulation,
        prepare=Equilibration(resets={"FLUID__Aurine_h2o": 0.0}),
    )
    ```

    Attributes:
        resets: the values at the start of the further states which grow
            without a steady state, e.g. the urine of a PD component.
        absolute_tolerance: the absolute tolerance of the rates (per minute).
        relative_tolerance: the relative tolerance of the rates.
        max_time: the longest equilibration (minutes).
    """

    resets: Mapping[str, Quantity | float] = field(
        default_factory=lambda: MappingProxyType({})
    )
    absolute_tolerance: float = EQUILIBRATION_TOLERANCES[0]
    relative_tolerance: float = EQUILIBRATION_TOLERANCES[1]
    max_time: float = EQUILIBRATION_TOLERANCES[2]

    def __post_init__(self) -> None:
        """Keep the resets read-only."""
        object.__setattr__(self, "resets", MappingProxyType(dict(self.resets)))

    def __call__(self, scan: Scan, model_path: Path) -> Scan:
        """Equilibrate the scan of a group on its build.

        Args:
            scan: the scan of the group over a design.
            model_path: the build of the group.

        Returns:
            The equilibrated scan.
        """
        return equilibrated_simulation(
            scan,
            model_path,
            resets=self.resets,
            absolute_tolerance=self.absolute_tolerance,
            relative_tolerance=self.relative_tolerance,
            max_time=self.max_time,
        )


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
