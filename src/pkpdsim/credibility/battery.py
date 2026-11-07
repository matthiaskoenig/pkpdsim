"""The battery of the checks of a PBPK model: the PK checks of its builds.

The battery of sbmlsim (`sbmlsim.credibility.battery`) checks every variant of
a project, a build of a model per sex (`pkpdsim.model.PKPDModel.variants`):
non-negativity, steady state, a second simulator, plausibility and the
rebuild through the factory. The PK battery `PK_BATTERY` (`PKBattery`, the
battery hook of every `PKPDModel`) adds what every PBPK model must fulfil, on
the simulations of every build with roadrunner without a dose and with a dose
of every route of every substance:

- `mass_balance`: the amount of the substances in the model, in its sinks and
  in its dose states changes by the infusion and by the flows of the reactions
  which do not read a dose (eliminations and productions) alone, so a dose
  enters with its molecular weight and nothing appears or disappears outside
  of the reactions (`MASS_LIMIT`); the fate of every dose is recorded (`Fate`,
  the section `fates` of the battery, `fates.tsv` and the table Fate of the
  doses of the report).
- `zero_dose`: without a dose the species of a substance which starts at
  zero (a drug, not an endogenous substance) stay zero.
- `superposition`: twice a small dose gives twice the response, and two small
  doses, at 0 and at `SECOND_DOSE`, the sum of the responses to each
  (`LIMIT`): the model is linear and time invariant in the limit of a small
  dose.

The hook gives the settings of the battery (`HORIZON` two days on a grid of
`STEPS`, the tolerances, a model starting at its steady state, AMICI as the
second simulator), the SBML with a parameter per flow of a reaction which
integrates it, the cases of the doses (`DOSE` of every dose state and infusion
rate, an infusion over `INFUSION`, repeated by the second simulator, and the
small doses of the superposition) and the order of the checks in the report.

The battery reads the conventions of the models of pkdb_models from the flat
SBML (`FlatBuild`): a substance has its molecular weight `Mr_<s>`, its doses
are the states `IVDOSE_<s>`, `PODOSE_<s>` and `SOLDOSE_<s>` (mg) and the
infusion rate `Ri_<s>` (mg/min), and a species belongs to the substance its
id (the part after the last `__`) ends with, `_<s>`, or else starts with,
`<s>_` (the lumen of the intestine of some models, `mid_lumen`). A check
without what it needs, e.g. a model without doses, is not applicable and
assesses nothing.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final, Self, override

import numpy as np
import pandas as pd
from sbmlsim.credibility.battery import (
    BatteryCase,
    BatteryCollector,
    BatteryContribution,
    BatteryHook,
    BatteryResult,
    BatteryRun,
    BatterySettings,
    not_applicable,
)
from sbmlsim.credibility.verification import (
    LIMIT,
    SCALE_FLOOR,
    Check,
    register_check,
)

from pkpdsim.credibility.blocks import render_macro

if TYPE_CHECKING:
    import libsbml

__all__ = [
    "ABSOLUTE_TOLERANCE",
    "DOSE",
    "FATE_COLUMNS",
    "HORIZON",
    "INFUSION",
    "MASS_BALANCE",
    "MASS_LIMIT",
    "PK_BATTERY",
    "PK_CHECKS",
    "RELATIVE_TOLERANCE",
    "SECOND_DOSE",
    "SMALL_DOSE",
    "STEADY_LIMIT",
    "STEPS",
    "SUPERPOSITION",
    "ZERO_DOSE",
    "Fate",
    "FlatBuild",
    "PKBattery",
]

logger = logging.getLogger(__name__)

#: the dose of a route (mg), an infusion gives it over `INFUSION`
DOSE: Final = 10.0
#: the small dose of the superposition (mg)
SMALL_DOSE: Final = 0.01
#: the duration of an infusion (min)
INFUSION: Final = 60.0
#: the end of every simulation (min), two days
HORIZON: Final = 2880.0
#: the time of the second dose of the superposition (min)
SECOND_DOSE: Final = 720.0
#: the intervals of the grid of the results, 10 min
STEPS: Final = 288
#: the relative tolerance of the simulations of the battery
RELATIVE_TOLERANCE: Final = 1e-10
#: the absolute tolerance of the simulations of the battery, also the floor of
#: a value which is zero
ABSOLUTE_TOLERANCE: Final = 1e-16
#: the largest imbalance of the amounts relative to the largest term of the balance
MASS_LIMIT: Final = 1e-6
#: the largest change of a state without a dose relative to its start
STEADY_LIMIT: Final = 1e-6

#: the prefixes of the dose states (mg) and of the infusion rate (mg/min)
_BOLUS: Final = ("IVDOSE", "PODOSE", "SOLDOSE")
_INFUSION: Final = "Ri"
#: the prefixes of the sinks of a substance, by the name of their part of a fate
_SINKS: Final = {"urine": "Aurine_", "feces": "Afeces_"}
#: the prefix of the parameters which integrate the flows of the reactions
_FLOW: Final = "battery_flow_"

MASS_BALANCE: Final = register_check(
    Check("mass_balance", "mass balance", "numerical_code_verification")
)
ZERO_DOSE: Final = register_check(
    Check("zero_dose", "zero dose", "numerical_code_verification")
)
SUPERPOSITION: Final = register_check(
    Check("superposition", "superposition", "numerical_code_verification")
)
#: the checks of the PK battery
PK_CHECKS: Final = (MASS_BALANCE, ZERO_DOSE, SUPERPOSITION)

#: the columns of `fates.tsv`
FATE_COLUMNS: Final = (
    "build",
    "dose",
    "substance",
    "amount",
    "part",
    "value",
    "share",
)


@dataclass(frozen=True, kw_only=True)
class Fate:
    """Where a dose went at the end of the horizon, the mass balance of a route.

    The parts are the differences to the simulation without a dose, so the
    fate of a dose of an endogenous substance is that of the dose and not of
    its production.

    Attributes:
        build: the build, e.g. `male`.
        dose: the dose state or infusion rate, e.g. `IVDOSE_tor`.
        substance: the dosed substance.
        amount: the dose (mmole).
        parts: the parts of the dose (mmole) by name: `body` (the species
            which are no sink), `urine` and `feces` (the sinks), `dose` (not
            yet applied) and a reaction which does not read a dose by its id
            (negative for an elimination, positive for a production).
        imbalance: the imbalance of the dose simulation relative to the
            largest term of its balance.
    """

    build: str
    dose: str
    substance: str
    amount: float
    parts: Mapping[str, float]
    imbalance: float

    def to_dict(self) -> dict[str, object]:
        """Get the fate as plain data."""
        return {
            "build": self.build,
            "dose": self.dose,
            "substance": self.substance,
            "amount": self.amount,
            "parts": dict(self.parts),
            "imbalance": self.imbalance,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the fate from plain data."""
        return cls(
            build=data["build"],
            dose=data["dose"],
            substance=data["substance"],
            amount=float(data["amount"]),
            parts={k: float(v) for k, v in data["parts"].items()},
            imbalance=float(data["imbalance"]),
        )


def _substance_of(sid: str, substances: Sequence[str]) -> str | None:
    """Get the substance of an id, by its suffix `_<s>` or its prefix `<s>_`.

    The local id (after the last `__`) is matched, the suffix before the prefix
    and the longest substance first, e.g. `Cve_mid1oh` and `mid_lumen`.
    """
    local = sid.rsplit("__", 1)[-1]
    for substance in substances:
        if local.endswith(f"_{substance}"):
            return substance
    for substance in substances:
        if local.startswith(f"{substance}_"):
            return substance
    return None


def _reads(ast: libsbml.ASTNode | None, assigned: Mapping[str, set[str]]) -> set[str]:
    """Get the names a formula reads, also through the assignment rules."""
    from pkpdsim.pbpk.flat_model import math_names

    found = math_names(ast)
    stack = list(found)
    while stack:
        for name in assigned.get(stack.pop(), set()) - found:
            found.add(name)
            stack.append(name)
    return found


@dataclass(frozen=True, kw_only=True)
class FlatBuild:
    """What the battery reads from a flat build of a model.

    Attributes:
        name: the build, e.g. `male`.
        path: the flat SBML file.
        substances: the substances with a molecular weight, the longest id
            first.
        species: the species which are states (not constant, no boundary
            condition, no assignment rule).
        tracked: the substance of every species which belongs to one.
        amounts: the species whose id is their amount (`hasOnlySubstanceUnits`),
            AMICI gives them as amounts and the others as concentrations.
        doses: the substance of every dose state and infusion rate, which a
            dose is given to.
        dose_states: the substance of every state of the dose (mg): the dose
            states and the parameters in mg with a rate rule which reads one, e.g.
            the drug in the stomach of a model, also through the assignment
            rules and transitively.
        flows: the net change of the amount of the tracked species per unit of
            the rate of every reaction which changes it, by reaction id.
        inputs: the reactions of `flows` which read a state of the dose or
            an infusion rate of a substance they form, also through assignment
            rules.
        outside: the tracked species which a rule or an event changes, not
            a reaction.
        states: the ids of the states without the species: the parameters
            and compartments with a rate rule.
        accumulating: the states which nothing drains (the sinks).
        sbml: the SBML of the build with a parameter per flow which
            integrates it, `battery_flow_<reaction>`.
    """

    name: str
    path: Path
    substances: tuple[str, ...]
    species: tuple[str, ...]
    tracked: Mapping[str, str]
    amounts: frozenset[str]
    doses: Mapping[str, str]
    dose_states: Mapping[str, str]
    flows: Mapping[str, float]
    inputs: frozenset[str]
    outside: tuple[str, ...]
    states: tuple[str, ...]
    accumulating: frozenset[str]
    sbml: str

    @classmethod
    def read(cls, name: str, path: Path) -> FlatBuild:
        """Read a flat build.

        Args:
            name: the build, e.g. `male`.
            path: the flat SBML file.

        Returns:
            The build with the SBML which integrates its flows.
        """
        import libsbml

        from pkpdsim.pbpk.flat_model import accumulating, math_names

        doc = libsbml.readSBMLFromFile(str(path))
        model = doc.getModel()
        parameters = {p.getId(): p for p in model.getListOfParameters()}
        substances = tuple(
            sorted(
                (pid[3:] for pid in parameters if pid.startswith("Mr_")),
                key=lambda s: (-len(s), s),
            )
        )
        # what an assignment rule and the rate of a reaction read: a formula
        # which reads a reaction id reads its rate
        assigned = {
            rule.getVariable(): math_names(rule.getMath())
            for rule in model.getListOfRules()
            if rule.isAssignment()
        }
        for reaction in model.getListOfReactions():
            law = reaction.getKineticLaw()
            if law is not None:
                local = {p.getId() for p in law.getListOfLocalParameters()}
                assigned[reaction.getId()] = math_names(law.getMath()) - local
        rates = {
            rule.getVariable(): _reads(rule.getMath(), assigned)
            for rule in model.getListOfRules()
            if rule.isRate()
        }
        rated = list(rates)
        species = tuple(
            sp.getId()
            for sp in model.getListOfSpecies()
            if not sp.getConstant()
            and not sp.getBoundaryCondition()
            and sp.getId() not in assigned
        )
        tracked = {
            sid: s
            for sid in species
            if (s := _substance_of(sid, substances)) is not None
        }
        doses = {
            pid: s
            for pid in parameters
            for prefix in (*_BOLUS, _INFUSION)
            for s in substances
            if pid == f"{prefix}_{s}"
        }
        dose_states = {
            pid: s for pid, s in doses.items() if pid.split("_")[0] in _BOLUS
        }
        # a parameter in mg whose rate rule reads the dose, e.g. the drug in the
        # stomach, is a state of the dose while a reaction or the rate rule of
        # another state of the dose reads it, not a sum of the given doses
        chain: dict[str, str] = {}
        grown = True
        while grown:
            grown = False
            for pid, reads in rates.items():
                sources = [d for d in (*doses, *chain) if d in reads and d != pid]
                if (
                    pid in parameters
                    and pid not in doses
                    and pid not in chain
                    and parameters[pid].getUnits().rsplit("__", 1)[-1] == "mg"
                    and sources
                ):
                    chain[pid] = {**doses, **chain}[sources[0]]
                    grown = True
        laws = set().union(
            *(
                _reads(r.getKineticLaw().getMath(), assigned)
                for r in model.getListOfReactions()
                if r.getKineticLaw() is not None
            )
        )
        pruned = True
        while pruned:
            pruned = False
            for pid in list(chain):
                others = {*dose_states, *chain} - {pid}
                if pid not in laws and not any(
                    pid in rates[o] for o in others if o in rates
                ):
                    del chain[pid]
                    pruned = True
        dose_states.update(chain)
        dosing = {**doses, **dose_states}
        flows: dict[str, float] = {}
        inputs: set[str] = set()
        for reaction in model.getListOfReactions():
            net = sum(
                p.getStoichiometry()
                for p in reaction.getListOfProducts()
                if p.getSpecies() in tracked
            ) - sum(
                r.getStoichiometry()
                for r in reaction.getListOfReactants()
                if r.getSpecies() in tracked
            )
            if net == 0.0:
                continue
            rid = reaction.getId()
            flows[rid] = net
            law = reaction.getKineticLaw()
            # an input reads a dose of a substance it forms, a secretion which a
            # dose of another substance stimulates is none
            formed = {
                tracked[x.getSpecies()]
                for x in reaction.getListOfProducts()
                if x.getSpecies() in tracked
            }
            if law is not None and formed & {
                dosing[d] for d in _reads(law.getMath(), assigned) & set(dosing)
            }:
                inputs.add(rid)
        outside = {r.getVariable() for r in model.getListOfRules() if r.isRate()}
        outside |= {
            a.getVariable()
            for event in model.getListOfEvents()
            for a in event.getListOfEventAssignments()
        }
        for rid, net in flows.items():
            flow = model.createParameter()
            flow.setId(f"{_FLOW}{rid}")
            flow.setValue(0.0)
            flow.setConstant(False)
            rule = model.createRateRule()
            rule.setVariable(f"{_FLOW}{rid}")
            rule.setMath(libsbml.parseL3Formula(f"{net!r} * {rid}"))
        return cls(
            name=name,
            path=path,
            substances=substances,
            species=species,
            tracked=tracked,
            amounts=frozenset(
                sp.getId()
                for sp in model.getListOfSpecies()
                if sp.getHasOnlySubstanceUnits()
            ),
            doses=doses,
            dose_states=dose_states,
            flows=flows,
            inputs=frozenset(inputs),
            outside=tuple(sorted(outside & set(tracked))),
            states=tuple(v for v in rated if v not in species and v in parameters),
            accumulating=frozenset(accumulating(model)),
            sbml=libsbml.writeSBMLToString(doc),
        )

    @property
    def selections(self) -> list[str]:
        """Get the selections of a simulation: the states and the flows.

        A species is selected as its amount, a parameter by its id.
        """
        return [
            *self.species,
            *self.states,
            *(f"{_FLOW}{rid}" for rid in self.flows),
        ]


def _dose_changes(dose: str, amount: float, at: float) -> dict[float, dict[str, float]]:
    """Get the changes of a dose given at a time, an infusion over `INFUSION`."""
    if dose.split("_")[0] == _INFUSION:
        return {at: {dose: amount / INFUSION}, at + INFUSION: {dose: 0.0}}
    return {at: {dose: amount}}


def _merge(
    *changes: Mapping[float, Mapping[str, float]],
) -> dict[float, dict[str, float]]:
    """Merge the changes of several doses, which change different times."""
    merged: dict[float, dict[str, float]] = {}
    for change in changes:
        for at, values in change.items():
            merged.setdefault(at, {}).update(values)
    return merged


def _case(
    name: str, changes: Mapping[float, Mapping[str, float]], *, repeated: bool
) -> BatteryCase:
    """Get a case of the battery, a dose state increased by its change."""
    added = frozenset(
        sid
        for values in changes.values()
        for sid in values
        if sid.split("_")[0] in _BOLUS
    )
    return BatteryCase(
        name=name, changes=changes, added=added, cross_simulator=repeated
    )


#: the small doses of the superposition, by the name of their case after the dose
_SMALL: Final = (
    "at 0",
    "twice at 0",
    f"at {SECOND_DOSE:g} min",
    f"at 0 and {SECOND_DOSE:g} min",
)


def dose_cases(build: FlatBuild) -> list[BatteryCase]:
    """Get the cases of the doses of a build.

    Every dose state and infusion rate gets `DOSE`, an infusion over
    `INFUSION`, the case named by the dose, e.g. `IVDOSE_d`, which the second
    simulator repeats; and the small doses of the superposition, e.g.
    `IVDOSE_d small at 0`, which it does not.

    Args:
        build: the build.

    Returns:
        The cases, the doses in the order of the build.
    """
    cases: list[BatteryCase] = []
    for dose in build.doses:
        cases.append(_case(dose, _dose_changes(dose, DOSE, 0.0), repeated=True))
        first = _dose_changes(dose, SMALL_DOSE, 0.0)
        second = _dose_changes(dose, SMALL_DOSE, SECOND_DOSE)
        small = (
            first,
            _dose_changes(dose, 2 * SMALL_DOSE, 0.0),
            second,
            _merge(first, second),
        )
        cases += [
            _case(f"{dose} small {at}", changes, repeated=False)
            for at, changes in zip(_SMALL, small, strict=True)
        ]
    return cases


@dataclass(frozen=True)
class _Runner:
    """What the checks read of the simulations of a build: the build, the columns and the values."""

    build: FlatBuild
    index: Mapping[str, int]
    values: Mapping[str, float]


@dataclass(frozen=True)
class _Simulations:
    """The simulations of a build: without a dose, with a dose and the superposition."""

    run: BatteryRun
    runner: _Runner
    zero: np.ndarray
    dosed: dict[str, np.ndarray]
    small: dict[str, dict[str, np.ndarray]]

    @classmethod
    def of(cls, run: BatteryRun) -> Self:
        """Get the simulations of the cases of the doses of a variant."""
        build = run.build.data
        if not isinstance(build, FlatBuild):
            raise TypeError(
                f"The variant '{run.build.variant}' was not read by the PK battery."
            )
        return cls(
            run=run,
            runner=_Runner(build=build, index=run.index, values=run.values),
            zero=run.zero,
            dosed={dose: run.cases[dose] for dose in build.doses},
            small={
                dose: {at: run.cases[f"{dose} small {at}"] for at in _SMALL}
                for dose in build.doses
            },
        )


def _mr(runner: _Runner, substance: str) -> float:
    """Get the molecular weight of a substance (g/mole)."""
    return runner.values[f"Mr_{substance}"]


def _amount(runner: _Runner, result: np.ndarray, row: int) -> dict[str, float]:
    """Get the amounts (mmole) of a row by part: body, sinks and doses."""
    build = runner.build
    parts: dict[str, float] = {"body": 0.0, "urine": 0.0, "feces": 0.0, "dose": 0.0}
    for sid in build.tracked:
        local = sid.rsplit("__", 1)[-1]
        part = next((k for k, p in _SINKS.items() if local.startswith(p)), "body")
        parts[part] += float(result[row, runner.index[sid]])
    for state, substance in build.dose_states.items():
        mg = (
            float(result[row, runner.index[state]])
            if state in runner.index
            else runner.values[state]
        )
        parts["dose"] += mg / _mr(runner, substance)
    return parts


def _balance(
    runner: _Runner, result: np.ndarray, dose: str | None
) -> tuple[dict[str, float], float, float]:
    """Get the parts of a simulation at its end and its balance.

    Returns:
        The parts at the end (mmole): the amounts by part and the integrated
        flows of the reactions which do not read a dose, by reaction; the
        imbalance (mmole) and the largest term of the balance.
    """
    build = runner.build
    begin = _amount(runner, result, 0)
    end = _amount(runner, result, -1)
    # a bolus is in the dose state of the first row, an infusion enters from outside
    infused = 0.0
    if dose is not None and dose.split("_")[0] == _INFUSION:
        infused = DOSE / _mr(runner, build.doses[dose])
    flows = {
        rid: float(result[-1, runner.index[f"{_FLOW}{rid}"]])
        for rid in build.flows
        if rid not in build.inputs
    }
    change = sum(end.values()) - sum(begin.values())
    imbalance = change - infused - sum(flows.values())
    terms = [
        abs(sum(begin.values())),
        abs(sum(end.values())),
        infused,
        *(abs(v) for v in flows.values()),
        *(
            abs(float(result[-1, runner.index[f"{_FLOW}{rid}"]]))
            for rid in build.inputs
        ),
    ]
    return {**end, **flows}, imbalance, max(terms, default=0.0)


def _mass_balance(
    simulations: Sequence[_Simulations], collector: BatteryCollector, fates: list[Fate]
) -> None:
    """Check the mass balance of every dose of every build."""
    for sims in simulations:
        runner = sims.runner
        build = runner.build
        for sid in build.outside:
            collector.problems.append(
                f"{build.name}: the species '{sid}' is changed by a rule or an "
                f"event, outside of the reactions"
            )
        # the amount of the smallest dose is the floor of the largest term
        smallest = min(
            (DOSE / _mr(runner, s) for s in build.doses.values()),
            default=ABSOLUTE_TOLERANCE,
        )
        zero_parts, imbalance, scale = _balance(runner, sims.zero, None)
        collector.cases += 1
        relative = abs(imbalance) / max(scale, smallest)
        collector.deviation(relative, f"{build.name} without a dose")
        if relative > MASS_LIMIT:
            collector.problems.append(
                f"{build.name}: the amounts without a dose change by "
                f"{imbalance:.3g} mmole more than the flows of the reactions "
                f"({relative:.2g} of the largest term)"
            )
        for dose, result in sims.dosed.items():
            substance = build.doses[dose]
            parts, imbalance, scale = _balance(runner, result, dose)
            collector.cases += 1
            relative = abs(imbalance) / max(scale, smallest)
            collector.deviation(relative, f"{build.name} {dose}")
            if relative > MASS_LIMIT:
                collector.problems.append(
                    f"{build.name} {dose}: the amounts change by {imbalance:.3g} "
                    f"mmole more than the dose and the flows of the reactions "
                    f"({relative:.2g} of the largest term)"
                )
            amount = DOSE / _mr(runner, substance)
            fates.append(
                Fate(
                    build=build.name,
                    dose=dose,
                    substance=substance,
                    amount=amount,
                    parts={
                        k: v - zero_parts.get(k, 0.0)
                        for k, v in parts.items()
                        if abs(v - zero_parts.get(k, 0.0)) > SCALE_FLOOR * amount
                    },
                    imbalance=relative,
                )
            )


def _scales(sims: _Simulations, ids: Sequence[str]) -> np.ndarray:
    """Get the largest absolute value of every id over the simulations of a build."""
    return sims.run.scales(ids)


def _floor(scale: np.ndarray) -> np.ndarray:
    """Get the floor of a deviation of an id: `SCALE_FLOOR` of its scale."""
    return np.maximum(SCALE_FLOOR * scale, ABSOLUTE_TOLERANCE)


def _zero_dose(
    simulations: Sequence[_Simulations], collector: BatteryCollector
) -> None:
    """Check that the species of a substance which starts at zero stay zero.

    A species stays zero within the floor of its largest value over the
    simulations of its build.
    """
    for sims in simulations:
        runner = sims.runner
        build = runner.build
        for substance in build.substances:
            species = [sid for sid, s in build.tracked.items() if s == substance]
            values = sims.zero[:, [runner.index[sid] for sid in species]]
            if not species or np.any(values[0] != 0.0):
                continue
            collector.cases += len(species)
            scale = _scales(sims, species)
            high = np.max(np.abs(values), axis=0)
            relative = high / np.maximum(scale, ABSOLUTE_TOLERANCE)
            k = int(np.argmax(relative))
            collector.deviation(float(relative[k]), f"{build.name}: {species[k]}")
            for k in np.flatnonzero(high > _floor(scale)):
                collector.problems.append(
                    f"{build.name}: '{species[k]}' of the substance '{substance}' "
                    f"reaches {high[k]:.3g} without a dose"
                )


def _superposition(
    simulations: Sequence[_Simulations], collector: BatteryCollector
) -> None:
    """Check that the response to small doses is linear and time invariant.

    The response is the difference to the simulation without a dose: twice
    the dose gives twice the response, and the doses at 0 and at
    `SECOND_DOSE` the sum of their responses; a difference is relative to the
    largest response of the species, at least `SCALE_FLOOR` of the largest
    response of all species.
    """
    later = f"at {SECOND_DOSE:g} min"
    for sims in simulations:
        runner = sims.runner
        build = runner.build
        species = list(build.tracked)
        columns = [runner.index[sid] for sid in species]
        zero = sims.zero[:, columns]
        for dose, results in sims.small.items():
            response = {k: v[:, columns] - zero for k, v in results.items()}
            comparisons = {
                "twice the dose": (response["twice at 0"], 2 * response["at 0"]),
                "two doses": (
                    response[f"at 0 and {SECOND_DOSE:g} min"],
                    response["at 0"] + response[later],
                ),
            }
            collector.cases += 1
            for name, (actual, expected) in comparisons.items():
                scale = np.max(np.abs(actual), axis=0)
                floor = SCALE_FLOOR * float(scale.max(initial=0.0))
                if floor == 0.0:
                    collector.problems.append(
                        f"{build.name} {dose}: the small doses change no species"
                    )
                    break
                difference = np.max(np.abs(actual - expected), axis=0) / np.maximum(
                    scale, floor
                )
                k = int(np.argmax(difference))
                collector.deviation(
                    float(difference[k]), f"{build.name} {dose} {name}: {species[k]}"
                )
                if difference[k] > LIMIT:
                    collector.problems.append(
                        f"{build.name} {dose}: the response of '{species[k]}' to "
                        f"{name} differs from the {'double' if name.startswith('twice') else 'sum'} "
                        f"of the single responses by {difference[k]:.2g} of its "
                        f"largest response"
                    )


_STATEMENTS: Final = {
    MASS_BALANCE.id: (
        "The amounts of the substances, of their sinks and doses change by the "
        "infusion and the flows of the reactions which read no dose alone, at "
        "most {mass:g} of the largest term."
    ),
    ZERO_DOSE.id: (
        "Without a dose the species of a substance which starts at zero stay zero."
    ),
    SUPERPOSITION.id: (
        "Twice a dose of {small:g} mg gives twice the response, and two doses "
        "at 0 and {second:g} min the sum of their responses, within {limit:g} "
        "of the largest response."
    ),
}


def _statement(check: Check) -> str:
    """Get what a check of the PK battery checks, in words."""
    return _STATEMENTS[check.id].format(
        mass=MASS_LIMIT, small=SMALL_DOSE, second=SECOND_DOSE, limit=LIMIT
    )


def fates_table(fates: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    """Get one row per build, dose and part of the fates of the doses of the battery.

    Args:
        fates: the fates as plain data, see `Fate.to_dict`.

    Returns:
        The columns of `FATE_COLUMNS`: the part in mmole and as a share of the
        dose.
    """
    rows = [
        {
            "build": f["build"],
            "dose": f["dose"],
            "substance": f["substance"],
            "amount": f["amount"],
            "part": part,
            "value": value,
            "share": value / f["amount"] if f["amount"] else math.nan,
        }
        for f in fates
        for part, value in f["parts"].items()
    ]
    return pd.DataFrame(rows, columns=list(FATE_COLUMNS))


class PKBattery(BatteryHook):
    """The battery hook of the PBPK models, see the module.

    Reads every build as a `FlatBuild`, simulates the SBML which integrates
    the flows of the reactions with the cases of the doses (`dose_cases`) and
    adds the checks `PK_CHECKS` with the fates of the doses (section `fates`).
    """

    name = "pk"
    checks = PK_CHECKS
    order = (
        "mass_balance",
        "non_negativity",
        "zero_dose",
        "steady_state",
        "superposition",
        "cross_simulator",
        "plausibility",
        "rebuild",
    )

    @override
    def settings(self) -> BatterySettings:
        """Get two days on a grid of 10 min, a model at its steady state and AMICI."""
        return BatterySettings(
            horizon=HORIZON,
            steps=STEPS,
            relative_tolerance=RELATIVE_TOLERANCE,
            absolute_tolerance=ABSOLUTE_TOLERANCE,
            steady_limit=STEADY_LIMIT,
            steady_state=True,
            engines=("amici",),
        )

    @override
    def read(self, variant: str, path: Path) -> object:
        """Read a build as a `FlatBuild`."""
        return FlatBuild.read(variant, path)

    @override
    def sbml(self, data: object) -> str | None:
        """Get the SBML of a build which integrates the flows of its reactions."""
        return data.sbml if isinstance(data, FlatBuild) else None

    @override
    def selections(self, data: object) -> Sequence[str]:
        """Select the flows of the reactions as well."""
        return data.selections if isinstance(data, FlatBuild) else ()

    @override
    def cases(self, data: object) -> Sequence[BatteryCase]:
        """Get the doses of a build, see `dose_cases`."""
        return dose_cases(data) if isinstance(data, FlatBuild) else ()

    @override
    def contribute(self, runs: Sequence[BatteryRun]) -> BatteryContribution:
        """Check the mass balance, the zero dose and the superposition of the builds."""
        builds = [run.build.data for run in runs]
        if not any(isinstance(b, FlatBuild) and b.substances for b in builds):
            reason = (
                "The model has no substance with a molecular weight `Mr_<s>`, the "
                "convention of the PBPK models the checks of the simulations read."
            )
            return BatteryContribution(
                results=tuple(not_applicable(check, reason) for check in PK_CHECKS),
                sections={"fates": []},
            )
        simulations = [_Simulations.of(run) for run in runs]
        results: dict[str, BatteryResult] = {}
        fates: list[Fate] = []
        if not any(sims.runner.build.doses for sims in simulations):
            reason = (
                "The model has no dose state or infusion rate of a substance, "
                "`IVDOSE_<s>`, `PODOSE_<s>`, `SOLDOSE_<s>` or `Ri_<s>`."
            )
            for check in (MASS_BALANCE, SUPERPOSITION):
                results[check.id] = not_applicable(check, reason)
        else:
            collector = BatteryCollector()
            _mass_balance(simulations, collector, fates)
            results[MASS_BALANCE.id] = collector.result(
                MASS_BALANCE, _statement(MASS_BALANCE)
            )
            collector = BatteryCollector()
            _superposition(simulations, collector)
            results[SUPERPOSITION.id] = collector.result(
                SUPERPOSITION, _statement(SUPERPOSITION)
            )
        collector = BatteryCollector()
        _zero_dose(simulations, collector)
        results[ZERO_DOSE.id] = (
            collector.result(ZERO_DOSE, _statement(ZERO_DOSE))
            if collector.cases
            else not_applicable(
                ZERO_DOSE,
                "Every substance of the model starts above zero, it is endogenous.",
            )
        )
        return BatteryContribution(
            results=tuple(results[check.id] for check in PK_CHECKS),
            sections={"fates": [fate.to_dict() for fate in fates]},
        )

    @override
    def describe(self) -> str:
        """Get the doses of the cases in words."""
        return (
            f"The time in min, {HORIZON / 60:g} h on a grid of {HORIZON / STEPS:g} "
            f"min; a dose of {DOSE:g} mg, an infusion over {INFUSION:g} min; the "
            f"superposition of {SMALL_DOSE:g} mg at 0 and {SECOND_DOSE:g} min."
        )

    @override
    def report_html(self, sections: Mapping[str, Any]) -> str:
        """Show the fates of the doses after the table of the battery, none without."""
        fates = list(sections.get("fates", ()))
        return render_macro("fates", fates) if fates else ""

    @override
    def report_tables(self, sections: Mapping[str, Any]) -> Mapping[str, pd.DataFrame]:
        """Write `fates.tsv`: one row per build, dose and part of the fate of a dose."""
        return {"fates": fates_table(list(sections.get("fates", ())))}


#: the battery hook of every `pkpdsim.model.PKPDModel`
PK_BATTERY: Final = PKBattery()
