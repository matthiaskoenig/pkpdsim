"""The battery of the checks of a model: its builds, not its predictions.

The verification of the predictions (`pkpdsim.credibility.verification`) asks
whether the numbers a context of use is decided on are the solution of the
model. The battery asks whether the model itself is sound, independent of any
data: every build is simulated with roadrunner without a dose and with a dose
of every route of every substance, and the results are checked against what
every PBPK model must fulfil (`Check`):

- `mass_balance`: the amount of the substances in the model, in its sinks and
  in its dose states changes by the infusion and by the flows of the reactions
  which do not read a dose (eliminations and productions) alone, so a dose
  enters with its molecular weight and nothing appears or disappears outside
  of the reactions (`MASS_LIMIT`); the fate of every dose is recorded (`Fate`).
- `non_negativity`: no species becomes negative beyond the floor of the
  integrator.
- `zero_dose`: without a dose the species of a substance which starts at
  zero (a drug, not an endogenous substance) stay zero.
- `steady_state`: without a dose every state but the accumulating ones (the
  sinks) stays at its start (`STEADY_LIMIT`): the physiology and the
  endogenous substances start at their steady state.
- `superposition`: twice a small dose gives twice the response, and two small
  doses, at 0 and at `SECOND_DOSE`, the sum of the responses to each
  (`LIMIT`): the model is linear and time invariant in the limit of a small
  dose.
- `cross_simulator`: AMICI (an optional dependency) gives the results of
  roadrunner without a dose and with every bolus dose and constant infusion
  (`LIMIT`).
- `plausibility`: the default individual has finite values and positive
  compartments, and no value of an adopted fit lies at a bound of its fit
  (`BOUND_MARGIN`).
- `rebuild`: the factory of the model writes the SBML files of the model
  again, equal but for the dates of the model history.

The battery reads the conventions of the models of pkdb_models from the flat
SBML: a substance has its molecular weight `Mr_<s>`, its doses are the states
`IVDOSE_<s>`, `PODOSE_<s>` and `SOLDOSE_<s>` (mg) and the infusion rate
`Ri_<s>` (mg/min), and a species belongs to the substance its id (the part
after the last `__`) ends with, `_<s>`, or else starts with, `<s>_` (the
lumen of the intestine of some models, `mid_lumen`). A check without what it needs, e.g. a
model without doses, is not applicable and assesses nothing.
"""

from __future__ import annotations

import dataclasses
import hashlib
import logging
import math
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum, unique
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final, Self

import numpy as np

from pkpdsim.credibility.verification import (
    BATTERY_CHECKS,
    LIMIT,
    LISTED_ISSUES,
    SCALE_FLOOR,
    Check,
)

if TYPE_CHECKING:
    import libsbml
    import roadrunner

    from pkpdsim.model import PKPDModel

__all__ = [
    "ABSOLUTE_TOLERANCE",
    "AMICI_TIMEOUT",
    "BATTERY_CHECKS",
    "BOUND_MARGIN",
    "DOSE",
    "HORIZON",
    "INFUSION",
    "MASS_LIMIT",
    "RELATIVE_TOLERANCE",
    "SECOND_DOSE",
    "SMALL_DOSE",
    "STEADY_LIMIT",
    "STEPS",
    "Battery",
    "BatteryResult",
    "BatteryStatus",
    "Fate",
    "FlatBuild",
    "run_battery",
    "same_sbml",
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
#: the longest compilation of a build by AMICI (s), see `_amici_model`
AMICI_TIMEOUT: Final = 900.0
#: the share of the range of a bound (in the scale of the fit) within which a
#: fitted value lies at the bound
BOUND_MARGIN: Final = 1e-3

#: the prefixes of the dose states (mg) and of the infusion rate (mg/min)
_BOLUS: Final = ("IVDOSE", "PODOSE", "SOLDOSE")
_INFUSION: Final = "Ri"
#: the prefixes of the sinks of a substance, by the name of their part of a fate
_SINKS: Final = {"urine": "Aurine_", "feces": "Afeces_"}
#: the prefix of the parameters which integrate the flows of the reactions
_FLOW: Final = "battery_flow_"
#: a line of the model history with a date, which changes with every build
_DATE: Final = re.compile(r"<dcterms:(created|modified)\b|<dcterms:W3CDTF>")


@unique
class BatteryStatus(StrEnum):
    """The result of a check of the battery."""

    PASSED = "passed"
    FAILED = "failed"
    NOT_APPLICABLE = "not_applicable"


@dataclass(frozen=True, kw_only=True)
class BatteryResult:
    """The result of a check of the battery over the builds of a model.

    Attributes:
        check: the check.
        status: whether it passed, failed or was not applicable.
        cases: what was checked, e.g. the doses of the mass balance over the
            builds or the species of the non-negativity.
        largest: the largest relative deviation of the check, `NaN` for a
            check without one.
        at: where the largest deviation is, e.g. `male IVDOSE_tor: Aurine_tor`.
        problems: what failed, a line each.
        statement: what was checked in words, the reason for a check which
            is not applicable.
    """

    check: Check
    status: BatteryStatus
    cases: int
    largest: float
    at: str
    problems: tuple[str, ...]
    statement: str

    @property
    def passed(self) -> bool | None:
        """Get whether the check passed, `None` if it was not applicable."""
        if self.status is BatteryStatus.NOT_APPLICABLE:
            return None
        return self.status is BatteryStatus.PASSED

    def to_dict(self) -> dict[str, object]:
        """Get the result as plain data, `NaN` as `None`."""
        return {
            "check": self.check.value,
            "status": self.status.value,
            "cases": self.cases,
            "largest": None if math.isnan(self.largest) else self.largest,
            "at": self.at,
            "problems": list(self.problems),
            "statement": self.statement,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the result from plain data, `None` as `NaN`."""
        largest = data["largest"]
        return cls(
            check=Check(data["check"]),
            status=BatteryStatus(data["status"]),
            cases=int(data["cases"]),
            largest=math.nan if largest is None else float(largest),
            at=data["at"],
            problems=tuple(data["problems"]),
            statement=data["statement"],
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


@dataclass(frozen=True, kw_only=True)
class Battery:
    """The results of the battery of a model.

    Attributes:
        results: the result of every check of `BATTERY_CHECKS`, in its order.
        fates: the fate of every dose of every build, see `Fate`.
        settings: the settings of the simulations in words.
    """

    results: tuple[BatteryResult, ...]
    fates: tuple[Fate, ...]
    settings: str

    def result(self, check: Check) -> BatteryResult:
        """Get the result of a check.

        Raises:
            KeyError: if the battery has no result of the check.
        """
        for result in self.results:
            if result.check is check:
                return result
        raise KeyError(f"The battery has no result of the check '{check}'.")

    @property
    def failed(self) -> list[BatteryResult]:
        """Get the checks which failed."""
        return [r for r in self.results if r.status is BatteryStatus.FAILED]

    def to_dict(self) -> dict[str, object]:
        """Get the battery as plain data."""
        return {
            "results": [r.to_dict() for r in self.results],
            "fates": [f.to_dict() for f in self.fates],
            "settings": self.settings,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the battery from plain data."""
        return cls(
            results=tuple(BatteryResult.from_dict(r) for r in data["results"]),
            fates=tuple(Fate.from_dict(f) for f in data["fates"]),
            settings=data["settings"],
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
        """Get the selections of a simulation: the time, the states and the flows.

        A species is selected as its amount, a parameter by its id.
        """
        return [
            "time",
            *self.species,
            *self.states,
            *(f"{_FLOW}{rid}" for rid in self.flows),
        ]


def same_sbml(a: str, b: str) -> bool:
    """Compare two SBML files without the dates of their model history.

    Args:
        a: the content of one file.
        b: the content of the other.

    Returns:
        Whether they are equal but for the lines with a date.
    """

    def lines(text: str) -> list[str]:
        return [line for line in text.splitlines() if not _DATE.search(line)]

    return lines(a) == lines(b)


class _Runner:
    """The simulations of a build with roadrunner on the grid of the battery."""

    def __init__(self, build: FlatBuild) -> None:
        """Load a build at `O0` with the tolerances of the battery."""
        from sbmlsim.model.model_roadrunner import load_roadrunner

        start = time.perf_counter()
        self.build = build
        self.r: roadrunner.RoadRunner = load_roadrunner(build.sbml, "O0")
        self.r.integrator.relative_tolerance = RELATIVE_TOLERANCE
        self.r.integrator.absolute_tolerance = ABSOLUTE_TOLERANCE
        self.columns = build.selections[1:]
        self.r.timeCourseSelections = build.selections
        self.index = {sid: k for k, sid in enumerate(self.columns)}
        self.r.resetAll()
        self.values = {
            pid: float(self.r[pid]) for pid in self.r.model.getGlobalParameterIds()
        }
        logger.info(
            "Load the %s build '%s' in %.1f s",
            build.name,
            build.path.name,
            time.perf_counter() - start,
        )

    def run(self, changes: Mapping[float, Mapping[str, float]]) -> np.ndarray:
        """Simulate from the start with changes at times of the grid.

        A dose state is increased by its change, any other id set to it. The
        row of a change time is the state before the change, except at 0.

        Args:
            changes: the changes by their time (min).

        Returns:
            The results on the grid, `STEPS + 1` rows of the selections
            without the time.
        """
        r = self.r
        r.resetAll()
        times = sorted({0.0, *changes})
        bounds = [*times, HORIZON]
        step = HORIZON / STEPS
        rows: list[np.ndarray] = []
        for k, start in enumerate(times):
            for sid, value in changes.get(start, {}).items():
                r[sid] = float(r[sid]) + value if sid.split("_")[0] in _BOLUS else value
            end = bounds[k + 1]
            points = round((end - start) / step) + 1
            result = np.asarray(r.simulate(start, end, points))[:, 1:]
            rows.append(result if k == 0 else result[1:])
        return np.vstack(rows)

    def column(self, result: np.ndarray, sid: str) -> np.ndarray:
        """Get the values of an id over the grid."""
        return result[:, self.index[sid]]


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


@dataclass
class _Collector:
    """The cases and problems of a check over the builds."""

    cases: int = 0
    largest: float = math.nan
    at: str = ""
    problems: list[str] = dataclasses.field(default_factory=list)

    def deviation(self, value: float, at: str) -> None:
        """Record a deviation, the largest one is kept, where only above zero."""
        if math.isnan(self.largest) or value > self.largest:
            self.largest, self.at = value, at if value > 0 else ""

    def result(self, check: Check, statement: str) -> BatteryResult:
        """Get the result of the check, failed with a problem."""
        problems = self.problems[:LISTED_ISSUES]
        if len(self.problems) > LISTED_ISSUES:
            problems.append(f"and {len(self.problems) - LISTED_ISSUES} more")
        return BatteryResult(
            check=check,
            status=BatteryStatus.FAILED if self.problems else BatteryStatus.PASSED,
            cases=self.cases,
            largest=self.largest,
            at=self.at,
            problems=tuple(problems),
            statement=statement,
        )


def _not_applicable(check: Check, reason: str) -> BatteryResult:
    """Get the result of a check which is not applicable."""
    return BatteryResult(
        check=check,
        status=BatteryStatus.NOT_APPLICABLE,
        cases=0,
        largest=math.nan,
        at="",
        problems=(),
        statement=reason,
    )


@dataclass
class _Simulations:
    """The simulations of a build: without a dose, with a dose and the superposition."""

    runner: _Runner
    zero: np.ndarray
    dosed: dict[str, np.ndarray]
    small: dict[str, dict[str, np.ndarray]]

    @classmethod
    def run(cls, runner: _Runner) -> _Simulations:
        """Run the simulations of a build."""
        start = time.perf_counter()
        zero = runner.run({})
        dosed: dict[str, np.ndarray] = {}
        small: dict[str, dict[str, np.ndarray]] = {}
        for dose in runner.build.doses:
            dosed[dose] = runner.run(_dose_changes(dose, DOSE, 0.0))
            first = _dose_changes(dose, SMALL_DOSE, 0.0)
            second = _dose_changes(dose, SMALL_DOSE, SECOND_DOSE)
            small[dose] = {
                "at 0": runner.run(first),
                "twice at 0": runner.run(_dose_changes(dose, 2 * SMALL_DOSE, 0.0)),
                f"at {SECOND_DOSE:g} min": runner.run(second),
                f"at 0 and {SECOND_DOSE:g} min": runner.run(_merge(first, second)),
            }
        logger.info(
            "Simulate the %s build '%s' in %.1f s",
            runner.build.name,
            runner.build.path.name,
            time.perf_counter() - start,
        )
        return cls(runner=runner, zero=zero, dosed=dosed, small=small)


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
    simulations: Sequence[_Simulations], collector: _Collector, fates: list[Fate]
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


def _results(sims: _Simulations) -> list[tuple[str, np.ndarray]]:
    """Get every simulation of a build by its name."""
    named = [("without a dose", sims.zero)]
    named += [(dose, result) for dose, result in sims.dosed.items()]
    for dose, results in sims.small.items():
        named += [(f"{dose} small {at}", result) for at, result in results.items()]
    return named


def _non_negativity(simulations: Sequence[_Simulations], collector: _Collector) -> None:
    """Check that no species becomes negative beyond the floor."""
    for sims in simulations:
        runner = sims.runner
        species = list(runner.build.species)
        columns = [runner.index[sid] for sid in species]
        named = _results(sims)
        stacked = np.vstack([result[:, columns] for _, result in named])
        scale = np.max(np.abs(stacked), axis=0)
        floor = np.maximum(SCALE_FLOOR * scale, ABSOLUTE_TOLERANCE)
        collector.cases += len(species)
        for name, result in named:
            low = result[:, columns].min(axis=0)
            relative = np.maximum(-low, 0.0) / np.maximum(scale, ABSOLUTE_TOLERANCE)
            if len(species):
                k = int(np.argmax(relative))
                collector.deviation(
                    float(relative[k]), f"{runner.build.name} {name}: {species[k]}"
                )
            for k in np.flatnonzero(low < -floor):
                collector.problems.append(
                    f"{runner.build.name} {name}: '{species[k]}' reaches "
                    f"{low[k]:.3g}, below the floor {-floor[k]:.3g}"
                )


def _scales(sims: _Simulations, ids: Sequence[str]) -> np.ndarray:
    """Get the largest absolute value of every id over the simulations of a build."""
    columns = [sims.runner.index[sid] for sid in ids]
    return np.max(
        np.abs(np.vstack([result[:, columns] for _, result in _results(sims)])),
        axis=0,
        initial=0.0,
    )


def _floor(scale: np.ndarray) -> np.ndarray:
    """Get the floor of a deviation of an id: `SCALE_FLOOR` of its scale."""
    return np.maximum(SCALE_FLOOR * scale, ABSOLUTE_TOLERANCE)


def _zero_dose(simulations: Sequence[_Simulations], collector: _Collector) -> None:
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


def _steady_state(simulations: Sequence[_Simulations], collector: _Collector) -> None:
    """Check that every state but the sinks stays at its start without a dose.

    A state may change by `STEADY_LIMIT` of its start plus the floor of its
    largest value over the simulations of its build, so a state which
    starts at zero may move by the noise of the integrator.
    """
    for sims in simulations:
        runner = sims.runner
        build = runner.build
        ids = [
            sid
            for sid in (*build.species, *build.states)
            if sid not in build.accumulating
        ]
        if not ids:
            continue
        values = sims.zero[:, [runner.index[sid] for sid in ids]]
        start = values[0]
        change = np.max(np.abs(values - start), axis=0)
        allowed = STEADY_LIMIT * np.abs(start) + _floor(_scales(sims, ids))
        relative = np.where(
            start != 0.0, change / np.maximum(np.abs(start), ABSOLUTE_TOLERANCE), 0.0
        )
        collector.cases += len(ids)
        k = int(np.argmax(relative))
        collector.deviation(float(relative[k]), f"{build.name}: {ids[k]}")
        for k in np.flatnonzero(change > allowed):
            collector.problems.append(
                f"{build.name}: '{ids[k]}' changes by {change[k]:.3g} from "
                f"{start[k]:.3g} without a dose"
            )


def _superposition(simulations: Sequence[_Simulations], collector: _Collector) -> None:
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


def _plausibility(
    model: PKPDModel, runners: Sequence[_Runner], collector: _Collector
) -> None:
    """Check the values of the default individual and of the adopted fits."""
    for runner in runners:
        r = runner.r
        r.resetAll()
        name = runner.build.name
        for cid in r.model.getCompartmentIds():
            size = float(r[cid])
            collector.cases += 1
            if not math.isfinite(size) or size <= 0.0:
                collector.problems.append(
                    f"{name}: the compartment '{cid}' has the size {size:g}"
                )
        for pid in r.model.getGlobalParameterIds():
            value = float(r[pid])
            collector.cases += 1
            if not math.isfinite(value):
                collector.problems.append(f"{name}: the parameter '{pid}' is {value}")
        for sid in r.model.getFloatingSpeciesIds():
            value = float(r[sid])
            collector.cases += 1
            if not math.isfinite(value) or value < 0.0:
                collector.problems.append(
                    f"{name}: the species '{sid}' starts at {value:g}"
                )
    if not model.adopted_ids:
        return
    try:
        records = model.fit_records()
    except (FileNotFoundError, ValueError) as err:
        collector.problems.append(str(err))
        return
    for pid, parameter in records.parameters.items():
        collector.cases += 1
        value = parameter._to_scale(parameter.value)
        for side, bound in (
            ("lower", parameter.lower_bound),
            ("upper", parameter.upper_bound),
        ):
            if not math.isfinite(bound):
                continue
            other = parameter.upper_bound if side == "lower" else parameter.lower_bound
            if not math.isfinite(other):
                continue
            span = abs(
                parameter._to_scale(parameter.upper_bound)
                - parameter._to_scale(parameter.lower_bound)
            )
            distance = (
                abs(value - parameter._to_scale(bound)) / span if span > 0 else 0.0
            )
            if distance < BOUND_MARGIN:
                collector.problems.append(
                    f"the value {parameter.value:g} of '{pid}' of the fit "
                    f"'{parameter.fit_id}' lies at its {side} bound {bound:g}"
                )


def _rebuild(model: PKPDModel) -> BatteryResult:
    """Check that the factory writes the SBML files of the model again."""
    statement = (
        "The factory writes the SBML files of the model again, equal but for "
        "the dates of the model history."
    )
    collector = _Collector()
    with tempfile.TemporaryDirectory(prefix="pkpdsim-rebuild-") as tmp:
        target = Path(tmp) / model.factory_dir.name
        start = time.perf_counter()
        try:
            model.load_factory()(target)
        except Exception as err:  # a failure of a factory of a model is the result
            collector.problems.append(
                f"the factory failed: {type(err).__name__}: {err}"
            )
            return collector.result(Check.REBUILD, statement)
        logger.info("Rebuild '%s' in %.1f s", model.name, time.perf_counter() - start)
        rebuilt = {p.relative_to(target): p for p in target.rglob("*.xml")}
        for path in model.model_paths.values():
            if path.relative_to(model.factory_dir) not in rebuilt:
                collector.problems.append(
                    f"the factory did not write the build '{path.name}'"
                )
        for relative, path in sorted(rebuilt.items()):
            existing = model.factory_dir / relative
            if not existing.exists():
                continue
            collector.cases += 1
            if not same_sbml(
                existing.read_text(encoding="utf-8"), path.read_text(encoding="utf-8")
            ):
                collector.problems.append(f"'{relative}' differs from its rebuild")
    return collector.result(Check.REBUILD, statement)


#: the compilation of a build by AMICI, in a process of its own
_COMPILE: Final = (
    "import sys\n"
    "from amici.importers.sbml import SbmlImporter\n"
    "SbmlImporter(sys.argv[1]).sbml2amici(sys.argv[2], sys.argv[3], verbose=False)\n"
)


class _CompileTimeout(Exception):
    """AMICI did not compile a build within `AMICI_TIMEOUT`."""


def _amici_model(build: FlatBuild) -> Any:  # noqa: ANN401  # amici is optional and untyped
    """Compile a build with AMICI once, cached by the SBML and the AMICI version.

    The compilation runs in a process of its own, stopped after
    `AMICI_TIMEOUT`: the symbolic code generation of AMICI takes a minute and a
    half for a PBPK model, but more than 45 min for the large nested formulas
    of the starts of the creatinine model of pkdb_models.

    Raises:
        _CompileTimeout: if the compilation takes longer.
        RuntimeError: if the compilation fails.
    """
    import amici

    content = build.path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()[:16]
    name = f"pkpdsim_{digest}"
    directory = Path(tempfile.gettempdir()) / "pkpdsim-amici" / amici.__version__ / name
    if not (directory / ".complete").exists():
        shutil.rmtree(directory, ignore_errors=True)
        start = time.perf_counter()
        try:
            result = subprocess.run(
                [sys.executable, "-c", _COMPILE, str(build.path), name, str(directory)],
                capture_output=True,
                text=True,
                timeout=AMICI_TIMEOUT,
                check=False,
            )
        except subprocess.TimeoutExpired as err:
            shutil.rmtree(directory, ignore_errors=True)
            raise _CompileTimeout from err
        if result.returncode != 0:
            last = (result.stderr.strip().splitlines() or ["no output"])[-1]
            raise RuntimeError(last)
        (directory / ".complete").touch()
        logger.info(
            "Compile the %s build '%s' with AMICI in %.1f s",
            build.name,
            build.path.name,
            time.perf_counter() - start,
        )
    return amici.import_model_module(name, directory).get_model()


def _cross_simulator(runners: Sequence[_Runner]) -> BatteryResult:
    """Compare the results of roadrunner with the results of AMICI."""
    statement = (
        "AMICI gives the results of roadrunner without a dose and with every "
        "bolus dose and constant infusion."
    )
    try:
        from amici.sim.sundials import run_simulation
    except ImportError:
        return _not_applicable(
            Check.CROSS_SIMULATOR,
            "AMICI is not installed, install pkpdsim with the extra 'amici'.",
        )
    collector = _Collector()
    skipped: list[str] = []
    times = np.linspace(0.0, HORIZON, STEPS + 1)
    for runner in runners:
        build = runner.build
        try:
            model = _amici_model(build)
        except _CompileTimeout:
            skipped.append(build.name)
            continue
        except Exception as err:  # a model AMICI cannot compile is the result
            collector.problems.append(
                f"{build.name}: AMICI cannot compile the build: {type(err).__name__}: {err}"
            )
            continue
        ids = list(model.get_state_ids())
        r = runner.r
        floating = set(r.model.getFloatingSpeciesIds())
        selections = [
            "time",
            *(f"[{i}]" if i in floating and i not in build.amounts else i for i in ids),
        ]
        model.set_timepoints(times)
        solver = model.create_solver()
        solver.set_relative_tolerance(RELATIVE_TOLERANCE)
        solver.set_absolute_tolerance(ABSOLUTE_TOLERANCE)
        solver.set_max_steps(10**6)
        default = np.asarray(model.get_initial_state(), dtype=float)
        cases: list[tuple[str, dict[str, float]]] = [("without a dose", {})]
        cases += [
            (dose, {dose: DOSE if dose.split("_")[0] in _BOLUS else DOSE / INFUSION})
            for dose in build.doses
        ]
        for case, changes in cases:
            collector.cases += 1
            initial = default.copy()
            parameters: dict[str, float] = {}
            for sid, value in changes.items():
                if sid in ids:
                    initial[ids.index(sid)] += value
                else:
                    parameters[sid] = value
            model.set_initial_state(initial)
            for pid, value in parameters.items():
                model.set_free_parameter_by_id(
                    pid, value
                ) if pid in model.get_free_parameter_ids() else model.set_fixed_parameter_by_id(
                    pid, value
                )
            rdata = run_simulation(model, solver)
            for pid in parameters:
                model.set_free_parameter_by_id(
                    pid, 0.0
                ) if pid in model.get_free_parameter_ids() else model.set_fixed_parameter_by_id(
                    pid, 0.0
                )
            model.set_initial_state(default)
            if rdata.status != 0:
                collector.problems.append(
                    f"{build.name} {case}: AMICI failed with the status {rdata.status}"
                )
                continue
            r.resetAll()
            for sid, value in changes.items():
                r[sid] = float(r[sid]) + value if sid.split("_")[0] in _BOLUS else value
            r.timeCourseSelections = selections
            expected = np.asarray(r.simulate(0.0, HORIZON, STEPS + 1))[:, 1:]
            r.timeCourseSelections = build.selections
            actual = np.asarray(rdata.x, dtype=float)
            scale = np.max(np.abs(expected), axis=0)
            floor = np.maximum(SCALE_FLOOR * scale, ABSOLUTE_TOLERANCE)
            with np.errstate(invalid="ignore", divide="ignore"):
                difference = np.abs(actual - expected) / np.maximum(
                    np.maximum(np.abs(expected), np.abs(actual)), floor
                )
            difference = np.where(actual == expected, 0.0, difference)
            largest = np.nanmax(
                np.where(np.isfinite(difference), difference, np.inf), axis=0
            )
            k = int(np.argmax(largest))
            collector.deviation(float(largest[k]), f"{build.name} {case}: {ids[k]}")
            if largest[k] > LIMIT:
                collector.problems.append(
                    f"{build.name} {case}: '{ids[k]}' differs between AMICI and "
                    f"roadrunner by {largest[k]:.2g}"
                )
    if skipped:
        late = (
            f"AMICI did not compile the {' and '.join(skipped)} "
            f"{'build' if len(skipped) == 1 else 'builds'} within "
            f"{AMICI_TIMEOUT / 60:g} min."
        )
        if not collector.cases and not collector.problems:
            return _not_applicable(Check.CROSS_SIMULATOR, late)
        statement = f"{statement} {late}"
    return collector.result(Check.CROSS_SIMULATOR, statement)


_STATEMENTS: Final = {
    Check.MASS_BALANCE: (
        "The amounts of the substances, of their sinks and doses change by the "
        "infusion and the flows of the reactions which read no dose alone, at "
        "most {mass:g} of the largest term."
    ),
    Check.NON_NEGATIVITY: (
        "No species becomes negative beyond {floor:g} of its largest value."
    ),
    Check.ZERO_DOSE: (
        "Without a dose the species of a substance which starts at zero stay zero."
    ),
    Check.STEADY_STATE: (
        "Without a dose every state but the sinks stays at its start, within "
        "{steady:g} of it."
    ),
    Check.SUPERPOSITION: (
        "Twice a dose of {small:g} mg gives twice the response, and two doses "
        "at 0 and {second:g} min the sum of their responses, within {limit:g} "
        "of the largest response."
    ),
    Check.PLAUSIBILITY: (
        "The default individual has finite values, positive compartments and "
        "species which are not negative, and no fitted value lies within "
        "{margin:g} of the range of a bound."
    ),
}


def _statement(check: Check) -> str:
    """Get what a check of the simulations checks, in words."""
    return _STATEMENTS[check].format(
        mass=MASS_LIMIT,
        floor=SCALE_FLOOR,
        steady=STEADY_LIMIT,
        small=SMALL_DOSE,
        second=SECOND_DOSE,
        limit=LIMIT,
        margin=BOUND_MARGIN,
    )


def settings() -> str:
    """Get the settings of the simulations of the battery in words."""
    return (
        f"roadrunner CVODE BDF, rtol {RELATIVE_TOLERANCE:g}, atol "
        f"{ABSOLUTE_TOLERANCE:g}; every simulation {HORIZON / 60:g} h on a grid of "
        f"{HORIZON / STEPS:g} min; a dose of {DOSE:g} mg, an infusion over "
        f"{INFUSION:g} min; the superposition of {SMALL_DOSE:g} mg at 0 and "
        f"{SECOND_DOSE:g} min."
    )


#: the checks of the simulations, which need a substance with a molecular weight
_SIMULATION_CHECKS: Final = (
    Check.MASS_BALANCE,
    Check.NON_NEGATIVITY,
    Check.ZERO_DOSE,
    Check.STEADY_STATE,
    Check.SUPERPOSITION,
    Check.CROSS_SIMULATOR,
)


def _simulation_checks(
    builds: Sequence[FlatBuild], runners: Sequence[_Runner], *, cross_simulator: bool
) -> tuple[dict[Check, BatteryResult], list[Fate]]:
    """Run the checks of the simulations of the builds."""
    simulations = [_Simulations.run(runner) for runner in runners]
    results: dict[Check, BatteryResult] = {}
    fates: list[Fate] = []
    if not any(build.doses for build in builds):
        reason = (
            "The model has no dose state or infusion rate of a substance, "
            "`IVDOSE_<s>`, `PODOSE_<s>`, `SOLDOSE_<s>` or `Ri_<s>`."
        )
        for check in (Check.MASS_BALANCE, Check.SUPERPOSITION):
            results[check] = _not_applicable(check, reason)
    else:
        collector = _Collector()
        _mass_balance(simulations, collector, fates)
        results[Check.MASS_BALANCE] = collector.result(
            Check.MASS_BALANCE, _statement(Check.MASS_BALANCE)
        )
        collector = _Collector()
        _superposition(simulations, collector)
        results[Check.SUPERPOSITION] = collector.result(
            Check.SUPERPOSITION, _statement(Check.SUPERPOSITION)
        )
    for check, function, nothing in (
        (Check.NON_NEGATIVITY, _non_negativity, "The model has no species."),
        (
            Check.ZERO_DOSE,
            _zero_dose,
            "Every substance of the model starts above zero, it is endogenous.",
        ),
        (Check.STEADY_STATE, _steady_state, "The model has no state."),
    ):
        collector = _Collector()
        function(simulations, collector)
        results[check] = (
            collector.result(check, _statement(check))
            if collector.cases
            else _not_applicable(check, nothing)
        )
    results[Check.CROSS_SIMULATOR] = (
        _cross_simulator(runners)
        if cross_simulator
        else _not_applicable(Check.CROSS_SIMULATOR, "Not run.")
    )
    return results, fates


def run_battery(
    model: PKPDModel, *, rebuild: bool = True, cross_simulator: bool = True
) -> Battery:
    """Run the battery of the checks of a model on its builds.

    The checks of the simulations need a substance with a molecular weight,
    `Mr_<s>`, the convention of the PBPK models; they are not applicable to a
    model without.

    Args:
        model: the model, every build of `model_paths`.
        rebuild: run the check `rebuild`, which runs the factory.
        cross_simulator: run the check `cross_simulator` when AMICI is
            installed, which compiles every build.

    Returns:
        The results of the checks and the fates of the doses.
    """
    builds = [
        FlatBuild.read(sex.value, path) for sex, path in model.model_paths.items()
    ]
    runners = [_Runner(build) for build in builds]
    if any(build.substances for build in builds):
        results, fates = _simulation_checks(
            builds, runners, cross_simulator=cross_simulator
        )
    else:
        reason = (
            "The model has no substance with a molecular weight `Mr_<s>`, the "
            "convention of the PBPK models the checks of the simulations read."
        )
        results = {
            check: _not_applicable(check, reason) for check in _SIMULATION_CHECKS
        }
        fates = []
    collector = _Collector()
    _plausibility(model, runners, collector)
    results[Check.PLAUSIBILITY] = collector.result(
        Check.PLAUSIBILITY, _statement(Check.PLAUSIBILITY)
    )
    results[Check.REBUILD] = (
        _rebuild(model) if rebuild else _not_applicable(Check.REBUILD, "Not run.")
    )
    return Battery(
        results=tuple(results[check] for check in BATTERY_CHECKS),
        fates=tuple(fates),
        settings=settings(),
    )
