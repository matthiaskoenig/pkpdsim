"""What the library reads from a flat model file, parsed and compiled once.

The equilibration of a condition (`pkpdsim.pbpk.equilibration`) and the
validity of the kidney outflow (`pkpdsim.pbpk.kidney_validity`) read the
same built model: its species and parameters, the states which nothing
drains, the reactions, the filtered substances and the ids the flows of the
kidney depend on. `flat_model` parses the file once per process (by its path
and modification time, so a file which is built again is read again) and
gives all of it as a `FlatModel`; `loaded_model` compiles the model once
(roadrunner at LLJIT `O0`, the fastest compilation, since its simulations
only reach a steady state or evaluate the initial assignments) with the
selections of both, so a model with two builds compiles each build once for
both checks. `math_names` is the one walker of the names a formula of the
model reads.

The module imports sbmlsim and libsbml; it is internal to the library.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import libsbml
from sbmlsim.model.model_roadrunner import RoadrunnerSBMLModel

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.organs import Organ, Part


def math_names(ast: libsbml.ASTNode | None) -> set[str]:
    """Get the names a formula of a model reads.

    Args:
        ast: the formula, `None` for none.

    Returns:
        The names of the formula (`AST_NAME`), not the time.
    """
    found: set[str] = set()
    stack = [] if ast is None else [ast]
    while stack:
        node = stack.pop()
        if node.getType() == libsbml.AST_NAME:
            found.add(node.getName())
        stack.extend(node.getChild(k) for k in range(node.getNumChildren()))
    return found


def accumulating(model: libsbml.Model) -> tuple[str, ...]:
    """Get the states of a parsed model which nothing drains, sorted.

    See `pkpdsim.pbpk.equilibration.accumulating_states`.
    """
    assigned: dict[str, set[str]] = {}
    rated: dict[str, set[str]] = {}
    for rule in model.getListOfRules():
        if rule.isAssignment():
            assigned[rule.getVariable()] = math_names(rule.getMath())
        elif rule.isRate():
            rated[rule.getVariable()] = math_names(rule.getMath())
    read = set().union(*rated.values())
    consumed: set[str] = set()
    changed: set[str] = set()
    for reaction in model.getListOfReactions():
        law = reaction.getKineticLaw()
        if law is not None:
            local = {p.getId() for p in law.getListOfLocalParameters()}
            read |= math_names(law.getMath()) - local
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


@dataclass(frozen=True)
class FilteredSubstance:
    """A filtered substance of a flat model.

    Attributes:
        sid: the substance.
        kinetic: its erythrocytes are kinetic: the plasma carries it, else
            `BP` does, which is `1 - HCT` without erythrocytes.
        threshold: it has a threshold reabsorption (`reab_threshold_<s>`),
            no `f_reabsorbed_<s>`.
    """

    sid: str
    kinetic: bool
    threshold: bool


@dataclass(frozen=True)
class FlatModel:
    """What the library reads from a flat model file.

    Attributes:
        species: the ids of the species.
        parameters: the ids of the parameters.
        accumulating: the states which nothing drains (`accumulating`).
        candidates: the species of the substances without a baseline which
            start at 0 (no initial assignment, no rule, not constant), see
            `pkpdsim.pbpk.equilibration.drug_free_zeros`.
        started: the species with a start: an initial value other than 0,
            an initial assignment or a rule.
        concentrations: the species whose id means their concentration
            (`hasOnlySubstanceUnits` false): a change sets them as `[sid]`.
        reactions: (reactants, products, reversible, the names its rate
            reads) of every reaction.
        filtered: the filtered substances of the kidney.
        water: the model has the GFR, the flow of the kidney and the
            hematocrit, whose filtered plasma water is bounded.
        quantities: the ids of the flows of the validity of the kidney
            outflow and `GFR_base` which the model has.
        inputs: the ids which the quantities depend on, through initial
            assignments and assignment rules (themselves included): the
            changes of an individual.
    """

    species: frozenset[str]
    parameters: frozenset[str]
    accumulating: tuple[str, ...]
    candidates: frozenset[str]
    started: frozenset[str]
    concentrations: frozenset[str]
    reactions: tuple[tuple[frozenset[str], frozenset[str], bool, frozenset[str]], ...]
    filtered: tuple[FilteredSubstance, ...]
    water: bool
    quantities: tuple[str, ...]
    inputs: frozenset[str]

    @property
    def selections(self) -> tuple[str, ...]:
        """The selections of the loaded model: the time, the states and the flows."""
        return ("time", *dict.fromkeys((*self.accumulating, *self.quantities)))


def flat_model(model_path: Path) -> FlatModel:
    """Get what the library reads from a flat model file, parsed once per process.

    Args:
        model_path: the flat SBML model.

    Returns:
        The parsed model, cached by the path and the modification time.
    """
    path = Path(model_path)
    return _read(path, path.stat().st_mtime_ns)


def loaded_model(model_path: Path) -> RoadrunnerSBMLModel:
    """Get a flat model file compiled once per process, at `O0`.

    The selections are those of the equilibration and of the validity of the
    kidney outflow (`FlatModel.selections`); the four models used last stay
    loaded.

    Args:
        model_path: the flat SBML model.

    Returns:
        The loaded model, cached by the path and the modification time.
    """
    path = Path(model_path)
    return _loaded(path, path.stat().st_mtime_ns)


@lru_cache(maxsize=4)
def _loaded(model_path: Path, mtime_ns: int) -> RoadrunnerSBMLModel:
    """Load a model, `mtime_ns` a key of the cache."""
    return RoadrunnerSBMLModel(
        source=model_path,
        selections=list(_read(model_path, mtime_ns).selections),
        jit_optimization="O0",
    )


@lru_cache(maxsize=8)
def _read(model_path: Path, mtime_ns: int) -> FlatModel:
    """Parse a flat model file, `mtime_ns` a key of the cache."""
    model = libsbml.readSBMLFromFile(str(model_path)).getModel()
    parameters = {p.getId() for p in model.getListOfParameters()}
    species = list(model.getListOfSpecies())
    all_ids = [sp.getId() for sp in species]
    candidates, started = _starts(model, parameters, species, all_ids)
    reactions = []
    for reaction in model.getListOfReactions():
        law = reaction.getKineticLaw()
        reactions.append(
            (
                frozenset(r.getSpecies() for r in reaction.getListOfReactants()),
                frozenset(r.getSpecies() for r in reaction.getListOfProducts()),
                reaction.getReversible(),
                frozenset(() if law is None else math_names(law.getMath())),
            )
        )
    filtered, water, quantities = _kidney(model, set(all_ids))
    return FlatModel(
        species=frozenset(all_ids),
        parameters=frozenset(parameters),
        accumulating=accumulating(model),
        candidates=frozenset(candidates),
        started=frozenset(started),
        concentrations=frozenset(
            sp.getId() for sp in species if not sp.getHasOnlySubstanceUnits()
        ),
        reactions=tuple(reactions),
        filtered=filtered,
        water=water,
        quantities=quantities,
        inputs=_dependencies(model, quantities),
    )


def _starts(
    model: libsbml.Model,
    parameters: set[str],
    species: list[libsbml.Species],
    all_ids: list[str],
) -> tuple[set[str], set[str]]:
    """The species of the drugs which start at 0, and the species with a start."""
    substances = {
        pid.removeprefix(ids.mr("")) for pid in parameters if pid.startswith(ids.mr(""))
    }
    exogenous = {s for s in substances if ids.baseline(s) not in parameters}
    assigned = {a.getSymbol() for a in model.getListOfInitialAssignments()}
    ruled = {r.getVariable() for r in model.getListOfRules()}
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
    return candidates, started


def _kidney(
    model: libsbml.Model, species: set[str]
) -> tuple[tuple[FilteredSubstance, ...], bool, tuple[str, ...]]:
    """The filtered substances, the water and the flows of the kidney of a model."""
    prefix = ids.flat(Organ.KIDNEY, ids.filtration_reaction(""))
    filtered = tuple(
        FilteredSubstance(
            s,
            ids.flat(Organ.KIDNEY, ids.species(Part.BLOOD_CELLS, s)) in species,
            model.getParameter(ids.reabsorption_threshold(s)) is not None,
        )
        for s in (
            reaction.getId().removeprefix(prefix)
            for reaction in model.getListOfReactions()
            if reaction.getId().startswith(prefix)
        )
    )
    quantities = [ids.GFR, ids.GFR_BASE, ids.organ_flow(Organ.KIDNEY), ids.HCT]
    for f in filtered:
        quantities += [ids.fu(f.sid), ids.bp(f.sid), ids.fraction_reabsorbed(f.sid)]
    present = tuple(q for q in quantities if model.getParameter(q) is not None)
    water = {ids.GFR, ids.organ_flow(Organ.KIDNEY), ids.HCT} <= set(present)
    return filtered, water, present


def _dependencies(model: libsbml.Model, roots: Iterable[str]) -> frozenset[str]:
    """The ids which the roots depend on through initial assignments and rules."""
    reads: dict[str, set[str]] = {}
    for assignment in model.getListOfInitialAssignments():
        reads.setdefault(assignment.getSymbol(), set()).update(
            math_names(assignment.getMath())
        )
    for rule in model.getListOfRules():
        if rule.isAssignment():
            reads.setdefault(rule.getVariable(), set()).update(
                math_names(rule.getMath())
            )
    seen: set[str] = set()
    stack = list(roots)
    while stack:
        sid = stack.pop()
        if sid not in seen:
            seen.add(sid)
            stack.extend(reads.get(sid, ()))
    return frozenset(seen)
