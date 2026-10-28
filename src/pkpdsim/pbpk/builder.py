"""The composition of a PBPK model from its definition, and its build.

The comp model has a top model and one submodel per organ:

- the top model holds the physiology (covariates, volumes, flows), the blood
  pools (arterial, venous, portal; plasma `C<pool>_<s>` and kinetic blood
  cells `C<pool>_bc_<s>`), the parameters of the substances (`Mr`, `fu`, `Kp`
  or the factor `fKp_<s>` and the `Kp_<o>_<s>` of a Kp method, `BP`, ...), the sinks (urine, feces, bile), the dose states of the routes
  and the parameters the intestine reads (`ka_<s>`, `FED`);
- every organ is a `ModelDefinition` with a `Submodel` (`LI`, `KI`, ...), the
  standard organ or a variant (`pkpdsim.pbpk.variants`); it reads the
  elements of the body through its ports, and its outflow concentrations
  replace the parameters `Cout_<o>_<s>` of the top model. The intestine
  (`segmented` by default, or `single_lumen`) replaces small and large
  intestine with one submodel `GU`, which owns the intestinal lumen and
  applies the oral doses.

Every port is declared by the organ (its `Boundary`) and every replacement is
generated here from it; the flat model has the ids of the scheme
(`LI__Cpl_<s>` in the organs, the ids of the top model unprefixed).

The sinks are those of `pkpdsim.pbpk.sinks`, the observers (unbound and blood concentrations, hepatic and arm vein) those
of `pkpdsim.pbpk.observers`.

Every quantity which depends on parameters only (the physiology, the
impairment, the compartments of the pools and the organs, `fu`, `BP`, the
coefficients and exchanges of the organs) is a constant with an initial
assignment, not an assignment rule: a change before the simulation (a change
of an experiment, which sbmlsim applies with the initial assignments which
read it) is followed, a change during the simulation is not. Only what reads
a species (the outflow concentrations `Cout_<o>_<s>`, the observers) and the
rate of an intravenous injection `Kiv_<s>` (of the dose protocol, which
changes it between doses) are assignment rules. A parameter which changes
during a simulation is named in `PBPKDefinition.dynamic`: every quantity
downstream of it, in the top model and in the organs, is an assignment rule
(`_dynamic`).

The top model holds the impaired organ function as well
(`pkpdsim.pbpk.impairment`): the severities and their factors on the
physiology, `fu` and the activity of the proteins of the liver.

A PD component (`pkpdsim.pbpk.pd`) is a further `ModelDefinition` with a
`Submodel` of its id: the elements of the body replace its inputs, its
outputs replace parameters of the body, which are dynamic.

Every element is annotated when it is created (`pkpdsim.pbpk.annotation`);
an input of a component takes the annotation of the element of the body
which replaces it, and the model the annotation of a PBPK model.
"""

import copy
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from pathlib import Path

from sbmlutils.comp import flatten_sbml
from sbmlutils.factory import (
    PORT_SUFFIX,
    AssignmentRule,
    Compartment,
    InitialAssignment,
    Model,
    ModelDefinition,
    Parameter,
    Port,
    PortType,
    RateRule,
    Reaction,
    ReplacedBy,
    ReplacedElement,
    Sbase,
    Species,
    Submodel,
    UnitDefinition,
    create_model,
)
from sbmlutils.factory.units import UnitType
from sbmlutils.metadata import SBO
from sbmlutils.validation import validate_doc

from pkpdsim.pbpk import distribution, ids, variants
from pkpdsim.pbpk.annotation import annotate, initial_assignment, model_annotations
from pkpdsim.pbpk.chemistry import Membrane, Metabolism, Transport
from pkpdsim.pbpk.data import cite
from pkpdsim.pbpk.definition import PBPKDefinition
from pkpdsim.pbpk.impairment import (
    impairment_objects,
    physiology_factors,
    unbound_fraction,
    unbound_fraction_notes,
)
from pkpdsim.pbpk.kp import kp_sources, kp_values
from pkpdsim.pbpk.observers import (
    TIME,
    definition_observer_objects,
    formula_symbols,
    observer_objects,
)
from pkpdsim.pbpk.organ import Boundary, GutComponent, OrganComponent, OrganContext
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.parts import parts_for
from pkpdsim.pbpk.pd import PDComponent, ports
from pkpdsim.pbpk.physiology import physiology_objects, reference_tissue_volumes
from pkpdsim.pbpk.placement import (
    FormulaParameter,
    SpeciesStart,
    formula_assignment,
    formula_ids,
    formula_problems,
    inhibition_parameters,
    kp_id,
    species_start,
    start_formula,
    value_parameter,
)
from pkpdsim.pbpk.routes import ORAL, route_objects
from pkpdsim.pbpk.sinks import sink_objects
from pkpdsim.pbpk.substance import Kp, KpMethod, Substance
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.pbpk.values import Formula, Q, Value

__all__ = [
    "BuildResult",
    "build",
    "built_organs",
    "compose",
    "observer_problems",
    "organ_component",
    "organ_context",
]

_ONE = "1 dimensionless"
_ANNOTATED = (Compartment, Species, Parameter, Reaction)
_MODEL_FIELDS = (
    "name",
    "sboTerm",
    "annotations",
    "notes",
    "creators",
    "model_units",
    "conversionFactor",
    "units",
    "functions",
    "compartments",
    "species",
    "parameters",
    "assignments",
    "rules",
    "rate_rules",
    "algebraic_rules",
    "reactions",
    "events",
    "constraints",
    "ports",
)
"""The fields of the model of a PD component which its model definition takes."""
_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_SOURCE = {
    "competitive": "Segel1975",
    "noncompetitive": "Segel1975",
    "uncompetitive": "Segel1975",
    "mixed": "Segel1975",
    "mechanism": "Mayhew2000",
}
"""The source of the rate law of an inhibition by its kind."""


@dataclass(frozen=True)
class BuildResult:
    """The files of a built model.

    Attributes:
        comp_path: the comp model, the organs as submodels.
        flat_path: the flattened model, which is simulated.
    """

    comp_path: Path
    flat_path: Path


def organ_component(definition: PBPKDefinition, organ: Organ) -> OrganComponent:
    """Get the component of an organ: its variant, the standard organ by default.

    Args:
        definition: the definition.
        organ: the organ.

    Returns:
        The component.

    Raises:
        ValueError: if the variant of the organ is a name of no variant.
    """
    variant = definition.variants.organs.get(
        organ, variants.DEFAULTS.get(organ, "standard")
    )
    if not isinstance(variant, str):
        return variant
    return variants.lookup(organ, variant)


def _gut(definition: PBPKDefinition) -> GutComponent | None:
    """The intestine which replaces several organs, `None` for the standard organs.

    Raises:
        ValueError: if an organ which the intestine replaces has a variant of
            its own.
    """
    component = organ_component(definition, Organ.SMALL_INTESTINE)
    if not isinstance(component, GutComponent):
        return None
    for organ in component.replaces[1:]:
        if organ in definition.variants.organs:
            raise ValueError(
                f"Organ '{organ.value}' has a variant, but the intestine "
                f"'{component.name}' replaces it."
            )
    return component


def built_organs(definition: PBPKDefinition) -> list[Organ]:
    """Get the organs which are built as a submodel.

    Args:
        definition: the definition.

    Returns:
        Every organ but those the intestine replaces besides the small
        intestine (the large intestine).
    """
    gut = _gut(definition)
    replaced = set(gut.replaces[1:]) if gut is not None else set()
    return [organ for organ in Organ if organ not in replaced]


def organ_context(
    definition: PBPKDefinition,
    organ: Organ,
    top: Mapping[str, tuple[str, UnitType]] | None = None,
) -> OrganContext:
    """Get the context of an organ: every substance flows through every organ.

    Args:
        definition: the definition.
        organ: the organ.
        top: the parameters and compartments of the top model with their
            kind and unit, which a `Formula` of the organ may read.

    Returns:
        The substances with their parts in the organ, the chemistry and the
        inhibitions and inductions of the substances.
    """
    chemistry = tuple(definition.chemistry)
    substances = tuple(definition.substances)
    return OrganContext(
        organ=organ,
        substances=substances,
        parts={s.sid: parts_for(s, organ, chemistry) for s in substances},
        chemistry=chemistry,
        inhibitions=tuple((s.sid, i) for s in substances for i in s.inhibits),
        inductions=tuple((s.sid, i) for s in substances for i in s.induces),
        routes=tuple(definition.routes),
        top=dict(top or {}),
    )


def compose(definition: PBPKDefinition) -> Model:
    """Compose the comp model of a definition.

    Args:
        definition: the definition.

    Returns:
        The comp model: the top model with one submodel per organ and one per
        PD component, the observers of the definition after those of the
        library, the creators of the definition as model history.

    Raises:
        ValueError: if the definition is invalid, an organ variant is unknown,
            an organ reads an element the body does not have, a route is
            given twice, an oral route or an apical transport of the small
            intestine has no intestine with a lumen, a PD component names
            a port or an element of the body which does not exist, or an
            observer of the definition is invalid in the top model
            (`observer_problems`).
    """
    definition.validate(symbols=False)
    top, model_definitions, submodels, replaced = _body(definition)
    problems = _observer_problems(definition, top, model_definitions, submodels)
    if problems:
        raise ValueError(
            f"Invalid PBPK definition '{definition.sid}':\n- " + "\n- ".join(problems)
        )
    top.add(definition_observer_objects(definition.observers))
    top.add(_delta_objects(definition, top))

    return Model(
        sid=definition.sid,
        name=definition.name or definition.sid,
        annotations=model_annotations(definition.substances),
        notes=definition.notes,
        creators=list(definition.creators),
        units=U,
        model_units=MODEL_UNITS,
        objects=top.objects,
        model_definitions=model_definitions,
        submodels=submodels,
        replaced_elements=replaced,
    )


def _delta_objects(definition: PBPKDefinition, top: "_Top") -> list[Sbase]:
    """The start values and the changes of the observers of `PBPKDefinition.deltas`.

    For every observer `<id>` (a parameter of the top model with an
    assignment rule: an observer of the library or of the definition, an
    output of an organ) the constant `<id>_0` with the initial assignment
    `<id>`, its value at the start, and the observer `<id>_delta = <id> -
    <id>_0`. The start value is evaluated after the changes of the
    individual before the simulation, so the change is 0 at the start of
    every simulation.

    Raises:
        ValueError: naming a delta which is no observer of the top model and
            an id `<id>_0` or `<id>_delta` which the model has.
    """
    if not definition.deltas:
        return []
    ruled = {str(o.variable) for o in top.objects if isinstance(o, AssignmentRule)}
    units = {str(o.sid): o.unit for o in top.objects if isinstance(o, Parameter)}
    known = {
        str(o.sid)
        for o in top.objects
        if isinstance(o, Parameter | Compartment | Species | Reaction)
    }
    unknown = [d for d in definition.deltas if d not in ruled or d not in units]
    taken = [
        x
        for d in definition.deltas
        for x in (ids.start_value(d), ids.delta(d))
        if x in known
    ]
    problems = []
    if unknown:
        problems.append(
            f"the deltas {unknown} are no observers of the top model (an observer "
            "of the library or the definition, an output of an organ)"
        )
    if taken:
        problems.append(f"the ids {taken} of the deltas are ids of the model")
    if problems:
        raise ValueError(
            f"Invalid PBPK definition '{definition.sid}': " + "; ".join(problems) + "."
        )
    objects: list[Sbase] = []
    for d in definition.deltas:
        start, change, unit = ids.start_value(d), ids.delta(d), units[d]
        objects += [
            annotate(
                Parameter(
                    start,
                    None,
                    unit,
                    constant=True,
                    name=f"{d} at the start",
                    sboTerm=_SBO_PARAMETER,
                    notes=f"The value of `{d}` at the start of the simulation, after "
                    "the changes before it (`PBPKDefinition.deltas`).\n\n"
                    f"`{start} = {d}` (initial assignment)",
                ),
                "observer",
            ),
            initial_assignment(start, d, unit, f"{d} at the start"),
            annotate(
                Parameter(
                    change,
                    None,
                    unit,
                    constant=False,
                    name=f"change of {d} over its start",
                    sboTerm=_SBO_PARAMETER,
                    notes=f"The change of `{d}` over its value at the start "
                    "(`PBPKDefinition.deltas`).\n\n"
                    f"`{change} = {d} - {start}`",
                ),
                "observer",
            ),
            AssignmentRule(change, f"{d} - {start}", unit, name=f"change of {d}"),
        ]
    return objects


def observer_problems(definition: PBPKDefinition) -> list[str]:
    """Check the observers of a definition against the top model.

    The model is composed without the observers of the definition; an
    observer id must be no id of the body (nor of a submodel or model
    definition), every symbol of its formula a value of the top model (a
    parameter, species, compartment or reaction, never a submodel or model
    definition) or another observer, the
    observers must not read each other in a cycle, and every observer must
    read a state (`_observer_problems`).

    Args:
        definition: the definition, otherwise valid.

    Returns:
        The problems, naming the observers and the ids.
    """
    if not definition.observers:
        return []
    top, model_definitions, submodels, _ = _body(definition)
    return _observer_problems(definition, top, model_definitions, submodels)


def _observer_problems(
    definition: PBPKDefinition,
    top: "_Top",
    model_definitions: list[ModelDefinition],
    submodels: list[Submodel],
) -> list[str]:
    """The problems of the observers of a definition in its composed body.

    A state is what changes during a simulation: a species, the variable of
    a rule (the observers of the library, the outputs of the organs and PD
    components, the quantities downstream of a dynamic parameter), a
    variable parameter (a dynamic one), a reaction and `time`. An observer
    which reaches none of them, directly or through other observers, reads
    constants only: its value is fixed before the simulation, and an
    assignment rule of constants slows the compile of the model (R53), so it
    is rejected.
    """
    # the ids an observer may read: the values of the top model
    values: set[str] = {
        str(o.sid)
        for o in top.objects
        if isinstance(o, Parameter | Species | Compartment | Reaction)
    }
    # the ids an observer must not take: the values, the submodels and the
    # model definitions (no values, a formula cannot read them)
    body = (
        values
        | {str(sm.sid) for sm in submodels}
        | {str(md.sid) for md in model_definitions}
    )
    states: set[str] = {TIME}
    for o in top.objects:
        if isinstance(o, Species | Reaction):
            states.add(str(o.sid))
        elif isinstance(o, AssignmentRule | RateRule):
            states.add(str(o.variable))
        elif isinstance(o, Parameter | Compartment) and not o.constant:
            states.add(str(o.sid))
    observers = {o.sid: formula_symbols(o.formula) for o in definition.observers}
    problems: list[str] = []
    for sid, symbols in observers.items():
        if sid in body:
            problems.append(f"observer id '{sid}' is an id of the body")
        unknown = sorted(symbols - values - set(observers) - {TIME})
        if unknown:
            problems.append(
                f"observer '{sid}' reads {unknown}, which are no values of the top "
                "model (the pools, sinks, observers, covariates and outputs of "
                "the organs; an id of an organ is not visible in the top model)"
            )
    if problems:
        return problems
    reads = {sid: symbols & set(observers) for sid, symbols in observers.items()}
    remaining = dict(reads)
    while True:
        done = [sid for sid, rs in remaining.items() if not rs & set(remaining)]
        if not done:
            break
        for sid in done:
            del remaining[sid]
    if remaining:
        return [f"observers {sorted(remaining)} read each other in a cycle"]
    stateful = set(states)
    changed = True
    while changed:
        changed = False
        for sid, symbols in observers.items():
            if sid not in stateful and symbols & stateful:
                stateful.add(sid)
                changed = True
    return [
        f"observer '{sid}' reads constants only ({sorted(symbols)}): its value is "
        "fixed before the simulation, and an assignment rule of constants slows "
        "the compile of the model (R53); read a state (a species, an observer, "
        "an output of an organ, a dynamic parameter, time) or compute it from "
        "the results"
        for sid, symbols in observers.items()
        if sid not in stateful
    ]


def _body(
    definition: PBPKDefinition,
) -> tuple["_Top", list[ModelDefinition], list[Submodel], list[ReplacedElement]]:
    """Compose the body of a definition: everything but its own observers.

    Returns:
        The top model, the model definitions and submodels of the organs and
        the PD components, and the replaced elements.
    """
    gut = _gut(definition)
    oral = sorted({r.substance for r in definition.routes if r.kind in ORAL})
    apical = sorted(
        p.sid
        for p in definition.chemistry
        if isinstance(p, Transport)
        and p.organ is Organ.SMALL_INTESTINE
        and p.membrane is Membrane.APICAL
    )
    if gut is None and (oral or apical):
        raise ValueError(
            f"Oral routes of {oral} and apical transport {apical} of the small "
            "intestine need an intestine with a lumen (variant 'segmented' or "
            "'single_lumen'), the small intestine is a standard organ."
        )
    top = _Top(definition)
    top.parameters()
    for s in definition.substances:
        top.substance(s)
    for pool in Pool:
        top.pool(pool)
    top.add(sink_objects(definition, lumen=gut is not None))
    if gut is not None:
        top.add(gut.body_objects(definition))
    top.routes()

    model_definitions: list[ModelDefinition] = []
    submodels: list[Submodel] = []
    replaced: list[ReplacedElement] = []
    organ_models: list[ModelDefinition] = []
    organs = built_organs(definition)
    reads = _organ_reads(definition, organs)
    exposed: dict[str, str] = {}
    started: set[str] = set(top.started)
    for organ in organs:
        inputs = top.inputs()
        md, boundary = organ_component(definition, organ).build(
            organ_context(definition, organ, inputs)
        )
        if any(md.sid == other.sid for other in model_definitions):
            raise ValueError(f"Two organs build the model definition '{md.sid}'.")
        boundary, applied = _apply_starts(
            md, boundary, organ, definition.starts, inputs
        )
        started |= applied
        species_outputs = _expose_species(md, organ, reads.get(organ, {}))
        if species_outputs:
            boundary = Boundary(
                boundary.replacements, {**boundary.outputs, **species_outputs}
            )
            exposed |= {
                ids.flat(organ, local): ids.organ_species(organ, local)
                for local in reads[organ]
            }
        _scope_port_metaids(md)
        model_definitions.append(md)
        organ_models.append(md)
        submodel = ids.submodel(organ)
        if submodel in top.known:
            raise ValueError(
                f"The submodel id '{submodel}' of organ '{organ.value}' is an id "
                "of the top model."
            )
        submodels.append(Submodel(submodel, modelRef=md.sid, name=organ.name.lower()))
        replaced += top.connect(
            organ, submodel, boundary, species=species_outputs.values()
        )
    unknown = sorted(set(definition.starts) - started)
    if unknown:
        raise ValueError(
            f"Invalid PBPK definition '{definition.sid}': the starts of {unknown} "
            "name no species of the model (a species of an organ by its flat id "
            "`<ORGAN>__<species>`, a species of the top model by its id)."
        )
    top.add(observer_objects(definition))
    outputs: list[str] = []
    for pd in definition.pd:
        coupled = replace(
            pd,
            inputs={
                local: exposed.get(body, body) for local, body in pd.inputs.items()
            },
        )
        md, submodel, pd_replaced = top.couple(
            coupled, (str(sm.sid) for sm in submodels)
        )
        if any(md.sid == other.sid for other in model_definitions):
            raise ValueError(
                f"The PD component '{md.sid}' has the id of another model definition."
            )
        model_definitions.append(md)
        submodels.append(submodel)
        replaced += pd_replaced
        outputs += pd.outputs.values()
    problems = _parameter_problems(
        definition, top, model_definitions, organ_models, submodels, replaced
    )
    if problems:
        raise ValueError(
            f"Invalid PBPK definition '{definition.sid}':\n- " + "\n- ".join(problems)
        )
    top.formulas()
    _annotate_inputs(top.objects, model_definitions, submodels, replaced)
    _dynamic(
        definition,
        (*definition.dynamic, *(o for o in outputs if o not in definition.dynamic)),
        top.objects,
        model_definitions,
        submodels,
        replaced,
    )
    return top, model_definitions, submodels, replaced


def _start_number(q: Q) -> float:
    """The number of a start which is no formula, a plain float or a `Value` (mM)."""
    return q.value if isinstance(q, Value) else float(q)  # ty: ignore[invalid-argument-type]


def _apply_starts(
    md: ModelDefinition,
    boundary: Boundary,
    organ: Organ,
    starts: Mapping[str, Q],
    top: Mapping[str, tuple[str, UnitType]],
) -> tuple[Boundary, set[str]]:
    """Give the species of an organ the starts of the definition.

    A start (`PBPKDefinition.starts`) names a species of the organ by its
    flat id `<ORGAN>__<species>`. It replaces the initial value of the
    species: a number is its initial concentration, a `Formula` its initial
    assignment (replacing the start at the equilibrium with the baseline).
    Every id of a formula is a constant of the top model, which the organ
    reads through an input port under the same id, as for a `Formula` of the
    organ: an input the organ has under that id is used, else the input is
    added to the model definition and the boundary.

    Args:
        md: the model definition of the organ, changed in place.
        boundary: the boundary of the organ.
        organ: the organ.
        starts: the starts of the definition by flat id.
        top: the elements of the top model, id -> (kind, unit).

    Returns:
        The boundary with the inputs the formulas add and the flat ids of the
        starts given to a species of the organ; a start of the prefix of the
        organ which names no species is left out, `compose` names every such
        start of all organs at once.

    Raises:
        ValueError: naming a formula which reads an id which is no constant
            of the top model or an element of the organ which is no input
            under that id.
    """
    prefix = f"{ids.submodel(organ)}__"
    mine = {
        sid.removeprefix(prefix): q
        for sid, q in starts.items()
        if sid.startswith(prefix)
    }
    species = {str(sp.sid): sp for sp in md.species or []}
    mine = {local: q for local, q in mine.items() if local in species}
    if not mine:
        return boundary, set()
    replacements = dict(boundary.replacements)
    elements = {str(e.sid) for e in [*(md.parameters or []), *(md.compartments or [])]}
    problems: list[str] = []
    assignments = list(md.assignments or [])
    for local, q in sorted(mine.items()):
        sp = species[local]
        assignments = [
            a
            for a in assignments
            if not (isinstance(a, InitialAssignment) and str(a.symbol) == local)
        ]
        if not isinstance(q, Formula):
            sp.initialConcentration = _start_number(q)
            continue
        where = f"The start of '{prefix}{local}' (`{q.expression}`)"
        problems += formula_problems(where, q.expression, top)
        for symbol in sorted(formula_ids(q.expression)):
            kind, unit = top.get(symbol, ("", None))
            if kind not in ("parameter", "compartment"):
                continue
            if symbol in elements or symbol in species:
                if replacements.get(symbol) != symbol:
                    problems.append(
                        f"{where} reads '{symbol}', a local element of the organ"
                    )
                continue
            port = Port(f"{symbol}{PORT_SUFFIX}", portType=PortType.INPUT_PORT)
            text = f"{symbol} (start of {local})"
            if kind == "compartment":
                md.compartments = [
                    *(md.compartments or []),
                    Compartment(
                        symbol, 1.0, U.litre, constant=True, name=text, port=port
                    ),
                ]
            else:
                md.parameters = [
                    *(md.parameters or []),
                    Parameter(
                        symbol,
                        1.0,
                        unit,
                        constant=True,
                        name=text,
                        sboTerm=_SBO_PARAMETER,
                        port=port,
                    ),
                ]
            elements.add(symbol)
            replacements[symbol] = symbol
        sp.initialConcentration = None
        assignments.append(species_start(local, q.expression, str(sp.name)))
    if problems:
        raise ValueError("; ".join(problems) + ".")
    md.assignments = assignments
    return Boundary(replacements, boundary.outputs), {f"{prefix}{x}" for x in mine}


def _organ_reads(
    definition: PBPKDefinition, organs: Iterable[Organ]
) -> dict[Organ, dict[str, list[str]]]:
    """The species of the organs which the PD components read.

    An input of a PD component which names an element of a built organ by its
    flat id (`LI__Cic_api`, the prefix of the submodel of the organ) reads a
    species of the organ, which the organ exposes (`_expose_species`).

    Returns:
        By organ, the local ids of the species with the descriptions of the
        inputs which read them (`'<component>' input '<port>'`).
    """
    prefixes = {f"{ids.submodel(organ)}__": organ for organ in organs}
    reads: dict[Organ, dict[str, list[str]]] = {}
    for pd in definition.pd:
        for local, body in pd.inputs.items():
            for prefix, organ in prefixes.items():
                if body.startswith(prefix):
                    reads.setdefault(organ, {}).setdefault(
                        body.removeprefix(prefix), []
                    ).append(f"PD component '{pd.model.sid}': input '{local}'")
    return reads


def _expose_species(
    md: ModelDefinition, organ: Organ, reads: Mapping[str, list[str]]
) -> dict[str, str]:
    """Expose species of an organ which PD components read as outputs of the organ.

    Every read species gets a parameter of the organ with its assignment rule
    and an output port, `C<o>_<part>_<s>` (`ids.organ_species`), which
    replaces the parameter of the same id in the top model (`_Top.connect`),
    which the component reads. The species keeps its id: the flat model has
    the species `LI__Cic_<s>` and the parameter `Cli_ic_<s>` of its value.

    Args:
        md: the model definition of the organ, to which the outputs are added.
        organ: the organ.
        reads: the local ids of the read species with the inputs which read
            them, see `_organ_reads`.

    Raises:
        ValueError: naming an input which reads no species of the organ, an
            amount species or a species whose id is no concentration `C...`.

    Returns:
        The outputs, local id -> id of the top model (the same id).
    """
    species = {str(sp.sid): sp for sp in md.species or []}

    def problem(local: str) -> str | None:
        """Why a local id cannot be read as a concentration, `None` if it can."""
        flat = ids.flat(organ, local)
        sp = species.get(local)
        if sp is None:
            return (
                f"'{flat}', which is no species of the organ '{organ.value}' (a PD "
                "component reads a concentration of an organ by the flat id of "
                "its species)"
            )
        if getattr(sp, "hasOnlySubstanceUnits", False):
            return (
                f"'{flat}', an amount (hasOnlySubstanceUnits) of the organ "
                f"'{organ.value}': a PD component reads a concentration (mM)"
            )
        if not local.startswith("C"):
            return (
                f"'{flat}', a species of the organ '{organ.value}' whose id is no "
                "concentration of the id scheme, `C<part>_<s>`: the organ exposes "
                "it as `C<o>_<part>_<s>` (ids.organ_species)"
            )
        return None

    problems = [
        f"{what} reads {why}"
        for local, inputs in sorted(reads.items())
        if (why := problem(local)) is not None
        for what in inputs
    ]
    if problems:
        raise ValueError("; ".join(problems) + ".")
    outputs: dict[str, str] = {}
    for local in sorted(reads):
        top_id = ids.organ_species(organ, local)
        md.parameters = [
            *(md.parameters or []),
            annotate(
                Parameter(
                    top_id,
                    None,
                    U.mM,
                    constant=False,
                    name=f"concentration {local} of the organ {organ.name.lower()} "
                    "[mM]",
                    sboTerm=_SBO_PARAMETER,
                    notes=f"The species `{local}` of the organ, an output which a "
                    f"PD component reads.\n\n`{top_id} = {local}`",
                    port=Port(f"{top_id}{PORT_SUFFIX}", portType=PortType.OUTPUT_PORT),
                ),
                "organ_concentration",
                organ=organ,
            ),
        ]
        md.rules = [*(md.rules or []), AssignmentRule(top_id, local, U.mM)]
        outputs[top_id] = top_id
    return outputs


def _formula_cycles(objects: Iterable[Sbase]) -> list[str]:
    """The formulas of the top model which read each other in a cycle.

    A cycle of initial assignments has no value. A formula which reads no
    remaining formula, or which no remaining formula reads, is no part of a
    cycle and is removed until none is left to remove; the formulas which
    only read a cycle are not named.
    """
    reads = {
        str(o.sid): formula_ids(o.formula.expression)
        for o in objects
        if isinstance(o, FormulaParameter)
    }
    remaining = {sid: symbols & set(reads) for sid, symbols in reads.items()}
    while True:
        read = {x for rs in remaining.values() for x in rs}
        done = [
            sid
            for sid, rs in remaining.items()
            if not rs & set(remaining) or sid not in read
        ]
        if not done:
            break
        for sid in done:
            del remaining[sid]
    if remaining:
        return [f"the formulas of {sorted(remaining)} read each other in a cycle"]
    return []


def _parameter_problems(
    definition: PBPKDefinition,
    top: "_Top",
    model_definitions: Iterable[ModelDefinition],
    organ_models: Iterable[ModelDefinition],
    submodels: Iterable[Submodel],
    replaced: Iterable[ReplacedElement],
) -> list[str]:
    """The parameters of a definition which take an id of the body or of the definition.

    An id of the top model (a parameter of a substance, a pool, a sink, an
    observer of the library, the submodel or the model definition of an
    organ or a PD component) would be defined twice, and so would the id of
    an observer of the definition (added after the body); a local id of an
    organ which is no input under the same id (`Vti`, `<process>_Vmax`, the
    liver's `fimp_<protein>`) is what a custom reaction of the organ reads
    instead of the parameter.

    Returns:
        The problems, naming the parameter and what has its id.
    """
    sids = {p.sid for p in definition.parameters}
    if not sids:
        return []
    local = [p for o in top.objects if isinstance(o, Reaction) for p in o.pars or []]
    counts: dict[str, int] = {}
    for o in [*top.objects, *local]:
        if isinstance(o, Parameter | Compartment | Species | Reaction):
            counts[str(o.sid)] = counts.get(str(o.sid), 0) + 1
    for x in [*submodels, *model_definitions]:
        counts[str(x.sid)] = counts.get(str(x.sid), 0) + 1
    observers = {o.sid for o in definition.observers}
    inputs = {
        (str(r.submodelRef), str(r.portRef).removesuffix(PORT_SUFFIX))
        for r in replaced
        if str(r.portRef).removesuffix(PORT_SUFFIX) == str(r.elementRef)
    }
    submodel_of = {str(sm.modelRef): str(sm.sid) for sm in submodels}
    organs: dict[str, list[str]] = {}
    for md in organ_models:
        sm = submodel_of[str(md.sid)]
        for name in ("compartments", "species", "parameters", "reactions"):
            for element in getattr(md, name) or []:
                sid = str(element.sid)
                if sid in sids and (sm, sid) not in inputs:
                    organs.setdefault(sid, []).append(sm)
        for r in md.reactions or []:
            for par in r.pars or []:
                if str(par.sid) in sids:
                    organs.setdefault(str(par.sid), []).append(sm)
    problems: list[str] = []
    for p in definition.parameters:
        if p.sid in observers:
            problems.append(
                f"parameter '{p.sid}' is the id of an observer of the definition"
            )
        if counts.get(p.sid, 0) > 1:
            problems.append(
                f"parameter '{p.sid}' is an id of the body (the top model of the "
                "library, a submodel or a model definition of an organ or a PD "
                "component)"
            )
        if p.sid in organs:
            problems.append(
                f"parameter '{p.sid}' is a local id of the organ submodels "
                f"{sorted(set(organs[p.sid]))} of the library, which a custom "
                "reaction there would read instead"
            )
    return problems


_TOKEN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
"""An identifier in a formula."""


def _symbols(formula: object) -> set[str]:
    """The identifiers of a formula (ids, units and functions alike)."""
    return set(_TOKEN.findall(str(formula)))


def _to_rules(
    assignments: list[InitialAssignment],
    elements: Iterable[Parameter | Compartment],
    dynamic: set[str],
) -> list[AssignmentRule]:
    """Turn the initial assignments which depend on dynamic ids into assignment rules.

    `dynamic` is extended by the converted ids; the converted assignments are
    removed from `assignments` and their elements become variable. The start
    of a species of an endogenous substance (`SpeciesStart`,
    `Substance.baseline`) stays an initial assignment: the species starts at
    the value of the dynamic quantities at the start, its reactions change it
    afterwards. Every other initial assignment, also of a species of a PD
    component, is converted as before.

    Args:
        assignments: the initial assignments.
        elements: the parameters and compartments.
        dynamic: the dynamic ids.

    Returns:
        The assignment rules.
    """
    rules: list[AssignmentRule] = []
    changed = True
    while changed:
        changed = False
        for ia in list(assignments):
            if (
                ia.symbol in dynamic
                or isinstance(ia, SpeciesStart)
                or not _symbols(ia.value) & dynamic
            ):
                continue
            assignments.remove(ia)
            rules.append(AssignmentRule(ia.symbol, ia.value, ia.unit, name=ia.name))
            dynamic.add(ia.symbol)
            changed = True
    for element in elements:
        if element.sid in dynamic:
            element.constant = False
    return rules


def _annotate_inputs(
    objects: list[Sbase],
    model_definitions: list[ModelDefinition],
    submodels: list[Submodel],
    replaced: list[ReplacedElement],
) -> None:
    """Give every element of a component which is linked to the body its annotation.

    An input of an organ is the element of the body (`V<o>`, `fu_<s>`,
    `Cve_<s>`), replaced by it in the comp model; the comp model shows the
    input with the SBO term and the annotations of the element of the body,
    so that every element of the comp model carries them as the flat model
    does. An output (of an organ `Cout_<s>`, of a PD component the parameter
    of the body it replaces) without annotations of its own takes those of
    the element of the body, which the flat model shows under its id.
    """
    by_sid = {str(o.sid): o for o in objects if isinstance(o, _ANNOTATED)}
    by_submodel = {str(sm.sid): str(sm.modelRef) for sm in submodels}
    by_md = {str(md.sid): md for md in model_definitions}
    links = [(r.elementRef, r.submodelRef, r.portRef) for r in replaced] + [
        (o.replacedBy.elementRef, o.replacedBy.submodelRef, o.replacedBy.portRef)
        for o in by_sid.values()
        if isinstance(o, Parameter) and o.replacedBy is not None
    ]
    for element_ref, submodel_ref, port_ref in links:
        top = by_sid.get(str(element_ref))
        md = by_md.get(by_submodel.get(str(submodel_ref), ""))
        if top is None or md is None:
            continue
        local = str(port_ref).removesuffix(PORT_SUFFIX)
        for name in ("compartments", "species", "parameters"):
            for element in getattr(md, name) or []:
                if str(element.sid) == local and not element.annotations:
                    element.sboTerm = top.sboTerm
                    element.annotations = list(top.annotations or [])
                    element.metaId = element.metaId or f"meta_{md.sid}__{local}"


def _dynamic(
    definition: PBPKDefinition,
    names: Iterable[str],
    objects: list[Sbase],
    model_definitions: list[ModelDefinition],
    submodels: list[Submodel],
    replaced: list[ReplacedElement],
) -> None:
    """Keep the quantities downstream of the dynamic parameters as assignment rules.

    The parameters `names` (`definition.dynamic` and the parameters which a
    PD component replaces) change during a simulation: every
    initial assignment of the top model which reads one of them (directly or
    through another such quantity) becomes an assignment rule, and so does
    every initial assignment of an organ which reads such an element through
    its port; these elements and the dynamic parameters are variable.
    Everything else stays an initial assignment (a change before the
    simulation is followed, a change during it is not).

    Raises:
        ValueError: if a dynamic name is no parameter of the top model.
    """
    dynamic = set(names)
    if not dynamic:
        return
    parameters = {str(o.sid) for o in objects if isinstance(o, Parameter)}
    unknown = sorted(dynamic - parameters)
    if unknown:
        raise ValueError(
            f"{definition.sid}: dynamic {unknown} are no parameters of the top model."
        )
    assignments = [o for o in objects if isinstance(o, InitialAssignment)]
    elements = [o for o in objects if isinstance(o, Parameter | Compartment)]
    rules = _to_rules(assignments, elements, dynamic)
    objects[:] = [
        o for o in objects if not isinstance(o, InitialAssignment) or o in assignments
    ]
    objects += rules
    by_submodel = {str(sm.sid): str(sm.modelRef) for sm in submodels}
    by_sid = {str(md.sid): md for md in model_definitions}
    seeds: dict[str, set[str]] = {}
    for r in replaced:
        if r.elementRef in dynamic:
            local = str(r.portRef).removesuffix(PORT_SUFFIX)
            seeds.setdefault(by_submodel[str(r.submodelRef)], set()).add(local)
    for md_sid, local in seeds.items():
        md = by_sid[md_sid]
        local_assignments = list(md.assignments or [])
        md.rules = [
            *(md.rules or []),
            *_to_rules(
                local_assignments,
                [*(md.parameters or []), *(md.compartments or [])],
                local,
            ),
        ]
        md.assignments = local_assignments


def build(definition: PBPKDefinition, output_dir: Path) -> BuildResult:
    """Build a definition: the comp model and the flat model.

    Args:
        definition: the definition.
        output_dir: the folder of the files `<sid>.xml` and `<sid>_flat.xml`.

    Returns:
        The paths of the comp and the flat model.

    Raises:
        ValueError: if the definition cannot be composed, or a model has SBML
            errors.
    """
    model = compose(definition)
    comp_path = output_dir / f"{definition.sid}.xml"
    flat_path = output_dir / f"{definition.sid}_flat.xml"
    # SBML L3V2, the level of the models of pkdb_models: the names of the rules
    # and initial assignments are written (L3V1 has no name on them)
    create_model(
        model,
        filepath=comp_path,
        sbml_level=3,
        sbml_version=2,
        validate=True,
        raise_on_error=True,
    )
    doc = flatten_sbml(comp_path, flat_path)
    validation = validate_doc(doc, title=flat_path.name)
    if validation.error_count:
        raise ValueError(
            f"The flat model '{flat_path}' has {validation.error_count} SBML errors."
        )
    return BuildResult(comp_path=comp_path, flat_path=flat_path)


class _Top:
    """The objects of the top model, collected while the definition is composed."""

    def __init__(self, definition: PBPKDefinition) -> None:
        self.definition = definition
        self.objects: list[Sbase] = [
            *physiology_objects(
                definition.sex, physiology_factors(), definition.population
            ),
            *impairment_objects(definition),
        ]
        self.known: set[str] = {
            str(o.sid) for o in self.objects if isinstance(o, Parameter)
        }
        #: the ids of `PBPKDefinition.starts` of species of the top model
        self.started: set[str] = set()

    def add(self, objects: Iterable[Sbase]) -> None:
        """Add objects, recording the ids of compartments, species and parameters."""
        for o in objects:
            if isinstance(o, Parameter | Compartment | Species):
                self.known.add(str(o.sid))
            self.objects.append(o)

    def parameters(self) -> None:
        """The parameters of the definition (`PBPKDefinition.parameters`).

        A number or a `Value` is a constant parameter with its provenance and
        uncertainty, a `Formula` a `FormulaParameter` whose initial
        assignment `formulas` adds.
        """
        self.add(
            annotate(
                value_parameter(
                    p.sid,
                    p.value,
                    p.unit,
                    p.name or p.sid,
                    _SBO_PARAMETER,
                    p.notes or "",
                ),
                "definition_parameter",
            )
            for p in self.definition.parameters
        )

    def inputs(self) -> dict[str, tuple[str, UnitType]]:
        """The elements of the top model with their kind and unit.

        The kinds: `parameter` and `compartment` for the constants, which a
        `Formula` may read; `species` and `variable` (an observer, an output
        of an organ, any element with a rule) for what changes during the
        simulation, which a formula would freeze at its initial value.
        """
        ruled = {
            str(o.variable)
            for o in self.objects
            if isinstance(o, AssignmentRule | RateRule)
        }
        local = [
            p for o in self.objects if isinstance(o, Reaction) for p in o.pars or []
        ]
        found: dict[str, tuple[str, UnitType]] = {}
        for o in [*self.objects, *local]:
            if isinstance(o, Species):
                found[str(o.sid)] = ("species", o.substanceUnits)
            elif isinstance(o, Parameter | Compartment):
                sid = str(o.sid)
                if not o.constant or sid in ruled:
                    found[sid] = ("variable", o.unit)
                else:
                    kind = "parameter" if isinstance(o, Parameter) else "compartment"
                    found[sid] = (kind, o.unit)
        return found

    def formulas(self) -> None:
        """Give every `Formula` of the top model its initial assignment.

        Raises:
            ValueError: if a formula reads an id which is no element of the
                top model, or a species or an observer, naming every such id.
        """
        top = self.inputs()
        local = [
            p for o in self.objects if isinstance(o, Reaction) for p in o.pars or []
        ]
        baselines = {
            ids.baseline(s.sid)
            for s in self.definition.substances
            if s.baseline is not None
        }
        downstream = self._pd_downstream(local) if baselines else set()
        problems: list[str] = []
        assignments: list[Sbase] = []
        for o in [*self.objects, *local]:
            if not isinstance(o, FormulaParameter):
                continue
            where = f"the formula of '{o.sid}' (`{o.formula.expression}`)"
            if str(o.sid) in baselines:
                problems += self._baseline_problems(
                    where, o.formula.expression, top, downstream
                )
            else:
                problems += formula_problems(where, o.formula.expression, top)
            assignments.append(formula_assignment(o))
        problems += _formula_cycles([*self.objects, *local])
        if problems:
            raise ValueError(
                f"Invalid PBPK definition '{self.definition.sid}':\n- "
                + "\n- ".join(problems)
            )
        self.objects += assignments

    def _pd_downstream(self, local: Iterable[Parameter]) -> set[str]:
        """The parameters replaced by a PD output and the quantities which read them.

        A PD output changes during the simulation and may read the species
        (an effect of the arm vein concentration); the initial assignments
        of the top model and the formulas which read it, directly or through
        each other, follow it (they become rules, `_dynamic`).
        """
        down = {o for pd in self.definition.pd for o in pd.outputs.values()}
        if not down:
            return down
        reads = {
            str(o.symbol): _symbols(o.value)
            for o in self.objects
            if isinstance(o, InitialAssignment)
        }
        reads |= {
            str(o.sid): formula_ids(o.formula.expression)
            for o in [*self.objects, *local]
            if isinstance(o, FormulaParameter)
        }
        changed = True
        while changed:
            changed = False
            for sid, symbols in reads.items():
                if sid not in down and symbols & down:
                    down.add(sid)
                    changed = True
        return down

    def _baseline_problems(
        self,
        where: str,
        expression: str,
        top: Mapping[str, tuple[str, UnitType]],
        downstream: set[str],
    ) -> list[str]:
        """Check the formula of a baseline: parameters, compartments and observers.

        The formula of `baseline_<s>` may read, besides the parameters and
        compartments of the top model, the observers of the definition (a
        production over a clearance): they are evaluated at the start. What
        it reads, directly or through the observers, must not change during
        the simulation (species, reactions, the observers of the library, the
        outputs of the organs, which read species, and the parameters which
        a PD output replaces with what reads them, `downstream`): the species
        start at the baseline, so it would read its own start.

        Returns:
            The problems: the unknown ids and the ids which change during the
            simulation, with the observer through which they are read.
        """
        observers = {
            o.sid: formula_symbols(o.formula) for o in self.definition.observers
        }
        symbols = formula_ids(expression)
        problems: list[str] = []
        unknown = sorted(x for x in symbols if x not in top and x not in observers)
        if unknown:
            problems.append(
                f"{where} reads ids which are unknown in the top model: {unknown}"
            )
        changing: list[str] = []
        seen: set[str] = set()
        stack: list[tuple[str, str | None]] = [(x, None) for x in sorted(symbols)]
        while stack:
            x, via = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            if x in observers:
                stack += [(y, x) for y in sorted(observers[x])]
                continue
            pd = x in downstream
            if pd or (x in top and top[x][0] not in ("parameter", "compartment")):
                label = f"{x} (a PD output or downstream of one)" if pd else x
                changing.append(
                    label if via is None else f"{label} (through observer {via})"
                )
        if changing:
            problems.append(
                f"{where} reads species, observers or PD outputs of the top model, "
                f"which change during the simulation: {sorted(changing)}; the "
                "species start at the baseline, which reads parameters and "
                "compartments of the top model and observers of the definition "
                "which read them"
            )
        return problems

    def _initial(
        self,
        pid: str,
        formula: str,
        unit: UnitDefinition,
        name: str,
        text: str,
        kind: str,
    ) -> None:
        """A constant parameter with its initial assignment, annotated as `kind`."""
        self.add(
            [
                annotate(
                    Parameter(
                        pid,
                        None,
                        unit,
                        constant=True,
                        name=name,
                        sboTerm=_SBO_PARAMETER,
                        notes=f"{text}\n\n`{pid} = {formula}`",
                    ),
                    kind,
                ),
                initial_assignment(pid, formula, unit, name),
            ]
        )

    def _compartment(
        self,
        cid: str,
        formula: str,
        name: str,
        text: str,
        kind: str,
        pool: Pool,
    ) -> None:
        """A constant compartment of a pool with its volume as initial assignment."""
        self.add(
            [
                annotate(
                    Compartment(
                        cid,
                        None,
                        U.litre,
                        constant=True,
                        name=name,
                        notes=f"{text}\n\n`{cid} = {formula}`",
                    ),
                    kind,
                    pool=pool,
                ),
                initial_assignment(cid, formula, U.litre, name),
            ]
        )

    def _species(
        self,
        sid: str,
        compartment: str,
        name: str,
        kind: str,
        s: Substance,
        start: str | None = None,
    ) -> None:
        """A species of a substance in the top model, initially zero or at `start`.

        `start` is the initial concentration as formula (the baseline of an
        endogenous substance, `start_formula`), the initial assignment of
        the species; a start of the definition (`PBPKDefinition.starts`)
        replaces it.

        Raises:
            ValueError: if the formula of a start of the definition reads an
                id which is no constant of the top model.
        """
        initial = 0.0 if start is None else None
        given = self.definition.starts.get(sid)
        if given is not None:
            self.started.add(sid)
            if isinstance(given, Formula):
                problems = formula_problems(
                    f"The start of '{sid}' (`{given.expression}`)",
                    given.expression,
                    self.inputs(),
                )
                if problems:
                    raise ValueError("; ".join(problems) + ".")
                initial, start = None, given.expression
            else:
                initial, start = float(_start_number(given)), None
        self.add(
            [
                annotate(
                    Species(
                        sid,
                        compartment,
                        initialConcentration=initial,
                        substanceUnit=U.mmole,
                        hasOnlySubstanceUnits=False,
                        name=name,
                    ),
                    kind,
                    substance=s,
                )
            ]
        )
        if start is not None:
            self.add([species_start(sid, start, name)])

    def _reaction(
        self, sid: str, equation: str, rate: str, name: str, pool: Pool
    ) -> None:
        """An exchange of a pool, the rate in mmole/min."""
        self.add(
            [
                annotate(
                    Reaction(
                        sid=sid,
                        equation=equation,
                        formula=(rate, U.mmole_per_min),
                        name=name,
                    ),
                    "exchange",
                    pool=pool,
                )
            ]
        )

    # --- substances ------------------------------------------------------

    def substance(self, s: Substance) -> None:
        """The parameters of a substance: Mr, fu, the distribution and the blood."""
        sid = s.sid
        self.add(
            [
                annotate(
                    value_parameter(
                        ids.mr(sid),
                        s.mr,
                        U.g_per_mole,
                        f"molecular weight {sid} [g/mole]",
                        SBO.MOLECULAR_MASS,
                    ),
                    "molecular_weight",
                ),
                annotate(
                    value_parameter(
                        ids.fu_ref(sid),
                        s.binding.fu if s.binding is not None else 1.0,
                        U.dimensionless,
                        f"reference fraction unbound {sid}",
                        _SBO_PARAMETER,
                        "Fraction unbound in plasma of the reference individual"
                        + (
                            f" (binding protein {s.binding.protein})."
                            if s.binding is not None
                            else "; the substance does not bind to plasma proteins."
                        ),
                    ),
                    "fu",
                ),
            ]
        )
        self._initial(
            ids.fu(sid),
            unbound_fraction(s),
            U.dimensionless,
            f"fraction unbound {sid}",
            unbound_fraction_notes(s),
            "fu",
        )
        if s.excretion.urine:
            self.add(
                [
                    annotate(
                        value_parameter(
                            ids.fraction_reabsorbed(sid),
                            s.excretion.fraction_reabsorbed,
                            U.dimensionless,
                            f"fraction reabsorbed {sid}",
                            _SBO_PARAMETER,
                            "Fraction of the glomerular filtrate which is reabsorbed "
                            "into the plasma of the kidney.",
                        ),
                        "fraction_reabsorbed",
                    )
                ]
            )
        if s.baseline is not None:
            self.add(
                [
                    annotate(
                        value_parameter(
                            ids.baseline(sid),
                            s.baseline,
                            U.mM,
                            f"baseline {sid} [mM]",
                            _SBO_PARAMETER,
                            "Plasma concentration of the endogenous substance at the "
                            "start: every plasma species starts at it, every other "
                            "species of the blood and the tissue at the equilibrium "
                            "of its exchange with plasma, the lumen and the sinks "
                            "at 0.",
                        ),
                        "baseline",
                    )
                ]
            )
        self._distribution(s)
        self._erythrocytes(s)
        self._interactions(s)

    def _interactions(self, s: Substance) -> None:
        """The constants of the inhibitions and inductions of a substance.

        Only of a protein which a process of the definition uses: the organs
        with such a process read them (`Ki`, `Kiu`, `kinact`, `Emax`,
        `EC50`), one value for every organ.
        """
        used = {
            p.protein
            for p in self.definition.chemistry
            if isinstance(p, Metabolism | Transport) and p.protein is not None
        }
        for inh in s.inhibits:
            if inh.protein not in used:
                continue
            for pid, field in inhibition_parameters(s.sid, inh):
                q = getattr(inh, field)
                assert q is not None
                if field == "kinact":
                    unit, sbo, text = (
                        U.per_min,
                        SBO.INACTIVATION_RATE_CONSTANT,
                        f"Maximal inactivation rate of {inh.protein} by {s.sid} "
                        "(mechanism-based inhibition).",
                    )
                else:
                    unit, sbo = U.mM, SBO.INHIBITORY_CONSTANT
                    text = {
                        "ki": f"Inhibition constant of {s.sid} on {inh.protein} "
                        f"({inh.kind}) on the unbound concentration"
                        + (
                            "; the KI of the inactivation, and the competitive "
                            "constant of the inhibitor while it binds."
                            if inh.kind == "mechanism"
                            else "."
                        ),
                        "kiu": f"Uncompetitive constant of the mixed inhibition of "
                        f"{inh.protein} by {s.sid} on the unbound concentration.",
                    }[field]
                self.add(
                    [
                        annotate(
                            value_parameter(
                                pid,
                                q,
                                unit,
                                f"{field} {inh.protein} {s.sid}",
                                sbo,
                                f"{text}\n\nSource of the rate law: "
                                f"{cite(_SOURCE[inh.kind])}",
                            ),
                            "kinact" if field == "kinact" else "ki",
                            protein=inh.protein,
                        )
                    ]
                )
        for ind in s.induces:
            if ind.protein not in used:
                continue
            self.add(
                [
                    annotate(
                        value_parameter(
                            ids.emax(ind.protein, s.sid),
                            ind.emax,
                            U.dimensionless,
                            f"Emax {ind.protein} {s.sid}",
                            _SBO_PARAMETER,
                            f"Maximal induction of {ind.protein} by {s.sid}, the "
                            "fold increase of its synthesis minus 1.",
                        ),
                        "emax",
                        protein=ind.protein,
                    ),
                    annotate(
                        value_parameter(
                            ids.ec50(ind.protein, s.sid),
                            ind.ec50,
                            U.mM,
                            f"EC50 {ind.protein} {s.sid}",
                            SBO.EC50,
                            f"Unbound concentration of {s.sid} of half-maximal "
                            f"induction of {ind.protein}.",
                        ),
                        "ec50",
                        protein=ind.protein,
                    ),
                ]
            )

    def _distribution(self, s: Substance) -> None:
        """The tissue to plasma coefficients, ftissue and the reference tissue volume."""
        dist = s.distribution
        if dist is None:
            return
        sid = s.sid
        objects = (
            self._kp_method(s, dist)
            if isinstance(dist, KpMethod)
            else self._kp_value(s, dist)
        )
        if dist.ftissue is not None:
            objects.append(
                annotate(
                    value_parameter(
                        ids.ftissue(sid),
                        dist.ftissue,
                        U.l_per_min,
                        f"whole-body exchange {sid} [l/min]",
                        _SBO_PARAMETER,
                        "Exchange between plasma and tissue of the whole body of the "
                        "reference individual; an organ has the share of its tissue "
                        "volume. It multiplies `fu*(Cpl - Cti/Kp)`, so its fitted "
                        "values differ from those of the hand-written models.",
                    ),
                    "permeability",
                )
            )
        shared = self._ftissue_organs(s)
        if shared:
            volumes = reference_tissue_volumes(self.definition.sex)
            total = sum(volumes[o] for o in shared)
            objects.append(
                annotate(
                    Parameter(
                        ids.tissue_volume_total_ref(sid),
                        total,
                        U.litre,
                        constant=True,
                        name=f"reference tissue volume of the organs {sid} enters [l]",
                        sboTerm=_SBO_PARAMETER,
                        notes="Tissue volume `V<o>*(1 - f_vas_<o>)` of the "
                        f"{self.definition.sex.value} reference individual summed "
                        "over the organs whose exchange is a share of "
                        f"`{ids.ftissue(sid)}`: "
                        + ", ".join(o.name.lower() for o in shared)
                        + ". A constant, so that the exchange scales with the "
                        "tissue volume of the individual.",
                    ),
                    "tissue_volume",
                )
            )
        self.add(objects)

    def _kp_value(self, s: Substance, dist: Kp) -> list[Sbase]:
        """The coefficient of a substance and the coefficients of single organs.

        With `interstitial` every organ the substance enters without an
        override has `Kp_<o>_<s> = fKp_<s> * share`, an initial assignment
        with the interstitial share of the organ as number
        (`distribution.interstitial_share`).
        """
        sid = s.sid
        objects: list[Sbase] = []
        if dist.value is not None:
            objects.append(
                annotate(
                    value_parameter(
                        ids.kp(sid),
                        dist.value,
                        U.dimensionless,
                        f"tissue to plasma coefficient {sid}",
                        _SBO_PARAMETER,
                        "Whole-tissue (interstitial and intracellular) to plasma "
                        "coefficient of the organs without an own value.",
                    ),
                    "kp",
                )
            )
        for organ, q in dist.overrides.items():
            objects.append(
                annotate(
                    value_parameter(
                        ids.kp_organ(organ, sid),
                        q,
                        U.dimensionless,
                        f"tissue to plasma coefficient {sid} {organ.name.lower()}",
                        _SBO_PARAMETER,
                        "Whole-tissue to plasma coefficient of the organ.",
                    ),
                    "kp_organ",
                    organ=organ,
                )
            )
        if dist.interstitial is not None:
            objects += self._kp_interstitial(s, dist)
        return objects

    def _kp_interstitial(self, s: Substance, dist: Kp) -> list[Sbase]:
        """The scale of the interstitial shares and the Kp of every other organ."""
        sid, fkp = s.sid, ids.fkp(s.sid)
        assert dist.interstitial is not None
        fu = distribution.number(s.binding.fu) if s.binding is not None else 1.0
        objects: list[Sbase] = [
            annotate(
                value_parameter(
                    fkp,
                    dist.interstitial,
                    U.dimensionless,
                    f"factor of the interstitial shares {sid}",
                    _SBO_PARAMETER,
                    "Factor on the interstitial share of every organ without an "
                    "own Kp, the extracellular distribution, so that a fit scales "
                    "all of them.",
                ),
                "kp_factor",
            )
        ]
        organs = [
            organ
            for organ in built_organs(self.definition)
            if (dist.organs is None or organ in dist.organs)
            and organ not in dist.overrides
        ]
        for organ in organs:
            pid = ids.kp_organ(organ, sid)
            share = distribution.interstitial_share(organ, fu)
            formula = f"{fkp} * {share!r} dimensionless"
            objects += [
                annotate(
                    Parameter(
                        pid,
                        None,
                        U.dimensionless,
                        constant=True,
                        name=f"tissue to plasma coefficient {sid} {organ.name.lower()}",
                        sboTerm=_SBO_PARAMETER,
                        notes="Whole-tissue to plasma coefficient of the organ, its "
                        "interstitial share times the factor: "
                        f"`f_int*K_int_pls(fu)/(1 - f_vas)` = {share:.6g} at fu "
                        f"{fu:.6g} with the fractions of the reference individuals, "
                        "the larger of both sexes; the substance reaches the "
                        f"interstitial fluid and not the cells.\n\n`{pid} = "
                        f"{formula}`",
                    ),
                    "kp_organ",
                    organ=organ,
                ),
                initial_assignment(
                    pid,
                    formula,
                    U.dimensionless,
                    f"tissue to plasma coefficient {sid} {organ.name.lower()}",
                ),
            ]
        return objects

    def _kp_method(self, s: Substance, method: KpMethod) -> list[Sbase]:
        """The factor of a Kp method and the coefficient of every organ it computes.

        Every `Kp_<o>_<s>` is an initial assignment `fKp_<s> * value` with the
        value of the method (`pkpdsim.pbpk.kp`) for the sex of the build.
        """
        sid, fkp = s.sid, ids.fkp(s.sid)
        objects: list[Sbase] = [
            annotate(
                value_parameter(
                    fkp,
                    method.factor,
                    U.dimensionless,
                    f"factor of the tissue to plasma coefficients {sid}",
                    _SBO_PARAMETER,
                    "Factor on every tissue to plasma coefficient of the Kp method "
                    f"({method.method}), so that a fit scales all of them.",
                ),
                "kp_factor",
            )
        ]
        for organ, v in kp_values(s, self.definition.sex).items():
            pid = ids.kp_organ(organ, sid)
            formula = f"{fkp} * {v.value!r} dimensionless"
            citations = "\n".join(f"- {cite(key)}" for key in kp_sources(organ))
            objects += [
                annotate(
                    Parameter(
                        pid,
                        None,
                        U.dimensionless,
                        constant=True,
                        name=f"tissue to plasma coefficient {sid} {organ.name.lower()}",
                        sboTerm=_SBO_PARAMETER,
                        notes=f"Whole-tissue to plasma coefficient of the organ, the "
                        f"Kp method times its factor. {v.notes}\n\n`{pid} = "
                        f"{formula}`\n\nSources:\n\n{citations}",
                    ),
                    "kp_organ",
                    organ=organ,
                ),
                initial_assignment(
                    pid,
                    formula,
                    U.dimensionless,
                    f"tissue to plasma coefficient {sid} {organ.name.lower()}",
                ),
            ]
        return objects

    def _ftissue_organs(self, s: Substance) -> list[Organ]:
        """The organs whose exchange of a substance is a share of `ftissue_<s>`.

        The merged and the full form of an organ the substance enters, unless
        the distribution gives the exchange of the organ (`ps`).
        """
        dist = s.distribution
        if dist is None:
            return []
        return [
            organ
            for organ in Organ
            if parts_for(s, organ, self.definition.chemistry).tissue
            in {"merged", "full"}
            and kp_id(s, organ) is not None
            and organ not in dist.ps
        ]

    def _erythrocytes(self, s: Substance) -> None:
        """The blood to plasma ratio and the parameters of the blood cells."""
        sid, ery = s.sid, s.erythrocytes
        bp, name = ids.bp(sid), f"blood to plasma ratio {sid}"
        if ery is None:
            self._initial(
                bp,
                f"{_ONE} - {ids.HCT}",
                U.dimensionless,
                name,
                "Blood to plasma ratio of a substance which does not enter the "
                "erythrocytes; an observer value, the flows carry plasma only.",
                "bp",
            )
            return
        if ery.bp is not None:
            self.add(
                [
                    annotate(
                        value_parameter(
                            bp,
                            ery.bp,
                            U.dimensionless,
                            name,
                            _SBO_PARAMETER,
                            "Blood to plasma ratio, erythrocytes in equilibrium.",
                        ),
                        "bp",
                    )
                ]
            )
            return
        assert ery.kbc is not None
        self.add(
            [
                annotate(
                    value_parameter(
                        ids.kbc(sid),
                        ery.kbc,
                        U.dimensionless,
                        f"blood cell to plasma coefficient {sid}",
                        _SBO_PARAMETER,
                    ),
                    "kbc",
                )
            ]
        )
        if ery.k_exchange is not None:
            self.add(
                [
                    annotate(
                        value_parameter(
                            ids.kbc_exchange(sid),
                            ery.k_exchange,
                            U.per_min,
                            f"exchange plasma and blood cells {sid} [1/min]",
                            _SBO_PARAMETER,
                            "Rate of the exchange between plasma and blood cells per "
                            "blood cell volume.",
                        ),
                        "exchange_rate",
                    )
                ]
            )
        self._initial(
            bp,
            f"{_ONE} - {ids.HCT} + {ids.HCT} * {ids.kbc(sid)}",
            U.dimensionless,
            name,
            "Blood to plasma ratio of the blood cell to plasma coefficient"
            + (
                " (at equilibrium of the kinetic blood cells)."
                if ery.mode == "kinetic"
                else "."
            ),
            "bp",
        )

    # --- pools -----------------------------------------------------------

    def pool(self, pool: Pool) -> None:
        """The compartments and species of a blood pool."""
        label = pool.name.lower()
        volume = ids.pool_volume(pool)
        plasma = ids.pool_compartment(pool, Part.PLASMA)
        cells = ids.pool_compartment(pool, Part.BLOOD_CELLS)
        self._compartment(
            plasma,
            f"{volume} * ({_ONE} - {ids.HCT})",
            f"{label} plasma",
            "Plasma of the blood pool.",
            "plasma",
            pool,
        )
        self._compartment(
            cells,
            f"{volume} * {ids.HCT}",
            f"{label} blood cells",
            "Blood cells of the blood pool.",
            "blood_cells",
            pool,
        )
        for s in self.definition.substances:
            ery = s.erythrocytes
            compartment = plasma
            if ery is not None and ery.mode == "equilibrium":
                compartment = ids.pool_capacity(pool, s.sid)
                self._compartment(
                    compartment,
                    f"{volume} * {ids.bp(s.sid)}",
                    f"plasma-equivalent capacity V*BP of {s.sid}, {label} pool",
                    "Not a physical volume: the plasma-equivalent capacity of the "
                    "blood of a substance with erythrocytes in equilibrium, the "
                    "pool volume times the blood to plasma ratio, so that the "
                    "plasma concentration in it gives the amount in the blood.",
                    "capacity",
                    pool,
                )
            cpl = ids.pool_species(pool, s.sid)
            self._species(
                cpl,
                compartment,
                f"{s.name} ({label} plasma)",
                "concentration_plasma",
                s,
                start_formula(s),
            )
            if ery is None or ery.mode != "kinetic":
                continue
            cbc = ids.pool_species(pool, s.sid, Part.BLOOD_CELLS)
            self._species(
                cbc,
                cells,
                f"{s.name} ({label} blood cells)",
                "concentration",
                s,
                start_formula(s, ids.kbc(s.sid)),
            )
            self._reaction(
                ids.pool_exchange_reaction(pool, s.sid),
                f"{cpl} <-> {cbc}",
                f"{ids.kbc_exchange(s.sid)} * {cells} * "
                f"({cpl} - {cbc} / {ids.kbc(s.sid)})",
                f"exchange {s.sid} plasma and blood cells, {label} pool",
                pool,
            )

    # --- routes ----------------------------------------------------------

    def routes(self) -> None:
        """The dose states of the routes.

        Raises:
            ValueError: if a route of a substance is given twice.
        """
        by_id = {s.sid: s for s in self.definition.substances}
        seen: set[tuple[str, str]] = set()
        for route in self.definition.routes:
            key = (route.kind, route.substance)
            if key in seen:
                raise ValueError(
                    f"Route '{route.kind}' of substance '{route.substance}' is "
                    "given twice."
                )
            seen.add(key)
            self.add(route_objects(route, by_id[route.substance]))

    # --- organs ----------------------------------------------------------

    def couple(
        self, pd: PDComponent, submodels: Iterable[str]
    ) -> tuple[ModelDefinition, Submodel, list[ReplacedElement]]:
        """The model definition, submodel and replacements of a PD component.

        Every input of the component is replaced by the element of the body;
        every parameter of the body which an output names is replaced by the
        element of the output (a `ReplacedBy` of the parameter). The target of
        an output is a plain parameter of a value: one which no initial
        assignment, rule or other component defines (a derived quantity
        such as `CO` is changed through the parameters it derives from). A
        parameter of `PBPKDefinition.dynamic` is such a parameter and may be
        the target; it is dynamic either way.

        Args:
            pd: the component.
            submodels: the ids of the submodels of the organs.

        Raises:
            ValueError: naming an input or output which is no port of the
                component, of the wrong type, an element of the body which
                does not exist (an input) or is no plain parameter of a value
                (an output), a component whose id is an id of the top model or
                of an organ submodel, or whose model units differ from those of
                the library.
        """
        model = copy.deepcopy(pd.model)
        sid = str(model.sid)
        if sid in {p.sid for p in self.definition.parameters}:
            raise ValueError(
                f"Invalid PBPK definition '{self.definition.sid}':\n- parameter "
                f"'{sid}' is an id of the PD component '{sid}' (its submodel and "
                "model definition)"
            )
        if sid in self.known:
            raise ValueError(
                f"The PD component '{sid}' has the id of an element of the body."
            )
        if sid in set(submodels):
            raise ValueError(
                f"The PD component '{sid}' has the id of the submodel of an organ."
            )
        _check_model_units(sid, model)
        if model.model_definitions or model.external_model_definitions:
            raise ValueError(
                f"The PD component '{sid}' has model definitions of its own; "
                "comp does not nest them."
            )
        declared = ports(model)
        parameters = {str(o.sid): o for o in self.objects if isinstance(o, Parameter)}
        derived = {
            str(o.symbol if isinstance(o, InitialAssignment) else o.variable)
            for o in self.objects
            if isinstance(o, InitialAssignment | AssignmentRule | RateRule)
        }
        problems: list[str] = []
        for what, mapping, wrong in (
            ("input", pd.inputs, PortType.OUTPUT_PORT),
            ("output", pd.outputs, PortType.INPUT_PORT),
        ):
            for local, body in mapping.items():
                if local not in declared:
                    problems.append(f"{what} '{local}' is no port of the component")
                elif declared[local][1] is wrong:
                    problems.append(
                        f"{what} '{local}' is declared as {wrong.value} of the "
                        "component"
                    )
                if what == "input" and body not in self.known:
                    problems.append(f"input '{local}' reads the unknown id '{body}'")
                if what == "output" and body not in parameters:
                    problems.append(
                        f"output '{local}' replaces '{body}', which is no "
                        "parameter of the body"
                    )
                elif what == "output" and (
                    body in derived or parameters[body].replacedBy is not None
                ):
                    problems.append(
                        f"output '{local}' replaces '{body}', which is derived "
                        "(an initial assignment, a rule or the output of another "
                        "component), not a parameter of a value"
                    )
        if problems:
            raise ValueError(f"PD component '{sid}': " + "; ".join(problems) + ".")
        md = ModelDefinition(
            sid=sid,
            **{f: getattr(model, f) for f in _MODEL_FIELDS if getattr(model, f)},
        )
        _scope_port_metaids(md)
        for local, body in pd.outputs.items():
            parameters[body].replacedBy = ReplacedBy(
                elementRef=body, submodelRef=sid, portRef=declared[local][0]
            )
        submodel = Submodel(sid, modelRef=sid, name=model.name or sid)
        replaced = [
            ReplacedElement(
                elementRef=body, submodelRef=sid, portRef=declared[local][0]
            )
            for local, body in pd.inputs.items()
        ]
        return md, submodel, replaced

    def connect(
        self,
        organ: Organ,
        submodel: str,
        boundary: Boundary,
        species: Iterable[str] = (),
    ) -> list[ReplacedElement]:
        """The replacements of the boundary of an organ.

        Every element of the body which the organ reads replaces the port of
        the organ; every output of the organ replaces a parameter of the top
        model, which is created here: the outflow concentration
        `Cout_<o>_<s>`, the urinary excretion rate `Rurine_<s>` of the
        kidney, a species of the organ which a PD component reads
        `C<o>_<part>_<s>` (`species`, see `_expose_species`).

        Raises:
            ValueError: if the organ reads an element the body does not have,
                or an exposed species has the id of an element of the body.
        """
        exposed = set(species)
        clash = sorted(exposed & self.known)
        if clash:
            raise ValueError(
                f"The species of organ '{organ.value}' which a PD component reads "
                f"have the ids of elements of the body: {clash}."
            )
        missing = sorted(
            t for t in boundary.replacements.values() if t not in self.known
        )
        if missing:
            raise ValueError(
                f"Organ '{organ.value}' reads elements which the body does not "
                f"have: {missing}."
            )
        cells = {ids.outflow_conc_bc(s.sid) for s in self.definition.substances}
        rates = {ids.urine_rate(s.sid) for s in self.definition.substances}
        for local, top in boundary.outputs.items():
            if top in rates:
                unit, name, kind = (
                    U.mmole_per_min,
                    f"urinary excretion rate {top} [mmole/min]",
                    "urine_rate",
                )
                text = "Rate of the urinary excretion, the sum of the fluxes of the "
            elif top in exposed:
                unit, name, kind = (
                    U.mM,
                    f"concentration {top} of an organ [mM]",
                    "organ_concentration",
                )
                text = "Concentration of a species, read by a PD component, of the "
            else:
                unit, name, kind = (
                    U.mM,
                    f"outflow concentration {top} [mM]",
                    "outflow_blood_cells" if local in cells else "outflow",
                )
                text = "Concentration (plasma or kinetic blood cells) leaving the "
            self.add(
                [
                    annotate(
                        Parameter(
                            top,
                            0.0,
                            unit,
                            constant=False,
                            name=name,
                            sboTerm=_SBO_PARAMETER,
                            notes=f"{text}organ {organ.name.lower()}, replaced by "
                            "the output of the organ; the value is a placeholder.",
                            replacedBy=ReplacedBy(
                                elementRef=top,
                                submodelRef=submodel,
                                portRef=f"{local}{PORT_SUFFIX}",
                            ),
                        ),
                        kind,
                        organ=organ,
                    )
                ]
            )
        return [
            ReplacedElement(
                elementRef=top,
                submodelRef=submodel,
                portRef=f"{local}{PORT_SUFFIX}",
            )
            for local, top in boundary.replacements.items()
        ]


_UNIT_FIELDS = ("time", "extent", "substance", "length", "area", "volume")


def _unit_id(unit: object) -> str | None:
    """The id of a model unit (a unit definition or its id)."""
    if unit is None:
        return None
    return str(getattr(unit, "sid", unit))


def _check_model_units(sid: str, model: Model) -> None:
    """Check that a PD component has the model units of the library.

    The flat model has one set of model units; a component in other units
    (hours) would be integrated in the units of the library without a
    conversion factor.

    Raises:
        ValueError: naming the model units which differ, or if the component
            has none.
    """
    units = model.model_units
    if units is None:
        raise ValueError(
            f"The PD component '{sid}' has no model units, the library has "
            "`pkpdsim.pbpk.MODEL_UNITS`."
        )
    differ = [
        f"{field} {_unit_id(getattr(units, field, None))} (library "
        f"{_unit_id(getattr(MODEL_UNITS, field))})"
        for field in _UNIT_FIELDS
        if _unit_id(getattr(units, field, None))
        != _unit_id(getattr(MODEL_UNITS, field))
    ]
    if differ:
        raise ValueError(
            f"The PD component '{sid}' has other model units than the library: "
            + ", ".join(differ)
            + "."
        )


def _scope_port_metaids(md: ModelDefinition) -> None:
    """Give the ports of a model definition metaids unique in the document.

    The port of an element is `<id>_port`, and sbmlutils gives a port the
    metaid of its id; a metaid is unique in the whole document, so the ports
    of two organs with the same local ids (`Vpl_port`) would collide. The
    metaid is prefixed with the id of the model definition, the port ids, to
    which the replacements refer, stay. The same holds for the metaid of an
    annotated element (`meta_<id>` by default), which is prefixed as well.
    """
    for name in ("compartments", "species", "parameters", "reactions"):
        for element in getattr(md, name) or []:
            if element.annotations and element.metaId is None:
                element.metaId = f"meta_{md.sid}__{element.sid}"
            port = getattr(element, "port", None)
            if port is None or port is False:
                continue
            if port is True:
                port = Port(
                    f"{element.sid}{PORT_SUFFIX}",
                    name=f"Port of {element.sid}",
                    sboTerm=SBO.PORT,
                )
            else:
                port = copy.copy(port)
            port.metaId = f"{md.sid}__{port.sid}"
            element.port = port
