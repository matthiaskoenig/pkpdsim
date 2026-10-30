"""The observers of a PBPK model: assignment rules of the top model, no states.

An observer reads the pool species and the outflow concentrations of the
organs and changes nothing, so it conserves mass. Per substance `<s>`:

- `C<pool>_free_<s> = fu*C<pool>_<s>`, the unbound plasma concentration of a
  pool (`ar`, `ve`, `po`);
- `C<pool>_blood_<s>`, the blood concentration of a pool,
  `(1-HCT)*C + HCT*C_bc`: `BP*C` unless the blood cells are kinetic;
- `Chv_<s>`, the hepatic vein: the outflow of the liver mixed with the
  portosystemic shunt by their flows,
  `((1-f_shunt)*Qh*Cout_li + f_shunt*(Qha*Car + Qpo*Cpo))/Qh`, with
  `Chv_free_<s>` and `Chv_blood_<s>`;
- `Carm_<s> = sum w*C`, the arm vein (peripheral venous) sampling site: the
  venous plasma outflows of organs and the arterial plasma mixed with the
  weights `w_arm_<site>` of a preset of `data/arm_vein.toml`
  (`PBPKDefinition.arm_vein`, by default the mix of skin and muscle of the
  PK-Sim documentation); `Carm_free_<s>` and `Carm_blood_<s>` as for the
  pools.

The tissue concentrations are the species of the organs themselves, their
unbound concentrations `Cu_<part>_<s>` observers of the organs, the amounts
in urine and feces the sink species, the urinary excretion rate
`Rurine_<s>` an output of the kidney.

A definition adds observers of its own (`Observer`, `PBPKDefinition.observers`,
e.g. the sum of two metabolites): parameters of the top model over the ids
of the top model, appended after the observers of the library. An observer
names a species or an unbound concentration of an organ by its flat id
(`MU__Cti_<s>`, `LI__Cu_ic_<s>`, an amount `MU__A<x>`), which the organ
exposes as a PD input is exposed (`pkpdsim.pbpk.builder`). An observer which
reads a state, directly or through other observers (a species, an observer
of the library, an output of an organ or a PD component, a dynamic
parameter, `time`), is an assignment rule; an observer of constants only is
a constant with an initial assignment, evaluated after the changes before
the simulation (an assignment rule of constants would slow the compile of
the model, R53).

The intestine adds readouts on request (`PBPKDefinition.readouts`,
`intestine_readouts`): the dissolved amount in its lumen `Alumen_<s>` and
the summed fluxes of the absorption `Rabs_<s>`, the efflux into the lumen
`Reff_<s>` and the loss to the feces `Rfeces_<s>`, outputs of the intestine
as `Rurine_<s>` is one of the kidney.
"""

from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import libsbml
from sbmlutils.factory import (
    PORT_SUFFIX,
    AssignmentRule,
    ModelDefinition,
    Parameter,
    Port,
    PortType,
    Sbase,
    UnitDefinition,
)
from sbmlutils.factory.units import UnitType
from sbmlutils.metadata import SBO

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.annotation import annotate, initial_assignment
from pkpdsim.pbpk.data import load, value
from pkpdsim.pbpk.organs import Organ, Part, Pool
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.substance import Substance
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value

if TYPE_CHECKING:
    from pkpdsim.pbpk.definition import PBPKDefinition

__all__ = [
    "Observer",
    "OrganOutput",
    "arm_vein_weights",
    "definition_observer_objects",
    "formula_functions",
    "formula_symbols",
    "intestine_readouts",
    "observer_objects",
    "readout_ids",
]

_ONE = "1 dimensionless"
_SBO_PARAMETER = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER
_OBSERVER = "An observer: an assignment rule which changes no state."


@dataclass(frozen=True)
class Observer:
    """An observer of a definition: a parameter of the top model with its formula.

    Attributes:
        sid: the id of the observer, an id of the top model.
        formula: the formula (SBML L3 infix) over the ids of the top model:
            the pools, the sinks, the observers of the library and of the
            definition, the covariates, the outputs of the organs, the
            readouts; and a species or an unbound concentration of an organ
            by its flat id (`MU__Cti_<s>`), which the organ exposes.
        unit: the unit of the formula; a unit which is no unit of `U` (the
            unit per body surface area `ML_PER_MIN_1P73M2`) is a unit
            definition of the model.
        name: the name of the observer.
        notes: the notes of the observer.
    """

    sid: str
    formula: str
    unit: UnitDefinition
    name: str
    notes: str | None = None


@dataclass(frozen=True)
class OrganOutput:
    """An output of an organ which replaces the parameter of the top model of its id.

    The parameter is created by the composition (`pkpdsim.pbpk.builder`):
    an element of an organ which a PD component or an observer reads, a
    readout of the intestine.

    Attributes:
        unit: the unit of the parameter.
        name: the name of the parameter.
        kind: the kind of its annotation (`pkpdsim.pbpk.annotation.KINDS`).
        text: the start of its notes, followed by `organ <organ>, replaced
            by the output of the organ`.
        constant: the output is a constant of the organ (a constant of the
            generated steady state), not an assignment rule.
    """

    unit: UnitType
    name: str
    kind: str
    text: str
    constant: bool = False


TIME = "time"
"""The symbol of the simulation time in the symbols of a formula."""


def formula_symbols(formula: str) -> set[str]:
    """Get the ids a formula reads, with libsbml's parser of SBML L3 infix.

    Args:
        formula: the formula; the units of its numbers (`1 dimensionless`)
            and the names of its functions are no ids.

    Returns:
        The ids of the formula, `time` for the simulation time.

    Raises:
        ValueError: if the formula cannot be parsed, with the message of
            libsbml.
    """
    found: set[str] = set()
    for node in _nodes(formula):
        if node.getType() == libsbml.AST_NAME_TIME:
            found.add(TIME)
        elif node.getType() == libsbml.AST_NAME:
            found.add(node.getName())
    return found


def formula_functions(formula: str) -> set[str]:
    """Get the functions a formula calls which are no functions of SBML L3 infix.

    Such a call (`foo(x)`) reads a function definition of the model; the
    models of the library define none.

    Args:
        formula: the formula.

    Returns:
        The names of the called functions which are no built-in functions
        (`exp`, `piecewise`, ...).

    Raises:
        ValueError: if the formula cannot be parsed, with the message of
            libsbml.
    """
    return {
        node.getName()
        for node in _nodes(formula)
        if node.getType() == libsbml.AST_FUNCTION
    }


def _nodes(formula: str) -> list[libsbml.ASTNode]:
    """Get the nodes of a formula, parsed with libsbml's parser of SBML L3 infix.

    Raises:
        ValueError: if the formula cannot be parsed, with the message of
            libsbml.
    """
    ast = libsbml.parseL3Formula(formula)
    if ast is None:
        raise ValueError(
            f"formula '{formula}' cannot be parsed: "
            f"{libsbml.getLastParseL3Error().strip()}"
        )
    nodes: list[libsbml.ASTNode] = []
    pending = [ast]
    while pending:
        node = pending.pop()
        nodes.append(node)
        pending.extend(node.getChild(i) for i in range(node.getNumChildren()))
    return nodes


def definition_observer_objects(
    observers: Iterable[Observer], constant: Collection[str] = ()
) -> list[Sbase]:
    """Get the parameters and the rules of the observers of a definition.

    The ids, the symbols of the formulas and whether an observer reads a
    state are checked when the model is composed (`pkpdsim.pbpk.builder`).

    Args:
        observers: the observers of the definition, their formulas over the
            ids of the top model.
        constant: the ids of the observers which read constants only: a
            constant with an initial assignment, every other observer an
            assignment rule.

    Returns:
        The parameters with their assignment rules or initial assignments,
        annotated as observers.
    """
    objects: list[Sbase] = []
    for o in observers:
        text = f"{o.notes} " if o.notes else ""
        fixed = o.sid in constant
        what = (
            "A constant of constants with an initial assignment, evaluated after "
            "the changes before the simulation."
            if fixed
            else _OBSERVER
        )
        objects += [
            annotate(
                Parameter(
                    o.sid,
                    None,
                    o.unit,
                    constant=fixed,
                    name=o.name,
                    sboTerm=_SBO_PARAMETER,
                    notes=f"{text}An observer of the definition. {what}\n\n"
                    f"`{o.sid} = {o.formula}`",
                ),
                "observer",
            ),
            initial_assignment(o.sid, o.formula, o.unit, o.name)
            if fixed
            else AssignmentRule(o.sid, o.formula, o.unit, name=o.name),
        ]
    return objects


_READOUT_KINDS = (
    (ids.lumen_amount, "lumen"),
    (ids.absorption_rate, "absorption"),
    (ids.efflux_rate, "efflux"),
    (ids.feces_rate, "feces"),
)
"""The readouts of the intestine by their id function and what they sum."""


def readout_ids(substances: Iterable[str]) -> dict[str, tuple[str, str]]:
    """Get the ids of the readouts of the intestine of substances.

    Args:
        substances: the substance ids.

    Returns:
        By readout id (`Alumen_<s>`, `Rabs_<s>`, `Reff_<s>`, `Rfeces_<s>`),
        what it sums (`lumen`, `absorption`, `efflux`, `feces`) and its
        substance.
    """
    return {f(s): (what, s) for s in substances for f, what in _READOUT_KINDS}


def intestine_readouts(
    md: ModelDefinition, readouts: Sequence[str], substances: Iterable[str]
) -> dict[str, OrganOutput]:
    """Add the readouts of a definition to the model definition of its intestine.

    Every readout is a parameter of the intestine with its assignment rule
    and an output port, which replaces the parameter of its id in the top
    model, as the urinary excretion rate of the kidney does. The sums are
    taken over the species and the reactions of the intestine, so that they
    hold for every intestine with a lumen (segmented, single lumen):

    - `Alumen_<s>`: the dissolved amount in the lumen, `sum V*C` of the
      lumen species `Clumen_<seg>_<s>` (mmole);
    - `Rabs_<s>`: the fluxes from a lumen species into the wall of the
      substance (the mucosa `Cmu_<seg>_<s>`, else a species of the wall), the
      absorption and an apical influx (mmole/min);
    - `Reff_<s>`: the fluxes from the wall into the lumen, the apical efflux
      (mmole/min);
    - `Rfeces_<s>`: the net flux of the intestine into the feces `Afeces_<s>`
      (mmole/min); with an intestine with a lumen every flux into the feces
      is one of the intestine, so it integrates to `Afeces_<s>`.

    Args:
        md: the model definition of the intestine, changed in place.
        readouts: the readouts of the definition (`PBPKDefinition.readouts`).
        substances: the substance ids of the definition.

    Returns:
        The outputs by id.

    Raises:
        ValueError: naming a readout of a substance without lumen species, or
            whose sum has no flux.
    """
    known = readout_ids(substances)
    species = {str(sp.sid): sp for sp in md.species or []}
    problems: list[str] = []
    outputs: dict[str, OrganOutput] = {}
    tokens = (None, *ids.GUT_SEGMENTS)
    for rid in readouts:
        what, s = known[rid]
        lumen = [ids.lumen(s, seg) for seg in tokens if ids.lumen(s, seg) in species]
        if what != "feces" and not lumen:
            problems.append(
                f"readout '{rid}': the substance '{s}' has no lumen species in the "
                "intestine (it reaches the lumen by an oral route, the bile or a "
                "transport of the intestine)"
            )
            continue
        if what == "lumen":
            formula = " + ".join(f"{species[c].compartment} * {c}" for c in lumen)
        else:
            formula = _flux_sum(md, what, s, set(lumen))
            if not formula:
                problems.append(
                    f"readout '{rid}': the intestine has no flux {_FLUXES[what]} "
                    f"of '{s}'"
                )
                continue
        output = _readout_output(what, s)
        md.parameters = [
            *(md.parameters or []),
            annotate(
                Parameter(
                    rid,
                    None,
                    output.unit,
                    constant=False,
                    name=output.name,
                    sboTerm=_SBO_PARAMETER,
                    notes=f"{_READOUT_NOTES[what]} A readout of the definition, an "
                    f"output of the intestine.\n\n`{rid} = {formula}`",
                    port=Port(f"{rid}{PORT_SUFFIX}", portType=PortType.OUTPUT_PORT),
                ),
                output.kind,
                organ=Organ.SMALL_INTESTINE,
            ),
        ]
        md.rules = [*(md.rules or []), AssignmentRule(rid, formula, output.unit)]
        outputs[rid] = output
    if problems:
        raise ValueError("; ".join(problems) + ".")
    return outputs


_FLUXES = {
    "absorption": "from the lumen into the wall",
    "efflux": "from the wall into the lumen (an apical efflux of the small intestine)",
    "feces": "into the feces",
}
"""What the fluxes of a readout are, in a message."""

_READOUT_NOTES = {
    "lumen": "The dissolved amount in the intestinal lumen, the sum over its "
    "segments (the undissolved drug of a tablet is not in it).",
    "absorption": "The rate of the absorption: the sum of the fluxes from the "
    "lumen into the intestinal wall (the absorption of every segment and an "
    "apical influx).",
    "efflux": "The rate of the efflux into the lumen: the sum of the fluxes from "
    "the intestinal wall into the lumen (the apical efflux of every segment).",
    "feces": "The rate of the fecal excretion: the net flux of the intestine "
    "into the feces.",
}
"""The notes of the parameter of a readout in the intestine."""


def _readout_output(what: str, s: str) -> OrganOutput:
    """The output of a readout: its unit, name, kind and the notes in the top model."""
    if what == "lumen":
        return OrganOutput(
            U.mmole,
            f"dissolved amount {s} in the intestinal lumen [mmole]",
            "lumen_amount",
            "Dissolved amount in the lumen, the sum over its segments, of the ",
        )
    process, kind = {
        "absorption": ("absorption", "intestine_absorption"),
        "efflux": ("efflux into the intestinal lumen", "intestine_efflux"),
        "feces": ("fecal excretion", "feces_rate"),
    }[what]
    return OrganOutput(
        U.mmole_per_min,
        f"rate of the {process} {s} [mmole/min]",
        kind,
        f"Rate of the {process}, a sum of the fluxes, of the ",
    )


def _stoichiometry(part: object) -> float:
    """The stoichiometry of a reactant or a product of a reaction, 1 if unset."""
    value = getattr(part, "stoichiometry", None)
    return 1.0 if value is None else float(value)


def _flux_sum(md: ModelDefinition, what: str, s: str, lumen: set[str]) -> str:
    """The sum of the fluxes of the intestine of a readout, empty without any.

    The amount which the reactions move into the targets (the wall of the
    absorption, the lumen of the efflux, the feces), each reaction weighted
    by the net stoichiometry of the targets: the reactions of the absorption
    take from the lumen, those of the efflux from the wall, any reaction
    counts for the feces (and one out of the feces less). The wall of a
    substance: its mucosa species and its species of the parts of the wall
    (`ids.species`).
    """
    tokens = (None, *ids.GUT_SEGMENTS)
    wall = {ids.mucosa(s, seg) for seg in tokens} | {ids.species(p, s) for p in Part}
    sources, targets = {
        "absorption": (lumen, wall),
        "efflux": (wall, lumen),
        "feces": (None, {ids.feces(s)}),
    }[what]
    terms: list[str] = []
    for r in md.reactions or []:
        reactants = r.equation.reactants
        if sources is not None and not {p.species for p in reactants} & sources:
            continue
        net = sum(
            _stoichiometry(p) for p in r.equation.products if p.species in targets
        )
        net -= sum(_stoichiometry(p) for p in reactants if p.species in targets)
        if net:
            sign = "-" if net < 0 else "+"
            factor = "" if abs(net) == 1 else f"{abs(net):g} dimensionless * "
            terms.append(f"{sign} {factor}{r.sid}")
    if not any(term.startswith("+") for term in terms):
        return ""
    terms.sort(key=lambda term: term.startswith("-"))
    return " ".join(terms).removeprefix("+ ")


def arm_vein_weights(preset: str) -> dict[Organ | Pool, Value]:
    """Get the weights of a preset of the arm vein.

    Args:
        preset: the name of the preset, e.g. `pksim_default`.

    Returns:
        The weight of every organ outflow or pool of the mix with its source.

    Raises:
        ValueError: if the preset is unknown or names neither an organ nor a
            pool.
    """
    data = load("arm_vein")
    if preset not in data:
        raise ValueError(
            f"Unknown arm_vein '{preset}', the presets are {sorted(data)}."
        )
    return {_site(preset, key): value("arm_vein", preset, key) for key in data[preset]}


def _site(preset: str, key: str) -> Organ | Pool:
    """Get the organ or the pool of a key of a preset, `skin`, `arterial`."""
    for enum in (Organ, Pool):
        if key.upper() in enum.__members__:
            return enum[key.upper()]
    raise ValueError(
        f"arm_vein.toml [{preset}] {key}: neither an organ nor a pool of the body."
    )


def observer_objects(definition: "PBPKDefinition") -> list[Sbase]:
    """Get the observers of a definition.

    Args:
        definition: the definition.

    Returns:
        The weights of the arm vein and the parameters with their assignment
        rules of every substance.
    """
    weights = arm_vein_weights(definition.arm_vein)
    objects: list[Sbase] = [
        annotate(
            value_parameter(
                ids.arm_weight(site),
                w,
                U.dimensionless,
                f"weight of {_site_label(site)} in the arm vein",
                _SBO_PARAMETER,
                f"Share of the {_site_label(site)} in the arm vein sample "
                f"(preset '{definition.arm_vein}').",
            ),
            "arm_weight",
            organ=site if isinstance(site, Organ) else None,
            pool=site if isinstance(site, Pool) else None,
        )
        for site, w in weights.items()
    ]
    for s in definition.substances:
        objects += _substance(s, weights)
    return objects


def _site_label(site: Organ | Pool) -> str:
    """The venous outflow of an organ or the plasma of a pool in a name."""
    if isinstance(site, Pool):
        return f"{site.name.lower()} plasma"
    return f"venous outflow of the {site.name.lower().replace('_', ' ')}"


def _observer(
    pid: str,
    formula: str,
    name: str,
    text: str,
    kind: str,
    pool: Pool | None = None,
    site: str | None = None,
) -> list[Sbase]:
    """A concentration (mM) with its assignment rule, annotated as `kind`."""
    return [
        annotate(
            Parameter(
                pid,
                None,
                U.mM,
                constant=False,
                name=name,
                sboTerm=_SBO_PARAMETER,
                notes=f"{text} {_OBSERVER}\n\n`{pid} = {formula}`",
            ),
            kind,
            pool=pool,
            site=site,
        ),
        AssignmentRule(pid, formula, U.mM),
    ]


def _blood(s: Substance, plasma: str, cells: str) -> str:
    """The blood concentration of the plasma and blood cell concentrations, as formula."""
    ery = s.erythrocytes
    if ery is not None and ery.mode == "kinetic":
        return f"({_ONE} - {ids.HCT}) * {plasma} + {ids.HCT} * {cells}"
    return f"{ids.bp(s.sid)} * {plasma}"


def _hepatic(s: Substance, cells: bool) -> str:
    """The hepatic vein mix of the liver outflow and the shunt, as formula.

    `((1-f_shunt)*Qh*Cout_li + f_shunt*(Qha*Car + Qpo*Cpo))/Qh` of the plasma,
    or of the kinetic blood cells.
    """
    sid = s.sid
    part = Part.BLOOD_CELLS if cells else Part.PLASMA
    liver = (
        ids.organ_outflow_bc(Organ.LIVER, sid)
        if cells
        else ids.organ_outflow(Organ.LIVER, sid)
    )
    car = ids.pool_species(Pool.ARTERIAL, sid, part)
    cpo = ids.pool_species(Pool.PORTAL, sid, part)
    return (
        f"(({_ONE} - {ids.F_SHUNT}) * {ids.QH} * {liver} + {ids.F_SHUNT} * "
        f"({ids.QHA} * {car} + {ids.QPO} * {cpo})) / {ids.QH}"
    )


def _substance(s: Substance, weights: dict[Organ | Pool, Value]) -> list[Sbase]:
    """The observers of a substance."""
    sid, fu = s.sid, ids.fu(s.sid)
    objects: list[Sbase] = []
    for pool in Pool:
        c = ids.pool_species(pool, sid)
        label = pool.name.lower()
        objects += _observer(
            ids.free(pool, sid),
            f"{fu} * {c}",
            f"unbound plasma concentration {sid} {label} [mM]",
            f"Unbound plasma concentration of the {label} pool.",
            "free_concentration",
            pool=pool,
        )
        objects += _observer(
            ids.blood(pool, sid),
            _blood(s, c, ids.pool_species(pool, sid, Part.BLOOD_CELLS)),
            f"blood concentration {sid} {label} [mM]",
            f"Blood concentration of the {label} pool.",
            "blood_concentration",
            pool=pool,
        )

    hv = ids.hepatic_vein(sid)
    objects += _observer(
        hv,
        _hepatic(s, cells=False),
        f"hepatic vein plasma concentration {sid} [mM]",
        "Plasma concentration of the hepatic vein: the outflow of the liver "
        "mixed with the portosystemic shunt by their flows.",
        "plasma_concentration",
        site="hepatic_vein",
    )
    objects += _observer(
        ids.hepatic_vein_free(sid),
        f"{fu} * {hv}",
        f"hepatic vein unbound plasma concentration {sid} [mM]",
        "Unbound plasma concentration of the hepatic vein.",
        "free_concentration",
        site="hepatic_vein",
    )
    ery = s.erythrocytes
    kinetic = ery is not None and ery.mode == "kinetic"
    objects += _observer(
        ids.hepatic_vein_blood(sid),
        _blood(s, hv, f"({_hepatic(s, cells=True)})" if kinetic else ""),
        f"hepatic vein blood concentration {sid} [mM]",
        "Blood concentration of the hepatic vein.",
        "blood_concentration",
        site="hepatic_vein",
    )

    def mix(cells: bool) -> str:
        terms = []
        for site in weights:
            if isinstance(site, Pool):
                c = ids.pool_species(
                    site, sid, Part.BLOOD_CELLS if cells else Part.PLASMA
                )
            elif cells:
                c = ids.organ_outflow_bc(site, sid)
            else:
                c = ids.organ_outflow(site, sid)
            terms.append(f"{ids.arm_weight(site)} * {c}")
        return " + ".join(terms)

    arm = ids.arm(sid)
    objects += _observer(
        arm,
        mix(cells=False),
        f"arm vein plasma concentration {sid} [mM]",
        "Plasma concentration of an arm vein sample: the venous outflows of the "
        "organs and the arterial plasma mixed with the weights of the preset. "
        "Plasma data without a stated sampling site map to it.",
        "plasma_concentration",
        site="arm_vein",
    )
    objects += _observer(
        ids.arm_free(sid),
        f"{fu} * {arm}",
        f"arm vein unbound plasma concentration {sid} [mM]",
        "Unbound plasma concentration of an arm vein sample.",
        "free_concentration",
        site="arm_vein",
    )
    cells = f"({mix(cells=True)})" if kinetic else ""
    objects += _observer(
        ids.arm_blood(sid),
        _blood(s, arm, cells),
        f"arm vein blood concentration {sid} [mM]",
        "Blood concentration of an arm vein sample.",
        "blood_concentration",
        site="arm_vein",
    )
    return objects
