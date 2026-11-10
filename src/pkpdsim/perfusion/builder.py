"""The build of a machine perfusion model: one flat SBML model of a definition.

`build` validates a `PerfusionDefinition`, composes the model (`compose`)
and writes `<sid>.xml`, SBML L3V2 (the level of the models of pkdb_models)
(the flat model of a liver, its comp package resolved): the circuit
(`pkpdsim.perfusion.circuit`), the zones of the organ with the cellular core
(`pkpdsim.perfusion.core`) and the observers (`pkpdsim.perfusion.observers`). Every element has a name and an
SBO term, every initial assignment its name; the written file is validated
with the units.

The organ is a chain of zones in series with the flow, the first zone first
(periportal in the liver). The zones of a liver are the zonation engine's
(`pkpdsim.zonation`, `_liver`): the well-stirred geometry with the vascular
space `V_vas_<zone>` and the cells `V_cell_<zone>` of the fractions of
the organ, the flow between the zones (`Q_li`, the pump flow of the circuit,
which is the boundary of the chain), the zonal factors of the profiles and the
cellular core (`pkpdsim.perfusion.cell_model`) as a submodel `M<kk>` of every
zone, composed through SBML comp and flattened (`build` writes both and the
cell model; the flat model is simulated). The states of the cells which the
liver observes are exposed under the ids of `pkpdsim.perfusion.ids`
(`atp_z1_cell`, `f_alive_z1`, `O2IM_z1`, see `compose(expose=)` of the engine),
the other elements of the cells are `M<kk>__<id>` in the flat model, the
constants of the core are shared by the zones. The zones of the kidney and of a generic organ are emitted
here, by the same core with the same ids, as are the circuit, the pump flow
and the kidney. An organ with the token `<o>` (`li` the liver, `ki` the
kidney, the key of a generic organ) has the mass `M_<o>` (kg), the density
`rho_<o>`, and the volume `V_<o> = M_<o>/rho_<o>` and the vascular and cell
fractions `f_vas_<o>` and `f_cell_<o>` as constants (the liver those of the
engine, `f_vas` and `f_cell`, and `V_li` for reference); the vascular space and
the cells of each of its `n` zones are the compartments `V_vas_<zone> =
f_vas_<o>*V_<o>/n` and `V_cell_<zone> = f_cell_<o>*V_<o>/n` (initial
assignments, so a change of the mass before a simulation sets the zones). A
profile of a process (`profiles` of the organ) is normalised to the zonal
factors `p_<process>_<zone> = n*p/sum(p)`, so the total capacity of the organ,
the rate per litre of cells times the cell volume of the zones, is that of one
zone of the whole organ, as in the zonated liver of the PBPK library; a
uniform profile is 1 in every zone. The cells release the quantities of the
organ (`PerfusionDefinition.released`), and the enzymes of the cells, the
perfusate and the observers are the released ones.

The liver secretes bile (`bile`, `bile_rate`). The kidney filters the
perfusate (`_kidney`): the GFR `GFR_ki` is that of an ample flow,
`GFR_max_ki = gfr_ki*M_ki` (`gfr_ki` in l/min/kg, numerically ml/min/g), at
most the maximal filtration fraction `FF_max_ki` of the pump flow (a rule,
so the kidney stops filtering with the pump), shared by the zones by their
volume; the tubules reabsorb the fraction `f_reab_ki` of the filtered water
and small solutes and none of a tracer, and the net urine flow `Q_urine =
(1 - f_reab_ki)*GFR_ki` (a rule) leaves the perfusate. The volumes of the
tubing and of the zones stay constant: the flow out of zone `k` of `n` is
`Q_pump - k/n*Q_urine`, out of the outflow tubing `Q_pump - Q_urine`, and
the reservoir loses the urine, its volume a rate rule `-Q_urine`, while the
urine volume `V_urine` gains it, so the perfusate and the urine keep their
water exactly. The urine holds the amounts `<q>_urine` (mmole, it starts
empty) and gives the observers `urine` (µl) and `urine_rate` (µl/min). A
generic organ has neither bile nor urine.
"""

import dataclasses
import math
import re
import textwrap
from collections.abc import Mapping
from pathlib import Path
from typing import NamedTuple

import libsbml
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    Model,
    Parameter,
    RateRule,
    Reaction,
    Species,
    create_model,
)
from sbmlutils.metadata import BQB, SBO
from sbmlutils.validation import ValidationOptions, validate_doc

from pkpdsim.pbpk.annotation import TAXA, initial_assignment
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.perfusion import circuit, core, ids, observers
from pkpdsim.perfusion.cell_model import ARRHENIUS, cell_model
from pkpdsim.perfusion.cell_model import constants as cell_constants
from pkpdsim.perfusion.definition import (
    KidneyOrgan,
    LiverOrgan,
    Organ,
    PerfusionDefinition,
)
from pkpdsim.perfusion.units import PU
from pkpdsim.zonation import composition
from pkpdsim.zonation import ids as zonation_ids
from pkpdsim.zonation.builder import ZonationBuilder
from pkpdsim.zonation.flow import ImposedFlow
from pkpdsim.zonation.geometry import Layer, geometry
from pkpdsim.zonation.profiles import Weights
from pkpdsim.zonation.settings import FlowBoundary, Geometry, ZonationSettings

__all__ = [
    "build",
    "compose",
    "model_ids",
    "tracer_ids",
    "zone_cell_model",
    "zone_factors",
]

_P = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


def zone_factors(organ: Organ) -> dict[str, dict[str, float]]:
    """Get the zonal factors of the processes with a profile, by zone.

    Args:
        organ: the organ.

    Returns:
        The factor `n*p/sum(p)` of every process with a profile in every
        zone (`p_<process>_<zone>`): the factors of a profile sum to `n`, so
        the zones of `1/n` of the volume have the capacity of the organ.
    """
    n = organ.zones
    out: dict[str, dict[str, float]] = {ids.zone_id(k): {} for k in range(1, n + 1)}
    for process, values in organ.weights().items():
        total = sum(values)
        for k, p in enumerate(values, start=1):
            out[ids.zone_id(k)][process] = n * p / total
    return out


_CELLS = {"liver": "hepatocytes", "kidney": "cells", "organ": "cells"}
_VASCULAR = {
    "liver": "vascular space (sinusoids)",
    "kidney": "vascular space (capillaries)",
    "organ": "vascular space",
}


class _Organ(NamedTuple):
    """The compartments and objects of the zones of an organ."""

    compartments: list[Compartment]
    objects: core.CoreObjects
    bile: list[str]
    """The formulas of the bile secretion of the zones (a liver)."""
    comp: Model | None = None
    """The composed model of the liver (its submodels), which the model keeps."""


def _zone_core(
    definition: PerfusionDefinition, zone: str, profile: Mapping[str, float]
) -> core.CoreObjects:
    """Get the cellular core of a zone which the builder emits itself (no liver)."""
    return core.core_objects(
        core.CoreNames(zone),
        ids.cells(zone),
        ids.vascular(zone),
        definition.core,
        profile=profile or None,
        composition=circuit.composition(definition.perfusate),
        bile=False,
        released=definition.released,
    )


def _compartment(
    sid: str, fraction: float, volume: float, n: int, text: str, k: int, what: str
) -> Compartment:
    """Get the compartment of the vascular space or the cells of a zone."""
    return Compartment(
        sid,
        fraction * volume / n,
        unit=U.litre,
        name=f"{text} of the {what} zone {k} of {n}",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
    )


def _organ(definition: PerfusionDefinition) -> _Organ:
    """Get the compartments and objects of the zones of the organ.

    The zones of a liver are the zonation engine's (`_liver`); the zones of
    the other organs are emitted here.

    Returns:
        The compartments of the zones, the objects (the constants of the
        organ, the cellular core of every zone, the tracers and the fresh
        perfusate of its vascular space) and the formulas of the bile
        secretion of the zones of a liver.
    """
    organ = definition.organ.resolved(definition.species)
    if isinstance(organ, LiverOrgan):
        return _liver(definition, organ)
    key, n, what = organ.key, organ.zones, organ.what
    m = circuit.converted(organ.mass, 1e-3, U.kg, "g")
    rho, f_vas, f_cell = organ.density, organ.vascular_fraction, organ.cell_fraction
    assert rho is not None and f_vas is not None and f_cell is not None
    volume = m.value / rho.value
    out = core.CoreObjects(
        parameters=[
            value_parameter(ids.mass(key), m, U.kg, f"mass of the {what}", _P),
            value_parameter(
                ids.density(key),
                rho,
                U.kg_per_l,
                f"density of the {what}",
                SBO.MASS_DENSITY_OF_AN_ENTITY,
            ),
            value_parameter(
                ids.vascular_fraction(key),
                f_vas,
                U.dimensionless,
                f"vascular fraction of the {what}",
                _P,
            ),
            value_parameter(
                ids.cell_fraction(key),
                f_cell,
                U.dimensionless,
                f"cell fraction of the {what}",
                _P,
            ),
            value_parameter(
                ids.organ_volume(key), volume, U.litre, f"volume of the {what}", _P
            ),
        ],
        assignments=[
            initial_assignment(
                ids.organ_volume(key),
                f"{ids.mass(key)}/{ids.density(key)}",
                U.litre,
                f"volume of the {what}",
            )
        ],
    )
    compartments: list[Compartment] = []
    factors = zone_factors(organ)
    for k in range(1, n + 1):
        zone = ids.zone_id(k)
        cells, vascular = ids.cells(zone), ids.vascular(zone)
        for sid, fraction, fid, where in (
            (cells, f_cell, ids.cell_fraction(key), _CELLS[what]),
            (vascular, f_vas, ids.vascular_fraction(key), _VASCULAR[what]),
        ):
            compartments.append(
                _compartment(sid, fraction.value, volume, n, where, k, what)
            )
            out.assignments.append(
                initial_assignment(
                    sid,
                    f"{fid} * {ids.organ_volume(key)}/{n} dimensionless",
                    U.litre,
                    str(compartments[-1].name),
                )
            )
        out += _zone_core(definition, zone, factors[zone])
        out.species += [
            circuit.tracer_species(t, zone, vascular)
            for t in definition.perfusate.tracers
        ]
        out += circuit.fresh_assignments(
            zone, definition.perfusate, definition.released
        )
    return _Organ(compartments, out, [])


class _Quantity(NamedTuple):
    """A quantity of the perfusate as a solute of the chain: an id and a name."""

    sid: str
    name: str
    mr: None = None
    diffusion: None = None


def zone_cell_model(definition: PerfusionDefinition) -> Model | None:
    """Get the cell model of the zones of the liver of a definition, `None` for no liver.

    The cell model has the zonal factor of every process with a profile as a
    port and no bile compartment (the liver of the builder has no bile of the
    engine), and is named after the model (`<sid>_core`, the stem of its file), so
    that the models of two definitions in one folder have a file each.
    """
    organ = definition.organ
    if not isinstance(organ, LiverOrgan):
        return None
    cell = cell_model(
        definition.core,
        bile=True,
        released=definition.released,
        processes=tuple(organ.profiles),
        bile_port=False,
    )
    cell.sid = f"{definition.sid}_core"
    return cell


def _cell_names(cell: Model, zone: str) -> dict[str, str]:
    """Get the names in the liver of the states of the cell model of a zone.

    The liver reaches a state of its cells through a proxy of its own model
    (`zonation.composition.compose(expose=)`) whose id survives the flattening,
    so it is the id of the flat model: the names of `pkpdsim.perfusion.ids`
    (`<quantity>_<zone>_cell`, `<quantity>_<zone>_dead`, `f_alive_<zone>`,
    `v_nec_<zone>`, `Vbile_<zone>`, `pH_<zone>`, `<reaction>_<zone>`).

    The states of the cell are its species of the cells and the necrotic cells,
    its reactions and its living fraction, necrosis, bile and pH. What the
    liver has of its own (the temperature factor `f_T`, the concentrations in
    the vascular space, the constants) is not exposed: the flat model has it
    once, and the copies in the submodels stay `M<kk>__<id>`.

    Returns:
        The id of the element in the cell model to its name in the liver.
    """
    names: dict[str, str] = {}
    for x in cell.species or []:
        sid = str(x.sid)
        for site in (ids.cell_site(zone), ids.dead_site(zone)):
            suffix = site[len(zone) + 1 :]
            if sid.endswith(f"_{suffix}"):
                names[sid] = ids.species(sid[: -len(suffix) - 1], site)
    names.update({str(r.sid): ids.reaction(str(r.sid), zone) for r in cell.reactions})
    names["BUF_ext"] = ids.buffer(zone)
    names.update(
        {
            "f_alive": ids.alive(zone),
            "v_nec": ids.necrosis_rate(zone),
            "Vbile": ids.bile(zone),
            "pH_ext": ids.ph(zone),
        }
    )
    return names


def _zone_formula(
    formula: str, names: Mapping[str, str], replaced: Mapping[str, str]
) -> str:
    """Get a formula of the cell model as it reads in the liver of a zone.

    Every element of the cell model is read under its name in the liver
    (`names`, the proxies of its states; `replaced`, the ports the zone
    replaces and what the liver has of its own, as the elements of the zone).
    """

    def rename(token: re.Match[str]) -> str:
        t = token.group(0)
        return replaced.get(t) or names.get(t) or t

    return re.sub(r"[A-Za-z_][A-Za-z0-9_]*", rename, formula)


def _liver(definition: PerfusionDefinition, organ: LiverOrgan) -> _Organ:
    """Get the zones of the liver from the zonation engine.

    The engine chains the zones (the vascular space and the cells of every
    zone, the flow between the zones, the zonal factors of the profiles) and
    places the cellular core as a cell model, a submodel of every zone; the
    circuit is the boundary of the chain (`FlowBoundary.EXTERNAL`) with the
    pump flow as the flow of the liver, and the organ keeps its own
    fractions of the vascular space and the cells (`ZonationSettings.fractions`).
    The perfusate species of the zones are those of the perfusion (the
    annotations, the start of the fresh perfusate), under the ids of the engine.

    Returns:
        The compartments, the objects of the engine and the fresh perfusate,
        the formulas of the bile secretion of the zones and the composed
        model.
    """
    n, released = organ.zones, definition.released
    cell = zone_cell_model(definition)
    assert cell is not None
    quantities = circuit.quantities(definition.perfusate, released)
    assert organ.density is not None
    assert organ.vascular_fraction is not None and organ.cell_fraction is not None
    settings = ZonationSettings(
        n_zones=n,
        geometry=Geometry.WELL_STIRRED,
        flow=ImposedFlow(
            circuit.converted(definition.circuit.flow, 1e-3, U.l_per_min, "ml/min")
        ),
        bile=None,
        profiles={process: Weights(w) for process, w in organ.weights().items()},
        species=definition.species,
        liver_mass=organ.mass,
        geometry_constants=dataclasses.replace(
            geometry(definition.species), f_par=1.0, density=organ.density.value
        ),
        fractions={
            Layer.PLASMA: organ.vascular_fraction.value,
            Layer.CELL: organ.cell_fraction.value,
        },
        boundary=FlowBoundary.EXTERNAL,
        cell_species=False,
    )
    builder = ZonationBuilder(
        settings,
        [_Quantity(q, circuit.quantity_label(q)) for q in quantities],
        cell,
        cell_processes=tuple(organ.profiles),
    )
    temperature = next(
        p
        for p in circuit.temperature_objects(definition.temperature).parameters
        if p.sid == ids.TEMPERATURE
    )
    shared = [
        *(
            str(p.sid)
            for p in cell_constants(definition.core, bile=True, released=released)
        ),
        *ARRHENIUS,
    ]
    (rule,) = [r for r in cell.rate_rules if r.variable == "Vbile"]

    def replaced_by(zone: str) -> dict[str, str]:
        """The ports of the cell model which the zone replaces, by its elements."""
        return {
            **{sid: sid for sid in shared},
            "Vli": ids.cells(zone),
            ids.F_T: ids.F_T,
            **{
                f"p_{process}": ids.profile(process, zone) for process in organ.profiles
            },
        }

    table = {ids.zone_id(k): _cell_names(cell, ids.zone_id(k)) for k in range(1, n + 1)}

    def names(zone: str) -> dict[str, str]:
        return table[zone]

    model = composition.compose(
        builder,
        cell.sid,
        composition.CellPorts(),
        exchange=ids.perfusate_quantities(released),
        temperature=temperature,
        shared=shared,
        expose={
            sid: (lambda zone, sid=sid: table[zone][sid]) for sid in sorted(table["z1"])
        },
        cell_model=cell,
    )
    # the species of the vascular space are those of the perfusion
    start = circuit.composition(definition.perfusate)
    c = start.concentrations(released)
    h, hco3 = c["h"], c["hco3"]
    assert h is not None and hco3 is not None
    ph = {ids.ph(ids.zone_id(k)): 3 - math.log10(h) for k in range(1, n + 1)}
    for p in model.parameters:
        if p.sid in ph and p.value is None:
            # The pH of the cells is an assignment rule of the cell model, `pH =
            # 3 - log10(h/1 mM)`, which the flat model keeps as the rule of the
            # element of the submodel (the proxy is merged with it): this value is
            # what the composed model carries, no simulation reads it.
            p.value = ph[str(p.sid)]
    species: dict[str, Species] = {}
    out = core.CoreObjects()
    bile: list[str] = []
    for k in range(1, n + 1):
        zone = ids.zone_id(k)
        vascular = ids.vascular(zone)
        perfusate = core.perfusate_objects(zone, vascular, start, released)
        zone_species = perfusate.species
        zone_species += [
            circuit.tracer_species(t, zone, vascular)
            for t in definition.perfusate.tracers
        ]
        for x in zone_species:
            if x.sid == ids.species("co2", zone):
                # The CO2 of the vascular space starts at the equilibrium of the
                # buffer, `h*hco3/Ka`. The flat model starts it with the initial
                # assignment of the cell model's `co2_ext`, which replaces this
                # species (it takes the assignment of the replacing element): this
                # value is what the composed model carries, no simulation reads it.
                x.initialConcentration = h * hco3 / core.Buffer.default().ka.value
        species.update({str(x.sid): x for x in zone_species})
        out += circuit.fresh_assignments(zone, definition.perfusate, released)
        bile.append(_zone_formula(str(rule.value), names(zone), replaced_by(zone)))
    out += core.CoreObjects(
        species=[species.get(str(x.sid), x) for x in model.species],
        reactions=list(model.reactions),
        parameters=list(model.parameters),
        assignments=list(model.assignments),
    )
    # the volume of the organ, `V_li = M_li/rho_li`, which `ids.organ_volume` names
    key, what = organ.key, organ.what
    volume = circuit.converted(organ.mass, 1e-3, U.kg, "g").value / organ.density.value
    out.parameters.append(
        value_parameter(
            ids.organ_volume(key), volume, U.litre, f"volume of the {what}", _P
        )
    )
    out.assignments.append(
        initial_assignment(
            ids.organ_volume(key),
            f"{ids.mass(key)}/{ids.density(key)}",
            U.litre,
            f"volume of the {what}",
        )
    )
    # the flow of the liver is the pump flow of the circuit
    (flow,) = [p for p in out.parameters if p.sid == zonation_ids.FLOW]
    flow.constant = False
    flow.value = None
    out.rules.append(
        AssignmentRule(
            zonation_ids.FLOW,
            ids.Q_PUMP,
            U.l_per_min,
            name="blood flow through the liver, the flow of the pump",
            sboTerm=SBO.MATHEMATICAL_EXPRESSION,
        )
    )
    return _Organ(list(model.compartments), out, bile, model)


def _kidney(
    definition: PerfusionDefinition, organ: KidneyOrgan
) -> tuple[Compartment, core.CoreObjects, dict[str, str]]:
    """Get the filtration and the urine of a kidney.

    The GFR of an ample flow `GFR_max_<organ> = gfr_<organ>*M_<organ>` (an
    initial assignment) is capped at the maximal filtration fraction of the
    pump flow, `GFR_<organ> = piecewise(GFR_max_<organ>, GFR_max_<organ> <=
    FF_max_<organ>*Q_pump, FF_max_<organ>*Q_pump)` (a rule, so the
    filtration follows a change of the flow and stops with the pump), and
    shared by the zones by their volume; each zone filters its small
    solutes (`ids.FILTERED`) at `(1 - f_reab_<organ>)*GFR_<organ>/n` times
    their concentration into the urine (the net of the filtration and the
    tubular reabsorption, which takes the same fraction of the water and the
    small solutes), its tracers at `GFR_<organ>/n` (not reabsorbed). The net urine flow
    `Q_urine = (1 - f_reab_<organ>)*GFR_<organ>` (a rule) leaves the
    perfusate, at most the fraction `FF_max_<organ>` of the flow: the
    flow out of zone `k` is `Q_pump - k/n*Q_urine`, out of the outflow tubing
    `Q_pump - Q_urine`, and the volume of the reservoir falls at `Q_urine`
    while the urine volume `V_urine` rises at it (rate rules), so the
    perfusate and the urine keep their volume.

    Returns:
        The urine compartment `V_urine`, the objects (the constants, the
        urine species, the filtration reactions, the outflows of the zones
        and the outflow tubing, the rate rules of the reservoir and the
        urine) and the outflow of every site whose flow is not the pump flow.
    """
    key, n = organ.key, organ.zones
    zones = [ids.zone_id(k) for k in range(1, n + 1)]
    gfr, f_reab, ff = organ.gfr, organ.reabsorption, organ.ff_max
    assert gfr is not None and f_reab is not None and ff is not None
    urine = ids.volume(ids.URINE_SITE)
    m, q_urine = ids.mass(key), ids.Q_URINE
    gfr_id, f_id = ids.gfr(key), ids.reabsorption(key)
    max_id, ff_id = ids.gfr_max(key), ids.ff_max(key)
    tracers = definition.perfusate.tracers
    out = core.CoreObjects(
        parameters=[
            value_parameter(
                ids.gfr_specific(key),
                circuit.converted(gfr, 1.0, PU.l_per_min_kg, "ml/min/g"),
                PU.l_per_min_kg,
                "glomerular filtration rate per mass of the kidney",
                _P,
            ),
            value_parameter(
                f_id,
                f_reab,
                U.dimensionless,
                "reabsorbed fraction of the filtrate of the kidney",
                _P,
            ),
            value_parameter(
                ff_id,
                ff,
                U.dimensionless,
                "maximal filtration fraction of the flow of the kidney",
                _P,
            ),
            Parameter(
                max_id,
                gfr.value * organ.mass.value * 1e-3,
                U.l_per_min,
                constant=True,
                name="glomerular filtration rate of the kidney at an ample flow",
                sboTerm=_P,
            ),
            Parameter(
                gfr_id,
                None,
                U.l_per_min,
                constant=False,
                name="glomerular filtration rate of the kidney",
                sboTerm=_P,
            ),
            Parameter(
                q_urine,
                None,
                U.l_per_min,
                constant=False,
                name="net urine flow of the kidney",
                sboTerm=_P,
            ),
        ],
        assignments=[
            initial_assignment(
                max_id,
                f"{ids.gfr_specific(key)} * {m}",
                U.l_per_min,
                "glomerular filtration rate of the kidney at an ample flow",
            ),
        ],
        rules=[
            AssignmentRule(
                gfr_id,
                f"piecewise({max_id}, {max_id} <= {ff_id} * {ids.Q_PUMP}, "
                f"{ff_id} * {ids.Q_PUMP})",
                U.l_per_min,
                name="glomerular filtration rate of the kidney",
                sboTerm=SBO.MATHEMATICAL_EXPRESSION,
                notes="The GFR of an ample flow, at most the maximal filtration "
                "fraction of the pump flow.",
            ),
            AssignmentRule(
                q_urine,
                f"(1 dimensionless - {f_id}) * {gfr_id}",
                U.l_per_min,
                name="net urine flow of the kidney",
                sboTerm=SBO.MATHEMATICAL_EXPRESSION,
            ),
        ],
        rate_rules=[
            RateRule(
                ids.volume(ids.RESERVOIR),
                f"-{q_urine}",
                U.l_per_min,
                name="water of the urine from the reservoir",
                sboTerm=SBO.MATHEMATICAL_EXPRESSION,
            ),
            RateRule(
                urine,
                q_urine,
                U.l_per_min,
                name="urine volume",
                sboTerm=SBO.MATHEMATICAL_EXPRESSION,
            ),
        ],
    )
    for q in (*ids.FILTERED, *tracers):
        out.species.append(
            Species(
                ids.species(q, ids.URINE_SITE),
                urine,
                initialAmount=0.0,
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=True,
                name=f"{circuit.quantity_label(q)} (urine)",
                sboTerm=SBO.SIMPLE_CHEMICAL,
                notes="An amount: the urine starts empty.",
            )
        )
        fraction = "" if q in tracers else f"(1 dimensionless - {f_id}) * "
        for zone in zones:
            a = ids.species(q, zone)
            out.reactions.append(
                Reaction(
                    ids.filtration(q, zone),
                    f"{a} -> {ids.species(q, ids.URINE_SITE)}",
                    formula=(
                        f"{fraction}{gfr_id}/{n} dimensionless * {a}",
                        U.mmole_per_min,
                    ),
                    compartment=ids.vascular(zone),
                    name=f"urine of {q} ({zone})",
                    sboTerm=SBO.TRANSPORT_REACTION,
                    notes=(
                        "The filtration, not reabsorbed (a tracer)."
                        if q in tracers
                        else "The filtration net of the tubular reabsorption."
                    ),
                )
            )
    outflows: dict[str, str] = {}
    for k, zone in enumerate(zones, start=1):
        outflows[zone] = ids.outflow(zone)
        _flow(out, zone, f"{ids.Q_PUMP} - {k}/{n} dimensionless * {q_urine}")
    outflows[ids.OUTFLOW] = ids.outflow(ids.OUTFLOW)
    _flow(out, ids.OUTFLOW, f"{ids.Q_PUMP} - {q_urine}")
    compartment = Compartment(
        urine,
        0.0,
        unit=U.litre,
        constant=False,
        name="urine of the kidney",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
    )
    return compartment, out, outflows


def _flow(objects: core.CoreObjects, site: str, formula: str) -> None:
    """Add the flow of the perfusate out of a site, `Q_<site>` (l/min)."""
    name = f"flow of the perfusate out of {site}"
    objects.parameters.append(
        Parameter(
            ids.outflow(site),
            None,
            U.l_per_min,
            constant=False,
            name=name,
            sboTerm=_P,
        )
    )
    objects.rules.append(
        AssignmentRule(
            ids.outflow(site),
            formula,
            U.l_per_min,
            name=name,
            sboTerm=SBO.MATHEMATICAL_EXPRESSION,
        )
    )


def _sites(definition: PerfusionDefinition) -> list[str]:
    """Get the sites of the perfusate in the order of the flow."""
    zones = [ids.zone_id(k) for k in range(1, definition.organ.zones + 1)]
    return [ids.RESERVOIR, ids.INFLOW, *zones, ids.OUTFLOW]


def tracer_ids(definition: PerfusionDefinition, tracer: str) -> set[str]:
    """Get the ids which a tracer adds to the model of a definition.

    Args:
        definition: the definition.
        tracer: the token of the tracer.

    Returns:
        The species at every site and mean, the fresh concentration and the
        flows; in a kidney the urine and the filtration of every zone.
    """
    sites = _sites(definition)
    out = {ids.species(tracer, s) for s in (*sites, ids.PERIPORTAL, ids.PERICENTRAL)}
    out.add(ids.fresh(tracer))
    out |= {
        ids.flow(tracer, a, b)
        for a, b in zip(sites, [*sites[1:], sites[0]], strict=True)
    }
    if isinstance(definition.organ, KidneyOrgan):
        out.add(ids.species(tracer, ids.URINE_SITE))
        out |= {ids.filtration(tracer, s) for s in sites if ids.is_zone(s)}
    return out


def model_ids(model: Model) -> set[str]:
    """Get the ids of the compartments, species, parameters and reactions of a model."""
    groups = (model.compartments, model.species, model.parameters, model.reactions)
    return {str(e.sid) for group in groups for e in group or []}


def compose(definition: PerfusionDefinition) -> Model:
    """Compose the model of a definition.

    Args:
        definition: the definition.

    Returns:
        The model: the circuit, the zones of the organ and the observers.

    Raises:
        ValueError: if the definition is invalid (`PerfusionDefinition.validate`).
    """
    problems = definition.validate()
    if problems:
        raise ValueError(
            f"Invalid perfusion definition '{definition.sid}':\n- "
            + "\n- ".join(problems)
        )
    organ = definition.organ.resolved(definition.species)
    perfusate = definition.perfusate
    released = definition.released
    enzymes = [e for e in ids.ENZYMES if e in released]
    zones = [ids.zone_id(k) for k in range(1, organ.zones + 1)]
    sites = _sites(definition)

    objects = circuit.temperature_objects(definition.temperature)
    objects += circuit.oxygen_objects()
    objects += core.CoreObjects(parameters=core.buffer_parameters())
    objects += circuit.perfusate_parameters(perfusate, released)
    objects += circuit.circuit_objects(definition.circuit, perfusate, released)
    zones_of = _organ(definition)
    objects += zones_of.objects
    compartments = [
        *circuit.circuit_compartments(
            definition.circuit, reservoir_constant=not isinstance(organ, KidneyOrgan)
        ),
        *zones_of.compartments,
    ]
    outflows: dict[str, str] = {}
    if isinstance(organ, KidneyOrgan):
        urine, kidney, outflows = _kidney(definition, organ)
        compartments.append(urine)
        objects += kidney
        objects += observers.urine_objects()
    objects += core.CoreObjects(
        parameters=core.core_parameters(
            definition.core,
            released=released,
            bile=isinstance(organ, LiverOrgan),
        )
    )
    objects += circuit.flow_objects(
        sites, perfusate, released, outflows, zone_flows=zones_of.comp is None
    )
    objects += core.CoreObjects(parameters=observers.enzyme_parameters(enzymes))
    for site in sites:
        objects += observers.site_objects(site, enzymes)
    objects += observers.mean_objects(
        zones,
        circuit.quantities(perfusate, released),
        ids.cell_quantities(released),
        enzymes,
    )
    if isinstance(organ, LiverOrgan):
        objects += observers.bile_objects(zones_of.bile, zones, ids.mass(organ.key))
    objects += observers.uptake_objects(zones)

    notes = (
        f"Machine perfusion of the {organ.what} ({organ.key}) of the "
        f"{definition.species} with {organ.zones} zone(s) in a circuit "
        "(pkpdsim.perfusion)."
    )
    if definition.notes:
        # the notes of a model package are an indented block, as its terms of use
        notes = f"{notes}\n\n{textwrap.dedent(definition.notes).strip()}"
    model = Model(
        sid=definition.sid,
        name=definition.name,
        notes=notes,
        annotations=[
            (BQB.HAS_TAXON, TAXA[definition.species]),
            (
                BQB.HAS_PROPERTY,
                "mamo/MAMO_0000046",
            ),  # ordinary differential equation model
        ],
        creators=list(definition.creators),
        units=PU,
        model_units=MODEL_UNITS,
        compartments=compartments,
        species=objects.species,
        parameters=_unique(objects.parameters),
        reactions=objects.reactions,
        rules=objects.rules,
        rate_rules=objects.rate_rules,
        assignments=objects.assignments,
    )
    if zones_of.comp is not None:
        # the zones of the liver are the submodels of the cell model
        comp = zones_of.comp
        model.external_model_definitions = comp.external_model_definitions
        model.submodels = comp.submodels
        model.replaced_elements = comp.replaced_elements
        model.deletions = comp.deletions
    return model


def _unique(parameters: list[Parameter]) -> list[Parameter]:
    """Get the parameters without a repeated id, the first of an id kept.

    The temperature of the circuit is the one of the zones of a liver, which
    the engine adds to its model as well.
    """
    seen: dict[str, Parameter] = {}
    for parameter in parameters:
        seen.setdefault(str(parameter.sid), parameter)
    return list(seen.values())


def build(definition: PerfusionDefinition, output_dir: Path) -> Path:
    """Build a definition: write the SBML model of the circuit.

    The model of an organ without submodels is written as `<sid>.xml`. The
    liver has the cellular core as a submodel of every zone: the composed model
    `<sid>.xml`, the cell model `<sid>_core.xml` and the flat model, the
    composition resolved, `<sid>_flat.xml` (the model which is simulated), by
    `pkpdsim.zonation.composition.write`.

    Args:
        definition: the definition.
        output_dir: the folder of the model.

    Returns:
        The path of the model which is simulated: `<sid>.xml`, `<sid>_flat.xml`
        for a liver.

    Raises:
        ValueError: if the definition is invalid or the model has SBML
            errors (the units included).
    """
    model = compose(definition)
    cell = zone_cell_model(definition)
    if cell is None:
        path = output_dir / f"{definition.sid}.xml"
        # SBML L3V2, the level of the models of pkdb_models: the names of the rules
        # and initial assignments are written (L3V1 has no name on them)
        create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    else:
        path = composition.write(model, output_dir, cell_model=cell)
    doc = libsbml.readSBMLFromFile(str(path))
    validation = validate_doc(doc, options=ValidationOptions(), title=str(path))
    if validation.error_count:
        raise ValueError(
            f"The model '{path}' has {validation.error_count} SBML errors."
        )
    return path
