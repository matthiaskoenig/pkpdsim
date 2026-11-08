"""The build of a machine perfusion model: one flat SBML model of a definition.

`build` validates a `PerfusionDefinition`, composes the model (`compose`)
and writes `<sid>.xml`, SBML L3V2 (the level of the models of pkdb_models)
without the comp package: the circuit (`pkpdsim.perfusion.circuit`), the
zones of the organ with the cellular core (`pkpdsim.perfusion.core`) and the
observers (`pkpdsim.perfusion.observers`). Every element has a name and an
SBO term, every initial assignment its name; the written file is validated
with the units.

The organ is a chain of zones in series with the flow, the first zone
first (periportal in the liver). An organ with the token `<o>` (`li` the
liver, `ki` the kidney, the key of a generic organ) has the mass `M_<o>`
(kg), the density `rho_<o>`, the volume `V_<o> = M_<o>/rho_<o>` and the
vascular and cell fractions `f_vas_<o>` and `f_cell_<o>` as constants; the
vascular space and the cells of each of its `n` zones are the compartments
`Vvas_<zone> = f_vas_<o>*V_<o>/n` and `Vcell_<zone> = f_cell_<o>*V_<o>/n`
(initial assignments, so a change of the mass before a simulation sets the
zones). A profile of a process (`profiles` of the organ) is normalised to
the zonal factors `p_<process>_<zone> = n*p/sum(p)`, so the total capacity
of the organ, the rate per litre of cells times the cell volume of the
zones, is that of one zone of the whole organ, as in the zonated liver of
the PBPK library; a uniform profile is 1 in every zone. The cells release
the quantities of the organ (`PerfusionDefinition.released`), and the
enzymes of the cells, the perfusate and the observers are the released ones.

The liver secretes bile (`bile`, `bile_rate`). The kidney filters the
perfusate (`_kidney`): the GFR `GFR_ki = gfr_ki*M_ki` (`gfr_ki` in l/min/kg,
numerically ml/min/g) is shared by the zones by their volume, the tubules
reabsorb the fraction `f_reab_ki` of the filtered water and small solutes
and none of a tracer, and the net urine flow `Q_urine = (1 - f_reab_ki)*
GFR_ki` leaves the perfusate. The volumes of the tubing and of the zones
stay constant: the flow out of zone `k` of `n` is `Q_pump - k/n*Q_urine`,
out of the outflow tubing `Q_pump - Q_urine`, and the reservoir loses the
urine, its volume a rate rule `-Q_urine`, while the urine volume `V_urine`
gains it, so the perfusate and the urine keep their water exactly. The
urine holds the amounts `<q>_urine` (mmole, it starts empty) and gives the
observers `urine` (µl) and `urine_rate` (µl/min). A generic organ has
neither bile nor urine.
"""

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
from pkpdsim.perfusion.definition import (
    KidneyOrgan,
    LiverOrgan,
    Organ,
    PerfusionDefinition,
)
from pkpdsim.perfusion.units import PU

__all__ = ["build", "compose", "model_ids", "tracer_ids", "zone_factors"]

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
    for process, values in organ.profiles.items():
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


def _organ(definition: PerfusionDefinition) -> _Organ:
    """Get the compartments and objects of the zones of the organ.

    Returns:
        The compartments of the zones, the objects (the constants of the
        organ, the cellular core of every zone, the tracers and the fresh
        perfusate of its vascular space) and the formulas of the bile
        secretion of the zones of a liver.
    """
    organ = definition.organ.resolved(definition.species)
    key, n, what = organ.key, organ.zones, organ.what
    released = definition.released
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
    start = circuit.composition(definition.perfusate)
    has_bile = isinstance(organ, LiverOrgan)
    bile: list[str] = []
    for k in range(1, n + 1):
        zone = ids.zone_id(k)
        cells, vascular = ids.cells(zone), ids.vascular(zone)
        for sid, fraction, fid, where in (
            (cells, f_cell, ids.cell_fraction(key), _CELLS[what]),
            (vascular, f_vas, ids.vascular_fraction(key), _VASCULAR[what]),
        ):
            name = f"{where} of the {what} zone {k} of {n}"
            compartments.append(
                Compartment(
                    sid,
                    fraction.value * volume / n,
                    unit=U.litre,
                    name=name,
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                )
            )
            out.assignments.append(
                initial_assignment(
                    sid,
                    f"{fid} * {ids.organ_volume(key)}/{n} dimensionless",
                    U.litre,
                    name,
                )
            )
        zone_objects = core.zone_objects(
            zone,
            cells,
            vascular,
            definition.core,
            profile=factors[zone] or None,
            composition=start,
            bile=has_bile,
            released=released,
        )
        bile += [
            str(r.value)
            for r in zone_objects.rate_rules
            if r.variable == ids.bile(zone)
        ]
        out += zone_objects
        out.species += [
            circuit.tracer_species(t, zone, vascular)
            for t in definition.perfusate.tracers
        ]
        out += circuit.fresh_assignments(zone, definition.perfusate, released)
    return _Organ(compartments, out, bile)


def _kidney(
    definition: PerfusionDefinition, organ: KidneyOrgan
) -> tuple[Compartment, core.CoreObjects, dict[str, str]]:
    """Get the filtration and the urine of a kidney.

    The GFR `GFR_<organ> = gfr_<organ>*M_<organ>` is shared by the zones by
    their volume; each zone filters its small solutes (`ids.FILTERED`) at
    `(1 - f_reab_<organ>)*GFR_<organ>/n` times their concentration into the
    urine (the net of the filtration and the tubular reabsorption, which
    takes the same fraction of the water and the small solutes), its
    tracers at `GFR_<organ>/n` (not reabsorbed). The net urine flow
    `Q_urine = (1 - f_reab_<organ>)*GFR_<organ>` leaves the perfusate: the
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
    gfr, f_reab = organ.gfr, organ.reabsorption
    assert gfr is not None and f_reab is not None
    urine = ids.volume(ids.URINE_SITE)
    m, q_urine = ids.mass(key), ids.Q_URINE
    gfr_id, f_id = ids.gfr(key), ids.reabsorption(key)
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
            Parameter(
                gfr_id,
                gfr.value * organ.mass.value * 1e-3,
                U.l_per_min,
                constant=True,
                name="glomerular filtration rate of the kidney",
                sboTerm=_P,
            ),
            Parameter(
                q_urine,
                (1 - f_reab.value) * gfr.value * organ.mass.value * 1e-3,
                U.l_per_min,
                constant=True,
                name="net urine flow of the kidney",
                sboTerm=_P,
            ),
        ],
        assignments=[
            initial_assignment(
                gfr_id,
                f"{ids.gfr_specific(key)} * {m}",
                U.l_per_min,
                "glomerular filtration rate of the kidney",
            ),
            initial_assignment(
                q_urine,
                f"(1 dimensionless - {f_id}) * {gfr_id}",
                U.l_per_min,
                "net urine flow of the kidney",
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
    objects += core.CoreObjects(parameters=core.core_parameters(definition.core))
    objects += circuit.flow_objects(sites, perfusate, released, outflows)
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
        objects += observers.bile_objects(zones_of.bile, zones)
    objects += observers.uptake_objects(zones)

    return Model(
        sid=definition.sid,
        name=definition.name,
        notes=(
            f"Machine perfusion of the {organ.what} ({organ.key}) of the "
            f"{definition.species} with {organ.zones} zone(s) in a circuit "
            "(pkpdsim.perfusion)."
        ),
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
        parameters=objects.parameters,
        reactions=objects.reactions,
        rules=objects.rules,
        rate_rules=objects.rate_rules,
        assignments=objects.assignments,
    )


def build(definition: PerfusionDefinition, output_dir: Path) -> Path:
    """Build a definition: write the SBML model `<sid>.xml`.

    Args:
        definition: the definition.
        output_dir: the folder of the model.

    Returns:
        The path of the model.

    Raises:
        ValueError: if the definition is invalid or the model has SBML
            errors (the units included).
    """
    model = compose(definition)
    path = output_dir / f"{definition.sid}.xml"
    # SBML L3V2, the level of the models of pkdb_models: the names of the rules
    # and initial assignments are written (L3V1 has no name on them)
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    doc = libsbml.readSBMLFromFile(str(path))
    validation = validate_doc(doc, options=ValidationOptions(), title=str(path))
    if validation.error_count:
        raise ValueError(
            f"The model '{path}' has {validation.error_count} SBML errors."
        )
    return path
