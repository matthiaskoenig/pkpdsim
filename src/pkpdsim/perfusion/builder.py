"""The build of a machine perfusion model: one flat SBML model of a definition.

`build` validates a `PerfusionDefinition`, composes the model (`compose`)
and writes `<sid>.xml`, SBML L3V2 (the level of the models of pkdb_models)
without the comp package: the circuit (`pkpdsim.perfusion.circuit`), the
zones of the organ with the cellular core (`pkpdsim.perfusion.core`) and the
observers (`pkpdsim.perfusion.observers`). Every element has a name and an
SBO term, every initial assignment its name; the written file is validated
with the units.

The organ is a chain of zones in series with the flow, periportal first.
The liver (`LiverOrgan`) has the mass `M_li` (kg), the density `rho_li`, the
volume `V_li = M_li/rho_li` and the vascular and cell fractions `f_vas_li`
and `f_cell_li` as constants; the vascular space and the cells of each of
its `n` zones are the compartments `Vvas_<zone> = f_vas_li*V_li/n` and
`Vcell_<zone> = f_cell_li*V_li/n` (initial assignments, so a change of the
liver mass before a simulation sets the zones). A profile of a process
(`LiverOrgan.profiles`) is normalised to the zonal factors `p_<process>_<zone>
= n*p/sum(p)`, so the total capacity of the organ, the rate per litre of
cells times the cell volume of the zones, is that of one zone of the whole
organ, as in the zonated liver of the PBPK library; a uniform profile is 1
in every zone.
"""

from pathlib import Path

import libsbml
from sbmlutils.factory import Compartment, Model, create_model
from sbmlutils.metadata import BQB, SBO
from sbmlutils.validation import ValidationOptions, validate_doc

from pkpdsim.pbpk.annotation import TAXA, initial_assignment
from pkpdsim.pbpk.placement import value_parameter
from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.perfusion import circuit, core, ids, observers
from pkpdsim.perfusion.definition import LiverOrgan, PerfusionDefinition
from pkpdsim.perfusion.units import PU

__all__ = ["build", "compose", "model_ids", "tracer_ids", "zone_factors"]

_P = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER


def zone_factors(organ: LiverOrgan) -> dict[str, dict[str, float]]:
    """Get the zonal factors of the processes with a profile, by zone.

    Args:
        organ: the liver.

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


def _liver(
    definition: PerfusionDefinition, organ: LiverOrgan
) -> tuple[list[Compartment], core.CoreObjects, list[str]]:
    """Get the compartments and objects of the zones of the liver.

    Returns:
        The compartments of the zones, the objects (the constants of the
        liver, the cellular core of every zone, the tracers and the fresh
        perfusate of its vascular space) and the formulas of the bile
        secretion of the zones.
    """
    key, n = organ.key, organ.zones
    organ = organ.resolved(definition.species)
    m = circuit.converted(organ.mass, 1e-3, U.kg, "g")
    rho, f_vas, f_cell = organ.density, organ.vascular_fraction, organ.cell_fraction
    assert rho is not None and f_vas is not None and f_cell is not None
    volume = m.value / rho.value
    out = core.CoreObjects(
        parameters=[
            value_parameter(ids.mass(key), m, U.kg, "mass of the liver", _P),
            value_parameter(
                ids.density(key),
                rho,
                U.kg_per_l,
                "density of the liver",
                SBO.MASS_DENSITY_OF_AN_ENTITY,
            ),
            value_parameter(
                ids.vascular_fraction(key),
                f_vas,
                U.dimensionless,
                "vascular fraction of the liver",
                _P,
            ),
            value_parameter(
                ids.cell_fraction(key),
                f_cell,
                U.dimensionless,
                "cell fraction of the liver",
                _P,
            ),
            value_parameter(
                ids.organ_volume(key), volume, U.litre, "volume of the liver", _P
            ),
        ],
        assignments=[
            initial_assignment(
                ids.organ_volume(key),
                f"{ids.mass(key)}/{ids.density(key)}",
                U.litre,
                "volume of the liver",
            )
        ],
    )
    compartments: list[Compartment] = []
    factors = zone_factors(organ)
    start = circuit.composition(definition.perfusate)
    bile: list[str] = []
    for k in range(1, n + 1):
        zone = ids.zone_id(k)
        cells, vascular = ids.cells(zone), ids.vascular(zone)
        for sid, fraction, fid, where in (
            (cells, f_cell, ids.cell_fraction(key), "hepatocytes"),
            (vascular, f_vas, ids.vascular_fraction(key), "vascular space (sinusoids)"),
        ):
            name = f"{where} of the liver zone {k} of {n}"
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
        out += circuit.fresh_assignments(zone, definition.perfusate)
    return compartments, out, bile


def _sites(definition: PerfusionDefinition) -> list[str]:
    """Get the sites of the perfusate in the order of the flow."""
    organ = definition.organ
    n = organ.zones if isinstance(organ, LiverOrgan) else 1
    zones = [ids.zone_id(k) for k in range(1, n + 1)]
    return [ids.RESERVOIR, ids.INFLOW, *zones, ids.OUTFLOW]


def tracer_ids(definition: PerfusionDefinition, tracer: str) -> set[str]:
    """Get the ids which a tracer adds to the model of a definition.

    Args:
        definition: the definition.
        tracer: the token of the tracer.

    Returns:
        The species at every site and mean, the fresh concentration and the
        flows.
    """
    sites = _sites(definition)
    out = {ids.species(tracer, s) for s in (*sites, ids.PERIPORTAL, ids.PERICENTRAL)}
    out.add(ids.fresh(tracer))
    out |= {
        ids.flow(tracer, a, b)
        for a, b in zip(sites, [*sites[1:], sites[0]], strict=True)
    }
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
    organ = definition.organ
    assert isinstance(organ, LiverOrgan)
    perfusate = definition.perfusate
    zones = [ids.zone_id(k) for k in range(1, organ.zones + 1)]
    sites = _sites(definition)

    objects = circuit.temperature_objects(definition.temperature)
    objects += circuit.oxygen_objects()
    objects += core.CoreObjects(parameters=core.buffer_parameters())
    objects += circuit.perfusate_parameters(perfusate)
    objects += circuit.circuit_objects(definition.circuit, perfusate)
    zone_compartments, organ_objects, bile = _liver(definition, organ)
    objects += organ_objects
    objects += core.CoreObjects(parameters=core.core_parameters(definition.core))
    objects += circuit.flow_objects(sites, perfusate)
    objects += core.CoreObjects(parameters=observers.enzyme_parameters())
    for site in sites:
        objects += observers.site_objects(site)
    objects += observers.mean_objects(zones, circuit.quantities(perfusate))
    objects += observers.bile_objects(bile, zones)
    objects += observers.uptake_objects(zones)

    return Model(
        sid=definition.sid,
        name=definition.name,
        notes=(
            f"Machine perfusion of the {organ.key} of the {definition.species} "
            f"with {organ.zones} zone(s) in a circuit (pkpdsim.perfusion)."
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
        compartments=[
            *circuit.circuit_compartments(definition.circuit),
            *zone_compartments,
        ],
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
