"""A cell model per zone, composed through SBML comp.

The cells of a zone (the hepatocytes of the liver) are a standalone cell model, which the engine
composes once per zone: one `ExternalModelDefinition`, one `Submodel` `M<kk>`
per zone (`ids.submodel`) and replaced elements which wire the compartments
and species of the zone onto the ports of the cell model. The cell model is
kept as its own model, because the individual cell models are coupled to
meshes and other architectures; the flat model is what is simulated.

The wiring is anatomical:

* The **external port** (`CellPorts.external`, `Vext`) is replaced by the
  **space of Disse** of the zone, not by the plasma: a hepatocyte faces the
  space of Disse, and the fenestrated barrier of the radial transport
  (`ZonationBuilder.radial_exchange`) is what separates the Disse space from
  the lumen. In the well-stirred geometry the Disse space is merged into the
  vascular site, which replaces the port instead.
* The **cell port** (`CellPorts.cell`, `Vli`) is replaced by the
  **hepatocyte** and the **bile port** (`CellPorts.bile`, `Vbi`) by the
  **bile canaliculus** (not without a bile).
* The concentration port `<sid>_ext` of every substance is replaced by the
  species of the substance in the Disse layer of the zone, so the cell model
  sees what the barrier delivers.
* `position` is replaced by a parameter `position_<zone>` per zone with the
  centre of the zone, from `profiles.zone_positions`, the function the
  profiles use, so the cell model and the profile factors agree on where a
  zone is.
* `temperature` is replaced by one global parameter of the liver, the same in
  every submodel: the liver is isothermal (a gradient would be a profile).
* Ids the caller replaces by its own model are deleted from every submodel
  (`deletions`).
* A **state of the cell which the liver observes** (`expose`) is reached through
  comp as well. The flat id of an element of a submodel, `M01__glc_cell`, exists
  in no model until the flattening has happened, so the math of the liver cannot
  name it. The engine adds to the liver one element per zone and exposed id, a
  proxy which the caller names (`expose` maps the id in the cell model to a
  function of the zone), of the same class as the element of the cell model (a
  species in the compartment of the cell model's species, replaced onto the
  same port compartment of the zone; a parameter; a reaction) with its unit,
  and wires it with `comp:ReplacedBy` onto the element of the submodel. The
  proxy is a name, not a state: it copies the class, units, compartment, flags and
  initial value of the element of the cell model when the model is built (so it cannot
  drift from it), has no rule or kinetic law, and the flattening merges it with the
  element of the submodel: **the id which survives in the flat model is the proxy's**
  (libsbml keeps the id of the parent element of a `ReplacedBy`), so the flat id of an
  exposed element is the name the caller gives, and every element which is not
  exposed stays `M<kk>__<id>`. The engine does not know how the caller names a
  state.

The source of the external model definition is `<cell_model_id>.xml`, relative
to the folder of the model that is written, so that the written pair (the comp
model and the cell model) is portable. `write` requires that file and raises a
`FileNotFoundError` naming the expected path. `compose` and `ZonationBuilder.build`
stay separate: `build` is the engine model without a cell model.
"""

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import libsbml
from sbmlutils.comp import flatten_sbml
from sbmlutils.factory import (
    PORT_SUFFIX,
    Deletion,
    ExternalModelDefinition,
    Model,
    Parameter,
    Reaction,
    ReplacedBy,
    ReplacedElement,
    Species,
    Submodel,
    create_model,
)
from sbmlutils.metadata import SBO

from pkpdsim.sbml.units import OrganUnits, U
from pkpdsim.zonation import ids
from pkpdsim.zonation.builder import ZonationBuilder
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.profiles import zone_positions
from pkpdsim.zonation.settings import Geometry

__all__ = ["CellPorts", "compose", "write"]

REFERENCE_TEMPERATURE = 37.0
"""The temperature (degC) of the liver."""


@dataclass(frozen=True)
class CellPorts:
    """The ids of the ports of a cell model.

    Attributes:
        external: the compartment of the surroundings of the cell, replaced by
            the space of Disse.
        cell: the compartment of the cell, replaced by the cell layer of the zone.
        bile: the compartment of the bile, replaced by the canaliculus.
        position: the parameter of the lobular position (0 periportal, 1
            pericentral).
        temperature: the parameter of the temperature (degC).
        concentration: the pattern of the concentration port of a substance
            in the external compartment, with the substance as `{sid}`.
        profile: the pattern of the port of the zonal factor of a process of
            the cell model, with the process as `{process}`.
    """

    external: str = "Vext"
    cell: str = "Vli"
    bile: str = "Vbi"
    position: str = "position"
    temperature: str = "temperature"
    concentration: str = "{sid}_ext"
    profile: str = "p_{process}"


def _port(sid: str) -> str:
    """The id of the port of an element."""
    return f"{sid}{PORT_SUFFIX}"


def _proxy(
    cell_model: Model,
    cell_id: str,
    proxy_id: str,
    sm: str,
    compartments: Mapping[str, str],
) -> Species | Parameter | Reaction:
    """Get the proxy in the liver of an element of the cell model.

    Args:
        cell_model: the cell model, which gives the class, the unit and the
            compartment of the element.
        cell_id: the id of the element in the cell model.
        proxy_id: the id of the proxy in the liver.
        sm: the submodel of the zone.
        compartments: the compartment of the liver of every compartment of the
            cell model which a zone replaces (the ports).

    Returns:
        The proxy, replaced by the element of the submodel.

    Raises:
        ValueError: if the cell model has no such element, or a species in a
            compartment which no port of the zone replaces.
    """
    by = ReplacedBy(submodelRef=sm, idRef=cell_id)
    for x in cell_model.species or []:
        if x.sid == cell_id:
            if x.compartment not in compartments:
                raise ValueError(
                    f"The exposed species '{cell_id}' is in the compartment "
                    f"'{x.compartment}' of the cell model, which no port of the "
                    f"zone replaces ({sorted(compartments)})."
                )
            return Species(
                proxy_id,
                compartments[x.compartment],
                initialAmount=x.initialAmount,
                initialConcentration=x.initialConcentration,
                substanceUnit=x.substanceUnits,
                hasOnlySubstanceUnits=x.hasOnlySubstanceUnits,
                constant=x.constant,
                boundaryCondition=x.boundaryCondition,
                name=x.name,
                sboTerm=x.sboTerm,
                replacedBy=by,
            )
    for p in cell_model.parameters or []:
        if p.sid == cell_id:
            return Parameter(
                proxy_id,
                p.value,
                unit=p.unit,
                constant=p.constant,
                name=p.name,
                sboTerm=p.sboTerm,
                replacedBy=by,
            )
    for r in cell_model.reactions or []:
        if r.sid == cell_id:
            return Reaction(
                proxy_id,
                equation="",
                reversible=r.reversible,
                fast=r.fast,
                name=r.name,
                sboTerm=r.sboTerm,
                replacedBy=by,
            )
    raise ValueError(f"The cell model '{cell_model.sid}' has no element '{cell_id}'.")


def compose(
    builder: ZonationBuilder,
    cell_model_id: str,
    ports: CellPorts,
    deletions: Sequence[str] = (),
    exchange: Sequence[str] | None = None,
    temperature: Parameter | None = None,
    shared: Sequence[str] = (),
    expose: Mapping[str, Callable[[str], str]] | None = None,
    cell_model: Model | None = None,
) -> Model:
    """Compose the model of a builder with a cell model per zone.

    Without a bile (`settings.bile is None`) the bile port is not replaced, so
    every submodel keeps its own local bile compartment (`M01__Vbi`, one more
    compartment per zone in the flat model); this is expected.

    Args:
        builder: the builder of the zonated liver; its `build()` is the model.
        cell_model_id: the id of the cell model, which is the id of its model
            and the stem of its file `<cell_model_id>.xml` next to the model.
        ports: the ports of the cell model. The cell model needs the port
            `ports.concentration` of every substance of `exchange`.
        deletions: the ids of ports of the cell model, deleted from every
            submodel (what the liver model replaces by its own).
        exchange: the ids of the solutes of the builder which the cell model
            has a concentration port of; `None` for all (a tracer which the
            cell does not see has none).
        temperature: the parameter of the temperature of the liver, which the
            caller defines (a circuit has one, which the zones follow); `None`
            adds the default one, `temperature` at 37 degC; the zones follow
            the parameter whatever its id.
        shared: the ids of parameters of the cell model (ports) which the
            liver replaces in every zone by its own parameter of the same id,
            which the caller adds to the model: a constant of the cell, one
            for all zones.
        expose: the elements of the cell model (species, parameters, reactions)
            which the math of the liver references: the id in the cell model
            to a function of the zone token (`z1`) which gives the id of the
            proxy in the liver. The engine adds one proxy per zone and id,
            replaced by the element of the submodel, so that the liver names
            the state of a cell by an id of its own model and not by the id
            the flattening gives it (`M01__<id>`, which is the id of the flat
            model, unchanged).
        cell_model: the cell model, which `expose` needs for the class, the unit
            and the compartment of the elements.

    Returns:
        The model of `builder.build()` with the external model definition,
        one submodel per zone, the replaced elements, the deletions and the
        proxies of the exposed elements.

    Raises:
        ValueError: if `expose` is given without the cell model, or names an
            element the cell model does not have.
    """
    if expose and cell_model is None:
        raise ValueError("`expose` needs the cell model (`cell_model`).")
    model = builder.build()
    s = builder.settings
    n = s.n_zones
    disse = ids.VASCULAR if s.geometry is Geometry.WELL_STIRRED else Layer.DISSE
    species = {x.sid: x for x in model.species}
    model.external_model_definitions = [
        ExternalModelDefinition(
            sid=cell_model_id,
            name="cell model of the cells of a zone",
            source=f"{cell_model_id}.xml",
            modelRef=cell_model_id,
        )
    ]
    model.parameters.append(
        temperature
        if temperature is not None
        else Parameter(
            ids.TEMPERATURE,
            REFERENCE_TEMPERATURE,
            unit=OrganUnits.degC,
            name="temperature of the liver, the same in all zones",
            constant=True,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        )
    )
    exchanged = (
        [x for x in builder.substances if x.sid in set(exchange)]
        if exchange is not None
        else list(builder.substances)
    )
    scaled = [p for p in builder.cell_processes if p in s.profiles]
    for k, x in zip(range(1, n + 1), zone_positions(n), strict=True):
        zone = ids.zone(k)
        sm = ids.submodel(k)
        model.submodels.append(
            Submodel(
                sid=sm,
                modelRef=cell_model_id,
                name=f"cell model of zone {k} of {n}",
            )
        )
        model.parameters.append(
            Parameter(
                ids.position(zone),
                x,
                unit=U.dimensionless,
                name=f"lobular position of zone {k} (0 periportal, 1 pericentral)",
                constant=True,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            )
        )
        wiring: list[tuple[str, str]] = [
            (ids.volume(disse, zone), ports.external),
            (ids.volume(Layer.CELL, zone), ports.cell),
            (ids.position(zone), ports.position),
            (
                temperature.sid if temperature is not None else ids.TEMPERATURE,
                ports.temperature,
            ),
        ]
        if s.bile is not None:
            wiring.append((ids.volume(Layer.BILE, zone), ports.bile))
        wiring.extend(
            (
                ids.species(substance.sid, disse, zone),
                ports.concentration.format(sid=substance.sid),
            )
            for substance in exchanged
        )
        wiring.extend((sid, sid) for sid in shared)
        wiring.extend(
            (ids.profile(process, zone), ports.profile.format(process=process))
            for process in scaled
        )
        if expose:
            assert cell_model is not None
            home = {
                ports.external: ids.volume(disse, zone),
                ports.cell: ids.volume(Layer.CELL, zone),
            }
            if s.bile is not None:
                home[ports.bile] = ids.volume(Layer.BILE, zone)
            for cell_id, name in expose.items():
                proxy = _proxy(cell_model, cell_id, name(zone), sm, home)
                if isinstance(proxy, Species):
                    model.species.append(proxy)
                elif isinstance(proxy, Parameter):
                    model.parameters.append(proxy)
                else:
                    model.reactions.append(proxy)
        for element, port in wiring:
            if element in species:
                # an annotated species of the cell model has a metaid, which
                # comp requires of the species which replaces it as well
                species[element].metaId = species[element].metaId or ids.meta_id(
                    element
                )
            rid = f"{sm}_{element}_{port}_RE"
            model.replaced_elements.append(
                ReplacedElement(
                    sid=rid,
                    metaId=rid,
                    elementRef=element,
                    submodelRef=sm,
                    portRef=_port(port),
                )
            )
        for sid in deletions:
            did = f"{sm}_{sid}_DE"
            model.deletions.append(
                Deletion(
                    sid=did,
                    metaId=did,
                    submodelRef=sm,
                    portRef=_port(sid),
                )
            )
    return model


def _model_id(path: Path) -> str:
    """The id of the model in a file (empty for a file without a model)."""
    document = libsbml.readSBMLFromFile(str(path))
    model = document.getModel()
    return model.getId() if model is not None else ""


def _write_cell_model(cell_model: Model, model: Model, output_dir: Path) -> None:
    """Write the cell model beside the composed model, never over another model.

    A file which holds a version of the same model (the id of its model is that
    of the cell model) is rewritten, as the composed model is.

    Raises:
        ValueError: if the composed model does not use the cell model.
        FileExistsError: if another model is in the file of the cell model.
    """
    used = {d.modelRef for d in model.external_model_definitions}
    if cell_model.sid not in used:
        raise ValueError(
            f"The composed model '{model.sid}' uses the cell models {sorted(used)}, "
            f"not '{cell_model.sid}'."
        )
    target = output_dir / f"{cell_model.sid}.xml"
    if target.is_file() and _model_id(target) != cell_model.sid:
        raise FileExistsError(
            f"'{target}' holds the model '{_model_id(target)}', not the cell "
            f"model '{cell_model.sid}': it is not overwritten."
        )
    create_model(
        cell_model, filepath=target, sbml_level=3, sbml_version=2, validate=False
    )


def write(model: Model, output_dir: Path, cell_model: Model | None = None) -> Path:
    """Write the comp model, flatten it and return the flat model.

    The comp model is written as `<model id>.xml` and the flat model as
    `<model id>_flat.xml` into `output_dir`, where the file of every external
    model definition has to be: the cell model is written there as
    `<cell model id>.xml` when it is given, over a file which holds a version
    of the same model but never over a file which holds another model.

    Args:
        model: the composed model (`compose`).
        output_dir: the folder.
        cell_model: the cell model of the composition; `None` requires its
            file in `output_dir` already.

    Returns:
        The path of the flat model.

    Raises:
        FileNotFoundError: if the file of an external model definition is not
            in `output_dir`.
        FileExistsError: if the file of the cell model holds another model.
        ValueError: if the composed model does not use the cell model.
    """
    if cell_model is not None:
        _write_cell_model(cell_model, model, output_dir)
    for definition in model.external_model_definitions:
        source = output_dir / definition.source
        if not source.is_file():
            raise FileNotFoundError(
                f"The cell model '{definition.modelRef}' is expected at "
                f"'{source}' (the source of its external model definition is "
                f"relative to the folder of the model)."
            )
    path = output_dir / f"{model.sid}.xml"
    flat_path = output_dir / f"{model.sid}_flat.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    flatten_sbml(path, flat_path, leave_ports=True, validate=False)
    return flat_path
