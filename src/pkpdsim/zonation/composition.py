"""A cell model per zone, composed through SBML comp.

The hepatocytes of a zone are a standalone cell model, which the engine
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

The source of the external model definition is `<cell_model_id>.xml`, relative
to the folder of the model that is written, so that the written pair (the comp
model and the cell model) is portable. `write` requires that file and raises a
`FileNotFoundError` naming the expected path. `compose` and `ZonationBuilder.build`
stay separate: `build` is the engine model without a cell model.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from sbmlutils.comp import flatten_sbml
from sbmlutils.factory import (
    PORT_SUFFIX,
    Deletion,
    ExternalModelDefinition,
    Model,
    Parameter,
    ReplacedElement,
    Submodel,
    create_model,
)
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.units import U
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
        cell: the compartment of the cell, replaced by the hepatocyte.
        bile: the compartment of the bile, replaced by the canaliculus.
        position: the parameter of the lobular position (0 periportal, 1
            pericentral).
        temperature: the parameter of the temperature (degC).
        concentration: the pattern of the concentration port of a substance
            in the external compartment, with the substance as `{sid}`.
    """

    external: str = "Vext"
    cell: str = "Vli"
    bile: str = "Vbi"
    position: str = "position"
    temperature: str = "temperature"
    concentration: str = "{sid}_ext"


def _port(sid: str) -> str:
    """The id of the port of an element."""
    return f"{sid}{PORT_SUFFIX}"


def compose(
    builder: ZonationBuilder,
    cell_model_id: str,
    ports: CellPorts,
    deletions: Sequence[str] = (),
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
            `ports.concentration` of every substance of the builder.
        deletions: the ids of ports of the cell model, deleted from every
            submodel (what the liver model replaces by its own).

    Returns:
        The model of `builder.build()` with the external model definition,
        one submodel per zone, the replaced elements and the deletions.
    """
    # the units of the perfusion library when the engine composes and not
    # when it is imported, as `zonation_units` of the builder takes them
    # (the three libraries import in any order)
    from pkpdsim.perfusion.units import PU

    model = builder.build()
    s = builder.settings
    n = s.n_zones
    disse = ids.VASCULAR if s.geometry is Geometry.WELL_STIRRED else Layer.DISSE
    model.external_model_definitions = [
        ExternalModelDefinition(
            sid=cell_model_id,
            name="cell model of the hepatocytes",
            source=f"{cell_model_id}.xml",
            modelRef=cell_model_id,
        )
    ]
    model.parameters.append(
        Parameter(
            ids.TEMPERATURE,
            REFERENCE_TEMPERATURE,
            unit=PU.degC,
            name="temperature of the liver, the same in all zones",
            constant=True,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        )
    )
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
            (ids.volume(Layer.HEPATOCYTE, zone), ports.cell),
            (ids.position(zone), ports.position),
            (ids.TEMPERATURE, ports.temperature),
        ]
        if s.bile is not None:
            wiring.append((ids.volume(Layer.BILE, zone), ports.bile))
        wiring.extend(
            (
                ids.species(substance.sid, disse, zone),
                ports.concentration.format(sid=substance.sid),
            )
            for substance in builder.substances
        )
        for element, port in wiring:
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


def write(model: Model, output_dir: Path) -> Path:
    """Write the comp model, flatten it and return the flat model.

    The comp model is written as `<model id>.xml` and the flat model as
    `<model id>_flat.xml` into `output_dir`, where the file of every external
    model definition has to be.

    Args:
        model: the composed model (`compose`).
        output_dir: the folder.

    Returns:
        The path of the flat model.

    Raises:
        FileNotFoundError: if the file of an external model definition is not
            in `output_dir`.
    """
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
