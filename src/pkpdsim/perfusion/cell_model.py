"""The cellular core as a standalone cell model with ports.

The cells of an organ zone of the machine perfusion library (`core`) as one
SBML model, which the zonation engine composes once per zone
(`pkpdsim.zonation.composition.compose`): the same objects as a zone of the
perfusion builder emits, with the same values, units, rate laws and SBO
terms, but one cell instead of a zone of an organ. There is one emitter of
the chemistry, `core.core_objects`; the zones of an organ without the engine
call it with the ids of a zone, `cell_model` with the ids of a cell.

The function is not among the lazy exports of `pkpdsim.perfusion`: it would
shadow this module of the same name; import it from here.

The ports are the contract of the cell models of the zonation
(`CellPorts`; "Cell model contract" of the zonation spec):

* the compartments `Vext` (the surroundings of the cells, which hold the
  perfusate species, the bicarbonate buffer and the pH), `Vli` (the cells)
  and `Vbi` (the bile),
* the parameters `position` (the lobular position, which the core does not
  read: its zonal factors are the profile of the organ) and `temperature`
  (degC), which sets the Arrhenius factor `f_T` of the rates, as the
  temperature of a perfusion model does, and
* one species `<q>_ext` per quantity of the perfusate
  (`ids.perfusate_quantities`), in `Vext`,
* the constants of the core and the buffer (`constants`), which the zones of a
  liver share (one parameter of the liver replaces them in every zone), and
* the zonal factor `p_<process>` of every process of `processes`.

The ids of the cell: `<q>_cell` and `<q>_dead` for the living and the necrotic
cells, `<q>_ext` for the perfusate, the reactions and the constants without a
zone (`GLCIM`, `f_alive`, `v_nec`, `Vbile`); the constants of the core, the
buffer and the Arrhenius factor are defined in the model, so that it
simulates on its own, and are ports. The compartments have placeholder sizes of a closed
example (the cells of a liver of 10 g, a perfusate of 1 ml), which a
composition replaces by the zone.
"""

from sbmlutils.factory import Compartment, Model, Parameter
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.units import MODEL_UNITS, U
from pkpdsim.perfusion import core, ids
from pkpdsim.perfusion.circuit import temperature_objects
from pkpdsim.perfusion.definition import table_value
from pkpdsim.perfusion.units import PU
from pkpdsim.zonation.composition import CellPorts

__all__ = ["CELL_MODEL_ID", "cell_model", "constants"]

CELL_MODEL_ID = "perfusion_core"
"""The id of the cell model, and the stem of its file."""

CELL_VOLUME = 0.0069
"""The cells of a liver of 10 g (l), the size of a closed example."""

EXTERNAL_VOLUME = 0.001
"""The perfusate around the cells of a closed example (l)."""

BILE_VOLUME = 1e-4
"""The bile compartment of a closed example (l)."""


def constants(
    params: core.CoreParameters,
    *,
    bile: bool = True,
    released: tuple[str, ...] = ids.RELEASABLE,
) -> list[Parameter]:
    """Get the constants of the cell model: those of the core and of the buffer.

    They are ports of the cell model: a liver which composes it into zones
    defines each once and the zones share it (`compose(shared=)`), as the zones
    of the perfusion model have always shared them.

    Args:
        params: the constants of the core.
        bile: the cell secretes bile.
        released: the quantities which the dying cells release.

    Returns:
        The parameters, `core.buffer_parameters` and `core.core_parameters`.
    """
    return [
        *core.buffer_parameters(),
        *core.core_parameters(params, released=released, bile=bile),
    ]


def cell_model(
    params: core.CoreParameters,
    *,
    bile: bool = True,
    released: tuple[str, ...] = ids.RELEASABLE,
    processes: tuple[str, ...] = (),
) -> Model:
    """Get the cellular core as a standalone model with the cell ports.

    Args:
        params: the constants of the core.
        bile: the cell secretes bile (a liver), `Vbile`; the compartment `Vbi`
            is a port in either case, so that a composition can wire it.
        released: the quantities of `ids.RELEASABLE` which the dying cells
            release; the enzymes of the cells and the perfusate are the
            released ones.
        processes: the processes of `core.PROCESSES` with a zonal factor,
            the parameter `p_<process>` (1, a port) in the rate of the
            process, which a composition replaces by the factor of the zone;
            a process without one has no factor.

    Returns:
        The model `CELL_MODEL_ID`: the objects of `core.core_objects`
        without a zone, the constants of the core, the buffer and the
        Arrhenius factor, and the ports.

    Raises:
        ValueError: if a released quantity is no quantity of
            `ids.RELEASABLE`, a process is no process of `core.PROCESSES`,
            or the constants lack what the core needs
            (`core.CoreParameters.problems`).
    """
    ports = CellPorts()
    names = core.CoreNames()
    objects = core.core_objects(
        names,
        ports.cell,
        ports.external,
        params,
        profile=dict.fromkeys(processes, 1.0),
        composition=None,
        bile=bile,
        released=released,
    )
    factors = {names.profile(process) for process in processes}
    for parameter in objects.parameters:
        if parameter.sid in factors:
            parameter.port = True
    for q in ids.perfusate_quantities(released):
        species = next(s for s in objects.species if s.sid == ids.species(q, "ext"))
        species.port = True
    temperature = temperature_objects(
        table_value("temperature", "reference"), ports.temperature
    )
    for parameter in temperature.parameters:
        if parameter.sid == ports.temperature:
            parameter.port = True
    shared = constants(params, bile=bile, released=released)
    for parameter in shared:
        parameter.port = True
    position = Parameter(
        ports.position,
        0.5,
        unit=U.dimensionless,
        name="lobular position (0 periportal, 1 pericentral)",
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        port=True,
    )
    return Model(
        sid=CELL_MODEL_ID,
        name="cellular core of the machine perfusion library",
        units=PU,
        model_units=MODEL_UNITS,
        compartments=[
            Compartment(
                sid,
                size,
                unit=U.litre,
                name=name,
                sboTerm=SBO.PHYSICAL_COMPARTMENT,
                port=True,
            )
            for sid, size, name in (
                (ports.external, EXTERNAL_VOLUME, "surroundings of the cells"),
                (ports.cell, CELL_VOLUME, "cells"),
                (ports.bile, BILE_VOLUME, "bile"),
            )
        ],
        parameters=[
            position,
            *shared,
            *temperature.parameters,
            *objects.parameters,
        ],
        species=objects.species,
        reactions=objects.reactions,
        rules=objects.rules,
        rate_rules=objects.rate_rules,
        assignments=[*temperature.assignments, *objects.assignments],
    )
