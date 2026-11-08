"""The simulation experiments of a machine perfusion: the condition and the medium changes.

`PerfusionSimulationExperiment` is the base of the experiments of a
perfusion model, a `pkpdsim.experiment.PKPDSimulationExperiment` with the
`PerfusionDefinition` of the model (`definition`), from which it takes the
ids of the condition and of the medium changes, and the labels and units of
the observables (`observable_labels`, `observable_units`), which a subclass
may replace.

The condition of a simulation is a change before the initialization
(`perfusion_changes`): the pump flow `Q_pump`, the mass of the organ
`M_<o>` (the volumes of its zones are initial assignments of it) and the
temperature `T` (the Arrhenius factor `f_T` and the solubility of oxygen are
initial assignments of it, so the temperature is a change before the
initialization). The temperature of the model is in `degC`, a unit of the
size of a kelvin without the offset (`pkpdsim.perfusion.units`), so the
change is the number of the temperature in degC and not a quantity, which
sbmlsim would convert to 310.15 for 37 degC.

A medium change replaces the fraction `f` of the perfusate of the reservoir
with fresh perfusate: every perfusate species of the reservoir, the
quantities of the perfusate which the organ has and the tracers, is set to
`(1 - f)*[X_res] + f*X_fresh` (`medium_change`), with the fresh perfusate
`X_fresh` of the model (the constants `<q>_fresh`); a `Change` of sbmlsim at
the time of the change (`PerfusionSimulationExperiment.medium_changes`), so
the model has no events and the times are a setting of the experiment. The
change keeps the volume of the reservoir: as much fresh perfusate comes in
as perfusate goes out. The reservoir of a kidney, which loses the urine
(its volume a state), keeps the volume it has at the change, the urine is
not made up. The tubing and the organ are not changed.

sbmlsim evaluates the values of a change with the state at its time and
assigns them at once; an output time at the time of a change has two rows,
the state before the change and the state after it, and a fit compares a
data point at that time with the state before
(`pkpdsim.fit.default_fit_settings`, `data_at_change`), the sample taken
before the medium is changed.
"""

import math
from collections.abc import Sequence
from types import MappingProxyType
from typing import ClassVar

from sbmlsim.simulation import Change
from sbmlsim.units import Quantity, ureg

from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.perfusion import circuit, core, ids
from pkpdsim.perfusion.definition import KidneyOrgan, LiverOrgan, PerfusionDefinition

__all__ = [
    "PerfusionSimulationExperiment",
    "medium_change",
    "observable_labels",
    "observable_units",
    "perfusion_changes",
    "site_label",
]

_SITES = {
    ids.RESERVOIR: "reservoir",
    ids.INFLOW: "inflow",
    ids.OUTFLOW: "outflow",
    ids.PERIPORTAL: "periportal",
    ids.PERICENTRAL: "pericentral",
}

_MM = "mM"


def site_label(site: str) -> str:
    """Get the label of a site of the perfusate.

    Args:
        site: a site of the circuit (`res`, `in`, `out`), a zone (`z3`) or a
            mean of the zones (`pp`, `pc`).

    Returns:
        The label, `reservoir`, `inflow`, `outflow`, `zone 3`, `periportal`,
        `pericentral`.

    Raises:
        KeyError: for another site.
    """
    if ids.is_zone(site):
        return f"zone {site[1:]}"
    return _SITES[site]


def medium_change(definition: PerfusionDefinition, fraction: float) -> dict[str, str]:
    """Get the values of a medium change of the reservoir.

    Args:
        definition: the definition of the model.
        fraction: the fraction of the reservoir which is replaced, in (0, 1].

    Returns:
        The formula `(1 - f)*[X_res] + f*X_fresh` by the concentration
        `[X_res]` of every perfusate species of the reservoir (the
        quantities of the perfusate of the organ and the tracers), the values
        of a `Change` of sbmlsim.

    Raises:
        ValueError: for a fraction which is not in (0, 1].
    """
    if not (math.isfinite(fraction) and 0.0 < fraction <= 1.0):
        raise ValueError(
            f"The fraction of a medium change is in (0, 1], not {fraction}."
        )
    f = repr(float(fraction))
    values: dict[str, str] = {}
    for q in circuit.quantities(definition.perfusate, definition.released):
        concentration = f"[{ids.species(q, ids.RESERVOIR)}]"
        values[concentration] = f"(1 - {f})*{concentration} + {f}*{ids.fresh(q)}"
    return values


def perfusion_changes(
    definition: PerfusionDefinition,
    flow: Quantity,
    mass: Quantity,
    temperature: Quantity | None = None,
) -> dict[str, Quantity | float]:
    """Get the condition of a perfusion, the changes before the initialization.

    Args:
        definition: the definition of the model.
        flow: the total flow of the pump, `Q_pump` (a volume per time, 0 or
            more).
        mass: the mass of the organ, `M_<o>` (positive).
        temperature: the temperature, `T`, in degC, kelvin or any other unit
            of temperature; `None` keeps the temperature of the model.

    Returns:
        The flow in l/min and the mass in kg (quantities in the units of the
        model) and the temperature as its number in degC (the unit of `T`,
        see the module).

    Raises:
        ValueError: for a quantity of another dimension, a negative flow, a
            mass which is not positive or a temperature below the absolute
            zero.
    """
    q_pump = _converted(flow, "l/min", "flow")
    m_organ = _converted(mass, "kg", "mass")
    if not (math.isfinite(q_pump.magnitude) and q_pump.magnitude >= 0):
        raise ValueError(f"The flow of the pump is 0 or more, not {flow}.")
    if not (math.isfinite(m_organ.magnitude) and m_organ.magnitude > 0):
        raise ValueError(f"The mass of the organ is positive, not {mass}.")
    out: dict[str, Quantity | float] = {
        ids.Q_PUMP: q_pump,
        ids.mass(definition.organ.key): m_organ,
    }
    if temperature is not None:
        celsius = float(_converted(temperature, "degC", "temperature").magnitude)
        if not (math.isfinite(celsius) and celsius > -circuit.ABSOLUTE_ZERO):
            raise ValueError(
                f"The temperature is above the absolute zero, not {temperature}."
            )
        out[ids.TEMPERATURE] = celsius
    return out


def _converted(value: Quantity, unit: str, what: str) -> Quantity:
    """Convert a quantity into a unit of the model.

    Raises:
        ValueError: for a quantity of another dimension.
    """
    try:
        return value.to(unit)
    except Exception as err:  # pint's DimensionalityError, also of an offset unit
        raise ValueError(f"The {what} '{value}' is no quantity of {unit}.") from err


def _observables(definition: PerfusionDefinition) -> dict[str, tuple[str, str]]:
    """Get the label and the unit of every observable of a definition."""
    released = definition.released
    zones = [ids.zone_id(k) for k in range(1, definition.organ.zones + 1)]
    sites = [ids.RESERVOIR, ids.INFLOW, *zones, ids.OUTFLOW]
    means = [ids.PERIPORTAL, ids.PERICENTRAL]
    enzymes = [e for e in ids.ENZYMES if e in released]
    out: dict[str, tuple[str, str]] = {}
    for site in [*sites, *means]:
        where = site_label(site)
        for q in circuit.quantities(definition.perfusate, released):
            out[ids.species(q, site)] = (f"{circuit.quantity_label(q)} ({where})", _MM)
        out[ids.po2(site)] = (f"pO2 ({where})", "kPa")
        out[ids.pco2(site)] = (f"pCO2 ({where})", "kPa")
        out[ids.ph(site)] = (f"pH ({where})", "dimensionless")
        for e in enzymes:
            out[ids.activity(e, site)] = (
                f"{core.quantity_name(e)} activity ({where})",
                "µkat/l",
            )
    for site in [*zones, *means]:
        where = site_label(site)
        for q in ids.cell_quantities(released):
            out[ids.species(q, ids.cell_site(site))] = (
                f"{core.quantity_name(q)} (cells, {where})",
                _MM,
            )
        out[ids.alive(site)] = (f"living fraction ({where})", "dimensionless")
    for zone in zones:
        out[ids.necrosis_rate(zone)] = (
            f"rate of necrosis ({site_label(zone)})",
            "1/min",
        )
    out[ids.VO2] = ("O2 uptake", "mmole/min")
    if isinstance(definition.organ, LiverOrgan):
        out[ids.BILE] = ("bile", "µl")
        out[ids.BILE_RATE] = ("bile rate", "µl/min")
    if isinstance(definition.organ, KidneyOrgan):
        out[ids.URINE] = ("urine", "µl")
        out[ids.URINE_RATE] = ("urine rate", "µl/min")
        for q in (*ids.FILTERED, *definition.perfusate.tracers):
            out[ids.species(q, ids.URINE_SITE)] = (
                f"{circuit.quantity_label(q)} (urine)",
                _MM,
            )
    return out


def observable_labels(definition: PerfusionDefinition) -> dict[str, str]:
    """Get the labels of the observables of a perfusion model, by id.

    Every perfusate quantity, pO2, pCO2, pH and enzyme activity at every
    site (`reservoir`, `inflow`, every zone, `outflow`, the means
    `periportal` and `pericentral`), the quantities of the cells and the
    living fraction of every zone and mean, the rate of necrosis of a zone,
    the oxygen uptake, the bile of a liver and the urine of a kidney.

    Args:
        definition: the definition of the model.
    """
    return {sid: label for sid, (label, _) in _observables(definition).items()}


def observable_units(definition: PerfusionDefinition) -> dict[str, str]:
    """Get the units of the observables of a perfusion model, by id.

    The units of the observables of `observable_labels` as sbmlsim reads
    them (concentrations `mM`, partial pressures `kPa`, activities
    `µkat/l`, the bile and the urine `µl` and `µl/min`).

    Args:
        definition: the definition of the model.
    """
    return {sid: unit for sid, (_, unit) in _observables(definition).items()}


class PerfusionSimulationExperiment(PKPDSimulationExperiment):
    """The base of the simulation experiments of a machine perfusion model.

    The base of the experiments of a model sets the model object
    (`pkpd_model`) and the definition of the model (`definition`); a class
    with a definition gets the labels and units of its observables
    (`observable_labels`, `observable_units`) unless it sets its own. The
    condition of a simulation is `perfusion_changes`, its medium changes
    `medium_changes`, see the module.

    Attributes:
        definition: the definition of the model, from which its build is
            made.
    """

    definition: ClassVar[PerfusionDefinition]

    def __init_subclass__(cls, **kwargs: object) -> None:
        """Give a class with a definition the labels and units of its observables."""
        super().__init_subclass__(**kwargs)
        definition = cls.__dict__.get("definition")
        if isinstance(definition, PerfusionDefinition):
            if "labels" not in cls.__dict__:
                cls.labels = MappingProxyType(observable_labels(definition))
            if "units" not in cls.__dict__:
                cls.units = MappingProxyType(observable_units(definition))

    @classmethod
    def _definition(cls) -> PerfusionDefinition:
        """Get the definition of the model.

        Raises:
            TypeError: if no base of the experiment sets `definition`.
        """
        definition = getattr(cls, "definition", None)
        if not isinstance(definition, PerfusionDefinition):
            raise TypeError(
                f"{cls.__name__} has no perfusion definition: set `definition` in "
                "the base of the experiments of the model."
            )
        return definition

    @classmethod
    def medium_changes(
        cls, times: Sequence[float], fraction: float, time_unit: str = "min"
    ) -> list[Change]:
        """Get the medium changes of the reservoir at times.

        Args:
            times: the times of the changes.
            fraction: the fraction of the reservoir which each change
                replaces with fresh perfusate, in (0, 1].
            time_unit: the unit of the times.

        Returns:
            One `Change` per time, every perfusate species of the reservoir
            set towards the fresh perfusate (`medium_change`).

        Raises:
            TypeError: if the experiment has no definition.
            ValueError: for a fraction which is not in (0, 1].
        """
        values = medium_change(cls._definition(), fraction)
        return [Change(ureg.Quantity(float(t), time_unit), dict(values)) for t in times]

    @classmethod
    def perfusion_changes(
        cls,
        flow: Quantity,
        mass: Quantity,
        temperature: Quantity | None = None,
    ) -> dict[str, Quantity | float]:
        """Get the condition of a perfusion, the `preinit_changes` of a simulation.

        See the function `perfusion_changes`: the flow of the pump, the mass
        of the organ of the definition and the temperature.

        Raises:
            TypeError: if the experiment has no definition.
            ValueError: for a value which `perfusion_changes` refuses.
        """
        return perfusion_changes(cls._definition(), flow, mass, temperature)
