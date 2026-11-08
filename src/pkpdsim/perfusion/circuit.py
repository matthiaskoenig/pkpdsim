"""The circuit of a machine perfusion: the sites, the flow, the oxygenator and the temperature.

The pump drives the perfusate at the flow `Q_pump` (l/min) from the
reservoir `res` through the inflow tubing `in`, the zones of the organ
(`z1` to `zn`) and the outflow tubing `out` back into the reservoir. Every
site is a well-stirred volume, and the flow from a site to the next is the
reaction `FLOW_<q>_<site>_<next>` of every quantity of the perfusate (the
quantities of `ids.PERFUSATE` and the tracers) at `Q_pump * <q>_<site>`, so
the amount of the perfusate is conserved and a flow of 0 is defined (no rate
divides by the flow). The bicarbonate buffer acts at every site
(`core.buffer_objects`).

The membrane oxygenator exchanges the gases of the inflow tubing with its
sweep gas, first order: the reaction `OXY_in` sets the oxygen towards the
target of the perfusate, `k_oxy * V_in * (alpha_o2 * pO2_target - o2_in)`,
and the reaction `OXY_co2_in` the dissolved CO2 towards the partial pressure
of CO2 of the gas, `k_oxy_co2 * V_in * (alpha_co2 * pCO2_gas - co2_in)`, so
the CO2 which the cells produce and the bicarbonate buffer releases is
vented and the bicarbonate buffer of the perfusate is open, as with a gas
of fixed pCO2. Without a pCO2 of the gas in the definition it is that of
the fresh perfusate (`co2_fresh/alpha_co2`), so the fresh perfusate is at
its steady state in the oxygenator.

The temperature `T` (degC) is a constant which an experiment changes; the
rates of the cellular core scale with the Arrhenius factor `f_T = exp(Ea/R
* (1/(T_ref + 273.15) - 1/(T + 273.15)))`, 1 at the reference temperature,
and the solubility of oxygen is `alpha_o2 = alpha_o2_ref * exp(k_alpha_o2 *
(T - T_alpha_o2))`, both initial assignments of `T`, so the free oxygen of a
pO2 rises in the cold; the solubility of CO2 `alpha_co2` is a constant (its
value at 37 C). `T` takes effect as a change before the initialization of a
simulation (the initial assignments of `f_T`, `alpha_o2` and `o2_fresh` are
evaluated then); a change of `T` during a timecourse does not evaluate them
again, so it leaves the rates and the solubility at the earlier temperature.

The fresh perfusate is the constants `<q>_fresh` (mM), the start of every
perfusate species at every site (an initial assignment); the free protons
follow from its pH, the CO2 from the buffer equilibrium and the oxygen from
the target at the temperature. A medium change of an experiment sets the
reservoir towards them.
"""

import dataclasses
from collections.abc import Sequence

from sbmlutils.factory import Compartment, Parameter, Reaction, Species
from sbmlutils.factory.units import UnitType
from sbmlutils.metadata import SBO

from pkpdsim.pbpk.annotation import initial_assignment
from pkpdsim.pbpk.placement import value_notes, value_parameter
from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion import core, ids
from pkpdsim.perfusion.definition import Circuit, Perfusate, table_value
from pkpdsim.perfusion.units import PU

__all__ = [
    "ABSOLUTE_ZERO",
    "circuit_compartments",
    "circuit_objects",
    "composition",
    "converted",
    "flow_objects",
    "fresh_assignments",
    "fresh_co2_pressure",
    "oxygen_objects",
    "perfusate_parameters",
    "quantities",
    "temperature_objects",
    "tracer_species",
]

ABSOLUTE_ZERO = 273.15
"""The absolute temperature of 0 degC (K), the definition of the scale."""

_P = SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER

_SITE_NAMES = {
    ids.RESERVOIR: "reservoir",
    ids.INFLOW: "inflow tubing",
    ids.OUTFLOW: "outflow tubing",
}


def quantities(perfusate: Perfusate) -> tuple[str, ...]:
    """Get the quantities of the perfusate: `ids.PERFUSATE` and the tracers."""
    return (*ids.PERFUSATE, *perfusate.tracers)


def _name(quantity: str) -> str:
    """Get the name of a quantity of the perfusate, a tracer by its token."""
    try:
        return core.quantity_name(quantity)
    except KeyError:
        return f"tracer {quantity}"


def converted(value: Value, factor: float, unit: UnitType, text: str) -> Value:
    """Get a value in the unit of the model, the conversion in the notes.

    Args:
        value: the value of the definition.
        factor: the factor of the unit of the definition in the unit of the
            model.
        unit: the unit of the model.
        text: the unit of the definition, for the notes.

    Returns:
        The value with its source, in the unit of the model.

    Raises:
        ValueError: if the value has an uncertainty, which would keep the
            unit of the definition (`PerfusionDefinition.validate` refuses it).
    """
    if value.uncertainties:
        raise ValueError(
            f"A value in {text} with an uncertainty is not converted; "
            "the definition gives no uncertainty of it."
        )
    notes = f"{value.value} {text} in the definition."
    if value.notes:
        notes = f"{value.notes}\n\n{notes}"
    return dataclasses.replace(
        value, value=value.value * factor, unit=unit, notes=notes
    )


def temperature_objects(temperature: Value) -> core.CoreObjects:
    """Get the temperature and the Arrhenius factor `f_T` of the rates.

    Args:
        temperature: the temperature of the perfusion (degC).

    Returns:
        The constants `T`, `T_ref`, `R_gas` and `Ea` (the activation energy
        `f_Ea*R*T_ref` of the database, a constant, so a change of the
        reference keeps it) and `f_T` with its initial assignment.
    """
    t_ref = table_value("temperature", "reference")
    r_gas = table_value("temperature", "gas_constant")
    f_ea = table_value("temperature", "activation_energy_factor")
    ea = Value(
        f_ea.value * r_gas.value * (t_ref.value + ABSOLUTE_ZERO),
        PU.J_per_mole,
        f_ea.source,
        f"derived: Ea = f_Ea*R*T_ref with [temperature] activation_energy_factor "
        f"f_Ea = {f_ea.value} and the reference {t_ref.value} C of the perfusion "
        f"table; {f_ea.notes}",
    )
    t_abs = f"{ABSOLUTE_ZERO} {PU.degC.sid}"
    return core.CoreObjects(
        parameters=[
            value_parameter(
                ids.TEMPERATURE,
                temperature,
                PU.degC,
                "temperature",
                SBO.THERMODYNAMIC_TEMPERATURE,
            ),
            value_parameter(
                ids.TEMPERATURE_REF,
                t_ref,
                PU.degC,
                "reference temperature of the rates",
                SBO.THERMODYNAMIC_TEMPERATURE,
            ),
            value_parameter(
                ids.GAS_CONSTANT, r_gas, PU.J_per_mole_K, "molar gas constant", _P
            ),
            value_parameter(
                ids.ACTIVATION_ENERGY,
                ea,
                PU.J_per_mole,
                "activation energy of the rates",
                _P,
            ),
            Parameter(
                ids.F_T,
                1.0,
                U.dimensionless,
                constant=True,
                name="temperature factor of the rates (Arrhenius)",
                sboTerm=_P,
            ),
        ],
        assignments=[
            initial_assignment(
                ids.F_T,
                f"exp({ids.ACTIVATION_ENERGY}/{ids.GAS_CONSTANT} * "
                f"(1 dimensionless/({ids.TEMPERATURE_REF} + {t_abs}) - "
                f"1 dimensionless/({ids.TEMPERATURE} + {t_abs})))",
                U.dimensionless,
                "temperature factor of the rates (Arrhenius)",
            )
        ],
    )


def oxygen_objects() -> core.CoreObjects:
    """Get the solubility of the gases in the perfusate at the temperature.

    Returns:
        The constants of the database, `alpha_o2` with its initial
        assignment and the solubility of CO2 `alpha_co2`.
    """
    solubility = table_value("oxygen", "solubility")
    return core.CoreObjects(
        parameters=[
            value_parameter(
                ids.CO2_SOLUBILITY,
                table_value("acid_base", "co2_solubility"),
                PU.mM_per_kPa,
                "solubility of CO2",
                _P,
            ),
            value_parameter(
                ids.SOLUBILITY_REF,
                solubility,
                PU.mM_per_kPa,
                "solubility of oxygen at its reference temperature",
                _P,
            ),
            value_parameter(
                ids.SOLUBILITY_TEMPERATURE,
                table_value("oxygen", "temperature"),
                PU.degC,
                "reference temperature of the solubility of oxygen",
                SBO.THERMODYNAMIC_TEMPERATURE,
            ),
            value_parameter(
                ids.SOLUBILITY_COEFFICIENT,
                table_value("oxygen", "temperature_coefficient"),
                PU.per_K,
                "temperature coefficient of the solubility of oxygen",
                _P,
            ),
            Parameter(
                ids.SOLUBILITY,
                solubility.value,
                PU.mM_per_kPa,
                constant=True,
                name="solubility of oxygen at the temperature",
                sboTerm=_P,
            ),
        ],
        assignments=[
            initial_assignment(
                ids.SOLUBILITY,
                f"{ids.SOLUBILITY_REF} * exp({ids.SOLUBILITY_COEFFICIENT} * "
                f"({ids.TEMPERATURE} - {ids.SOLUBILITY_TEMPERATURE}))",
                PU.mM_per_kPa,
                "solubility of oxygen at the temperature",
            )
        ],
    )


def composition(perfusate: Perfusate) -> core.Composition:
    """Get the composition of the core for the initial values of the species.

    The initial assignments of the fresh perfusate (`fresh_assignments`) set
    the starts; the composition gives the initial values the same numbers at
    the reference temperature of the solubility.
    """
    alpha = table_value("oxygen", "solubility")
    target = perfusate.o2_target_kpa
    return core.Composition(
        glc=perfusate.glucose,
        o2=Value(alpha.value * target.value, U.mM, target.source, target.notes),
        lac=perfusate.lactate,
        hco3=perfusate.hco3,
        k=perfusate.k,
        ph=perfusate.ph,
    )


def perfusate_parameters(perfusate: Perfusate) -> core.CoreObjects:
    """Get the fresh perfusate, the constants `<q>_fresh`, and the oxygen target.

    Args:
        perfusate: the perfusate.

    Returns:
        The constants `pO2_target`, `pH_fresh` and `<q>_fresh` of glucose,
        lactate, bicarbonate, potassium, the enzymes (0) and the tracers (0);
        the free protons, the CO2 and the oxygen of the fresh perfusate as
        initial assignments.
    """
    given = {
        "glc": perfusate.glucose,
        "lac": perfusate.lactate,
        "hco3": perfusate.hco3,
        "k": perfusate.k,
    }
    none = Value(0.0, U.mM, "ModelChoice", "a fresh perfusate has none")
    for q in (*ids.ENZYMES, *perfusate.tracers):
        given[q] = none
    names = {q: _name(q) for q in quantities(perfusate)}
    out = core.CoreObjects(
        parameters=[
            value_parameter(
                ids.O2_TARGET,
                perfusate.o2_target_kpa,
                PU.kPa,
                "oxygen partial pressure of the oxygenator",
                _P,
            ),
            value_parameter(
                ids.fresh("pH"), perfusate.ph, U.dimensionless, "pH (fresh)", SBO.PH
            ),
            *(
                value_parameter(
                    ids.fresh(q),
                    v,
                    U.mM,
                    f"{names[q]} (fresh)",
                    SBO.MOLAR_CONCENTRATION_OF_AN_ENTITY,
                )
                for q, v in given.items()
            ),
        ]
    )
    derived = {
        "h": f"power(10 dimensionless, 3 dimensionless - {ids.fresh('pH')}) * 1 mM",
        "co2": f"{ids.fresh('h')} * {ids.fresh('hco3')}/{ids.BUFFER_KA}",
        "o2": f"{ids.SOLUBILITY} * {ids.O2_TARGET}",
    }
    for q, formula in derived.items():
        name = f"{names[q]} (fresh)"
        out.parameters.append(
            Parameter(
                ids.fresh(q),
                0.0,
                U.mM,
                constant=True,
                name=name,
                sboTerm=SBO.MOLAR_CONCENTRATION_OF_AN_ENTITY,
            )
        )
        out.assignments.append(initial_assignment(ids.fresh(q), formula, U.mM, name))
    return out


def fresh_assignments(site: str, perfusate: Perfusate) -> core.CoreObjects:
    """Get the start of the perfusate species of a site, the fresh perfusate.

    The CO2 keeps the initial assignment of the buffer equilibrium of its
    site (`core.perfusate_objects`), which is the fresh CO2.
    """
    return core.CoreObjects(
        assignments=[
            initial_assignment(
                ids.species(q, site),
                ids.fresh(q),
                U.mM,
                f"{_name(q)} (perfusate {site}) from the fresh perfusate",
            )
            for q in quantities(perfusate)
            if q != "co2"
        ]
    )


def tracer_species(tracer: str, site: str, compartment: str) -> Species:
    """Get the species of an inert tracer of the perfusate at a site."""
    return Species(
        ids.species(tracer, site),
        compartment,
        initialConcentration=0.0,
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        name=f"tracer {tracer} (perfusate {site})",
        sboTerm=SBO.SIMPLE_CHEMICAL,
    )


def circuit_compartments(circuit: Circuit) -> list[Compartment]:
    """Get the compartments of the reservoir and the tubing (l)."""
    volumes = {
        ids.RESERVOIR: circuit.reservoir,
        ids.INFLOW: circuit.tubing_in,
        ids.OUTFLOW: circuit.tubing_out,
    }
    out = []
    for site, value in volumes.items():
        v = converted(value, 1e-3, U.litre, "ml")
        out.append(
            Compartment(
                ids.volume(site),
                v.value,
                unit=U.litre,
                name=f"perfusate in the {_SITE_NAMES[site]}",
                sboTerm=SBO.PHYSICAL_COMPARTMENT,
                notes=value_notes(v),
            )
        )
    return out


def fresh_co2_pressure(perfusate: Perfusate) -> Value:
    """Get the partial pressure of CO2 of the fresh perfusate (kPa).

    The CO2 of the buffer equilibrium at the pH and the bicarbonate of the
    fresh perfusate, `10^(3 - pH) mM * HCO3/Ka_hco3`, over the solubility of
    CO2; the default pCO2 of the sweep gas, at which the fresh perfusate
    stays at its pH.
    """
    ka = core.Buffer.default().ka.value
    alpha = table_value("acid_base", "co2_solubility").value
    co2 = 10 ** (3 - perfusate.ph.value) * perfusate.hco3.value / ka
    return Value(
        co2 / alpha,
        PU.kPa,
        "ModelChoice",
        "derived: the pCO2 of the fresh perfusate at the equilibrium of its "
        "bicarbonate buffer, CO2 = 10^(3 - pH_fresh) mM * HCO3_fresh/Ka_hco3, "
        "over the solubility of CO2, so the fresh perfusate is at its steady "
        "state in the oxygenator; the sweep gas of the study protocol (about "
        "5 % CO2 is usual) is a TODO",
    )


def circuit_objects(circuit: Circuit, perfusate: Perfusate) -> core.CoreObjects:
    """Get the perfusate of the reservoir and the tubing, the pump and the oxygenator.

    Args:
        circuit: the circuit.
        perfusate: the perfusate.

    Returns:
        The perfusate species, the buffer and the start of the fresh
        perfusate at the sites of the circuit, the constants `Q_pump` (l/min),
        `k_oxy`, `k_oxy_co2` and `pCO2_gas`, the reactions of the oxygenator
        (oxygen and CO2).
    """
    out = core.CoreObjects(
        parameters=[
            value_parameter(
                ids.Q_PUMP,
                converted(circuit.flow, 1e-3, U.l_per_min, "ml/min"),
                U.l_per_min,
                "total flow of the pump",
                _P,
            ),
            value_parameter(
                ids.OXYGENATOR_K,
                circuit.oxygenator_k,
                U.per_min,
                "rate constant of the oxygenator (oxygen)",
                SBO.KINETIC_CONSTANT,
            ),
            value_parameter(
                ids.OXYGENATOR_K_CO2,
                _co2_rate(circuit),
                U.per_min,
                "rate constant of the oxygenator (CO2)",
                SBO.KINETIC_CONSTANT,
            ),
            value_parameter(
                ids.CO2_GAS,
                circuit.pco2_gas or fresh_co2_pressure(perfusate),
                PU.kPa,
                "partial pressure of CO2 of the sweep gas",
                SBO.PRESSURE,
            ),
        ]
    )
    start = composition(perfusate)
    for site in (ids.RESERVOIR, ids.INFLOW, ids.OUTFLOW):
        compartment = ids.volume(site)
        out += core.perfusate_objects(site, compartment, start)
        out += core.buffer_objects(site, compartment)
        out.species += [tracer_species(t, site, compartment) for t in perfusate.tracers]
        out += fresh_assignments(site, perfusate)
    o2, v_in = ids.species("o2", ids.INFLOW), ids.volume(ids.INFLOW)
    out.reactions.append(
        Reaction(
            ids.OXYGENATOR,
            f"<-> {o2}",
            formula=(
                f"{ids.OXYGENATOR_K} * {v_in} * "
                f"({ids.SOLUBILITY} * {ids.O2_TARGET} - {o2})",
                U.mmole_per_min,
            ),
            compartment=v_in,
            name="oxygenator",
            sboTerm=SBO.TRANSPORT_REACTION,
            notes="Sets the oxygen of the inflow tubing towards the target, "
            "first order; the gas side is the parameter pO2_target.",
        )
    )
    co2 = ids.species("co2", ids.INFLOW)
    out.reactions.append(
        Reaction(
            ids.OXYGENATOR_CO2,
            f"{co2} <->",
            formula=(
                f"{ids.OXYGENATOR_K_CO2} * {v_in} * "
                f"({co2} - {ids.CO2_SOLUBILITY} * {ids.CO2_GAS})",
                U.mmole_per_min,
            ),
            compartment=v_in,
            name="oxygenator (CO2)",
            sboTerm=SBO.TRANSPORT_REACTION,
            notes="Vents the CO2 of the inflow tubing towards the pCO2 of the "
            "sweep gas, first order; the gas side is the parameter pCO2_gas.",
        )
    )
    return out


def _co2_rate(circuit: Circuit) -> Value:
    """Get the rate constant of the CO2 exchange, that of the oxygen without one."""
    if circuit.oxygenator_k_co2 is not None:
        return circuit.oxygenator_k_co2
    k = circuit.oxygenator_k
    notes = "the rate constant of the oxygen (oxygenator_k)"
    return dataclasses.replace(k, notes=f"{notes}; {k.notes}" if k.notes else notes)


def flow_objects(sites: Sequence[str], perfusate: Perfusate) -> core.CoreObjects:
    """Get the flow of the perfusate through the sites of the circuit.

    Args:
        sites: the sites in the order of the flow, the last flows into the
            first (the reservoir, the inflow, the zones, the outflow).
        perfusate: the perfusate.

    Returns:
        The reaction `FLOW_<q>_<site>_<next>` of every quantity and site,
        `Q_pump * <q>_<site>`.
    """
    out = core.CoreObjects()
    for source, target in zip(sites, [*sites[1:], sites[0]], strict=True):
        for q in quantities(perfusate):
            a, b = ids.species(q, source), ids.species(q, target)
            out.reactions.append(
                Reaction(
                    ids.flow(q, source, target),
                    f"{a} -> {b}",
                    formula=(f"{ids.Q_PUMP} * {a}", U.mmole_per_min),
                    compartment=ids.volume(source),
                    name=f"flow of {q} from {source} to {target}",
                    sboTerm=SBO.TRANSPORT_REACTION,
                )
            )
    return out
