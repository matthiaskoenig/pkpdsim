# Machine perfusion

Ex vivo machine perfusion keeps a donor organ alive in a circuit: a pump drives an oxygenated perfusate from a reservoir through the organ and back. The organ consumes oxygen and substrates, makes bile or urine and, when it is damaged, releases its enzymes into the perfusate. `pkpdsim.perfusion` builds the SBML of such a circuit with one mounted organ from a definition, as the [PBPK library](pbpk/index.md) builds a body: the reservoir, the tubing, the pump, the oxygenator, the temperature, the perfusate and the organ, a chain of zones with the cellular core. The organs are the liver with its bile, the kidney with its filtration and urine and a generic organ without excretion, so another organ is a definition and not a new model (pkpdsim#146).

## A model

A `PerfusionDefinition` has the id and name of the model, the species of the organ, the organ, the circuit, the fresh perfusate, the temperature and the constants of the cellular core, and the creators and notes of the model (the notes, e.g. the terms of use of a model package, follow the description of the circuit). The defaults of `Circuit`, `Perfusate`, `LiverOrgan`, `KidneyOrgan` and `CoreParameters` are the numbers of the parameter database (`tables/perfusion.toml`, see [Parameter database](database.md)), in the units an experimenter gives them (ml, ml/min, g, kPa, degC); `build` validates the definition and writes `<sid>.xml`, SBML L3V2 with units, every element named and annotated with SBO. The rat liver of the study of pkdb_models, in four zones from periportal to pericentral:

```python
from pathlib import Path

from pkpdsim.pbpk.values import Value
from pkpdsim.perfusion import (
    Circuit,
    CoreParameters,
    LiverOrgan,
    Perfusate,
    PerfusionDefinition,
    build,
)
from pkpdsim.perfusion.units import PU
from pkpdsim.species import Species

definition = PerfusionDefinition(
    sid="rat_liver_nmp",
    name="rat liver in a normothermic machine perfusion",
    species=Species.RAT,
    organ=LiverOrgan.default(zones=4),
    circuit=Circuit.default(),
    perfusate=Perfusate.default(),
    temperature=Value(37.0, PU.degC, "ModelChoice", "a normothermic perfusion"),
    core=CoreParameters.default(),
)
assert definition.validate() == []
path = build(definition, Path("models"))
print(path.name)
```

```text
rat_liver_nmp.xml
```

<!-- generated: figures-perfusion-1 (scripts/visual/docs_figures.py), do not edit -->
The content of the model `rat_liver_nmp`: the organ in its circuit (A) and the schematic of the circuit (B), hover them for details ([figures](pbpk/figures.md)).

<iframe class="pkpd-figure" src="img/figures/models/perfusion-1-anatomy.html" title="The anatomical figure of the model rat_liver_nmp" style="aspect-ratio:1044 / 608" loading="lazy"></iframe>

<iframe class="pkpd-figure" src="img/figures/models/perfusion-1-schematic.html" title="The schematic of the model rat_liver_nmp" style="aspect-ratio:1273 / 634" loading="lazy"></iframe>
<!-- end generated -->

`create_archive(definition, output_dir, model=, copyright=)` (`pkpdsim.perfusion.archive`) writes the model with its ODEs, its antimony, the two figures and the COMBINE archive, as the factory of a model package does. The archive names the anatomical figure A `<sid>_anatomy.svg` (key `<sid>_anatomy_svg`; the PBPK library names its anatomical figure `<sid>_body.svg`) and the schematic B `<sid>.svg`. The legend of an interactive figure credits pkpdsim, the machine perfusion library with this page (`CREDIT` and `ICON_CREDIT` of `pkpdsim.perfusion.figure`) and the anatomograms of the Expression Atlas, whose organ shapes (CC BY 4.0) the figures draw.

## The circuit

The pump drives the perfusate at the flow `Q_pump` (l/min, the total flow of the pump) from the reservoir `res` through the inflow tubing `in`, the zones of the organ `z1` to `zn` (their vascular space is the site `vas_<z>` of the zonation engine, `<q>_vas_<z>`) and the outflow tubing `out` back into the reservoir. Every site is a well-stirred volume (`V_res`, `V_in`, `V_out`, `V_vas_<z>`), and every quantity of the perfusate flows from a site to the next at `Q_pump` times its concentration (in the liver the flow between the zones is `Q_li`, the flow through the liver, which the pump sets: the zonation engine chains the zones), so the amount of the perfusate is conserved and a flow of 0 is defined (no rate divides by the flow). The membrane oxygenator in the inflow tubing exchanges the gases with its sweep gas, first order: the oxygen towards the target of the perfusate, the CO2 towards the pCO2 of the gas, so the CO2 which the cells produce is vented and the bicarbonate buffer of the perfusate is open. Without a pCO2 of the gas in the definition it is that of the fresh perfusate, the initial assignment `pCO2_gas = co2_fresh/alpha_co2`, at whose buffer equilibrium the oxygenator then leaves it, also after a change of the fresh perfusate before the initialization.

The equations below are those of a toy (`<q>` a quantity of the perfusate, `<z>` a zone, `<m>` the submodel of a zone, `M01` for `z1` (the cell elements which the liver does not expose keep `<m>__<id>` in the flat model, the states it observes have the ids of the zone), `<o>` the organ, `li` or `ki`); `tests/docs/test_perfusion_equations.py` builds a liver and a kidney of two zones and checks every equation of this page against them.

```equations title="the flow through the circuit (a liver of two zones)"
FLOW_<q>_res_in: <q>_res -> <q>_in = Q_pump * <q>_res
FLOW_<q>_in_vas_z1: <q>_in -> <q>_vas_z1 = Q_pump * <q>_in
FLOW_<q>_vas_z1_vas_z2: <q>_vas_z1 -> <q>_vas_z2 = Q_li * <q>_vas_z1
FLOW_<q>_vas_z2_out: <q>_vas_z2 -> <q>_out = Q_pump * <q>_vas_z2
FLOW_<q>_out_res: <q>_out -> <q>_res = Q_pump * <q>_out
```

```equations title="the oxygenator"
OXY_in: -> o2_in = k_oxy * V_in * (alpha_o2 * pO2_target - o2_in)
OXY_co2_in: co2_in -> = k_oxy_co2 * V_in * (co2_in - alpha_co2 * pCO2_gas)
pCO2_gas = co2_fresh / alpha_co2
```

The fresh perfusate is the constants `<q>_fresh` (mM), the start of every perfusate species at every site: glucose, lactate, bicarbonate, potassium and the enzymes (none) as given, the free protons from its pH (between 6 and 8, `PH_RANGE` of `pkpdsim.perfusion.definition`, `validate` refuses another), the dissolved CO2 from the buffer equilibrium and the free oxygen from the target at the temperature. The bicarbonate buffer acts at every site, `H+ + HCO3- <-> CO2` by mass action (the hydration of CO2 sets the time), and the pH follows from the free protons.

```equations title="the fresh perfusate and the buffer"
h_fresh = 10^(3 - pH_fresh) * 1
co2_fresh = h_fresh * hco3_fresh / Ka_hco3
o2_fresh = alpha_o2 * pO2_target
glc_res = glc_fresh
BUF_res: h_res + hco3_res -> co2_res = V_res * k_co2 * (h_res * hco3_res / Ka_hco3 - co2_res)
pH_res = 3 - log10(h_res / 1)
```

## The zonated organ

The organ is a chain of `n` zones in series with the flow, the first zone first (periportal in the liver, as the zonated liver of the PBPK library): each zone a vascular space through which the perfusate flows and its cells with the cellular core, so that oxygen, glucose, lactate, ATP and the damage form gradients along the sinusoid; `n = 1` is the single compartment. An organ with the token `<o>` has the mass `M_<o>` (kg), the density `rho_<o>`, the volume `V_<o>` and the vascular and cell fractions `f_vas_<o>` and `f_cell_<o>`; the zones share the volume equally. The zones of the liver are those of the zonation engine (`pkpdsim.zonation`, in the well-stirred geometry with the fractions of the organ), so the liver (the last two equations) has the engine's ids `f_vas`, `f_cell`, `f_par` (1) and `V_vas_<z>`, `V_cell_<z>` and `V_li = M_li/rho_li`; the cells of a zone are a submodel of the cell model `<sid>_core`, `M<kk>__<id>` in the flat model, and the constants of the core are shared by the zones. The volumes are initial assignments of the mass, so a change of the mass before a simulation (the liver weight of an animal) sets the zones.

```equations title="the volumes of the zones (an organ of two zones)"
V_<o> = M_<o> / rho_<o>
V_vas_<z> = f_vas_<o> * V_<o> / 2
V_cell_<z> = f_cell_<o> * V_<o> / 2
V_vas_<z> = f_vas * f_par * M_li / rho_li / 2
V_cell_<z> = f_cell * f_par * M_li / rho_li / 2
```

A profile of a process (`LiverOrgan(profiles={"OXPHOS": [2, 1]})`) weighs it per zone with the factor `p_<process>_<z> = n*p/sum(p)` in its rate, so the factors sum to `n` and the total capacity of the organ is that of one zone of the whole organ; a process without a profile has none, a uniform profile is 1 in every zone.

## The cellular core

The cells of a zone take up glucose and oxygen, store glycogen, run glycolysis, the lactate dehydrogenase, the TCA cycle and the oxidative phosphorylation, use ATP, produce reactive oxygen species (ROS) at low oxygen and detoxify them with glutathione, regenerated with the NADPH of the pentose phosphate pathway, export lactate with a proton, die below an ATP threshold or in an acidic perfusate (necrosis) and release their enzymes and potassium into the vascular space, and the liver secretes bile. The processes and values are ported from the liver model of `pkdb_models.misc.machine_perfusion` (from the PVL liver model) with its defects fixed, the rates per litre of cells; the module docstring of `pkpdsim.perfusion.core` gives every rate law with the reason for it. The constants are shared by the zones, so a fit has one value per constant. The starts of the cells are shared as well: the constants `<q>_start` (mM, `CoreParameters.starts`) of every quantity of the cells but the released ones, whose start is their content (`ALT_content`, `K_content`, ...), and every zone starts from them by an initial assignment, so a covariate on the starting ATP is one change of `atp_start`. A constant in another unit than that of the model, a quantity of the cells without a start and a start of no quantity are refused by `validate`, and an organ without bile has no bile constants, an organ which does not release an enzyme no content of it.

```equations title="the starts of the cells"
atp_<z>_cell = atp_start
alt_<z>_cell = ALT_content
```

Every biological rate scales with the temperature factor `f_T` (the liver has its own `f_T`, which its bile secretion uses, and every cell model has one, `<m>__f_T`, which the equations below use; both are the same function of the temperature), the living fraction `f_alive_<z>` and the cell volume, only the living cells metabolise; a Hill term reads the absolute value of its species, so it stays defined for the small negative values of the integrator, and no rate divides by a species which can reach 0 (the inhibition of the glycolysis by ATP is defined at ATP 0).

```equations title="transport between the perfusate and the cells"
GLCIM_<z>: glc_vas_<z> -> glc_<z>_cell = (<m>__f_T * f_alive_<z> * V_cell_<z> * GLCIM_Vmax / GLCIM_Km_glc) * (glc_vas_<z> - glc_<z>_cell) / (1 + glc_vas_<z> / GLCIM_Km_glc + glc_<z>_cell / GLCIM_Km_glc)
O2IM_<z>: o2_vas_<z> -> o2_<z>_cell = (<m>__f_T * f_alive_<z> * V_cell_<z> * O2IM_Vmax / O2IM_Km_o2) * (o2_vas_<z> - o2_<z>_cell) / (1 + o2_vas_<z> / O2IM_Km_o2 + o2_<z>_cell / O2IM_Km_o2)
LACEX_<z>: lac_<z>_cell -> lac_vas_<z> + h_vas_<z> = (<m>__f_T * f_alive_<z> * V_cell_<z> * LACEX_Vmax / LACEX_Km_lac) * (lac_<z>_cell - lac_vas_<z> * h_vas_<z> / LACEX_h_cell) / (1 + lac_<z>_cell / LACEX_Km_lac + lac_vas_<z> / LACEX_Km_lac)
```

```equations title="the energy metabolism"
GLGP_<z>: glyc_<z>_cell -> glc_<z>_cell = <m>__f_T * f_alive_<z> * V_cell_<z> * GLGP_Vmax * glyc_<z>_cell / (GLGP_Km_glyc + glyc_<z>_cell)
GLC2PYR_<z>: glc_<z>_cell + adp_<z>_cell + nad_<z>_cell -> pyr_<z>_cell + atp_<z>_cell + nadh_<z>_cell = <m>__f_T * f_alive_<z> * V_cell_<z> * GLC2PYR_Vmax * glc_<z>_cell / (GLC2PYR_Km_glc + glc_<z>_cell) * adp_<z>_cell / (GLC2PYR_Km_adp + adp_<z>_cell) * nad_<z>_cell / (GLC2PYR_Km_nad + nad_<z>_cell) * 1 / (1 + (abs(atp_<z>_cell) / GLC2PYR_Ki_atp)^GLC2PYR_n_atp)
LDH_<z>: pyr_<z>_cell + nadh_<z>_cell -> lac_<z>_cell + nad_<z>_cell = (<m>__f_T * f_alive_<z> * V_cell_<z> * LDH_Vmax / (LDH_Km_pyr * LDH_Km_nadh)) * (pyr_<z>_cell * nadh_<z>_cell - lac_<z>_cell * nad_<z>_cell / LDH_Keq) / ((1 + pyr_<z>_cell / LDH_Km_pyr + lac_<z>_cell / LDH_Km_lac) * (1 + nadh_<z>_cell / LDH_Km_nadh + nad_<z>_cell / LDH_Km_nad))
TCA_<z>: pyr_<z>_cell + nad_<z>_cell -> nadh_<z>_cell + co2_vas_<z> = <m>__f_T * f_alive_<z> * V_cell_<z> * TCA_Vmax * pyr_<z>_cell / (TCA_Km_pyr + pyr_<z>_cell) * nad_<z>_cell / (TCA_Km_nad + nad_<z>_cell)
OXPHOS_<z>: nadh_<z>_cell + o2_<z>_cell + adp_<z>_cell -> nad_<z>_cell + atp_<z>_cell = <m>__f_T * f_alive_<z> * V_cell_<z> * OXPHOS_Vmax * nadh_<z>_cell / (OXPHOS_Km_nadh + nadh_<z>_cell) * o2_<z>_cell / (OXPHOS_Km_o2 + o2_<z>_cell) * adp_<z>_cell / (OXPHOS_Km_adp + adp_<z>_cell)
ATPU_<z>: atp_<z>_cell -> adp_<z>_cell = <m>__f_T * f_alive_<z> * V_cell_<z> * ATPU_Vmax * atp_<z>_cell / (ATPU_Km_atp + atp_<z>_cell)
```

The stoichiometry of a reaction is in its SBML (the equations list the species): the glycolysis `glc + 2 ADP + 2 NAD -> 2 pyr + 2 ATP + 2 NADH`, the TCA cycle `pyr + 5 NAD -> 5 NADH + 3 CO2` (the CO2 into the perfusate), the oxidative phosphorylation `NADH + 0.5 O2 + 3 ADP -> NAD + 3 ATP`.

```equations title="the redox"
ROSP_<z>: nadh_<z>_cell + o2_<z>_cell -> nad_<z>_cell + ros_<z>_cell = <m>__f_T * f_alive_<z> * V_cell_<z> * ROSP_Vmax * nadh_<z>_cell / (ROSP_Km_nadh + nadh_<z>_cell) * o2_<z>_cell / (ROSP_Km_o2 + o2_<z>_cell) * 1 / (1 + (abs(o2_<z>_cell) / ROSP_Ki_o2)^ROSP_n_o2)
ROSDETOX_<z>: ros_<z>_cell + gsh_<z>_cell -> gssg_<z>_cell = <m>__f_T * f_alive_<z> * V_cell_<z> * ROSDETOX_k * gsh_<z>_cell * ros_<z>_cell / (ROSDETOX_Km_ros + ros_<z>_cell)
GR_<z>: gssg_<z>_cell + nadph_<z>_cell -> gsh_<z>_cell + nadp_<z>_cell = <m>__f_T * f_alive_<z> * V_cell_<z> * GR_k * gssg_<z>_cell * nadph_<z>_cell / (GR_Km_nadph + nadph_<z>_cell)
PPP_<z>: glc_<z>_cell + nadp_<z>_cell -> nadph_<z>_cell + co2_vas_<z> = <m>__f_T * f_alive_<z> * V_cell_<z> * PPP_Vmax * glc_<z>_cell / (PPP_Km_glc + glc_<z>_cell) * nadp_<z>_cell / (PPP_Km_nadp + nadp_<z>_cell)
```

The cells die at the rate `v_nec_<z>`, either drive alone gives the full rate: the energy failure, a smooth threshold of ATP, and the acidosis of the perfusate, a Hill term of the free protons. The acidosis drives the damage and does not inhibit the energy metabolism, since the medium changes of the study protect the liver far beyond the dilution of the enzymes. A dying cell moves its enzymes and potassium into the dead cells of the zone (`<q>_<z>_dead`), which release them first order into the vascular space, so the totals are conserved; the release is the leak of a ruptured membrane and does not scale with the temperature. The liver releases ALT, AST, LDH and potassium, the kidney and a generic organ LDH and potassium.

```equations title="necrosis and release"
v_nec_<z> = NEC_k * <m>__f_T * (1 - (abs(atp_<z>_cell) / NEC_atp)^NEC_n_atp / ((1 + (abs(atp_<z>_cell) / NEC_atp)^NEC_n_atp) * (1 + (abs(h_vas_<z>) / NEC_h)^NEC_n_h)))
d(f_alive_<z>)/dt = -v_nec_<z> * f_alive_<z>
ALTD_<z>: alt_<z>_cell -> alt_<z>_dead = v_nec_<z> * alt_<z>_cell * V_cell_<z>
ALTR_<z>: alt_<z>_dead -> alt_vas_<z> = PROTEINR_k * alt_<z>_dead * V_cell_<z>
KD_<z>: k_<z>_cell -> k_<z>_dead = v_nec_<z> * k_<z>_cell * V_cell_<z>
KR_<z>: k_<z>_dead -> k_vas_<z> = KR_k * k_<z>_dead * V_cell_<z>
```

```equations title="bile (a liver of two zones)"
d(Vbile_<z>)/dt = <m>__f_T * f_alive_<z> * V_cell_<z> * BILE_Vmax * atp_<z>_cell / (BILE_Km_atp + atp_<z>_cell)
bile = 1000000 * (Vbile_z1 + Vbile_z2)
```

### An extension of a package

A package adds a mechanism which the library does not have with a `PerfusionExtension` (`build(definition, output_dir, extension)`, `create_archive(..., extension=...)`): `cell` extends the cell model of the zones, `composed` the composed model, and `shared` names the elements of the cell model which the liver holds once.

The cells of a mechanism are written into the cell model once and not per zone, so the engine exposes them as it does the cells of the library: a species `<q>_cell` of the cell model is `<q>_<z>_cell` in the flat model, a reaction `R` is `R_<z>`, and every other element of the cell model stays `M<kk>__<id>`. A constant of the mechanism belongs in `shared`, a port of the cell model and an element of the composed model, so the fit has one value of it instead of one per zone; `compose` refuses a shared id which the composed model does not hold. What the circuit or the zones hold (a quantity of the perfusate, an observer over the zones) is written into the composed model, after the build has placed the circuit, the zones and the observers.

The reperfusion injury of the machine perfusion model of pkdb_models (`pkdb_models.misc.machine_perfusion.models.reperfusion`) is such an extension: the hypoxic and the dying cells release a signal which recirculates with the perfusate and drives the ROS, the loss of nuclear HMGB1 and the necrosis.

## The temperature

The temperature `T` (degC) is a constant of the model; the rates of the cells scale with the Arrhenius factor `f_T`, 1 at the reference temperature `T_ref` (37 degC), and the solubility of oxygen rises in the cold, both initial assignments of `T`, so a pO2 holds more free oxygen in a hypothermic perfusion. SBML has no degree Celsius: `T` is in `degC`, a unit of the size of a kelvin, and the formulas write out the absolute temperature. The temperature takes effect as a change before the initialization (`perfusion_changes` below); a change of `T` during a timecourse does not evaluate the initial assignments again.

```equations title="the temperature"
f_T = exp((Ea / R_gas) * (1 / (T_ref + 273.15) - 1 / (T + 273.15)))
alpha_o2 = alpha_o2_ref * exp(k_alpha_o2 * (T - T_alpha_o2))
```

## The experiments and the medium change

`PerfusionSimulationExperiment` (`pkpdsim.perfusion.experiment`) is the base of the simulation experiments of a perfusion model, a `PKPDSimulationExperiment` with the definition of the model as `definition`; a class with a definition gets the labels and units of its observables (`observable_labels`, `observable_units`) unless it sets its own. The condition of a simulation is a change before the initialization, `perfusion_changes(flow, mass, temperature=None)`: the pump flow `Q_pump`, the mass of the organ `M_<o>` and the temperature `T` as its number in degC. Never give `T` as a quantity: sbmlsim converts a quantity of a temperature to kelvin, `Q(37, "degC")` to 310.15, while `T` counts kelvin steps from 0 degC, so the model would run at 310 degC (`f_T` about 276); the initialization of a `PerfusionSimulationExperiment` refuses a simulation which sets `T` to a quantity (before the initialization, in its presimulation, at a time or in a dimension of a scan). The groups of a study differ by these and by their medium changes, the parameters of the model are shared.

A medium change replaces the fraction `f` of the reservoir with fresh perfusate: every perfusate species of the reservoir, the quantities of the perfusate of the organ and the tracers, is set to `(1 - f)*[X_res] + f*X_fresh` with the fresh perfusate of the model. It is a `Change` of sbmlsim at its time (`medium_changes(times, fraction, time_unit="min")`, one per time), so the model has no events and the times are a setting of the experiment. The change keeps the volume of the reservoir, as much fresh perfusate comes in as goes out; the reservoir of a kidney, which loses the urine, keeps the volume it has at the change, the urine is not made up. The tubing and the organ are not changed.

sbmlsim evaluates the values of a change with the state at its time and assigns them at once. An output time at a change has two rows, the state before and the state after the change, and a fit compares a data point at the time of a change with the state before it (`data_at_change` of `pkpdsim.fit.default_fit_settings`), a sample taken before the medium is changed; `tests/perfusion/test_experiment.py` pins this order.

```python
import numpy as np
from sbmlsim import Q
from sbmlsim.simulation import Change, Simulation
from sbmlsim.simulator import Simulator

from pkpdsim.perfusion import medium_change, perfusion_changes

values = medium_change(definition, 0.5)
print(values["[glc_res]"])

simulator = Simulator(n_workers=1)
model = simulator.load(path)
model.selections = ["time", "[glc_res]", "pO2_out"]
condition = perfusion_changes(
    definition, Q(3.0, "ml/min"), Q(10.0, "g"), Q(37.0, "degC")
)
for n, text in ((0, "no change"), (1, "one change"), (2, "two changes")):
    simulation = Simulation(
        end=6.0,
        time_unit="hr",
        preinit_changes=condition,
        changes=[Change(Q(t, "hr"), values) for t in (2.0, 4.0)[:n]],
        times=[0.0, 6.0],
    )
    result = simulator.simulate(model, simulation)
    glc = np.asarray(result["[glc_res]"]).ravel()[-1]
    po2 = np.asarray(result["pO2_out"]).ravel()[-1]
    print(f"{text}: glucose {glc:.3f} mM, outflow pO2 {po2:.1f} kPa at 6 h")
```

```text
(1 - 0.5)*[glc_res] + 0.5*glc_fresh
no change: glucose 6.921 mM, outflow pO2 46.1 kPa at 6 h
one change: glucose 6.858 mM, outflow pO2 46.2 kPa at 6 h
two changes: glucose 6.755 mM, outflow pO2 46.3 kPa at 6 h
```

The experiments of a study subclass the base with the model object and the definition of the package, each group a simulation:

```python
# not executed
from sbmlsim import Q
from sbmlsim.simulation import Scan, Simulation

from pkpdsim.perfusion.experiment import PerfusionSimulationExperiment

from .. import DEFINITION, MODEL


class FlowRate(PerfusionSimulationExperiment):
    pkpd_model = MODEL
    definition = DEFINITION

    def simulations(self) -> dict[str, Simulation | Scan]:
        return {
            "q3_mc2": Simulation(
                end=6.0,
                time_unit="hr",
                preinit_changes=self.perfusion_changes(Q(3.0, "ml/min"), Q(10.0, "g")),
                changes=self.medium_changes([2.0, 4.0], 0.5, time_unit="hr"),
                steps=360,
            ),
        }
```

The fit mappings of a perfusion describe their data with the vocabulary of `pkpdsim.metadata` shared by all models: the tissue `Tissue.PERFUSATE`, the route `Route.EX_VIVO` (an organ outside the body), the dosing `Dosing.PERFUSION` (no dose, the organ is kept in the perfusate of the circuit) and the sampling site of the circuit, `SamplingSite.RESERVOIR` (the default of the perfusate), `INFLOW` or `OUTFLOW`, the sites `res`, `in` and `out` of the observers below. A site of the circuit on another tissue and another site on the perfusate are a `ValueError`.

## The observers and the sites

Every quantity can be read at every site: the reservoir `res` (where the study samples the perfusate), the inflow `in` after the oxygenator, every zone `z<k>`, the outflow `out` after the organ, and the periportal and pericentral means `pp` and `pc`, the first and the last third of the zones, each zone weighted by its overlap with the third (with one zone both means are the zone). The cells of a zone are the site `<z>_cell` (the species `<q>_<z>_cell`) (`pp_cell`, `pc_cell` for the means). A species is an amount in roadrunner and sbmlsim and its concentration is selected as `[<q>_<site>]` (the 10 ml of the inflow tubing hold 0.065 mmole of glucose, `glc_in`, at 6.5 mM, `[glc_in]`); the means and the other observers are assignment rules, selected by their id. The labels and units of `observable_labels` and `observable_units` are keyed by these selections.

| observer | unit | what |
|---|---|---|
| `[<q>_<site>]` | mM | a quantity of the perfusate (glucose, oxygen, lactate, CO2, bicarbonate, protons, potassium, the released enzymes, the tracers) at the reservoir, the inflow, a zone or the outflow, a species |
| `<q>_pp`, `<q>_pc` | mM | the periportal and pericentral mean of a quantity of the perfusate, an assignment rule |
| `[<q>_<z>_cell]`, `<q>_pp_cell`, `<q>_pc_cell` | mM | a quantity of the cells of a zone (a species) and its means |
| `<q>_urine` | mmole | the amount of a filtered quantity in the urine of a kidney: the urine volume starts at 0, so its concentration `[<q>_urine]` is 0/0 at the start and the urine is selected by its amount |
| `pO2_<site>` | kPa | the oxygen partial pressure, the free oxygen over its solubility at the temperature |
| `pCO2_<site>` | kPa | the CO2 partial pressure |
| `pH_<site>` | | the pH of the free protons |
| `<e>_act_<site>` | µkat/l | the catalytic activity of an enzyme (ALT, AST, LDH), the unit of the data (1 µkat/l = 60 U/l) |
| `f_alive_<z>`, `f_alive_pp`, `f_alive_pc` | | the living fraction of the cells |
| `VO2` | mmole/min | the oxygen uptake of the organ |
| `bile`, `bile_rate` | µl, µl/min | the bile of a liver |
| `bile_per_mass` | µl/g | the bile per gram of the liver, `bile/(1000 g/kg * M_li)`, the unit of the data |
| `urine`, `urine_rate` | µl, µl/min | the urine of a kidney |

```equations title="observers"
pO2_<z> = o2_vas_<z> / alpha_o2
pCO2_<z> = co2_vas_<z> / alpha_co2
alt_act_res = alt_res * a_alt
VO2 = O2IM_z1 + O2IM_z2
bile_per_mass = bile / (1000 * M_li)
o2_pp = o2_vas_z1
glc_pc = glc_vas_z2
```

The specific activities `a_<e>` (µkat per mmole of enzyme) are placeholders of 1 until the enzyme contents are fitted, see the notes of `[enzymes]` of the parameter database.

## The kidney and the generic organ

The kidney filters the perfusate: the GFR `GFR_ki` is that of an ample flow, `gfr_ki*M_ki` (ml/min/g), at most the maximal filtration fraction `FF_max_ki` of the pump flow, so the kidney stops filtering with the pump; the zones share it by their volume. The tubules reabsorb the fraction `f_reab_ki` of the filtered water and small solutes and none of a tracer (a filtration marker), and the net urine flow `Q_urine` leaves the perfusate. The volumes of the tubing and the zones stay constant: the flow out of a zone is the pump flow less the urine of the zones up to it, and the reservoir, whose volume is a state, loses the urine. The enzymes are not filtered. A generic organ (`GenericOrgan(key, mass)`, the fractions of the organ of the physiology table of the species) has the cellular core and neither bile nor urine.

```equations title="the kidney (two zones)"
GFR_max_ki = gfr_ki * M_ki
GFR_ki = piecewise(GFR_max_ki, GFR_max_ki <= FF_max_ki * Q_pump, FF_max_ki * Q_pump)
Q_urine = (1 - f_reab_ki) * GFR_ki
Q_z1 = Q_pump - (1 / 2) * Q_urine
FLOW_<q>_vas_z1_vas_z2: <q>_vas_z1 -> <q>_vas_z2 = Q_z1 * <q>_vas_z1
FILT_glc_<z>: glc_vas_<z> -> glc_urine = ((1 - f_reab_ki) * GFR_ki / 2) * glc_vas_<z>
FILT_trc_<z>: trc_vas_<z> -> trc_urine = (GFR_ki / 2) * trc_vas_<z>
d(V_res)/dt = -Q_urine
d(V_urine)/dt = Q_urine
urine_rate = 1000000 * Q_urine
```

## Limitations

The open points of the first version, each a decision of the implementation (pkpdsim#146) and a TODO of the model package:

- The oxygen carrier `rbc` (red blood cells) is declared and refused by `validate` until a sourced hemoglobin binding exists; the perfusate is acellular, its oxygen the free oxygen only.
- The liver has a single portal inflow, the rat setup of the study; the split between the hepatic artery and the portal vein is not built.
- The solubility of CO2 `alpha_co2` is held at its value at 37 degC (a `ModelChoice`), the solubility of oxygen follows the temperature; the bicarbonate buffer and the release of the dead cells do not scale with the temperature, so a hypothermic perfusion needs the temperature of the buffer as well.
- The figures draw the organ of a laboratory animal with the shape of the human organ of the anatomy resources (`pkpdsim.perfusion.anatomy`).
- The maximal filtration fraction `FF_max_ki` of the kidney (0.5) is a `ModelChoice` without a source; the kidney has no data and is not validated.
- The enzyme activities are placeholders and most constants of the core are values of the ported model, to be fitted with the study (`ModelChoice`, see the notes of the parameter database).
- The rates of the cells are those of the living cells (`f_T*f_alive_<z>*V_cell_<z>` times a rate per litre of cells), while the concentrations of the cells are per litre of all cells: the dead cells keep their metabolites (glucose, ATP, glycogen, ...), only the enzymes and the potassium move into the dead pool, and the time constants of the cells grow as `1/f_alive_<z>` while the cells die, so a zone with few living cells changes its metabolites slowly.
- A Michaelis-Menten term `s/(km + s)` is singular at `s = -km`, which the integrator reaches only at an undershoot of a species far below 0 (the Hill terms read `abs` and stay defined); the tolerances of the tests keep the species far from it.
