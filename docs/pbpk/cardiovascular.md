---
hide:
  - navigation
  - toc
---

# Pressure-driven cardiovascular circulation

The optional circulation variants replace the prescribed organ flows with a 0D cardiovascular network. The organs, chemistry and administration routes are those of the whole-body PBPK model. `Variants(circulation=PressureFlow())` solves the steady resistive network between prescribed arterial and venous pressures; `Variants(circulation=DynamicHeart())` adds a closed systemic and pulmonary loop with ventricular pumping and vascular storage. With no circulation choice the existing fractional-perfusion model is unchanged.

```python
from pkpdsim.pbpk import DynamicHeart, PressureFlow, Variants

averaged = Variants(circulation=PressureFlow())
heartbeat_resolved = Variants(circulation=DynamicHeart())
```

Both constructors accept `parameters={"P_ar": 100.0}` to replace numerical defaults by their short database keys, `resistance={Organ.LIVER: 1.5}` for organ resistance multipliers, and `systemic_resistance`, `hepatic_resistance` and `portal_resistance` for the shared systemic, hepatic outflow and portal-to-sinusoidal multipliers. The liver organ multiplier changes its hepatic arterial resistance; the separate hepatic outflow multiplier changes the resistance after the supplies merge. The lung multiplier applies to the pulmonary resistance of `DynamicHeart`. `DynamicHeart` also accepts `contractility`. Unspecified defaults are read from the parameter database.

Both variants emit SBML algebraic rules and use the **`core` engine of sbmlsim**, which selects SUNDIALS IDA for these models. The pressure-flow equations and the pharmacokinetics are solved together. A second cardiovascular simulator or synchronization step is unnecessary.

## Interactive examples

Change hepatic, systemic, portal and renal resistance with the same sliders in both views. **Pressure network** solves the steady junction balances. **Beating heart** runs the compiled SBML model with sbmlsim-core and SUNDIALS IDA in a browser worker; it also provides heart-rate and contractility controls. Resistance settings carry over when switching views.

<iframe class="pkpd-figure cardiovascular-dashboard" src="examples/cardiovascular.html" title="Interactive cardiovascular model: resistance sliders and animated flow schematic" style="height:780px" allowfullscreen loading="lazy"></iframe>

The compact dashboard uses the male and female whole-body anatomy of the physiology figures. **Physiology** switches both the body and the underlying sex-specific SBML circulation; resistance and heart controls carry over. **Full screen** uses the available display.

Edge widths show instantaneous blood flow during heartbeat playback. Adjust **Edge-width scale** to magnify the shared flow-to-width mapping without changing the model. The orange dots, halos, value labels and cursor lines highlight the displayed instant in every heart plot. Click or drag a time course, or use the time slider, to inspect the final beats. Pause or scrub the simulation to compare pressures, valve opening and regional perfusion at the same time. The [interactive examples walkthrough](cardiovascular_examples.md) explains the controls, the shared flow-width scale and how to reproduce the simulations. [Open the explorer in its own window](examples/cardiovascular.html){ target="_blank" rel="noopener" } for a dashboard that fills its window.

## The network and its reference state

The network follows the PBPK topology. The systemic organs run in parallel from artery to vein; stomach, small intestine, large intestine, spleen and pancreas drain into the portal vein. The liver receives the portal blood and a separate hepatic arterial supply. A sinusoidal node connects these supplies to the hepatic venous outflow. The pulmonary path is in series with the body.

```text
                        systemic organs
                     /------------------\
systemic artery ----+                    +---- systemic vein
                     \-- portal organs -- portal --\
                      \-- hepatic artery ---------- sinusoid --/
```

Each forward-flow resistance obeys `Q = max(P_up - P_down, 0)/R`. This ideal-diode convention also applies to organ beds: the PBPK transport interface assumes the established flow direction and does not represent reversed circulation. The reference resistance is derived from the individual's un-impaired PBPK flow and the chosen reference pressure difference, rather than from a second catalogue of organ perfusion fractions. For a systemic organ `i`, for example,

\[
R_{i,0} = \frac{P_{ar,0}-P_{ve,0}}{Q_{i,0}},
\qquad
Q_i = \frac{\max(P_{ar}-P_{ve},0)\,fimp_{Q_i}}
{R_{i,0}\,fR_i\,fR_{systemic}}.
\]

The reference flows are the PBPK fractions times `CO_base`. Existing organ impairment factors act as conductance multipliers in the pressure-flow equation. At healthy reference conditions the pressure variant reproduces the original organ perfusion at its reference pressures and unit resistance factors. Impairment of the portal organs or liver can shift the portal pressure balance, so the impaired solution need not equal the old independently prescribed flows. Once a resistance changes, the flows are determined by the network equations; the fractional perfusions are reference information, not constraints on the new solution.

Changing a reference pressure through the constructor's `parameters` recalibrates the reference resistances as well as the default boundary. To study a pressure perturbation of the same calibrated network, change `CV_P_ar` or `CV_P_ve` with sbmlsim `preinit_changes`, a timed change or a PD output instead.

The portal and sinusoidal connections are

\[
Q_{po}=\frac{\max(P_{po}-P_{li},0)}{R_{po}},\qquad
Q_{ha}=\frac{\max(P_{ar}-P_{li},0)\,fimp_{Q_{li}}}{R_{ha}},\qquad
Q_h=\frac{\max(P_{li}-P_{ve},0)}{R_{hepatic}}.
\]

The sinusoidal node has zero storage, so `Qha + Qpo - Qh = 0`. In `PressureFlow` the portal node also has zero storage, so `sum(Q_portal_organs) - Qpo = 0`. These two conservation equations determine the portal and sinusoidal pressures. In particular, portal pressure is an output, not a prescribed boundary.

| quantity | unit | convention |
| --- | --- | --- |
| pressure | mmHg | all nodes use the same external-pressure reference |
| blood flow | l/min | the units of the existing PBPK organ flows |
| blood volume | l | whole blood, with PBPK hematocrit used for plasma and cells |
| resistance | mmHg min/l | pressure difference divided by blood flow |
| compliance | l/mmHg | change in blood volume per pressure change |
| ventricular elastance | mmHg/l | pressure divided by stressed ventricular volume |
| heart rate | 1/min | cycles per minute; model time remains minutes |

The numerical defaults and their provenance are in `database/tables/cardiovascular.toml`. They describe an illustrative adult human state: systemic arterial/venous pressures 95/5 mmHg, portal/sinusoidal pressures 10/8 mmHg, pulmonary arterial/venous pressures 15/8 mmHg and heart rate 75/min. These numbers are explicitly labelled `ModelChoice`. They are not a reproduced patient dataset or a validated parameter set of an external model.

## A resistance study in portal hypertension

An increased intrahepatic outflow resistance raises sinusoidal and portal pressures and reduces hepatic flow. This is one mechanical component of cirrhosis; fibrosis, collateral development, altered binding, enzyme abundance and systemic vasodilation require their own mechanisms. The distinction between portal, sinusoidal and venous pressures is also central to the hepatic resistance model of [Levitt and Levitt (2012)](https://doi.org/10.1186/1471-230X-12-26).

The following example builds a complete PBPK model and compares the network before and after increasing hepatic resistance. `preinit_changes` changes the resistance factor before consistent initialization. It preserves the baseline resistance calibration, so the factor represents a perturbation of the same individual. The corresponding constructor option is `PressureFlow(hepatic_resistance=4.0)`.

```python
from pathlib import Path

from sbmlsim.model import AbstractModel
from sbmlsim.simulation import Simulation
from sbmlsim.simulator import Simulator

from pkpdsim.pbpk import PBPKDefinition, PressureFlow, Substance, Variants, build

definition = PBPKDefinition(
    sid="portal_resistance",
    substances=[Substance(sid="tracer", name="tracer", mr=300.0)],
    variants=Variants(circulation=PressureFlow()),
)
built = build(definition, Path("models"))
model = AbstractModel(
    source=built.flat_path,
    selections=["time", "CV_P_po", "CV_P_li", "Qpo", "Qha", "Qh", "CO"],
)
simulator = Simulator(
    engine="core",
    n_workers=1,
    absolute_tolerance=1e-10,
    relative_tolerance=1e-8,
)
for factor in (1.0, 2.0, 4.0):
    result = simulator.simulate(
        model,
        Simulation(
            end=1.0,
            steps=10,
            preinit_changes={"CV_fR_hepatic": factor},
        ),
    )
    print(factor, result["CV_P_po"][0], result["Qpo"][0], result["Qh"][0])
```

<!-- generated: figures-cardiovascular-2 (scripts/visual/docs_figures.py), do not edit -->
The content of the model `portal_resistance`: the anatomical figure (A) and the schematic (B), hover them for details ([figures](figures.md)).

<iframe class="pkpd-figure" src="img/figures/models/cardiovascular-2-anatomy.html" title="The anatomical figure of the model portal_resistance" style="aspect-ratio:1460 / 885" loading="lazy"></iframe>

<iframe class="pkpd-figure" src="img/figures/models/cardiovascular-2-schematic.html" title="The schematic of the model portal_resistance" style="aspect-ratio:1041 / 716" loading="lazy"></iframe>
<!-- end generated -->

For a constant boundary condition the steady network flows remain constant during each simulation. A timed sbmlsim `Change`, or a PD component, can change a factor during the simulation; IDA then restores the algebraic balance at the change. Increasing `CV_fR_portal`, the portal-to-sinusoidal resistance factor, instead isolates a pre-sinusoidal mechanism. Increasing the hepatic outflow resistance isolates a post-sinusoidal mechanism. A pressure difference here is a model-node difference, not automatically a clinical measurement of hepatic venous pressure gradient.

The pre-existing cirrhosis shunt fraction `f_shunt` controls how much drug bypasses functional hepatic tissue. It does not create a new pressure-dependent collateral vessel. Do not interpret a resistance scan alone as a complete cirrhosis progression model or a prediction of TIPS flow.

## The dynamic heart

`DynamicHeart()` retains the organ resistance network and adds left and right ventricles, four ideal one-way valves, pulmonary arterial and venous reservoirs and compliant systemic arterial, systemic venous and portal pools. The dynamic circuit is

```text
vein -> tricuspid -> RV -> pulmonary valve -> pulmonary artery
  ^                                             |
  |                                        lung organ
  |                                             |
  +-- body <- artery <- aortic valve <- LV <- mitral <- pulmonary vein
```

The ventricle uses the time-varying pressure-volume relationship introduced experimentally by [Suga and Sagawa (1974)](https://doi.org/10.1161/01.RES.35.1.117):

\[
P_v = E_v(t)(V_v - V_{0,v}),
\qquad
E_v(t) = \frac{E_{min,v} + a(t)\,CV_{contractility}\,(E_{max,v}-E_{min,v})}
{CV_{scale}}.
\]

The activation `a(t)` is a raised cosine during the systolic fraction of each cycle and zero during diastole. A continuously integrated cardiac phase determines the cycle, so a change in heart rate preserves the phase already reached. The illustrative activation is an implementation choice, not a reproduced experimental waveform.

An ideal valve has `Q = max((P_up - P_down)/R_valve, 0)`. Valve closure is continuous in flow but introduces a derivative discontinuity. Vascular pressures follow a linear compliance relationship, `P = P_reference + (V - V_reference)/C`; each reservoir obeys `dV/dt = sum(Q_in) - sum(Q_out)`. The portal reservoir therefore stores the difference between the splanchnic organ inflows and the portal flow into the liver.

The starting left ventricular blood volume is reallocated from the original arterial pool, and the right ventricular and pulmonary volumes from the original venous pool. They are not added on top of PBPK blood volume. Compliances and ventricular volumes scale by `CV_scale = Vblood/CV_Vblood_ref`, with a reference blood volume of 5.3 l. The pulmonary arterial and venous reservoirs reserve 0.08 and 0.22 of the original venous pool; the remaining systemic arterial and venous volumes must stay positive. The portal compliance is small (0.0003 l/mmHg before scaling), because the explicit PBPK portal pool is a small conduit rather than the entire splanchnic venous bed. The lung's drug transport connects the pulmonary arterial reservoir to the pulmonary venous reservoir; the ventricular and valve paths transport the same substances as the rest of the circulation. The myocardial PBPK organ remains a systemic tissue with its own perfusion.

The changing blood-pool volumes also change their plasma and blood-cell compartments. SBML stores species amounts and computes concentrations using those volumes, so concentration changes include dilution as well as transport. Equilibrium and kinetic erythrocyte models use the PBPK blood/plasma convention. Each internal transfer removes and adds the same drug amount.

To simulate a heartbeat-resolved model, change the variant in the example above to `DynamicHeart()` and request a sufficiently dense output grid. For example, `Simulation(end=1.0, steps=6000)` returns 10 ms output samples for one minute. Output sampling is separate from the adaptive solver's internal steps. Inspect multiple final beats for convergence before interpreting stroke volume, mean flow or pressure-volume loops. Initialization at the reference pressures is consistent with the DAE constraints; it is not a guarantee of a periodic cardiac state.

The principal cardiac outputs are `CV_P_lv`, `CV_P_rv`, `CV_V_lv` and `CV_V_rv`; valve flows are `CV_Q_av`, `CV_Q_tv`, `CV_Q_pv` and `CV_Q_mv` (aortic, tricuspid, pulmonary and mitral). `CV_phase` is the integrated cardiac phase, `CV_HR` the heart rate and `CV_contractility` the active-elastance multiplier. `CV_V_total` sums the central storage volumes; `CV_V_blood_total` includes the fixed organ blood volumes and should remain equal to `Vblood`. The original organ flows `Q<organ>`, `Qha`, `Qpo`, `Qh` and `CO` remain available. `CO` is instantaneous aortic valve flow in the dynamic model, so it is zero during aortic closure and is not instantaneously equal to the lung flow or the sum of the organ flows. Compare cardiac output by averaging complete beats after settling. Cardiovascular unit definitions are available as `CVU` from `pkpdsim.pbpk.cardiovascular` for external PD components which manipulate dimensional cardiovascular quantities.

## Pharmacodynamic feedback

Resistance factors are ordinary parameters of the top model and can be replaced by a `PDComponent` output. This uses the existing composition interface: a concentration is the component's input, and the dimensionless factor is its output. An input may also read an [element of an organ](observers.md#an-element-of-an-organ-for-a-pd-component-or-an-observer). The following component produces an Emax reduction of systemic resistance while keeping the factor positive:

```python
from sbmlutils.factory import AssignmentRule, Model, Parameter, Port, PortType

from pkpdsim.pbpk import MODEL_UNITS, PDComponent, U

vascular_effect = PDComponent(
    model=Model(
        sid="vascular_effect",
        units=U,
        model_units=MODEL_UNITS,
        objects=[
            Parameter(
                "C_drug",
                0.0,
                U.mM,
                constant=False,
                port=Port("C_drug_port", portType=PortType.INPUT_PORT),
            ),
            Parameter("Emax", 0.5, U.dimensionless),
            Parameter("EC50", 0.001, U.mM),
            Parameter(
                "factor",
                None,
                U.dimensionless,
                constant=False,
                port=Port("factor_port", portType=PortType.OUTPUT_PORT),
            ),
            AssignmentRule(
                "factor",
                "1 dimensionless - Emax*C_drug/(EC50 + C_drug)",
                U.dimensionless,
            ),
        ],
    ),
    inputs={"C_drug": "Carm_free_tracer"},
    outputs={"factor": "CV_fR_systemic"},
)
```

Add this component and an intravenous route to the definition above, then simulate a dose:

```python
from dataclasses import replace

from sbmlsim.simulation import Change

from pkpdsim.pbpk import Route

feedback = replace(
    definition,
    sid="vascular_feedback",
    routes=[Route.iv("tracer")],
    pd=[vascular_effect],
)
feedback_built = build(feedback, Path("models"))
feedback_model = AbstractModel(
    source=feedback_built.flat_path,
    selections=["time", "Carm_tracer", "CV_fR_systemic", "CO", "Qki"],
)
response = simulator.simulate(
    feedback_model,
    Simulation(
        end=2.0,
        steps=200,
        changes=[Change(0.0, {"IVDOSE_tracer": 10.0})],  # mg
    ),
)
print(response["CV_fR_systemic"].min(), response["CO"].max())
```

A changing concentration changes resistance, pressure and organ perfusion in the same model that determines subsequent concentrations. The example is a generic vasodilator mechanism. It does not supply a calibrated losartan or antihypertensive clinical model. Individual organ resistance factors allow regional effects instead of a shared systemic effect; the heart variant additionally exposes rate and contractility controls.

## Settling, a resistance step and a dose

A finite drug-free interval lets the heart approach its periodic state. In the following example the hepatic resistance changes after one minute, and the drug is administered after another half minute. This is a single time-domain simulation with timed changes: its state continues across both changes. A stationary `presimulation=SteadyState()` would ask the beating heart to stop changing and is unsuitable.

```python
import numpy as np

from pkpdsim.pbpk import DynamicHeart

heart_definition = replace(
    feedback,
    sid="dynamic_heart_feedback",
    variants=Variants(circulation=DynamicHeart()),
)
heart_built = build(heart_definition, Path("models"))
heart_model = AbstractModel(
    source=heart_built.flat_path,
    selections=[
        "time",
        "CO",
        "Qlu",
        "Qh",
        "CV_P_ar",
        "CV_P_po",
        "CV_V_lv",
        "CV_V_blood_total",
        "Vblood",
        "CV_fR_systemic",
    ],
)
heart_rate = 75.0
cycle = 1.0 / heart_rate  # minutes
heart_solver = Simulator(
    engine="core",
    n_workers=1,
    absolute_tolerance=1e-10,
    relative_tolerance=1e-8,
    maximum_time_step=cycle / 50,
)
heart_result = heart_solver.simulate(
    heart_model,
    Simulation(
        end=2.0,
        steps=12000,  # 10 ms output
        changes=[
            Change(1.0, {"CV_fR_hepatic": 3.0}),
            Change(1.5, {"IVDOSE_tracer": 10.0}),
        ],
    ),
)


def beat_mean(column, first_cycle, last_cycle):
    time = heart_result["time"]
    first, last = first_cycle * cycle, last_cycle * cycle
    mask = (time >= first - 1e-10) & (time <= last + 1e-10)
    return np.trapezoid(heart_result[column][mask], time[mask]) / (last - first)


# Compare complete beats before the resistance step, then after its transient.
print(
    "mean CO before/after resistance step [l/min]:",
    beat_mean("CO", 60, 70),
    beat_mean("CO", 100, 110),
)
print("mean lung flow before step [l/min]:", beat_mean("Qlu", 60, 70))
print("stroke volume before step [ml]:", 1000 * beat_mean("CO", 60, 70) / heart_rate)
print(
    "adjacent beat mean CO [l/min]:", beat_mean("CO", 68, 69), beat_mean("CO", 69, 70)
)
print(
    "maximum total blood volume drift [l]:",
    np.max(np.abs(heart_result["CV_V_blood_total"] - heart_result["Vblood"])),
)
```

The adjacent beat means help assess settling; compare the complete waveforms as well before choosing a settling interval for a new condition. The beat mean of `CO` is the pumped output and its integral over one beat is the stroke volume. The mean lung and aortic flows approach each other in a periodic state, while their instantaneous waveforms differ. The interval used for the resistance comparison ends before the drug is given, so the mechanical resistance step can be distinguished from the subsequent pharmacodynamic response.

## Solver, initialization and scope

The algebraic unknowns are nonconstant SBML parameters with conservation residuals. Ordinary assignment rules define derived flows and pressure-volume relations; rate rules define vascular and chamber storage. The residual system is index 1 in the physiological forward-flow regime with positive resistances and compliances. A disconnected network in which all paths around a node close can lose a unique algebraic pressure solution. sbmlsim's `core` engine performs consistent initialization and integrates it with IDA. An initial pressure guess is not an extra pressure boundary. RoadRunner cannot simulate these algebraic-rule models.

For the commands and experiments of a [model package](../models.md), choose the same engine in the model's `simulator` settings:

```python
from sbmlsim.project import SimulatorSettings

cardiovascular_solver = SimulatorSettings(
    engine="core",
    absolute_tolerance=1e-10,
    relative_tolerance=1e-8,
)
```

Pass these settings as `PKPDModel(..., simulator=cardiovascular_solver)` so runs and fits use the DAE-capable engine.

For the dynamic heart, `maximum_time_step=1/(HR*50)` minutes is a useful starting bound for the fast valve transitions. Tighten it and the solver tolerances to check convergence of pressure extrema, stroke volume and exposure for the actual application. The equidistant output grid does not constrain the integrator's internal steps, and sparse output can miss pressure peaks even when integration is accurate.

Explicit parameter values, timed changes and PD coupling are supported. sbmlsim's reference-based generic parameter sampling currently reloads models through RoadRunner (`sampling.references()` and `parameters_of()`), so that path is unavailable for these DAE models even with the `core` engine selected. The verification here does not establish support for every sensitivity, sampling or optimization API.

All resistance and compliance values must be positive. The pressure variant represents steady forward circulation between its chosen boundaries. The heart variant is a lumped model with fixed heart-rate/contractility defaults and no baroreflex, autonomic regulation, renal autoregulation, hepatic arterial buffer response, vessel collapse, valve leaflet dynamics, inertance, wave propagation or remodeling. Clinical calibration is needed before drawing patient-specific conclusions.

`DynamicHeart` currently supports adult human models with static anatomy. Age-dependent physiology, laboratory-animal physiology, pregnancy, dynamically changing body composition and generated stationary custom chemistry are outside that contract. Pregnancy is also unsupported by `PressureFlow`, because the cardiovascular network has no fetoplacental branch. A periodic heart is not a biochemical steady state: use a time-domain settling simulation rather than asking a constant-state steady-state solver to converge a beating heart.

`PressureFlow` supports the library's generated stationary endogenous concentrations. These use the consistently initialized pressure network, including an altered hepatic resistance, and retain their initial baseline after a later intervention. A generated baseline whose perfusion depends on concentration-driven PD is refused because its initialization would be circular; that condition requires an explicit initialization strategy.

## Existing implementations and provenance

The equations are independently encoded in the SBML model builder using established resistance, compliance, elastance and diode-valve components. No external implementation is vendored and none of the following libraries is a runtime dependency. Their published implementations provide useful resources for more detailed components and future comparisons.

| resource | interface and relevant content | license | relationship to this implementation |
| --- | --- | --- | --- |
| [svZeroDSolver](https://github.com/SimVascular/svZeroDSolver), [documentation](https://simvascular.github.io/svZeroDSolver/) | C++ solver and Python `svzerod`; JSON networks of 0D hemodynamic components | BSD-3-Clause | reference for modular pressure-flow networks; the runtime solver here is sbmlsim/IDA |
| [ModularCirc](https://github.com/alan-turing-institute/ModularCirc) | Python chambers, valves and vessels; published Naghavi and Korakianitis models | MIT | reference for modular circulation composition; its model parameter sets are not reproduced here |
| [circulation](https://github.com/ComputationalPhysiology/circulation) | Python implementations including `Regazzoni2020`; pressure-volume and whole-circulation simulations | MIT | candidate independent reference for a fuller cardiac model; this is not a Regazzoni-model port |

The physiological references establish the component relationships, not the clinical validity of the new PBPK network. The parameter database keeps that distinction explicit. Validation checks algebraic conservation, reference-flow recovery, pressure responses to resistance, blood-volume and drug-mass conservation, and heartbeat dynamics. `tests/pbpk/test_cardiovascular_reference.py` also compares the generated heart against an independently assembled hydraulic circuit integrated with SciPy DOP853, including a changed hepatic resistance: the explicit reference eliminates hepatic pressure with a scalar conservation solve and checks all seven storage volumes and eight pressures. This verifies the SBML translation and integration independently of IDA; it is not a reproduction of a complete external model parameter set. A comparison of a fitted whole-circulation parameter set against an external implementation remains a separate validation task.
