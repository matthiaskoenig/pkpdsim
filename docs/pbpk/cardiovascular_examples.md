# Interactive cardiovascular examples

Explore how vascular resistance changes pressures and perfusion in the whole-body circulation. This explorer is part of the [cardiovascular model documentation](cardiovascular.md#interactive-examples). Both model views have the same hepatic, systemic, portal and renal resistance sliders, and retain those settings when you switch views. The heart view also provides heart-rate and contractility sliders.

[Open the dashboard in the cardiovascular model documentation](cardiovascular.md#interactive-examples). The app is embedded only there. Use **Full screen** or open its own window to fill your display. It runs locally in the browser; CSV downloads contain the displayed values and units.

Choose **Male** or **Female** under **Physiology**. This selects the corresponding whole-body anatomy, PBPK organ-flow calibration and compiled heart model. The two versions use the same illustrative cardiovascular pressure and elastance assumptions; selecting female physiology does not imply a separately clinically calibrated female heart. Resistance, rate and contractility settings carry over between versions, while the results and reference comparison remain specific to the selected sex.


## Read the flow schematic

The dashboard overlays flow on the male and female bodies from the [physiology figures](figures.md). It uses the same theme-aware inline SVG mechanism as the [sbmlsim-core illustrations](https://matthiaskoenig.github.io/sbmlsim-core/javascript/): explicit compartments, directed connections and labels drawn from the simulated values. The heart view adds ventricular chambers, the pulmonary circuit and the four valve paths to the whole-body organ branches.

Connections describe the model topology rather than anatomical vessel trajectories. An edge becomes wider as its instantaneous flow increases. All edges use the same scale, with width growing with the square root of flow; the legend gives examples in l/min. **Edge-width scale** adjusts this mapping from 0.25× to 3× for all vessels and the legend together. It changes only the display. A closed valve retains a thin path to show the connection while its flow is zero. Node labels show pressures, and edge labels show the instantaneous flow at the displayed time.

The orange point and halo, value label, shaded band and vertical cursor show the displayed instant on the heart plots. Click or drag a time course to move it, or use the keyboard-accessible time slider. The **Final 3 beats** view is selected initially.

Use the playback controls to play or pause, change replay speed, or scrub to a specific time. Replay advances through an already computed IDA trajectory; changing a model slider computes a new trajectory. This lets you distinguish a change in the physical model from a change in playback speed.

## Increased hepatic resistance

Choose **Pressure network**, then set **Hepatic outflow resistance** to **4×**, or use its preset. The resistance after the hepatic arterial and portal supplies merge rises fourfold. Portal and sinusoidal pressures increase, and portal and total hepatic perfusion decrease. The arterial and systemic venous pressure boundaries remain fixed.

```python
from pkpdsim.pbpk import PressureFlow, Variants

hepatic_resistance = Variants(circulation=PressureFlow(hepatic_resistance=4.0))
```

The network diagram locates the pressure rise upstream of the resistance. The pressure-response plot shows the full resistance scan while keeping the other controls at their current values. Open **All organ flows** to inspect the regional perfusion changes. Switch to **Beating heart** to apply the same resistance setting in the closed loop and replay the resulting flow through the schematic.

This isolates a mechanical component of portal hypertension. It does not introduce collateral vessels, a pressure-dependent shunt, hepatic arterial buffering or a complete cirrhosis phenotype. The displayed portal-to-systemic venous gradient is a difference between model nodes, not a simulated clinical HVPG measurement.

## Systemic and regional resistance

Reset the pressure network, then set **Systemic resistance** to **0.7×**, or use its preset. At prescribed systemic pressure boundaries, reducing arterial organ resistance increases blood flow. The total flow can change because the averaged model supplies the pressure boundaries externally.

```python
vasodilation = Variants(circulation=PressureFlow(systemic_resistance=0.7))
```

Now change **Renal resistance** alone. Renal perfusion changes independently of the other systemic branches. Change **Portal resistance** to distinguish resistance before the sinusoidal node from hepatic outflow resistance. In **Pressure network**, **Pressure boundaries** lets you investigate venous congestion or a different arterial pressure without recalibrating the reference resistance network. Those prescribed boundaries apply only to the steady model; the closed-loop heart predicts its own arterial and venous pressures.

The sliders calculate the exact two-junction steady network in JavaScript using resistance values exported from the SBML model. The browser result is checked against independent sbmlsim/IDA simulations. It does not simulate vascular storage or beating-heart dynamics.

## A beating heart and a resistance intervention

Choose **Beating heart** and set **Hepatic outflow resistance** to **4×**. The same continuous resistance sliders are available here as in the pressure network, so hepatic, systemic, portal and renal perturbations can be combined. Each calculation starts from the same initial conditions. The circulation settles for 24 seconds, then the selected intervention is applied. Live runs end on a complete post-intervention beat near 48 seconds, and the final values use the last three complete beats.

```python
from pkpdsim.pbpk import DynamicHeart

heartbeat = Variants(circulation=DynamicHeart())
# In a Simulation with time in minutes:
# Change(0.4, {"CV_fR_hepatic": 4.0})
```

The arterial and venous pressures are now outputs of a closed loop. The heart, compliant vessels and organ resistance network respond together. Compare the full response with **Final 3 beats**, and inspect the left ventricular pressure-volume loop. Pause the schematic during systole and diastole to see the changes in valve flow. Aortic flow is pulsatile and falls to zero while the aortic valve is closed. **Mean cardiac output** averages that flow over complete beats; **stroke volume** is its per-beat integral.

Eight presets configure the controls for common examples. Each preset can be changed further with the sliders:

| example | parameter changed at 24 s | purpose |
| --- | --- | --- |
| Reference circulation | none | settling, normal waveform and baseline comparison |
| Hepatic resistance 2x, 4x, 8x | `CV_fR_hepatic` | increasing hepatic outflow resistance |
| Lower systemic resistance | `CV_fR_systemic = 0.7` | a vasodilation mechanism in the closed loop |
| Higher systemic resistance | `CV_fR_systemic = 1.4` | increased afterload |
| Reduced contractility | `CV_contractility = 0.65` | reduced active ventricular elastance |
| Faster heart rate | `CV_HR = 100` | changing rate while preserving the integrated cardiac phase |

The heart view integrates the compiled SBML model in the browser with **sbmlsim-core WebAssembly and SUNDIALS IDA**, using a worker so that simulation does not block the controls. Changing a slider reruns the model with that intervention. Saved reference cases provide the initial comparison and numerical checks; arbitrary combinations are solved directly rather than interpolated between saved examples. Readouts use complete beats of the solution, and the total blood-volume error accompanies each calculation.

## Reproduce the examples

The source of the numerical examples is [the reference-data generator](https://github.com/matthiaskoenig/pkpdsim/blob/develop/scripts/cardiovascular_examples.py). From a checkout with the development environment installed, run:

```bash
uv run python scripts/cardiovascular_examples.py
```

This builds the pressure-driven and dynamic-heart SBML models, solves the eight reference heartbeat cases with IDA, and writes the male reference data to `docs/pbpk/examples/cardiovascular-data.json` and the female data under `docs/pbpk/examples/female/`. The anatomy export `cardiovascular-anatomy.json` reuses `reference_body` from the physiology figures, including sex-specific organ anchors and the EMBL-EBI anatomogram attribution. It also exports the compiled model and browser runtime under `docs/pbpk/examples/solver/`, with a separate compiled female model under `female/solver/`. The common runtime is shared, including the [runtime license](examples/solver/LICENSE) and [third-party notices](examples/solver/THIRD_PARTY.md). The app uses these assets to run new slider combinations locally. The reference data records solver tolerances, package versions, units and a hash of the SBML models. To regenerate somewhere else, pass `--output /path/to/cardiovascular-data.json`.

The numerical defaults describe illustrative male and female adult humans. They establish example behavior and numerical conservation, rather than clinical calibration. For a coupled drug example, a complete build and a finite settling-and-dosing simulation, see [Cardiovascular circulation](cardiovascular.md#pharmacodynamic-feedback).
