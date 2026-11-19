# Unreleased

## Features

- Optional whole-body cardiovascular circulation ([#208](https://github.com/matthiaskoenig/pkpdsim/issues/208)): `Variants(circulation=PressureFlow())` derives steady organ perfusion from pressures and vascular resistance; `DynamicHeart()` adds compliant reservoirs, time-varying ventricular elastance and valve flows. Both retain the portal circulation and hepatic dual supply, support resistance interventions and pharmacodynamic feedback, and use sbmlsim's core engine with SUNDIALS IDA. Fractional perfusion remains the default.
- Interactive cardiovascular examples in the model documentation: shared resistance sliders, heart-rate and contractility controls, live browser IDA integration, flow-width schematics, playback and CSV downloads. The examples include hepatic resistance changes relevant to portal hypertension. Their adult-human parameters are illustrative and require calibration for a specific application.

- Make the cardiovascular explorer a compact dashboard with the physiology figures’ male/female bodies and sex-specific SBML/IDA simulations, adjustable flow-width scaling and prominent synchronized plot markers. Embed the app only in the cardiovascular model documentation.

## Fixes

- Use the core engine for PBPK initialization, renal-flow validation and endogenous balance checks when a model contains algebraic rules. Pressure-driven models can initialize generated stationary concentrations consistently and retain their baseline after a later intervention.
