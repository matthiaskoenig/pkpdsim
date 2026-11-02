# pbpk equilibration

The drug-free equilibration of a condition: `equilibrated_simulation`, which `PBPKSimulationExperiment.equilibrated` calls with the file of its loaded model and the groups of a sensitivity analysis use outside an experiment, `drug_free_zeros` (the species of the drugs which start exactly at 0), `accumulating_states` and `growing_states`, see [Impairment](../pbpk/impairment.md#equilibration-of-a-condition). `flat_model` is what the equilibration and the validity of the kidney outflow read from a flat model file, parsed and compiled once per process.

::: pkpdsim.pbpk.equilibration

## flat_model

::: pkpdsim.pbpk.flat_model
