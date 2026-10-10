# API reference

The API reference is generated from the docstrings of the package.

## pkpdsim

| module | description |
| --- | --- |
| [pkpdsim](pkpdsim.md) | the package and its version |
| [templates](templates.md) | `creators` and `terms_of_use` of a model, see [Templates and annotations](../templates.md) |
| [annotations](annotations.md) | `compartments`, the annotations of the recurring compartments, see [Templates and annotations](../templates.md) |
| [data](data.md) | `studies_dir`, `collect_tsv_files`: the study data of a model from the study folders of PK-DB, see [Study data](../data.md) |
| [initial_values](initial_values.md) | `init_concentrations_changes`, `init_amounts_changes`, `get_species_keys`: the initial values of a substance in every tissue, see [Initial values](../initial_values.md) |
| [licenses](licenses.md) | `License`, `license_text`, `write_license`: the texts of CC BY-SA 4.0 and MIT, see [Licenses](../licenses.md) |
| [model](model.md) | `PKPDModel`, `SimulatorSettings`, `FigureSettings`, `ExtraAction`, `resolve`: the description of a model, see [Models and their commands](../models.md) |
| [run](run.md) | `run_experiments`, `run_simulation_experiments`: the run of the simulation experiments of a model |
| [cli](cli.md) | `run_command`, `fit_command`, `report_command`, `identifiability_command`: the commands of a model |
| [fit](fit.md) | `default_fit_settings`, `mapping_metadata`, `observable_id`: the shared fit settings and the helpers of the mapping collections |
| [testing](testing.md) | `check_model`, `check_registry`, `initialize_experiments`, `run_all_experiments`, `start_fits`: the tests of a model as functions |
| [verify](verify.md) | `pkpdsim-verify`: the verification of a change of a model |
