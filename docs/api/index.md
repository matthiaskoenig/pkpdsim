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
| [model](model.md) | `PKPDModel`, `SimulatorSettings`, `FigureSettings`, `CheckSettings`, `ExtraAction`, `resolve`: the description of a model, see [Models and their commands](../models.md) |
| [run](run.md) | `run_experiments`, `run_simulation_experiments`: the run of the simulation experiments of a model |
| [cli](cli.md) | `run_command`, `fit_command`, `report_command`, `identifiability_command`: the commands of a model |
| [fit](fit.md) | `default_fit_settings`, `mapping_metadata`, `observable_id`, `metadata_filter`, `observable_filter`, `select_mappings`: the shared fit settings, the filters and the selection of the mappings of a fit, see [Experiments and their metadata](../experiments.md) |
| [metadata](metadata.md) | `Tissue`, `Route`, `Dosing`, `Fasting`, `ApplicationForm`, `DataType`, `PKPDData`, `PKPDMappingMetaData`: the vocabulary of the metadata of a fit mapping, see [Experiments and their metadata](../experiments.md) |
| [experiment](experiment.md) | `PKPDSimulationExperiment`, `to_sid`, `define_unit`, `rename_units`: the base of the simulation experiments of a model, the SIds of its keys and the units of its data, see [Experiments and their metadata](../experiments.md) |
| [impairment](impairment.md) | the presets of the renal function, the cirrhosis and the cardiac function with their colors |
| [pk](pk.md) | `PK_LABELS`, `PK_UNITS`: the labels and units of the pharmacokinetic parameters |
| [scan](scan.md) | `single_dimension`, `point_timecourse`: the dimension of a scan and the timecourse of a point of its result, see [Scans and sensitivity analyses](../scans.md) |
| [sensitivity](sensitivity.md) | `parameters_from_sbml`, `set_bounds`, `parameters_to_df`, `write_parameters`, `distributions`, `free_parameters`, `sample_size`: the parameters of the sensitivity analyses of a model; `SensitivityRunner` (`run`, `parse_command_line`, `main`) with `LocalAnalysis`, `SamplingAnalysis`, `FastAnalysis`, `SobolAnalysis`, `SensitivityGroup`, `SensitivityOutput`, `CommandLine`: the analyses of its script, see [Scans and sensitivity analyses](../scans.md) |
| [testing](testing.md) | `check_model`, `check_registry`, `initialize_experiments`, `run_all_experiments`, `start_fits`: the tests of a model as functions |
| [verify](verify.md) | `pkpdsim-verify`: the verification of a change of a model |
