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
| [timecourse](timecourse.md) | `SimPart`, `Simulation`, `run_simulation`: the timecourse of a roadrunner model which is simulated in parts, see [Timecourses in parts](../timecourse.md) |
| [licenses](licenses.md) | `License`, `license_text`, `write_license`: the texts of CC BY-SA 4.0 and MIT, see [Licenses](../licenses.md) |
| [model](model.md) | `PKPDModel`, `SimulatorSettings`, `FigureSettings`, `CheckSettings`, `ExtraAction`, `resolve`: the description of a model, see [Models and their commands](../models.md) |
| [run](run.md) | `run_experiments`, `run_simulation_experiments`: the run of the simulation experiments of a model |
| [cli](cli.md) | `run_command`, `fit_command`, `report_command`, `identifiability_command`: the commands of a model |
| [profiling](profiling.md) | `run_profiled`, `profiler_available`, `profile_dir`: a command under the sampling profiler of python 3.15 (`--profile`), see [Commands](../commands.md) |
| [fit](fit.md) | `default_fit_settings`, `mapping_metadata`, `observable_id`, `metadata_filter`, `observable_filter`, `select_mappings`: the shared fit settings, the filters and the selection of the mappings of a fit; `fit_subset`, `current_fixed_fits`, `recorded_fixed_fits`, `report_context`, `records_context`: the subset of a fit id and the context of the report of a fit, see [Experiments and their metadata](../experiments.md) |
| [metadata](metadata.md) | `Tissue`, `Route`, `Dosing`, `Fasting`, `ApplicationForm`, `Health`, `Coadministration`, `Genotype`, `Phenotype`, `AdministrationSite`, `DataType`, `PKPDData`, `SamplingSite`, `PKPDMappingMetaData`: the vocabulary of the metadata of a fit mapping, see [Experiments and their metadata](../experiments.md) |
| [experiment](experiment.md) | `PKPDSimulationExperiment`, `to_sid`, `define_unit`, `rename_units`: the base of the simulation experiments of a model, the SIds of its keys and the units of its data, see [Experiments and their metadata](../experiments.md) |
| [impairment](impairment.md) | the presets of the renal function, the cirrhosis and the cardiac function with their colors |
| [pk](pk.md) | `scan_nca`, `nca_frame`: the non-compartmental analysis of a timecourse at every point of a scan, vectorized, and its rows; `PK_LABELS`, `PK_UNITS`: the labels and units of the pharmacokinetic parameters |
| [scan](scan.md) | `single_dimension`, `scan_dimensions`, `point_timecourse`, `point_band`, `start_values`, `end_values`: the dimensions of a scan, the timecourse of a point of its result and the band of its draws, see [Scans and sensitivity analyses](../scans.md) |
| [sensitivity](sensitivity.md) | `parameters_from_sbml`, `set_bounds`, `parameters_to_df`, `write_parameters`, `distributions`, `free_parameters`, `sample_size`: the parameters of the sensitivity analyses of a model; `SensitivityRunner` (`run`, `parse_command_line`, `main`), `run_subfolders` with `LocalAnalysis`, `SamplingAnalysis`, `FastAnalysis`, `SobolAnalysis`, `SensitivityGroup`, `SensitivityOutput`, `CommandLine`: the analyses of its script, see [Scans and sensitivity analyses](../scans.md) |
| [uncertainty](uncertainty.md) | `adopt_fit`, `FitRecord`, `FitRecords`, `ParameterUncertainty`: the adopted fits of a model and the distrib uncertainty of its fitted parameters; `parameter_distributions`, `uncertainty_dimension`, `UncertaintySettings`, `UncertaintyRun`, `apply_uncertainty`, `check_fit_parameters`: the run of the experiments over draws of the fitted parameters, see [Uncertainty of the fitted parameters](../uncertainty.md) |
| [reports](reports.md) | `adopted_record`, `adopted_reports`: the reports of the adopted fits of a model, each written from its record, see [Commands](../commands.md) |
| [release](release.md) | `ReleaseSettings`, `Author`: the release metadata of a model (`PKPDModel.release`); `release`, `main`: the release of a model repository, the command `pkpd_release`, see [Release of a model repository](../release.md) |
| [testing](testing.md) | `check_model`, `check_registry`, `initialize_experiments`, `run_all_experiments`, `start_fits`: the tests of a model as functions |
| [verify](verify.md) | `pkpdsim-verify`: the verification of a change of a model |
| [pbpk](pbpk.md) | `PBPKDefinition`, `Substance`, `Binding`, `Erythrocytes`, `Kp`, `KpMethod`, `Excretion`, `Metabolism`, `Transport`, `CustomReaction`, `Route`, `Variants`, `PDComponent`, `build`, `compose`: the PBPK model library and its modules, see [PBPK library](../pbpk/index.md) |
