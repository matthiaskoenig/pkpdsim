# Commands

pkpdsim installs six commands which work with every model package. The model is the first argument of a command, the dotted name of its package; everything after it are the options of the command for that model:

- `run_simulation` creates the models with the factory, runs the simulation experiments and the extra actions of the model.
- `run_fit` fits the parameters of a fit problem, analyses their identifiability, or adopts a finished fit.
- `run_report` writes the report of a fit for stored parameter sets.
- `run_identifiability` analyses the identifiability of a stored parameter set.
- `run_credibility` assesses the contexts of use of a model against its adopted fits.
- `pkpd_release` makes the complete release of the repository of a model.

```bash
run_simulation pkdb_models.models.hctz --action simulate --experiments all
run_fit pkdb_models.models.hctz --subset PK --runs 10 --cores 10
run_report pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_credibility pkdb_models.models.hctz
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
```

A model package registers no commands of its own: the commands are defined once, in pkpdsim, and everything which differs between two models comes from the model object `MODEL` of the package (its actions, fit problems, settings and folders) or from the options of the command line.

## The run layer of sbmlsim

`run_simulation`, `run_fit`, `run_report`, `run_identifiability` and `run_credibility` are wrappers of the run layer of sbmlsim, `python -m sbmlsim.run` (see the documentation of sbmlsim, `sbmlsim.run`): they keep their options and defaults, translate them into the spec of a kind and run it with `sbmlsim.run.execute` on the model object of the package, the project `<package>:MODEL` (`pkpdsim.cli.project_reference`). A kind runs into an empty output directory and writes the records of the run next to its results: `run.json` (the kind, the versions, the host, the status), `spec.json` (the planned spec), `summary.json` and `progress.jsonl`. A command therefore runs the kind into the folder of its result below the folder it always wrote into. The run writes into a hidden staging directory next to that folder, `.<name>.<id>`, and replaces an earlier result only after it succeeded (the earlier result is renamed aside, the new one renamed into place, the earlier one removed); a run or a plan which fails leaves the earlier result as it was, and a run which was killed leaves its staging directory, which the next run into the folder removes (it says so). An error of the run itself keeps its traceback; a result which is not moved into its place exits with 1 and names where the earlier result is when it could not be restored. An output which is an input of the command (the records of the fit it reads, the adopted fits `MODEL.fits_dir`, the factory directory or a build), contains one or lies inside one, or which contains the package or a folder of its data, is a usage error before anything runs, e.g. `run_report --name <fit id> -o results/fit results/fit/<fit id>/parameters.json`. A spec which does not plan is a usage error (exit 2), a configuration error of the run itself exits with 1:

| command | kind | output of the run | the spec |
|---------|------|-------------------|----------|
| `run_simulation -a simulate -e misc` | `simulate` | `results/simulation/misc` | `{project, experiments: [misc], figure_formats, uncertainty}` |
| `run_simulation -a factory --rebuild` | `factory` | `results/factory` | `{project}` |
| `run_fit --subset PK` | `fit` | `results/fit/<fit id>` | `{project, subset, opid, repeats, algorithm, seed, studies, timeout, start, record_adopted: true, identifiability, profile, fisher, report: true, mapping_figures}` |
| `run_fit --adopt <fit id>` | `adopt` | a temporary folder, the record copied into `MODEL.fits_dir` | `{project, force}`, the input `fit` the folder of the fit |
| `run_report <parameters.json>` | `report` | `results/report/<name>` | `{project, subset, mapping_figures, group_by}`, the input `fit` the folder of the file |
| `run_report --adopted` | `report` | `results/report/<fit id>` | `{project, subset, adopted: true, mapping_figures}` |
| `run_identifiability <parameters.json>` | `identifiability`, `report` | `results/identifiability/<name>/analysis`, `.../report` | the profile settings, `pids`, `fisher`; the report of the analysis |
| `run_credibility` | `credibility` | `results/credibility` | `{project, verification, lock, uncertainty, draws, seed, variability, subjects, studies, source_commit}`, the input `record` the file of `--from-json` |

The same run from the command line of sbmlsim, e.g. for a fit:

```bash
python -m sbmlsim.run fit spec.yaml --out results/fit/run --workers 10
```

with `spec.yaml` the spec of the table, `project: pkdb_models.models.hctz:MODEL`, `subset: PK`, `repeats: 10`. The kind `fit` writes the records of the fit into the folder `<fit id>` of its output; `run_fit` moves them up into `results/fit/<fit id>`, so the records, the report and `parameters.json` are where they always were, next to the records of the run. The report of the kind `report` (its folder `report`) and of an adopted fit (`<fit id>`) are moved up the same way; `run_identifiability` writes the report of the analysis into `results/identifiability/<name>` itself, the run of the analysis into `analysis/` and the records of the run of the report into `report/` of it. The kind `credibility` writes its outputs into `credibility/` of its output and the lock of `lock: true` into the output itself; `run_credibility` moves the outputs up into `results/credibility` and the lock into the folder of the model, where the next assessment compares with it.

Three parts use the library of sbmlsim directly instead: `--action factory` builds the models of the package in place (`create_models`), the extra actions write into the results of the model, which are no empty output (they are the `actions` of the project as well, so `python -m sbmlsim.run action` runs them, see [Extra actions](#extra-actions)), and `run_fit --strategy single`, `run_report` of several parameter files and a parameter file without the records of its fit next to it run with `sbmlsim.fit.cli` (`run_fit`, `report_cli`, `identifiability_cli`), which no kind expresses.

## The model argument

The model is given as the dotted name of its package, e.g. `pkdb_models.models.hctz`. The command imports the package and takes its model object `MODEL`, the `PKPDModel` of its `__init__.py` (see [Models and their commands](models.md)). The package has to be importable in the environment the command runs in, i.e. installed, e.g. the checkout of a repository of models with `uv sync` or `pip install -e .`; a console script does not import from the working directory.

The model comes first, since the options of a command depend on the model: the actions of `run_simulation` include the extra actions of the model, and the choices and defaults of the fit commands are its fit problems and its folders. An option before the model is a usage error.

`--help` without a model prints the usage of the command, with a model the options of that model:

```bash
run_fit --help                                  # the usage of the command
run_fit pkdb_models.models.hctz --help          # the options for hctz
```

A model which cannot be resolved stops the command with a usage error and the exit code 2:

| cause | error |
|-------|-------|
| no argument | `the model package is required` |
| an option first | `the model package is the first argument, not '--subset'` |
| not a dotted name, e.g. a path | `'pkdb_models/models/hctz' is no dotted package.` |
| the package or one of its parents is not installed | `The model package 'pkdb_models.models.hctzz' is not found.` |
| no `MODEL` in the package | `The package '...' has no model object MODEL.` |
| `MODEL` is not a `PKPDModel` | `MODEL of the package '...' is no PKPDModel but a ...` |
| `MODEL` of another package, e.g. a module which imports it | `MODEL of the package '...' belongs to the package '...'.` |
| a fit command for a model without fit definitions | `the model 'name' (...) has no fit definitions` |

Any other error of the import of the package, e.g. a module which the model imports and which is not installed, is raised with its traceback, so that a broken model shows where it breaks.

`pkpdsim.cli.resolve_model(package)` is the resolution as a function: it returns the model object or raises `ModelNotFoundError` with the messages of the table.

## run_simulation

`run_simulation` performs one `--action` (`-a`):

- `factory` creates the models with the factory of the model: the SBML, the ODEs and the COMBINE archive, in the package. With `--rebuild` the kind `factory` builds them into `factory` of the results instead and reports the outputs of the factory which differ from the ones of the package; the command then exits with 1.
- `list_experiments` lists the groups and experiments of the registry of the model.
- `simulate` runs the experiments of `--experiments` (required).
- `all` runs the factory, all experiments and every extra action with `in_all`.
- an extra action of the model, see below.

```bash
run_simulation pkdb_models.models.hctz --action factory
run_simulation pkdb_models.models.hctz --action list_experiments
run_simulation pkdb_models.models.hctz --action simulate --experiments "misc,Anderson1961" --open
run_simulation pkdb_models.models.hctz --action all --results-dir /path/to/results --formats svg,png
run_simulation pkdb_models.models.hctz --action simulate --experiments all --workers 4
run_simulation pkdb_models.models.hctz --action simulate --experiments all --uncertainty --draws 100
```

The options:

- `--experiments` (`-e`): comma separated groups and experiments of the registry, e.g. `all`, `misc,Anderson1961`. A single group is written into the folder of the group, any other selection into `custom_selection`, both in the `simulation` folder of the results; an earlier result in the folder is replaced. A name which is neither a group nor an experiment is reported and left out, a selection without any experiment exits with 1.
- `--results-dir` (`-r`): the folder of the results, `MODEL.results_dir` by default.
- `--formats` (`-f`): comma separated formats of the figures, the static images of matplotlib (`svg`, `png`, ...) and the interactive pages of plotly (`html`); the formats of the `FigureSettings` of the model by default. The report shows the SVG and links the interactive page; without `svg` the interactive page is embedded in the report.
- `--open` (`-o`): opens the report of the simulations in the web browser.
- `--workers` (`-w`): the processes of the run: a number (the kind `simulate` with the experiments in a pool of that many processes, at most one per experiment, every experiment serially in its worker, a single experiment with its scans in the pool), `model` (the default of a model and of the option: the kind with a worker per four experiments, at most one per CPU, serially below, every experiment with its scans serially in its worker) or `auto` (sbmlsim's choice: a worker per four experiments, or the scans in a pool of every CPU when one has 256 points or more; the kind resolves its workers, so `auto` runs the experiments with `Project.run_experiments` of sbmlsim itself, without the records of a run); an extra action gets the number, 1 for `model` (`pkpdsim.model.RUN_WORKERS`) and `None` for `auto`, see [Models and their commands](models.md).
- `--rebuild`: the factory of `factory` and `all` builds into the results and compares, see above.
- `--uncertainty` (`-u`) with `--draws` (100), `--seed` (1234) and `--quantiles` (0.05 0.95): runs the experiments over draws of the fitted parameters of the model and draws the curves of the simulations as bands, into the folder `simulation_uncertainty` of the results, see [Uncertainty of the fitted parameters](uncertainty.md).

A run which fails is not the end of the command. An experiment which cannot be initialized or whose figure fails is reported, the other experiments run, the report and the figures are written, and the command then runs its remaining actions, prints every failure and exits with 1.

A model without experiments (`experiments=None`) offers only `--action factory`, `--action all` (which runs only the factory and the extra actions with `in_all`) and its extra actions, and none of the options of the experiments; `simulate` and `list_experiments` are a usage error for it.

### Extra actions

A model adds actions of its own with `ExtraAction` and their options with `extra_arguments` of the model object (see [Models and their commands](models.md#commands)), e.g. analyses which are not simulation experiments:

```bash
run_simulation pkdb_models.models.glimepiride --action analyze_cyp2c9
```

`run_simulation` calls an extra action as before, with the folder of the results. The extra actions are also the `actions` of the project of the model (`pkpdsim.model.ExtraActionFunction`), which the kind `action` of sbmlsim runs into its output, `python -m sbmlsim.run action spec.yaml --out DIR` with the spec `{project: <package>:MODEL, name: analyze_cyp2c9, arguments: {analyses: pk}}`: the arguments are the options of the action over their defaults, `options.n_workers` the processes of the run and `options.uncertainty` `None`.

An extra action gets the options of the command line, the settings of the simulator of the run as `options.simulator` (the `SimulatorSettings` of the model), the processes of `--workers` as `options.n_workers` and the settings of `--uncertainty` as `options.uncertainty` (`None` without). An action which simulates passes them on, `model.run_experiments(..., simulator=options.simulator, n_workers=options.n_workers)` or `options.simulator.simulator(options.n_workers)` (`sbmlsim.project.SimulatorSettings.simulator`) for a runner of its own, otherwise `--workers` does not reach it. An option of the model named `--simulator` or `--n-workers`, or `--uncertainty` for a model without experiments, is refused with a `ValueError` before the action runs, because it would be replaced.

## run_fit

`run_fit` fits the parameters of a fit problem of the model (`--subset`, the first fit definition by default) and, with `--identifiability`, analyses the identifiability of the fitted parameters: the profile likelihood of every fitted parameter and the Fisher information around the best parameter set, which end up in the identifiability section of the report of the fit. The analysis is the default of a model with `fit_identifiability`, which `--no-identifiability` skips.

```bash
run_fit pkdb_models.models.hctz --subset PK --runs 10 --cores 10 --max-points 30
run_fit pkdb_models.models.hctz --subset PK --method DE --strategy SINGLE
run_fit pkdb_models.models.hctz --adopt PK_<fit id>
run_fit pkdb_models.models.inulin --subset PK --runs 1 --start profile_point.json
```

The options of the fit are those of sbmlsim (`sbmlsim.fit.cli.fit_cli`): `--subset` (`-x`), `--name` (`-n`, the name of the reports of `--strategy SINGLE`; ignored with a warning otherwise, the report of a fit is in the folder of its id), `--runs` (`-r`), `--cores` (`-c`), `--seed` (`-s`), `--method` (`-m`, `LSQ` least squares with latin hypercube starts or `DE` differential evolution), `--strategy` (`-t`, `ALL` the experiments together or `SINGLE` every experiment on its own), `--experiments` (`-e`), `--timeout`, `--start`, `--no-mapping-figures` and `--output_dir` (`-o`, the `fit` folder of the results by default). The options of the identifiability are `--parameter` (`-p`, repeatable), `--alpha` (`-a`), `--df`, `--max-points`, `--no-reoptimize` and `--no-fisher`.

Least squares starts every run from a sample of the bounds of the parameters (log-uniform latin hypercube sampling, `sbmlsim.fit.sampling`); the start value of a parameter (`FitParameter.start_value`) is neither the first start nor the center of the samples, it is used only for a parameter with an infinite bound, which is not sampled. `--start` starts least squares from a given point instead, e.g. a point of a profile likelihood below the cost of the fit: a json file of the values by parameter id in the units of the fit parameters (`{"kslow_inu": 0.002}`), or the id of a fit, an adopted fit of the model (`MODEL.fits_dir`) or else a finished fit in the output folder, whose best parameters are converted into the units of the fit parameters. The point replaces the start values of the parameters it gives (`sbmlsim.fit.project_fits.start_point`, `with_start`, on copies of the parameters), and every run starts from the start values (the sampling `START` of sbmlsim), so a fit with `--start` has one run (`--runs 1`); a parameter which the start does not give starts from its start value, a parameter of the file which is not fitted, a value which is no number or outside of the bounds, `--runs` other than 1 and differential evolution, which draws its population from the bounds, are usage errors.

The problem of the fit, of its identifiability analysis and of its report loads the model at the optimization level of the model (`SimulatorSettings.jit_optimization`), as `run_report` and `run_identifiability` do: sbmlsim builds a fit problem without a level, so the commands hand sbmlsim the fit definitions at the level (`PKPDModel.load_fit_definitions`, `FitDefinition.jit_optimization` of sbmlsim), whose problems carry it into every worker of a parallel fit; the environment of the process and of a forkserver which starts the workers does not matter, and nothing is set in it. The selection of the mappings (`select_mappings`, which instantiates the experiments once in the process) is no problem of the fit and keeps the level of sbmlsim.

The fit runs the kind `fit` (`--strategy ALL`, the default) into `results/fit/<fit id>`: the repeats in `runs/`, `optimization_result.json`, `parameters.json`, `fits.tsv`, with `--identifiability` `identifiability.json` and `fisher.json`, the report and the records of the run. A fit which fails keeps the records it wrote in its folder, so a later `--start <fit id>` finds them. The plan of the kind refuses a fit whose bounds let a remainder of a model become negative, an invalid start and an unknown study before anything runs, a usage error with the exit code 2. The report of the fit carries the title of the model (`PKPDModel.title`), the description of the fit (`FitDefinition.description`) and of its parameters (`FitParameter.description`), and a note on the adopted fits of the other subsets of the model, whose parameters were fixed at their fitted values. The provenance of the fit records `start` and `fixed_fits` (the ids of the fixed adopted fits), the start values of the parameters and the fixed values, i.e. the values of the parameters of the fixed adopted fits which the fit does not optimize itself; the options of the command are the spec of the run, `spec.json`. A fit of `--strategy SINGLE` records the arguments of the command (`subset`, `runs`, `cores`, `method`, `command`, `start` and `fixed_fits`) as before.

`--cross-validation <fit id>` cross-validates an adopted fit instead of fitting: a refit per training study without it, stored in the record of the fit, see [Cross-validation](credibility.md#cross-validation). `--adopt <fit id>` checks that `<fit id>` is the id of a fit with its `parameters.json` in the output folder, runs the kind `adopt` on it into a temporary folder, which checks the fit and writes its record, and copies the record into the adopted fits of the model (`MODEL.fits_dir`) instead of fitting, `--force` replaces an adopted record of the same id; the id then goes into `fits` of the model object, see [Uncertainty of the fitted parameters](uncertainty.md).

## run_report and run_identifiability

Both work on the `parameters.json` of an earlier fit, the kinds `report` and `identifiability` with the records of the fit next to it as the input `fit` (several parameter files of `run_report`, or a file of another name, with `sbmlsim.fit.cli.report_cli` and `identifiability_cli`):

```bash
run_report pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.hctz --subset PK --max-points 30 results/fit/PK_<fit id>/parameters.json
```

`run_report` writes the report of the fit problem for one or more files of parameter sets, `run_identifiability` the profile likelihood and the Fisher information of the first parameter set of a file. `--subset` is the fit problem, `--output_dir` the folder of the reports and `--name` the folder of the report in it, `report` and `identifiability` by default.

The report of `run_report` has the context of the model as the report of `run_fit` has it: the title of the model, the description of the fit and of its parameters and the note on the fixed adopted fits, which are the fits the provenance of the record next to the first file names, not the adopted fits of the model now (a record without provenance names none).

`run_report <package> --adopted` writes the report of every adopted fit of the model (`MODEL.fits`) instead, each from its record in `MODEL.fits_dir`, without a parameters file and without running a fit:

```bash
run_report pkdb_models.models.hctz --adopted
```

Every fit has a report of its own in `MODEL.report_dir/<fit id>` (`--output_dir` replaces `MODEL.report_dir`). It is the report which `run_report <package> --subset <subset> --name <fit id> <record>/parameters.json` writes for the record: the complete report of sbmlsim with the title of the model, the description of the fit and its parameters, the adopted fits whose parameters were fixed as the provenance of the record names them, and the provenance of the record (a record without it states that it was not recorded and names no fixed fits). The subset of a fit is the part of its id before the date, `PK` for `PK_20261016_185018__3efe3`, and must have a fit definition. There is no report which merges the fits of a model.

## run_credibility

Assesses the contexts of use a model declares (`PKPDModel.credibility`, see [Context of use and credibility](credibility.md)) against its adopted fits. The adopted fits are evaluated on the current model and data, the kinds of their mappings are those of their records; the command needs a declaration and at least one adopted fit:

```bash
run_credibility pkdb_models.models.hctz
run_credibility pkdb_models.models.hctz --output_dir results/credibility
run_credibility pkdb_models.models.hctz --no-verification
run_credibility pkdb_models.models.hctz --uncertainty --draws 128
run_credibility pkdb_models.models.hctz --variability --subjects 64
run_credibility pkdb_models.models.hctz --from-json results/credibility/credibility.json
```

It writes `credibility.json` (the assessment), `contexts.tsv`, `evidence.tsv`, `domain.tsv`, `quantities.tsv`, `ratios.tsv`, `temporal.tsv`, `studies.tsv`, `inputs.tsv`, `input_ranges.tsv`, `cross_validation.tsv`, `uncertainty.tsv`, `uncertainty_quantities.tsv`, `uncertainty_sensitivity.tsv`, `variability.tsv`, `variability_quantities.tsv`, `vpc.tsv`, `verification.tsv`, `battery.tsv`, `fates.tsv` and the report `index.html` into `MODEL.credibility_dir` (`results/credibility`), next to the records of the run of the kind `credibility` (see [The run layer of sbmlsim](#the-run-layer-of-sbmlsim)), and prints a table of the contexts of use with their verdicts, the role of the evidence which decided them (decided on: calibration for a description, C1, else validation) and the evaluable points with their minimum (`evaluable/minimum`), and how many of them are verified. The decisive predictions are verified (the fits evaluated again with tightened tolerances and with a second integration method, and the builds checked, see [Verification](credibility.md#verification)) and the builds run the battery of their checks (mass balance, non-negativity, zero dose, steady state, superposition, AMICI as a second simulator when installed, plausibility and rebuild, the battery of sbmlsim with the PK battery of pkpdsim, see [Battery of the builds](credibility.md#battery-of-the-builds)); `--no-verification` skips both, the factors of the verification are then not assessed. `--uncertainty` propagates the parameter uncertainty of the adopted fits to the decisive predictions: every fit is evaluated at `--draws` draws of its parameters (128 by default, with the seed `--seed`, 1234), which costs as many evaluations of every fit, and the factor model inputs is assessed (see [Parameter uncertainty](credibility.md#parameter-uncertainty)); `--draws` or `--seed` without `--uncertainty` and fewer than two draws are usage errors. `--variability` checks the decisive points with the variability of the subjects the declaration has (`Credibility.variability`): every fit is evaluated for `--subjects` subjects (64 by default) at every parameter vector, the adopted one or the draws of `--uncertainty`, with a variability between the studies in `--studies` studies (16 by default) of `--subjects` subjects each, and the coverage of the bands of the predictive check is recorded (see [Variability and the predictive check](credibility.md#variability-and-the-predictive-check)); `--subjects` or `--studies` without `--variability`, fewer than two subjects, no study and `--variability` with `--from-json` are usage errors, a declaration without a variability is an error. `--output_dir` (`-o`) replaces the directory; an output which is or contains an input of the model (the adopted fits, the factory, a build) or contains the package or its data is a usage error. `--from-json` writes the outputs again from a written `credibility.json` without evaluating the fits. `--lock` writes the lock of the evidence (`credibility.lock.json` in the folder of the model) instead of comparing the evidence with it; an evidence which differs from an existing lock is an error (see [Locked splits](credibility.md#locked-splits)). `--source-commit SHA` records the commit of the source repository in the provenance instead of the commit of the working directory, with the dirty state unknown; `pkpd_release` passes it, since it runs the assessment in the release repository. Without it the provenance is the commit of the repository of the model and whether the folder of the model has changes. A directory that holds a `credibility.json` is replaced, so no file of an earlier assessment survives (e.g. of a context of use that was removed); a non-empty directory without one is refused and left as it is. An error of the declaration, of the evidence (a validation on fitted data, an evidence without a point), of the records (a split that differs from the kinds the adopted fits recorded, a record without the kinds of its mappings), of a `credibility.json` that is no assessment or of the directory is printed and ends the command with exit code 1. `--no-verification`, `--lock`, `--uncertainty`, `--variability` or `--source-commit` with `--from-json` is a usage error (exit code 2), and so are settings of the uncertainty or the variability which the kind refuses.

## pkpd_release

`pkpd_release` makes the complete release of the repository of a model (`hctz-model`) from its package with the release of sbmlsim: the build (the kind `release` of the run layer of sbmlsim: export, models, simulations, fit reports, credibility, release notes, docs) into a clone of the release repository in an output, then the publishing (commit, tag, push, Docker image and GitHub release). `--no-publish` stops before the commit, so the build can be reviewed; a later run without the flag publishes the reviewed build:

```bash
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
pkpd_release pkdb_models.models.hctz --version 0.8.0
pkpd_release pkdb_models.models.hctz --version 0.8.0 --repo ../hctz-model
```

`--repo` is the checkout of the release repository (by default a sibling of the monorepo named like the repository), which the build clones and only reads (the repository on GitHub is cloned when the checkout is missing) and which publishing fast-forwards to the release, `--out` the output of the build (by default `<repo>-release-<version>` next to the checkout, the release repository in its folder `repository`), `--pkpdsim-path` a local pkpdsim for a run with `--no-publish`, `--workers` the workers of the simulations. A refused check or a failed step ends the command with exit code 1 and the reason. The release metadata of a model, the one-time setup and the steps are in [Release of a model repository](release.md).

## Profiling a command

`--profile` anywhere after the model package runs a command under the statistical sampling profiler of python 3.15 (`profiling.sampling`, [PEP 799](https://peps.python.org/pep-0799/)), which samples the stack of the command from a process of its own, so the command runs at nearly its normal speed:

```bash
run_simulation pkdb_models.models.hctz --profile --action simulate --experiments all
run_fit pkdb_models.models.hctz --profile --subset PK
```

The profiler is the one of sbmlsim (`sbmlsim.run.profiling.run_profiled`, `--profile` of `python -m sbmlsim.run`), which runs `python -m pkpdsim.cli`. The profile is written into `profiling/<command>_<time>/` of the results of the model: `flamegraph.html` (the interactive flame graph), `heatmap/` (the sampled lines of the source code), `profile.pstats` (the statistics of the functions, for `pstats` and its viewers), `summary.txt` (the functions with the most time) and `profile.bin` (the samples; `python -m profiling.sampling replay --diff-flamegraph old/profile.bin new/profile.bin` compares two profiles). The command keeps its exit code. With an older python `--profile` is a usage error and the command does not run. The profiler runs the command as `python -m pkpdsim.cli <command> <package> ...`, which works without the profiler as well. A command which runs its simulations in a pool of processes (`--workers`) is profiled in its main process only.

## The folders of the results

Every command writes into the folders of the model object unless an option names another one; the results are not tracked, the adopted fits are:

- `run_simulation`: `MODEL.simulation_dir` (`results/simulation`), with `--uncertainty` `MODEL.uncertainty_dir` (`results/simulation_uncertainty`), `--rebuild` `results/factory`; `--results-dir` replaces `results`.
- `run_fit`: `MODEL.fit_dir` (`results/fit/<fit id>`), `--output_dir` replaces it.
- `run_fit --adopt`: `MODEL.fits_dir` (`fits` next to the output of the factory).
- `run_report`: `MODEL.report_dir` (`results/report/report`, `--name` names the last folder), `--output_dir` replaces it; with `--adopted` one folder per adopted fit below it.
- `run_identifiability`: `MODEL.identifiability_dir` (`results/identifiability/identifiability`, `--name` names the last folder), `--output_dir` replaces it.
- `run_credibility`: `MODEL.credibility_dir` (`results/credibility`), `--output_dir` replaces it.

`results` is `MODEL.results_dir`, the folder `results` of the model package unless the model object sets `results_path`; so every result stays in the folder of its model.

## A model package of its own

Since the commands come from pkpdsim, a model which is packaged on its own needs no `[project.scripts]`: it depends on pkpdsim and defines its model object, and once it is installed the commands run it by the name of its package:

```bash
pip install -e .                                # the package hctz_model with MODEL
run_simulation hctz_model --action simulate --experiments all
run_fit hctz_model --subset PK
```

## From python

The commands are functions of a model object as well, which take the arguments after the model: `run_command(model, args)`, `fit_command` (which returns the folder of every fit by its id), `report_command` and `identifiability_command` of `pkpdsim.cli`. The model is the model object `MODEL` of its package, which the run layer loads by its reference; a copy of it with other fields is a `ValueError`. The console scripts are `run_simulation_main(argv)`, `run_fit_main`, `run_report_main` and `run_identifiability_main`, which take the model package first. In a script or the `__main__` block of an experiment `MODEL.run_experiments` runs experiments directly:

```python
# not executed
from pkdb_models.models.hctz import MODEL
from pkdb_models.models.hctz.experiments.studies import Anderson1961

MODEL.run_experiments(Anderson1961, output_dir=Anderson1961.__name__)
```
