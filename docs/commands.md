# Commands

pkpdsim installs five commands which work with every model package. The model is the first argument of a command, the dotted name of its package; everything after it are the options of the command for that model:

- `run_simulation` creates the models with the factory, runs the simulation experiments and the extra actions of the model.
- `run_fit` fits the parameters of a fit problem, analyses their identifiability, or adopts a finished fit.
- `run_report` writes the report of a fit for stored parameter sets.
- `run_identifiability` analyses the identifiability of a stored parameter set.
- `pkpd_release` makes the complete release of the repository of a model.

```bash
run_simulation pkdb_models.models.hctz --action simulate --experiments all
run_fit pkdb_models.models.hctz --subset PK --runs 10 --cores 10
run_report pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
```

A model package registers no commands of its own: the commands are defined once, in pkpdsim, and everything which differs between two models comes from the model object `MODEL` of the package (its actions, fit problems, settings and folders) or from the options of the command line.

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

- `factory` creates the models with the factory of the model: the SBML, the ODEs and the COMBINE archive.
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

- `--experiments` (`-e`): comma separated groups and experiments of the registry, e.g. `all`, `misc,Anderson1961`. A single group is written into the folder of the group, any other selection into `custom_selection`, both in the `simulation` folder of the results.
- `--results-dir` (`-r`): the folder of the results, `MODEL.results_dir` by default.
- `--formats` (`-f`): comma separated formats of the figures, the static images of matplotlib (`svg`, `png`, ...) and the interactive pages of plotly (`html`); the formats of the `FigureSettings` of the model by default. The report shows the SVG and links the interactive page; without `svg` the interactive page is embedded in the report.
- `--open` (`-o`): opens the report of the simulations in the web browser.
- `--workers` (`-w`): the processes of the scans for this run: a number, `auto` (sbmlsim's choice, a pool of every CPU for a scan of 256 points or more) or `model` (the `n_workers` of the model object, the default). With serial scans (`--workers 1`, the default of a model) sbmlsim runs the experiments in a pool, a worker per four experiments, at most one per CPU; with a pool of the scans the experiments run one after the other.
- `--uncertainty` (`-u`) with `--draws` (100), `--seed` (1234) and `--quantiles` (0.05 0.95): runs the experiments over draws of the fitted parameters of the model and draws the curves of the simulations as bands, into the folder `simulation_uncertainty` of the results, see [Uncertainty of the fitted parameters](uncertainty.md).

A run which fails is not the end of the command. An experiment which cannot be initialized or whose figure fails is reported, the other experiments run, the report and the figures are written, and the command then runs its remaining actions, prints every failure and exits with 1.

A model without experiments (`experiments=None`) offers only `--action factory`, `--action all` (which runs only the factory and the extra actions with `in_all`) and its extra actions, and none of the options of the experiments; `simulate` and `list_experiments` are a usage error for it.

### Extra actions

A model adds actions of its own with `ExtraAction` and their options with `extra_arguments` of the model object (see [Models and their commands](models.md#commands)), e.g. analyses which are not simulation experiments:

```bash
run_simulation pkdb_models.models.glimepiride --action analyze_cyp2c9
```

An extra action gets the options of the command line, the settings of the simulator of the run as `options.simulator` (the `SimulatorSettings` of the model with the processes of `--workers`) and the settings of `--uncertainty` as `options.uncertainty` (`None` without). An action which simulates passes them on, `model.run_experiments(..., simulator=options.simulator)` or `pkpdsim.run.create_simulator(options.simulator)` for a runner of its own, otherwise `--workers` does not reach it. An option of the model named `--simulator`, or `--uncertainty` for a model without experiments, is refused with a `ValueError` before the action runs, because it would be replaced.

## run_fit

`run_fit` fits the parameters of a fit problem of the model (`--subset`, the first fit definition by default) and, with `--identifiability`, analyses the identifiability of the fitted parameters: the profile likelihood of every fitted parameter and the Fisher information around the best parameter set, which end up in the identifiability section of the report of the fit. The analysis is the default of a model with `fit_identifiability`, which `--no-identifiability` skips.

```bash
run_fit pkdb_models.models.hctz --subset PK --runs 10 --cores 10 --max-points 30
run_fit pkdb_models.models.hctz --subset PK --method DE --strategy SINGLE
run_fit pkdb_models.models.hctz --adopt PK_<fit id>
run_fit pkdb_models.models.inulin --subset PK --runs 1 --start profile_point.json
```

The options of the fit are those of sbmlsim (`sbmlsim.fit.cli.fit_cli`): `--subset` (`-x`), `--name` (`-n`), `--runs` (`-r`), `--cores` (`-c`), `--seed` (`-s`), `--method` (`-m`, `LSQ` least squares with latin hypercube starts or `DE` differential evolution), `--strategy` (`-t`, `ALL` the experiments together or `SINGLE` every experiment on its own), `--experiments` (`-e`), `--timeout`, `--start`, `--no-mapping-figures` and `--output_dir` (`-o`, the `fit` folder of the results by default). The options of the identifiability are `--parameter` (`-p`, repeatable), `--alpha` (`-a`), `--df`, `--max-points`, `--no-reoptimize` and `--no-fisher`.

Least squares starts every run from a sample of the bounds of the parameters (log-uniform latin hypercube sampling, `sbmlsim.fit.sampling`); the start value of a parameter (`FitParameter.start_value`) is neither the first start nor the center of the samples, it is used only for a parameter with an infinite bound, which is not sampled. `--start` starts least squares from a given point instead, e.g. a point of a profile likelihood below the cost of the fit: a json file of the values by parameter id in the units of the fit parameters (`{"kslow_inu": 0.002}`), or the id of a fit, an adopted fit of the model (`MODEL.fits_dir`) or else a finished fit in the output folder, whose best parameters are converted into the units of the fit parameters. The point replaces the start values of the parameters it gives (`pkpdsim.fit.start_point`, `with_start`, on copies of the parameters), and every run starts from the start values (the sampling `START` of sbmlsim), so a fit with `--start` has one run (`--runs 1`); a parameter which the start does not give starts from its start value, a parameter of the file which is not fitted, a value which is no number or outside of the bounds, `--runs` other than 1 and differential evolution, which draws its population from the bounds, are usage errors.

The problem of the fit, of its identifiability analysis and of its report loads the model at the optimization level of the model (`SimulatorSettings.jit_optimization`), as `run_report` and `run_identifiability` do: sbmlsim builds a fit problem without a level, so the commands hand sbmlsim the fit definitions at the level (`pkpdsim.fit.at_level`, `leveled_definitions`), whose problems carry it into every worker of a parallel fit; the environment of the process and of a forkserver which starts the workers does not matter, and nothing is set in it. The selection of the mappings (`select_mappings`, which instantiates the experiments once in the process) is no problem of the fit and keeps the level of sbmlsim.

The report of the fit carries the title of the model (`PKPDModel.title`), the description of the fit (`FitDefinition.description`) and of its parameters (`FitParameter.description`), and a note on the adopted fits of the other subsets of the model, whose parameters were fixed at their fitted values. The provenance of the fit records the arguments of the command (`subset`, `runs`, `cores`, `method`, `command`, `start` and `fixed_fits`, the ids of the fixed adopted fits), the start values of the parameters and the fixed values, i.e. the values of the parameters of the fixed adopted fits which the fit does not optimize itself.

`--adopt <fit id>` copies the record of the finished fit `<fit id>` from the output folder into the adopted fits of the model (`MODEL.fits_dir`) instead of fitting, `--force` replaces an adopted record of the same id; the id then goes into `fits` of the model object, see [Uncertainty of the fitted parameters](uncertainty.md).

## run_report and run_identifiability

Both work on the `parameters.json` of an earlier fit (`sbmlsim.fit.cli.report_cli`, `identifiability_cli`):

```bash
run_report pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.hctz --subset PK --max-points 30 results/fit/PK_<fit id>/parameters.json
```

`run_report` writes the report of the fit problem for one or more files of parameter sets, `run_identifiability` the profile likelihood and the Fisher information of the first parameter set of a file. `--subset` is the fit problem, `--output_dir` the folder of the report.

The report of `run_report` has the context of the model as the report of `run_fit` has it: the title of the model, the description of the fit and of its parameters and the note on the fixed adopted fits, which are the fits the provenance of the record next to the first file names, not the adopted fits of the model now (a record without provenance names none).

`run_report <package> --adopted` writes the report of every adopted fit of the model (`MODEL.fits`) instead, each from its record in `MODEL.fits_dir`, without a parameters file and without running a fit:

```bash
run_report pkdb_models.models.hctz --adopted
```

Every fit has a report of its own in `MODEL.report_dir/<fit id>` (`--output_dir` replaces `MODEL.report_dir`). It is the report which `run_report <package> --subset <subset> --name <fit id> <record>/parameters.json` writes for the record: the complete report of sbmlsim with the title of the model, the description of the fit and its parameters, the adopted fits whose parameters were fixed as the provenance of the record names them, and the provenance of the record (a record without it states that it was not recorded and names no fixed fits). The subset of a fit is the part of its id before the date, `PK` for `PK_20261016_185018__3efe3`, and must have a fit definition. There is no report which merges the fits of a model.

## pkpd_release

`pkpd_release` makes the complete release of the repository of a model (`hctz-model`) from its package: export, build of the models, simulations and fit reports, release notes, docs, commit, tag, push, Docker image and GitHub release. `--no-publish` stops before the commit, so the build can be reviewed; a later run without the flag publishes the reviewed build:

```bash
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
pkpd_release pkdb_models.models.hctz --version 0.8.0
pkpd_release pkdb_models.models.hctz --version 0.8.0 --repo ../hctz-model
```

`--repo` is the checkout of the release repository (by default a sibling of the monorepo named like the repository, cloned when missing), `--pkpdsim-path` a local pkpdsim for a run with `--no-publish`, `--workers` the workers of the simulations. A refused check or a failed step ends the command with exit code 1 and the reason. The release metadata of a model, the one-time setup and the steps are in [Release of a model repository](release.md).

## Profiling a command

`--profile` anywhere after the model package runs a command under the statistical sampling profiler of python 3.15 (`profiling.sampling`, [PEP 799](https://peps.python.org/pep-0799/)), which samples the stack of the command from a process of its own, so the command runs at nearly its normal speed:

```bash
run_simulation pkdb_models.models.hctz --profile --action simulate --experiments all
run_fit pkdb_models.models.hctz --profile --subset PK
```

The profile is written into `profiling/<command>_<time>/` of the results of the model: `flamegraph.html` (the interactive flame graph), `heatmap/` (the sampled lines of the source code), `profile.pstats` (the statistics of the functions, for `pstats` and its viewers), `summary.txt` (the functions with the most time) and `profile.bin` (the samples; `python -m profiling.sampling replay --diff-flamegraph old/profile.bin new/profile.bin` compares two profiles). The command keeps its exit code. With an older python `--profile` is a usage error and the command does not run. The profiler runs the command as `python -m pkpdsim.cli <command> <package> ...`, which works without the profiler as well. A command which runs its simulations in a pool of processes (`--workers`) is profiled in its main process only.

## The folders of the results

Every command writes into the folders of the model object unless an option names another one; the results are not tracked, the adopted fits are:

- `run_simulation`: `MODEL.simulation_dir` (`results/simulation`), with `--uncertainty` `MODEL.uncertainty_dir` (`results/simulation_uncertainty`); `--results-dir` replaces `results`.
- `run_fit`: `MODEL.fit_dir` (`results/fit/<fit id>`), `--output_dir` replaces it.
- `run_fit --adopt`: `MODEL.fits_dir` (`fits` next to the output of the factory).
- `run_report`: `MODEL.report_dir` (`results/report`), `--output_dir` replaces it; with `--adopted` one folder per adopted fit below it.
- `run_identifiability`: `MODEL.identifiability_dir` (`results/identifiability`), `--output_dir` replaces it.

`results` is `MODEL.results_dir`, the folder `results` of the model package unless the model object sets `results_path`; so every result stays in the folder of its model.

## A model package of its own

Since the commands come from pkpdsim, a model which is packaged on its own needs no `[project.scripts]`: it depends on pkpdsim and defines its model object, and once it is installed the commands run it by the name of its package:

```bash
pip install -e .                                # the package hctz_model with MODEL
run_simulation hctz_model --action simulate --experiments all
run_fit hctz_model --subset PK
```

## From python

The commands are functions of a model object as well, which take the arguments after the model: `run_command(model, args)`, `fit_command`, `report_command` and `identifiability_command` of `pkpdsim.cli`. The console scripts are `run_simulation_main(argv)`, `run_fit_main`, `run_report_main` and `run_identifiability_main`, which take the model package first. In a script or the `__main__` block of an experiment `MODEL.run_experiments` runs experiments directly:

```python
# not executed
from pkdb_models.models.hctz import MODEL
from pkdb_models.models.hctz.experiments.studies import Anderson1961

MODEL.run_experiments(Anderson1961, output_dir=Anderson1961.__name__)
```
