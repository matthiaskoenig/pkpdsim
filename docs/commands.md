# Commands

pkpdsim installs four commands which work with every model package. The model
is the first argument of a command, the dotted name of its package; everything
after it are the options of the command for that model:

- `run_simulation` creates the models with the factory, runs the simulation
  experiments and the extra actions of the model.
- `run_fit` fits the parameters of a fit problem, analyses their
  identifiability, or adopts a finished fit.
- `run_report` writes the report of a fit for stored parameter sets.
- `run_identifiability` analyses the identifiability of a stored parameter set.

```bash
run_simulation pkdb_models.models.hctz --action simulate --experiments all
run_fit pkdb_models.models.hctz --subset PK --runs 10 --cores 10
run_report pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
```

A model package registers no commands of its own: the commands are defined once,
in pkpdsim, and everything which differs between two models comes from the
model object `MODEL` of the package (its actions, fit problems, settings and
folders) or from the options of the command line.

## The model argument

The model is given as the dotted name of its package, e.g.
`pkdb_models.models.hctz`. The command imports the package and takes its model
object `MODEL`, the `PKPDModel` of its `__init__.py` (see
[Models and their commands](models.md)). The package has to be importable in
the environment the command runs in, i.e. installed, e.g. the checkout of a
repository of models with `uv sync` or `pip install -e .`; a console script
does not import from the working directory.

The model comes first, since the options of a command depend on the model:
the actions of `run_simulation` include the extra actions of the model, and
the choices and defaults of the fit commands are its fit problems and its
folders. An option before the model is a usage error.

`--help` without a model prints the usage of the command, with a model the
options of that model:

```bash
run_fit --help                                  # the usage of the command
run_fit pkdb_models.models.hctz --help          # the options for hctz
```

A model which cannot be resolved stops the command with a usage error and the
exit code 2:

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

Any other error of the import of the package, e.g. a module which the model
imports and which is not installed, is raised with its traceback, so that a
broken model shows where it breaks.

`pkpdsim.cli.resolve_model(package)` is the resolution as a function: it
returns the model object or raises `ModelNotFoundError` with the messages of
the table.

## run_simulation

`run_simulation` performs one `--action` (`-a`):

- `factory` creates the models with the factory of the model: the SBML, the
  ODEs and the COMBINE archive.
- `list_experiments` lists the groups and experiments of the registry of the
  model.
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

- `--experiments` (`-e`): comma separated groups and experiments of the
  registry, e.g. `all`, `misc,Anderson1961`. A single group is written into
  the folder of the group, any other selection into `custom_selection`, both in
  the `simulation` folder of the results.
- `--results-dir` (`-r`): the folder of the results, `MODEL.results_dir` by
  default.
- `--formats` (`-f`): comma separated formats of the figures, the static images
  of matplotlib (`svg`, `png`, ...) and the interactive pages of plotly
  (`html`); the formats of the `FigureSettings` of the model by default. The
  report shows the SVG and links the interactive page; without `svg` the
  interactive page is embedded in the report.
- `--open` (`-o`): opens the report of the simulations in the web browser.
- `--workers` (`-w`): the processes of the scans for this run: a number, `auto`
  (sbmlsim's choice, a pool of every CPU for a scan of 256 points or more) or
  `model` (the `n_workers` of the model object, the default).
- `--uncertainty` (`-u`) with `--draws` (100), `--seed` (1234) and
  `--quantiles` (0.05 0.95): runs the experiments over draws of the fitted
  parameters of the model and draws the curves of the simulations as bands,
  into the folder `simulation_uncertainty` of the results, see
  [Uncertainty of the fitted parameters](uncertainty.md).

A run which fails is not the end of the command. An experiment which cannot be
initialized or whose figure fails is reported, the other experiments run, the
report and the figures are written, and the command then runs its remaining
actions, prints every failure and exits with 1.

A model without experiments (`experiments=None`) offers only `--action
factory`, `--action all` (which runs only the factory and the extra actions
with `in_all`) and its extra actions, and none of the options of the
experiments; `simulate` and `list_experiments` are a usage error for it.

### Extra actions

A model adds actions of its own with `ExtraAction` and their options with
`extra_arguments` of the model object (see
[Models and their commands](models.md#commands)),
e.g. analyses which are not simulation experiments:

```bash
run_simulation pkdb_models.models.glimepiride --action analyze_cyp2c9
```

An extra action gets the options of the command line, the settings of the
simulator of the run as `options.simulator` (the `SimulatorSettings` of the
model with the processes of `--workers`) and the settings of `--uncertainty`
as `options.uncertainty` (`None` without). An action which simulates passes
them on, `model.run_experiments(..., simulator=options.simulator)` or
`pkpdsim.run.create_simulator(options.simulator)` for a runner of its own,
otherwise `--workers` does not reach it. An option of the model named
`--simulator`, or `--uncertainty` for a model without experiments, is refused
with a `ValueError` before the action runs, because it would be replaced.

## run_fit

`run_fit` fits the parameters of a fit problem of the model (`--subset`, the
first fit definition by default) and, with `--identifiability`, analyses the
identifiability of the fitted parameters: the profile likelihood of every
fitted parameter and the Fisher information around the best parameter set,
which end up in the identifiability section of the report of the fit. The
analysis is the default of a model with `fit_identifiability`, which
`--no-identifiability` skips.

```bash
run_fit pkdb_models.models.hctz --subset PK --runs 10 --cores 10 --max-points 30
run_fit pkdb_models.models.hctz --subset PK --method DE --strategy SINGLE
run_fit pkdb_models.models.hctz --adopt PK_<fit id>
```

The options of the fit are those of sbmlsim (`sbmlsim.fit.cli.fit_cli`):
`--subset` (`-x`), `--name` (`-n`), `--runs` (`-r`), `--cores` (`-c`),
`--seed` (`-s`), `--method` (`-m`, `LSQ` least squares with latin hypercube
starts or `DE` differential evolution), `--strategy` (`-t`, `ALL` the
experiments together or `SINGLE` every experiment on its own), `--experiments`
(`-e`), `--timeout`, `--no-mapping-figures` and `--output_dir` (`-o`, the
`fit` folder of the results by default). The options of the identifiability
are `--parameter` (`-p`, repeatable), `--alpha` (`-a`), `--df`,
`--max-points`, `--no-reoptimize` and `--no-fisher`.

`--adopt <fit id>` copies the record of the finished fit `<fit id>` from the
output folder into the adopted fits of the model (`MODEL.fits_dir`) instead of
fitting, `--force` replaces an adopted record of the same id; the id then goes
into `fits` of the model object, see
[Uncertainty of the fitted parameters](uncertainty.md).

## run_report and run_identifiability

Both work on the `parameters.json` of an earlier fit
(`sbmlsim.fit.cli.report_cli`, `identifiability_cli`):

```bash
run_report pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.hctz --subset PK --max-points 30 results/fit/PK_<fit id>/parameters.json
```

`run_report` writes the report of the fit problem for one or more files of
parameter sets, `run_identifiability` the profile likelihood and the Fisher
information of the first parameter set of a file. `--subset` is the fit
problem, `--output_dir` the folder of the report.

## The folders of the results

Every command writes into the folders of the model object unless an option
names another one; the results are not tracked, the adopted fits are:

- `run_simulation`: `MODEL.simulation_dir` (`results/simulation`), with
  `--uncertainty` `MODEL.uncertainty_dir` (`results/simulation_uncertainty`);
  `--results-dir` replaces `results`.
- `run_fit`: `MODEL.fit_dir` (`results/fit/<fit id>`), `--output_dir` replaces
  it.
- `run_fit --adopt`: `MODEL.fits_dir` (`fits` next to the output of the
  factory).
- `run_report`: `MODEL.report_dir` (`results/report`), `--output_dir` replaces
  it.
- `run_identifiability`: `MODEL.identifiability_dir`
  (`results/identifiability`), `--output_dir` replaces it.

`results` is `MODEL.results_dir`, the folder `results` of the model package
unless the model object sets `results_path`; so every result stays in the
folder of its model.

## A model package of its own

Since the commands come from pkpdsim, a model which is packaged on its own needs
no `[project.scripts]`: it depends on pkpdsim and defines its model object, and
once it is installed the commands run it by the name of its package:

```bash
pip install -e .                                # the package hctz_model with MODEL
run_simulation hctz_model --action simulate --experiments all
run_fit hctz_model --subset PK
```

## From python

The commands are functions of a model object as well, which take the arguments
after the model: `run_command(model, args)`, `fit_command`, `report_command`
and `identifiability_command` of `pkpdsim.cli`. The console scripts are
`run_simulation_main(argv)`, `run_fit_main`, `run_report_main` and
`run_identifiability_main`, which take the model package first. In a script or
the `__main__` block of an experiment `MODEL.run_experiments` runs experiments
directly:

```python
# not executed
from pkdb_models.models.hctz import MODEL
from pkdb_models.models.hctz.experiments.studies import Anderson1961

MODEL.run_experiments(Anderson1961, output_dir=Anderson1961.__name__)
```
