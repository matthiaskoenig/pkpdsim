# Models and their commands

A model package describes itself once, in its `__init__.py`, as a `PKPDModel`.
The commands, the verification and the tests of pkpdsim work with this object;
what stays in the package is what is the model's own: its factory, the registry
of its experiments, its fit definitions and parameters, the filters of its
mapping collections and its tolerances. The base class of its experiments, the
vocabulary of their metadata and the filters of the fits are described in
[Experiments and their metadata](experiments.md).

```python
# not executed
from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings

HCTZ_PATH = Path(__file__).parent
MODEL_BASE_PATH = HCTZ_PATH / "models" / "results" / "models"
MODEL_PATH = MODEL_BASE_PATH / "hctz_body_flat.xml"
DATA_PATHS = [HCTZ_PATH / "data" / "hydrochlorothiazide"]

MODEL = PKPDModel(
    name="hctz",
    package=__name__,
    title="HCTZ PBPK/PD model",
    path=HCTZ_PATH,
    model_path=MODEL_PATH,
    data_paths=DATA_PATHS,
    archive_path=MODEL_BASE_PATH.parent / "hctz_model.omex",
    factory=".models.factory:create_hctz_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
    simulator=SimulatorSettings(absolute_tolerance=1e-10, relative_tolerance=1e-10),
    figures=FigureSettings(formats=("svg",), savefig_dpi=300, legend_fontsize=10),
    fit_identifiability=True,
)
```

`PKPDModel`, `FigureSettings`, `SimulatorSettings` and `ExtraAction` take
keyword arguments only, so the order of their fields is not API. A list of
`data_paths` is stored as a tuple; a model is hashable unless its absolute
tolerance is an `AbsoluteTolerance` of sbmlsim.

## References

`factory`, `experiments` and `fit_definitions` are references
`"module:attribute"`, resolved when a command uses them. The experiments of a
model import its `__init__.py` for the paths, so the model object cannot import
them; a module which starts with a dot is relative to `package`.

| field | the reference names |
| --- | --- |
| `factory` | a function which creates the models in a directory |
| `experiments` | the registry: groups of experiment classes by name, with the group `all` |
| `fit_definitions` | the fit problems by name (`sbmlsim.fit.cli.FitDefinition`) |

## Settings

- `SimulatorSettings`: the absolute tolerance (one value, or an `AbsoluteTolerance` of sbmlsim per kind and state), the relative tolerance and the processes of a scan, `n_workers`. A run creates an sbmlsim `Simulator` from them.
- `FigureSettings`: sbmlsim's `FigureSettings` (the resolution, the sizes of the panels, the fonts and the legends, with the defaults of sbmlsim) and the two fields of a model, the `formats` of the figures (`svg`, `png`, `html`) and the resolution of the static images (`savefig_dpi`). A run hands the settings of sbmlsim to every figure; an experiment (`figure_settings`) and a figure (`Figure(settings=...)`) override them. Nothing is applied to the process, so nothing has to be restored after a run.
- `reduced_selections`, `save_results`, `fit_identifiability`: what a run simulates, whether it stores its results, whether a fit analyses the identifiability of its parameters by default.
- `fits`, `fits_path`: the ids of the adopted fits, whose records in `fits_dir` (`fits` next to the output directory of the factory by default) give the fitted values and their uncertainty to the factory (`fit_records()`), see [Uncertainty of the fitted parameters](uncertainty.md).
- `release`: the release metadata of the model (`ReleaseSettings`: the repository, the Docker image, the authors, ...), which `pkpd_release` needs to make the release of its model repository, see [Release of a model repository](release.md); `None` (the default) for a model without a release repository.

The simulations run serially by default (`n_workers=1`): a pool needs the guard `if __name__ == "__main__":` in the script which runs the model (the processes of the pool import it again), and the memory of a run grows with every process. `n_workers=4` runs a pool of four processes for a scan, `n_workers=None` leaves the choice to sbmlsim (a pool of every CPU for a scan of 256 points or more, serially below):

```python
from pkpdsim.model import FigureSettings, SimulatorSettings

serial = SimulatorSettings(absolute_tolerance=1e-10, relative_tolerance=1e-10)
assert serial.n_workers == 1
pool = SimulatorSettings(n_workers=4)
sbmlsim_chooses = SimulatorSettings(n_workers=None)

try:
    SimulatorSettings(n_workers=0)
except ValueError as err:
    print(err)

figures = FigureSettings(formats=("svg", "png"), savefig_dpi=300, legend_fontsize=10)
print(figures.sbmlsim_settings().legend_fontsize)
```

## Commands

The commands of pkpdsim run every model by the dotted name of its package; the
package registers no commands of its own:

```bash
run_simulation pkdb_models.models.hctz --action simulate --experiments all
run_fit pkdb_models.models.hctz --subset PK --identifiability
run_report pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.hctz --subset PK results/fit/PK_<fit id>/parameters.json
```

The model object gives everything which differs between two models: the
actions of `run_simulation` (the extra actions of `extra_actions`, an
`ExtraAction` with `name`, `help`, the `function` called with the model, the
parsed arguments and the folder of the results, and `in_all` to run it with
`--action all`; their options are added by the function of `extra_arguments`),
the fit problems of the fit commands, the settings of the simulator and the
figures and the folders of the results. See [Commands](commands.md) for the
commands, their options and their errors.

## Tests

`pkpdsim.testing` holds the tests of a model as functions: the registry is
complete (`check_registry`), the model object belongs to its package and its
references resolve (`check_model`), the experiments initialize
(`initialize_experiments`), all experiments run (`run_all_experiments`) and every
fit starts (`start_fits`). A function which checks several things raises one
`AssertionError` which names every failure, not only the first.

The known exceptions of a model (experiments which are not registered, modules
which do not import, experiments which run one at a time) are its `checks`, see
[Experiments and their metadata](experiments.md#exceptions-of-the-tests).

`check_model` also checks that the name of the model is the last component of its
package (or that component ends with `_<name>`), that `model_path`, the factory
directory and `archive_path` lie below the path of the model, and that the data
paths exist.

## Verification of a change

`pkpdsim-verify` dumps the experiments of a model in two checkouts and compares
them, writes the fit definitions and the settings of a run as json, and compares
figures and SBML:

```bash
pkpdsim-verify dump --root ../base --model pkdb_models.models.hctz --out base.pkl --atol model
pkpdsim-verify dump --root . --model pkdb_models.models.hctz --out new.pkl --order base.pkl --atol model
pkpdsim-verify compare base.pkl new.pkl
pkpdsim-verify definitions --root . --model pkdb_models.models.hctz --out new.json
pkpdsim-verify settings --root ../base --model pkdb_models.models.hctz --out base.json
pkpdsim-verify settings --root . --model pkdb_models.models.hctz --out new.json
diff base.json new.json
pkpdsim-verify png old/results/simulation new/results/simulation
pkpdsim-verify sbml old.xml new.xml
```

`--atol model` and `--rtol model` take the tolerances of the model object; any
other value which is not a number exits with 2. A package without a model object
exits with 2, and `compare --renames FILE` applies the expected renames of a
migration to the old dump (see
[Experiments and their metadata](experiments.md#renames-in-a-verification) and
[Verification across versions of sbmlsim](#verification-across-versions-of-sbmlsim)).
An experiment which cannot be initialized is recorded in the dump as an error and
the command exits with 1 after the other experiments are stored.

`settings` captures what a run of the experiments hands to sbmlsim, bound to the
signatures: the arguments of `Simulator` (its tolerances and `n_workers`), of
`ExperimentRunner` (without the simulator) and of the run, with its
`figure_settings` (without the output path), the matplotlib `rc` settings a run
changes and the counts of the calls. It exits with 2 for a run which did not
complete. The json has `run.figure_settings`; the part `figure` of pkpdsim 0.4 (the
class attributes of sbmlsim's `Figure`) is gone.

## Verification across versions of sbmlsim

A dump records the versions of sbmlsim and pkpdsim it was made with, and `compare`
prints them. A dump of pkpdsim 0.4 (sbmlsim 0.8) and one of pkpdsim 0.5 (sbmlsim
0.9) differ in their layout, which `compare` takes into account:

- sbmlsim 0.8 puts the time first, `(_point, *dims)`, sbmlsim 0.9 last,
  `(*dims, time)`. `compare` moves the time of every array first and matches the
  other dimensions by name. A dimension of size 1 which only one array has is
  squeezed (a scan of one point is a simulation now), different names of the
  dimensions are an infinite difference.
- sbmlsim 0.9 keeps the state before and after a change at the time of the change
  (a second dose), sbmlsim 0.8 only the state after it. Two dumps with the same
  time vectors are compared row by row. Otherwise the rows of a repeated time are
  matched from the last one, the k-th last old row with the k-th last new row, so
  the state after a change is compared with the state after it, and only the rows
  which one dump has more of at a time are left out (the states before a change
  of a dump of sbmlsim 0.9 compared with one of sbmlsim 0.8; two dumps of sbmlsim
  0.9 with two rows at a time compare both): the summary counts them and names
  the experiments with `changes: <experiment> <rows> rows`.
- The rows of every point of a scan are matched by the times of that point. The
  points of a grid share their times; in the ragged layout of sbmlsim 0.9
  (`_point`, every point with its own times, padded with NaN) and in a scan of
  sbmlsim 0.8 a point is matched and its times are compared on its own.
- A row at a time which only the old dump has is always a difference, and no
  pattern of `expected` hides it. A row at a time which only the new dump has (a
  new output time of a task) is a difference as well, which is printed per task
  (`times of <experiment>/<task>: <n> rows only in new, at <times>`) and which
  the pattern `experiment/task/time` of `expected` names.

The specification of a figure (its attributes as json) changed as well. A figure which the old dump has from sbmlsim 0.8 and the new dump from sbmlsim 0.9 is compared in the form of sbmlsim 0.9, and `compare` prints how many (`figures across versions of sbmlsim: <n>`):

- The settings. sbmlsim 0.8 kept a setting which a figure set (`fig.legend_position = "outside"`) as an attribute of the figure, sbmlsim 0.9 keeps it in `_settings` (`fig.settings = {"legend_position": "outside"}`). The attributes of the old figure which are fields of `FigureSettings` move into its `_settings`, which are compared with the settings the new figure is rendered with apart from the settings of the run: the `figure_settings` of its experiment overridden by its own `settings`. Another value, a setting which only one figure has and one which moved into the experiment are compared by value.
- The size. sbmlsim 0.8 computed the size of a figure when it was created (`_width` and `_height`: the explicit size, or the number of panels times the panel size), sbmlsim 0.9 keeps only the explicit size (`width` and `height`, `None` for the size of the panels) and computes the size when it renders the figure. The old size is compared with the size the new figure is rendered with apart from the settings of the run, which a dump of pkpdsim 0.5 stores with the settings (`__resolved__` of the specification, `pkpdsim.verify.figure_spec`), so an explicit size and the size of the panels are compared by value. A plot which is added to a figure gets its share of the size of the figure, from sbmlsim 0.8 always and from sbmlsim 0.9 only for an explicit size; such a size of a plot is left out on both sides, a size of its own is compared. A figure whose size was set after its plots were added gave them its earlier size in sbmlsim 0.8, which then counts as a size of their own and differs. A size computed from a setting which sbmlsim 0.8 read from the class `Figure` (`Figure.panel_width`; a dump of pkpdsim 0.4.1 ran without the settings of the model and saw sbmlsim's default) and which sbmlsim 0.9 reads from the figure settings of the model differs whenever the model changes that setting; the figure is the same.
- The fields which only sbmlsim 0.9 has, the selection of a `Data` (`sel`), the dimensions of the lines of a `Curve` (`over`) and the bands of a `Plot` (`_bands`), have no counterpart in sbmlsim 0.8: they are left out of the comparison, and every one which is not empty is listed for the review of the migration, `new fields: <experiment>/<figure>: <path> = <value>` (the path as in a difference, e.g. `subplots[0].plot._curves[1].y.sel`).

Two figures of one version of sbmlsim are compared as they are, `sel`, `over` and `_bands` included. Neither dump has the settings of the run (the `FigureSettings` of the model, which a run hands to every figure); `settings` compares them.

The ids of sbmlsim 0.9 are SIds (see [Keys](experiments.md#keys)). The migration says so in the renames file instead of weakening the comparison, which then passes with exit code 0 or fails:

| part | content |
| --- | --- |
| `ids` | an old id to its new id: the task of a result, the dataset, the figure and mapping keys and the references to them in a figure |
| `sids` | `true`: every id of the old dump becomes `pkpdsim.experiment.to_sid(id)` (after `ids`) |
| `figure_fields` | the keys of the specification of a figure which are left out of both dumps, at the level of the figure only (not in its plots, curves or axes), in the form in which the figures are compared: a field of a figure which `compare` cannot map. The specification is the attributes of sbmlsim's `Figure` and its objects as json (`pkpdsim.verify.serialize`, `figure_spec`), so the fields are those of the version of sbmlsim which made the dump. The port from sbmlsim 0.8 to 0.9 needs none |
| `expected` | a part (`results`, `datasets`, `mappings`, `figures`, `ranges`) to the patterns of the keys of the differences which are expected |
| `metadata`, `strings` | as in [Experiments and their metadata](experiments.md#renames-in-a-verification) |

```json
{
  "sids": true,
  "expected": {
    "results": ["Anderson1961/task_po_2p5/time"],
    "mappings": ["Anderson1961/Fig1_hctz"]
  }
}
```

The patterns of `expected` are `fnmatch` patterns of the keys:

| part | key |
| --- | --- |
| `results` | `experiment/task/variable` |
| `datasets` | `experiment/dataset/column` |
| `mappings`, `ranges` | `experiment/mapping` |
| `figures` | `experiment/figure` |

`*` and `?` are wildcards, a bracket is literal: `Anderson1961/task_po_2p5/[Cve_hctz]` names the variable `[Cve_hctz]` of the task `task_po_2p5` (the base names the task of a simulation `task_<simulation>`). A difference which matches is printed as expected and not counted, a pattern which matches no difference is counted, so the file says exactly what changes. Two old ids with one new id, a `strings` rename which merges two keys of a part and two mappings of one experiment with one key are errors (exit code 2). The renames are applied to the content of the parts of the old dump, never to the names of the parts, and the loaded dump is not changed.

`png` compares the figures of two directories by their file names. After the
renames to SIds the old file `A_fig 2.5.png` is `A_fig_2p5.png`: `png --sids`
pairs an old figure with the new figure of the SId of its name (the stem through
`to_sid`) and exits with 2 when two old figures have one SId:

```bash
pkpdsim-verify png --sids old/results/simulation new/results/simulation
```
