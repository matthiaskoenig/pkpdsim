# Models and their commands

A model package describes itself once, in its `__init__.py`, as a `PKPDModel`.
The commands, the verification and the tests of pkpdsim work with this object;
what stays in the package is what is the model's own: its factory, the registry
of its experiments, its fit definitions and parameters, the filters of its
mapping collections and its tolerances.

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
    data_paths=tuple(DATA_PATHS),
    archive_path=MODEL_BASE_PATH.parent / "hctz_model.omex",
    factory=".models.factory:create_hctz_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
    simulator=SimulatorSettings(absolute_tolerance=1e-10, relative_tolerance=1e-10),
    figures=FigureSettings(formats=("svg",), fig_dpi=300, legend_fontsize=10),
    fit_identifiability=True,
)
```

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

- `SimulatorSettings`: the absolute tolerance (one value, or an
  `AbsoluteTolerance` of sbmlsim per kind and state) and the relative tolerance
  of the experiments.
- `FigureSettings`: the formats of the figures (`svg`, `png`, `html`), the
  resolution of the static images (`savefig_dpi`) and the sizes and positions
  of sbmlsim's `Figure`. They hold for the process, so they are applied for a
  run and restored after it.
- `reduced_selections`, `save_results`, `fit_identifiability`: what a run
  simulates, whether it stores its results, whether a fit analyses the
  identifiability of its parameters by default.

## Commands

`pyproject.toml` of the repository of the model registers the commands with the
model object:

```toml
[project.scripts]
run_hctz = "pkdb_models.models.hctz:MODEL.run_main"
fit_hctz = "pkdb_models.models.hctz:MODEL.fit_main"
report_hctz = "pkdb_models.models.hctz:MODEL.report_main"
identifiability_hctz = "pkdb_models.models.hctz:MODEL.identifiability_main"
```

```bash
run_hctz --action factory                       # the SBML, the ODEs, the archive
run_hctz --action list_experiments              # the groups and experiments
run_hctz --action simulate --experiments "misc,Anderson1961" --open
run_hctz --action all --results-dir /path/to/results --formats svg,png
fit_hctz --subset=PK --runs=10 --cores=10 --max-points=30
report_hctz --subset=PK results/fit/PK_<fit id>/parameters.json
identifiability_hctz --subset=PK results/fit/PK_<fit id>/parameters.json
```

A model adds actions of its own to `run_<model>` with `ExtraAction` and their
options with `extra_arguments`, e.g. analyses which are not simulation
experiments.

In a script or the `__main__` block of an experiment:

```python
# not executed
from pkdb_models.models.hctz import MODEL
from pkdb_models.models.hctz.experiments.studies import Anderson1961

MODEL.run_experiments(Anderson1961, output_dir=Anderson1961.__name__)
```

## Tests

`pkpdsim.testing` holds the tests of a model as functions: the registry is
complete (`check_registry`), the experiments initialize
(`initialize_experiments`), all experiments run (`run_all_experiments`) and every
fit starts (`start_fits`).

## Verification of a change

`pkpdsim-verify` dumps the experiments of a model in two checkouts and compares
them, writes the fit definitions as json, and compares figures and SBML:

```bash
pkpdsim-verify dump --root ../base --model pkdb_models.models.hctz --out base.pkl --atol model
pkpdsim-verify dump --root . --model pkdb_models.models.hctz --out new.pkl --order base.pkl --atol model
pkpdsim-verify compare base.pkl new.pkl
pkpdsim-verify definitions --root . --model pkdb_models.models.hctz --out new.json
pkpdsim-verify png old/figures new/figures
pkpdsim-verify sbml old.xml new.xml
```
