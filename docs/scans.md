# Scans and sensitivity analyses

Two modules hold what the scripts of the models around scans and sensitivity analyses share, up to the runner of the sensitivity analyses of a script. Both build on the API of sbmlsim 0.9, where a design of `sbmlsim.simulation.sampling` is a dimension of a `Scan` and its result is a `ScanResult`.

## The timecourse of a point of a scan

sbmlsim 0.9 returns a scan as a `ScanResult` with the layout `(*dims, time)`, and at the time of a change after the start (a second dose) it keeps the state before and after the change. The pharmacokinetic parameters of a scan (`<model>_pk.py`) and the figures of a model read the timecourse of one point as sbmlsim 0.8 gave it, with the state after the change only, so that their numbers and figures stay the same. `pkpdsim.scan.single_dimension(res)` gives the id of the dimension of a scan of one dimension and `point_timecourse(res, key, /, **point)` the time points and the values of a variable at the point `point` (the position along every dimension, e.g. `dim_dose=2`; `res` and `key` are positional only, so a dimension may have any name); a result without dimension needs none:

```python
# not executed
from pkpdsim.scan import point_timecourse, single_dimension

dim = single_dimension(res)
for k in range(res.ds.sizes[dim]):
    time, values = point_timecourse(res, "[Cve_hctz]", **{dim: k})
```

A scan with several dimensions, a point which misses a dimension or names an unknown one, and a variable which is no timecourse raise a `ValueError`. The padding of a result whose simulations have different times is removed.

## The parameters of a sensitivity analysis

sbmlsim 0.9 analyses the sensitivity of a model on the result of a scan, with a design (`local`, `sobol`, `fast`, `lhs`) as a dimension of the scan; it has no table of the parameters of an analysis any more, which the sensitivity scripts of the models write next to their results. `pkpdsim.sensitivity` is that table with the behavior of sbmlsim 0.8, so the tables of the scripts stay the same, and the distributions the global designs need:

```python
# not executed
from pkpdsim.sensitivity import (
    ParameterType,
    distributions,
    free_parameters,
    parameters_from_sbml,
    sample_size,
    set_bounds,
    set_relative_bounds,
    write_parameters,
)
from sbmlsim.simulation.sampling import sobol

parameters = parameters_from_sbml(MODEL_PATH, exclude_ids={"Mr_hctz"})
set_bounds(parameters, [("f_renal_function", 0.85, 1.15, ParameterType.SCALING)])
set_relative_bounds(parameters, 0.15)
write_parameters(parameters, MODEL.results_dir / "sensitivity")
# a design must not vary what the simulation of a group sets
free = free_parameters(parameters, [simulation])
design = sobol(distributions(free), n=sample_size(4096), seed=1234, model=MODEL_PATH)
```

- `parameters_from_sbml` takes the constant parameters, compartments and species of a model in its order, with their values after the initialization, and leaves out those without a value, those close to zero and `exclude_ids`.
- `set_bounds` sets the bounds and the type of parameters by id (an unknown id is logged and skipped), `set_relative_bounds` lets every parameter without bounds vary by `value * (1 +- fraction)`. For a negative value these bounds are reversed (the lower bound `value * (1 - fraction)` is above the upper one), as in the table of sbmlsim 0.8, which stays the same; such a parameter gets its bounds with `set_bounds` before a global analysis.
- `write_parameters` writes `parameters.tsv` and the LaTeX table `parameters.tex` (`jinja2` is needed for the table), `parameters_to_df` gives the table as a `DataFrame` with the columns `uid`, `name`, `value`, `lower_bound`, `upper_bound`, `unit`, `type` and `reference`, sorted by type and id.
- `distributions` gives the uniform distribution of every parameter in its bounds, the input of a global design; a parameter without bounds or with a lower bound above the upper one raises a `ValueError` naming it.
- `free_parameters(parameters, simulations)` gives the parameters which none of the simulations sets, in their order, in its `preinit_changes`, in the `preinit_changes` of its presimulation or in a `Change`. A design replaces its targets wherever a simulation sets them, so a parameter of a group (e.g. `f_renal_function` of the degree of renal impairment) would take the values of the design in every group and the groups would be one analysis; sbmlsim 0.8.5 left such a parameter out ("Removing fixed parameter"), sbmlsim 0.9's designs do not. The ids are compared as written, plain ids as in sbmlsim 0.8.5: a target given as a selection (`[S]`) does not match the parameter `S`. Every parameter left out is logged at level info.
- `sample_size(n)` is `n`, or the value of the environment variable `PKPDSIM_SENSITIVITY_N` for a quick run of the scripts in a verification.

The distributions do not need a model:

```python
from pkpdsim.sensitivity import (
    ParameterType,
    SensitivityParameter,
    distributions,
    parameters_to_df,
    set_bounds,
    set_relative_bounds,
)

parameters = [
    SensitivityParameter(uid="Ka", name="absorption", value=2.0, unit="1/hr"),
    SensitivityParameter(uid="f_renal_function", name="renal function", value=1.0),
]
set_bounds(parameters, [("f_renal_function", 0.85, 1.15, ParameterType.SCALING)])
set_relative_bounds(parameters, 0.5)
print(parameters_to_df(parameters)[["uid", "lower_bound", "upper_bound", "type"]])
print(distributions(parameters))
```

A parameter which a simulation sets is left out of the design:

```python
from sbmlsim import Q
from sbmlsim.simulation import Simulation

from pkpdsim.sensitivity import SensitivityParameter, free_parameters

parameters = [
    SensitivityParameter(uid="Ka", name="absorption", value=2.0),
    SensitivityParameter(uid="f_renal_function", name="renal function", value=1.0),
]
group = Simulation(
    end=24,
    steps=100,
    preinit_changes={"f_renal_function": 0.5, "PODOSE": Q(10, "mg")},
)
print([p.uid for p in free_parameters(parameters, [group])])
```

## The sensitivity analyses of a script

The sensitivity script of a model (`<model>/sensitivity/sensitivity_analysis.py`) runs the analyses of sbmlsim 0.9 over the simulations of its groups and writes the tables and figures of its outputs. What every script does is `SensitivityRunner` of `pkpdsim.sensitivity`; a script gives only what is the model's own: the parameters with their bounds, the outputs with their observables and the groups of its analyses.

```python
# not executed
from pkpdsim.sensitivity import (
    FastAnalysis,
    LocalAnalysis,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityRunner,
    SobolAnalysis,
)
from sbmlsim import Q
from sbmlsim.simulation import PK, Formula, Simulation

from pkdb_models.models.hctz import MODEL
from pkdb_models.models.hctz.experiments.base_experiment import (
    HCTZSimulationExperiment,
)

OUTPUTS = [
    SensitivityOutput("cmax", "HCTZ Cmax"),
    SensitivityOutput("max_diuresis", "max diuresis"),
]
OBSERVABLES = [
    PK("pk_hctz", "[Cve_hctz]", dose="PODOSE_hctz", route="oral", parameters=["cmax"]),
    Formula("cmax", "pk_hctz.cmax", unit="mM"),
    Formula("max_diuresis", "max(diuresis)", unit="l/min"),
]


def sensitivity_simulation(changes: dict[str, float]) -> Simulation:
    return Simulation(
        end=5 * 24 * 60,
        steps=10000,
        time_unit="min",
        preinit_changes={"PODOSE_hctz": Q(25, "mg"), **changes},
    )


renal_map = HCTZSimulationExperiment.renal_map
GROUPS = [
    SensitivityGroup("control", "Control", sensitivity_simulation({})),
    *(
        SensitivityGroup(
            uid, name, sensitivity_simulation({"f_renal_function": renal_map[name]})
        )
        for uid, name in [
            ("mildRI", "Mild renal impairment"),
            ("modRI", "Moderate renal impairment"),
            ("sevRI", "Severe renal impairment"),
        ]
    ),
]

SENSITIVITY = SensitivityRunner(
    model=MODEL,
    parameters=sensitivity_parameters,  # the table of the section above
    outputs=OUTPUTS,
    observables=OBSERVABLES,
    analyses={
        "local": LocalAnalysis(groups=GROUPS[:1]),
        "sampling": SamplingAnalysis(groups=GROUPS),
        "fast": FastAnalysis(groups=GROUPS[:1]),
        "sobol": SobolAnalysis(groups=GROUPS[:1]),
    },
)

if __name__ == "__main__":
    SENSITIVITY.main()
```

```bash
python -m pkdb_models.models.hctz.sensitivity.sensitivity_analysis             # the default analyses, local
python -m pkdb_models.models.hctz.sensitivity.sensitivity_analysis local sobol
PKPDSIM_SENSITIVITY_N=65 python -m pkdb_models.models.hctz.sensitivity.sensitivity_analysis fast
```

- **The run.** `main(argv)` reads the names of the analyses from the command line (`default`, `("local",)` unless the script gives others, without one; an unknown name exits with 2) and calls `run(names)`. `run` writes the table of the parameters (`parameters.tsv`, `parameters.tex` of `write_parameters`) into `model.results_dir / "sensitivity"` and runs the analyses in their order, every one into `<analysis>/<group>` below it; `path` replaces the folder (e.g. a folder per simulation of a model with several, which makes a runner per simulation), `simulator` the settings of the simulator. The simulator has the tolerances of the model and `n_workers=None`: sbmlsim chooses the processes of the pool, so the script needs the guard `if __name__ == "__main__":`. An unknown name is a `ValueError` before anything runs.
- **The parameters** are a function, `parameters`, which a run calls (it reads the model). An analysis varies the parameters which no simulation of its groups sets, `free_parameters`: `sampling` over the four groups of hctz leaves out `f_renal_function`, which the groups of renal impairment set.
- **The outputs.** A `SensitivityOutput` is an observable of the simulation of a group: `uid` names it in the tables and figures, `name` in their titles, and `observable` is the id of the observable in the result, the uid when not given. An id which is no identifier (`[Cve_api]_cmax`, the id of the tables of sbmlsim 0.8) names an observable `Cve_api_cmax`. The ids and the observables of the outputs are unique, and every observable of an output is one of `observables`.
- **The groups.** A `SensitivityGroup` is a simulation with the changes of the group in its `preinit_changes`; `uid` is the folder, `name` in the titles of the figures. Every analysis has its own groups, at least one (the local and global analyses of most models run for the control only, `GROUPS[:1]`, the sampling for every group).
- **The seed** of the designs is `seed`, 1234 by default.

| analysis | design and indices | per group | defaults |
| --- | --- | --- | --- |
| `LocalAnalysis` | `sampling.local` with the values of the simulation of the group, `sensitivity.local` | `<output>.normalized.tsv`, the heatmap `normalized.png` | `delta=0.01` |
| `SamplingAnalysis` | `sampling.lhs`, one design for every group | `samples.tsv` (the outputs and the parameters of every point), `summary.tsv` (a row per output: `output`, `name`, the statistics of `ScanResult.summary` with the `quantiles`, `unit`), `distributions.png` | `n=1000`, `quantiles=(0.025, 0.05, 0.5, 0.95, 0.975)` |
| `FastAnalysis` | `sampling.fast`, `sensitivity.fast` | `<output>.S1.tsv`, `<output>.ST.tsv`, `<output>.png` | `n=1000`, more than 64 |
| `SobolAnalysis` | `sampling.sobol`, `sensitivity.sobol` | `<output>.S1.tsv`, `<output>.ST.tsv`, `<output>.png` | `n=4096`, a power of 2 |

The number of samples `n` goes through `sample_size`, so `PKPDSIM_SENSITIVITY_N` replaces it in a smoke run. An analysis without a group or with two groups of one id, a number of samples or a `delta` which is not positive, and a runner without outputs or analyses, with two outputs of one id or one observable, with an output whose observable is none of `observables` or with a default which is none of its analyses are a `ValueError` when they are created, not after the first analysis ran.

The outputs and the analyses do not need a model:

```python
from sbmlsim.simulation import Simulation

from pkpdsim.sensitivity import LocalAnalysis, SensitivityGroup, SensitivityOutput

print(SensitivityOutput("cmax", "Cmax").observable)
print(SensitivityOutput("[Cve_api]_cmax", "API Cmax", observable="Cve_api_cmax"))
control = SensitivityGroup("control", "Control", Simulation(end=24, steps=100))
print(LocalAnalysis(groups=[control]).delta)
```
