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

### The state after a change of the whole result

A figure which reads the arrays of a result (a band over the points of a scan, the result of a task of an experiment) needs the state after a change for all points at once. `state_after_changes(res)` gives the result on its time points without the repeated ones, with the state after a change at its time (`res.interpolate(np.unique(res.time_points()))`), as a `ScanResult` whose arrays keep their dimensions; for one point it has the time points and values of `point_timecourse`. A result without a change is unchanged.

```python
# not executed
from pkpdsim.scan import state_after_changes

res = state_after_changes(self.results[self.to_sid("task_sim")])
```

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

- **The run.** `main(argv)` reads the names of the analyses from the command line (`default`, `("local",)` unless the script gives others, without one; an unknown or repeated name exits with 2) and calls `run(names)`. `run` writes the table of the parameters (`parameters.tsv`, `parameters.tex` of `write_parameters`) into `model.results_dir / "sensitivity"` and runs the analyses in their order, every one into `<analysis>` below it (see the table); `path` replaces the folder, `subfolder` puts the analyses into a folder below it while the table stays in `path` (see below), `simulator` replaces the settings of the simulator. The simulator has the tolerances of the model and `n_workers=None`: sbmlsim chooses the processes of the pool, so the script needs the guard `if __name__ == "__main__":`.
- **The workers.** A pool needs about 300 to 340 MB per process for a PBPK model, so every CPU of a 20-CPU machine exceeds a memory budget of 8 GB. `--workers N` (an int of at least 1; `0` or absent leaves the setting of the simulator, `n_workers=None`) caps the processes: `python -m ...sensitivity_analysis sobol --workers 8`. `main` passes it as the settings of the model with `n_workers=N`; a script with its own `main` reads `parse_command_line(...).workers` and passes `simulator=runner.simulator_settings(words.workers)` to `run` or `run_subfolders`.
- **The parameters** are a function, `parameters`, which a run calls (it reads the model). An analysis varies the parameters which no simulation of its groups sets, `free_parameters`: `sampling` over the four groups of hctz leaves out `f_renal_function`, which the groups of renal impairment set.
- **The outputs.** A `SensitivityOutput` is an observable of the simulation of a group: `uid` names it in the tables and the names of the files, `name` in the titles and the heatmaps of the figures and in the column `name` of the tables, and `observable` is the id of the observable in the result, the uid when not given. An id which is no identifier (`[Cve_api]_cmax`, the id of the tables of sbmlsim 0.8) names an observable `Cve_api_cmax`. The ids, the names and the observables of the outputs are unique, and every observable of an output is one of `observables`. The unit of an output in the tables of sbmlsim 0.8.5 is the one its observable was given (`Formula(..., unit="mM*min")`), else the one of the result.
- **The groups.** A `SensitivityGroup` is a simulation with the changes of the group in its `preinit_changes`; `uid` is its name in the files and folders, `name` in the titles of the figures, `color` the color of the group in the figures of the sampling which compare the groups (the property cycle of matplotlib when not given). Every analysis has its own groups, at least one (the local and global analyses of most models run for the control only, `GROUPS[:1]`, the sampling for every group).
- **The seed** of the designs is `seed`, 1234 by default.

| analysis | design and indices | the tables and figures of sbmlsim 0.8.5, `<analysis>/` | per group, `<analysis>/<group>/` | defaults |
| --- | --- | --- | --- | --- |
| `LocalAnalysis` | `sampling.local` of `delta * k / n_var`, `k = 1, ..., n_var`, with the values of the simulation of the group, `sensitivity.local` | `local_d<delta>_nvar<n_var>_<k>_<group>.tsv` (a row per parameter and output: `parameter`, `output`, `effect`, `Sij`, `|Sij|`, `classification`, `method`, `difference`, `n_var`, sorted by `|Sij|`; `Sij` the mean of `normalized` over the designs, as sbmlsim 0.8.5 averaged its steps) and its heatmap `.png` | none | `delta=0.01`, in (0, 1); `n_var=3` |
| `SamplingAnalysis` | `sampling.lhs`, one design for every group | `sampling_N<n>.tsv` (a row per group and output: `gid`, `gname`, `uid`, `name`, `N`, `mean`, `median`, `std`, `cv`, `min`, `q005`, `q095`, `max`, the uncertainty `U = (P97.5 - P2.5) / P50`, `symbol`, `classification`, `unit`), `sampling_N<n>_compact.tsv` and `sampling_N<n>.tex` (an output per row, `mean (cv in %)` per group), the box, violin and points of every group, `sampling_N<n>_sensitivity_samples.png` (the parameters) and `..._outputs.png` | `samples.tsv` (the outputs and the parameters of every point) | `n=1000` |
| `FastAnalysis` | `sampling.fast`, `sensitivity.fast` | `fast_M4_N<n>_sensitivity_<k>_<group>_S1.png`, `..._ST.png` (heatmaps) and `..._<o>_<output>.png` (the indices of an output with their intervals) | `<output>.S1.tsv`, `<output>.ST.tsv` | `n=1000`, more than 64 |
| `SobolAnalysis` | `sampling.sobol`, `sensitivity.sobol` | `sobol_N<n>_sensitivity_<k>_<group>_S1.png`, `..._ST.png` and `..._<o>_<output>.png` | `<output>.S1.tsv`, `<output>.ST.tsv` | `n=4096`, a power of two |

The number of samples `n` goes through `sample_size`, so `PKPDSIM_SENSITIVITY_N` replaces it in a smoke run (a smoke run of `fast` needs more than 64 samples and one of `sobol` a power of two, so they run with different values). What can go wrong is found before the first analysis runs, not hours into a run:

- when the runner is created: an analysis without a group or with two groups of one id, a `delta` outside of (0, 1), a number of samples which the design refuses (a positive int; for `FastAnalysis` more than 64, for `SobolAnalysis` a power of two); a runner without outputs, analyses or default analyses, with two outputs of one id or one observable, with an output whose observable is none of `observables`, with a default which is none of its analyses or is given twice (a `ValueError`; a string as `default` is a `TypeError`);
- when a run starts, before anything is written: the names (an unknown or repeated name, a string instead of a sequence), the subfolder, and for every named analysis the parameters it varies (none, or one without valid bounds, see `distributions`) and the number of samples which `PKPDSIM_SENSITIVITY_N` gives; the error names the analysis in a note.

The analyses are these four: `Analysis` has no public hook for an analysis of a script's own, a script with its own analysis (simvastatin's central differences) runs it itself.

### The tables and figures of sbmlsim 0.8.5

The sensitivity analyses of the models keep their parameters and outputs across the port onto sbmlsim 0.9: the numbers differ by construction (other designs and samplers), the tables and figures a manuscript uses do not. The files of sbmlsim 0.8.5 are therefore the results of an analysis, with their names, columns and content (see the table above), one result per quantity: the `<k>` of a name is the position of the group in the analysis (`00`), `<o>` the position of the output (`000`). Only what sbmlsim 0.8.5 wrote no file of is written in another form, per group: the samples of the sampling and the tables of the global indices. The heatmaps label the parameters `uid: name` and the outputs by name and leave out the parameters below 0.05 without clustering them (local: `seismic` from -2 to 2; `S1`, `ST`: `viridis` from 0 to 1); a heatmap of which no parameter reaches the cutoff shows every parameter, where `plot_heatmap` would refuse it. An undefined index (an output which is zero at the reference) has no class, a negative uncertainty (a negative median) the class of its magnitude, where the classifications of sbmlsim refuse them. What sbmlsim 0.9 does not give:

- the caches of sbmlsim 0.8.5, `<prefix>_results.pkl` and `<prefix>_sensitivity.pkl`: pickles of its own objects, written with `cache_results=False` too, which no script read;
- the steps of the local analysis: sbmlsim 0.9 has one step per side, so the table averages the designs of `delta / 3`, `2 delta / 3` and `delta` (three scans instead of one), which is the mean over the steps of sbmlsim 0.8.5;
- separate samples per group: the sampling shares one design for every group, so the samples of a parameter are the same in every group;
- the bars of the indices of an output in the form of sbmlsim 0.8.5: they are the figure of `plot_indices`, which draws `S1` and `ST` with their intervals side by side and scales every output, where sbmlsim 0.8.5 overlaid them on one scale for every output of the group.

### Several simulations

A model with several simulations (dulaglutide: a single dose with the pharmacokinetic and the bodyweight outputs, weekly doses with the bodyweight outputs) has a runner per simulation and its own `main`, which reads the simulations with the analyses from the command line (`parse_command_line(argv, choices=...)`, in any order: `local sampling multi`) and runs them with `run_subfolders`: the analyses of every simulation into a subfolder, `sensitivity/<simulation>/<analysis>`, and the one table of the parameters into `sensitivity`, written and printed once:

```python
# not executed
SIMULATIONS = {"single": SINGLE, "multi": MULTI}  # a SensitivityRunner each


def main(argv: Sequence[str] | None = None) -> None:
    words = SINGLE.parse_command_line(argv, choices=list(SIMULATIONS))
    run_subfolders(SIMULATIONS, words.analyses, subfolders=words.choices or ("single",))


if __name__ == "__main__":
    main()
```

`parse_command_line` returns a `CommandLine` of the names of the analyses (the default analyses without one) the words of the choices in their order (none without one) and the `workers` of `--workers` (`None` if absent or `0`); a word which is neither, or is given twice, exits with 2, `--help` lists both. A choice which is the name of an analysis is a `ValueError`. `run_subfolders(runners, names, *, subfolders, path, simulator)` checks everything of every runner before anything is written (the names, the parameters and the sizes, as `run`), and that the runners share one model and give the same parameters, which share one table; `run(names, subfolder=...)` runs one runner into a subfolder.

The outputs and the analyses do not need a model:

```python
from sbmlsim.simulation import Simulation

from pkpdsim.sensitivity import LocalAnalysis, SensitivityGroup, SensitivityOutput

print(SensitivityOutput("cmax", "Cmax").observable)
print(SensitivityOutput("[Cve_api]_cmax", "API Cmax", observable="Cve_api_cmax"))
control = SensitivityGroup("control", "Control", Simulation(end=24, steps=100))
print(LocalAnalysis(groups=[control]).delta)
```
