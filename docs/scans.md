# Scans and sensitivity analyses

The scans and the sensitivity analyses of a model are sbmlsim's: a design of `sbmlsim.simulation.sampling` is a dimension of a `Scan`, its result a `ScanResult`, and the runner of the sensitivity analyses of a model with the table of their parameters is `sbmlsim.sensitivity.runner` and `sbmlsim.sensitivity.parameters`. `pkpdsim.sensitivity` holds what is pharmacokinetic or of the PBPK library: the observables of the non-compartmental analysis of a substance, the parameters of an analysis of a model of the library and its local analysis of the nonzero parameters.

## The points of a scan

What the figures and the pharmacokinetic parameters of a model read from the points of a scan are the functions of [`sbmlsim.result.scan`](https://matthiaskoenig.github.io/sbmlsim/api/result.scan/) (`point_timecourse`, `single_dimension`, `state_after_changes`, `start_values`, `end_values`, `scan_dimensions`, `point_band`), described with examples in [Scans](https://matthiaskoenig.github.io/sbmlsim/scans/) of sbmlsim.

## The sensitivity analyses of a model

The sensitivity script of a model (`<model>/sensitivity/sensitivity_analysis.py`) runs its analyses with the runner of sbmlsim, see [Sensitivity analyses of a project](https://matthiaskoenig.github.io/sbmlsim/sensitivity/#sensitivity-analyses-of-a-project) and the API of [`sbmlsim.sensitivity.runner`](https://matthiaskoenig.github.io/sbmlsim/api/sensitivity.runner/) and [`sbmlsim.sensitivity.parameters`](https://matthiaskoenig.github.io/sbmlsim/api/sensitivity.parameters/): `SensitivityRunner`, `SensitivityGroup`, `SensitivityOutput`, the analyses `LocalAnalysis`, `SamplingAnalysis`, `FastAnalysis` and `SobolAnalysis`, `run_subfolders`, `sample_size` and the table of the parameters, `parameters_from_sbml`, `set_bounds`, `set_relative_bounds`, `write_parameters`, `distributions` and `free_parameters`. The runner writes the tables and figures of sbmlsim 0.8.5, which the manuscripts read. A `PKPDModel` is the project of a runner:

- **The project.** `SensitivityRunner(project=MODEL, ..., out=MODEL.results_dir / "sensitivity")`: the runner simulates with the tolerances of the model (`MODEL.simulator`) and writes into `out`, the folder of the sensitivity analyses of the model package, unless the command line gives `--out`. A model may list its runners in `sensitivities` of its `PKPDModel` (a field of the project, a runner, a mapping of runners by subfolder or a reference `"module:attribute"`), which runs them with the kind `sensitivity` of [`python -m sbmlsim.run`](https://matthiaskoenig.github.io/sbmlsim/run/).
- **The builds per sex** are the variants of the project (`male`, the default, and `female`, see [Builds per sex](models.md#builds-per-sex)): a group of women runs the female build with `SensitivityGroup(..., variant=Sex.FEMALE)`. A global design is drawn once and carried to the female build relative to the values of the table, which are those of the male build; the table gets the column `value_female`.
- **The equilibration.** A group of a model of the PBPK library whose condition moves a slow state (the fluid and the blood pressure of a PD component in renal impairment) runs from the drug-free steady state of its condition: `SensitivityGroup(..., prepare=Equilibration(resets={"FLUID__Vurine": 0.0}))` with `pkpdsim.pbpk.Equilibration`, the equilibration of the experiments (see [Equilibration of a condition](pbpk/impairment.md#equilibration-of-a-condition)). The simulation of the group starts at the dose (time 0), so an output which compares with the start (a change of the blood pressure against its value at the dose) reads time 0.
- **The command line.** `SENSITIVITY.main()` reads the analyses, `--workers N` (a pool of a PBPK model needs about 300 to 340 MB per process) and `--out DIR`; `SBMLSIM_SENSITIVITY_N` replaces the numbers of samples for a smoke run.

```python
# not executed
from pkpdsim.sensitivity import pk_observables
from sbmlsim import Q
from sbmlsim.sensitivity.runner import (
    LocalAnalysis,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityRunner,
)
from sbmlsim.simulation import Simulation

from pkdb_models.models.hctz import MODEL

OUTPUTS = [
    SensitivityOutput("cmax", "HCTZ Cmax", observable="pk_hctz.cmax"),
    SensitivityOutput("hctz_cl_f", "HCTZ CL/F"),
]


def observables() -> list:
    return pk_observables(
        "hctz",
        "[Cve_hctz]",
        {"cmax": "mM", "cl_f": "ml/min"},
        dose="PODOSE_hctz",
        route="oral",
        per_dose_mass={"cl_f"},
        mr="Mr_hctz",
    )


def sensitivity_simulation(changes: dict[str, float]) -> Simulation:
    return Simulation(
        end=5 * 24 * 60,
        time_unit="min",
        preinit_changes={"PODOSE_hctz": Q(25, "mg"), **changes},
    )


GROUPS = [
    SensitivityGroup("control", "Control", sensitivity_simulation({})),
    SensitivityGroup(
        "sevRI",
        "Severe renal impairment",
        sensitivity_simulation({"f_renal_function": 0.1}),
    ),
]

SENSITIVITY = SensitivityRunner(
    project=MODEL,
    parameters=sensitivity_parameters,  # the table, see below
    outputs=OUTPUTS,
    observables=observables,
    analyses={
        "local": LocalAnalysis(groups=GROUPS[:1]),
        "sampling": SamplingAnalysis(groups=GROUPS),
    },
    out=MODEL.results_dir / "sensitivity",
)

if __name__ == "__main__":
    SENSITIVITY.main()
```

```bash
python -m pkdb_models.models.hctz.sensitivity.sensitivity_analysis             # the default analyses, local
python -m pkdb_models.models.hctz.sensitivity.sensitivity_analysis local sampling --workers 8
SBMLSIM_SENSITIVITY_N=65 python -m pkdb_models.models.hctz.sensitivity.sensitivity_analysis fast
```

A model with several simulations (dulaglutide) has a runner per simulation and runs them with `run_subfolders(runners, words.analyses, subfolders=..., workers=words.workers)` of sbmlsim, the words read with `parse_command_line(argv, choices=...)` of one of its runners.

## The observables of the pharmacokinetics

`pk_observables(prefix, selection, outputs, *, dose, route, options, per_dose_mass, mr, pk_id)` gives the non-compartmental analysis of a concentration, the `PK` observable `pk_<prefix>` (or `pk_id`, since an id of the model, e.g. the species `hctz`, is no id of an observable) with the parameters of the outputs, and a `Formula` `<prefix>_<key>` per output in its unit. The dose of an analysis is a mass (`PODOSE_<s>` in mg), so an output whose unit has the amount of the dose (`per_dose_mass`, e.g. `cl_f` and `vz_f`) is divided by the molecular weight of the dosed substance at the start, `at(<mr>, 0)`, which it then requires. The observables are the `observables` of a runner, a list or a function which builds them when a run starts; an output names the formula `<prefix>_<key>` or the parameter `<pk_id>.<key>` of the analysis. `max_after_dose(time, values)` is the function of a `Custom` observable of one selection, its maximum from the dose at time 0 on (a simulation which starts before the dose).

```python
from pkpdutils import NCAOptions
from sbmlsim.simulation import Custom

from pkpdsim.sensitivity import max_after_dose, pk_observables

observables = [
    *pk_observables(
        "hctz",
        "[Cve_hctz]",
        {"auc_inf_obs": "mM*min", "cl_f": "ml/min"},
        dose="PODOSE_hctz",
        route="oral",
        options=NCAOptions(lloq=4.5e-10),
        per_dose_mass={"cl_f"},
        mr="Mr_hctz",
        pk_id="pk_hctz",
    ),
    Custom("max_diuresis", max_after_dose, "l/min", symbols=["diuresis"]),
]
for observable in observables:
    print(observable)
```

## The parameters of a model of the PBPK library

A model of the PBPK library has hundreds of constants, most of them initial assignments of its covariates, so its analysis names the parameters it varies: `library_parameters(model_path, ids, fit_ids=..., scaling_ids=..., references={}, relative=0.15)` is the table of the `ids` in their order (`SensitivityParameter`s of `sbmlsim.sensitivity.parameters`). The fitted parameters vary by `value * (1 +- relative)`, the scaling factors (1 in the healthy individual) by `1 +- relative`, the severities `SEVERITIES` (`CP_severity` 0 to 3, `HF_severity` 0 to 4) from 0 to `relative` of their top, and any other parameter by `value * (1 +- relative)`. The reference of a parameter is the `Sources:` list of its notes (`notes_sources`), replaced by the text of `references`. An id which is no parameter of the model, and an id of `fit_ids`, `scaling_ids` or `references` which is not in `ids`, is a `ValueError` naming it. `parameters_from_sbml(..., exclude_derived=True)` of sbmlsim is the table of every constant without the initial assignments. `NonzeroLocalAnalysis` is the `LocalAnalysis` of sbmlsim which leaves out the parameters whose value is 0 (a severity of the healthy individual), whose local sensitivity is 0/0; the sampled analyses vary the severities within their bounds.

```python
# not executed
from pkpdsim.sensitivity import NonzeroLocalAnalysis, library_parameters


def sensitivity_parameters() -> list:
    return library_parameters(
        MODEL_PATH,
        ["BW", "f_renal_function", "CP_severity", "fu_ref_api"],
        fit_ids={"fu_ref_api"},
        scaling_ids={"f_renal_function"},
        references={"f_renal_function": r"\cite{Stevens2024}"},
    )


LOCAL = NonzeroLocalAnalysis(groups=GROUPS[:1])
```
