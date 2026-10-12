# Scans and sensitivity analyses

Two modules hold what the scripts of the models around scans and sensitivity analyses share. Both build on the API of sbmlsim 0.9, where a design of `sbmlsim.simulation.sampling` is a dimension of a `Scan` and its result is a `ScanResult`.

## The timecourse of a point of a scan

sbmlsim 0.9 returns a scan as a `ScanResult` with the layout `(*dims, time)`, and at the time of a change after the start (a second dose) it keeps the state before and after the change. The pharmacokinetic parameters of a scan (`<model>_pk.py`) and the figures of a model read the timecourse of one point as sbmlsim 0.8 gave it, with the state after the change only, so that their numbers and figures stay the same. `pkpdsim.scan.single_dimension(res)` gives the id of the dimension of a scan of one dimension and `point_timecourse(res, key, **point)` the time points and the values of a variable at the point `point` (the position along every dimension, e.g. `dim_dose=2`); a result without dimension needs none:

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
design = sobol(
    distributions(parameters), n=sample_size(4096), seed=1234, model=MODEL_PATH
)
```

- `parameters_from_sbml` takes the constant parameters, compartments and species of a model in its order, with their values after the initialization, and leaves out those without a value, those close to zero and `exclude_ids`.
- `set_bounds` sets the bounds and the type of parameters by id (an unknown id is logged and skipped), `set_relative_bounds` lets every parameter without bounds vary by `value * (1 +- fraction)`.
- `write_parameters` writes `parameters.tsv` and the LaTeX table `parameters.tex` (`jinja2` is needed for the table), `parameters_to_df` gives the table as a `DataFrame` with the columns `uid`, `name`, `value`, `lower_bound`, `upper_bound`, `unit`, `type` and `reference`, sorted by type and id.
- `distributions` gives the uniform distribution of every parameter in its bounds, the input of a global design; a parameter without bounds raises a `ValueError` naming it.
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
