# Uncertainty of the fitted parameters

The fitted parameters of a model are known only within their confidence intervals, and so are its predictions. `pkpdsim.uncertainty` makes the uncertainty of the fitted parameters part of the model and draws the predictions of the simulation experiments as bands. A fit becomes part of a model in three steps: it is adopted, the factory takes its values with their uncertainty, and a run draws them.

## Adopting a fit

A fit writes its results into the results of the model (`results/fit/<fit id>`), which are not tracked. `run_fit <package> --adopt <fit id>` copies the record of a finished fit into the folder of the adopted fits of the model, `PKPDModel.fits_dir` (`fits` next to the output directory of the factory, e.g. `models/results/fits`), which is tracked with the model:

```bash
run_fit pkdb_models.models.hctz --subset PK --runs 40 --cores 10 --identifiability
run_fit pkdb_models.models.hctz --adopt PK_20261015_235422__7a62a
```

The record is the part of the results which the model needs (`RECORD_FILES`): the best parameter sets (`parameters.json`), the problem with its settings and bounds (`optimization_result.json`), the profile likelihood (`identifiability.json`, `identifiability.tsv`), the Fisher information (`fisher.json`, `fisher.tsv`), the metrics and the text report, about 300 kB. A fit without a profile likelihood and without a Fisher information cannot be adopted (run it with `--identifiability`); an adopted fit is only replaced with `--force`. The model object lists its adopted fits:

```python
# not executed
MODEL = PKPDModel(
    name="hctz",
    ...,
    fits=("PK_20261015_235422__7a62a", "PD_20261016_000200__28b30"),
)
```

## The fitted values in the factory

A model file takes the value of a fitted parameter from the record of its fit, so the value, its uncertainty and its provenance come from one place and cannot drift apart. `PKPDModel.fit_records()` reads the records (`FitRecords`), and `FitRecords.parameter(pid, unit=...)` gives the keyword arguments of an `sbmlutils.factory.Parameter`: its `value`, its `uncertainties` and its `unit`. The pid is the id of the parameter in the flat model, as the fit names it:

```python
# not executed
FITS = MODEL.fit_records()

Parameter(
    "HCTZEX_f",
    **FITS.parameter("KI__HCTZEX_f", unit=U.dimensionless),
    name="renal clearance of hydrochlorothiazide relative to the GFR",
)
```

A parameter which no adopted fit has, and a parameter which two adopted fits have, raise. `pkpdsim.testing.check_model` checks with `check_fit_parameters` that every parameter of an adopted fit is a parameter of the flat model with the value of the fit and a distribution.

## The uncertainty in the SBML

The uncertainty is a `distrib:uncertainty` of the SBML parameter (the distrib package of SBML), which comp flattening keeps. It has three children and the provenance in its notes (`ParameterUncertainty`):

| child | value |
| --- | --- |
| `uncertSpan` of the type `confidenceInterval` | the confidence interval of the profile likelihood at the level of the fit (0.95); a side which the profile does not close is the bound of the fit |
| `uncertParameter` of the type `standardError` | the standard error of the Fisher information in the unit of the model (delta method), if the fit has one |
| `uncertParameter` of the type `distribution` | the distribution the draws of the parameter come from |

The distribution follows the confidence interval in the space the fit searched the parameter in, the `parameter_scale` of the fit:

| the profile | the distribution |
| --- | --- |
| closed, logarithmic scale | `lognormal(ln(value), sigma, lb, ub)`: sigma is the mean of the closed half widths of the interval in the logarithm divided by `z = sqrt(chi2(alpha, df))`, truncated to the bounds of the fit |
| open on one side | the closed side gives sigma, the bound cuts off the open side (a value at its bound gives a half lognormal) |
| open on both sides (not identifiable) | `uniform(lb, ub)` |
| linear scale | `normal(value, sigma, lb, ub)` with the same rule |
| no profile (the Fisher fallback) | sigma is the standard error of the Fisher information in the space of the fit |

The notes name the fit, the record, the cost and the number of data points, the bounds and the scale, the method, the class of the profile, the confidence intervals of the profile and of the Fisher information and the distribution. The draws are independent per parameter: distrib states one distribution per parameter, the correlations of the Fisher information stay in the record.

## Running over the draws

`run_simulation <package> --uncertainty` runs the experiments over draws of the fitted parameters of the model into the folder `simulation_uncertainty` of the results; without the option the experiments run once into `simulation`, as before:

```bash
run_simulation pkdb_models.models.hctz --action simulate --experiments all                 # fast
run_simulation pkdb_models.models.hctz --action simulate --experiments all --uncertainty   # with bands
run_simulation pkdb_models.models.hctz --action simulate --experiments all -u --draws 200 --seed 1 --quantiles 0.025 0.975
```

The run reads the distribution of every parameter of `model_path` (`parameter_distributions`) and draws them once in a Latin hypercube into the dimension `uncertainty` (`uncertainty_dimension`, `UncertaintySettings`: 100 draws, seed 1234, the quantiles 0.05 and 0.95). Every experiment gets the draws as its setting `uncertainty` (`UncertaintyRun`), and `PKPDSimulationExperiment.initialize` applies them (`apply_uncertainty`):

- every simulation becomes a scan over the draws, a scan gets the dimension of the draws appended; a fitted parameter which a simulation, or a dimension of its scan, sets itself keeps its value and is left out of the draws of that simulation, and a simulation which sets all of them runs once;
- every curve of a figure which reads a task over the draws becomes the band of the quantiles of the draws with the median, in the colour and with the name of the curve; the curves of the data stay.

An experiment therefore needs no code of its own for the bands. A figure which an experiment draws itself with matplotlib (`figures_mpl`) reads `self.uncertainty` (`None` without draws) and draws the bands itself: `pkpdsim.scan.scan_dimensions(res, across)` gives the dimensions of a scan without the draws and `point_band(res, key, across, quantiles, **point)` the time points, the lower quantile, the median and the upper quantile of a variable at a point:

```python
# not executed
from pkpdsim.scan import point_band, point_timecourse

run = self.uncertainty
if run is None:
    time, values = point_timecourse(res, "[Cve_hctz]", dim_scan=k)
    ax.plot(time, values, color=color)
else:
    time, lower, median, upper = point_band(
        res, "[Cve_hctz]", across=run.id, quantiles=run.quantiles, dim_scan=k
    )
    ax.fill_between(time, lower, upper, color=color, alpha=0.2)
    ax.plot(time, median, color=color)
```

An extra action of `run_simulation` gets the settings of `--uncertainty` as `options.uncertainty`, and `PKPDModel.run_experiments(..., uncertainty=UncertaintySettings())` runs experiments over the draws from a script.
