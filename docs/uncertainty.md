# Uncertainty of the fitted parameters

The fitted parameters of a model are known only within their confidence intervals, and so are its predictions. The records of the adopted fits, the `distrib` uncertainty they give the parameters of the SBML and the run of the experiments over their draws are part of sbmlsim, see [The uncertainty of the fitted parameters](https://matthiaskoenig.github.io/sbmlsim/uncertainty/#the-uncertainty-of-the-fitted-parameters); a model of pkpdsim uses them to make the uncertainty of the fitted parameters part of the model and to draw the predictions of the simulation experiments as bands. This page keeps the commands and the model object of pkpdsim.

## Adopting a fit

A fit writes its results into the results of the model (`results/fit/<fit id>`), which are not tracked. `run_fit <package> --adopt <fit id>` copies the record of a finished fit into the folder of the adopted fits of the model, `PKPDModel.fits_dir` (`fits` next to the output directory of the factory, e.g. `models/results/fits`), which is tracked with the model:

```bash
run_fit pkdb_models.models.hctz --subset PK --runs 40 --cores 10 --identifiability
run_fit pkdb_models.models.hctz --adopt PK_20261015_235422__7a62a
```

The record is the part of the results which the model needs (`sbmlsim.fit.records.RECORD_FILES`, `adopt_fit`, see [Adopting a fit](https://matthiaskoenig.github.io/sbmlsim/uncertainty/#adopting-a-fit)). A fit without a profile likelihood and without a Fisher information cannot be adopted (run it with `--identifiability`); an adopted fit is only replaced with `--force`. The model object lists its adopted fits:

```python
# not executed
MODEL = PKPDModel(
    name="hctz",
    ...,
    fits=("PK_20261015_235422__7a62a", "PD_20261016_000200__28b30"),
)
```

`run_report <package> --adopted` writes the report of every adopted fit from its record, one folder per fit in the reports of the model, see [Commands](commands.md).

## The fitted values in the factory

A model file takes the value of a fitted parameter from the record of its fit, so the value, its uncertainty and its provenance come from one place and cannot drift apart. `PKPDModel.fit_records()` reads the records (`sbmlsim.fit.records.AdoptedRecords`), and `AdoptedRecords.parameter(pid, unit=...)` gives the keyword arguments of an `sbmlutils.factory.Parameter`: its `value`, its `uncertainties` and its `unit`. The pid is the id of the parameter in the flat model, as the fit names it:

```python
# not executed
FITS = MODEL.fit_records()

Parameter(
    "HCTZEX_f",
    **FITS.parameter("KI__HCTZEX_f", unit=U.dimensionless),
    name="renal clearance of hydrochlorothiazide relative to the GFR",
)
```

A parameter which no adopted fit has, and a parameter which two adopted fits have, raise. `sbmlsim.project_checks.check_model` checks with `sbmlsim.experiment.draws.check_fit_parameters` that every parameter of an adopted fit is a parameter of every build of the model (its variants, the male and the female build) with the value of the fit and a distribution, and that the female build has the distribution of the male build: the draws are those of the male build, the default variant, and run on every build.

## The uncertainty in the SBML

The uncertainty is a `distrib:uncertainty` of the SBML parameter with the confidence interval, the standard error and the distribution the draws come from, see [The uncertainty in the SBML](https://matthiaskoenig.github.io/sbmlsim/uncertainty/#the-uncertainty-in-the-sbml) of sbmlsim for the children and the rule by which the distribution follows from the profile. The draws are independent per parameter; the correlations of the Fisher information stay in the record.

## Running over the draws

`run_simulation <package> --uncertainty` runs the experiments over draws of the fitted parameters of the model into the folder `simulation_uncertainty` of the results; without the option the experiments run once into `simulation`, as before:

```bash
run_simulation pkdb_models.models.hctz --action simulate --experiments all                 # fast
run_simulation pkdb_models.models.hctz --action simulate --experiments all --uncertainty   # with bands
run_simulation pkdb_models.models.hctz --action simulate --experiments all -u --draws 200 --seed 1 --quantiles 0.025 0.975
```

The run reads the distribution of every parameter of `model_path`, the default variant, and draws them once in a Latin hypercube (100 draws, seed 1234, the quantiles 0.05 and 0.95 by default, `sbmlsim.experiment.draws.UncertaintySettings`). Every experiment gets the draws as its setting `uncertainty` and the `initialize` of `SimulationExperiment`, which `PKPDSimulationExperiment` derives from, turns its simulations into scans over the draws and its curves into bands, see [Running over the draws](https://matthiaskoenig.github.io/sbmlsim/uncertainty/#running-over-the-draws) of sbmlsim; an experiment therefore needs no code of its own for the bands. A figure which an experiment draws itself with matplotlib (`figures_mpl`) reads `self.uncertainty` (`None` without draws) and draws the bands with `point_band`:

```python
# not executed
from sbmlsim.result.scan import point_band, point_timecourse

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

An extra action of `run_simulation` gets the settings of `--uncertainty` as `options.uncertainty`, and `PKPDModel.run_experiments(..., uncertainty=UncertaintySettings())` (`from sbmlsim.experiment.draws import UncertaintySettings`) runs experiments over the draws from a script, `sbmlsim.project.Project.run_experiments` with the settings and the result directories of the model.
