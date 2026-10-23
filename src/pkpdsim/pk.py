"""The pharmacokinetic parameters of the scans of a model.

`scan_nca` is the non-compartmental analysis of pkpdutils of a timecourse at
every point of a scan at once (`pkpdutils.nca` of a batch), instead of one
analysis per point (`nca_single`):

```python
# not executed
dose = Q(start_values(res, "PODOSE_hctz").values, res.units["PODOSE_hctz"]) / Mr
result = scan_nca(
    res, "[Cve_hctz]", dose=dose, route=Route.ORAL, substance="hctz",
    lloq_fraction=1e-6,
)
df = result.to_dataframe()  # a row per point, a column per dimension
```

The labels and the units of the parameters in the figures are `PK_LABELS` and
`PK_UNITS`; the keys are the parameters of the analysis, the base of the
experiments of a model has them as `pk_labels` and `pk_units`.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

from pkpdutils import NCAOptions, NCAResult, Route, Timecourses, nca

if TYPE_CHECKING:
    from sbmlsim.result import ScanResult
    from sbmlsim.units import Quantity

#: the label of a pharmacokinetic parameter
PK_LABELS: Mapping[str, str] = {
    "auc_last": "AUCend",
    "auc_inf_obs": "AUC",
    "cl_f": "Apparent clearance",
    "cl_renal": "Renal clearance",
    "cl_fecal": "Fecal clearance",
    "cmax": "Cmax",
    "thalf": "Half-life",
    "lambda_z": "Elimination rate",
    "vz_f": "Apparent volume",
}

#: the unit of a pharmacokinetic parameter in the figures
PK_UNITS: Mapping[str, str] = {
    "auc_last": "µmole/l*hr",
    "auc_inf_obs": "µmole/l*hr",
    "cl_f": "ml/min",
    "cl_renal": "ml/min",
    "cl_fecal": "ml/min",
    "cmax": "µmole/l",
    "thalf": "hr",
    "lambda_z": "1/hr",
    "vz_f": "l",
}


def scan_nca(
    res: ScanResult,
    key: str,
    *,
    dose: Quantity,
    route: Route,
    substance: str,
    lloq_fraction: float | None = None,
    options: NCAOptions | None = None,
) -> NCAResult:
    """Get the non-compartmental analysis of a timecourse at every point of a scan.

    The timecourses of all points are one batch of pkpdutils
    (`ScanResult.to_timecourses`: the state after a change at its time, the
    padding of a ragged result left out), which `pkpdutils.nca` analyses
    vectorized, with the parameters of `nca_single` of every point.

    Args:
        res: the result of the scan, of any number of dimensions.
        key: the variable, a concentration.
        dose: the dose of every point as an amount of substance, a quantity
            whose magnitude has the shape of the dimensions of the scan, in
            their order (`start_values(res, dose_id)` converted with the
            molecular weight), or one dose of all points.
        route: the route of the dose.
        substance: the name of the substance.
        lloq_fraction: the limit of quantification of every point relative to
            its maximum, `None` without; long simulations decay into the
            numerical noise of the integrator, which would be the terminal
            phase.
        options: the options of the analysis, the defaults of pkpdutils by
            default.

    Returns:
        The analysis over the dimensions of the scan; `to_dataframe()` has a
        row per point with its position along every dimension.
    """
    amount = dose.to("mmole").magnitude
    timecourses = res.to_timecourses(
        key,
        dose={"amount": amount, "unit": "mmol"},
        route=route,
        substance=substance,
    )
    if lloq_fraction is not None:
        values = timecourses.ds["value"]
        sample_dims = set(timecourses.sample_dims)
        cmax = values.max(dim=[d for d in values.dims if d not in sample_dims])
        lloq = (lloq_fraction * cmax).where(cmax > 0)
        timecourses = Timecourses(timecourses.ds.assign(lloq=lloq))
    return nca(timecourses, options=options)
