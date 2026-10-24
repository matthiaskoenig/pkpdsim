"""The results of the scans of a model: the dimension, the timecourse of a point, the state after a change.

sbmlsim 0.9 returns a scan as a `ScanResult` with the layout `(*dims, time)`
and, at the time of a change after the start, the states before and after it.
The PK of the scans of a model (`<model>_pk.py`) and its own figures read the
timecourse of one point of a scan as sbmlsim 0.8 gave it, with the state after
a change only, so that their numbers and figures stay the same:

```python
# not executed
dim = single_dimension(res)
for k in range(res.ds.sizes[dim]):
    time, values = point_timecourse(res, "[Cve_hctz]", **{dim: k})
```

`state_after_changes(res)` gives the whole result in that form, for the figures
which read the arrays of a result instead of one point; a ragged result of
several simulations on the time points of a figure (`plot_times`), not on the
union of their steps.

A simulation without `steps` or `times` outputs the steps of the integrator, a
scan of such simulations is ragged: every point keeps its own time points
(`(*dims, _point)`, padded with `NaN`). `start_values(res, key)` and
`end_values(res, key)` give the first and the last value of a variable at every
point of a scan in both layouts, e.g. the dose or the scanned parameter at the
start and the amount in the urine at the end.

A run over the draws of the fitted parameters adds the dimension of the draws to
every result (`pkpdsim.uncertainty`); `scan_dimensions(res, across)` gives the
dimensions of the scan without it, and `point_band(res, key, across, quantiles,
**point)` the quantiles and the median of the draws at a point of the scan.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
from sbmlsim.result.scan import CHANGES, POINT, TIME
from sbmlsim.result.timecourse import plot_times, without_states_before

if TYPE_CHECKING:
    import xarray as xr
    from sbmlsim.result import ScanResult


def single_dimension(res: ScanResult) -> str:
    """Get the dimension of a scan of one dimension.

    Args:
        res: the result of the scan.

    Returns:
        The id of the dimension.

    Raises:
        ValueError: if the scan has no dimension or several.
    """
    if len(res.dims) != 1:
        raise ValueError(
            f"A scan of one dimension is required, not of {list(res.dims)}."
        )
    return res.dims[0]


def scan_dimensions(res: ScanResult, across: str | None = None) -> list[str]:
    """Get the dimensions of a scan without the dimension of the draws.

    Args:
        res: the result of the scan.
        across: the dimension of the draws, `None` for a result without.

    Returns:
        The ids of the dimensions, in the order of the result.
    """
    return [d for d in res.dims if d != across]


def state_after_changes(res: ScanResult) -> ScanResult:
    """Get a result with the state after a change at its time, every time once.

    The output of sbmlsim 0.9 has the states before and after a change at its
    time; the figures and analyses of a model read the state after it, as
    sbmlsim 0.8 gave it. The result is interpolated onto its time points
    without the repeated ones (`ScanResult.interpolate`, which gives the state
    after a change at its time). It is the result of `point_timecourse` for all
    points of a scan at once, as a `ScanResult`; a result without a change is
    unchanged.

    A ragged result of several simulations (the steps of the integrator, a
    simulation without `steps` or `times`) is interpolated onto the time points
    of a figure instead (`sbmlsim.result.timecourse.plot_times`: a few hundred,
    dense where the simulations stepped densely, with the times of the
    changes), not onto the union of the steps of all its simulations, whose
    size grows with the square of the number of simulations: a scan of 11
    points over 100 draws took more than 40 GB.

    Args:
        res: the result of a simulation or a scan.

    Returns:
        The result on its time points, each time once.
    """
    times = res.time_points()
    simulations = res.ds[TIME].size // res.ds.sizes[POINT] if res.ragged else 1
    if simulations <= 1:
        return res.interpolate(np.unique(times))
    changes = np.asarray(res.ds.attrs.get(CHANGES, []), dtype=float)
    return res.interpolate(plot_times(times, keep=changes.tolist()))


def _timecourse_variable(res: ScanResult, key: str) -> xr.DataArray:
    """Get a variable which is a timecourse of the result.

    Raises:
        ValueError: if the variable is no timecourse.
    """
    tdim = POINT if res.ragged else TIME
    if key not in res.ds.data_vars or tdim not in res.ds[key].dims:
        raise ValueError(f"'{key}' is no timecourse of the result.")
    return res.ds[key]


def start_values(res: ScanResult, key: str, /) -> xr.DataArray:
    """Get the value of a variable at the start of every point of a scan.

    The first output of every simulation, on a grid (`time`) and in the ragged
    layout of the steps of the integrator (`_point`) alike.

    Args:
        res: the result of the scan, positional only.
        key: the variable, a timecourse, positional only.

    Returns:
        The values over the dimensions of the scan, with their coordinates.

    Raises:
        ValueError: if the variable is no timecourse.
    """
    values = _timecourse_variable(res, key)
    tdim = POINT if res.ragged else TIME
    return values.isel({tdim: 0}).drop_vars(TIME, errors="ignore")


def end_values(res: ScanResult, key: str, /) -> xr.DataArray:
    """Get the value of a variable at the end of every point of a scan.

    The last output of every simulation: the last time of a grid (`time`), the
    last time before the padding of every point of a ragged result (`_point`);
    a point whose simulation failed has `NaN`.

    Args:
        res: the result of the scan, positional only.
        key: the variable, a timecourse, positional only.

    Returns:
        The values over the dimensions of the scan, with their coordinates.

    Raises:
        ValueError: if the variable is no timecourse.
    """
    values = _timecourse_variable(res, key)
    if not res.ragged:
        return values.isel({TIME: -1}).drop_vars(TIME, errors="ignore")
    # the padding is at the end: the last time of a point is its last finite one
    last = (res.ds[TIME].notnull().sum(POINT) - 1).clip(min=0)
    return values.isel({POINT: last})


def point_timecourse(
    res: ScanResult, key: str, /, **point: int
) -> tuple[np.ndarray, np.ndarray]:
    """Get the timecourse of a variable at one point of a scan, as sbmlsim 0.8 gave it.

    At the time of a change after the start the timecourse has the state after
    the change only (`sbmlsim.result.timecourse.without_states_before`), and a
    result whose simulations have different times loses its padding.

    Args:
        res: the result of the scan, positional only.
        key: the variable, a timecourse, positional only, so that a dimension
            of any name, also `res` or `key`, is a position.
        **point: the position of the point along every dimension of the scan,
            e.g. `dim_dose=2`; the result of a simulation needs none.

    Returns:
        The time points and the values.

    Raises:
        ValueError: if a dimension of the scan has no position, a position
            names no dimension of the scan, or the variable is no timecourse.
    """
    missing = sorted(set(res.dims) - set(point))
    unknown = sorted(set(point) - set(res.dims))
    if missing or unknown:
        raise ValueError(
            f"A point needs a position along every dimension of the scan "
            f"{list(res.dims)}: missing {missing}, unknown {unknown}."
        )
    selected = res.isel(**point) if point else res
    tdim = POINT if selected.ragged else TIME
    if key not in selected.ds.data_vars or tdim not in selected.ds[key].dims:
        raise ValueError(f"'{key}' is no timecourse of the result.")
    values = np.asarray(selected.ds[key].values, dtype=float)
    time = np.asarray(selected.ds[TIME].values, dtype=float)
    changes = np.asarray(selected.ds.attrs.get(CHANGES, []), dtype=float)
    kept_time, kept_values = without_states_before(
        time[np.newaxis, :], values[np.newaxis, :], [changes]
    )
    finite = np.isfinite(kept_time[0])
    return kept_time[0][finite], kept_values[0][finite]


def point_band(
    res: ScanResult,
    key: str,
    /,
    across: str,
    quantiles: tuple[float, float] = (0.05, 0.95),
    **point: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Get the band of the draws of a variable at one point of a scan.

    The timecourse of every draw is the one of `point_timecourse`. Draws
    with different time points (the steps of the integrator, a simulation
    without `steps` or `times`) are interpolated linearly onto the time points
    of a figure (`sbmlsim.result.timecourse.plot_times`: a few hundred, dense
    where the draws stepped densely, with the times of the changes), `NaN`
    outside of the times of a draw; a draw whose simulation failed (`NaN`) is
    ignored.

    Args:
        res: the result of the scan, positional only.
        key: the variable, a timecourse, positional only.
        across: the dimension of the draws.
        quantiles: the lower and the upper quantile.
        **point: the position of the point along every other dimension.

    Returns:
        The time points, the lower quantile, the median and the upper quantile.

    Raises:
        ValueError: if the result has no dimension `across` or a dimension
            has no position.
    """
    if across not in res.dims:
        raise ValueError(
            f"The result has no dimension '{across}' of draws: {list(res.dims)}."
        )
    if across in point:
        raise ValueError(
            f"The dimension '{across}' of the draws is reduced, not a point."
        )
    draws = [
        point_timecourse(res, key, **point, **{across: k})
        for k in range(res.ds.sizes[across])
    ]
    if not draws:
        raise ValueError(f"The result has no draws along '{across}'.")
    first = draws[0][0]
    if all(t.shape == first.shape and np.allclose(t, first) for t, _ in draws):
        times = first
        stacked = np.vstack([values for _, values in draws])
    else:
        # the steps of the integrator: every draw on the time points of a
        # figure, not on the union of the steps of all draws
        changes = np.asarray(res.ds.attrs.get(CHANGES, []), dtype=float)
        times = plot_times(np.concatenate([t for t, _ in draws]), keep=changes.tolist())
        stacked = np.vstack(
            [
                np.interp(times, t, values, left=np.nan, right=np.nan)
                if t.size
                else np.full(times.shape, np.nan)
                for t, values in draws
            ]
        )
    low, high = quantiles
    # a time at which every draw failed has no quantile, it is NaN
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lower, median, upper = np.nanquantile(stacked, [low, 0.5, high], axis=0)
    return times, lower, median, upper
