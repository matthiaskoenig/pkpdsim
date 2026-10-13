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
which read the arrays of a result instead of one point.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from sbmlsim.result.scan import CHANGES, POINT, TIME
from sbmlsim.result.timecourse import without_states_before

if TYPE_CHECKING:
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


def state_after_changes(res: ScanResult) -> ScanResult:
    """Get a result with the state after a change at its time, every time once.

    The output of sbmlsim 0.9 has the states before and after a change at its
    time; the figures and analyses of a model read the state after it, as
    sbmlsim 0.8 gave it. The result is interpolated onto its own time points
    without the repeated ones (`ScanResult.interpolate`, which gives the state
    after a change at its time). It is the result of `point_timecourse` for all
    points of a scan at once, as a `ScanResult`; a result without a change is
    unchanged.

    Args:
        res: the result of a simulation or a scan.

    Returns:
        The result on the grid of its time points, each time once.
    """
    return res.interpolate(np.unique(res.time_points()))


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
