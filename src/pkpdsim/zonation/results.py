"""The bridge from a simulation result to the arrays the zone figures draw.

`pkpdsim.zonation.figures` draws plain numpy arrays and knows no simulator;
this module reads them out of a `sbmlsim.result.scan.ScanResult`, with the
units of the result and with the draws of an adopted fit as the band of a
figure: `zone_values` the zones at a time, `zone_band` the same over the
draws, `series_band` a time course over the draws, `zone_column` the
`ZoneColumn` of a condition of `zone_profiles` and `band_text` the text a
title names the band with. `Band` is a drawing input of plain arrays and lives
with the figures, next to `ZoneColumn`; it is re-exported here, so a caller
imports it from either module.

A time is matched and not interpolated: of the points of the result the one
nearest the time is taken, and of several points at the same distance the
last, because a change of the simulation (a medium change, a second dose)
gives two points at its time and the state after the change is the one to
draw. Every time is read and returned in `time_unit`, minutes by default, the
convention of the engine.

A band reads every selection once over all its draws (`_rows`, the draws of
the result as the first axis) and converts the time of the result once per
call, not once per draw and selection: a figure of three conditions, five
rows, three times, ten zones and a hundred draws is a few hundred reads of
the dataset instead of tens of thousands. Taking the last of the points
nearest the time over the raw rows of the dataset is what
`sbmlsim.result.scan.point_timecourse` of a single draw gives: it drops the
state before a change, which is the earlier of the two rows at its time, and
the later row is the one taken here either way;
`tests/zonation/test_results.py` pins the two against each other, on a grid
of times and on a ragged result whose draws keep their own time points.
"""

from collections.abc import Callable, Sequence

import numpy as np
from sbmlsim import Q
from sbmlsim.experiment.draws import UncertaintyRun
from sbmlsim.result.scan import POINT, TIME, ScanResult, point_timecourse

from pkpdsim.zonation.figures import (
    DEFAULT_TIME_UNIT,
    Band,
    ZoneColumn,
    ZoneProfile,
)

__all__ = [
    "Band",
    "band_text",
    "series_band",
    "zone_band",
    "zone_column",
    "zone_values",
]


def _rows(res: ScanResult, key: str, across: str | None) -> np.ndarray:
    """Get a variable of the result as `(rows, points)`, a row per draw.

    Args:
        res: the result.
        key: the variable or the coordinate.
        across: the dimension of the draws, `None` for a single simulation.

    Returns:
        The values, one row per draw of `across` and one row for a variable
        which does not carry the dimension (the time of a grid of times, which
        every draw shares).

    Raises:
        KeyError: naming the variables of the result, if it has no variable of
            the key (`ScanResult.__getitem__`).
        ValueError: if the result has a dimension of a scan besides `across`,
            which leaves no one value per time, or if the key is no time
            course of the result.
    """
    extra = [d for d in res.dims if d != across]
    if extra:
        raise ValueError(
            f"The result has the dimensions {extra} besides the draws: read "
            f"one point of the scan (`res.isel(...)`), it has no one value "
            f"per time."
        )
    array = res[key]
    tdim = POINT if res.ragged else TIME
    if tdim not in array.dims:
        raise ValueError(f"'{key}' is no timecourse of the result.")
    if across is not None and across in array.dims:
        return np.asarray(array.transpose(across, tdim).values, dtype=float)
    return np.asarray(array.transpose(tdim).values, dtype=float)[np.newaxis, :]


def _nearest(times: np.ndarray, time: float) -> np.ndarray:
    """Get the index of the last point nearest a time, per row of the times.

    The padding of a ragged result is `NaN` and is no point: a row of padding
    alone (a draw whose simulation failed) gives the index 0, whose value is
    the `NaN` the median and the quantiles of the draws leave out.

    Args:
        times: the time points `(rows, points)`, padded with `NaN`.
        time: the time, in the unit of `times`.

    Returns:
        The index of the point per row.
    """
    distance = np.abs(times - time)
    with np.errstate(invalid="ignore"):
        nearest = np.nanmin(distance, axis=1, keepdims=True)
    at = distance == nearest
    # the last of the nearest points: a change gives two points at its time,
    # and the second is the state after it
    return at.shape[1] - 1 - np.argmax(at[:, ::-1], axis=1)


def _values_at(
    res: ScanResult,
    selections: Sequence[str],
    time: float,
    unit: str,
    *,
    across: str | None,
    time_unit: str,
) -> np.ndarray:
    """Get the values of the selections at a time, a row per draw.

    Every selection is read once over all its draws and the time of the
    result is converted once, see the module.

    Args:
        res: the result of one point of a scan.
        selections: the selections.
        time: the time, in `time_unit`.
        unit: the unit of the values.
        across: the dimension of the draws, `None` for a single simulation.
        time_unit: the unit of `time`.

    Returns:
        The values `(rows, len(selections))`, one row per draw and one row
        without `across`.
    """
    at = Q(_rows(res, TIME, across), res.units[TIME]).m_as(time_unit)
    index = _nearest(at, time)
    rows = int(res.ds.sizes[across]) if across is not None else 1
    if index.size != rows:  # a grid of times, which every draw shares
        index = np.broadcast_to(index, (rows,))
    draw = np.arange(rows)
    values = np.empty((rows, len(selections)), dtype=float)
    for column, selection in enumerate(selections):
        series = _rows(res, selection, across)
        picked = series[draw if series.shape[0] == rows else 0, index]
        values[:, column] = Q(picked, res.units[selection]).m_as(unit)
    return values


def zone_values(
    res: ScanResult,
    selections: Sequence[str],
    time: float,
    unit: str,
    *,
    time_unit: str = DEFAULT_TIME_UNIT,
) -> np.ndarray:
    """Get the values of the zones at a time, in a unit (the state after a change).

    Args:
        res: the result of the simulation, of one point of a scan.
        selections: the selections of the zones, periportal first.
        time: the time, in `time_unit`.
        unit: the unit of the values.
        time_unit: the unit of `time`, the convention of the engine (min).

    Returns:
        The value of every selection, in `unit`.

    Raises:
        KeyError: naming the variables of the result, if it has no variable of
            a selection.
        ValueError: if a selection is no time course of the result, or if the
            result is more than one point of a scan.
    """
    return _values_at(res, selections, time, unit, across=None, time_unit=time_unit)[0]


def _draws(res: ScanResult, run: UncertaintyRun | None) -> int:
    """Get the number of draws of a result, 0 for a result without draws."""
    return int(res.ds.sizes[run.id]) if run is not None and run.id in res.dims else 0


def _band(values: np.ndarray, run: UncertaintyRun) -> Band:
    """Get the median and the quantiles of values by draw (the first axis)."""
    lower, upper = np.nanquantile(values, run.quantiles, axis=0)
    return Band(np.nanmedian(values, axis=0), lower, upper)


def zone_band(
    res: ScanResult,
    selections: Sequence[str],
    time: float,
    unit: str,
    *,
    run: UncertaintyRun | None = None,
    time_unit: str = DEFAULT_TIME_UNIT,
) -> Band:
    """Get the values of the zones at a time, over the draws of a run.

    Args:
        res: the result of the simulation, of one point of a scan.
        selections: the selections of the zones, periportal first.
        time: the time, in `time_unit`.
        unit: the unit of the values.
        run: the run of the draws; `None`, or a result which carries no draws
            of the run, gives the values as the median and no band.
        time_unit: the unit of `time`, the convention of the engine (min).

    Returns:
        The band: the median over the draws and the quantiles of the run, or
        the values alone for a result without draws.

    Raises:
        KeyError: naming the variables of the result, if it has no variable of
            a selection.
        ValueError: if a selection is no time course of the result, or if the
            result is more than one point of a scan.
    """
    draws = _draws(res, run)
    if not draws or run is None:
        return Band(zone_values(res, selections, time, unit, time_unit=time_unit))
    return _band(
        _values_at(res, selections, time, unit, across=run.id, time_unit=time_unit),
        run,
    )


def series_band(
    res: ScanResult,
    selection: str,
    unit: str,
    *,
    run: UncertaintyRun | None = None,
    time_unit: str = DEFAULT_TIME_UNIT,
) -> tuple[np.ndarray, Band]:
    """Get the time and the values of one selection, over the draws of a run.

    Args:
        res: the result of the simulation, every draw on the same times.
        selection: the selection.
        unit: the unit of the values.
        run: the run of the draws; `None`, or a result which carries no draws
            of the run, gives the values as the median and no band.
        time_unit: the unit of the returned time, the convention of the engine
            (min).

    Returns:
        The time in `time_unit` and the band of the values in `unit`; every
        draw is on the same times, so the times of the first draw are the
        times of the band.

    Raises:
        ValueError: if the selection is no time course of the result.
    """

    def series(result: ScanResult) -> tuple[np.ndarray, np.ndarray]:
        """The time and the values of one point of the scan, as they are read."""
        return point_timecourse(result, selection)

    def value(values: np.ndarray) -> np.ndarray:
        """The values of a draw in the unit of the call."""
        return Q(values, res.units[selection]).m_as(unit)

    draws = _draws(res, run)
    if not draws or run is None:
        points, values = series(res)
        return Q(points, res.units[TIME]).m_as(time_unit), Band(value(values))
    # a whole time course per draw: the states before the changes are dropped
    # per draw (`point_timecourse`), so this keeps the loop, and only the time
    # of the first draw is read and converted, as every draw has the same
    points, first = series(res.isel(**{run.id: 0}))
    values = np.stack(
        [first, *(series(res.isel(**{run.id: k}))[1] for k in range(1, draws))]
    )
    return Q(points, res.units[TIME]).m_as(time_unit), _band(value(values), run)


def band_text(run: UncertaintyRun | None) -> str:
    """Get the text of the band of a run for a title, empty without a run.

    Args:
        run: the run of the draws, `None` for a simulation without draws.

    Returns:
        The text, ` (median, 5 to 95 % band of the fit)` for the default
        quantiles of a run.
    """
    if run is None:
        return ""
    lower, upper = run.quantiles
    return f" (median, {100 * lower:g} to {100 * upper:g} % band of the fit)"


def _transformed(
    band: Band, transform: Callable[[np.ndarray], np.ndarray] | None
) -> Band:
    """Transform the median and both ends of a band, nothing without a transform."""
    if transform is None:
        return band
    return Band(
        transform(band.median),
        None if band.lower is None else transform(band.lower),
        None if band.upper is None else transform(band.upper),
    )


def zone_column(
    res: ScanResult,
    profiles: Sequence[ZoneProfile],
    *,
    times: Sequence[float],
    run: UncertaintyRun | None = None,
    time_unit: str = DEFAULT_TIME_UNIT,
    label: str = "",
    color: str | None = None,
) -> ZoneColumn:
    """Get the column of a condition of `zone_profiles` from its result.

    Every profile is read at its own times (`ZoneProfile.times`) or at the
    times of the call, in `ZoneProfile.source_unit or ZoneProfile.unit`, and
    `ZoneProfile.transform` is applied to the median and to each end of the
    band: a living fraction read as dimensionless becomes the percentage of
    dead cells.

    Args:
        res: the result of the condition, of one point of a scan.
        profiles: the profiles of the figure, a row each.
        times: the times of the figure, of a profile which has none of its own.
        run: the run of the draws, `None` for a simulation without draws.
        time_unit: the unit of the times, the convention of the engine (min).
        label: the name of the condition, the title of its panel.
        color: the colour of the condition, `None` the colour of the time.

    Returns:
        The column, the band of every profile and time it is read at.
    """
    values: dict[tuple[str, float], Band] = {}
    for profile in profiles:
        unit = profile.source_unit or profile.unit
        at = times if profile.times is None else profile.times
        for time in at:
            band = zone_band(
                res,
                profile.selections,
                float(time),
                unit,
                run=run,
                time_unit=time_unit,
            )
            values[(profile.label, float(time))] = _transformed(band, profile.transform)
    return ZoneColumn(label=label, values=values, color=color)
