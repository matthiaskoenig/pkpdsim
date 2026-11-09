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
"""

from collections.abc import Callable, Sequence

import numpy as np
from sbmlsim import Q
from sbmlsim.experiment.draws import UncertaintyRun
from sbmlsim.result.scan import ScanResult, point_timecourse

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
        ValueError: if a selection is no time course of the result.
    """
    values: list[float] = []
    for selection in selections:
        points, series = point_timecourse(res, selection)
        at = Q(points, res.units["time"]).m_as(time_unit)
        distance = np.abs(at - time)
        # the last of the nearest points: a change gives two points at its time
        k = int(np.flatnonzero(distance == distance.min())[-1])
        values.append(Q(series[k], res.units[selection]).m_as(unit))
    return np.asarray(values, dtype=float)


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
    """
    draws = _draws(res, run)
    if not draws or run is None:
        return Band(zone_values(res, selections, time, unit, time_unit=time_unit))
    return _band(
        np.stack(
            [
                zone_values(
                    res.isel(**{run.id: k}),
                    selections,
                    time,
                    unit,
                    time_unit=time_unit,
                )
                for k in range(draws)
            ]
        ),
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
        """The time and the values of one point of the scan, in their units."""
        points, values = point_timecourse(result, selection)
        return (
            Q(points, res.units["time"]).m_as(time_unit),
            Q(values, res.units[selection]).m_as(unit),
        )

    draws = _draws(res, run)
    if not draws or run is None:
        time, values = series(res)
        return time, Band(values)
    time, _ = series(res.isel(**{run.id: 0}))
    return time, _band(
        np.stack([series(res.isel(**{run.id: k}))[1] for k in range(draws)]), run
    )


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
