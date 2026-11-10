"""The bridge from a scan result of sbmlsim to the arrays the zone figures draw."""

from collections.abc import Mapping, Sequence

import matplotlib
import numpy as np
import pytest
import xarray as xr

matplotlib.use("Agg")

from sbmlsim import Q
from sbmlsim.experiment.draws import UncertaintyRun, UncertaintySettings
from sbmlsim.result.scan import ScanResult, point_timecourse
from sbmlsim.simulation import Dimension
from sbmlsim.simulation.sampling import UNCERTAINTY_DIMENSION

from pkpdsim.zonation.figures import Band, ZoneProfile
from pkpdsim.zonation.results import (
    band_text,
    series_band,
    zone_band,
    zone_column,
    zone_values,
)

ZONES = 3
CELL = tuple(f"drg_cell_z{k}" for k in range(1, ZONES + 1))
ALIVE = tuple(f"alive_z{k}" for k in range(1, ZONES + 1))


def _scan(
    series: Mapping[str, np.ndarray],
    time: Sequence[float] | np.ndarray,
    units: Mapping[str, str],
    *,
    draws: int = 0,
    time_unit: str = "min",
) -> ScanResult:
    """A result over the time, or over the draws of a run and the time."""
    t = np.asarray(time, dtype=float)
    dims = (UNCERTAINTY_DIMENSION, "time") if draws else ("time",)
    coords: dict[str, np.ndarray] = {"time": t}
    if draws:
        coords[UNCERTAINTY_DIMENSION] = np.arange(draws)
    return ScanResult(
        xr.Dataset(
            {k: (dims, np.asarray(v, dtype=float)) for k, v in series.items()},
            coords=coords,
            attrs={
                "dims": [UNCERTAINTY_DIMENSION] if draws else [],
                "units": {"time": time_unit, **dict(units)},
            },
        )
    )


def _flat(
    levels: Sequence[float] = (1.0, 2.0, 3.0),
    time: Sequence[float] = (0.0, 30.0, 60.0, 120.0),
    unit: str = "mM",
) -> ScanResult:
    """A result of constant zones, the level of the zone `k` in `levels[k - 1]`."""
    t = np.asarray(time, dtype=float)
    return _scan(
        {sid: np.full_like(t, level) for sid, level in zip(CELL, levels, strict=True)},
        t,
        dict.fromkeys(CELL, unit),
    )


def _run(draws: int, quantiles: tuple[float, float] = (0.05, 0.95)) -> UncertaintyRun:
    """The run of a number of draws, the dimension sbmlsim gives them."""
    return UncertaintyRun(
        settings=UncertaintySettings(draws=draws, quantiles=quantiles),
        dimension=Dimension(
            UNCERTAINTY_DIMENSION, values={"p": np.arange(float(draws))}
        ),
    )


def _draws(draws: int = 5) -> ScanResult:
    """A result over draws: the zone `k` is `k + d` in the draw `d`, constant in time."""
    t = np.asarray([0.0, 60.0], dtype=float)
    offsets = np.arange(float(draws))[:, np.newaxis]
    return _scan(
        {
            sid: np.full((draws, t.size), float(k)) + offsets
            for k, sid in enumerate(CELL, start=1)
        },
        t,
        dict.fromkeys(CELL, "mM"),
        draws=draws,
    )


def test_zone_values_converts_the_unit() -> None:
    res = _flat()
    assert zone_values(res, CELL, 30.0, "mM") == pytest.approx([1.0, 2.0, 3.0])
    assert zone_values(res, CELL, 30.0, "uM") == pytest.approx([1e3, 2e3, 3e3])


def test_zone_values_takes_the_state_after_a_change() -> None:
    # two points at t = 60 min, the second after the change: the second is read
    res = _scan(
        {CELL[0]: np.array([0.0, 1.0, 2.0, 5.0, 6.0])},
        [0.0, 30.0, 60.0, 60.0, 120.0],
        {CELL[0]: "mM"},
    )
    assert zone_values(res, CELL[:1], 60.0, "mM") == pytest.approx([5.0])
    # the nearest point, not an interpolation, and the last of the nearest
    assert zone_values(res, CELL[:1], 58.0, "mM") == pytest.approx([5.0])
    assert zone_values(res, CELL[:1], 31.0, "mM") == pytest.approx([1.0])


def test_zone_values_reads_the_time_in_the_time_unit() -> None:
    res = _scan(
        {CELL[0]: np.array([0.0, 1.0, 2.0])},
        [0.0, 60.0, 120.0],
        {CELL[0]: "mM"},
    )
    assert zone_values(res, CELL[:1], 1.0, "mM", time_unit="hr") == pytest.approx([1.0])
    assert zone_values(res, CELL[:1], 1.0, "mM", time_unit="min") == pytest.approx(
        [0.0]
    )


def test_zone_values_converts_the_time_of_the_result() -> None:
    res = _scan(
        {CELL[0]: np.array([0.0, 1.0, 2.0])},
        [0.0, 1.0, 2.0],
        {CELL[0]: "mM"},
        time_unit="hr",
    )
    assert zone_values(res, CELL[:1], 60.0, "mM") == pytest.approx([1.0])
    assert zone_values(res, CELL[:1], 2.0, "mM", time_unit="hr") == pytest.approx([2.0])


def test_zone_band_without_a_run_has_no_band() -> None:
    band = zone_band(_flat(), CELL, 30.0, "mM")
    assert isinstance(band, Band)
    assert band.median == pytest.approx([1.0, 2.0, 3.0])
    assert band.lower is None
    assert band.upper is None


def test_zone_band_of_a_result_without_the_draws_of_the_run_has_no_band() -> None:
    band = zone_band(_flat(), CELL, 30.0, "mM", run=_run(5))
    assert band.median == pytest.approx([1.0, 2.0, 3.0])
    assert (band.lower, band.upper) == (None, None)


def test_zone_band_is_the_median_and_the_quantiles_of_the_draws() -> None:
    band = zone_band(_draws(5), CELL, 60.0, "mM", run=_run(5))
    # the zone k over the draws is k, k + 1, ... k + 4
    assert band.median == pytest.approx([3.0, 4.0, 5.0])
    assert band.lower is not None
    assert band.upper is not None
    assert band.lower == pytest.approx([1.2, 2.2, 3.2])
    assert band.upper == pytest.approx([4.8, 5.8, 6.8])
    assert band.lower == pytest.approx(
        [float(np.nanquantile(np.arange(float(k), k + 5.0), 0.05)) for k in (1, 2, 3)]
    )


def test_series_band_returns_the_times_in_the_unit_of_the_call() -> None:
    res = _scan(
        {CELL[0]: np.array([0.0, 1.0, 2.0])},
        [0.0, 60.0, 120.0],
        {CELL[0]: "mM"},
    )
    time, band = series_band(res, CELL[0], "uM", time_unit="hr")
    assert time == pytest.approx([0.0, 1.0, 2.0])
    assert band.median == pytest.approx([0.0, 1e3, 2e3])
    assert (band.lower, band.upper) == (None, None)
    minutes, _ = series_band(res, CELL[0], "mM")
    assert minutes == pytest.approx([0.0, 60.0, 120.0])


def test_series_band_over_the_draws_is_the_median_and_the_quantiles() -> None:
    time, band = series_band(_draws(5), CELL[1], "mM", run=_run(5))
    assert time == pytest.approx([0.0, 60.0])
    assert band.median == pytest.approx([4.0, 4.0])
    assert band.lower is not None
    assert band.lower == pytest.approx([2.2, 2.2])


def test_band_text_names_the_quantiles() -> None:
    assert band_text(_run(5)) == " (median, 5 to 95 % band of the fit)"
    assert band_text(_run(5, (0.1, 0.9))) == " (median, 10 to 90 % band of the fit)"


def test_band_text_of_no_run_is_empty() -> None:
    assert band_text(None) == ""


def _alive(level: float = 0.8, draws: int = 0) -> ScanResult:
    """A result of the living fraction of the zones, dimensionless."""
    t = np.asarray([0.0, 60.0, 120.0], dtype=float)
    if not draws:
        return _scan(
            {sid: np.full_like(t, level) for sid in ALIVE},
            t,
            dict.fromkeys(ALIVE, "dimensionless"),
        )
    spread = np.linspace(-0.1, 0.1, draws)[:, np.newaxis]
    return _scan(
        {sid: np.full((draws, t.size), level) + spread for sid in ALIVE},
        t,
        dict.fromkeys(ALIVE, "dimensionless"),
        draws=draws,
    )


def test_zone_column_reads_the_source_unit_and_applies_the_transform() -> None:
    # a living fraction of 0.8 read as dimensionless with transform 100 * (1 - v) is 20.0
    profile = ZoneProfile(
        "Dead cells",
        "%",
        ALIVE,
        transform=lambda v: 100.0 * (1.0 - v),
        source_unit="dimensionless",
    )
    column = zone_column(
        _alive(), [profile], times=(60.0,), label="1 ml/min", color="C1"
    )
    assert column.label == "1 ml/min"
    assert column.color == "C1"
    assert column.values[("Dead cells", 60.0)].median == pytest.approx([20.0] * ZONES)


def test_zone_column_transforms_both_ends_of_a_band() -> None:
    profile = ZoneProfile(
        "Dead cells",
        "%",
        ALIVE,
        transform=lambda v: 100.0 * (1.0 - v),
        source_unit="dimensionless",
    )
    run = _run(5)
    column = zone_column(_alive(draws=5), [profile], times=(60.0,), run=run)
    band = column.values[("Dead cells", 60.0)]
    assert band.median == pytest.approx([20.0] * ZONES)
    assert band.lower is not None
    assert band.upper is not None
    # the draws are 0.7 to 0.9 alive, so the 5 % quantile is 0.71 and the 95 % 0.89;
    # the transform decreases, so the end of the lower draws is the larger one
    assert band.lower == pytest.approx([29.0] * ZONES)
    assert band.upper == pytest.approx([11.0] * ZONES)


def test_zone_column_covers_the_times_of_every_profile() -> None:
    # a profile with its own times contributes those keys, not the times of the call
    plain = ZoneProfile("Drug", "mM", CELL)
    own = ZoneProfile("Drug late", "mM", CELL, times=(60.0, 120.0))
    column = zone_column(_flat(), [plain, own], times=(0.0, 30.0))
    assert sorted(column.values) == [
        ("Drug", 0.0),
        ("Drug", 30.0),
        ("Drug late", 60.0),
        ("Drug late", 120.0),
    ]
    assert column.label == ""
    assert column.color is None


def test_zone_column_reads_the_unit_of_the_profile_without_a_source_unit() -> None:
    profile = ZoneProfile("Drug", "uM", CELL)
    column = zone_column(_flat(), [profile], times=(30.0,))
    assert column.values[("Drug", 30.0)].median == pytest.approx([1e3, 2e3, 3e3])


def _reference(
    res: ScanResult,
    selections: Sequence[str],
    time: float,
    unit: str,
    run: UncertaintyRun,
    time_unit: str = "min",
) -> Band:
    """The band read one draw at a time, the implementation the fast read replaces."""

    def values(point: ScanResult) -> np.ndarray:
        out = []
        for selection in selections:
            points, series = point_timecourse(point, selection)
            at = Q(points, res.units["time"]).m_as(time_unit)
            distance = np.abs(at - time)
            k = int(np.flatnonzero(distance == distance.min())[-1])
            out.append(Q(series[k], res.units[selection]).m_as(unit))
        return np.asarray(out, dtype=float)

    draws = int(res.ds.sizes[run.id])
    stacked = np.stack([values(res.isel(**{run.id: k})) for k in range(draws)])
    lower, upper = np.nanquantile(stacked, run.quantiles, axis=0)
    return Band(np.nanmedian(stacked, axis=0), lower, upper)


def _changing(draws: int = 7) -> ScanResult:
    """Draws on a grid of times with two points at the time of a change (60 min)."""
    t = np.array([0.0, 30.0, 60.0, 60.0, 90.0, 120.0])
    rng = np.random.default_rng(11)
    series = {}
    for k, sid in enumerate(CELL, start=1):
        # the state after the change jumps, so the two rows at 60 differ
        shape = np.array([0.0, 0.5, 1.0, 4.0, 4.5, 5.0]) * k
        series[sid] = shape[np.newaxis, :] + rng.random((draws, t.size))
    return _scan(series, t, dict.fromkeys(CELL, "mM"), draws=draws)


def _ragged(draws: int = 5) -> ScanResult:
    """Draws which keep their own time points, padded with NaN (the steps of a run)."""
    points = 6
    rng = np.random.default_rng(3)
    time = np.full((draws, points), np.nan)
    series = {sid: np.full((draws, points), np.nan) for sid in CELL}
    for d in range(draws):
        kept = points - (d % 2)  # the last draw of every pair is padded
        time[d, :kept] = np.linspace(0.0, 120.0, kept)
        for k, sid in enumerate(CELL, start=1):
            series[sid][d, :kept] = k * np.linspace(0.0, 1.0, kept) + rng.random(kept)
    return ScanResult(
        xr.Dataset(
            {
                "time": ((UNCERTAINTY_DIMENSION, "_point"), time),
                **{
                    sid: ((UNCERTAINTY_DIMENSION, "_point"), v)
                    for sid, v in series.items()
                },
            },
            coords={UNCERTAINTY_DIMENSION: np.arange(draws)},
            attrs={
                "dims": [UNCERTAINTY_DIMENSION],
                "units": {"time": "min", **dict.fromkeys(CELL, "mM")},
            },
        )
    )


@pytest.mark.parametrize("at", [0.0, 45.0, 60.0, 61.0, 120.0])
def test_the_band_of_the_fast_read_is_the_band_of_the_draws_one_at_a_time(
    at: float,
) -> None:
    run = _run(7)
    res = _changing(7)
    band = zone_band(res, CELL, at, "uM", run=run)
    expected = _reference(res, CELL, at, "uM", run)
    assert band.median == pytest.approx(expected.median)
    assert band.lower == pytest.approx(expected.lower)
    assert band.upper == pytest.approx(expected.upper)


def test_the_fast_read_matches_the_reference_on_a_ragged_result() -> None:
    run = _run(5)
    res = _ragged(5)
    assert res.ragged
    for at in (0.0, 50.0, 120.0):
        band = zone_band(res, CELL, at, "mM", run=run)
        expected = _reference(res, CELL, at, "mM", run)
        assert band.median == pytest.approx(expected.median), at
        assert band.lower == pytest.approx(expected.lower), at
        assert band.upper == pytest.approx(expected.upper), at


def test_the_band_of_a_scan_of_several_points_is_refused() -> None:
    res = _draws(5)
    extra = ScanResult(
        res.ds.expand_dims({"dose": [1.0, 2.0]}).assign_attrs(
            dims=[UNCERTAINTY_DIMENSION, "dose"], units=res.units
        )
    )
    with pytest.raises(ValueError, match="dose"):
        zone_band(extra, CELL, 60.0, "mM", run=_run(5))


def test_a_selection_the_result_does_not_have_is_named() -> None:
    with pytest.raises(KeyError, match="nope_z1"):
        zone_values(_flat(), ["nope_z1"], 30.0, "mM")
    with pytest.raises(KeyError, match="nope_z1"):
        zone_band(_draws(5), ["nope_z1"], 60.0, "mM", run=_run(5))
