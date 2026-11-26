"""The zone figures: over a synthetic result, and over the result of a built chain.

A figure is a `sbmlsim.plot.Chart`; the tests read its spec, the plots, their
axes, colour bars and elements of the contract `figure/1`, and render it only
where the drawing itself is the question (no pixels are compared).
"""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import sbmlsim_core
from pbpk.conftest import individual
from sbmlsim.plot.chart import Chart
from sbmlutils.factory import create_model

from pkpdsim.zonation import ids, observers
from pkpdsim.zonation.figures import (
    _PT,
    LOBULE_X,
    LOBULE_Y,
    LOBULE_Y_BARE,
    MISSING,
    POSITION,
    Band,
    ColorScale,
    ZoneColumn,
    ZoneProfile,
    _limits,
    _scale,
    colorbar,
    data_axes,
    draw_lobule,
    fill_band,
    layer_profile,
    lobule_axes,
    lobule_map,
    panel_aspects,
    result_frame,
    space_time,
    zone_profiles,
)
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.observers import add_observers, observer_objects
from pkpdsim.zonation.render import lobule_bands, lobule_isoline

from .test_chain import TRACER, _Diffusing, _settings

N = 10
LITERAL_LAYERS = ["rbc", "plasma", "endothelium", "disse", "cell", "bile"]
LAYERS = (
    Layer.RBC,
    Layer.PLASMA,
    Layer.ENDOTHELIUM,
    Layer.DISSE,
    Layer.CELL,
    Layer.BILE,
)
NEUTRAL = "#4d4d4d"
"""The grey 0.3 of matplotlib, the colour of an entry which stands for a time."""
TRIAD_COLORS = ("#1f77b4", "#d62728", "#2ca02c")
POSITION_LINES = "Position\n(0 periportal, 1 pericentral)"


# --- reading a chart ---------------------------------------------------------


def _plots(chart: Chart) -> list[dict[str, Any]]:
    """The plots of the spec of a chart, through its JSON as the renderer reads it."""
    assert isinstance(chart, Chart)
    return json.loads(chart.to_json())["plots"]


def _elements(plot: dict[str, Any], kind: str) -> list[dict[str, Any]]:
    """The elements of a kind of a plot, in the order of drawing."""
    return [e for e in plot["elements"] if e["kind"] == kind]


def _lines(plot: dict[str, Any]) -> list[dict[str, Any]]:
    """The lines of the simulation of a plot, one per time."""
    return [
        e
        for e in _elements(plot, "curve")
        if str(e.get("legend", "")).startswith("Sim ")
    ]


def _markers(plot: dict[str, Any]) -> list[dict[str, Any]]:
    """The markers of the data of a plot, one curve per time."""
    return [
        e
        for e in _elements(plot, "curve")
        if str(e.get("legend", "")).startswith("Data ")
    ]


def _numbers(values: Sequence[float | None]) -> np.ndarray:
    """The numbers of an array of the spec, a gap (`null`) NaN."""
    return np.array([np.nan if v is None else v for v in values], dtype=float)


def _x(element: dict[str, Any]) -> np.ndarray:
    """The x values of an element."""
    return _numbers(element["x"])


def _y(element: dict[str, Any]) -> np.ndarray:
    """The y values of an element."""
    return _numbers(element["y"])


def _colorbars(plot: dict[str, Any]) -> list[dict[str, Any]]:
    """The colour bars of a plot, one or several."""
    if "colorbar" in plot:
        return [plot["colorbar"]]
    return list(plot.get("colorbars", []))


def _time_colours(count: int) -> list[str]:
    """The colours of the times: the part 0.1 to 0.75 of viridis."""
    return sbmlsim_core.color_scale(
        "viridis", 0.0, 1.0, [float(v) for v in np.linspace(0.1, 0.75, count)]
    )


def _rendered(chart: Chart) -> str:
    """The SVG of a chart, which must be vector."""
    svg = chart.svg()
    assert svg.lstrip().startswith("<svg")
    assert "<image" not in svg
    return svg


# --- the fixtures ------------------------------------------------------------


def _result(layers: tuple[Layer | str, ...] = LAYERS, n: int = N) -> pd.DataFrame:
    """A synthetic result of the substance `drg` with a gradient in space and time."""
    t = np.linspace(0.0, 120.0, 25)
    cols: dict[str, np.ndarray] = {"time": t}
    for i, layer in enumerate(layers):
        for k, sid in enumerate(observers.profile_ids("drg", layer, n), start=1):
            cols[sid] = (i + 1) * k * t / 120.0
    cols[observers.outflow_concentration("drg")] = t / 120.0
    return pd.DataFrame(cols)


def _profiles(count: int, n: int = N) -> list[ZoneProfile]:
    """Profiles `p0`, `p1`... over the plasma ids of `drg`."""
    return [
        ZoneProfile(f"p{i}", "mM", tuple(observers.profile_ids("drg", Layer.PLASMA, n)))
        for i in range(count)
    ]


def _data(**columns: list[object]) -> pd.DataFrame:
    """A long frame of data, by default two points of `p0` and one of `p1` at 60."""
    base: dict[str, list[object]] = {
        "label": ["p0", "p0", "p1"],
        "time": [60.0, 60.0, 60.0],
        "position": [0.25, 0.75, 0.5],
        "value": [1.0, 2.0, 3.0],
        "sd": [0.1, 0.2, 0.3],
    }
    base.update(columns)
    return pd.DataFrame(base)


def _flat(plasma: float = 3.0, cell: float = 0.0) -> pd.DataFrame:
    """A result of constant layers: the plasma at a value, the cells at another."""
    t = np.linspace(0.0, 10.0, 11)
    cols: dict[str, np.ndarray] = {"time": t}
    for layer, value in ((Layer.PLASMA, plasma), (Layer.CELL, cell)):
        for sid in observers.profile_ids("drg", layer, 4):
            cols[sid] = np.full_like(t, value)
    cols[observers.outflow_concentration("drg")] = np.full_like(t, plasma)
    return pd.DataFrame(cols)


def _plasma_cell() -> list[ZoneProfile]:
    """Two profiles over four zones, the plasma and the cells."""
    return [
        ZoneProfile("Plasma", "mM", tuple(observers.profile_ids("drg", "plasma", 4))),
        ZoneProfile("Cell", "mM", tuple(observers.profile_ids("drg", "cell", 4))),
    ]


def _band(values: np.ndarray, spread: float | None = None) -> Band:
    """A band of values, optionally with ends at a distance from the median."""
    median = np.asarray(values, dtype=float)
    if spread is None:
        return Band(median)
    return Band(median, median - spread, median + spread)


def _medians(
    profiles: Sequence[ZoneProfile], times: Sequence[float], level: float
) -> dict[tuple[str, float], np.ndarray]:
    """A value per zone, per profile and per time: `level * (i + 1) * t * k / n`."""
    out: dict[tuple[str, float], np.ndarray] = {}
    for i, profile in enumerate(profiles):
        n = len(profile.selections)
        for t in profile.times if profile.times is not None else times:
            out[(profile.label, float(t))] = (
                level * (i + 1) * float(t) * np.arange(1, n + 1) / n
            )
    return out


def _column(
    label: str,
    *,
    profiles: Sequence[ZoneProfile],
    times: Sequence[float],
    color: str | None = None,
    level: float = 1.0,
    spread: float | None = None,
) -> ZoneColumn:
    """A column of a condition, the values distinct per row, time and level."""
    return ZoneColumn(
        label,
        {
            key: _band(median, spread)
            for key, median in _medians(profiles, times, level).items()
        },
        color,
    )


# --- the profiles over the zones ---------------------------------------------


def test_the_profile_grid_has_a_row_per_profile_with_label_unit_and_a_line_per_time() -> (
    None
):
    chart = zone_profiles(_result(), _profiles(3), times=(0, 60, 120))
    plots = _plots(chart)
    assert len(plots) == 3
    for i, plot in enumerate(plots):
        assert (plot["row"], plot["column"]) == (i, 0)
        assert plot["y"]["label"] == f"p{i} [mM]"
        assert len(_lines(plot)) == 3


def test_the_x_axis_is_the_zone_midpoint() -> None:
    (plot,) = _plots(zone_profiles(_result(), _profiles(1), times=(0,)))
    assert _x(_lines(plot)[0]) == pytest.approx([(k - 0.5) / 10 for k in range(1, 11)])
    assert (plot["x"]["min"], plot["x"]["max"]) == (0.0, 1.0)


def test_the_profile_values_are_the_result_at_the_time() -> None:
    (plot,) = _plots(zone_profiles(_result(), _profiles(1), times=(60,)))
    expected = [(2 * k * 0.5) for k in range(1, N + 1)]  # plasma is the second layer
    assert _y(_lines(plot)[0]) == pytest.approx(expected)


def test_the_data_is_overlaid_where_given() -> None:
    chart = zone_profiles(_result(), _profiles(2), times=(0, 60), data=_data())
    first, second = _plots(chart)
    (markers,) = _markers(first)
    assert list(_x(markers)) == [0.25, 0.75]
    assert list(_y(markers)) == [1.0, 2.0]
    # the error bars are the sd, symmetric about the value, per row
    assert _numbers(markers["yerr"]) == pytest.approx([0.1, 0.2])
    assert markers["style"]["line"]["type"] == "none"
    assert markers["style"]["marker"]["edge"] == "#000000"
    (markers,) = _markers(second)
    assert (list(_x(markers)), list(_y(markers))) == ([0.5], [3.0])
    assert _numbers(markers["yerr"]) == pytest.approx([0.3])
    # the data are drawn over the lines
    assert first["elements"].index(markers_of(first)) > first["elements"].index(
        _lines(first)[-1]
    )


def markers_of(plot: dict[str, Any]) -> dict[str, Any]:
    """The one curve of the markers of the data of a plot."""
    (markers,) = _markers(plot)
    return markers


def test_the_data_of_another_time_is_skipped() -> None:
    data = _data(time=[60.0, 30.0, 60.0])
    first, second = _plots(
        zone_profiles(_result(), _profiles(2), times=(0, 60), data=data)
    )
    assert list(_x(markers_of(first))) == [0.25]
    assert list(_x(markers_of(second))) == [0.5]


def test_data_with_a_label_no_profile_has_is_refused() -> None:
    with pytest.raises(ValueError) as info:
        zone_profiles(
            _result(), _profiles(2), times=(60,), data=_data(label=["p0", "p9", "p1"])
        )
    assert str(info.value) == (
        "The data has the labels ['p9'] which no profile has; "
        "the labels of the profiles are ['p0', 'p1']."
    )


def test_data_of_which_no_row_is_drawn_is_refused() -> None:
    with pytest.raises(ValueError) as info:
        zone_profiles(_result(), _profiles(2), times=(0, 120), data=_data())
    assert str(info.value) == (
        "No row of the data is drawn: its times [60.0] include none of the "
        "times [0, 120] of the figure."
    )


def test_nothing_to_draw_is_refused() -> None:
    result = _result()
    with pytest.raises(ValueError, match="profiles are empty"):
        zone_profiles(result, [], times=(0,))
    with pytest.raises(ValueError, match="times are empty"):
        zone_profiles(result, _profiles(1), times=())
    with pytest.raises(ValueError, match="columns are empty"):
        zone_profiles([], _profiles(1), times=(0,))
    with pytest.raises(ValueError, match="layers are empty"):
        space_time(result, "drg", layers=(), outflow=False)


def test_the_unit_of_a_label_is_an_argument() -> None:
    heat, out = _plots(space_time(_result(), "drg", layers=(Layer.PLASMA,), unit="uM"))
    assert out["y"]["label"] == "drg outflow [uM]"
    assert _colorbars(heat)[0]["label"] == "drg plasma [uM]"
    # the unit as the caller wrote it, not the short form of pint
    assert "unit" not in out["y"]
    plots = _plots(layer_profile(_result(), "drg", time=0, unit="uM"))
    assert plots[0]["y"]["label"] == "drg [uM]"
    assert _plots(layer_profile(_result(), "drg", time=0))[0]["y"]["label"] == (
        "drg [mM]"
    )


def test_a_missing_selection_names_the_column() -> None:
    bad = [ZoneProfile("x", "mM", ("nope_z1",))]
    with pytest.raises(KeyError, match="nope_z1"):
        zone_profiles(_result(), bad, times=(0,))


def test_every_position_label_names_the_ends_of_the_sinusoid() -> None:
    end = "(0 periportal, 1 pericentral)"
    plots = _plots(zone_profiles(_result(), _profiles(2), times=(0,)))
    assert plots[-1]["x"]["label"] == f"Position {end}" == POSITION
    heat = _plots(space_time(_result(), "drg", layers=(Layer.PLASMA,)))[0]
    assert heat["y"]["label"].replace("\n", " ") == f"Position {end}"
    bottom = _plots(layer_profile(_result(), "drg", time=0))[-1]
    assert bottom["x"]["label"].replace("\n", " ") == f"Position {end}"


def test_the_time_carries_its_unit_in_every_figure() -> None:
    result = _result()
    (plot,) = _plots(zone_profiles(result, _profiles(1), times=(0, 60)))
    assert [e["label"] for e in _lines(plot)] == ["Sim t = 0 min", "Sim t = 60 min"]
    (plot,) = _plots(zone_profiles(result, _profiles(1), times=(60,), time_unit="h"))
    assert _lines(plot)[0]["label"] == "Sim t = 60 h"
    heat = space_time(result, "drg", layers=(Layer.PLASMA,))
    assert _plots(heat)[-1]["x"]["label"] == "Time [min]"
    heat = space_time(result, "drg", layers=(Layer.PLASMA,), time_unit="h")
    assert _plots(heat)[-1]["x"]["label"] == "Time [h]"
    assert layer_profile(result, "drg", time=60).spec()["title"] == "drg at t = 60 min"
    assert layer_profile(result, "drg", time=2, time_unit="h").spec()["title"] == (
        "drg at t = 2 h"
    )


def test_a_layer_of_zeros_is_drawn_from_zero_to_one() -> None:
    result = _flat()
    plasma, cell = _plots(zone_profiles(result, _plasma_cell(), times=(5,)))
    assert plasma["y"]["min"] == 0.0
    assert (cell["y"]["min"], cell["y"]["max"]) == (0.0, 1.0)
    heat_plasma, heat_cells, out = _plots(
        space_time(result, "drg", layers=(Layer.PLASMA, Layer.CELL))
    )
    scale = _colorbars(heat_cells)[0]["scale"]
    assert (scale["low"], scale["high"]) == (0.0, 1.0)
    scale = _colorbars(heat_plasma)[0]["scale"]
    assert (scale["low"], scale["high"]) == (0.0, 3.0)
    assert (out["y"]["min"], out["y"]["max"]) == pytest.approx((0.0, 3.0 * 1.05))
    own = {
        p["title"]: p["y"]
        for p in _plots(layer_profile(result, "drg", time=5, share_y=False))
    }
    assert (own["cell"]["min"], own["cell"]["max"]) == (0.0, 1.0)
    assert (own["plasma"]["min"], own["plasma"]["max"]) == pytest.approx(
        (0.0, 3.0 * 1.05)
    )


def test_a_constant_layer_is_drawn_from_zero_to_its_value() -> None:
    """Zero to one would hide a constant above one, so the scale runs to the value."""
    heat, _ = _plots(space_time(_flat(plasma=0.25), "drg", layers=(Layer.PLASMA,)))
    scale = _colorbars(heat)[0]["scale"]
    assert (scale["low"], scale["high"]) == (0.0, 0.25)


def test_a_negative_value_leaves_the_axis_to_the_renderer() -> None:
    result = _flat()
    result[observers.profile_ids("drg", "cell", 4)] = -1.0
    _, cell = _plots(zone_profiles(result, _plasma_cell(), times=(5,)))
    assert "min" not in cell["y"]
    assert "max" not in cell["y"]


def test_every_time_has_a_colour_and_its_data_has_the_colour_of_the_time() -> None:
    data = _data(time=[0.0, 60.0, 60.0], label=["p0", "p0", "p0"])
    chart = zone_profiles(
        _result(), _profiles(1), times=(0, 60, 120), data=data.iloc[:2]
    )
    (plot,) = _plots(chart)
    lines = _lines(plot)
    colours = [line["style"]["color"] for line in lines]
    assert colours == _time_colours(3)
    assert [line["style"]["line"]["type"] for line in lines] == [
        "dotted",
        "dashed",
        "solid",
    ]
    # the markers of 0 and of 60 are drawn in the colour of their curve
    markers = _markers(plot)
    assert [m["style"]["color"] for m in markers] == colours[:2]
    assert [m["style"]["marker"]["type"] for m in markers] == ["square", "diamond"]


def test_every_row_has_one_legend_right_of_the_grid() -> None:
    chart = zone_profiles(_result(), _profiles(2), times=(0, 60))
    assert chart.spec()["legend"] == {"scope": "row", "position": "outside"}
    data = _data(label=["p0", "p0", "p0"])
    first, second = _plots(
        zone_profiles(_result(), _profiles(2), times=(0, 60), data=data)
    )
    # the entries of a row are its lines and the data of its times, in that order
    assert [e["legend"] for e in first["elements"] if "legend" in e] == [
        "Sim t = 0 min",
        "Sim t = 60 min",
        "Data t = 60 min",
    ]
    assert first["legend_order"][:4] == [
        "Sim t = 0 min",
        "Data t = 0 min",
        "Sim t = 60 min",
        "Data t = 60 min",
    ]
    assert [e["legend"] for e in second["elements"] if "legend" in e] == [
        "Sim t = 0 min",
        "Sim t = 60 min",
    ]
    # the colours are the times, so no entry takes another style than its line
    assert "legend_swatches" not in first


def test_a_frame_draws_one_column_as_before() -> None:
    chart = zone_profiles(_result(), _profiles(2), times=(0, 60))
    plots = _plots(chart)
    assert len(plots) == 2
    assert (chart.spec()["grid"]["rows"], chart.spec()["grid"]["columns"]) == (2, 1)
    assert ["title" in p for p in plots] == [False, False]
    assert [p["y"]["label"] for p in plots] == ["p0 [mM]", "p1 [mM]"]
    assert plots[-1]["x"]["label"] == POSITION
    assert all(not _elements(p, "area") for p in plots)


def test_a_column_per_condition_has_its_label_as_the_title_of_the_top_row() -> None:
    profiles = _profiles(2, n=4)
    times = (0.0, 60.0)
    columns = [
        _column("1 ml/min", profiles=profiles, times=times, level=1.0),
        _column("3 ml/min", profiles=profiles, times=times, level=2.0),
    ]
    chart = zone_profiles(columns, profiles, times=times)
    plots = _plots(chart)
    assert len(plots) == 4
    assert (chart.spec()["grid"]["rows"], chart.spec()["grid"]["columns"]) == (2, 2)
    assert [p.get("title") for p in plots] == ["1 ml/min", "3 ml/min", None, None]
    # the y label on the first column of every row, the x label on the bottom row
    assert [p["y"]["label_visible"] for p in plots] == [True, False, True, False]
    assert [p["y"].get("label") for p in plots] == ["p0 [mM]", None, "p1 [mM]", None]
    assert [p["x"]["label_visible"] for p in plots] == [False, False, True, True]
    assert [p["x"].get("label") for p in plots] == [None, None, POSITION, POSITION]
    # the inner panels keep the marks of the ticks and lose their numbers
    assert [p["y"]["ticks_visible"] for p in plots] == [True, False, True, False]
    assert all(p["y"]["tick_marks_visible"] for p in plots)
    # the values of the column are drawn, the level of the second twice the first
    assert _y(_lines(plots[1])[1]) == pytest.approx(2.0 * _y(_lines(plots[0])[1]))


def test_the_lines_of_a_column_take_its_colour_and_the_style_of_their_time() -> None:
    profiles = _profiles(1, n=4)
    times = (0.0, 60.0, 120.0)
    columns = [
        _column("A", profiles=profiles, times=times, color="#ff0000"),
        _column("B", profiles=profiles, times=times, color="blue", level=2.0),
    ]
    first, second = _plots(zone_profiles(columns, profiles, times=times))
    assert [line["style"]["color"] for line in _lines(first)] == ["#ff0000"] * 3
    assert [line["style"]["color"] for line in _lines(second)] == ["#0000ff"] * 3
    # the style is the time, so the same styles in both columns
    styles = ["dotted", "dashed", "solid"]
    assert [line["style"]["line"]["type"] for line in _lines(first)] == styles
    assert [line["style"]["line"]["type"] for line in _lines(second)] == styles


def test_a_column_without_a_colour_takes_the_colour_of_the_time() -> None:
    profiles = _profiles(1, n=4)
    times = (0.0, 60.0, 120.0)
    columns = [_column("A", profiles=profiles, times=times)]
    (plot,) = _plots(zone_profiles(columns, profiles, times=times))
    (frame,) = _plots(zone_profiles(_result(n=4), profiles, times=times))
    drawn = [line["style"]["color"] for line in _lines(plot)]
    assert len(set(drawn)) == 3
    assert drawn == [line["style"]["color"] for line in _lines(frame)]
    assert drawn == _time_colours(3)


def test_a_band_is_filled_between_its_ends() -> None:
    profiles = _profiles(1, n=4)
    columns = [
        _column("A", profiles=profiles, times=(60.0,), color="#ff0000", spread=0.5)
    ]
    (plot,) = _plots(zone_profiles(columns, profiles, times=(60.0,)))
    (fill,) = _elements(plot, "area")
    assert fill["style"]["color"] == "#ff0000"
    assert fill["style"]["fill"]["opacity"] == pytest.approx(0.25)
    assert "label" not in fill  # the line of the median names it
    (line,) = _lines(plot)
    median = _y(line)
    assert _numbers(fill["from"]) == pytest.approx(median - 0.5)
    assert _numbers(fill["to"]) == pytest.approx(median + 0.5)
    # the fill is drawn under its line, which stays on top
    assert plot["elements"].index(fill) < plot["elements"].index(line)
    # the band is in the limits of the row, the median alone would cut it off
    assert plot["y"]["max"] > median.max() + 0.5


def test_a_band_without_ends_fills_nothing() -> None:
    profiles = _profiles(1, n=4)
    columns = [_column("A", profiles=profiles, times=(60.0,), color="#ff0000")]
    (plot,) = _plots(zone_profiles(columns, profiles, times=(60.0,)))
    assert [e["kind"] for e in plot["elements"]] == ["curve"]


def test_a_row_draws_its_own_times() -> None:
    plain = ZoneProfile("Drug", "mM", tuple(observers.profile_ids("drg", "plasma", 4)))
    own = ZoneProfile(
        "Dead",
        "%",
        tuple(observers.profile_ids("drg", "cell", 4)),
        times=(30.0, 90.0),
    )
    profiles = [plain, own]
    times = (0.0, 60.0)
    # the column covers every time of the figure and of a row, so a row which
    # drew the times of the other one would draw a band it has, not raise
    every = (0.0, 30.0, 60.0, 90.0)
    covering = [ZoneProfile(p.label, p.unit, p.selections) for p in profiles]
    columns = [_column("A", profiles=covering, times=every, color="#ff0000")]
    top, bottom = _plots(zone_profiles(columns, profiles, times=times))
    expected = _medians(covering, every, 1.0)
    assert [_y(line) for line in _lines(top)] == [
        pytest.approx(expected[("Drug", t)]) for t in times
    ]
    assert [_y(line) for line in _lines(bottom)] == [
        pytest.approx(expected[("Dead", t)]) for t in (30.0, 90.0)
    ]
    # the legend of a row names the times of that row
    assert [line["label"] for line in _lines(top)] == [
        "Sim t = 0 min",
        "Sim t = 60 min",
    ]
    assert [line["label"] for line in _lines(bottom)] == [
        "Sim t = 30 min",
        "Sim t = 90 min",
    ]


def test_the_transform_of_a_profile_is_applied() -> None:
    """A frame carries its values in the unit of the profile; the transform applies."""
    result = _flat(plasma=0.8)
    profile = ZoneProfile(
        "Dead cells",
        "%",
        tuple(observers.profile_ids("drg", "plasma", 4)),
        transform=lambda v: 100.0 * (1.0 - v),
    )
    (plot,) = _plots(zone_profiles(result, [profile], times=(5.0,)))
    assert _y(_lines(plot)[0]) == pytest.approx([20.0] * 4)
    assert plot["y"]["label"] == "Dead cells [%]"


def test_the_data_of_a_column_is_drawn_in_that_column_only() -> None:
    profiles = _profiles(1, n=4)
    times = (60.0,)
    columns = [
        _column("A", profiles=profiles, times=times, color="#ff0000"),
        _column("B", profiles=profiles, times=times, color="#0000ff"),
    ]
    data = pd.DataFrame(
        {
            "label": ["p0", "p0"],
            "time": [60.0, 60.0],
            "position": [0.25, 0.75],
            "value": [1.0, 2.0],
            "column": ["A", None],
        }
    )
    first, second = _plots(zone_profiles(columns, profiles, times=times, data=data))
    # the row of the column A and the row of no column in A, the row of no column in B
    assert list(_x(markers_of(first))) == [0.25, 0.75]
    assert list(_x(markers_of(second))) == [0.75]
    # without an sd the markers have no error bars
    assert "yerr" not in markers_of(first)
    # the markers take the colour of their column
    assert markers_of(first)["style"]["color"] == "#ff0000"
    assert markers_of(second)["style"]["color"] == "#0000ff"


def test_data_of_a_column_which_no_column_has_raises() -> None:
    profiles = _profiles(1, n=4)
    times = (60.0,)
    columns = [_column("A", profiles=profiles, times=times)]
    data = pd.DataFrame(
        {
            "label": ["p0"],
            "time": [60.0],
            "position": [0.25],
            "value": [1.0],
            "column": ["Z"],
        }
    )
    with pytest.raises(ValueError, match="Z") as info:
        zone_profiles(columns, profiles, times=times, data=data)
    assert "['A']" in str(info.value)


def test_a_missing_band_raises_and_names_the_pair() -> None:
    profiles = _profiles(1, n=4)
    columns = [_column("A", profiles=profiles, times=(60.0,))]
    with pytest.raises(KeyError) as info:
        zone_profiles(columns, profiles, times=(60.0, 120.0))
    message = str(info.value)
    assert "p0" in message
    assert "120" in message
    assert "A" in message


def test_the_rows_of_a_grid_share_their_y_axis() -> None:
    profiles = _profiles(2, n=4)
    times = (60.0,)
    columns = [
        _column("A", profiles=profiles, times=times, level=1.0),
        _column("B", profiles=profiles, times=times, level=3.0),
    ]
    first, second, third, fourth = _plots(zone_profiles(columns, profiles, times=times))

    def limits(plot: dict[str, Any]) -> tuple[float, float]:
        return (plot["y"]["min"], plot["y"]["max"])

    assert limits(first) == limits(second)
    assert limits(third) == limits(fourth)
    # the limit of a row is the largest value of its columns, the second here
    expected = _medians(profiles, times, 3.0)
    assert limits(first)[1] == pytest.approx(expected[("p0", 60.0)].max() * 1.05)
    assert limits(third)[1] == pytest.approx(expected[("p1", 60.0)].max() * 1.05)
    assert limits(first) != limits(third)


def test_the_legend_entries_are_neutral_where_the_columns_carry_the_colour() -> None:
    profiles = _profiles(1, n=4)
    times = (0.0, 60.0)
    columns = [
        _column("A", profiles=profiles, times=times, color="#ff0000"),
        _column("B", profiles=profiles, times=times, color="#0000ff"),
    ]
    for plot in _plots(zone_profiles(columns, profiles, times=times)):
        # the titles name the columns, so an entry stands for its time, not a column
        swatches = plot["legend_swatches"]
        assert list(swatches) == ["Sim t = 0 min", "Sim t = 60 min"]
        assert [s["color"] for s in swatches.values()] == [NEUTRAL] * 2
        assert [s["line"]["type"] for s in swatches.values()] == ["dotted", "dashed"]


def test_the_data_entry_of_a_row_follows_the_colour_rule() -> None:
    profiles = _profiles(1, n=4)
    times = (60.0,)
    data = _data(label=["p0", "p0", "p0"])
    plain = [_column("A", profiles=profiles, times=times)]
    coloured = [_column("A", profiles=profiles, times=times, color="#ff0000")]
    (plot,) = _plots(zone_profiles(coloured, profiles, times=times, data=data))
    swatch = plot["legend_swatches"]["Data t = 60 min"]
    assert swatch["color"] == NEUTRAL
    assert swatch["marker"]["type"] == "square"
    assert swatch["line"]["type"] == "none"
    (plot,) = _plots(zone_profiles(plain, profiles, times=times, data=data))
    assert "legend_swatches" not in plot
    assert markers_of(plot)["style"]["color"] == _time_colours(1)[0]


def test_a_swatch_names_only_a_group_its_plot_has() -> None:
    """The renderer refuses a swatch of a group of no element of its plot."""
    profiles = _profiles(1, n=4)
    times = (60.0,)
    columns = [
        _column("A", profiles=profiles, times=times, color="#ff0000"),
        _column("B", profiles=profiles, times=times, color="#0000ff"),
    ]
    data = pd.DataFrame(
        {"label": ["p0"], "time": [60.0], "position": [0.5], "value": [1.0]}
    ).assign(column=["B"])
    chart = zone_profiles(columns, profiles, times=times, data=data)
    first, second = _plots(chart)
    assert "Data t = 60 min" not in first["legend_swatches"]
    assert "Data t = 60 min" in second["legend_swatches"]
    _rendered(chart)


# --- the edge cases of a profile: gaps ---------------------------------------


def test_a_profile_with_nan_data_or_a_missing_layer_renders_with_gaps() -> None:
    """A NaN of a band (a layer a condition lacks) or of the data is a gap, no error."""
    profiles = _profiles(1, n=4)
    times = (60.0,)
    gap = _column("A", profiles=profiles, times=times, color="#ff0000", spread=0.5)
    nothing = {key: Band(np.full(4, np.nan)) for key in gap.values}
    lacking = ZoneColumn("B", nothing, "#0000ff")
    band = gap.values[("p0", 60.0)]
    holed = Band(
        np.where(np.arange(4) == 1, np.nan, band.median), band.lower, band.upper
    )
    columns = [ZoneColumn("A", {("p0", 60.0): holed}, "#ff0000"), lacking]
    data = _data(
        label=["p0", "p0"],
        time=[60.0, 60.0],
        position=[0.25, 0.75],
        value=[1.0, np.nan],
        sd=[0.1, 0.2],
    )
    chart = zone_profiles(columns, profiles, times=times, data=data)
    first, second = _plots(chart)
    # the gap of the median is `null` in the spec, a break of its line
    assert _lines(first)[0]["y"][1] is None
    assert np.isfinite(_y(_lines(first)[0])).sum() == 3
    # a layer the condition lacks is a line of gaps, drawn as nothing
    assert all(v is None for v in _lines(second)[0]["y"])
    # the NaN value of the data is a gap of its marker, the other one is drawn
    assert markers_of(first)["y"] == [1.0, None]
    # the scale is the one of the values the figure has
    assert band.upper is not None
    assert first["y"]["max"] == pytest.approx(band.upper.max() * 1.05)
    _rendered(chart)


def test_a_row_of_nothing_but_gaps_is_refused() -> None:
    """No data is no reason to draw a plausible row."""
    profiles = _profiles(1, n=4)
    column = ZoneColumn("A", {("p0", 60.0): Band(np.full(4, np.nan))})
    with pytest.raises(ValueError, match="None of the 4 values is finite"):
        zone_profiles([column], profiles, times=(60.0,))


def test_a_profile_of_a_single_zone_is_drawn_as_a_point() -> None:
    """One zone is the well-stirred limit and must draw as easily as ten.

    A line of one point has no segment, so the panel was empty; the median
    gets a marker and the band of the draws, which a fill between one x does
    not show either, is drawn as the segment it is.
    """
    profiles = _profiles(1, n=1)
    column = _column("one", profiles=profiles, times=(60.0,), spread=0.5)
    (plot,) = _plots(zone_profiles([column], profiles, times=(60.0,)))
    segment, line = _elements(plot, "curve")
    assert line["style"]["marker"]["type"] == "circle"
    assert list(_x(segment)) == [0.5, 0.5]
    median = _y(line)[0]
    assert list(_y(segment)) == pytest.approx([median - 0.5, median + 0.5])
    assert segment["style"]["opacity"] == pytest.approx(0.25)
    assert segment["style"]["line"]["width"] == pytest.approx(6.0)
    (frame,) = _plots(zone_profiles(_result(n=1), _profiles(1, n=1), times=(60,)))
    assert _lines(frame)[0]["style"]["marker"]["type"] == "circle"
    (panel,) = _plots(layer_profile(_result((Layer.PLASMA,), n=1), "drg", time=60))
    assert _elements(panel, "curve")[0]["style"]["marker"]["type"] == "circle"


def test_a_profile_of_many_zones_is_a_line_without_markers() -> None:
    """The marker is for the point which would otherwise not be drawn."""
    (plot,) = _plots(zone_profiles(_result(), _profiles(1), times=(60,)))
    assert _lines(plot)[0]["style"]["marker"]["type"] == "none"


def test_a_row_without_times_of_its_own_says_so() -> None:
    """`times=()` of a profile drew no line and failed on a bare `min()`."""
    profile = ZoneProfile(
        "empty", "mM", tuple(observers.profile_ids("drg", Layer.PLASMA, N)), times=()
    )
    with pytest.raises(ValueError, match="times of the profile 'empty' are empty"):
        zone_profiles(_result(), [profile], times=(60,))


def test_fill_band_fills_between_the_ends_and_a_band_without_ends_nothing() -> None:
    chart = Chart()
    plot = chart.plot()
    x = np.linspace(0.0, 1.0, 5)
    fill_band(plot, x, Band(x), "C0")
    assert plot.spec()["elements"] == []
    fill_band(plot, x, Band(x, x - 0.1, x + 0.1), "C0", alpha=0.4)
    (fill,) = plot.spec()["elements"]
    assert fill["kind"] == "area"
    assert fill["style"]["color"] == "#1f77b4"
    assert fill["style"]["fill"]["opacity"] == pytest.approx(0.4)
    assert fill["from"] == pytest.approx(list(x - 0.1))


# --- the space and the time --------------------------------------------------


def _cells(plot: dict[str, Any]) -> list[tuple[np.ndarray, str]]:
    """The cells of a heatmap: the ring of every cell and its colour."""
    return [
        (np.asarray(ring, dtype=float), polygon["style"]["color"])
        for polygon in _elements(plot, "polygon")
        for ring in polygon["paths"]
    ]


def _cell_at(plot: dict[str, Any], x: float, y: float) -> str:
    """The colour of the cell of a heatmap which holds a point."""
    held = [
        colour
        for ring, colour in _cells(plot)
        if ring[:, 0].min() < x < ring[:, 0].max()
        and ring[:, 1].min() < y < ring[:, 1].max()
    ]
    assert len(held) == 1, f"({x}, {y}) lies in {len(held)} cells"
    return held[0]


def _scale_of(plot: dict[str, Any]) -> ColorScale:
    """The colour scale of the one colour bar of a plot."""
    (bar,) = _colorbars(plot)
    scale = bar["scale"]
    return ColorScale(
        scale["low"],
        scale["high"],
        scale["name"],
        log=scale["log"],
        start=scale.get("start", 0.0),
        end=scale.get("end", 1.0),
    )


def test_the_space_time_panel_carries_the_outflow_below() -> None:
    chart = space_time(_result(), "drg", layers=(Layer.PLASMA,))
    heat, out = _plots(chart)
    assert len(_colorbars(heat)) == 1
    assert _colorbars(out) == []
    assert len(_cells(heat)) == 25 * N
    assert _y(_elements(out, "curve")[0])[-1] == pytest.approx(1.0)
    assert out["x"]["label"] == "Time [min]"
    assert heat["row"] < out["row"]
    # the time is shared: the heatmap keeps the marks of its ticks, the outflow
    # carries the numbers and the label
    assert heat["x"]["ticks_visible"] is False
    assert heat["x"]["tick_marks_visible"] is True
    assert (heat["x"]["min"], heat["x"]["max"]) == (out["x"]["min"], out["x"]["max"])
    _rendered(chart)


def test_the_heatmap_has_the_zones_periportal_first_and_the_time_to_the_right() -> None:
    """The plasma of the fixture is `2 k t / 120` in zone `k` at time `t`."""
    heat = _plots(space_time(_result(), "drg", layers=(Layer.PLASMA,)))[0]
    scale = _scale_of(heat)
    assert (scale.low, scale.high) == (0.0, 2.0 * N)

    def colour(value: float) -> str | None:
        return scale.colors([value])[0]

    # z1 at 120, z10 at 120, z10 at 60 and z1 at 0: the periportal zone is lowest
    assert _cell_at(heat, 119.0, 0.05) == colour(2.0)
    assert _cell_at(heat, 119.0, 0.95) == colour(2.0 * N)
    assert _cell_at(heat, 60.0, 0.95) == colour(N)
    assert _cell_at(heat, 1.0, 0.05) == colour(0.0)
    # the cells of a zone span its tenth, the cells of a time reach half way to
    # its neighbours, and the cells of the ends end at them
    rings = np.vstack([ring for ring, _ in _cells(heat)])
    assert sorted(set(np.round(rings[:, 1], 12))) == pytest.approx(
        list(np.linspace(0.0, 1.0, N + 1))
    )
    edges = sorted(set(np.round(rings[:, 0], 12)))
    assert edges == pytest.approx([0.0, *np.arange(2.5, 120.0, 5.0), 120.0])
    assert (heat["y"]["min"], heat["y"]["max"]) == (0.0, 1.0)
    assert (heat["x"]["min"], heat["x"]["max"]) == (0.0, 120.0)
    # the cells of one colour are one element, so a heatmap stays light
    colours = [p["style"]["color"] for p in _elements(heat, "polygon")]
    assert len(colours) == len(set(colours))
    # a cell is opaque and outlined in its own colour, which closes the seams
    for polygon in _elements(heat, "polygon"):
        assert polygon["style"]["fill"]["opacity"] == 1.0
        assert "color" not in polygon["style"]["fill"]


def test_the_space_time_without_the_outflow_has_a_panel_per_layer() -> None:
    plots = _plots(
        space_time(_result(), "drg", layers=(Layer.PLASMA, Layer.CELL), outflow=False)
    )
    assert len(plots) == 2
    assert [len(_colorbars(p)) for p in plots] == [1, 1]
    # the last heatmap carries the time
    assert plots[-1]["x"]["label"] == "Time [min]"
    assert plots[-1]["x"]["ticks_visible"] is True


def test_the_outflow_panel_is_lower_than_a_heatmap_and_the_name_is_chosen() -> None:
    chart = space_time(_result(), "drg", layers=(Layer.PLASMA,), name="Drug")
    heat, out = _plots(chart)
    assert out["row_span"] < heat["row_span"]
    assert out["y"]["label"] == "Drug outflow [mM]"
    assert _colorbars(heat)[0]["label"] == "Drug plasma [mM]"
    default = _plots(space_time(_result(), "drg", layers=(Layer.PLASMA,)))[0]
    assert _colorbars(default)[0]["label"] == "drg plasma [mM]"
    # as the renderer lays it out: the outflow is about two thirds of a heatmap,
    # and the panels of the column share their x, the colour bar included
    geometry = sbmlsim_core.figure_geometry(chart.to_json())
    rects = {p["id"]: p["rect"] for p in geometry["plots"]}
    heat_rect, out_rect = rects[heat["id"]], rects[out["id"]]
    assert out_rect["height"] / heat_rect["height"] == pytest.approx(0.65, abs=0.1)
    assert heat_rect["x"] == pytest.approx(out_rect["x"])
    assert heat_rect["width"] == pytest.approx(out_rect["width"])
    # the label of the time lies in the figure
    assert out_rect["y"] + out_rect["height"] + 30.0 < geometry["height"]


def test_a_heatmap_of_a_result_with_a_gap_keeps_the_scale_of_its_values() -> None:
    """The scale of a figure of a result with a gap is the one of what it has."""
    result = _result((Layer.PLASMA,), n=4)
    first = observers.profile_ids("drg", Layer.PLASMA, 4)[0]
    result.loc[3, first] = np.nan
    values = result[observers.profile_ids("drg", Layer.PLASMA, 4)].to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    chart = space_time(result, "drg", layers=(Layer.PLASMA,), outflow=False)
    (heat,) = _plots(chart)
    scale = _scale_of(heat)
    assert (scale.low, scale.high) == (float(finite.min()), float(finite.max()))
    # the cell of the gap is left out, the background shows through it
    assert len(_cells(heat)) == 25 * 4 - 1
    _rendered(chart)


def test_a_space_time_heatmap_on_a_log_scale_leaves_its_zeros_out() -> None:
    """A log colour scale spans the positive values; a zero is a gap, not an end."""
    result = _result((Layer.PLASMA, Layer.CELL), n=4)
    for sid in observers.profile_ids("drg", Layer.CELL, 4):
        result[sid] = np.where(
            result["time"] >= 60.0, 10.0 ** (result["time"] / 30.0), 0.0
        )
    chart = space_time(result, "drg", layers=(Layer.CELL,), log=True)
    heat, out = _plots(chart)
    scale = _scale_of(heat)
    assert scale.log
    assert (scale.low, scale.high) == pytest.approx((10.0**2, 10.0**4))
    # the zeros of the times before 60 are no cells, the 13 times from 60 are
    assert len(_cells(heat)) == 13 * 4
    assert min(ring[:, 0].min() for ring, _ in _cells(heat)) == pytest.approx(57.5)
    # a cell has the colour of its value on the log scale
    assert _cell_at(heat, 90.0, 0.6) == scale.colors([10.0**3])[0]
    # the outflow is a line axis, which the scale of the heatmap does not touch
    assert out["y"].get("scale", "linear") == "linear"
    _rendered(chart)


def test_space_time_needs_two_time_points() -> None:
    """One time point gives cells of no width and an axis of min == max (#c)."""
    one = _result().iloc[:1]
    with pytest.raises(ValueError, match="at least two time points"):
        space_time(one, "drg", layers=LAYERS)


def test_a_log_scale_needs_a_positive_value_and_a_positive_end() -> None:
    result = _flat(plasma=0.0)
    with pytest.raises(ValueError, match="No value is positive"):
        space_time(result, "drg", layers=(Layer.PLASMA,), log=True)
    with pytest.raises(ValueError, match="positive lower end"):
        space_time(_result(), "drg", layers=(Layer.PLASMA,), log=True, vmin=0.0)
    # a constant on a log scale spans the decade below it
    (heat, _) = _plots(
        space_time(_flat(plasma=5.0), "drg", layers=(Layer.PLASMA,), log=True)
    )
    scale = _scale_of(heat)
    assert (scale.low, scale.high) == pytest.approx((0.5, 5.0))


# --- the layers --------------------------------------------------------------


def test_the_layer_profile_shows_every_layer_that_carries_substance() -> None:
    plots = _plots(layer_profile(_result(), "drg", time=60))
    assert len(plots) == 6
    assert [p["title"] for p in plots] == [str(layer) for layer in LAYERS]
    for i, plot in enumerate(plots):
        assert _y(_elements(plot, "curve")[0]) == pytest.approx(
            [(i + 1) * k * 0.5 for k in range(1, N + 1)]
        )


def test_the_layer_profile_shows_only_the_layers_present() -> None:
    plots = _plots(layer_profile(_result((ids.VASCULAR, Layer.CELL)), "drg", time=0))
    assert [p["title"] for p in plots] == ["cell", "vas"]


def test_the_layer_profile_of_an_absent_substance_raises() -> None:
    with pytest.raises(KeyError, match="other"):
        layer_profile(_result(), "other", time=0)


def test_the_layer_grid_carries_its_labels_on_the_outer_panels() -> None:
    for layers in (LAYERS, LAYERS[:5]):
        plots = _plots(layer_profile(_result(layers), "drg", time=0))
        # the label of the zones under the last panel of every column, in two
        # lines so that it is no wider than the panel; the concentration left
        # of the first panel of every row
        bottom = [i + 3 >= len(layers) for i in range(len(layers))]
        assert [p["x"].get("label") for p in plots] == [
            POSITION_LINES if b else None for b in bottom
        ]
        assert [p["y"].get("label") for p in plots] == [
            "drg [mM]" if p["column"] == 0 else None for p in plots
        ]


def test_the_layers_share_one_concentration_scale_unless_told_otherwise() -> None:
    result = _flat()

    def limits(chart: Chart) -> dict[str, tuple[float, float]]:
        return {p["title"]: (p["y"]["min"], p["y"]["max"]) for p in _plots(chart)}

    shared = limits(layer_profile(result, "drg", time=5))
    assert all(v == pytest.approx((0.0, 3.15)) for v in shared.values())
    own = limits(layer_profile(result, "drg", time=5, share_y=False))
    assert own["plasma"] == pytest.approx((0.0, 3.15))
    assert own["cell"] == (0.0, 1.0)
    zero = limits(layer_profile(_flat(plasma=0.0), "drg", time=5))
    assert set(zero.values()) == {(0.0, 1.0)}


def test_a_shared_scale_labels_the_ticks_of_the_outer_panels_only() -> None:
    def ticks(plot: dict[str, Any]) -> tuple[bool, bool]:
        return (plot["y"]["ticks_visible"], plot["x"]["ticks_visible"])

    shared = [ticks(p) for p in _plots(layer_profile(_result(), "drg", time=0))]
    assert shared == [
        (True, False),
        (False, False),
        (False, False),
        (True, True),
        (False, True),
        (False, True),
    ]
    five = [ticks(p) for p in _plots(layer_profile(_result(LAYERS[:5]), "drg", time=0))]
    assert five[2] == (False, True)  # above a gap, so on its own bottom edge
    own = _plots(layer_profile(_result(), "drg", time=0, share_y=False))
    assert all(ticks(p) == (True, True) for p in own)


# --- the figures as a whole --------------------------------------------------


def test_the_figures_are_charts_and_save_as_svg_png_and_json(tmp_path: Path) -> None:
    result = _result()
    figures = {
        "profiles": zone_profiles(result, _profiles(2), times=(0, 60)),
        "space_time": space_time(result, "drg", layers=(Layer.PLASMA, Layer.CELL)),
        "layers": layer_profile(result, "drg", time=60),
        "lobule": lobule_map(result, "drg", time=60, zones=N, layer=Layer.CELL),
    }
    for name, chart in figures.items():
        assert isinstance(chart, Chart), name
        for suffix in (".svg", ".png", ".json"):
            path = chart.save(tmp_path / f"{name}{suffix}")
            assert path.stat().st_size > 0, path
        assert "<image" not in (tmp_path / f"{name}.svg").read_text(encoding="utf-8")
        assert (tmp_path / f"{name}.png").read_bytes()[:4] == b"\x89PNG"
        spec = json.loads((tmp_path / f"{name}.json").read_text(encoding="utf-8"))
        assert spec["contract"] == "figure/1"


def test_data_axes_are_the_plots_of_a_chart() -> None:
    """A colour bar is a field of its plot, so every plot is a panel."""
    chart = space_time(_result(), "drg", layers=(Layer.PLASMA, Layer.CELL))
    assert data_axes(chart) == chart.plots
    assert len(data_axes(chart)) == 3


def test_the_figures_draw_the_result_of_a_built_chain(tmp_path: Path) -> None:
    """The ids the figures look for are those the engine emits, over a real chain."""
    n = 4
    settings = _settings(n, 0.01)
    model = _Diffusing(settings, [TRACER], None).build()
    add_observers(model, observer_objects(settings, [TRACER]))
    path = tmp_path / "chain.xml"
    create_model(model, filepath=path, sbml_level=3, sbml_version=2, validate=False)
    r = individual(path)
    r[f"[{ids.species('tr', Layer.PLASMA, ids.zone(1))}]"] = 1.0
    species = list(r.model.getFloatingSpeciesIds())
    selections = ["time", *(f"[{s}]" for s in species)]
    out = observers.outflow_concentration("tr")
    selections.append(out)
    result = result_frame(r.simulate(0, 5, 51, selections=selections))
    assert "time" in result.columns
    assert "tr_plasma_z1" in result.columns

    plasma = observers.profile_ids("tr", Layer.PLASMA, n)
    profile = ZoneProfile("tracer", "mM", tuple(plasma))
    (plot,) = _plots(zone_profiles(result, [profile], times=(0.5, 5)))
    lines = _lines(plot)
    assert len(_x(lines[0])) == n
    assert all(np.isfinite(_y(line)).all() for line in lines)
    assert np.abs(_y(lines[0])).max() > 0

    plots = _plots(layer_profile(result, "tr", time=1.0))
    assert [p["title"] for p in plots] == LITERAL_LAYERS
    curves = [_elements(p, "curve")[0] for p in plots]
    assert all(len(_x(c)) == n and np.isfinite(_y(c)).all() for c in curves)
    assert any(np.abs(_y(c)).max() > 0 for c in curves)

    chart = space_time(result, "tr", layers=(Layer.PLASMA,))
    heat, flow_out = _plots(chart)
    assert len(_cells(heat)) == n * 51
    assert _scale_of(heat).high > 0
    assert np.isfinite(_y(_elements(flow_out, "curve")[0])).all()
    chart.save(tmp_path / "chain.svg")
    assert "<image" not in (tmp_path / "chain.svg").read_text(encoding="utf-8")


# --- the scales --------------------------------------------------------------


def test_a_gap_in_the_data_does_not_collapse_a_scale() -> None:
    """One NaN made every scale the fallback zero to one.

    `np.min` and `np.max` of an array with a NaN are NaN and every comparison
    against NaN is false, so both scales fell through to their no-range
    branch: data from 200 to 1000 was drawn from 0 to 1, every zone saturated
    and the colorbar stated a range the figure did not have.
    """
    whole = np.array([200.0, 500.0, 1000.0])
    gap = np.array([200.0, np.nan, 1000.0])
    assert _limits([gap]) == _limits([whole]) == (0.0, 1000.0)
    assert _scale([gap]) == _scale([whole]) == (200.0, 1000.0)


def test_a_scale_of_nothing_but_nan_is_refused() -> None:
    """No data is no reason to draw a plausible figure."""
    nans = np.array([np.nan, np.nan])
    for scale in (_limits, _scale):
        with pytest.raises(ValueError, match="None of the 2 values is finite"):
            scale([nans])


def test_a_scale_of_one_value_is_never_degenerate() -> None:
    """The renderer needs `low < high`, and a bar of one value states nothing."""
    assert _scale([np.array([2.5, 2.5])]) == (0.0, 2.5)
    assert _scale([np.array([0.0, 0.0])]) == (0.0, 1.0)
    assert _scale([np.array([-2.0, -2.0])]) == (-2.0, 0.0)
    assert _scale([np.array([4.0, 4.0])], log=True) == (0.4, 4.0)
    for low, high in (
        _scale([np.array([2.5])]),
        _scale([np.array([-2.0])]),
        _scale([np.array([0.0])]),
    ):
        assert low < high


def test_the_colour_scale_refuses_what_has_no_range() -> None:
    with pytest.raises(ValueError, match="no range"):
        ColorScale(1.0, 1.0)
    with pytest.raises(ValueError, match="not finite"):
        ColorScale(0.0, float("inf"))
    with pytest.raises(ValueError, match="positive lower end"):
        ColorScale(0.0, 1.0, log=True)
    with pytest.raises(ValueError, match="0 <= start < end <= 1"):
        ColorScale(0.0, 1.0, start=0.5, end=0.5)


def test_the_colour_scale_is_the_one_of_sbmlsim_core_and_its_gaps_none() -> None:
    scale = ColorScale(1.0, 100.0, log=True, start=0.1, end=0.9)
    values = [1.0, 10.0, 100.0]
    assert scale.colors([0.0, *values, np.nan, -1.0]) == [
        None,
        *sbmlsim_core.color_scale(
            "viridis", 1.0, 100.0, values, log=True, start=0.1, end=0.9
        ),
        None,
        None,
    ]
    plot = Chart().plot()
    colorbar(plot, scale, "x [mM]")
    bar = plot.spec()["colorbar"]
    assert bar["label"] == "x [mM]"
    assert bar["scale"] == {
        "name": "viridis",
        "low": 1.0,
        "high": 100.0,
        "log": True,
        "start": 0.1,
        "end": 0.9,
    }


# --- the lobule --------------------------------------------------------------


def _values() -> list[float]:
    """A distinct value per zone, the gradient the maps draw."""
    return [float(k) for k in range(1, N + 1)]


def _scale_n() -> ColorScale:
    """A viridis scale of the values 1 to N."""
    return ColorScale(1.0, float(N))


def _lobule(values: Sequence[float] | None = None, **arguments: Any) -> dict[str, Any]:  # noqa: ANN401
    """The spec of a plot of the lobule a caller drew with `draw_lobule`."""
    chart = Chart()
    x, y = lobule_axes()
    plot = chart.plot(x=x, y=y, legend=False)
    arguments.setdefault("scale", _scale_n())
    draw_lobule(plot, _values() if values is None else values, **arguments)
    (spec,) = _plots(chart)
    return spec


def _bands(plot: dict[str, Any]) -> list[dict[str, Any]]:
    """The polygons of the bands of a lobule, one per zone."""
    return _elements(plot, "polygon")


def _zone_face(plot: dict[str, Any], k: int) -> str:
    """The colour of the band of zone `k`."""
    (band,) = [b for b in _bands(plot) if b["id"] == f"{plot['id']}-{ids.zone(k)}"]
    return band["style"]["color"]


def _inside(ring: np.ndarray, point: tuple[float, float]) -> bool:
    """Whether a point lies inside a ring, by the crossings of a ray to the right."""
    x, y = point
    xs, ys = ring[:, 0], ring[:, 1]
    xn, yn = np.roll(xs, -1), np.roll(ys, -1)
    crossing = (ys > y) != (yn > y)
    at = xs + (y - ys) * (xn - xs) / np.where(yn == ys, 1.0, yn - ys)
    return bool(np.count_nonzero(crossing & (x < at)) % 2)


def _face_at(plot: dict[str, Any], point: tuple[float, float]) -> str:
    """The colour of the band which holds a point."""
    colours = {
        band["style"]["color"]
        for band in _bands(plot)
        if any(_inside(np.asarray(ring), point) for ring in band["paths"])
    }
    assert len(colours) == 1, f"{point} lies in {len(colours)} colours"
    return colours.pop()


def _curves(plot: dict[str, Any]) -> list[dict[str, Any]]:
    """The curves of a lobule: the isolines, the outline and the triads."""
    return _elements(plot, "curve")


def test_draw_lobule_colours_a_band_per_zone() -> None:
    values = _values()
    plot = _lobule(values)
    assert [b["id"] for b in _bands(plot)] == [
        f"{plot['id']}-{ids.zone(k)}" for k in range(1, N + 1)
    ]
    # every zone is one polygon of the twelve wedges of `lobule_bands`
    for band, wedges in zip(_bands(plot), lobule_bands(N), strict=True):
        assert len(band["paths"]) == 12
        for ring, wedge in zip(band["paths"], wedges, strict=True):
            assert np.asarray(ring) == pytest.approx(wedge)
        # opaque, outlined in its own colour against the seams of the wedges
        assert band["style"]["fill"]["opacity"] == 1.0
        assert band["style"]["line"]["type"] == "solid"
    colours = _scale_n().colors(values)
    for k in (1, 4, N):
        assert _zone_face(plot, k) == colours[k - 1]
    # the band of a zone lies where the zone lies: 1 at a corner, N at the centre
    assert _face_at(plot, (0.98, 0.0)) == colours[0]
    assert _face_at(plot, (0.0, 0.0)) == colours[-1]


def test_the_axes_of_a_lobule_hold_it_at_one_scale_without_ticks() -> None:
    for aspect in (0.5, LOBULE_X / LOBULE_Y_BARE, 2.0):
        x, y = lobule_axes(aspect)
        assert x.min is not None and x.max is not None
        assert y.min is not None and y.max is not None
        assert x.min == -x.max
        assert y.min == -y.max
        # a unit of x and a unit of y are one length on a panel of the aspect
        assert (x.max - x.min) / (y.max - y.min) == pytest.approx(aspect)
        # the hexagon and its triads are inside
        assert x.max >= LOBULE_X
        assert y.max >= LOBULE_Y_BARE
        assert min(x.max - LOBULE_X, y.max - LOBULE_Y_BARE) == pytest.approx(
            0.0, abs=1e-12
        )
        for axis in (x, y):
            assert axis.ticks_visible is False
            assert axis.tick_marks_visible is False
    # the triads stand outside the corners, the name of the portal triad below
    assert LOBULE_X > 1.0 + 0.045
    assert abs(_PT[1]) < LOBULE_Y
    with pytest.raises(ValueError, match="positive"):
        lobule_axes(0.0)


def test_the_axes_of_a_lobule_reserve_room_for_the_names_only_when_it_has_them() -> (
    None
):
    """The hexagon fills the panel height unless `PT` stands below it (#a)."""
    _, bare = lobule_axes()
    _, named = lobule_axes(labels=True)
    assert bare.max == pytest.approx(LOBULE_Y_BARE)
    assert named.max == pytest.approx(LOBULE_Y)
    # the triad dots reach 0.91, the hexagon 0.866: inside the bare axes
    assert 0.91 < LOBULE_Y_BARE < 0.95
    # the map has no names, so its hexagon fills its panel
    (plot,) = _plots(lobule_map(_result(), "drg", time=60.0, zones=N))
    assert plot["y"]["max"] < LOBULE_Y_BARE * 1.5  # widened only for the aspect
    outline = _curves(plot)[0]
    assert np.abs(_y(outline)).max() / plot["y"]["max"] > 0.85


def test_the_name_of_the_portal_triad_is_clear_of_its_triad_dots() -> None:
    """The `PT` glyph of the labelled lobule must not touch the green dot (#b)."""
    spec = _lobule(labels=True, anatomy=True)
    (text,) = [e for e in _elements(spec, "text") if e.get("text") == "PT"]
    lowest = min(_y(e).min() for e in _curves(spec)[1:])
    assert text["y"] < lowest - 0.08
    assert abs(text["y"]) < LOBULE_Y
    x, _ = lobule_axes(label="below")
    assert (x.label, x.label_visible) == ("below", True)


def test_draw_lobule_draws_a_zone_without_data_in_the_missing_colour() -> None:
    values = _values()
    for k in (4, 5, 6, 7):
        values[k - 1] = float("nan")
    plot = _lobule(values)
    colours = _scale_n().colors(_values())
    for k in (4, 7):
        # neither left out, which is a hole in the hexagon, nor a value
        assert _zone_face(plot, k) == MISSING
    # the zones with data keep the colour of their own value, not of a neighbour
    for k in (1, 3, 8, N):
        assert _zone_face(plot, k) == colours[k - 1]


def test_the_colour_of_a_zone_without_data_is_the_callers_to_choose() -> None:
    values = _values()
    values[2] = float("nan")
    plot = _lobule(values, missing="blue")
    assert _zone_face(plot, 3) == "#0000ff"
    assert _zone_face(plot, 4) == _scale_n().colors([4.0])[0]


def test_draw_lobule_draws_the_thirds_only_when_asked() -> None:
    assert _curves(_lobule()) == []
    thirds = _curves(_lobule(thirds=True))
    assert len(thirds) == 2
    for line, position in zip(thirds, (1.0 / 3.0, 2.0 / 3.0), strict=True):
        assert np.column_stack([_x(line), _y(line)]) == pytest.approx(
            lobule_isoline(position)
        )
        assert line["style"]["line"] == {"type": "dotted", "width": 0.7}
        assert line["style"]["color"] == "#404040"
    # the ring of 1/3 lies wholly outside the ring of 2/3, in that order
    radii = [np.hypot(_x(line), _y(line)) for line in thirds]
    assert radii[0].min() == pytest.approx(2.0 / 3.0)  # towards a portal triad
    assert radii[1].min() == pytest.approx(1.0 / 3.0)
    assert radii[0].min() > radii[1].max()


def test_draw_lobule_draws_the_triads_only_when_asked() -> None:
    outline, *triads = _curves(_lobule(anatomy=True))
    assert np.column_stack([_x(outline), _y(outline)]) == pytest.approx(
        lobule_isoline(0.0)
    )
    assert outline["style"]["color"] == "#262626"
    assert outline["style"]["line"]["width"] == pytest.approx(1.1)
    assert [line["style"]["color"] for line in triads] == list(TRIAD_COLORS)
    points = np.vstack([np.column_stack([_x(t), _y(t)]) for t in triads])
    assert len(points) == 18  # three dots at every one of the six corners
    angles = np.deg2rad(60.0 * np.arange(6))
    corners = np.column_stack([np.cos(angles), np.sin(angles)])
    distance = np.linalg.norm(points[:, None, :] - corners[None, :, :], axis=2)
    assert distance.min(axis=1) == pytest.approx(0.045, abs=1e-12)
    nearest = distance.argmin(axis=1)
    assert sorted(np.bincount(nearest, minlength=6)) == [3] * 6
    # the three vessels of a triad stand 120 degrees apart around their corner
    for corner in range(6):
        triad = points[nearest == corner] - corners[corner]
        around = np.sort(np.rad2deg(np.arctan2(triad[:, 1], triad[:, 0])) % 360.0)
        assert np.diff(np.append(around, around[0] + 360.0)) == pytest.approx(120.0)
    for line in triads:
        assert line["style"]["line"]["type"] == "none"
        marker = line["style"]["marker"]
        assert marker["type"] == "circle"
        assert marker["edge"] == "#262626"
        assert marker["edge_width"] == pytest.approx(0.3)
    # every dot lies in the limits of the axes, so none is cut in half
    assert np.abs(points[:, 0]).max() + 0.02 < LOBULE_X
    assert np.abs(points[:, 1]).max() + 0.02 < LOBULE_Y


def test_draw_lobule_names_the_vessels_only_when_asked() -> None:
    assert _elements(_lobule(anatomy=True), "text") == []
    centre, triad = _elements(_lobule(anatomy=True, labels=True), "text")
    assert (centre["text"], triad["text"]) == ("CV", "PT")
    assert (centre["x"], centre["y"]) == (0.0, 0.0)
    assert centre["anchor"] == triad["anchor"] == "middle"
    # below the corner at 300 degrees, clear of its three dots
    assert triad["x"] == 0.5
    assert triad["y"] < -np.cos(np.pi / 6.0) - 0.045 - 0.08
    # the name of the central vein reads on the map: dark on a light centre
    # (the high end of viridis, yellow), light on a dark one
    assert centre["style"]["color"] == "#000000"
    dark = _elements(_lobule(_values()[::-1], labels=True), "text")[0]
    assert dark["style"]["color"] == "#ffffff"


def test_the_points_of_a_panel_are_the_callers_to_choose() -> None:
    def size(plot: dict[str, Any]) -> int:
        return sum(len(ring) for band in _bands(plot) for ring in band["paths"])

    coarse, fine = _lobule(points=4), _lobule(points=64)
    assert size(coarse) < size(fine)
    assert len(_bands(coarse)) == len(_bands(fine)) == N


def test_lobule_map_draws_its_panel_with_the_renderer() -> None:
    """The map draws its panel with `draw_lobule`: the bands, the colours, the title."""
    chart = lobule_map(_result(), "drg", time=60, zones=N)
    (plot,) = _plots(chart)
    assert len(_bands(plot)) == N
    # the cell layer is the fifth: the value of zone k at the time 60 is 2.5 k
    scale = _scale_of(plot)
    assert (scale.low, scale.high) == (2.5, 2.5 * N)
    low, high = sbmlsim_core.color_scale("viridis", 0.0, 1.0, [0.0, 1.0])
    assert _zone_face(plot, 1) == low
    assert _zone_face(plot, N) == high
    assert _face_at(plot, (0.98, 0.0)) == _zone_face(plot, 1)
    assert _face_at(plot, (0.0, 0.0)) == _zone_face(plot, N)
    assert plot["title"] == "drg cell at t = 60 min"
    assert "portal triads (0, periportal)" in plot["x"]["label"]
    assert _colorbars(plot)[0]["label"] == "drg cell [mM]"
    # and it draws the anatomy: the outline and the three dots of every triad
    assert len(_curves(plot)) == 4
    assert _elements(plot, "text") == []
    assert plot["legend"] is False


def test_lobule_map_draws_a_regular_hexagon() -> None:
    """The renderer has no equal aspect: the limits make a unit one length both ways."""
    chart = lobule_map(_result(), "drg", time=60, zones=N, name="Drug")
    (plot,) = _plots(chart)
    (geometry,) = sbmlsim_core.figure_geometry(chart.to_json())["plots"]
    rect = geometry["rect"]
    per_x = rect["width"] / (plot["x"]["max"] - plot["x"]["min"])
    per_y = rect["height"] / (plot["y"]["max"] - plot["y"]["min"])
    assert per_x == pytest.approx(per_y, rel=1e-6)
    assert panel_aspects(chart) == {
        plot["id"]: pytest.approx(rect["width"] / rect["height"])
    }
    _rendered(chart)


def test_lobule_map_names_its_ends_and_the_unit_and_has_one_colorbar() -> None:
    (plot,) = _plots(
        lobule_map(_result(), "drg", time=60, zones=N, unit="uM", name="Drug")
    )
    text = plot["x"]["label"] + plot["title"]
    assert "periportal" in text
    assert "pericentral" in text
    assert "Drug" in text
    (bar,) = _colorbars(plot)
    assert bar["label"] == "Drug cell [uM]"


def test_the_map_of_a_layer_reads_that_layer() -> None:
    result = _result((Layer.PLASMA,))
    (plot,) = _plots(lobule_map(result, "drg", time=60, zones=N, layer=Layer.PLASMA))
    assert len(_bands(plot)) == N
    with pytest.raises(KeyError):
        lobule_map(result, "drg", time=60, zones=N)


def test_a_lobule_of_one_value_has_a_colour_bar_with_a_range() -> None:
    """One zone, or a constant, is one colour on a scale which is no point."""
    one = _result((Layer.CELL,), n=1)
    chart = lobule_map(one, "drg", time=60, zones=1)
    (plot,) = _plots(chart)
    scale = _scale_of(plot)
    value = 1.0 * 60.0 / 120.0  # the first layer, zone 1, at 60
    assert (scale.low, scale.high) == (0.0, value)
    assert _zone_face(plot, 1) == scale.colors([value])[0]
    _rendered(chart)
    flat = _flat(cell=-2.0)
    (plot,) = _plots(lobule_map(flat, "drg", time=5, zones=4))
    scale = _scale_of(plot)
    assert (scale.low, scale.high) == (-2.0, 0.0)
    (plot,) = _plots(lobule_map(_flat(cell=0.0), "drg", time=5, zones=4))
    scale = _scale_of(plot)
    assert (scale.low, scale.high) == (0.0, 1.0)


def test_vmin_and_vmax_override_the_data_range_and_default_to_it() -> None:
    result = _result((Layer.CELL,))
    (plot,) = _plots(lobule_map(result, "drg", time=60, zones=N))
    assert (_scale_of(plot).low, _scale_of(plot).high) == (0.5, 0.5 * N)
    (plot,) = _plots(lobule_map(result, "drg", time=60, zones=N, vmin=1.0, vmax=4.0))
    assert (_scale_of(plot).low, _scale_of(plot).high) == (1.0, 4.0)
    # a value beyond an end takes the colour of that end
    assert _zone_face(plot, N) == _scale_of(plot).colors([4.0])[0]
    assert _zone_face(plot, 1) == _scale_of(plot).colors([1.0])[0]
    (heat,) = _plots(
        space_time(
            result, "drg", layers=(Layer.CELL,), outflow=False, vmin=0.0, vmax=3.0
        )
    )
    assert (_scale_of(heat).low, _scale_of(heat).high) == (0.0, 3.0)
    with pytest.raises(ValueError, match="no range"):
        lobule_map(result, "drg", time=60, zones=N, vmin=5.0, vmax=5.0)


def test_a_lobule_of_another_number_of_zones_than_the_result_is_refused() -> None:
    """`zones=5` on a ten zone result drew five periportal zones as a lobule."""
    with pytest.raises(ValueError, match="holds 10 zones"):
        lobule_map(_result(), "drg", time=60, zones=5)
    assert isinstance(lobule_map(_result(), "drg", time=60, zones=N), Chart)


def test_a_grid_of_lobules_takes_one_colour_bar_and_stays_light(tmp_path: Path) -> None:
    """A model which composes its own grid of maps: the pattern of the docs."""
    scale = _scale_n()

    def lobules(aspects: dict[str, float]) -> Chart:
        chart = Chart(460.0, 230.0, columns=2)
        for column, values in enumerate((_values(), _values()[::-1])):
            x, y = lobule_axes(aspects.get(f"p{column}"))
            plot = chart.plot(0, column, x=x, y=y, legend=False)
            draw_lobule(plot, values, scale=scale, anatomy=True)
        colorbar(plot, scale, "Drug [mM]")
        return chart

    first = lobules({})
    aspects = panel_aspects(first)
    chart = lobules(aspects)
    # the layout does not depend on the limits of the lobules, so the second
    # chart has the panels the first one measured
    assert panel_aspects(chart) == pytest.approx(aspects)
    left, right = _plots(chart)
    assert _colorbars(left) == []
    assert _colorbars(right)[0]["label"] == "Drug [mM]"
    # the element ids of the bands are unique over the panels
    ids_ = [e["id"] for p in (left, right) for e in p["elements"]]
    assert len(ids_) == len(set(ids_))
    path = chart.save(tmp_path / "grid.svg")
    assert "<image" not in path.read_text(encoding="utf-8")
    # the budget of this figure, two lobules of ten zones with one colour bar
    # at the default sampling of the bands; `points=256` would be a regression
    assert path.stat().st_size < 260 * 1024


def test_a_missing_colour_which_draws_nothing_is_refused() -> None:
    with pytest.raises(ValueError, match="draws nothing"):
        _lobule([float("nan")] * N, missing="none")
