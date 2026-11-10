"""The figures of a zonated model: the profiles over the zones, in space and time.

The functions draw a simulation result of any zonated model, whatever built it,
and know nothing but the id scheme of the engine (`pkpdsim.zonation.ids`,
`pkpdsim.zonation.observers`). The result is a `pandas.DataFrame` with a `time`
column and a column per selection named by its plain id (`glc_cell_z2`);
`result_frame` converts a roadrunner result. The number of zones and the layers
which carry a substance are read from the ids present in the result.

`zone_profiles` draws a grid of conditions as well: a `ZoneColumn` per
condition carries plain arrays, the `Band` of every profile and time, so the
figures stay free of the result objects of a simulator.
`pkpdsim.zonation.results` is the bridge which reads them out of a
`sbmlsim.result.scan.ScanResult`, with the draws of an adopted fit as the band.

A figure is built directly as a `matplotlib.figure.Figure` and owned by the
caller; `matplotlib.pyplot` and its global state are not used. Every figure is
vector: the heatmap is not rasterized and neither is the solid of its colorbar,
so the SVG holds no embedded image (a viewer which drops embedded raster data
would show an empty panel). `colorbar` is the one way to that guarantee, here
and in a model, and `draw_lobule` the renderer of the lobule, which `lobule_map`
draws its one panel with and a model its own grid of lobules.
"""

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import Colorbar
from matplotlib.colors import Normalize, to_rgba
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.patheffects import withStroke
from matplotlib.typing import LineStyleType

from pkpdsim.zonation import ids, observers
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.render import lobule_bands, lobule_isoline

__all__ = [
    "MISSING",
    "Band",
    "ZoneColumn",
    "ZoneProfile",
    "colorbar",
    "data_axes",
    "draw_lobule",
    "layer_profile",
    "lobule_map",
    "result_frame",
    "space_time",
    "zone_profiles",
]

#: the group id of the axes of a colorbar
COLORBAR = "colorbar"

_LINE_STYLES: tuple[LineStyleType, ...] = (":", "--", "-", "-.")
_MARKERS: tuple[str, ...] = ("s", "D", "o", "^")
_LAYER_ORDER: tuple[str, ...] = (*(str(layer) for layer in Layer), ids.VASCULAR)
POSITION = "Position (0 periportal, 1 pericentral)"
"""The label of the axis of the zones, which names its ends: a flip must be visible."""
DEFAULT_TIME_UNIT = "min"
"""The unit of the time of the engine (the rates of the engine are per minute)."""
HEATMAP_HEIGHT = 1.0
OUTFLOW_HEIGHT = 0.6
"""The height of the outflow panel of `space_time` relative to a heatmap."""
_HEADROOM = 0.05
DEFAULT_UNIT = "mM"
"""The unit of the concentrations of the profile observers (`observers.profile_ids`)."""
_BAND_ALPHA = 0.25
"""The opacity of the fill of a band of the draws."""
_NEUTRAL = "0.3"
"""The colour of a legend handle which stands for a time in every column."""
_TRIAD_OFFSET = 0.045
"""The distance of a dot of a portal triad from its corner, of the circumradius."""
_TRIAD_COLORS: tuple[str, ...] = ("#1f77b4", "#d62728", "#2ca02c")
"""The colours of the portal vein, the hepatic artery and the bile duct of a triad."""
_THIRDS: tuple[float, ...] = (1.0 / 3.0, 2.0 / 3.0)
"""The lobular positions which bound the thirds a histology counts its zones in."""
_ANATOMY_GREY = "0.15"
"""The grey of the outline of the hexagon and of the edge of a dot of a triad."""
_THIRDS_GREY = "0.25"
"""The grey of the dotted isolines of the thirds."""
_VESSEL_SIZE = 7
"""The font size of the names of the vessels (`CV`, `PT`) on a lobule."""
MISSING = "#f0f0f0"
"""The colour of a zone without data on a lobule, the default of `draw_lobule`."""
NO_LEGEND = "_nolegend_"
"""The label matplotlib leaves out of a legend, the lines of a grid of columns.

The grid draws the handles of its legend itself (`_row_legend`), one set per
row and not one per column, so the lines and markers of a panel carry this
label: a caller who calls `ax.legend()` on a panel gets nothing, not a second
legend of the colours of one column.
"""


@dataclass(frozen=True)
class ZoneProfile:
    """A profile over the zones, a row of `zone_profiles`.

    `transform` and `source_unit` belong to the reading of a result
    (`pkpdsim.zonation.results.zone_column`), which reads the selections in
    `source_unit or unit` and transforms the values into `unit`: a living
    fraction read as dimensionless becomes the percentage of dead cells with
    `transform=lambda v: 100.0 * (1.0 - v)`. A frame carries its values in the
    unit of the profile already, so `zone_profiles` of a frame ignores
    `source_unit` and applies `transform` to the interpolated values.

    Attributes:
        label: the name of the quantity, the label of the y axis.
        unit: the unit of the quantity, the unit of the values which are drawn.
        selections: the ids of the result by zone, periportal first; their
            number is the number of zones of the row.
        times: the times of this row, `None` the times of the figure.
        transform: the transformation of the values into `unit`, `None` none.
        source_unit: the unit the selections are read in, `None` is `unit`.
    """

    label: str
    unit: str
    selections: tuple[str, ...]
    times: tuple[float, ...] | None = None
    transform: Callable[[np.ndarray], np.ndarray] | None = None
    source_unit: str | None = None


@dataclass(frozen=True, eq=False)
class Band:
    """The values over the zones or over the time: the median and its band.

    A simulation without draws carries its values as the median and has no
    ends, so it draws a line and no fill. A transformation which decreases
    swaps the ends (the lower draws give the larger value), which fills the
    same area.

    The fields are arrays, which have no one truth value, so the comparison is
    the identity of the object (`eq=False`) and not a field by field one,
    which would raise `ValueError: The truth value of an array ... is
    ambiguous`: compare the fields with `numpy.allclose`.

    Attributes:
        median: the values, the median over the draws of a run over them.
        lower: the lower quantile of the draws, `None` without draws.
        upper: the upper quantile of the draws, `None` without draws.
    """

    median: np.ndarray
    lower: np.ndarray | None = None
    upper: np.ndarray | None = None


@dataclass(frozen=True)
class ZoneColumn:
    """A column of `zone_profiles`: one condition, its values by profile and time.

    Attributes:
        label: the name of the condition, the title of its panel in the top
            row; empty draws no title.
        values: the band of every profile and time of the figure, by the label
            of the profile and the time.
        color: the colour of every line of the column, `None` the colour of
            the time.
    """

    label: str
    values: Mapping[tuple[str, float], Band]
    color: str | None = None


def result_frame(result: object) -> pd.DataFrame:
    """Convert a roadrunner result into the frame the figures draw.

    Args:
        result: the result of `RoadRunner.simulate` with explicit selections
            (a `NamedArray`, which carries `colnames`); `time` must be selected.

    Returns:
        The frame with a column per selection, `[glc]` named `glc`.
    """
    names = [str(c).removeprefix("[").removesuffix("]") for c in result.colnames]  # ty: ignore[unresolved-attribute]
    return pd.DataFrame(np.asarray(result), columns=names)


def data_axes(figure: Figure) -> list[Axes]:
    """Get the axes of the panels of a figure, the axes of a colorbar left out.

    A caller which post-processes a returned figure (a title, a limit, a
    legend) needs to know which axes are panels: the colorbars of the
    heatmaps carry the group id `COLORBAR`, and this function is the
    documented way to tell them apart.
    """
    return [ax for ax in figure.axes if ax.get_gid() != COLORBAR]


def _positions(n: int) -> np.ndarray:
    """Get the midpoints of `n` zones from 0 (periportal) to 1 (pericentral)."""
    return (np.arange(1, n + 1) - 0.5) / n


def _columns(result: pd.DataFrame, selections: Sequence[str]) -> None:
    """Raise a KeyError which names the selections the result does not hold."""
    missing = [s for s in selections if s not in result.columns]
    if missing:
        raise KeyError(f"The result has no column {missing}.")


def _at(result: pd.DataFrame, selections: Sequence[str], time: float) -> np.ndarray:
    """Get the values of the selections at a time, interpolated in time."""
    _columns(result, ["time", *selections])
    t = result["time"].to_numpy(dtype=float)
    return np.array(
        [np.interp(time, t, result[s].to_numpy(dtype=float)) for s in selections]
    )


def _label(label: str, unit: str) -> str:
    """Get the label of an axis with its unit."""
    return label if unit in {"", "dimensionless"} else f"{label} [{unit}]"


def colorbar(
    figure: Figure,
    mappable: ScalarMappable,
    ax: Axes | Sequence[Axes],
    label: str | None = None,
    **kwargs: Any,  # noqa: ANN401
) -> Colorbar:
    """Add a colorbar which is vector and is no data axes.

    The one way a figure of this package or of a model gets a colorbar: the
    solids of a colorbar are a rasterized image by default, which an SVG viewer
    that drops embedded raster data shows as an empty bar, so they are
    un-rasterized here (`tests/zonation/test_figures.py` looks for `<image` in
    the saved figures). The axes of the bar carry the group id `COLORBAR`, so
    `data_axes` tells them from the panels.

    Args:
        figure: the figure the colorbar belongs to.
        mappable: the scale which is drawn, a `ScalarMappable` of a norm and a
            colormap or an artist which carries one (a heatmap).
        ax: the panel the colorbar is taken from, or the panels of a grid it
            spans.
        label: the label of the bar, `None` none.
        **kwargs: further arguments of `matplotlib.figure.Figure.colorbar`
            (`location`, `shrink`, `aspect`, ... of a bar beside a grid).

    Returns:
        The colorbar.
    """
    if label is not None:
        kwargs["label"] = label
    bar = figure.colorbar(mappable, ax=ax, **kwargs)
    if bar.solids is not None:
        bar.solids.set_rasterized(False)
    bar.ax.set_gid(COLORBAR)
    return bar


def _limits(
    values: Sequence[np.ndarray], headroom: float = 0.0
) -> tuple[float, float] | None:
    """Get the limits of a concentration axis, from zero if nothing is negative.

    A concentration is never negative, so an axis of values at or above zero
    starts at zero; a layer which is zero everywhere gets zero to one. Values
    below zero leave the limits to matplotlib (`None`). A colour scale does not
    start at zero, it spans its values (`_scale`).
    """
    low = min(float(np.min(v)) for v in values)
    high = max(float(np.max(v)) for v in values)
    if low < 0.0:
        return None
    return (0.0, high * (1.0 + headroom) if high > 0.0 else 1.0)


def _scale(
    values: Sequence[np.ndarray], vmin: float | None = None, vmax: float | None = None
) -> tuple[float, float] | None:
    """Get the limits of a colour scale: the range of the values it draws.

    A colour scale spans the range of its values, a line axis keeps its zero
    baseline (`_limits`): the zero is explicit on the ticks of a colorbar, and
    a scale from zero would draw a gradient of a few percent as one shade. A
    scale of no range (a constant, zeros) takes the limits of `_limits`, from
    zero to the value or to one; values below zero with no range leave it to
    matplotlib (`None`). An explicit `vmin` or `vmax` replaces that end of the
    range, the way to clip a scale which a spike dominates; a `vmin` equal to
    the value of a constant therefore raises (no range), while the default
    call for the same data returns the fallback (0, 1) or (0, value).

    Raises:
        ValueError: if the limits are not `vmin < vmax`.
    """
    low = min(float(np.min(v)) for v in values)
    high = max(float(np.max(v)) for v in values)
    if vmin is not None or vmax is not None:
        low = low if vmin is None else vmin
        high = high if vmax is None else vmax
        if not low < high:
            raise ValueError(f"The colour scale {low} to {high} has no range.")
        return (low, high)
    if high > low:
        return (low, high)
    return _limits(values)


def _time_text(time: float, time_unit: str) -> str:
    """Get a time with its unit, `t = 0.5 min`."""
    return f"t = {time:g} {time_unit}"


def _colors(count: int) -> np.ndarray:
    """Get the colour of every time, a sequential map from early to late."""
    return colormaps["viridis"](np.linspace(0.1, 0.75, count))


def _row_times(profile: ZoneProfile, times: Sequence[float]) -> list[float]:
    """Get the times of a row: its own times, else the times of the figure."""
    return [float(t) for t in (times if profile.times is None else profile.times)]


def _frame_column(
    result: pd.DataFrame, profiles: Sequence[ZoneProfile], times: Sequence[float]
) -> ZoneColumn:
    """Get the one column of a frame: the values at the times, no label, no colour.

    A frame carries its values in the unit of the profile, so `source_unit` is
    ignored; `transform` is applied to the interpolated values.

    Raises:
        KeyError: if the result has no column of a selection.
    """
    values: dict[tuple[str, float], Band] = {}
    for profile in profiles:
        for time in _row_times(profile, times):
            y = _at(result, profile.selections, time)
            values[(profile.label, time)] = Band(
                y if profile.transform is None else profile.transform(y)
            )
    return ZoneColumn("", values)


def _band_of(column: ZoneColumn, label: str, time: float) -> Band:
    """Get the band of a profile at a time in a column.

    Raises:
        KeyError: naming the pair, if the column has no band of it.
    """
    try:
        return column.values[(label, time)]
    except KeyError:
        raise KeyError(
            f"The column '{column.label}' has no band of the profile '{label}' "
            f"at the time {time:g}."
        ) from None


def _fill_band(ax: Axes, x: np.ndarray, band: Band, color: str | np.ndarray) -> None:
    """Fill a band between its ends; a band without ends fills nothing."""
    if band.lower is None or band.upper is None:
        return
    ax.fill_between(x, band.lower, band.upper, color=color, alpha=_BAND_ALPHA, lw=0.0)


def _check_data(
    data: pd.DataFrame, profiles: Sequence[ZoneProfile], columns: Sequence[ZoneColumn]
) -> None:
    """Check the labels of the data against the profiles and the columns.

    Raises:
        ValueError: naming them, if a row has a label which no profile has or
            a column which no column of the figure has.
    """
    labels = [profile.label for profile in profiles]
    unknown = sorted(set(data["label"]) - set(labels))
    if unknown:
        raise ValueError(
            f"The data has the labels {unknown} which no profile has; "
            f"the labels of the profiles are {labels}."
        )
    if "column" not in data:
        return
    names = [column.label for column in columns]
    unknown = sorted(set(data["column"].dropna()) - set(names))
    if unknown:
        raise ValueError(
            f"The data has the column labels {unknown} which no column has; "
            f"the labels of the columns are {names}."
        )


def _data_rows(
    data: pd.DataFrame, profile: ZoneProfile, column: ZoneColumn
) -> pd.DataFrame:
    """Get the rows of the data of a profile which a column draws.

    A row with a `column` is drawn in the column of that label only, a row
    without one (no `column` at all, or no label in it) in every column.
    """
    rows = data[data["label"] == profile.label]
    if "column" in rows:
        rows = rows[rows["column"].isna() | (rows["column"] == column.label)]
    return rows


def _draw_data(
    ax: Axes,
    rows: pd.DataFrame,
    times: Sequence[float],
    colors: np.ndarray,
    column: ZoneColumn,
    time_unit: str,
    *,
    labelled: bool,
) -> dict[float, np.ndarray]:
    """Draw the data of a row in a column, the markers of a time.

    Args:
        ax: the panel.
        rows: the rows of the data of this profile and column.
        times: the times of the row.
        colors: the colour of every time.
        column: the column, whose colour the markers take where it has one.
        time_unit: the unit of the times.
        labelled: name the markers in the legend of the panel; the grid draws
            its own handles (`_row_legend`), see `zone_profiles`.

    Returns:
        The drawn values plus their error, by the time they are drawn at; a
        time of which no row is drawn is absent.
    """
    drawn: dict[float, np.ndarray] = {}
    for k, time in enumerate(times):
        at = rows[np.isclose(rows["time"], time)]
        if at.empty:
            continue
        err = at["sd"].to_numpy(dtype=float) if "sd" in at else 0.0
        drawn[time] = at["value"].to_numpy(dtype=float) + err
        text = f"Data {_time_text(time, time_unit)}"
        ax.errorbar(
            at["position"],
            at["value"],
            yerr=at.get("sd"),
            linestyle="",
            marker=_MARKERS[k % len(_MARKERS)],
            markeredgecolor="black",
            color=colors[k] if column.color is None else column.color,
            label=text if labelled else NO_LEGEND,
        )
    return drawn


def _row_legend(
    ax: Axes,
    times: Sequence[float],
    colors: np.ndarray,
    *,
    by_time: bool,
    data_at: Sequence[float],
    time_unit: str,
) -> None:
    """Draw the legend of a row right of it: its times and the data of its times.

    The columns are named by the titles of the top row, so the handles stand
    for the times: in the colour of the time where the colours are the times
    (no column has a colour), neutral where a column has one.
    """
    handles: list[Line2D] = []
    labels: list[str] = []
    for k, time in enumerate(times):
        color = colors[k] if by_time else _NEUTRAL
        handles.append(
            Line2D(
                [],
                [],
                color=color,
                linestyle=_LINE_STYLES[k % len(_LINE_STYLES)],
                linewidth=2.0,
            )
        )
        labels.append(f"Sim {_time_text(time, time_unit)}")
        if time in data_at:
            handles.append(
                Line2D(
                    [],
                    [],
                    color=color,
                    linestyle="",
                    marker=_MARKERS[k % len(_MARKERS)],
                    markeredgecolor="black",
                )
            )
            labels.append(f"Data {_time_text(time, time_unit)}")
    ax.legend(
        handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0
    )


def zone_profiles(
    result: pd.DataFrame | Sequence[ZoneColumn],
    profiles: Sequence[ZoneProfile],
    *,
    times: Sequence[float],
    data: pd.DataFrame | None = None,
    time_unit: str = DEFAULT_TIME_UNIT,
) -> Figure:
    """Draw a row per profile and a column per condition, a line per time.

    A `pandas.DataFrame` is one condition: one column, no title, a colour and
    a line style per time (a sequential map, early to late), the data of a time
    in the colour of that time, and one legend for the figure, right of the
    grid, if every row has the same entries, else one per row. A sequence of
    `ZoneColumn` is a grid of conditions: the label of a column is the title of
    its panel in the top row, the colour of a line is the colour of its column
    where it has one and the line style stays the style of the time, so a
    reader tells the times apart inside a column and the columns apart by
    colour. A band of the draws of a fit is filled between its ends; the
    legend is then one per row, right of the last column, and names the times
    of that row. `pkpdsim.zonation.results.zone_column` builds a column from a
    result of sbmlsim.

    The x axis is shared and the y axis is shared per row, over the values of
    every column of that row. A line axis starts at zero when no value is
    negative; a colour scale (`space_time`, `lobule_map`) spans the range of
    its values instead.

    Args:
        result: the simulation result as a frame, or a column per condition.
        profiles: the profiles, a row each; the number of zones is the number
            of selections of a profile. A profile with `times` draws those
            times in its row, including its data and its legend.
        times: the times of the lines, of a row which has none of its own;
            interpolated in a frame, read at the nearest point of a column
            (`results.zone_values`).
        data: optional data, a long frame with the columns `label` (of a
            profile), `time`, `position` (0 periportal to 1 pericentral),
            `value` (in the unit of the profile) and optionally `sd` (a
            symmetric error in the same unit) and `column` (the label of a
            column, a row without one is drawn in every column), drawn over
            the lines as markers by time. The rows of a time which the row
            does not draw are skipped, as plotting a subset of the times of a
            dataset is legitimate.
        time_unit: the unit of the times, the convention of the engine (min).

    Returns:
        The figure.

    Raises:
        KeyError: if the result has no column of a selection, or if a column
            has no band of a profile and a time it is asked for. A frame is
            read before the data is checked, so a missing selection is named
            before a label of the data which no profile has.
        ValueError: if `profiles`, `times` or `result` is empty, if a row of
            `data` has a label which no profile or no column has, or if no row
            of `data` is drawn.
    """
    if not profiles:
        raise ValueError("The profiles are empty, there is no row to draw.")
    if not len(times):
        raise ValueError("The times are empty, there is no line to draw.")
    single = isinstance(result, pd.DataFrame)
    columns = [_frame_column(result, profiles, times)] if single else list(result)
    if not columns:
        raise ValueError("The columns are empty, there is no condition to draw.")
    if data is not None:
        _check_data(data, profiles, columns)
    n_rows, n_columns = len(profiles), len(columns)
    by_time = all(column.color is None for column in columns)
    figure = Figure(
        figsize=(7.0 if single else 3.6 * n_columns + 2.4, 3.0 * n_rows),
        layout="constrained",
    )
    axes = figure.subplots(n_rows, n_columns, sharex=True, sharey="row", squeeze=False)
    drawn = 0
    for row, profile in enumerate(profiles):
        row_times = _row_times(profile, times)
        colors = _colors(len(row_times))
        x = _positions(len(profile.selections))
        values: list[np.ndarray] = []
        data_at: set[float] = set()
        for ax, column in zip(axes[row], columns, strict=True):
            for k, time in enumerate(row_times):
                band = _band_of(column, profile.label, time)
                color = colors[k] if column.color is None else column.color
                _fill_band(ax, x, band, color)
                text = f"Sim {_time_text(time, time_unit)}"
                ax.plot(
                    x,
                    band.median,
                    linestyle=_LINE_STYLES[k % len(_LINE_STYLES)],
                    color=color,
                    linewidth=2.0,
                    label=text if single else NO_LEGEND,
                )
                values += [
                    v for v in (band.median, band.lower, band.upper) if v is not None
                ]
            if data is not None:
                points = _draw_data(
                    ax,
                    _data_rows(data, profile, column),
                    row_times,
                    colors,
                    column,
                    time_unit,
                    labelled=single,
                )
                drawn += sum(v.size for v in points.values())
                values += list(points.values())
                data_at |= set(points)
        limits = _limits(values, _HEADROOM)
        for ax in axes[row]:
            if limits is not None:
                ax.set_ylim(*limits)
            ax.set_xlim(0.0, 1.0)
        axes[row, 0].set_ylabel(_label(profile.label, profile.unit))
        if not single:
            _row_legend(
                axes[row, -1],
                row_times,
                colors,
                by_time=by_time,
                data_at=sorted(data_at),
                time_unit=time_unit,
            )
    if data is not None and not drawn:
        raise ValueError(
            f"No row of the data is drawn: its times {sorted(set(data['time']))} "
            f"include none of the times {list(times)} of the figure."
        )
    for ax, column in zip(axes[0], columns, strict=True):
        if column.label:
            ax.set_title(column.label, fontweight="bold")
    if single:
        entries = [ax.get_legend_handles_labels() for ax in axes[:, 0]]
        if all(e[1] == entries[0][1] for e in entries):
            figure.legend(*entries[0], loc="outside right center")
        else:
            for ax in axes[:, 0]:
                ax.legend(
                    loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0
                )
    for ax in axes[-1]:
        ax.set_xlabel(POSITION)
    return figure


def _zones(result: pd.DataFrame, quantity: str, layer: Layer | str) -> int:
    """Get the number of zones of a layer in the result, 0 if it has none."""
    n = 0
    while observers.profile_ids(quantity, layer, n + 1)[-1] in result.columns:
        n += 1
    return n


def space_time(
    result: pd.DataFrame,
    quantity: str,
    *,
    layers: Sequence[Layer | str],
    outflow: bool = True,
    unit: str = DEFAULT_UNIT,
    time_unit: str = DEFAULT_TIME_UNIT,
    name: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Figure:
    """Draw the concentration over the zones and the time, a heatmap per layer.

    The heatmaps are as tall as the figure allows and the outflow panel below
    them is 0.6 of one. The axis of the outflow starts at zero when no value is
    negative, but a colour scale spans the range of the values it draws: the
    zero is explicit on the ticks of a colorbar, and a scale from zero would
    show a gradient of a few percent as one shade.

    Args:
        result: the simulation result.
        quantity: the id of the substance.
        layers: the layers, a heatmap each, from the top.
        outflow: add the concentration of the blood leaving the organ below.
        unit: the unit of the concentrations in the labels; the default is the
            unit of the profile observers (`observers.profile_ids`, mM), pass
            the unit of the result if it is another.
        time_unit: the unit of the time, the convention of the engine (min).
        name: the display name of the substance in the labels, the id by
            default.
        vmin: the lower end of the colour scale, the smallest value drawn by
            default; set it (with `vmax`) to clip a scale which a spike, e.g.
            the initial condition in one zone, dominates. It must be below
            `vmax`, so a `vmin` equal to the value of constant data raises,
            where the default scale of such data is the fallback zero to the
            value.
        vmax: the upper end of the colour scale, the largest value drawn by
            default.

    Returns:
        The figure.

    Raises:
        KeyError: if the result has no zone of a layer or no outflow.
        ValueError: if there is no layer and no outflow to draw.
    """
    if not layers and not outflow:
        raise ValueError(
            "The layers are empty and there is no outflow: nothing to draw."
        )
    _columns(result, ["time"])
    name = quantity if name is None else name
    time = result["time"].to_numpy(dtype=float)
    rows = len(layers) + int(outflow)
    heights = [HEATMAP_HEIGHT] * len(layers) + [OUTFLOW_HEIGHT] * int(outflow)
    figure = Figure(figsize=(7.0, 2.6 * sum(heights)), layout="constrained")
    axes = figure.subplots(
        rows, 1, sharex=True, squeeze=False, gridspec_kw={"height_ratios": heights}
    )[:, 0]
    for ax, layer in zip(axes[: len(layers)], layers, strict=True):
        n = _zones(result, quantity, layer)
        if n == 0:
            raise KeyError(f"The result has no zone of '{quantity}' in '{layer}'.")
        values = result[observers.profile_ids(quantity, layer, n)].to_numpy(dtype=float)
        limits = _scale([values], vmin, vmax) or (None, None)
        # `rasterized=False` is the matplotlib default and is stated for the
        # reader only: what keeps the SVG free of images is the colorbar
        # (`colorbar` un-rasterizes its solids, the default is an image) and
        # the test which looks for `<image` in the saved figures.
        mesh = ax.pcolormesh(
            time,
            _positions(n),
            values.T,
            shading="nearest",
            rasterized=False,
            vmin=limits[0],
            vmax=limits[1],
        )
        colorbar(figure, mesh, ax, _label(f"{name} {layer}", unit))
        ax.set_ylabel(POSITION.replace(" (", "\n("))
        ax.set_ylim(0.0, 1.0)
    if outflow:
        sid = observers.outflow_concentration(quantity)
        _columns(result, [sid])
        ax = axes[-1]
        y = result[sid].to_numpy(dtype=float)
        ax.plot(time, y, color="C0", linewidth=2.0)
        limits = _limits([y], _HEADROOM)
        if limits is not None:
            ax.set_ylim(*limits)
        ax.set_ylabel(_label(f"{name} outflow", unit))
    axes[-1].set_xlabel(f"Time [{time_unit}]")
    axes[-1].set_xlim(float(time.min()), float(time.max()))
    # no `align_ylabels`: constrained layout puts the labels of the column at one
    # x already (the wrapped label of a heatmap is wider, its right edge is the one)
    return figure


def layer_profile(
    result: pd.DataFrame,
    substance: str,
    *,
    time: float,
    unit: str = DEFAULT_UNIT,
    time_unit: str = DEFAULT_TIME_UNIT,
    share_y: bool = True,
) -> Figure:
    """Draw the profile over the zones in every layer which carries the substance.

    The layers are those of which the result holds the substance in zones
    from `z1` on, so a well-stirred build shows its vascular site and a build of
    the radial geometry its layers. The grid carries the label of the zones
    and the label of the concentration once each, under and left of it.

    The axes start at zero when no value is negative; a colour scale
    (`space_time`, `lobule_map`) spans the range of its values instead.

    Args:
        result: the simulation result.
        substance: the id of the substance.
        time: the time of the profile, interpolated in the result.
        unit: the unit of the concentration in the labels; the default is the
            unit of the profile observers (`observers.profile_ids`, mM).
        time_unit: the unit of the time, the convention of the engine (min).
        share_y: draw every panel on one concentration scale, from zero to the
            maximum over all layers (the default): the figure compares the
            layers at one time, and a layer of which the substance is nearly
            absent must look so. Pass `False` for a scale per panel, to see the
            detail of a layer of low concentration; then every panel keeps its
            own tick labels, as the numbers differ.

    Returns:
        The figure, a panel per layer.

    Raises:
        KeyError: if the result holds the substance in no layer.
    """
    present = [
        (layer, n)
        for layer in _LAYER_ORDER
        if (n := _zones(result, substance, layer)) > 0
    ]
    if not present:
        raise KeyError(f"The result has no zones of '{substance}' in any layer.")
    ncols = min(3, len(present))
    nrows = -(-len(present) // ncols)
    figure = Figure(figsize=(3.4 * ncols, 2.8 * nrows), layout="constrained")
    axes = figure.subplots(nrows, ncols, squeeze=False, sharey=share_y)
    flat = list(axes.flat)
    for ax in flat[len(present) :]:
        ax.remove()
    profiles = [
        _at(result, observers.profile_ids(substance, layer, n), time)
        for layer, n in present
    ]
    shared = _limits(profiles, _HEADROOM) if share_y else None
    for i, (ax, (layer, n), y) in enumerate(
        zip(flat[: len(present)], present, profiles, strict=True)
    ):
        ax.plot(_positions(n), y, color="C0", linewidth=2.0)
        limits = shared if share_y else _limits([y], _HEADROOM)
        if limits is not None:
            ax.set_ylim(*limits)
        ax.set_title(str(layer))
        ax.set_xlim(0.0, 1.0)
        if share_y:
            # one scale: its numbers on the left column (`sharey`) and the bottom row
            ax.tick_params(labelbottom=i + ncols >= len(present))
    figure.supxlabel(POSITION)
    figure.supylabel(_label(substance, unit))
    figure.suptitle(f"{substance} at {_time_text(time, time_unit)}")
    return figure


def _portal_triads(ax: Axes) -> None:
    """Draw the portal triad of every corner of the hexagon.

    Three dots around the corner, the portal vein, the hepatic artery and the
    bile duct, the convention of the figures of the machine perfusion model.
    The corners are on the angles 0, 60, ... degrees (`lobule_bands`).
    """
    offsets = _TRIAD_OFFSET * np.array([[1.0, 0.0], [-0.5, 0.87], [-0.5, -0.87]])
    angles = np.deg2rad(60.0 * np.arange(6))
    corners = np.column_stack([np.cos(angles), np.sin(angles)])
    for offset, color in zip(offsets, _TRIAD_COLORS, strict=True):
        x, y = (corners + offset).T
        # a dot lies a little beyond the limits, which hold the hexagon and not
        # the triads around its corners: clipped, it would be half a dot
        ax.plot(
            x,
            y,
            "o",
            ms=2.6,
            color=color,
            mec=_ANATOMY_GREY,
            mew=0.3,
            ls="",
            lw=0.0,
            clip_on=False,
        )


def draw_lobule(
    ax: Axes,
    values: Sequence[float] | np.ndarray,
    *,
    mappable: ScalarMappable,
    missing: str = MISSING,
    thirds: bool = False,
    anatomy: bool = False,
    labels: bool = False,
) -> None:
    """Draw a value per zone on the idealized hexagonal lobule, into given axes.

    The renderer of the lobule of this package, which `lobule_map` draws its
    one panel with and a model composes its own grid of maps from (its
    conditions, its measured rings, its legend stay in the model; a grid takes
    one colorbar over its panels through `colorbar`). The hexagon has the
    circumradius 1, the central vein at the origin and the portal triads at the
    corners on the angles 0, 60, ... degrees; a value is drawn as the band of
    its zone (`pkpdsim.zonation.render.lobule_bands`), twelve patches bounded
    by the isolines of the lobular position at `(k - 1) / n` and `k / n`, so
    the area of a band is the area share of the zone, `zone_areas(n)`. The
    centre of this idealization is the central vein itself, so no disc of a
    vein lumen covers the pericentral bands.

    The axes are set up for the lobule: the limits of the hexagon, an equal
    aspect, no ticks and no spines. Everything is a vector patch or line, no
    image.

    Args:
        ax: the axes, which the caller made and owns.
        values: the value of every zone, periportal first; their number is the
            number of zones. A value which is not finite is a zone without
            data and is drawn in `missing`.
        mappable: the colour of a value, a `ScalarMappable` of a norm and a
            colormap; the same one goes to `colorbar`.
        missing: the colour of a zone without data, a grey by default. It is
            taken here and not from the `bad` colour of the colormap, whose
            default is transparent: a caller who forgets to set it would get
            white holes in the hexagon, which read as a defect of the geometry
            and not as a zone without data.
        thirds: draw the isolines of 1/3 and 2/3 as dotted lines, the thirds of
            the lobulus the histology of a study counts its zones in.
        anatomy: draw the outline of the hexagon and the three dots of the
            portal triad of every corner (portal vein, hepatic artery, bile
            duct).
        labels: name the vessels, `CV` at the centre and `PT` below the lower
            right corner, for the first panel of a grid of several; both are
            written over a white stroke, so they read on a dark map and on a
            light one.

    Raises:
        ValueError: if there is no value, so no zone to draw.
    """
    per_zone = np.asarray(values, dtype=float)
    bands = lobule_bands(len(per_zone))
    # the mask keeps the norm off the values which are not finite; their colour
    # is `missing` and not what the colormap happens to give a masked value
    colours = mappable.to_rgba(np.ma.masked_invalid(per_zone))
    colours[~np.isfinite(per_zone)] = to_rgba(missing)
    for k, polygons in enumerate(bands, start=1):
        colour = colours[k - 1]
        for polygon in polygons:
            ax.add_patch(
                Polygon(
                    polygon,
                    closed=True,
                    facecolor=colour,
                    edgecolor=colour,
                    linewidth=0.3,
                    gid=ids.zone(k),
                )
            )
    if thirds:
        for position in _THIRDS:
            ax.plot(*lobule_isoline(position).T, ":", color=_THIRDS_GREY, lw=0.7)
    if anatomy:
        ax.plot(*lobule_isoline(0.0).T, color=_ANATOMY_GREY, lw=1.1)
        _portal_triads(ax)
    if labels:
        # the stroke carries the name over whatever the map is: the centre is
        # the dark end of the scale whenever the gradient falls to the vein
        halo = [withStroke(linewidth=2.0, foreground="white")]
        ax.text(
            0.0,
            0.0,
            "CV",
            ha="center",
            va="center",
            fontsize=_VESSEL_SIZE,
            path_effects=halo,
        )
        # below the corner at 300 degrees, under the panel and clear of the
        # dots of that triad: no band is covered and no dot is touched
        ax.text(
            0.5,
            -0.95,
            "PT",
            ha="center",
            va="top",
            fontsize=_VESSEL_SIZE,
            path_effects=halo,
        )
    ax.set_xlim(-1.02, 1.02)
    ax.set_ylim(-0.9, 0.9)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def lobule_map(
    result: pd.DataFrame,
    quantity: str,
    *,
    time: float,
    zones: int,
    layer: Layer | str = Layer.CELL,
    unit: str = DEFAULT_UNIT,
    time_unit: str = DEFAULT_TIME_UNIT,
    name: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Figure:
    """Draw the value of every zone on the idealized hexagonal lobule.

    The central vein is the centre of the hexagon, the six portal triads are
    its corners, and a point has the normalized periportal to pericentral
    distance `d_pp / (d_pp + d_pc)` of the nearest triad and the central
    vein: the rim and the corners are zone 1 (periportal), the centre is zone
    `zones` (pericentral). The bands of the zones are bounded by the isolines
    of that distance at `k / zones`, an interval of the position per zone, so
    their areas are the area shares of the zones, `zone_areas(zones)`. The
    figure is vector: a patch per band and wedge, no image, and the solids of
    its colorbar are not rasterized either.

    The panel is drawn by `draw_lobule`, with the outline of the hexagon and
    the portal triads of its corners; a model which composes its own grid of
    lobules calls that renderer and `colorbar` directly.

    Args:
        result: the simulation result.
        quantity: the id of the substance.
        time: the time of the map, interpolated in the result.
        zones: the number of zones of the result.
        layer: the layer whose zones are drawn, the hepatocytes by default.
        unit: the unit of the concentrations in the labels (the default is the
            unit of the profile observers, mM).
        time_unit: the unit of the time, the convention of the engine (min).
        name: the display name of the substance in the labels, the id by
            default.
        vmin: the lower end of the colour scale, the smallest value drawn by
            default; set it (with `vmax`) to clip a scale which a spike, e.g.
            the initial condition in one zone, dominates. It must be below
            `vmax`, so a `vmin` equal to the value of constant data raises,
            where the default scale of such data is the fallback zero to the
            value.
        vmax: the upper end of the colour scale, the largest value drawn by
            default.

    Returns:
        The figure, one panel of the hexagon and its colorbar. The colour
        scale spans the range of the values: a line axis starts at zero, a
        colour scale does not, because the zero is explicit on the ticks of a
        colorbar and a scale from zero would show a small gradient as one shade.

    Raises:
        KeyError: if the result has no column of a zone.
        ValueError: for less than one zone.
    """
    if zones < 1:
        raise ValueError(f"A lobule has at least one zone, not {zones}.")
    name = quantity if name is None else name
    values = _at(result, observers.profile_ids(quantity, layer, zones), time)
    limits = _scale([values], vmin, vmax) or (float(values.min()), float(values.max()))
    mappable = ScalarMappable(Normalize(*limits), colormaps["viridis"])
    figure = Figure(figsize=(5.0, 4.4), layout="constrained")
    ax = figure.subplots()
    draw_lobule(ax, values, mappable=mappable, anatomy=True)
    ax.set_xlabel(
        "Corners: portal triads (0, periportal)\nCentre: central vein (1, pericentral)"
    )
    ax.set_title(f"{name} {layer} at {_time_text(time, time_unit)}")
    colorbar(figure, mappable, ax, _label(f"{name} {layer}", unit))
    return figure
