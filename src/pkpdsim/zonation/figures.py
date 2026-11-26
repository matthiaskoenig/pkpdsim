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

A figure is a `sbmlsim.plot.Chart`, the spec `figure/1` which the renderer of
sbmlsim-core draws as SVG or PNG (`Chart.svg`, `Chart.png`, `Chart.save`), and
which is owned by the caller. Every figure is vector: the cells of a heatmap
and the bands of a lobule are polygons and the colour bar a gradient. A colour
of a value comes from `ColorScale`, which colours with the colour maps of
sbmlsim-core (`sbmlsim_core.color_scale`, imported when a colour is computed),
and `colorbar` draws the bar of the same scale, so a cell and the bar agree.
`draw_lobule` is the renderer of the lobule, which `lobule_map` draws its one
panel with and a model its own grid of lobules (`lobule_axes`,
`panel_aspects`).
"""

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd
from sbmlsim.engine import import_core
from sbmlsim.plot.chart import Chart, ChartAxis, ChartPlot
from sbmlsim.plot.colors import normalize_color

from pkpdsim.zonation import ids, observers
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.render import lobule_bands, lobule_isoline

__all__ = [
    "MISSING",
    "Band",
    "ColorScale",
    "ZoneColumn",
    "ZoneProfile",
    "colorbar",
    "data_axes",
    "draw_lobule",
    "fill_band",
    "layer_profile",
    "lobule_axes",
    "lobule_map",
    "panel_aspects",
    "result_frame",
    "space_time",
    "zone_profiles",
]

type _LineType = Literal["solid", "dashed", "dotted", "dashdot", "none"]

_LINE_STYLES: tuple[_LineType, ...] = ("dotted", "dashed", "solid", "dashdot")
"""The line style of every time, in the order of the times."""
_MARKERS: tuple[str, ...] = ("square", "diamond", "circle", "triangle")
"""The marker of the data of every time, in the order of the times."""
_LAYER_ORDER: tuple[str, ...] = (*(str(layer) for layer in Layer), ids.VASCULAR)
POSITION = "Position (0 periportal, 1 pericentral)"
"""The label of the axis of the zones, which names its ends: a flip must be visible."""
_POSITION_LINES = POSITION.replace(" (", "\n(")
"""The label of the axis of the zones in two lines, for a narrow panel or a y axis."""
DEFAULT_TIME_UNIT = "min"
"""The unit of the time of the engine (the rates of the engine are per minute)."""
HEATMAP_SPAN = 3
"""The rows of the grid of `space_time` a heatmap spans."""
OUTFLOW_SPAN = 2
"""The rows of the grid of `space_time` the outflow panel spans, about 0.65 of a heatmap.

A span of the grid is its rows and the space between them: with the narrow
space of `_SPACE_TIME_HSPACE` the outflow panel is a little more than two
thirds of three rows. The rows are few and tall on purpose, the renderer gives
every row of the grid a height of at least 20 pixels and a figure of many thin
rows overflows its bottom, where the label of the time is."""
_HEADROOM = 0.05
DEFAULT_UNIT = "mM"
"""The unit of the concentrations of the profile observers (`observers.profile_ids`)."""
_BAND_ALPHA = 0.25
"""The opacity of the fill of a band of the draws."""
_LINE_WIDTH = 2.0
"""The width of a line of a profile, in pixels."""
_DATA_SIZE = 6.0
"""The size of a marker of the data, in pixels."""
_SEGMENT_WIDTH = 6.0
"""The width of the segment which draws the band of a single zone, in pixels."""
_ONE_ZONE_MARKER = "circle"
"""The marker of a profile of a single zone.

One zone is the well-stirred limit of the engine and must draw as easily as a
sinusoid of ten: a line of one point has no segment to draw, so the panel is
empty without a marker."""
_NEUTRAL = "#4d4d4d"
"""The colour of a legend entry which stands for a time in every column (grey 0.3)."""
_EDGE = "#000000"
"""The colour of the edge of a marker of the data."""
_TRIAD_OFFSET = 0.045
"""The distance of a dot of a portal triad from its corner, of the circumradius."""
_TRIAD_COLORS: tuple[str, ...] = ("#1f77b4", "#d62728", "#2ca02c")
"""The colours of the portal vein, the hepatic artery and the bile duct of a triad."""
_TRIAD_SIZE = 4.0
"""The size of a dot of a portal triad, in pixels."""
_THIRDS: tuple[float, ...] = (1.0 / 3.0, 2.0 / 3.0)
"""The lobular positions which bound the thirds a histology counts its zones in."""
_ANATOMY_GREY = "#262626"
"""The grey of the outline of the hexagon and of the edge of a dot of a triad (0.15)."""
_THIRDS_GREY = "#404040"
"""The grey of the dotted isolines of the thirds (0.25)."""
_SEAM_WIDTH = 0.6
"""The width of the outline of a cell or a band in its own colour, in pixels.

Two polygons which share an edge leave a seam of the background where their
antialiased edges meet; an outline in the colour of the fill closes it."""
MISSING = "#f0f0f0"
"""The colour of a zone without data on a lobule, the default of `draw_lobule`."""
LOBULE_X = 1.1
"""The half width the axes of a lobule hold: the hexagon and the dots of its triads."""
LOBULE_Y = 1.08
"""The half height the axes of a lobule with the names hold: the hexagon, its triads and `PT` below."""
LOBULE_Y_BARE = 0.93
"""The half height the axes of a lobule without names hold: the hexagon and its triads, which reach 0.91."""
_CV = (0.0, 0.0)
"""The position of the name of the central vein, the centre of the hexagon."""
_PT = (0.5, -1.03)
"""The position of the name of the portal triad, below the triad dots of the corner at 300 degrees."""
_LABEL_MARGIN = 0.04
"""The room a name needs between its position and the limit of its axis."""
_LIGHT = 0.5
"""The luminance above which a colour is light and a text on it is dark."""
_PX_PER_INCH = 72.0
"""The pixels of a chart per inch of the figures this module drew with matplotlib."""
_ROW_HEIGHT = 3.0 * _PX_PER_INCH
"""The height of a row of `zone_profiles`."""
_COLUMN_WIDTH = 3.6 * _PX_PER_INCH
"""The width of a column of a grid of `zone_profiles`, without its legend."""
_SINGLE_WIDTH = 5.0 * _PX_PER_INCH
"""The width of `zone_profiles` of a frame, without its legend."""
_SPACE_TIME_WIDTH = 7.0 * _PX_PER_INCH
"""The width of `space_time`, its colour bars included."""
_HEATMAP_HEIGHT = 3.2 * _PX_PER_INCH
"""The height of a heatmap of `space_time` with its share of the margins.

Tall enough for its panel to hold the label of its position axis, which runs
along it in two lines; a lower panel lets the label stick out of the figure."""
_SPACE_TIME_HSPACE = 0.1
"""The space between the rows of the grid of `space_time`, of a row.

Small, so a span is its rows; the space between two panels grows to what
their axes need."""
_LAYER_WIDTH = 3.4 * _PX_PER_INCH
"""The width of a column of `layer_profile`."""
_LAYER_HEIGHT = 2.8 * _PX_PER_INCH
"""The height of a row of `layer_profile`."""
_LOBULE_SIZE = (5.0 * _PX_PER_INCH, 4.4 * _PX_PER_INCH)
"""The size of `lobule_map`, its colour bar and its caption included."""


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
            the time; any colour `sbmlsim.plot.colors.normalize_color` reads.
    """

    label: str
    values: Mapping[tuple[str, float], Band]
    color: str | None = None


@dataclass(frozen=True)
class ColorScale:
    """A colour scale: the values from `low` to `high` on a colour map of sbmlsim-core.

    The one way a figure of this package or of a model colours a value: the
    colours of the cells or the bands come from `colors`, and `colorbar`
    draws the bar of the same scale, so a value has the colour the bar shows
    for it. The colour maps and the part of a map are those of
    `sbmlsim_core.color_scale`.

    Attributes:
        low: the value of the lower end of the scale.
        high: the value of the upper end of the scale, above `low`.
        name: the colour map, `viridis`, `magma`, `rdbu`, `coolwarm` or `greys`.
        log: whether the scale is logarithmic; `low` must then be positive.
        start: the part of the colour map `low` takes, from 0 to 1.
        end: the part of the colour map `high` takes, above `start`.

    Raises:
        ValueError: for ends which are not finite or not `low < high`, a log
            scale whose `low` is not positive and a part which is not
            `0 <= start < end <= 1`.
    """

    low: float
    high: float
    name: str = "viridis"
    log: bool = False
    start: float = 0.0
    end: float = 1.0

    def __post_init__(self) -> None:
        """Check the ends and the part of the colour map."""
        if not (math.isfinite(self.low) and math.isfinite(self.high)):
            raise ValueError(
                f"The colour scale {self.low} to {self.high} has an end which is "
                "not finite."
            )
        if not self.low < self.high:
            raise ValueError(
                f"The colour scale {self.low} to {self.high} has no range."
            )
        if self.log and self.low <= 0.0:
            raise ValueError(
                f"The log colour scale {self.low} to {self.high} needs a positive "
                "lower end."
            )
        if not 0.0 <= self.start < self.end <= 1.0:
            raise ValueError(
                f"The part {self.start} to {self.end} of the colour map is not "
                "0 <= start < end <= 1."
            )

    def colors(self, values: Sequence[float] | np.ndarray) -> list[str | None]:
        """Get the colour of every value as `#rrggbb`.

        A value beyond an end takes the colour of that end. A value which is
        not finite, and on a log scale one which is not positive, has no colour
        (`None`): it is a gap, which a figure leaves out or draws in a colour of
        its own.

        Args:
            values: the values, flattened.

        Returns:
            The colour of every value, `None` for a gap.
        """
        flat = np.asarray(values, dtype=float).ravel()
        drawn = np.isfinite(flat)
        if self.log:
            drawn &= flat > 0.0
        colours: list[str | None] = [None] * flat.size
        if drawn.any():
            computed = import_core().color_scale(
                self.name,
                float(self.low),
                float(self.high),
                [float(v) for v in flat[drawn]],
                log=self.log,
                start=float(self.start),
                end=float(self.end),
            )
            for index, colour in zip(np.flatnonzero(drawn), computed, strict=True):
                colours[int(index)] = str(colour)
        return colours


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


def data_axes(figure: Chart) -> list[ChartPlot]:
    """Get the plots of the panels of a figure.

    A caller which post-processes a returned figure (an element, a colour
    bar) adds to its plots: a colour bar is a field of its plot and no plot
    of its own, so every plot of a chart is a panel.

    Args:
        figure: the chart of a function of this module.

    Returns:
        The plots, in the order they were added.
    """
    return list(figure.plots)


def colorbar(plot: ChartPlot, scale: ColorScale, label: str | None = None) -> None:
    """Add the colour bar of a scale to a plot.

    The one way a figure of this package or of a model gets a colour bar: it
    shows exactly the scale `ColorScale.colors` colours the values with, the
    same map, ends, part of the map and log scale. Several bars of one plot
    stack right of its panel.

    Args:
        plot: the plot the bar is drawn right of.
        scale: the scale which is drawn.
        label: the label of the bar, `None` none.
    """
    plot.colorbar(
        scale.name,
        scale.low,
        scale.high,
        label=label,
        log=scale.log,
        start=scale.start,
        end=scale.end,
    )


def _color(color: str) -> str:
    """Get a colour as `#rrggbb`, through `sbmlsim.plot.colors.normalize_color`.

    Raises:
        ValueError: for a colour which is no colour or is `none`, which draws
            nothing.
    """
    normalized = normalize_color(color)
    if normalized is None:
        raise ValueError(
            f"The colour '{color}' draws nothing, a figure needs a colour."
        )
    return normalized


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
    """Get the label of an axis with its unit, the unit as the caller wrote it."""
    return label if unit in {"", "dimensionless"} else f"{label} [{unit}]"


def _finite(values: Sequence[np.ndarray]) -> np.ndarray:
    """Get the finite values of a sequence of arrays, flat.

    A single NaN makes `np.min` and `np.max` NaN and every comparison against
    it false, so a scale computed from them falls through to its no-range
    branch: data from 200 to 1000 with one gap was drawn from 0 to 1, every
    zone saturated and the colorbar stated a range the figure did not have. A
    gap of the data is not a reason to draw a wrong figure, and no data at all
    is not a reason to draw a plausible one.

    Args:
        values: the arrays of the axis or the colour scale.

    Returns:
        The finite values of all of them.

    Raises:
        ValueError: if no value is finite, which is no scale (and if there is
            no value at all).
    """
    flat = [np.asarray(v, dtype=float).ravel() for v in values]
    finite = np.concatenate(flat) if flat else np.empty(0)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        total = sum(int(v.size) for v in flat)
        raise ValueError(
            f"None of the {total} values is finite (NaN or infinite), so there "
            "is no scale to draw them on: the figure would state a range it "
            "does not have."
        )
    return finite


def _limits(
    values: Sequence[np.ndarray], headroom: float = 0.0
) -> tuple[float, float] | None:
    """Get the limits of a concentration axis, from zero if nothing is negative.

    A concentration is never negative, so an axis of values at or above zero
    starts at zero; a layer which is zero everywhere gets zero to one. Values
    below zero leave the limits to the renderer (`None`). A colour scale does
    not start at zero, it spans its values (`_scale`). NaN is ignored, and
    nothing but NaN raises (`_finite`).

    Raises:
        ValueError: if no value is finite.
    """
    finite = _finite(values)
    low, high = float(finite.min()), float(finite.max())
    if low < 0.0:
        return None
    return (0.0, high * (1.0 + headroom) if high > 0.0 else 1.0)


def _scale(
    values: Sequence[np.ndarray],
    vmin: float | None = None,
    vmax: float | None = None,
    *,
    log: bool = False,
) -> tuple[float, float]:
    """Get the limits of a colour scale: the range of the values it draws.

    A colour scale spans the range of its values, a line axis keeps its zero
    baseline (`_limits`): the zero is explicit on the ticks of a colour bar,
    and a scale from zero would draw a gradient of a few percent as one shade.
    A scale of no range is never degenerate, the renderer needs `low < high`
    and a bar of one value states nothing: a constant at or above zero takes
    zero to the value (zero to one for zeros), a negative constant the value
    to zero, and on a log scale a constant the decade below it to the value.
    An explicit `vmin` or `vmax` replaces that end of the range, the way to
    clip a scale which a spike dominates; a `vmin` equal to the value of a
    constant therefore raises (no range), while the default call for the same
    data returns the fallback.

    A log scale spans the positive values: a zero or a negative value has no
    place on it and is a gap of the figure, not an end of its scale. NaN is
    ignored, and nothing but NaN raises (`_finite`).

    Raises:
        ValueError: if the limits are not `vmin < vmax`, if no value is finite,
            on a log scale if no value is positive or an end is not positive.
    """
    finite = _finite(values)
    if log:
        finite = finite[finite > 0.0]
        if finite.size == 0:
            raise ValueError(
                "No value is positive, so a log colour scale has nothing to span."
            )
    low, high = float(finite.min()), float(finite.max())
    if vmin is not None or vmax is not None:
        low = low if vmin is None else vmin
        high = high if vmax is None else vmax
        if not low < high:
            raise ValueError(f"The colour scale {low} to {high} has no range.")
        if log and low <= 0.0:
            raise ValueError(
                f"The log colour scale {low} to {high} needs a positive lower end."
            )
        return (low, high)
    if high > low:
        return (low, high)
    if log:
        return (high / 10.0, high)
    if high < 0.0:
        return (high, 0.0)
    return (0.0, high if high > 0.0 else 1.0)


def _time_text(time: float, time_unit: str) -> str:
    """Get a time with its unit, `t = 0.5 min`."""
    return f"t = {time:g} {time_unit}"


def _time_colors(count: int) -> list[str]:
    """Get the colour of every time, a sequential map from early to late.

    The part 0.1 to 0.75 of viridis, which leaves out its darkest and its
    lightest colours, so every line reads on white.
    """
    scale = ColorScale(0.0, 1.0, start=0.1, end=0.75)
    return [str(c) for c in scale.colors(np.linspace(0.0, 1.0, count))]


def _row_times(profile: ZoneProfile, times: Sequence[float]) -> list[float]:
    """Get the times of a row: its own times, else the times of the figure.

    Raises:
        ValueError: for a profile whose own times are empty, which draws
            nothing; the figure says so instead of failing on a bare
            `min() iterable argument is empty` when the row has no values to
            scale.
    """
    own = profile.times is not None
    row = [float(t) for t in (profile.times if own else times)]
    if not row:
        raise ValueError(
            f"The times of the profile '{profile.label}' are empty, so its row "
            "draws no line; leave them out to take the times of the figure."
        )
    return row


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


def fill_band(
    plot: ChartPlot,
    x: np.ndarray,
    band: Band,
    color: str,
    *,
    alpha: float = _BAND_ALPHA,
) -> None:
    """Fill a band between its ends; a band without ends fills nothing.

    The look of the bands of `zone_profiles`, for a model which draws a profile
    of its own beside a map and fills the band of a run the same way. The fill
    is in no legend, the line of the median names it. A gap of an end (NaN) is
    a gap of the fill.

    Args:
        plot: the plot to fill into.
        x: the positions of the values of the band.
        band: the band, of which the lower and the upper end are filled.
        color: the colour of the fill, any colour
            `sbmlsim.plot.colors.normalize_color` reads.
        alpha: the opacity of the fill.
    """
    if band.lower is None or band.upper is None:
        return
    colour = _color(color)
    if x.size == 1:
        # a fill between one x is empty, so the band of a single zone is drawn
        # as the segment it is; without it the panel would show the median
        # alone and state no uncertainty
        plot.curve(
            np.repeat(x, 2),
            [float(band.lower[0]), float(band.upper[0])],
            width=_SEGMENT_WIDTH,
            color=colour,
            opacity=alpha,
        )
        return
    plot.area(x, band.lower, band.upper, color=colour, fill_opacity=alpha)


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


def _data_at(rows: pd.DataFrame, times: Sequence[float]) -> dict[int, pd.DataFrame]:
    """Get the rows of the data of every time of a row, by the index of the time.

    A time of which no row is drawn is absent; the rows of a time which the
    row does not draw are skipped, as plotting a subset of the times of a
    dataset is legitimate.
    """
    at: dict[int, pd.DataFrame] = {}
    for k, time in enumerate(times):
        matched = rows[np.isclose(rows["time"], time)]
        if not matched.empty:
            at[k] = matched
    return at


def _sim_text(time: float, time_unit: str) -> str:
    """Get the legend entry of the line of a time."""
    return f"Sim {_time_text(time, time_unit)}"


def _data_text(time: float, time_unit: str) -> str:
    """Get the legend entry of the data of a time."""
    return f"Data {_time_text(time, time_unit)}"


def _swatches(
    times: Sequence[float], data_at: Sequence[int], time_unit: str
) -> dict[str, dict[str, Any]]:
    """Get the neutral swatches of the legend of a grid whose columns carry the colour.

    The columns are named by the titles of the top row, so an entry stands
    for a time in every column: its line style or its marker, in a neutral
    grey and not in the colour of one column.

    Args:
        times: the times of the row.
        data_at: the indices of the times of which the plot draws data.
        time_unit: the unit of the times.

    Returns:
        The style of every legend group of the plot.
    """
    swatches: dict[str, dict[str, Any]] = {}
    for k, time in enumerate(times):
        swatches[_sim_text(time, time_unit)] = {
            "color": _NEUTRAL,
            "line": _LINE_STYLES[k % len(_LINE_STYLES)],
            "width": _LINE_WIDTH,
        }
        if k in data_at:
            swatches[_data_text(time, time_unit)] = {
                "color": _NEUTRAL,
                "line": "none",
                "marker": _MARKERS[k % len(_MARKERS)],
                "size": _DATA_SIZE,
                "marker_edge": _EDGE,
            }
    return swatches


def _legend_order(times: Sequence[float], time_unit: str) -> list[str]:
    """Get the order of the legend of a row: every time, its line and its data."""
    order: list[str] = []
    for time in times:
        order += [_sim_text(time, time_unit), _data_text(time, time_unit)]
    return order


def _shared_axis(
    *, label: str | None, outer: bool, limits: tuple[float, float] | None
) -> ChartAxis:
    """Get an axis shared by the panels of a row or a column.

    The outer panel carries the label and the numbers of the ticks, an inner
    one the marks of the ticks alone.
    """
    return ChartAxis(
        label=label if outer else None,
        label_visible=outer,
        ticks_visible=outer,
        tick_marks_visible=True,
        min=None if limits is None else limits[0],
        max=None if limits is None else limits[1],
    )


def _draw_data(
    plot: ChartPlot,
    at: Mapping[int, pd.DataFrame],
    times: Sequence[float],
    colors: Sequence[str],
    column: ZoneColumn,
    time_unit: str,
) -> list[np.ndarray]:
    """Draw the data of a row in a column, the markers of a time.

    A value which is NaN is a gap: its marker is left out, the others of its
    time are drawn.

    Args:
        plot: the panel.
        at: the rows of the data of every time drawn, by the index of the time.
        times: the times of the row.
        colors: the colour of every time.
        column: the column, whose colour the markers take where it has one.
        time_unit: the unit of the times.

    Returns:
        The drawn values plus their error, which the axis must hold.
    """
    drawn: list[np.ndarray] = []
    for k, rows in at.items():
        time = times[k]
        value = rows["value"].to_numpy(dtype=float)
        sd = rows["sd"].to_numpy(dtype=float) if "sd" in rows else None
        drawn.append(value if sd is None else value + np.nan_to_num(sd))
        text = _data_text(time, time_unit)
        plot.curve(
            rows["position"].to_numpy(dtype=float),
            value,
            yerr=sd,
            line="none",
            marker=_MARKERS[k % len(_MARKERS)],
            size=_DATA_SIZE,
            marker_edge=_EDGE,
            color=colors[k] if column.color is None else _color(column.color),
            label=text,
            legend=text,
        )
    return drawn


def zone_profiles(
    result: pd.DataFrame | Sequence[ZoneColumn],
    profiles: Sequence[ZoneProfile],
    *,
    times: Sequence[float],
    data: pd.DataFrame | None = None,
    time_unit: str = DEFAULT_TIME_UNIT,
) -> Chart:
    """Draw a row per profile and a column per condition, a line per time.

    A `pandas.DataFrame` is one condition: one column, no title, a colour and
    a line style per time (a sequential map, early to late) and the data of a
    time in the colour of that time. A sequence of `ZoneColumn` is a grid of
    conditions: the label of a column is the title of its panel in the top
    row, the colour of a line is the colour of its column where it has one and
    the line style stays the style of the time, so a reader tells the times
    apart inside a column and the columns apart by colour. A band of the draws
    of a fit is filled between its ends. Every row has one legend right of its
    last column, which names the times of that row and their data: in the
    colour of the time where the colours are the times (no column has a
    colour), in a neutral grey where a column has one.
    `pkpdsim.zonation.results.zone_column` builds a column from a result of
    sbmlsim.

    The x axis is shared and the y axis is shared per row, over the values of
    every column of that row: the outer panels carry the labels and the
    numbers. A line axis starts at zero when no value is negative; a colour
    scale (`space_time`, `lobule_map`) spans the range of its values instead.
    A value which is NaN, of a band (a layer a condition does not have) or of
    the data, is a gap of its line or its markers.

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
        The chart, a plot per row and column in the order of the rows.

    Raises:
        KeyError: if the result has no column of a selection, or if a column
            has no band of a profile and a time it is asked for. A frame is
            read before the data is checked, so a missing selection is named
            before a label of the data which no profile has.
        ValueError: if `profiles`, `times`, `result` or the own times of a
            profile are empty, if a row of `data` has a label which no profile
            or no column has, if no row of `data` is drawn, or if no value of a
            row is finite.
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
    chart = Chart(
        _SINGLE_WIDTH if single else _COLUMN_WIDTH * n_columns,
        _ROW_HEIGHT * n_rows,
        rows=n_rows,
        columns=n_columns,
        legend_scope="row",
        legend_position="outside",
    )
    drawn = 0
    for row, profile in enumerate(profiles):
        row_times = _row_times(profile, times)
        colors = _time_colors(len(row_times))
        x = _positions(len(profile.selections))
        bands = [
            [_band_of(column, profile.label, time) for time in row_times]
            for column in columns
        ]
        at = [
            {}
            if data is None
            else _data_at(_data_rows(data, profile, column), row_times)
            for column in columns
        ]
        values = [
            v
            for column_bands in bands
            for band in column_bands
            for v in (band.median, band.lower, band.upper)
            if v is not None
        ]
        values += [
            rows["value"].to_numpy(dtype=float)
            + (np.nan_to_num(rows["sd"].to_numpy(dtype=float)) if "sd" in rows else 0.0)
            for column_at in at
            for rows in column_at.values()
        ]
        drawn += sum(len(rows) for column_at in at for rows in column_at.values())
        limits = _limits(values, _HEADROOM)
        for col, column in enumerate(columns):
            plot = chart.plot(
                row,
                col,
                title=column.label if row == 0 and column.label else None,
                x=_shared_axis(
                    label=POSITION, outer=row == n_rows - 1, limits=(0.0, 1.0)
                ),
                y=_shared_axis(
                    label=_label(profile.label, profile.unit),
                    outer=col == 0,
                    limits=limits,
                ),
                legend_order=_legend_order(row_times, time_unit),
                legend_swatches=None
                if by_time
                else _swatches(row_times, list(at[col]), time_unit),
            )
            for k, (time, band) in enumerate(zip(row_times, bands[col], strict=True)):
                color = colors[k] if column.color is None else _color(column.color)
                fill_band(plot, x, band, color)
                text = _sim_text(time, time_unit)
                plot.curve(
                    x,
                    band.median,
                    line=_LINE_STYLES[k % len(_LINE_STYLES)],
                    width=_LINE_WIDTH,
                    marker=_ONE_ZONE_MARKER if x.size == 1 else "none",
                    color=color,
                    label=text,
                    legend=text,
                )
            _draw_data(plot, at[col], row_times, colors, column, time_unit)
    if data is not None and not drawn:
        raise ValueError(
            f"No row of the data is drawn: its times {sorted(set(data['time']))} "
            f"include none of the times {list(times)} of the figure."
        )
    return chart


def _zones(result: pd.DataFrame, quantity: str, layer: Layer | str) -> int:
    """Get the number of zones of a layer in the result, 0 if it has none."""
    n = 0
    while observers.profile_ids(quantity, layer, n + 1)[-1] in result.columns:
        n += 1
    return n


def _edges(centres: np.ndarray) -> np.ndarray:
    """Get the edges of the cells around their centres, the ends at the first and last.

    The edge between two centres is their midpoint, the cell of a centre
    reaches half way to its neighbours, and the outer cells end at their
    centre, the ends of the axis of the time.
    """
    mids = 0.5 * (centres[1:] + centres[:-1])
    return np.concatenate([centres[:1], mids, centres[-1:]])


def _draw_cells(
    plot: ChartPlot,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    values: np.ndarray,
    scale: ColorScale,
) -> None:
    """Draw the cells of a heatmap as polygons on numeric axes.

    The cell `(i, j)` spans `x_edges[i]..x_edges[i + 1]` and
    `y_edges[j]..y_edges[j + 1]` and has the colour of `values[i, j]`; the
    cells of one colour are one polygon element of a ring each, so a heatmap
    of thousands of cells is a few hundred elements. A cell without a colour
    (NaN, a value which is not positive on a log scale) and a cell of no width
    (a time twice, at a change) are left out: a gap shows the background.

    Args:
        plot: the plot of the heatmap.
        x_edges: the edges of the cells along x, one more than the rows of
            `values`.
        y_edges: the edges along y, one more than the columns of `values`.
        values: the values, `(x, y)`.
        scale: the colour scale.
    """
    colours = scale.colors(values)
    rings: dict[str, list[np.ndarray]] = {}
    n_y = values.shape[1]
    for index, colour in enumerate(colours):
        if colour is None:
            continue
        i, j = divmod(index, n_y)
        x0, x1 = x_edges[i], x_edges[i + 1]
        if not x1 > x0:
            continue
        y0, y1 = y_edges[j], y_edges[j + 1]
        rings.setdefault(colour, []).append(
            np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
        )
    for colour, cells in rings.items():
        plot.polygon(cells, color=colour, width=_SEAM_WIDTH, fill_opacity=1.0)


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
    log: bool = False,
) -> Chart:
    """Draw the concentration over the zones and the time, a heatmap per layer.

    A heatmap is the cells of the zones and the time points on the numeric
    axes of the time and the position, which the outflow panel below them
    shares: the cell of a time point reaches half way to its neighbours, as
    the time points of a simulation are not equidistant. The outflow panel is
    about two thirds of a heatmap. The axis of the outflow starts at zero when
    no value is negative, but a colour scale spans the range of the values it
    draws: the
    zero is explicit on the ticks of a colour bar, and a scale from zero would
    show a gradient of a few percent as one shade. A value which is NaN, and
    on a log scale a value which is not positive, is a gap of the heatmap.

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
        log: a logarithmic colour scale, over the positive values: a scale of
            values over decades, of which a linear scale shows the largest
            only. A zero, e.g. the zones the substance has not reached yet, is
            a gap.

    Returns:
        The chart, a plot per layer from the top and the outflow last.

    Raises:
        KeyError: if the result has no zone of a layer or no outflow.
        ValueError: if there is no layer and no outflow to draw, if the time
            has fewer than two distinct points (the cells of the time would
            have no width), a non-finite point or is not increasing, if a
            layer has no finite value, and on a log scale if a layer has no
            positive value.
    """
    if not layers and not outflow:
        raise ValueError(
            "The layers are empty and there is no outflow: nothing to draw."
        )
    _columns(result, ["time"])
    name = quantity if name is None else name
    time = result["time"].to_numpy(dtype=float)
    if not np.isfinite(time).all():
        raise ValueError("space_time needs finite time points.")
    if np.unique(time).size < 2:
        raise ValueError("space_time needs at least two time points.")
    if (np.diff(time) < 0.0).any():
        raise ValueError("space_time needs the time points in increasing order.")
    rows = HEATMAP_SPAN * len(layers) + OUTFLOW_SPAN * int(outflow)
    height = _HEATMAP_HEIGHT * rows / HEATMAP_SPAN
    chart = Chart(
        _SPACE_TIME_WIDTH, height, rows=rows, columns=1, hspace=_SPACE_TIME_HSPACE
    )
    limits_time = (float(time.min()), float(time.max()))
    x_edges = _edges(time)
    for index, layer in enumerate(layers):
        n = _zones(result, quantity, layer)
        if n == 0:
            raise KeyError(f"The result has no zone of '{quantity}' in '{layer}'.")
        values = result[observers.profile_ids(quantity, layer, n)].to_numpy(dtype=float)
        scale = ColorScale(*_scale([values], vmin, vmax, log=log), log=log)
        last = index == len(layers) - 1 and not outflow
        plot = chart.plot(
            HEATMAP_SPAN * index,
            0,
            row_span=HEATMAP_SPAN,
            x=_shared_axis(label=f"Time [{time_unit}]", outer=last, limits=limits_time),
            y=ChartAxis(label=_POSITION_LINES, min=0.0, max=1.0),
            legend=False,
        )
        _draw_cells(plot, x_edges, np.linspace(0.0, 1.0, n + 1), values, scale)
        colorbar(plot, scale, _label(f"{name} {layer}", unit))
    if outflow:
        sid = observers.outflow_concentration(quantity)
        _columns(result, [sid])
        y = result[sid].to_numpy(dtype=float)
        plot = chart.plot(
            HEATMAP_SPAN * len(layers),
            0,
            row_span=OUTFLOW_SPAN,
            x=ChartAxis(
                label=f"Time [{time_unit}]", min=limits_time[0], max=limits_time[1]
            ),
            y=_shared_axis(
                label=_label(f"{name} outflow", unit),
                outer=True,
                limits=_limits([y], _HEADROOM),
            ),
            legend=False,
        )
        plot.curve(time, y, width=_LINE_WIDTH, color=_color("C0"))
    return chart


def layer_profile(
    result: pd.DataFrame,
    substance: str,
    *,
    time: float,
    unit: str = DEFAULT_UNIT,
    time_unit: str = DEFAULT_TIME_UNIT,
    share_y: bool = True,
) -> Chart:
    """Draw the profile over the zones in every layer which carries the substance.

    The layers are those of which the result holds the substance in zones
    from `z1` on, so a well-stirred build shows its vascular site and a build of
    the radial geometry its layers. The label of the zones is under the last
    panel of every column and the label of the concentration left of the first
    panel of every row; the title names the substance and the time.

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
            absent must look so; the numbers of the scale are on the first
            panel of every row. Pass `False` for a scale per panel, to see the
            detail of a layer of low concentration; then every panel keeps its
            own numbers, as they differ.

    Returns:
        The chart, a plot per layer.

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
    chart = Chart(
        _LAYER_WIDTH * ncols,
        _LAYER_HEIGHT * nrows,
        title=f"{substance} at {_time_text(time, time_unit)}",
        rows=nrows,
        columns=ncols,
    )
    profiles = [
        _at(result, observers.profile_ids(substance, layer, n), time)
        for layer, n in present
    ]
    shared = _limits(profiles, _HEADROOM) if share_y else None
    label = _label(substance, unit)
    for i, ((layer, n), y) in enumerate(zip(present, profiles, strict=True)):
        row, col = divmod(i, ncols)
        bottom = i + ncols >= len(present)
        limits = shared if share_y else _limits([y], _HEADROOM)
        plot = chart.plot(
            row,
            col,
            title=str(layer),
            x=ChartAxis(
                label=_POSITION_LINES if bottom else None,
                label_visible=bottom,
                ticks_visible=bottom or not share_y,
                tick_marks_visible=True,
                min=0.0,
                max=1.0,
            ),
            y=ChartAxis(
                label=label if col == 0 else None,
                label_visible=col == 0,
                ticks_visible=col == 0 or not share_y,
                tick_marks_visible=True,
                min=None if limits is None else limits[0],
                max=None if limits is None else limits[1],
            ),
            legend=False,
        )
        plot.curve(
            _positions(n),
            y,
            width=_LINE_WIDTH,
            marker=_ONE_ZONE_MARKER if n == 1 else "none",
            color=_color("C0"),
        )
    return chart


def lobule_axes(
    aspect: float | None = None, *, label: str | None = None, labels: bool = False
) -> tuple[ChartAxis, ChartAxis]:
    """Get the axes of a panel of the lobule: the hexagon at one scale, no ticks.

    The renderer has no equal aspect of its own, so the limits do it: they
    hold the hexagon with its triads (`LOBULE_X`, and `LOBULE_Y_BARE` high, or
    `LOBULE_Y` with the room of the names below it) and widen the
    direction in which the panel has room to spare, so a unit of x and a
    unit of y are one length and the hexagon is regular. `panel_aspects`
    measures the aspect of the panels of a chart.

    Args:
        aspect: the width over the height of the panel in pixels, `None` the
            aspect of the hexagon itself.
        label: the label under the panel, `None` none.
        labels: whether the panel carries the names of the vessels
            (`draw_lobule(labels=True)`), which need room below the hexagon.

    Returns:
        The x and the y axis.

    Raises:
        ValueError: for an aspect which is not positive.
    """
    height = LOBULE_Y if labels else LOBULE_Y_BARE
    own = LOBULE_X / height
    aspect = own if aspect is None else float(aspect)
    if not aspect > 0.0:
        raise ValueError(f"The aspect of a panel is positive, not {aspect}.")
    half_x, half_y = (
        (LOBULE_X, LOBULE_X / aspect)
        if aspect < own
        else (
            height * aspect,
            height,
        )
    )
    hidden: dict[str, Any] = {"ticks_visible": False, "tick_marks_visible": False}
    return (
        ChartAxis(
            label=label,
            label_visible=label is not None,
            min=-half_x,
            max=half_x,
            **hidden,
        ),
        ChartAxis(label_visible=False, min=-half_y, max=half_y, **hidden),
    )


def panel_aspects(chart: Chart) -> dict[str, float]:
    """Get the width over the height of the panel of every plot, as the renderer lays them out.

    The aspect a panel of the lobule needs for `lobule_axes`: a caller builds
    the chart once with any aspect, measures it and builds it again with the
    measured ones. The layout does not depend on the limits of axes without
    ticks, so the second chart has the panels the first one had.

    Args:
        chart: the chart.

    Returns:
        The aspect of the panel of every plot, by the id of the plot.
    """
    geometry = import_core().figure_geometry(chart.to_json())
    return {
        str(plot["id"]): float(plot["rect"]["width"]) / float(plot["rect"]["height"])
        for plot in geometry["plots"]
    }


def _portal_triads(plot: ChartPlot) -> None:
    """Draw the portal triad of every corner of the hexagon.

    Three dots around the corner, the portal vein, the hepatic artery and the
    bile duct, the convention of the figures of the machine perfusion model.
    The corners are on the angles 0, 60, ... degrees (`lobule_bands`).
    """
    # the three dots stand at the same distance from the corner, 120 degrees apart
    around = np.deg2rad(120.0 * np.arange(3))
    offsets = _TRIAD_OFFSET * np.column_stack([np.cos(around), np.sin(around)])
    angles = np.deg2rad(60.0 * np.arange(6))
    corners = np.column_stack([np.cos(angles), np.sin(angles)])
    for offset, color in zip(offsets, _TRIAD_COLORS, strict=True):
        x, y = (corners + offset).T
        plot.curve(
            x,
            y,
            line="none",
            marker="circle",
            size=_TRIAD_SIZE,
            marker_edge=_ANATOMY_GREY,
            marker_edge_width=0.3,
            color=color,
        )


def _luminance(colour: str) -> float:
    """Get the relative luminance of a `#rrggbb` colour, 0 black and 1 white."""
    r, g, b = (int(colour[i : i + 2], 16) / 255.0 for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _check_label_room(plot: ChartPlot) -> None:
    """Raise if the limits of a plot do not hold the names `CV` and `PT`.

    Axes without limits (automatic) are not checked.

    Args:
        plot: the plot the names are drawn into.

    Raises:
        ValueError: if a name lies outside the limits, with the margin of its
            glyph, e.g. the axes of `lobule_axes` without `labels=True`.
    """
    spec = plot.spec()
    for name, (x, y) in (("CV", _CV), ("PT", _PT)):
        for axis, value in (("x", x), ("y", y)):
            low, high = spec[axis].get("min"), spec[axis].get("max")
            if (low is not None and value < low + _LABEL_MARGIN) or (
                high is not None and value > high - _LABEL_MARGIN
            ):
                raise ValueError(
                    f"The name {name} at {axis}={value} is outside the limits "
                    f"({low}, {high}) of the plot {plot.id!r}: make its axes "
                    "with lobule_axes(..., labels=True), which reserves the room."
                )


def draw_lobule(
    plot: ChartPlot,
    values: Sequence[float] | np.ndarray,
    *,
    scale: ColorScale,
    missing: str = MISSING,
    thirds: bool = False,
    anatomy: bool = False,
    labels: bool = False,
    points: int | None = None,
) -> None:
    """Draw a value per zone on the idealized hexagonal lobule, into a plot.

    The renderer of the lobule of this package, which `lobule_map` draws its
    one panel with and a model composes its own grid of maps from (its
    conditions, its measured rings, its legend stay in the model; a grid takes
    one colour bar of its scale through `colorbar`). The plot is made with the
    axes of `lobule_axes`, which hold the hexagon at one scale. The hexagon
    has the circumradius 1, the central vein at the origin and the portal
    triads at the corners on the angles 0, 60, ... degrees; a value is drawn
    as the band of its zone (`pkpdsim.zonation.render.lobule_bands`), the
    twelve wedges of a band as the rings of one polygon element whose id is
    `<plot id>-<zone id>` (`ids.zone`), bounded by the isolines of the lobular
    position at `(k - 1) / n` and `k / n`, so the area of a band is the area
    share of the zone, `zone_areas(n)`. The centre of this idealization is
    the central vein itself, so no disc of a vein lumen covers the pericentral
    bands. Everything is a vector polygon, line or marker.

    Args:
        plot: the plot, which the caller made with `lobule_axes` and owns.
        values: the value of every zone, periportal first; their number is the
            number of zones. A value which is not finite is a zone without
            data and is drawn in `missing`.
        scale: the colour of a value; the same scale goes to `colorbar`.
        missing: the colour of a zone without data, a grey by default (any
            colour `sbmlsim.plot.colors.normalize_color` reads). A gap of the
            hexagon would read as a defect of the geometry and not as a zone
            without data, so it is drawn.
        thirds: draw the isolines of 1/3 and 2/3 as dotted lines, the thirds of
            the lobulus the histology of a study counts its zones in.
        anatomy: draw the outline of the hexagon and the three dots of the
            portal triad of every corner (portal vein, hepatic artery, bile
            duct).
        labels: name the vessels, `CV` at the centre and `PT` below the lower
            right corner, for the first panel of a grid of several; `CV` is
            dark on a light centre and light on a dark one, so it reads on any
            map.
        points: the sampling of the isolines of the bands, the directions over
            the arc of an Apollonius circle (`lobule_bands`); a caller which
            draws many small panels samples more coarsely. `None` is the default
            of `lobule_bands`.

    Raises:
        ValueError: if there is no value, so no zone to draw, and for
            `labels` on a plot whose limits do not hold the names (the axes of
            `lobule_axes` without `labels=True`).
    """
    if labels:
        _check_label_room(plot)
    per_zone = np.asarray(values, dtype=float)
    bands = lobule_bands(len(per_zone), points=points)
    grey = _color(missing)
    colours = [grey if c is None else c for c in scale.colors(per_zone)]
    for k, polygons in enumerate(bands, start=1):
        plot.polygon(
            polygons,
            color=colours[k - 1],
            width=_SEAM_WIDTH,
            fill_opacity=1.0,
            id=f"{plot.id}-{ids.zone(k)}",
        )
    if thirds:
        for position in _THIRDS:
            plot.curve(
                *lobule_isoline(position).T,
                line="dotted",
                width=0.7,
                color=_THIRDS_GREY,
            )
    if anatomy:
        plot.curve(*lobule_isoline(0.0).T, width=1.1, color=_ANATOMY_GREY)
        _portal_triads(plot)
    if labels:
        light = _luminance(colours[-1]) > _LIGHT
        plot.text(*_CV, "CV", anchor="middle", color="#000000" if light else "#ffffff")
        plot.text(*_PT, "PT", anchor="middle", color="#000000")


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
) -> Chart:
    """Draw the value of every zone on the idealized hexagonal lobule.

    The central vein is the centre of the hexagon, the six portal triads are
    its corners, and a point has the normalized periportal to pericentral
    distance `d_pp / (d_pp + d_pc)` of the nearest triad and the central
    vein: the rim and the corners are zone 1 (periportal), the centre is zone
    `zones` (pericentral). The bands of the zones are bounded by the isolines
    of that distance at `k / zones`, an interval of the position per zone, so
    their areas are the area shares of the zones, `zone_areas(zones)`. The
    figure is vector: a polygon of the wedges per band, and a gradient as the
    colour bar.

    The panel is drawn by `draw_lobule`, with the outline of the hexagon and
    the portal triads of its corners, on the axes of `lobule_axes` at the
    aspect of its panel (`panel_aspects`); a model which composes its own grid
    of lobules calls these and `colorbar` directly.

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
        The chart, one plot of the hexagon with its colour bar. The colour
        scale spans the range of the values: a line axis starts at zero, a
        colour scale does not, because the zero is explicit on the ticks of a
        colour bar and a scale from zero would show a small gradient as one
        shade. A scale of one value is never degenerate (`_scale`).

    Raises:
        KeyError: if the result has no column of a zone.
        ValueError: for less than one zone, and for a number of zones which is
            not the one of the result: the bands of the map are the area
            shares of `zones`, so `zones=5` on a ten zone result would draw a
            lobule of five zones out of the five periportal ones of the
            result, a figure of nothing.
    """
    if zones < 1:
        raise ValueError(f"A lobule has at least one zone, not {zones}.")
    found = _zones(result, quantity, layer)
    if found and found != zones:
        raise ValueError(
            f"The result holds {found} zones of '{quantity}' in '{layer}' and "
            f"not the {zones} of the lobule: the map would draw the first "
            f"{zones} of them as the whole lobule."
        )
    name = quantity if name is None else name
    values = _at(result, observers.profile_ids(quantity, layer, zones), time)
    scale = ColorScale(*_scale([values], vmin, vmax))

    def build(aspect: float | None) -> Chart:
        """Build the map for the aspect of its panel."""
        chart = Chart(*_LOBULE_SIZE)
        x, y = lobule_axes(
            aspect,
            labels=False,
            label="Corners: portal triads (0, periportal)\n"
            "Centre: central vein (1, pericentral)",
        )
        plot = chart.plot(
            title=f"{name} {layer} at {_time_text(time, time_unit)}",
            x=x,
            y=y,
            legend=False,
        )
        draw_lobule(plot, values, scale=scale, anatomy=True)
        colorbar(plot, scale, _label(f"{name} {layer}", unit))
        return chart

    (aspect,) = panel_aspects(build(None)).values()
    return build(aspect)
