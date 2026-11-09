"""The figures of a zonated model: the profiles over the zones, in space and time.

The functions draw a simulation result of any zonated model, whatever built it,
and know nothing but the id scheme of the engine (`pkpdsim.zonation.ids`,
`pkpdsim.zonation.observers`). The result is a `pandas.DataFrame` with a `time`
column and a column per selection named by its plain id (`glc_cell_z2`);
`result_frame` converts a roadrunner result. The number of zones and the layers
which carry a substance are read from the ids present in the result.

A figure is built directly as a `matplotlib.figure.Figure` and owned by the
caller; `matplotlib.pyplot` and its global state are not used. Every figure is
vector: the heatmap is not rasterized and neither is the solid of its colorbar,
so the SVG holds no embedded image (a viewer which drops embedded raster data
would show an empty panel).
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import Colorbar
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.patches import Polygon

from pkpdsim.zonation import ids, observers
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.render import lobule_bands

__all__ = [
    "ZoneProfile",
    "data_axes",
    "layer_profile",
    "lobule_map",
    "result_frame",
    "space_time",
    "zone_profiles",
]

#: the group id of the axes of a colorbar
COLORBAR = "colorbar"

_LINE_STYLES: tuple[str, ...] = (":", "--", "-", "-.")
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


@dataclass(frozen=True)
class ZoneProfile:
    """A profile over the zones, a row of `zone_profiles`.

    Attributes:
        label: the name of the quantity, the label of the y axis.
        unit: the unit of the quantity.
        selections: the ids of the result by zone, periportal first.
    """

    label: str
    unit: str
    selections: tuple[str, ...]


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


def _colorbar(figure: Figure, mesh: object, ax: Axes, label: str) -> Colorbar:
    """Add a colorbar which is vector and is no data axes."""
    bar = figure.colorbar(mesh, ax=ax, label=label)  # ty: ignore[invalid-argument-type]
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


def zone_profiles(
    result: pd.DataFrame,
    profiles: Sequence[ZoneProfile],
    *,
    times: Sequence[float],
    data: pd.DataFrame | None = None,
    time_unit: str = DEFAULT_TIME_UNIT,
) -> Figure:
    """Draw a row per profile, with a line per time over the zones.

    Every time has its colour (a sequential map, early to late) and its line
    style; the data of a time are markers of the colour of that time. The
    legend is one for the figure, right of the grid, if every row has the same
    entries, else one per row.

    A line axis starts at zero when no value is negative; a colour scale
    (`space_time`, `lobule_map`) spans the range of its values instead.

    Args:
        result: the simulation result.
        profiles: the profiles, a row each; the number of zones is the number
            of selections of a profile.
        times: the times of the lines, interpolated in the result.
        data: optional data, a long frame with the columns `label` (of a
            profile), `time`, `position` (0 periportal to 1 pericentral),
            `value` (in the unit of the profile) and optionally `sd` (a
            symmetric error in the same unit), drawn over the lines as
            markers by time. The rows of a time which is not in `times` are
            skipped, as plotting a subset of the times of a dataset is
            legitimate.
        time_unit: the unit of the times, the convention of the engine (min).

    Returns:
        The figure.

    Raises:
        KeyError: if the result has no column of a selection.
        ValueError: if `profiles` or `times` is empty, if a row of `data` has a
            label which no profile has, or if no row of `data` is drawn.
    """
    if not profiles:
        raise ValueError("The profiles are empty, there is no row to draw.")
    if not len(times):
        raise ValueError("The times are empty, there is no line to draw.")
    if data is not None:
        labels = [profile.label for profile in profiles]
        unknown = sorted(set(data["label"]) - set(labels))
        if unknown:
            raise ValueError(
                f"The data has the labels {unknown} which no profile has; "
                f"the labels of the profiles are {labels}."
            )
    drawn = 0
    colors = colormaps["viridis"](np.linspace(0.1, 0.75, len(times)))
    n_rows = len(profiles)
    figure = Figure(figsize=(7.0, 3.0 * n_rows), layout="constrained")
    axes = figure.subplots(n_rows, 1, sharex=True, squeeze=False)[:, 0]
    for ax, profile in zip(axes, profiles, strict=True):
        x = _positions(len(profile.selections))
        values: list[np.ndarray] = []
        for k, t in enumerate(times):
            y = _at(result, profile.selections, t)
            values.append(y)
            ax.plot(
                x,
                y,
                linestyle=_LINE_STYLES[k % len(_LINE_STYLES)],
                color=colors[k],
                linewidth=2.0,
                label=f"Sim {_time_text(t, time_unit)}",
            )
        if data is not None:
            rows = data[data["label"] == profile.label]
            for k, t in enumerate(times):
                at = rows[np.isclose(rows["time"], t)]
                if at.empty:
                    continue
                drawn += len(at)
                err = at["sd"].to_numpy(dtype=float) if "sd" in at else 0.0
                values.append(at["value"].to_numpy(dtype=float) + err)
                ax.errorbar(
                    at["position"],
                    at["value"],
                    yerr=at.get("sd"),
                    linestyle="",
                    marker=_MARKERS[k % len(_MARKERS)],
                    markeredgecolor="black",
                    color=colors[k],
                    label=f"Data {_time_text(t, time_unit)}",
                )
        limits = _limits(values, _HEADROOM)
        if limits is not None:
            ax.set_ylim(*limits)
        ax.set_ylabel(_label(profile.label, profile.unit))
        ax.set_xlim(0.0, 1.0)
    if data is not None and not drawn:
        raise ValueError(
            f"No row of the data is drawn: its times {sorted(set(data['time']))} "
            f"include none of the times {list(times)} of the figure."
        )
    entries = [ax.get_legend_handles_labels() for ax in axes]
    if all(e[1] == entries[0][1] for e in entries):
        figure.legend(*entries[0], loc="outside right center")
    else:
        for ax in axes:
            ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    axes[-1].set_xlabel(POSITION)
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
            the initial condition in one zone, dominates. It must be below `vmax`, so a `vmin` equal
            to the value of constant data raises, where the default scale of
            such data is the fallback zero to the value.
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
        # (`_colorbar` un-rasterizes its solids, the default is an image) and
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
        _colorbar(figure, mesh, ax, _label(f"{name} {layer}", unit))
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
            the initial condition in one zone, dominates. It must be below `vmax`, so a `vmin` equal
            to the value of constant data raises, where the default scale of
            such data is the fallback zero to the value.
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
    for k, polygons in enumerate(lobule_bands(zones), start=1):
        colour = mappable.to_rgba(values[k - 1])
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
    ax.set_xlim(-1.02, 1.02)
    ax.set_ylim(-0.9, 0.9)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel(
        "Corners: portal triads (0, periportal)\nCentre: central vein (1, pericentral)"
    )
    ax.set_title(f"{name} {layer} at {_time_text(time, time_unit)}")
    _colorbar(figure, mappable, ax, _label(f"{name} {layer}", unit))
    return figure
