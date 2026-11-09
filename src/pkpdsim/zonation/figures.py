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
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.figure import Figure

from pkpdsim.zonation import ids, observers
from pkpdsim.zonation.geometry import Layer

__all__ = [
    "ZoneProfile",
    "data_axes",
    "layer_profile",
    "result_frame",
    "space_time",
    "zone_profiles",
]

#: the group id of the axes of a colorbar
COLORBAR = "colorbar"

_LINE_STYLES: tuple[str, ...] = (":", "--", "-", "-.")
_MARKERS: tuple[str, ...] = ("s", "D", "o", "^")
_LAYER_ORDER: tuple[str, ...] = (*(str(layer) for layer in Layer), ids.VASCULAR)
_UNIT = "mM"
"""The unit of the concentrations of the profile observers."""


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
    """Get the axes of the panels of a figure, the axes of a colorbar left out."""
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


def zone_profiles(
    result: pd.DataFrame,
    profiles: Sequence[ZoneProfile],
    *,
    times: Sequence[float],
    data: pd.DataFrame | None = None,
) -> Figure:
    """Draw a row per profile, with a line per time over the zones.

    Args:
        result: the simulation result.
        profiles: the profiles, a row each; the number of zones is the number
            of selections of a profile.
        times: the times of the lines, interpolated in the result; the line
            style is by time.
        data: optional data, a long frame with the columns `label` (of a
            profile), `time`, `position` (0 periportal to 1 pericentral),
            `value` and optionally `sd`, drawn over the lines as markers by
            time; the rows of a time which is not in `times` are not drawn.

    Returns:
        The figure.

    Raises:
        KeyError: if the result has no column of a selection.
    """
    n_rows = len(profiles)
    figure = Figure(figsize=(6.0, 3.0 * n_rows), layout="constrained")
    axes = figure.subplots(n_rows, 1, sharex=True, squeeze=False)[:, 0]
    for ax, profile in zip(axes, profiles, strict=True):
        x = _positions(len(profile.selections))
        for k, t in enumerate(times):
            ax.plot(
                x,
                _at(result, profile.selections, t),
                linestyle=_LINE_STYLES[k % len(_LINE_STYLES)],
                color="C0",
                linewidth=2.0,
                label=f"Sim {t:g}",
            )
        if data is not None:
            rows = data[data["label"] == profile.label]
            for k, t in enumerate(times):
                at = rows[np.isclose(rows["time"], t)]
                if at.empty:
                    continue
                ax.errorbar(
                    at["position"],
                    at["value"],
                    yerr=at.get("sd"),
                    linestyle="",
                    marker=_MARKERS[k % len(_MARKERS)],
                    markeredgecolor="black",
                    color="C0",
                    label=f"Data {t:g}",
                )
        ax.set_ylabel(_label(profile.label, profile.unit))
        ax.set_xlim(0.0, 1.0)
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    axes[-1].set_xlabel("Position along the sinusoid (0 periportal, 1 pericentral)")
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
) -> Figure:
    """Draw the concentration over the zones and the time, a heatmap per layer.

    Args:
        result: the simulation result.
        quantity: the id of the substance.
        layers: the layers, a heatmap each, from the top.
        outflow: add the concentration of the blood leaving the organ below.

    Returns:
        The figure.

    Raises:
        KeyError: if the result has no zone of a layer or no outflow.
    """
    _columns(result, ["time"])
    time = result["time"].to_numpy(dtype=float)
    rows = len(layers) + int(outflow)
    figure = Figure(figsize=(7.0, 2.6 * rows), layout="constrained")
    axes = figure.subplots(rows, 1, sharex=True, squeeze=False)[:, 0]
    for ax, layer in zip(axes, layers, strict=False):
        n = _zones(result, quantity, layer)
        if n == 0:
            raise KeyError(f"The result has no zone of '{quantity}' in '{layer}'.")
        values = result[observers.profile_ids(quantity, layer, n)].to_numpy(dtype=float)
        mesh = ax.pcolormesh(
            time, _positions(n), values.T, shading="nearest", rasterized=False
        )
        _colorbar(figure, mesh, ax, _label(f"{quantity} {layer}", _UNIT))
        ax.set_ylabel("Position")
        ax.set_ylim(0.0, 1.0)
    if outflow:
        sid = observers.outflow_concentration(quantity)
        _columns(result, [sid])
        ax = axes[-1]
        ax.plot(time, result[sid].to_numpy(dtype=float), color="C0", linewidth=2.0)
        ax.set_ylabel(_label(f"{quantity} outflow", _UNIT))
    axes[-1].set_xlabel("Time")
    axes[-1].set_xlim(float(time.min()), float(time.max()))
    return figure


def layer_profile(result: pd.DataFrame, substance: str, *, time: float) -> Figure:
    """Draw the profile over the zones in every layer which carries the substance.

    The layers are those of which the result holds the substance in zones
    from `z1` on, so a well-stirred build shows its vascular site and a build of
    the radial geometry its layers.

    Args:
        result: the simulation result.
        substance: the id of the substance.
        time: the time of the profile, interpolated in the result.

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
    axes = figure.subplots(nrows, ncols, squeeze=False)
    flat = list(axes.flat)
    for ax in flat[len(present) :]:
        ax.remove()
    for ax, (layer, n) in zip(flat, present, strict=False):
        selections = observers.profile_ids(substance, layer, n)
        ax.plot(_positions(n), _at(result, selections, time), color="C0", linewidth=2.0)
        ax.set_title(str(layer))
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel("Position")
        ax.set_ylabel(_label(substance, _UNIT))
    figure.suptitle(f"{substance} at {time:g}")
    return figure
