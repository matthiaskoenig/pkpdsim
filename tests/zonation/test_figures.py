"""The zone figures: over a synthetic result, and over the result of a built chain."""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection, PolyCollection, QuadMesh
from matplotlib.colors import Normalize, to_rgba
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.patheffects import withStroke
from matplotlib.text import Text
from pbpk.conftest import individual
from sbmlutils.factory import create_model

from pkpdsim.zonation import ids, observers
from pkpdsim.zonation.figures import (
    COLORBAR,
    MISSING,
    POSITION,
    Band,
    ZoneColumn,
    ZoneProfile,
    colorbar,
    data_axes,
    draw_lobule,
    layer_profile,
    lobule_map,
    result_frame,
    space_time,
    zone_profiles,
)
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.observers import add_observers, observer_objects
from pkpdsim.zonation.render import lobule_isoline

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


def _x(line: Line2D) -> np.ndarray:
    """The x data of a line as an array."""
    return np.asarray(line.get_xdata(), dtype=float)


def _y(line: Line2D) -> np.ndarray:
    """The y data of a line as an array."""
    return np.asarray(line.get_ydata(), dtype=float)


def _cells(mesh: object) -> np.ndarray:
    """The values of a heatmap as an array."""
    assert isinstance(mesh, QuadMesh)
    return np.asarray(mesh.get_array(), dtype=float)


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


def test_the_profile_grid_has_a_row_per_profile_with_label_unit_and_a_line_per_time() -> (
    None
):
    fig = zone_profiles(_result(), _profiles(3), times=(0, 60, 120))
    axes = data_axes(fig)
    assert len(axes) == 3
    for i, ax in enumerate(axes):
        assert ax.get_ylabel() == f"p{i} [mM]"
        assert len(ax.lines) == 3


def test_the_x_axis_is_the_zone_midpoint() -> None:
    fig = zone_profiles(_result(), _profiles(1), times=(0,))
    assert fig.axes[0].lines[0].get_xdata() == pytest.approx(
        [(k - 0.5) / 10 for k in range(1, 11)]
    )


def test_the_profile_values_are_the_result_at_the_time() -> None:
    fig = zone_profiles(_result(), _profiles(1), times=(60,))
    expected = [(2 * k * 0.5) for k in range(1, N + 1)]  # plasma is the second layer
    assert fig.axes[0].lines[0].get_ydata() == pytest.approx(expected)


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


def _error_bars(ax: Axes) -> list[np.ndarray]:
    """The segments of the error bars of an axes as arrays."""
    out = []
    for collection in ax.collections:
        assert isinstance(collection, LineCollection)
        out.extend(np.asarray(seg, dtype=float) for seg in collection.get_segments())
    return out


def test_the_data_is_overlaid_where_given() -> None:
    fig = zone_profiles(_result(), _profiles(2), times=(0, 60), data=_data())
    first, second = data_axes(fig)
    # the two lines of the simulation, then the markers of the data
    assert len(first.lines) == 3
    assert list(_x(first.lines[2])) == [0.25, 0.75]
    assert list(_y(first.lines[2])) == [1.0, 2.0]
    assert list(_x(second.lines[2])) == [0.5]
    assert list(_y(second.lines[2])) == [3.0]
    # the error bars are value -/+ sd, per row
    bars = _error_bars(first)
    assert [float(v) for seg in bars for v in (seg[0, 1], seg[1, 1])] == pytest.approx(
        [0.9, 1.1, 1.8, 2.2]
    )
    bars = _error_bars(second)
    assert [float(v) for seg in bars for v in (seg[0, 1], seg[1, 1])] == pytest.approx(
        [2.7, 3.3]
    )


def test_the_data_of_another_time_is_skipped() -> None:
    data = _data(time=[60.0, 30.0, 60.0])
    fig = zone_profiles(_result(), _profiles(2), times=(0, 60), data=data)
    first, second = data_axes(fig)
    assert list(_x(first.lines[2])) == [0.25]
    assert list(_x(second.lines[2])) == [0.5]


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
    fig = space_time(_result(), "drg", layers=(Layer.PLASMA,), unit="uM")
    assert data_axes(fig)[1].get_ylabel() == "drg outflow [uM]"
    assert fig.axes[-1].get_ylabel() == "drg plasma [uM]"
    assert layer_profile(_result(), "drg", time=0, unit="uM").get_supylabel() == (
        "drg [uM]"
    )
    assert layer_profile(_result(), "drg", time=0).get_supylabel() == "drg [mM]"


def test_a_missing_selection_names_the_column() -> None:
    bad = [ZoneProfile("x", "mM", ("nope_z1",))]
    with pytest.raises(KeyError, match="nope_z1"):
        zone_profiles(_result(), bad, times=(0,))


def test_the_space_time_panel_carries_the_outflow_below() -> None:
    fig = space_time(_result(), "drg", layers=(Layer.PLASMA,))
    heat, out = data_axes(fig)
    assert len(fig.axes) == 3  # and the colorbar
    assert len(heat.collections) == 1
    assert _cells(heat.collections[0]).size == 25 * N
    assert _y(out.lines[0])[-1] == pytest.approx(1.0)
    assert out.get_xlabel() == "Time [min]"


def test_the_heatmap_has_the_zones_periportal_first_and_the_time_to_the_right() -> None:
    """The plasma of the fixture is `2 k t / 120` in zone `k` at time `t`."""
    fig = space_time(_result(), "drg", layers=(Layer.PLASMA,))
    heat = data_axes(fig)[0]
    mesh = heat.collections[0]
    cells = _cells(mesh).reshape(N, 25)  # a row per zone, a column per time
    assert cells[0, -1] == pytest.approx(2.0)  # z1 at 120
    assert cells[N - 1, -1] == pytest.approx(2.0 * N)  # z10 at 120
    assert cells[N - 1, 12] == pytest.approx(N)  # z10 at 60
    assert cells[0, 0] == pytest.approx(0.0)  # z1 at 0
    # the periportal zone is lowest, as the axis (0 periportal) says
    assert cells[N - 1, -1] > cells[0, -1]
    assert isinstance(mesh, QuadMesh)
    corners = np.asarray(mesh.get_coordinates(), dtype=float)
    assert corners.shape == (N + 1, 26, 2)
    assert corners[0, 0, 1] == pytest.approx(0.0)
    assert corners[1, 0, 1] == pytest.approx(0.1)
    assert corners[-1, 0, 1] == pytest.approx(1.0)
    centers = (corners[0, :-1, 0] + corners[0, 1:, 0]) / 2
    assert centers == pytest.approx(np.linspace(0.0, 120.0, 25))
    assert heat.get_ylim() == (0.0, 1.0)
    assert heat.get_xlim() == pytest.approx((0.0, 120.0))


def test_the_space_time_without_the_outflow_has_a_panel_per_layer() -> None:
    fig = space_time(_result(), "drg", layers=(Layer.PLASMA, Layer.CELL), outflow=False)
    assert len(data_axes(fig)) == 2
    assert len(fig.axes) == 4


def test_the_layer_profile_shows_every_layer_that_carries_substance() -> None:
    fig = layer_profile(_result(), "drg", time=60)
    axes = data_axes(fig)
    assert len(axes) == 6
    assert [ax.get_title() for ax in axes] == [str(layer) for layer in LAYERS]
    for i, ax in enumerate(axes):
        assert ax.lines[0].get_ydata() == pytest.approx(
            [(i + 1) * k * 0.5 for k in range(1, N + 1)]
        )


def test_the_layer_profile_shows_only_the_layers_present() -> None:
    fig = layer_profile(_result((ids.VASCULAR, Layer.CELL)), "drg", time=0)
    assert [ax.get_title() for ax in data_axes(fig)] == ["cell", "vas"]


def test_the_layer_profile_of_an_absent_substance_raises() -> None:
    with pytest.raises(KeyError, match="other"):
        layer_profile(_result(), "other", time=0)


def test_the_figures_are_owned_not_pyplot_state() -> None:
    import matplotlib.pyplot as plt

    before = plt.get_fignums()
    fig = zone_profiles(_result(), _profiles(1), times=(0,))
    assert isinstance(fig, Figure)
    assert plt.get_fignums() == before


def test_the_figures_save_as_svg(tmp_path: Path) -> None:
    zone_profiles(_result(), _profiles(1), times=(0,)).savefig(tmp_path / "f.svg")
    assert (tmp_path / "f.svg").read_text().startswith("<?xml")


def test_no_figure_embeds_a_raster_image(tmp_path: Path) -> None:
    """The heatmap and its colorbar are vector, the maintainer's viewer drops images."""
    result = _result()
    figures = {
        "profiles": zone_profiles(result, _profiles(2), times=(0, 60)),
        "space_time": space_time(result, "drg", layers=(Layer.PLASMA, Layer.CELL)),
        "layers": layer_profile(result, "drg", time=60),
        "lobule": lobule_map(result, "drg", time=60, zones=N, layer=Layer.CELL),
    }
    for name, fig in figures.items():
        path = tmp_path / f"{name}.svg"
        fig.savefig(path)
        text = path.read_text()
        assert "<image" not in text, name
    # the heatmap is drawn: a path per cell
    assert (tmp_path / "space_time.svg").read_text().count("<path") > 25 * N


def test_a_rasterized_colorbar_would_embed_an_image(tmp_path: Path) -> None:
    """The test above can fail: the matplotlib default of the colorbar is an image."""
    fig = Figure()
    ax = fig.subplots()
    mesh = ax.pcolormesh(np.arange(12.0).reshape(3, 4))
    fig.colorbar(mesh, ax=ax)
    fig.savefig(tmp_path / "default.svg")
    assert "<image" in (tmp_path / "default.svg").read_text()


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
    fig = zone_profiles(result, [profile], times=(0.5, 5))
    ax = data_axes(fig)[0]
    assert len(_x(ax.lines[0])) == n
    for line in ax.lines:
        y = np.asarray(line.get_ydata())
        assert np.isfinite(y).all()
    assert np.abs(ax.lines[0].get_ydata()).max() > 0

    fig = layer_profile(result, "tr", time=1.0)
    axes = data_axes(fig)
    assert [ax.get_title() for ax in axes] == LITERAL_LAYERS
    for ax in axes:
        assert len(_x(ax.lines[0])) == n
        y = np.asarray(ax.lines[0].get_ydata())
        assert np.isfinite(y).all()
    assert any(np.abs(ax.lines[0].get_ydata()).max() > 0 for ax in axes)

    fig = space_time(result, "tr", layers=(Layer.PLASMA,))
    heat, flow_out = data_axes(fig)
    values = _cells(heat.collections[0])
    assert values.size == n * 51
    assert np.isfinite(values).all()
    assert values.max() > 0
    assert np.isfinite(flow_out.lines[0].get_ydata()).all()
    fig.savefig(tmp_path / "chain.svg")
    assert "<image" not in (tmp_path / "chain.svg").read_text()


def _drawn(fig: Figure) -> Figure:
    """Lay a figure out on a canvas which draws, as the tests of the label and legend positions need."""
    FigureCanvasAgg(fig).draw()
    return fig


def _legend(ax: Axes) -> Legend:
    """The legend of an axes, which must exist."""
    legend = ax.get_legend()
    assert legend is not None
    return legend


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


def test_every_position_label_names_the_ends_of_the_sinusoid() -> None:
    end = "(0 periportal, 1 pericentral)"
    fig = zone_profiles(_result(), _profiles(2), times=(0,))
    assert data_axes(fig)[-1].get_xlabel() == f"Position {end}"
    heat = data_axes(space_time(_result(), "drg", layers=(Layer.PLASMA,)))[0]
    assert heat.get_ylabel().replace("\n", " ") == f"Position {end}"
    assert layer_profile(_result(), "drg", time=0).get_supxlabel() == (
        f"Position {end}"
    )


def test_the_time_carries_its_unit_in_every_figure() -> None:
    result = _result()
    fig = zone_profiles(result, _profiles(1), times=(0, 60), data=None)
    legend = fig.legends[0]
    assert [t.get_text() for t in legend.get_texts()] == [
        "Sim t = 0 min",
        "Sim t = 60 min",
    ]
    fig = zone_profiles(result, _profiles(1), times=(60,), time_unit="h")
    assert fig.legends[0].get_texts()[0].get_text() == "Sim t = 60 h"
    heat = space_time(result, "drg", layers=(Layer.PLASMA,))
    assert data_axes(heat)[-1].get_xlabel() == "Time [min]"
    heat = space_time(result, "drg", layers=(Layer.PLASMA,), time_unit="h")
    assert data_axes(heat)[-1].get_xlabel() == "Time [h]"
    layers = layer_profile(result, "drg", time=60)
    assert layers.get_suptitle() == "drg at t = 60 min"
    assert layer_profile(result, "drg", time=2, time_unit="h").get_suptitle() == (
        "drg at t = 2 h"
    )


def test_a_layer_of_zeros_is_drawn_from_zero_to_one() -> None:
    result = _flat()
    for profile, ax in zip(
        _plasma_cell(),
        data_axes(zone_profiles(result, _plasma_cell(), times=(5,))),
        strict=True,
    ):
        assert ax.get_ylim()[0] == 0.0, profile.label
    cell = data_axes(zone_profiles(result, _plasma_cell(), times=(5,)))[1]
    assert cell.get_ylim() == (0.0, 1.0)
    fig = space_time(result, "drg", layers=(Layer.PLASMA, Layer.CELL))
    plasma, cells, out = data_axes(fig)
    assert cells.collections[0].get_clim() == (0.0, 1.0)
    assert plasma.collections[0].get_clim() == (0.0, 3.0)
    assert out.get_ylim() == pytest.approx((0.0, 3.0 * 1.05))
    axes = {
        ax.get_title(): ax
        for ax in data_axes(layer_profile(result, "drg", time=5, share_y=False))
    }
    assert axes["cell"].get_ylim() == (0.0, 1.0)
    assert axes["plasma"].get_ylim() == pytest.approx((0.0, 3.0 * 1.05))


def test_a_constant_layer_is_drawn_from_zero_to_its_value() -> None:
    """Zero to one would hide a constant above one, so the scale runs to the value."""
    fig = space_time(_flat(plasma=0.25), "drg", layers=(Layer.PLASMA,))
    assert data_axes(fig)[0].collections[0].get_clim() == (0.0, 0.25)


def test_a_negative_value_leaves_the_axis_to_matplotlib() -> None:
    result = _flat()
    result[observers.profile_ids("drg", "cell", 4)] = -1.0
    ax = data_axes(zone_profiles(result, _plasma_cell(), times=(5,)))[1]
    assert ax.get_ylim()[0] < -1.0


def test_every_time_has_a_colour_and_its_data_has_the_colour_of_the_time() -> None:
    data = _data(time=[0.0, 60.0, 60.0], label=["p0", "p0", "p0"])
    fig = zone_profiles(_result(), _profiles(1), times=(0, 60, 120), data=data.iloc[:2])
    ax = data_axes(fig)[0]
    lines = ax.lines[:3]
    colors = [line.get_color() for line in lines]
    assert len({tuple(c) for c in colors}) == 3
    assert [line.get_linestyle() for line in lines] == [":", "--", "-"]
    # the markers of 0 and of 60 are drawn in the colour of their curve
    assert tuple(ax.lines[3].get_color()) == tuple(colors[0])
    assert tuple(ax.lines[4].get_color()) == tuple(colors[1])


def test_the_legend_is_one_when_the_rows_agree_and_one_per_row_when_they_differ() -> (
    None
):
    same = zone_profiles(_result(), _profiles(2), times=(0, 60))
    assert len(same.legends) == 1
    assert all(ax.get_legend() is None for ax in data_axes(same))
    data = _data(label=["p0", "p0", "p0"])
    differ = zone_profiles(_result(), _profiles(2), times=(0, 60), data=data)
    assert differ.legends == []
    first, second = data_axes(differ)
    assert len(_legend(first).get_texts()) == 3
    assert len(_legend(second).get_texts()) == 2


def test_the_layer_grid_carries_each_label_once() -> None:
    for layers in (LAYERS, LAYERS[:5]):
        fig = layer_profile(_result(layers), "drg", time=0)
        assert fig.get_supxlabel() == "Position (0 periportal, 1 pericentral)"
        assert fig.get_supylabel() == "drg [mM]"
        assert [ax.get_xlabel() for ax in data_axes(fig)] == [""] * len(layers)
        assert [ax.get_ylabel() for ax in data_axes(fig)] == [""] * len(layers)


def test_the_layers_share_one_concentration_scale_unless_told_otherwise() -> None:
    result = _flat()
    shared = data_axes(layer_profile(result, "drg", time=5))
    assert {tuple(round(v, 6) for v in ax.get_ylim()) for ax in shared} == {(0.0, 3.15)}
    own = {
        ax.get_title(): ax.get_ylim()
        for ax in data_axes(layer_profile(result, "drg", time=5, share_y=False))
    }
    assert own["plasma"] == pytest.approx((0.0, 3.15))
    assert own["cell"] == (0.0, 1.0)
    zero = _flat(plasma=0.0)
    limits = {ax.get_ylim() for ax in data_axes(layer_profile(zero, "drg", time=5))}
    assert limits == {(0.0, 1.0)}


def test_the_single_legend_is_centred_on_the_grid() -> None:
    fig = zone_profiles(_result(), _profiles(2), times=(0, 60))
    _drawn(fig)
    box = fig.legends[0].get_window_extent()
    axes = data_axes(fig)
    top = max(ax.get_window_extent().y1 for ax in axes)
    bottom = min(ax.get_window_extent().y0 for ax in axes)
    # centred on the figure, which is within the span of the grid
    assert (box.y0 + box.y1) / 2 == pytest.approx(fig.bbox.height / 2, rel=0.01)
    assert bottom < box.y0 < box.y1 < top
    assert box.x0 > max(ax.get_window_extent().x1 for ax in axes)


def test_the_outflow_panel_is_lower_than_a_heatmap_and_the_name_is_chosen() -> None:
    fig = space_time(_result(), "drg", layers=(Layer.PLASMA,), name="Drug")
    heat, out = data_axes(fig)
    spec = heat.get_subplotspec()
    assert spec is not None
    assert spec.get_gridspec().get_height_ratios() == [1.0, 0.6]
    assert out.get_ylabel() == "Drug outflow [mM]"
    assert fig.axes[-1].get_ylabel() == "Drug plasma [mM]"
    default = space_time(_result(), "drg", layers=(Layer.PLASMA,))
    assert default.axes[-1].get_ylabel() == "drg plasma [mM]"


def test_the_outflow_panel_is_drawn_at_0p6_of_a_heatmap() -> None:
    fig = _drawn(space_time(_result(), "drg", layers=(Layer.PLASMA,)))
    heat, out = data_axes(fig)
    ratio = out.get_window_extent().height / heat.get_window_extent().height
    assert ratio == pytest.approx(0.6, rel=0.05)


def test_a_shared_scale_labels_the_ticks_of_the_outer_panels_only() -> None:
    def ticks(ax: Axes) -> tuple[bool, bool]:
        return (
            bool(ax.yaxis.get_tick_params()["labelleft"]),
            bool(ax.xaxis.get_tick_params()["labelbottom"]),
        )

    shared = [ticks(ax) for ax in data_axes(layer_profile(_result(), "drg", time=0))]
    assert shared == [
        (True, False),
        (False, False),
        (False, False),
        (True, True),
        (False, True),
        (False, True),
    ]
    five = [
        ticks(ax) for ax in data_axes(layer_profile(_result(LAYERS[:5]), "drg", time=0))
    ]
    assert five[2] == (False, True)  # above a gap, so on its own bottom edge
    own = data_axes(layer_profile(_result(), "drg", time=0, share_y=False))
    assert all(ticks(ax) == (True, True) for ax in own)


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


def _fills(ax: Axes) -> list[PolyCollection]:
    """The filled bands of an axes."""
    return [c for c in ax.collections if isinstance(c, PolyCollection)]


def test_a_frame_draws_one_column_as_before() -> None:
    fig = zone_profiles(_result(), _profiles(2), times=(0, 60))
    axes = data_axes(fig)
    assert len(axes) == 2
    assert tuple(fig.get_size_inches()) == pytest.approx((7.0, 6.0))
    assert [ax.get_title() for ax in axes] == ["", ""]
    assert len(fig.legends) == 1
    assert [ax.get_ylabel() for ax in axes] == ["p0 [mM]", "p1 [mM]"]
    assert axes[-1].get_xlabel() == POSITION
    assert _fills(axes[0]) == []


def test_a_column_per_condition_has_its_label_as_the_title_of_the_top_row() -> None:
    profiles = _profiles(2, n=4)
    times = (0.0, 60.0)
    columns = [
        _column("1 ml/min", profiles=profiles, times=times, level=1.0),
        _column("3 ml/min", profiles=profiles, times=times, level=2.0),
    ]
    fig = zone_profiles(columns, profiles, times=times)
    axes = data_axes(fig)
    assert len(axes) == 4
    assert tuple(fig.get_size_inches()) == pytest.approx((3.6 * 2 + 2.4, 6.0))
    assert [ax.get_title() for ax in axes] == ["1 ml/min", "3 ml/min", "", ""]
    assert axes[0].title.get_fontweight() == "bold"
    # the y label on the first column of every row, the x label on the bottom row
    assert [ax.get_ylabel() for ax in axes] == ["p0 [mM]", "", "p1 [mM]", ""]
    assert [ax.get_xlabel() for ax in axes] == ["", "", POSITION, POSITION]
    # the values of the column are drawn, the level of the second twice the first
    assert _y(axes[1].lines[1]) == pytest.approx(2.0 * _y(axes[0].lines[1]))


def test_the_lines_of_a_column_take_its_colour_and_the_style_of_their_time() -> None:
    profiles = _profiles(1, n=4)
    times = (0.0, 60.0, 120.0)
    columns = [
        _column("A", profiles=profiles, times=times, color="#ff0000"),
        _column("B", profiles=profiles, times=times, color="#0000ff", level=2.0),
    ]
    first, second = data_axes(zone_profiles(columns, profiles, times=times))
    assert [line.get_color() for line in first.lines] == ["#ff0000"] * 3
    assert [line.get_color() for line in second.lines] == ["#0000ff"] * 3
    # the style is the time, so the same styles in both columns
    assert [line.get_linestyle() for line in first.lines] == [":", "--", "-"]
    assert [line.get_linestyle() for line in second.lines] == [":", "--", "-"]


def test_a_column_without_a_colour_takes_the_colour_of_the_time() -> None:
    profiles = _profiles(1, n=4)
    times = (0.0, 60.0, 120.0)
    columns = [_column("A", profiles=profiles, times=times)]
    ax = data_axes(zone_profiles(columns, profiles, times=times))[0]
    frame = data_axes(zone_profiles(_result(n=4), profiles, times=times))[0]
    drawn = [tuple(line.get_color()) for line in ax.lines]
    assert len(set(drawn)) == 3
    assert drawn == [tuple(line.get_color()) for line in frame.lines]
    expected = colormaps["viridis"](np.linspace(0.1, 0.75, 3))
    assert np.asarray(drawn) == pytest.approx(expected)


def test_a_band_is_filled_between_its_ends() -> None:
    profiles = _profiles(1, n=4)
    columns = [
        _column("A", profiles=profiles, times=(60.0,), color="#ff0000", spread=0.5)
    ]
    ax = data_axes(zone_profiles(columns, profiles, times=(60.0,)))[0]
    fills = _fills(ax)
    assert len(fills) == 1
    fill = fills[0]
    assert fill.get_alpha() == pytest.approx(0.25)
    assert np.asarray(fill.get_facecolor())[0] == pytest.approx((1.0, 0.0, 0.0, 0.25))
    assert fill.get_linewidth() == pytest.approx([0.0])
    median = _y(ax.lines[0])
    y = np.asarray(fill.get_paths()[0].vertices, dtype=float)[:, 1]
    assert y.min() == pytest.approx(median.min() - 0.5)
    assert y.max() == pytest.approx(median.max() + 0.5)
    # the band is in the limits of the row, the median alone would cut it off
    assert ax.get_ylim()[1] > median.max() + 0.5


def test_a_band_without_ends_fills_nothing() -> None:
    profiles = _profiles(1, n=4)
    columns = [_column("A", profiles=profiles, times=(60.0,), color="#ff0000")]
    ax = data_axes(zone_profiles(columns, profiles, times=(60.0,)))[0]
    assert len(ax.lines) == 1
    assert _fills(ax) == []


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
    fig = zone_profiles(columns, profiles, times=times)
    top, bottom = data_axes(fig)
    expected = _medians(covering, every, 1.0)
    assert [_y(line) for line in top.lines] == [
        pytest.approx(expected[("Drug", t)]) for t in times
    ]
    assert [_y(line) for line in bottom.lines] == [
        pytest.approx(expected[("Dead", t)]) for t in (30.0, 90.0)
    ]
    # the legend of a row names the times of that row
    assert [t.get_text() for t in _legend(top).get_texts()] == [
        "Sim t = 0 min",
        "Sim t = 60 min",
    ]
    assert [t.get_text() for t in _legend(bottom).get_texts()] == [
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
    ax = data_axes(zone_profiles(result, [profile], times=(5.0,)))[0]
    assert _y(ax.lines[0]) == pytest.approx([20.0] * 4)
    assert ax.get_ylabel() == "Dead cells [%]"


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
    first, second = data_axes(zone_profiles(columns, profiles, times=times, data=data))
    # the row of the column A and the row of no column in A, the row of no column in B
    assert list(_x(first.lines[1])) == [0.25, 0.75]
    assert list(_x(second.lines[1])) == [0.75]
    # the markers take the colour of their column
    assert first.lines[1].get_color() == "#ff0000"
    assert second.lines[1].get_color() == "#0000ff"


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
    fig = zone_profiles(columns, profiles, times=times)
    first, second, third, fourth = data_axes(fig)
    assert first.get_ylim() == second.get_ylim()
    assert third.get_ylim() == fourth.get_ylim()
    # the limit of a row is the largest value of its columns, the second here
    expected = _medians(profiles, times, 3.0)
    assert first.get_ylim()[1] == pytest.approx(expected[("p0", 60.0)].max() * 1.05)
    assert third.get_ylim()[1] == pytest.approx(expected[("p1", 60.0)].max() * 1.05)
    assert first.get_ylim() != third.get_ylim()


def test_the_svg_of_a_grid_holds_no_image(tmp_path: Path) -> None:
    profiles = _profiles(2, n=4)
    times = (0.0, 60.0)
    columns = [
        _column("A", profiles=profiles, times=times, color="#ff0000", spread=0.5),
        _column("B", profiles=profiles, times=times, color="#0000ff", spread=0.5),
    ]
    path = tmp_path / "grid.svg"
    zone_profiles(columns, profiles, times=times).savefig(path)
    text = path.read_text()
    assert text.startswith("<?xml")
    assert "<image" not in text


def _handles(ax: Axes) -> list[Line2D]:
    """The handles of the legend of an axes, which must exist."""
    handles = _legend(ax).legend_handles
    assert all(isinstance(h, Line2D) for h in handles)
    return [h for h in handles if isinstance(h, Line2D)]


def _handle_colors(ax: Axes) -> np.ndarray:
    """The colours of the handles of the legend of an axes, as RGBA."""
    return np.asarray([to_rgba(h.get_color()) for h in _handles(ax)])


def test_the_legend_handles_are_neutral_where_the_columns_carry_the_colour() -> None:
    profiles = _profiles(1, n=4)
    times = (0.0, 60.0)
    columns = [
        _column("A", profiles=profiles, times=times, color="#ff0000"),
        _column("B", profiles=profiles, times=times, color="#0000ff"),
    ]
    ax = data_axes(zone_profiles(columns, profiles, times=times))[-1]
    # the titles name the columns, so a handle stands for its time, not a column
    assert _handle_colors(ax) == pytest.approx(np.vstack([to_rgba("0.3")] * 2))
    assert [h.get_linestyle() for h in _handles(ax)] == [":", "--"]


def test_the_legend_handles_take_the_colour_of_the_time_where_no_column_has_one() -> (
    None
):
    profiles = _profiles(1, n=4)
    times = (0.0, 60.0)
    columns = [
        _column("A", profiles=profiles, times=times),
        _column("B", profiles=profiles, times=times, level=2.0),
    ]
    axes = data_axes(zone_profiles(columns, profiles, times=times))
    expected = colormaps["viridis"](np.linspace(0.1, 0.75, 2))
    assert _handle_colors(axes[-1]) == pytest.approx(expected)
    # and the handles are the colours the lines of the panel have
    assert _handle_colors(axes[-1]) == pytest.approx(
        np.asarray([to_rgba(line.get_color()) for line in axes[0].lines])
    )


def test_the_data_handle_of_a_row_follows_the_colour_rule() -> None:
    profiles = _profiles(1, n=4)
    times = (60.0,)
    data = _data(label=["p0", "p0", "p0"])
    plain = [_column("A", profiles=profiles, times=times)]
    coloured = [_column("A", profiles=profiles, times=times, color="#ff0000")]
    ax = data_axes(zone_profiles(coloured, profiles, times=times, data=data))[0]
    assert _handle_colors(ax) == pytest.approx(np.vstack([to_rgba("0.3")] * 2))
    assert [h.get_marker() for h in _handles(ax)] == ["None", "s"]
    ax = data_axes(zone_profiles(plain, profiles, times=times, data=data))[0]
    expected = colormaps["viridis"](np.linspace(0.1, 0.75, 1))[0]
    assert _handle_colors(ax) == pytest.approx(np.vstack([expected, expected]))


def test_a_panel_of_a_grid_draws_no_legend_of_its_own() -> None:
    """The grid draws the handles of a row itself, so the artists carry no label."""
    profiles = _profiles(1, n=4)
    times = (0.0, 60.0)
    data = _data(label=["p0", "p0", "p0"], time=[60.0, 60.0, 60.0])
    columns = [
        _column("A", profiles=profiles, times=times, color="#ff0000"),
        _column("B", profiles=profiles, times=times, color="#0000ff"),
    ]
    first, _ = data_axes(zone_profiles(columns, profiles, times=times, data=data))
    assert first.get_legend_handles_labels() == ([], [])
    frame = data_axes(zone_profiles(_result(n=4), profiles, times=times))[0]
    assert frame.get_legend_handles_labels()[1] == ["Sim t = 0 min", "Sim t = 60 min"]


BAD = "#f0f0f0"
TRIAD_COLORS = ("#1f77b4", "#d62728", "#2ca02c")


def _mappable() -> ScalarMappable:
    """A viridis scale of the values 1 to N, the grey of a zone without data."""
    return ScalarMappable(
        Normalize(1.0, float(N)), colormaps["viridis"].with_extremes(bad=BAD)
    )


def _colour(mappable: ScalarMappable, value: float) -> np.ndarray:
    """The rgba of one value on a scale (`to_rgba` is annotated for arrays)."""
    return np.asarray(mappable.to_rgba(np.array([value])))[0]


def _values() -> list[float]:
    """A distinct value per zone, the gradient the maps draw."""
    return [float(k) for k in range(1, N + 1)]


def _axes() -> Axes:
    """One axes of a figure of its own, the panel a caller hands to the renderer."""
    return Figure().subplots()


def _polygons(ax: Axes) -> list[Polygon]:
    """The patches of the bands of a panel."""
    return [p for p in ax.patches if isinstance(p, Polygon)]


def _zone_face(ax: Axes, k: int) -> tuple[float, ...]:
    """The one face colour of the patches of zone `k`."""
    colours = {
        to_rgba(p.get_facecolor()) for p in _polygons(ax) if p.get_gid() == ids.zone(k)
    }
    assert len(colours) == 1
    return colours.pop()


def _face_at(ax: Axes, point: tuple[float, float]) -> tuple[float, ...]:
    """The face colour of the patch which holds a point (data coordinates)."""
    colours = {
        to_rgba(p.get_facecolor())
        for p in _polygons(ax)
        if p.get_path().contains_point(point, radius=0.0)
    }
    assert len(colours) == 1, f"{point} lies in {len(colours)} colours"
    return colours.pop()


def _stroke(text: Text) -> dict[str, Any]:
    """The arguments of the one stroke behind a text (`withStroke` keeps them private)."""
    (effect,) = text.get_path_effects()
    assert isinstance(effect, withStroke)
    arguments = vars(effect)["_gc"]
    assert isinstance(arguments, dict)
    return arguments


def _xy(line: Line2D) -> np.ndarray:
    """The points of a line as an array of shape (n, 2)."""
    return np.column_stack([_x(line), _y(line)])


def test_draw_lobule_colours_a_band_per_zone() -> None:
    values = _values()
    mappable = _mappable()
    ax = _axes()
    draw_lobule(ax, values, mappable=mappable)
    assert len(_polygons(ax)) == 12 * N
    for k in (1, 4, N):
        assert _zone_face(ax, k) == pytest.approx(_colour(mappable, values[k - 1]))
    # the band of a zone lies where the zone lies: 1 at a corner, N at the centre
    assert _face_at(ax, (0.98, 0.0)) == pytest.approx(_colour(mappable, values[0]))
    assert _face_at(ax, (0.0, 0.0)) == pytest.approx(_colour(mappable, values[-1]))
    assert ax.get_xlim() == pytest.approx((-1.02, 1.02))
    assert ax.get_ylim() == pytest.approx((-0.9, 0.9))
    assert ax.get_aspect() == 1.0
    assert list(ax.get_xticks()) == []
    assert list(ax.get_yticks()) == []
    assert not any(spine.get_visible() for spine in ax.spines.values())


def test_draw_lobule_draws_a_zone_without_data_in_the_missing_colour() -> None:
    """A plain colormap: the renderer takes the colour, the caller sets no `bad`."""
    values = _values()
    for k in (4, 5, 6, 7):
        values[k - 1] = float("nan")
    mappable = ScalarMappable(Normalize(1.0, float(N)), colormaps["viridis"])
    assert np.asarray(colormaps["viridis"].get_bad()) == pytest.approx(to_rgba("none"))
    ax = _axes()
    draw_lobule(ax, values, mappable=mappable)
    for k in (4, 7):
        assert _zone_face(ax, k) == pytest.approx(to_rgba(MISSING))
        # neither transparent, which is a hole in the hexagon, nor a value
        assert _zone_face(ax, k)[3] == pytest.approx(1.0)
        assert _zone_face(ax, k) != pytest.approx(_colour(mappable, 1.0))
    # the zones with data keep the colour of their own value, not of a neighbour
    for k in (1, 3, 8, N):
        assert _zone_face(ax, k) == pytest.approx(_colour(mappable, float(k)))


def test_the_colour_of_a_zone_without_data_is_the_callers_to_choose() -> None:
    values = _values()
    values[2] = float("nan")
    ax = _axes()
    draw_lobule(ax, values, mappable=_mappable(), missing="#ff00ff")
    assert _zone_face(ax, 3) == pytest.approx(to_rgba("#ff00ff"))
    assert _zone_face(ax, 4) == pytest.approx(_colour(_mappable(), 4.0))


def test_draw_lobule_draws_the_thirds_only_when_asked() -> None:
    plain = _axes()
    draw_lobule(plain, _values(), mappable=_mappable())
    assert len(plain.lines) == 0
    ax = _axes()
    draw_lobule(ax, _values(), mappable=_mappable(), thirds=True)
    assert len(ax.lines) == 2
    for line, position in zip(ax.lines, (1.0 / 3.0, 2.0 / 3.0), strict=True):
        assert _xy(line) == pytest.approx(lobule_isoline(position))
        assert line.get_linestyle() == ":"
        assert to_rgba(line.get_color()) == to_rgba("0.25")
        assert line.get_linewidth() == pytest.approx(0.7)
    # the ring of 1/3 lies wholly outside the ring of 2/3, in that order
    radii = [np.hypot(*_xy(line).T) for line in ax.lines]
    assert radii[0].min() == pytest.approx(2.0 / 3.0)  # towards a portal triad
    assert radii[1].min() == pytest.approx(1.0 / 3.0)
    assert radii[0].min() > radii[1].max()


def test_draw_lobule_draws_the_triads_only_when_asked() -> None:
    plain = _axes()
    draw_lobule(plain, _values(), mappable=_mappable())
    assert len(plain.lines) == 0
    ax = _axes()
    draw_lobule(ax, _values(), mappable=_mappable(), anatomy=True)
    outline, *triads = ax.lines
    assert _xy(outline) == pytest.approx(lobule_isoline(0.0))
    assert to_rgba(outline.get_color()) == to_rgba("0.15")
    assert outline.get_linewidth() == pytest.approx(1.1)
    assert [to_rgba(line.get_color()) for line in triads] == [
        to_rgba(colour) for colour in TRIAD_COLORS
    ]
    points = np.vstack([_xy(line) for line in triads])
    assert len(points) == 18  # three dots at every one of the six corners
    angles = np.deg2rad(60.0 * np.arange(6))
    corners = np.column_stack([np.cos(angles), np.sin(angles)])
    distance = np.linalg.norm(points[:, None, :] - corners[None, :, :], axis=2)
    assert distance.min(axis=1) == pytest.approx(0.045, abs=5e-4)
    assert sorted(np.bincount(distance.argmin(axis=1), minlength=6)) == [3] * 6
    for line in triads:
        assert line.get_marker() == "o"
        assert line.get_linestyle() == "None"
        assert line.get_markersize() == pytest.approx(2.6)
        assert to_rgba(line.get_markeredgecolor()) == to_rgba("0.15")
        assert line.get_markeredgewidth() == pytest.approx(0.3)


def test_draw_lobule_names_the_vessels_only_when_asked() -> None:
    plain = _axes()
    draw_lobule(plain, _values(), mappable=_mappable(), anatomy=True)
    assert [t.get_text() for t in plain.texts] == []
    ax = _axes()
    draw_lobule(ax, _values(), mappable=_mappable(), anatomy=True, labels=True)
    assert [t.get_text() for t in ax.texts] == ["CV", "PT"]
    centre, triad = ax.texts
    assert centre.get_position() == pytest.approx((0.0, 0.0))
    assert centre.get_fontsize() == pytest.approx(7.0)
    assert triad.get_fontsize() == pytest.approx(7.0)
    # below the corner at 300 degrees, under the panel, clear of its three dots
    x, y = triad.get_position()
    assert (x, y) == pytest.approx((0.5, -0.95))
    assert triad.get_verticalalignment() == "top"
    assert y < -np.cos(np.pi / 6.0) - 0.045 - 0.02
    # both read on a dark map and on a light one: a white stroke behind them
    for text in (centre, triad):
        stroke = _stroke(text)
        assert to_rgba(stroke["foreground"]) == to_rgba("white")
        assert stroke["linewidth"] == pytest.approx(2.0)


def test_lobule_map_is_unchanged_by_the_renderer() -> None:
    """The map draws its panel with the renderer: the patches, the colours, the title."""
    fig = lobule_map(_result(), "drg", time=60, zones=N)
    (ax,) = data_axes(fig)
    assert len(_polygons(ax)) == 12 * N
    # the cell layer is the fifth: the value of zone k at the time 60 is 2.5 k
    assert _zone_face(ax, 1) == pytest.approx(colormaps["viridis"](0.0))
    assert _zone_face(ax, N) == pytest.approx(colormaps["viridis"](1.0))
    assert _face_at(ax, (0.98, 0.0)) == pytest.approx(_zone_face(ax, 1))
    assert _face_at(ax, (0.0, 0.0)) == pytest.approx(_zone_face(ax, N))
    assert ax.get_title() == "drg cell at t = 60 min"
    assert "portal triads (0, periportal)" in ax.get_xlabel()
    assert len([a for a in fig.axes if a.get_gid() == COLORBAR]) == 1
    # and it draws the anatomy: the outline and the three dots of every triad
    assert len(ax.lines) == 4
    assert len(ax.texts) == 0


def test_the_colorbar_of_several_axes_is_vector(tmp_path: Path) -> None:
    """A model which composes its own grid of maps needs the guarantee of the module."""
    fig = Figure(figsize=(6.0, 3.0), layout="constrained")
    left, right = fig.subplots(1, 2)
    mappable = _mappable()
    draw_lobule(left, _values(), mappable=mappable, anatomy=True)
    draw_lobule(right, _values()[::-1], mappable=mappable, anatomy=True)
    bar = colorbar(
        fig,
        mappable,
        [left, right],
        label="Drug [mM]",
        location="right",
        shrink=0.8,
        aspect=30,
    )
    assert bar.solids is not None
    assert not bar.solids.get_rasterized()
    assert bar.ax.get_gid() == COLORBAR
    assert bar.ax.get_ylabel() == "Drug [mM]"
    assert data_axes(fig) == [left, right]
    fig.savefig(tmp_path / "grid.svg")
    assert "<image" not in (tmp_path / "grid.svg").read_text()
