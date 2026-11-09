"""The zone figures: over a synthetic result, and over the result of a built chain."""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from pathlib import Path

from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import LineCollection, QuadMesh
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from pbpk.conftest import individual
from sbmlutils.factory import create_model

from pkpdsim.zonation import ids, observers
from pkpdsim.zonation.figures import (
    ZoneProfile,
    data_axes,
    layer_profile,
    result_frame,
    space_time,
    zone_profiles,
)
from pkpdsim.zonation.geometry import Layer
from pkpdsim.zonation.observers import add_observers, observer_objects

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
    """Lay a figure out: a figure of the library has no canvas which draws."""
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


def test_the_y_labels_are_aligned() -> None:
    fig = space_time(_result(), "drg", layers=(Layer.PLASMA, Layer.CELL))
    _drawn(fig)
    heat, other, out = data_axes(fig)
    xs = {round(ax.yaxis.label.get_position()[0], 3) for ax in (heat, other, out)}
    assert len(xs) == 1


def test_the_outflow_panel_is_lower_than_a_heatmap_and_the_name_is_chosen() -> None:
    fig = space_time(_result(), "drg", layers=(Layer.PLASMA,), name="Drug")
    heat, out = data_axes(fig)
    spec = heat.get_subplotspec()
    assert spec is not None
    assert spec.get_gridspec().get_height_ratios() == [1.0, 0.6]
    _drawn(fig)
    assert out.get_position().height < 0.7 * heat.get_position().height
    assert fig.axes[-1].get_ylabel() == "Drug plasma [mM]"
    assert out.get_ylabel() == "Drug outflow [mM]"
    default = space_time(_result(), "drg", layers=(Layer.PLASMA,))
    assert default.axes[-1].get_ylabel() == "drg plasma [mM]"
