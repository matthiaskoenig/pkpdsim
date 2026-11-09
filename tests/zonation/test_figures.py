"""The zone figures: over a synthetic result, and over the result of a built chain."""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from pathlib import Path

from matplotlib.axes import Axes
from matplotlib.collections import LineCollection, QuadMesh
from matplotlib.figure import Figure
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
    assert (
        data_axes(layer_profile(_result(), "drg", time=0, unit="uM"))[0].get_ylabel()
        == "drg [uM]"
    )
    assert data_axes(layer_profile(_result(), "drg", time=0))[0].get_ylabel() == (
        "drg [mM]"
    )


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
    assert out.get_xlabel() == "Time"


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
