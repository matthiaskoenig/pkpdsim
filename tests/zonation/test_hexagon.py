"""The lobule hexagon and the tube mesh, measured from the artists and not from pixels."""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from itertools import pairwise
from pathlib import Path

from matplotlib import colormaps
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch

from pkpdsim.species import Species
from pkpdsim.zonation import ids, observers
from pkpdsim.zonation.figures import COLORBAR, data_axes, lobule_map, space_time
from pkpdsim.zonation.geometry import Layer, geometry
from pkpdsim.zonation.profiles import lobular_area, zone_areas
from pkpdsim.zonation.render import lobule_bands, lobule_isoline, tube_mesh

N = 10
HEXAGON_AREA = 3.0 * np.sqrt(3.0) / 2.0
APOTHEM = np.cos(np.pi / 6.0)


def _result(layer: Layer = Layer.CELL, n: int = N) -> pd.DataFrame:
    """A synthetic result of `drg`, the value of a zone `k` is `k` at the time 60."""
    t = np.linspace(0.0, 120.0, 25)
    cols: dict[str, np.ndarray] = {"time": t}
    for k, sid in enumerate(observers.profile_ids("drg", layer, n), start=1):
        cols[sid] = k * t / 60.0
    cols[observers.outflow_concentration("drg")] = t / 120.0
    return pd.DataFrame(cols)


def _patches(fig: Figure) -> list[PathPatch]:
    """The compound patches of the hexagon, one per zone, from the artists."""
    (ax,) = data_axes(fig)
    return [p for p in ax.patches if isinstance(p, PathPatch)]


def _colour_at(fig: Figure, point: tuple[float, float]) -> tuple[float, ...]:
    """The face colour of the patch which holds a point (data coordinates)."""
    colours = {
        to_rgba(p.get_facecolor())
        for p in _patches(fig)
        if p.get_path().contains_point(point, radius=0.0)
    }
    assert len(colours) == 1, f"{point} lies in {len(colours)} colours"
    return colours.pop()


def _rim_point() -> tuple[float, float]:
    """A point at a corner of the hexagon (the portal triad at angle 0)."""
    return (0.98, 0.0)


def _zone_colour(fig: Figure, k: int) -> tuple[float, ...]:
    """The face colour of the patches of zone `k`."""
    colours = {
        to_rgba(p.get_facecolor()) for p in _patches(fig) if p.get_gid() == ids.zone(k)
    }
    assert len(colours) == 1
    return colours.pop()


def _area(patch: PathPatch) -> float:
    """The area of a compound patch by the shoelace formula over its wedges."""
    total = 0.0
    for wedge in patch.get_path().to_polygons(closed_only=True):
        v = np.asarray(wedge)
        x, y = v[:, 0], v[:, 1]
        total += 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))
    return total


def _areas_per_zone(fig: Figure) -> list[float]:
    """The share of the area of the hexagon of every zone."""
    total = dict.fromkeys(range(1, N + 1), 0.0)
    for p in _patches(fig):
        total[int(str(p.get_gid())[1:])] += _area(p)
    s = sum(total.values())
    return [a / s for a in total.values()]


def test_the_hexagon_maps_zone_one_to_the_rim_and_zone_n_to_the_centre() -> None:
    fig = lobule_map(_result(), "drg", time=60, zones=N)
    assert _colour_at(fig, _rim_point()) == pytest.approx(_zone_colour(fig, 1))
    assert _colour_at(fig, (0.0, 0.0)) == pytest.approx(_zone_colour(fig, N))
    assert _zone_colour(fig, 1) != pytest.approx(_zone_colour(fig, N))


def test_the_colours_are_the_scale_of_the_values() -> None:
    """The value of zone `k` is `k`, the scale runs from the smallest value to the largest."""
    fig = lobule_map(_result(), "drg", time=60, zones=N)
    cmap = colormaps["viridis"]
    for k in (1, 5, N):
        assert _zone_colour(fig, k) == pytest.approx(cmap((k - 1) / (N - 1)))


def _band_areas(**sampling: int) -> list[float]:
    """The share of the area of the hexagon of every zone, from the polygons."""
    bands = lobule_bands(N, **sampling)
    return [
        sum(_shoelace(np.asarray(wedge)) for wedge in wedges) / HEXAGON_AREA
        for wedges in bands
    ]


def test_the_hexagon_area_per_zone_is_the_exact_zone_area() -> None:
    """A zone is the interval `[(k - 1) / n, k / n]`, its area `zone_areas(n)`.

    The polygons are the chords of the isolines, so the accurate sampling is
    asked for here; the areas of the engine come from `lobular_area` and never
    from these polygons, and the default of a drawing is the coarser one below.
    """
    areas = _band_areas(points=256)
    assert areas == pytest.approx(list(zone_areas(N)), rel=1e-5)
    assert sum(areas) == pytest.approx(1.0, rel=1e-9)


def test_the_sampling_of_a_drawing_keeps_the_areas_within_a_thousandth() -> None:
    """What the figures draw is coarse, and a zone is still its area."""
    assert _band_areas() == pytest.approx(list(zone_areas(N)), rel=1e-3)
    drawn = _areas_per_zone(lobule_map(_result(), "drg", time=60, zones=N))
    assert drawn == pytest.approx(list(zone_areas(N)), rel=1e-3)


def test_the_bands_refuse_a_wedge_of_one_point() -> None:
    with pytest.raises(ValueError, match="two points"):
        lobule_bands(4, points=1)


def _wedge_arc(position: float, points: int) -> np.ndarray:
    """The first wedge of the isoline of a position, the arc two bands share."""
    ring = lobule_isoline(position, points=points)
    return ring[: (len(ring) - 1) // 12 + 1]


def test_a_band_and_the_band_inside_it_share_their_isoline() -> None:
    """No seam: the inner arc of a band is the outer arc of the band inside it."""
    n, points = 4, 8
    bands = lobule_bands(n, points=points)
    for k, (outside, inside) in enumerate(pairwise(bands), start=1):
        arc = _wedge_arc(k / n, points)
        assert outside[0][-len(arc) :] == pytest.approx(arc[::-1])
        assert inside[0][: len(arc)] == pytest.approx(arc)


def test_every_point_of_the_hexagon_lies_in_exactly_one_patch() -> None:
    fig = lobule_map(_result(), "drg", time=60, zones=N)
    rng = np.random.default_rng(0)
    checked = 0
    for x, y in rng.uniform(-1.0, 1.0, (1500, 2)):
        if abs(y) >= np.sqrt(3.0) / 2.0 - 1e-9:
            continue
        if abs(y) * 0.5 + abs(x) * np.sqrt(3.0) / 2.0 >= np.sqrt(3.0) / 2.0 - 1e-9:
            continue
        held = [p for p in _patches(fig) if p.get_path().contains_point((x, y))]
        assert len(held) == 1
        checked += 1
    assert checked > 500


def test_the_hexagon_is_the_unit_hexagon_with_its_corners_on_the_axis() -> None:
    fig = lobule_map(_result(), "drg", time=60, zones=N)
    v = np.vstack([p.get_path().vertices for p in _patches(fig)])
    assert v[:, 0].max() == pytest.approx(1.0)
    assert np.hypot(v[:, 0], v[:, 1]).max() == pytest.approx(1.0)
    total = sum(_area(p) for p in _patches(fig))
    assert total == pytest.approx(3.0 * np.sqrt(3.0) / 2.0, rel=1e-6)


def test_the_map_names_its_ends_and_the_unit_and_has_one_colorbar() -> None:
    fig = lobule_map(_result(), "drg", time=60, zones=N, unit="uM", name="Drug")
    (ax,) = data_axes(fig)
    text = ax.get_xlabel() + ax.get_title()
    assert "periportal" in text
    assert "pericentral" in text
    assert "t = 60 min" in text
    bars = [a for a in fig.axes if a.get_gid() == COLORBAR]
    assert len(bars) == 1
    assert "[uM]" in bars[0].get_ylabel()
    assert "Drug" in bars[0].get_ylabel()


def test_the_map_of_a_layer_reads_that_layer() -> None:
    result = _result(Layer.PLASMA)
    fig = lobule_map(result, "drg", time=60, zones=N, layer=Layer.PLASMA)
    assert len(_patches(fig)) > 0
    with pytest.raises(KeyError, match="column"):
        lobule_map(result, "drg", time=60, zones=N)


def test_the_map_is_vector_and_uses_no_pyplot(tmp_path: Path) -> None:
    import matplotlib.pyplot as plt

    before = plt.get_fignums()
    fig = lobule_map(_result(), "drg", time=60, zones=N)
    assert isinstance(fig, Figure)
    assert plt.get_fignums() == before
    fig.savefig(tmp_path / "lobule.svg")
    text = (tmp_path / "lobule.svg").read_text()
    assert "<image" not in text
    assert text.count("<path") >= N


def test_one_zone_is_one_colour() -> None:
    fig = lobule_map(_result(n=1), "drg", time=60, zones=1)
    assert {to_rgba(p.get_facecolor()) for p in _patches(fig)} == {
        _zone_colour_one(fig)
    }


def _zone_colour_one(fig: Figure) -> tuple[float, ...]:
    return _zone_colour(fig, 1)


def _shoelace(ring: np.ndarray) -> float:
    """The area enclosed by a closed polyline, by the shoelace formula."""
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _angles(ring: np.ndarray) -> np.ndarray:
    """The direction of every point of a ring around the centre (degrees)."""
    return np.rad2deg(np.arctan2(ring[:, 1], ring[:, 0])) % 360.0


def _radius_at(ring: np.ndarray, degrees: float) -> float:
    """The radius of an isoline in a direction, interpolated along its points."""
    angle, radius = _angles(ring), np.hypot(ring[:, 0], ring[:, 1])
    order = np.argsort(angle)
    return float(np.interp(degrees, angle[order], radius[order]))


def test_the_isoline_of_zero_is_the_hexagon() -> None:
    ring = lobule_isoline(0.0)
    radius = np.hypot(ring[:, 0], ring[:, 1])
    # radius 1 towards a corner, cos(pi / 6) towards the middle of an edge
    assert radius.max() == pytest.approx(1.0)
    assert radius.min() == pytest.approx(APOTHEM)
    assert _radius_at(ring, 0.0) == pytest.approx(1.0)
    assert _radius_at(ring, 30.0) == pytest.approx(APOTHEM)
    corners = np.unique(np.round(_angles(ring)[np.isclose(radius, 1.0, atol=1e-9)]))
    assert sorted(corners % 360.0) == [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
    assert _shoelace(ring) == pytest.approx(HEXAGON_AREA, rel=1e-9)


def test_the_isoline_of_one_is_the_centre() -> None:
    ring = lobule_isoline(1.0)
    assert np.hypot(ring[:, 0], ring[:, 1]) == pytest.approx(0.0)
    assert _shoelace(ring) == pytest.approx(0.0)
    # the centre is the same point in every direction: two per wedge, not more
    assert len(ring) == len(lobule_isoline(1.0, points=256)) == 13


def test_the_isoline_shrinks_with_the_position() -> None:
    """The radius in a direction never grows with the position, and the area shrinks."""
    rings = [lobule_isoline(x) for x in (0.0, 0.25, 0.5, 0.75, 0.9)]
    for degrees in (0.0, 5.0, 30.0, 60.0, 95.0, 200.0, 330.0):
        radii = [_radius_at(ring, degrees) for ring in rings]
        assert all(b <= a + 1e-9 for a, b in pairwise(radii)), (degrees, radii)
    # towards a portal triad it decreases strictly, and so does the area inside
    for degrees in (0.0, 60.0, 180.0):
        radii = [_radius_at(ring, degrees) for ring in rings]
        assert all(b < a for a, b in pairwise(radii)), (degrees, radii)
    areas = [_shoelace(ring) for ring in rings]
    assert all(b < a for a, b in pairwise(areas)), areas


@pytest.mark.parametrize("x", [0.0, 0.2, 0.5, 0.8, 1.0])
def test_the_radius_towards_a_triad_is_the_distance_to_the_central_vein(
    x: float,
) -> None:
    """On the line triad to centre the position is the distance from the triad."""
    ring = lobule_isoline(x)
    for degrees in (0.0, 120.0, 240.0):
        assert _radius_at(ring, degrees) == pytest.approx(1.0 - x)


@pytest.mark.parametrize("x", [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0])
def test_the_area_inside_an_isoline_is_the_pericentral_area(x: float) -> None:
    """The shoelace area of the isoline of `x` is `1 - lobular_area(x)` of the hexagon."""
    share = _shoelace(lobule_isoline(x)) / HEXAGON_AREA
    assert share == pytest.approx(1.0 - lobular_area(x), abs=1e-3)


def test_the_isoline_is_closed_and_ordered() -> None:
    ring = lobule_isoline(0.4)
    assert ring.ndim == 2
    assert ring.shape[1] == 2
    assert ring[-1] == pytest.approx(ring[0])
    angle = np.unwrap(np.arctan2(ring[:, 1], ring[:, 0]))
    assert np.all(np.diff(angle) >= -1e-12)
    assert angle[-1] - angle[0] == pytest.approx(2.0 * np.pi)


def test_the_isoline_samples_its_arc_as_asked_and_refuses_one_point() -> None:
    assert len(lobule_isoline(0.5, points=8)) < len(lobule_isoline(0.5, points=64))
    assert _shoelace(lobule_isoline(0.5, points=16)) == pytest.approx(
        _shoelace(lobule_isoline(0.5, points=256)), rel=1e-3
    )
    # the isoline of 0 is the hexagon, twelve straight pieces and no arc at all
    assert len(lobule_isoline(0.0, points=64)) == 13
    with pytest.raises(ValueError, match="two points"):
        lobule_isoline(0.5, points=1)


def test_the_isoline_refuses_a_position_outside_the_lobule() -> None:
    with pytest.raises(ValueError, match="outside"):
        lobule_isoline(1.5)


def test_the_tube_mesh_has_a_segment_per_zone() -> None:
    _, cells = tube_mesh(geometry(Species.RAT), N)
    assert len(cells) == N


def test_the_tube_mesh_is_the_tube_of_the_geometry() -> None:
    geo = geometry(Species.RAT)
    points, cells = tube_mesh(geo, N, sides=12)
    assert points.shape == ((N + 1) * 12, 3)
    assert cells.shape == (N, 24)
    assert cells.max() == len(points) - 1
    assert points[:, 0].min() == pytest.approx(0.0)
    assert points[:, 0].max() == pytest.approx(geo.L_sin.value)
    assert np.hypot(points[:, 1], points[:, 2]) == pytest.approx(geo.r_sin.value)
    # the segments follow each other from the periportal end
    x = points[cells][:, :, 0]
    assert x[:, :12].max(axis=1) == pytest.approx(np.arange(N) * geo.L_sin.value / N)
    assert x[:, 12:].max(axis=1) == pytest.approx(
        np.arange(1, N + 1) * geo.L_sin.value / N
    )


def test_the_tube_mesh_refuses_no_zone() -> None:
    with pytest.raises(ValueError, match="zone"):
        tube_mesh(geometry(Species.RAT), 0)


def test_one_zone_has_one_peclet_number_of_two_and_constant_profiles() -> None:
    """The claims of the documentation page for `n_zones=1`."""
    from pkpdsim.pbpk.physiology import reference_flow
    from pkpdsim.zonation import effective_peclet
    from pkpdsim.zonation.profiles import Exponential, Linear, Sigmoid
    from pkpdsim.zonation.settings import ZonationSettings

    settings = ZonationSettings(n_zones=1, reference_flow=reference_flow(Species.RAT))
    assert effective_peclet(settings, 0.0) == pytest.approx(2.0)
    for profile in (Linear(1.0, 4.0), Sigmoid(0.5, 6.0), Exponential(1.0, 3.0)):
        assert profile.values(1) == pytest.approx((1.0,))


def test_the_caption_is_whole_and_clear_of_the_colorbar() -> None:
    fig = lobule_map(_result(), "drg", time=60, zones=N)
    FigureCanvasAgg(fig).draw()
    (ax,) = data_axes(fig)
    renderer = fig.canvas.get_renderer()  # ty: ignore[unresolved-attribute]
    caption = ax.xaxis.label.get_window_extent(renderer)
    bar = next(a for a in fig.axes if a.get_gid() == COLORBAR)
    bar_box = bar.get_tightbbox(renderer)
    assert ax.get_xlabel().count("\n") == 1
    assert "portal triads (0, periportal)" in ax.get_xlabel()
    assert "central vein (1, pericentral)" in ax.get_xlabel()
    assert fig.bbox.x0 <= caption.x0
    assert caption.x1 <= fig.bbox.x1
    assert caption.y0 >= fig.bbox.y0
    assert bar_box is not None
    assert not caption.overlaps(bar_box)


def test_the_colour_scale_spans_the_range_of_the_values() -> None:
    """A small gradient far from zero must not be one shade; a line axis keeps its zero."""
    result = _result()
    for k, sid in enumerate(observers.profile_ids("drg", Layer.CELL, N)):
        result[sid] = 0.078 + 0.012 * k / (N - 1)
    fig = lobule_map(result, "drg", time=60, zones=N)
    bar = next(a for a in fig.axes if a.get_gid() == COLORBAR)
    assert bar.get_ylim() == pytest.approx((0.078, 0.090))
    heat = space_time(result, "drg", layers=(Layer.CELL,))
    assert data_axes(heat)[0].collections[0].get_clim() == pytest.approx((0.078, 0.090))
    assert min(bar.get_yticks()) >= 0.0


def test_vmin_and_vmax_override_the_data_range_and_default_to_it() -> None:
    """A caller clips a scale which a spike dominates; leaving them out keeps the data range."""
    result = _result()
    values = (1.0, 10.0)  # zone 1 and zone 10 at time 60
    fig = lobule_map(result, "drg", time=60, zones=N)
    bar = next(a for a in fig.axes if a.get_gid() == COLORBAR)
    assert bar.get_ylim() == pytest.approx(values)
    clipped = lobule_map(result, "drg", time=60, zones=N, vmin=2.0, vmax=6.0)
    bar = next(a for a in clipped.axes if a.get_gid() == COLORBAR)
    assert bar.get_ylim() == pytest.approx((2.0, 6.0))
    half = lobule_map(result, "drg", time=60, zones=N, vmax=5.0)
    bar = next(a for a in half.axes if a.get_gid() == COLORBAR)
    assert bar.get_ylim() == pytest.approx((1.0, 5.0))
    heat = space_time(result, "drg", layers=(Layer.CELL,), vmin=0.0, vmax=3.0)
    assert data_axes(heat)[0].collections[0].get_clim() == (0.0, 3.0)
    default = space_time(result, "drg", layers=(Layer.CELL,))
    assert data_axes(default)[0].collections[0].get_clim()[1] == pytest.approx(20.0)
    with pytest.raises(ValueError, match="no range"):
        lobule_map(result, "drg", time=60, zones=N, vmin=5.0, vmax=5.0)


def test_a_vmin_at_a_constant_raises_where_the_default_has_a_fallback() -> None:
    result = _result()
    for sid in observers.profile_ids("drg", Layer.CELL, N):
        result[sid] = 0.25
    fig = lobule_map(result, "drg", time=60, zones=N)
    bar = next(a for a in fig.axes if a.get_gid() == COLORBAR)
    assert bar.get_ylim() == pytest.approx((0.0, 0.25))
    with pytest.raises(ValueError, match="no range"):
        lobule_map(result, "drg", time=60, zones=N, vmin=0.25)
