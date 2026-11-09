"""The geometry of the zonated liver for renderers: the tube and the lobule hexagon.

Both are plain numpy arrays and know no figure library. `lobule_bands` is the
polygons of the idealized hexagonal lobule which `pkpdsim.zonation.figures.lobule_map`
draws; `tube_mesh` is the sinusoid as points and cells.

The tube is the geometry of the engine (`SinusoidGeometry`): a cylinder of the
sinusoid radius and length along the x axis, from the portal end (x = 0) to
the central end (x = `L_sin`), with one segment per zone. The arrays are what a
renderer of unstructured meshes takes; the engine neither depends on one nor
renders anything: a script which draws (the videos of gdeobdtpa) installs
pyvista for its run, builds a grid of these arrays and colours the segments by
the zone values, `values[k]` for the segment `k`.
"""

import math
from itertools import pairwise

import numpy as np

from pkpdsim.zonation.geometry import SinusoidGeometry
from pkpdsim.zonation.profiles import lobular_radius

__all__ = ["lobule_bands", "tube_mesh"]


_WEDGE = math.pi / 6.0
"""The angle of one twelfth of the hexagon, the wedge of one portal triad."""
_ARC_POINTS = 256
"""The points of the isoline of a band in a wedge; the error of the area is below 1e-4 of a zone."""


def _zone_edges(n: int) -> list[float]:
    """Get the positions which bound the `n` bands of the hexagon, `k / n`.

    A zone is an interval of the lobular position, `[(k - 1) / n, k / n]`:
    zone 1 starts at 0 (the portal triads, the rim) and zone `n` ends at 1
    (the central vein, the centre). The area of a band is then the difference
    of the areas of its ends, `zone_areas(n)` exactly; it is not the weight at
    the centre of the zone, which overstates the zone of the peak by a few
    percent (see `pkpdsim.zonation.profiles`).
    """
    return [k / n for k in range(n + 1)]


def _band(low: float, high: float) -> np.ndarray:
    """Get the polygon of the band `[low, high]` of a position in the wedge of a triad.

    The outer arc is the isoline of `low`, the inner arc that of `high`, both
    from the direction of the triad (theta 0) to that of the edge midpoint.
    The directions in which an isoline leaves the hexagon are sampled too, so
    the corner where it meets the edge is exact.
    """
    theta = np.linspace(0.0, _WEDGE, _ARC_POINTS)
    # the direction in which an isoline leaves through the edge of the hexagon
    # (`sin(theta) = sqrt(3) / 2 * x / (1 - x)`) is a corner of the band
    kinks = [
        math.asin(min(1.0, math.sqrt(3.0) / 2.0 * x / (1.0 - x)))
        for x in (low, high)
        if 0.0 < x < 1.0
    ]
    theta = np.unique(np.concatenate([theta, np.clip(kinks, 0.0, _WEDGE)]))
    outer = [lobular_radius(low, t) for t in theta]
    inner = [lobular_radius(high, t) for t in theta]
    return np.vstack(
        [
            np.column_stack(
                [np.array(outer) * np.cos(theta), np.array(outer) * np.sin(theta)]
            ),
            np.column_stack(
                [
                    np.array(inner[::-1]) * np.cos(theta[::-1]),
                    np.array(inner[::-1]) * np.sin(theta[::-1]),
                ]
            ),
        ]
    )


def _twelve(wedge: np.ndarray) -> list[np.ndarray]:
    """Get the twelve copies of a polygon of the first wedge over the hexagon.

    The hexagon has the corners on the angles 0, 60, ... degrees; the wedge
    `[0, 30]` degrees is mirrored about the 30 degree line and rotated, the
    copies of the odd wedges are the mirror images of the even ones.
    """
    copies: list[np.ndarray] = []
    for j in range(6):
        c, s = math.cos(j * math.pi / 3.0), math.sin(j * math.pi / 3.0)
        rotation = np.array([[c, -s], [s, c]])
        copies.append(wedge @ rotation.T)
        # the mirror image about the line at 30 degrees: (x, y) -> (cos 60 x + sin 60 y, sin 60 x - cos 60 y)
        mirror = np.array(
            [
                [math.cos(math.pi / 3.0), math.sin(math.pi / 3.0)],
                [math.sin(math.pi / 3.0), -math.cos(math.pi / 3.0)],
            ]
        )
        copies.append(wedge @ mirror.T @ rotation.T)
    return copies


def lobule_bands(n: int) -> list[list[np.ndarray]]:
    """Get the polygons of the zones of the idealized hexagonal lobule.

    The hexagon has the circumradius 1, the central vein (zone `n`) at the
    origin and the portal triads (zone 1) at the corners on the angles 0, 60,
    ... degrees. The band of zone `k` is `[(k - 1) / n, k / n]` of the
    normalized periportal to pericentral distance, in twelve wedges.

    Args:
        n: the number of zones.

    Returns:
        For every zone from periportal, its twelve polygons, arrays of shape
        `(points, 2)`; the areas of the zones are `zone_areas(n)` to the chord
        error of the isolines (8.7e-5 relative at 10 zones).

    Raises:
        ValueError: for less than one zone.
    """
    if n < 1:
        raise ValueError(f"A lobule has at least one zone, not {n}.")
    return [_twelve(_band(low, high)) for low, high in pairwise(_zone_edges(n))]


def tube_mesh(
    geo: SinusoidGeometry, n_zones: int, *, sides: int = 24
) -> tuple[np.ndarray, np.ndarray]:
    """Get the mesh of the sinusoid, a segment per zone.

    Args:
        geo: the geometry; its sinusoid radius and length (m) give the tube.
        n_zones: the number of zones, the number of segments.
        sides: the number of points of a ring (at least 3).

    Returns:
        The points, an array of shape `((n_zones + 1) * sides, 3)` in m, the
        rings from the periportal end (x = 0) to the pericentral end, and the
        cells, an integer array of shape `(n_zones, 2 * sides)`: the points of
        segment `k` (zone `k + 1`) are the `sides` points of the ring `k`
        followed by the `sides` points of the ring `k + 1`, each in ring
        order, counterclockwise about the x axis. A side face joins the ring
        point `j` to `j + 1` modulo `sides`, so the quad of the side `j` is
        `(a[j], a[j + 1], b[j + 1], b[j])` of the two rings `a` and `b`, and
        the two caps are the rings themselves. VTK has no cell type of a
        prism of `sides` corners: a consumer builds a polyhedron from these
        faces (see the documentation page).

    Raises:
        ValueError: for less than one zone or less than three sides.
    """
    if n_zones < 1:
        raise ValueError(f"A tube has at least one zone, not {n_zones}.")
    if sides < 3:
        raise ValueError(f"A ring has at least three sides, not {sides}.")
    angle = np.linspace(0.0, 2.0 * np.pi, sides, endpoint=False)
    x = np.linspace(0.0, geo.L_sin, n_zones + 1)
    points = np.array(
        [(xi, geo.r_sin * np.cos(a), geo.r_sin * np.sin(a)) for xi in x for a in angle]
    )
    ring = np.arange(sides)
    cells = np.array(
        [
            np.concatenate([k * sides + ring, (k + 1) * sides + ring])
            for k in range(n_zones)
        ]
    )
    return points, cells
