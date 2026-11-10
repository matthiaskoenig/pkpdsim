"""The geometry of the zonated liver for renderers: the tube and the lobule hexagon.

Both are plain numpy arrays and know no figure library. `lobule_bands` is the
polygons of the idealized hexagonal lobule which
`pkpdsim.zonation.figures.draw_lobule` colours and `lobule_isoline` one contour
of the lobular position over the same hexagon, the outline of the lobule (the
position 0), the thirds of a histology or any other level; `tube_mesh` is the
sinusoid as points and cells.

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
from pkpdsim.zonation.profiles import lobular_edge, lobular_radius

__all__ = ["lobule_bands", "lobule_isoline", "tube_mesh"]


_WEDGE = math.pi / 6.0
"""The angle of one twelfth of the hexagon, the wedge of one portal triad."""
_ARC_POINTS = 16
"""The directions sampled over the arc of an isoline by default, what a drawing needs.

A polygon is the chords of its isolines, so its area falls a little short of
the area of the zone; the areas of the engine come from `lobular_area`
analytically and never from these polygons, so the sampling is a question of
the picture alone. At 16 directions the area of a zone is within 8.6e-4 of
`zone_areas` (the smallest zone, the one at the triads; the others are nearer)
and the sagitta of a chord is below 1.2e-4 of the circumradius, a twentieth of
a pixel in a panel of 3 inches at 200 dpi, while the SVG of a panel of ten
zones is about 113 kB. A caller who measures the polygons themselves asks for
more: `points=256` is within 3e-6, at 1.4 MB a panel.
"""


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


def _directions(position: float, points: int) -> np.ndarray:
    """Get the directions which the isoline of a position is sampled in.

    The isoline is the arc of an Apollonius circle up to the direction in which
    it leaves through the edge of the hexagon and the straight edge beyond it,
    so `points` directions cover the arc and the midpoint of the edge (theta
    `pi / 6`) closes it: a chord of the straight part is the part itself, and
    no point is spent where the polygon is exact already. The direction is
    `lobular_edge` of the profiles, the one the engine integrates the areas of
    the zones with, so the polygons and `zone_areas` cut the lobule in the same
    place.

    The isoline of 0 is therefore two directions, the corner and the midpoint
    of the edge, and the isoline of 1 is the central vein, the radius 0 in
    every direction, which two directions give as well.
    """
    if position >= 1.0:
        return np.array([0.0, _WEDGE])
    theta = np.linspace(0.0, lobular_edge(position), points)
    return np.unique(np.append(theta, _WEDGE))


def _arc(position: float, points: int) -> np.ndarray:
    """Get the polyline of one position over the wedge of a triad, `(n, 2)`."""
    theta = _directions(position, points)
    radius = np.array([lobular_radius(position, t) for t in theta])
    return np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])


def _band(low: float, high: float, points: int) -> np.ndarray:
    """Get the polygon of the band `[low, high]` of a position in the wedge of a triad.

    The outer arc is the isoline of `low` and the inner arc that of `high`,
    reversed, each sampled in its own directions (`_directions`): the isoline
    of a position is one polyline wherever it is used, so the inner arc of a
    band is point for point the outer arc of the band inside it and the two
    leave no seam, and `lobule_isoline` draws the same line.
    """
    return np.vstack([_arc(low, points), _arc(high, points)[::-1]])


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


def lobule_bands(n: int, *, points: int = _ARC_POINTS) -> list[list[np.ndarray]]:
    """Get the polygons of the zones of the idealized hexagonal lobule.

    The hexagon has the circumradius 1, the central vein (zone `n`) at the
    origin and the portal triads (zone 1) at the corners on the angles 0, 60,
    ... degrees. The band of zone `k` is `[(k - 1) / n, k / n]` of the
    normalized periportal to pericentral distance, in twelve wedges.

    Args:
        n: the number of zones.
        points: the directions sampled over one wedge, the resolution of the
            polygons; the default is what a drawing needs (`_ARC_POINTS`), a
            caller who measures the polygons themselves asks for more.

    Returns:
        For every zone from periportal, its twelve polygons, arrays of the
        points of the two isolines which bound the band; the areas of the
        zones are `zone_areas(n)` to the chord error of the isolines (at 10
        zones within 8.6e-4 of a zone at the default, 3e-6 at `points=256`).

    Raises:
        ValueError: for less than one zone or less than two points.
    """
    if n < 1:
        raise ValueError(f"A lobule has at least one zone, not {n}.")
    if points < 2:
        raise ValueError(f"A wedge has at least two points, not {points}.")
    return [_twelve(_band(low, high, points)) for low, high in pairwise(_zone_edges(n))]


def lobule_isoline(position: float, *, points: int = _ARC_POINTS) -> np.ndarray:
    """Get the closed polyline of the points of one lobular position.

    The isoline over the whole hexagon of `lobule_bands`: the twelve copies of
    the arc of the wedge `[0, pi / 6]` (`lobular_radius`), joined into one ring
    around the centre, counterclockwise from the triad at the angle 0. A plot
    of the ring is the contour of the position and a `Polygon` of it the area
    at or beyond the position, the pericentral side, whose share of the hexagon
    is `1 - lobular_area(position)`.

    The isoline of 0 is the hexagon itself, of the circumradius 1 towards a
    portal triad and of the apothem `cos(pi / 6)` towards the middle of an
    edge; the isoline of 1 is the centre, the central vein, every radius 0. The
    radius in a direction does not grow with the position, and where an isoline
    has left the hexagon it follows the edge, so two positions share that part
    of their rings.

    Args:
        position: the lobular position, 0 periportal (the portal triads) and 1
            pericentral (the central vein).
        points: the directions sampled over one wedge; the kink where the
            isoline leaves through an edge is always one of them. The default
            is what a drawing needs (`_ARC_POINTS`), as for `lobule_bands`.

    Returns:
        The points of the ring, an array of shape `(n, 2)` in the hexagon of
        the circumradius 1 with its centre at the origin, in order around the
        centre; the last point is the first, so the ring is closed.

    Raises:
        ValueError: for a position outside the lobule, or for less than two
            points per wedge.
    """
    if points < 2:
        raise ValueError(f"A wedge has at least two points, not {points}.")
    wedge = _arc(position, points)
    # the odd wedges are the mirror images, which run from the edge back to the
    # triad: reversed, so that the ring goes around the centre in one direction
    pieces = [
        copy if k % 2 == 0 else copy[::-1] for k, copy in enumerate(_twelve(wedge))
    ]
    # two wedges share the direction where they meet, so a copy drops its first
    # point; the last point of the last copy is the first point of the first one
    ring = np.vstack([pieces[0], *(piece[1:] for piece in pieces[1:])])
    ring[-1] = ring[0]  # the rotation of the last copy closes it to rounding only
    return ring


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
