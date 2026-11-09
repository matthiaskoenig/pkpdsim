"""The sinusoid as a tube mesh: points and cells of plain arrays.

The tube is the geometry of the engine (`SinusoidGeometry`): a cylinder of the
sinusoid radius and length along the x axis, from the portal end (x = 0) to
the central end (x = `L_sin`), with one segment per zone. The arrays are what a
renderer of unstructured meshes takes; the engine neither depends on one nor
renders anything: a script which draws (the videos of gdeobdtpa) installs
pyvista for its run, builds a grid of these arrays and colours the segments by
the zone values, `values[k]` for the segment `k`.
"""

import numpy as np

from pkpdsim.zonation.geometry import SinusoidGeometry

__all__ = ["tube_mesh"]


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
        segment `k` are the `sides` points of the ring `k` followed by those of
        the ring `k + 1`, counterclockwise about the x axis, the prism of the
        zone `k + 1`.

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
