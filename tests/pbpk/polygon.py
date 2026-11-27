"""Whether points lie inside a polygon, for the tests of the outlines of the bodies.

`contains_points` is the crossing test of `matplotlib.path.Path.contains_points`
for a polygon without codes: a point is inside when a ray from it to the right
crosses the edges of the closed polygon an odd number of times (the even-odd
rule), so the inner pentagon of a pentagram is outside. A point on an edge has
no defined answer, as in matplotlib.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def contains_points(polygon: ArrayLike, points: ArrayLike) -> NDArray[np.bool_]:
    """Get whether each point lies inside the closed polygon, by the even-odd rule.

    Args:
        polygon: the vertices, `(n, 2)`; the last one is joined to the first,
            which may also be repeated at the end.
        points: the points, `(m, 2)`.

    Returns:
        One boolean per point.
    """
    vertices = np.asarray(polygon, dtype=float)
    xy = np.asarray(points, dtype=float).reshape(-1, 2)
    x, y = xy[:, :1], xy[:, 1:]
    x0, y0 = vertices[:, 0], vertices[:, 1]
    x1, y1 = np.roll(x0, -1), np.roll(y0, -1)
    # an edge crosses the horizontal line of a point when its ends lie on either
    # side of it, the lower end included: a vertex on the line counts once
    spans = (y0 > y) != (y1 > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        x_line = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
    crossings = np.count_nonzero(spans & (x < x_line), axis=1)
    return crossings % 2 == 1
