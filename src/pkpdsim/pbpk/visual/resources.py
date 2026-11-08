"""The 2D anatomy resources: the drawings of the bodies and organs of the figures.

One JSON file per species and sex in `pkpdsim/resources/anatomy/`, built by
`scripts/visual/anatomy_assets.py` from the anatomograms of the EBI Expression
Atlas (CC BY 4.0, see the `NOTICE.md` of the folder): the outline of the body
and its details, the organs of the library by their `Organ` key with their
UBERON ids, further anatomy, the landmarks of the human drawings and the body
cavities, all in the coordinates of the drawing (units of its view box, y
downwards).
"""

import json
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from types import MappingProxyType

import numpy as np

__all__ = ["RESOURCES", "Anatomy", "Shape", "anatomy", "parse_d"]

RESOURCES = Path(__file__).resolve().parents[2] / "resources" / "anatomy"
"""The folder of the anatomy resources."""


def parse_d(d: str) -> tuple[np.ndarray, ...]:
    """The polylines of a path of the resources (`M x,y x,y ... Z`, absolute).

    Args:
        d: the path data, subpaths of absolute points.

    Returns:
        One array of shape (n, 2) per subpath.
    """
    out = []
    for sub in d.split("M")[1:]:
        points = sub.strip().rstrip("Z").split()
        out.append(np.array([[float(v) for v in p.split(",")] for p in points]))
    return tuple(out)


@dataclass(frozen=True)
class Shape:
    """A drawn tissue: an organ of the library or further anatomy.

    Attributes:
        key: the key, an `Organ` value or the name of the anatomy (`bladder`).
        name: the name of the tissue.
        uberon: the UBERON terms of the drawing which form it.
        cavity: the body cavity it lies in (`thorax`, `abdomen`), `head`,
            `limb` or `None`.
        d: the path data.
        paths: the polylines of the path.
    """

    key: str
    name: str
    uberon: tuple[str, ...]
    cavity: str | None
    d: str
    paths: tuple[np.ndarray, ...]


@dataclass(frozen=True)
class Anatomy:
    """The drawing of a species and sex.

    Attributes:
        key: the resource (`human_male`).
        species: the species (`human`, `mouse`, `rat`).
        sex: the sex (`male`, `female`).
        viewbox: the box of the outline, x, y, width, height.
        midline: the x of the midline of the body.
        landmarks: the heights of the human drawings (`top`, `neck`,
            `shoulder`, `armpit`, `waist`, `crotch`, `bottom`), empty for the
            animals.
        outline: the outline of the body, one closed polyline.
        details: the small shapes of the outline (eyes, nose, mouth, hands).
        organs: the organs of the library by their `Organ` value.
        anatomy: further anatomy by name.
        cavities: the polylines of the cavities `thorax` and `abdomen`.
        attribution: the source, license and changes.
    """

    key: str
    species: str
    sex: str
    viewbox: tuple[float, float, float, float]
    midline: float
    landmarks: Mapping[str, float]
    outline: np.ndarray
    details: tuple[np.ndarray, ...]
    organs: Mapping[str, Shape]
    anatomy: Mapping[str, Shape]
    cavities: Mapping[str, tuple[np.ndarray, ...]]
    attribution: str

    def shape(self, key: str) -> Shape | None:
        """The organ or the anatomy of a key, `None` if the drawing lacks it."""
        return self.organs.get(key) or self.anatomy.get(key)


def _shapes(entries: Mapping[str, Mapping]) -> Mapping[str, Shape]:
    return MappingProxyType(
        {
            k: Shape(
                key=k,
                name=v["name"],
                uberon=tuple(v["uberon"]),
                cavity=v["cavity"],
                d=v["d"],
                paths=parse_d(v["d"]),
            )
            for k, v in entries.items()
        }
    )


@cache
def anatomy(species: str = "human", sex: str = "male") -> Anatomy:
    """The anatomy resource of a species and sex.

    Args:
        species: `human`, `mouse` or `rat` (male only).
        sex: `male` or `female`.

    Returns:
        The drawing, read once per process.

    Raises:
        ValueError: for a species and sex without a resource.
    """
    path = RESOURCES / f"{species}_{sex}.json"
    if not path.exists():
        known = sorted(p.stem for p in RESOURCES.glob("*.json"))
        raise ValueError(f"No anatomy for {species} {sex}; the resources are {known}.")
    r = json.loads(path.read_text(encoding="utf-8"))
    return Anatomy(
        key=r["id"],
        species=r["species"],
        sex=r["sex"],
        viewbox=tuple(r["viewbox"]),
        midline=float(r["midline"]),
        landmarks=MappingProxyType(dict(r.get("landmarks", {}))),
        outline=parse_d(r["outline"])[0],
        details=parse_d(r["details"]) if r["details"] else (),
        organs=_shapes(r["organs"]),
        anatomy=_shapes(r["anatomy"]),
        cavities=MappingProxyType(
            {c: parse_d(d) for c, d in r["cavities"].items() if d}
        ),
        attribution=r["source"]["attribution"],
    )
