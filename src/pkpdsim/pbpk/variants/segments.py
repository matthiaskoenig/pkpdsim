"""The lumen segments of the segmented intestine: their data and geometry.

The segments of Thelen 2011 in the order of the transit, each with its table
in `data/intestine_human.toml`; the share of an intestinal segment in its
region (small or large intestine) is its share of the transit time of the
region, by which the mucosa gets its tissue volume and blood flow. The
single lumen is one segment, the fluid of the small intestine.
"""

from dataclasses import dataclass

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import value
from pkpdsim.pbpk.organs import Organ

__all__ = ["DATA", "SEGMENTS", "SMALL", "Segment", "segments", "single_lumen"]

DATA = "intestine_human"
"""The data file of the intestine."""

SEGMENTS: dict[str, str] = dict(
    zip(
        ids.GUT_SEGMENTS,
        (
            "stomach",
            "duodenum",
            "upper_jejunum",
            "lower_jejunum",
            "upper_ileum",
            "lower_ileum",
            "caecum",
            "colon_ascendens",
            "colon_transversum",
            "colon_descendens",
            "colon_sigmoid",
            "rectum",
        ),
        strict=True,
    )
)
"""The data table of every lumen segment by its id token."""

SMALL = ("du", "jej1", "jej2", "ile1", "ile2")
"""The segments of the small intestine."""


@dataclass(frozen=True)
class Segment:
    """A lumen segment of the intestine.

    Attributes:
        token: the id token (`du`), `None` for the single lumen.
        key: the table of its data in `data/intestine_human.toml`.
        region: the organ whose tissue and blood flow its mucosa shares,
            `None` without mucosa (the stomach).
        share: the share of the segment in its region (the transit share).
        absorbs: the segment absorbs (with the area factor of the data).
    """

    token: str | None
    key: str
    region: Organ | None
    share: float
    absorbs: bool


def segments() -> tuple[Segment, ...]:
    """Get the lumen segments of the segmented intestine.

    The share of a segment in its region is its share of the transit time of
    the region (small intestine, large intestine).

    Returns:
        The segments in the order of the transit.
    """
    transit = {t: value(DATA, key, "transit_time").value for t, key in SEGMENTS.items()}
    small = sum(transit[t] for t in SMALL)
    large = sum(transit[t] for t in SEGMENTS if t not in SMALL and t != "st")
    result = []
    for token, key in SEGMENTS.items():
        if token == "st":
            result.append(Segment(token, key, None, 0.0, absorbs=False))
        elif token in SMALL:
            result.append(
                Segment(token, key, Organ.SMALL_INTESTINE, transit[token] / small, True)
            )
        else:
            result.append(
                Segment(token, key, Organ.LARGE_INTESTINE, transit[token] / large, True)
            )
    return tuple(result)


def single_lumen() -> tuple[Segment, ...]:
    """The one segment of the single lumen."""
    return (Segment(None, "single_lumen", Organ.SMALL_INTESTINE, 1.0, absorbs=False),)
