"""The lumen segments of the segmented intestine: their data and geometry.

The segments of Thelen 2011 in the order of the transit, each with its table
in `database/tables/intestine_human.toml`; the share of an intestinal segment in its
region (small or large intestine) is its share of the transit time of the
region, by which the mucosa gets its tissue volume and blood flow. The
single lumen is one segment, the fluid of the small intestine.
"""

from dataclasses import dataclass

from pkpdsim.pbpk import ids
from pkpdsim.pbpk.data import value
from pkpdsim.pbpk.organs import Organ
from pkpdsim.species import Species

__all__ = [
    "DATA",
    "SEGMENTS",
    "SMALL",
    "Segment",
    "intestine_table",
    "segments",
    "single_lumen",
]

DATA = "intestine_human"
"""The data file of the intestine of the human."""


def intestine_table(species: Species) -> str:
    """The data file of the intestine of a species, `intestine_rat`."""
    return f"intestine_{species.value}"


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
        key: the table of its data in `database/tables/intestine_human.toml`.
        region: the organ whose tissue and blood flow its mucosa shares,
            `None` without mucosa (the stomach).
        share: the share of the segment in its region (the transit share).
        absorbs: the segment absorbs (with the area factor of the data).
        table: the data file of the intestine of the species.
    """

    token: str | None
    key: str
    region: Organ | None
    share: float
    absorbs: bool
    table: str = DATA

    @property
    def label(self) -> str:
        """The segment in names, `upper jejunum`."""
        return self.key.replace("_", " ")


def segments(table: str = DATA) -> tuple[Segment, ...]:
    """Get the lumen segments of the segmented intestine.

    The share of a segment in its region is its share of the transit time of
    the region (small intestine, large intestine).

    Args:
        table: the data file of the intestine of the species
            (`intestine_table`), the human by default.

    Returns:
        The segments in the order of the transit.
    """
    transit = {
        t: value(table, key, "transit_time").value for t, key in SEGMENTS.items()
    }
    small = sum(transit[t] for t in SMALL)
    large = sum(transit[t] for t in SEGMENTS if t not in SMALL and t != "st")
    result = []
    for token, key in SEGMENTS.items():
        if token == "st":
            result.append(Segment(token, key, None, 0.0, absorbs=False, table=table))
        elif token in SMALL:
            share = transit[token] / small
            result.append(
                Segment(token, key, Organ.SMALL_INTESTINE, share, True, table)
            )
        else:
            share = transit[token] / large
            result.append(
                Segment(token, key, Organ.LARGE_INTESTINE, share, True, table)
            )
    return tuple(result)


def single_lumen(table: str = DATA) -> tuple[Segment, ...]:
    """The one segment of the single lumen, of the data file of the species."""
    return (
        Segment(
            None, "single_lumen", Organ.SMALL_INTESTINE, 1.0, absorbs=False, table=table
        ),
    )
