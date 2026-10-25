"""The data files of the PBPK library.

The data are TOML files in `pkpdsim/pbpk/data`. A number is a table
`{ value = 1.8, unit = "kg", source = "ICRP89" }` whose `source` is a key of
`sources.toml`, so every number states where it comes from.
"""

import tomllib
from dataclasses import dataclass
from functools import cache
from importlib import resources
from typing import Any

from sbmlutils.factory import UnitDefinition

from pkpdsim.pbpk.units import U
from pkpdsim.pbpk.values import Value

__all__ = ["Source", "cite", "load", "sources", "value"]

_UNITS: dict[str, UnitDefinition] = {
    "dimensionless": U.dimensionless,
    "l": U.litre,
    "s": U.second,
    "min": U.min,
    "1/min": U.per_min,
    "1/hr": U.per_hr,
}
"""The units of the model for the units of the data files read by `value`."""


@dataclass(frozen=True)
class Source:
    """A publication a number is taken from.

    Attributes:
        key: the key of the source, e.g. `ICRP89`.
        citation: the citation.
        doi: the DOI of the publication.
    """

    key: str
    citation: str
    doi: str | None = None


def load(name: str) -> dict[str, Any]:
    """Load a data file.

    Args:
        name: the name of the file without extension, e.g. `sources`.

    Returns:
        The content of `data/<name>.toml`.
    """
    text = (
        resources.files("pkpdsim.pbpk")
        .joinpath("data", f"{name}.toml")
        .read_text("utf-8")
    )
    return tomllib.loads(text)


def value(name: str, section: str, key: str) -> Value:
    """Get a number of a data file with its unit, source and notes.

    Args:
        name: the name of the file without extension, e.g. `routes`.
        section: the table of the number.
        key: the key of the number in the table.

    Returns:
        The value.

    Raises:
        ValueError: if the unit of the number is not a unit of the model.
    """
    entry = load(name)[section][key]
    try:
        unit = _UNITS[entry["unit"]]
    except KeyError:
        raise ValueError(
            f"{name}.toml [{section}] {key}: unit '{entry['unit']}' is not a unit "
            "of the model."
        ) from None
    return Value(entry["value"], unit, entry["source"], entry.get("notes"))


@cache
def sources() -> dict[str, Source]:
    """Get the sources by their key."""
    return {
        key: Source(key=key, citation=entry["citation"], doi=entry.get("doi"))
        for key, entry in load("sources").items()
    }


def cite(key: str) -> str:
    """Get the citation of a source with its DOI.

    Args:
        key: the key of the source.

    Returns:
        The citation, followed by the DOI if the source has one.

    Raises:
        KeyError: if no source has the key.
    """
    try:
        source = sources()[key]
    except KeyError:
        raise KeyError(f"Unknown source '{key}', see data/sources.toml.") from None
    if source.doi:
        return f"{source.citation} doi:{source.doi}"
    return source.citation
