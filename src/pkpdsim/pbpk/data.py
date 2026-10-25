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

__all__ = ["Source", "cite", "load", "sources"]


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
