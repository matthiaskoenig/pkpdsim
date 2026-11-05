"""The tables and the sources of the database.

A table is a TOML file in `pkpdsim/database/tables` (a subfolder, e.g.
`populations`, is a group of tables). A number is a TOML table
`{ value = 1.8, unit = "kg", source = "ICRP89", notes = "..." }` whose
`source` is a key of `sources.toml`, so every number states where it comes
from. A number may state its uncertainty as its source reports it: `sd` (the
standard deviation of the individuals), `se` (the standard error of the
value), `cv` (the coefficient of variation, a fraction), `n` (the number of
individuals) and `ci` (`[lower, upper]`, the 95 % confidence interval of the
value); `entry` reads it into `Uncertainty`.
"""

import copy
import tomllib
from collections.abc import Iterator
from dataclasses import dataclass
from functools import cache
from importlib import resources
from importlib.resources.abc import Traversable
from typing import Any

__all__ = [
    "Entry",
    "Source",
    "Uncertainty",
    "cite",
    "entry",
    "load",
    "load_folder",
    "sources",
]

_UNCERTAINTY = ("sd", "se", "cv", "n", "ci")
"""The keys of the uncertainty of a number of a table."""


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


@dataclass(frozen=True)
class Uncertainty:
    """The uncertainty of a number as its source states it.

    Attributes:
        sd: the standard deviation of the individuals, in the unit of the number.
        se: the standard error of the number.
        cv: the coefficient of variation of the individuals, a fraction.
        n: the number of individuals.
        ci: the 95 % confidence interval of the number, `(lower, upper)`.
    """

    sd: float | None = None
    se: float | None = None
    cv: float | None = None
    n: int | None = None
    ci: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        """Check the values.

        Raises:
            ValueError: for a negative spread, an `n` below 1 or an interval
                whose lower bound exceeds the upper one.
        """
        for name in ("sd", "se", "cv"):
            v = getattr(self, name)
            if v is not None and v < 0:
                raise ValueError(f"The uncertainty '{name}' is negative: {v}.")
        if self.n is not None and self.n < 1:
            raise ValueError(f"The uncertainty 'n' is below 1: {self.n}.")
        if self.ci is not None and self.ci[0] > self.ci[1]:
            raise ValueError(f"The interval {self.ci} is not ordered.")


@dataclass(frozen=True)
class Entry:
    """A number of a table with its unit, provenance and uncertainty.

    Attributes:
        value: the number (or the numbers of a vector).
        unit: the unit as the table writes it, e.g. `1/hr`.
        source: the key of the source.
        notes: a remark on the number.
        uncertainty: the uncertainty the source states, `None` without.
    """

    value: float | tuple[float, ...]
    unit: str
    source: str
    notes: str | None = None
    uncertainty: Uncertainty | None = None


def _root() -> Traversable:
    """The folder of the database."""
    return resources.files("pkpdsim.database")


def load(name: str) -> dict[str, Any]:
    """Load a table.

    Args:
        name: the name of the table without extension, e.g. `routes`.

    Returns:
        The content of `tables/<name>.toml`, a copy which the caller may change.
    """
    return copy.deepcopy(_parsed(name))


@cache
def _parsed(name: str) -> dict[str, Any]:
    """Parse a table once per process; the content must not be changed.

    Args:
        name: the name of the table without extension, `sources` for the
            sources of the database.

    Returns:
        The content of the file, shared by every caller.
    """
    path = (
        _root().joinpath("sources.toml")
        if name == "sources"
        else _root().joinpath("tables", f"{name}.toml")
    )
    return tomllib.loads(path.read_text("utf-8"))


def load_folder(name: str) -> Iterator[tuple[str, dict[str, Any]]]:
    """Load the tables of a folder of the tables, by name in order.

    Args:
        name: the folder of `tables`, e.g. `populations`.

    Yields:
        The name of every table of the folder without extension and its content.
    """
    folder = _root().joinpath("tables", name)
    for path in sorted(folder.iterdir(), key=lambda p: p.name):
        if path.name.endswith(".toml"):
            yield (
                path.name.removesuffix(".toml"),
                tomllib.loads(path.read_text("utf-8")),
            )


def _uncertainty(data: dict[str, Any], where: str) -> Uncertainty | None:
    """The uncertainty of a number of a table, `None` without."""
    if not any(key in data for key in _UNCERTAINTY):
        return None
    ci = data.get("ci")
    if ci is not None and len(ci) != 2:
        raise ValueError(f"{where}: the interval 'ci' is not [lower, upper].")
    try:
        return Uncertainty(
            sd=data.get("sd"),
            se=data.get("se"),
            cv=data.get("cv"),
            n=data.get("n"),
            ci=(float(ci[0]), float(ci[1])) if ci is not None else None,
        )
    except ValueError as err:
        raise ValueError(f"{where}: {err}") from None


def entry(name: str, section: str, key: str) -> Entry:
    """Get a number of a table with its unit, source, notes and uncertainty.

    Args:
        name: the name of the table without extension, e.g. `routes`.
        section: the table of the number, dotted for a nested table.
        key: the key of the number in the table.

    Returns:
        The number.

    Raises:
        KeyError: if the table has no such number.
        ValueError: for a malformed uncertainty or a source which
            `sources.toml` does not have.
    """
    data: Any = _parsed(name)
    for part in section.split("."):
        data = data[part]
    data = data[key]
    where = f"tables/{name}.toml [{section}] {key}"
    if data["source"] not in sources():
        raise ValueError(f"{where}: unknown source '{data['source']}'.")
    raw = data["value"]
    value = tuple(float(v) for v in raw) if isinstance(raw, list) else float(raw)
    return Entry(
        value=value,
        unit=data["unit"],
        source=data["source"],
        notes=data.get("notes"),
        uncertainty=_uncertainty(data, where),
    )


@cache
def sources() -> dict[str, Source]:
    """Get the sources by their key."""
    return {
        key: Source(key=key, citation=entry["citation"], doi=entry.get("doi"))
        for key, entry in _parsed("sources").items()
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
        raise KeyError(
            f"Unknown source '{key}', see pkpdsim/database/sources.toml."
        ) from None
    if source.doi:
        return f"{source.citation} doi:{source.doi}"
    return source.citation
