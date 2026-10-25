"""Tests of the data files and their sources."""

import tomllib
from importlib import resources
from typing import Any

import pytest

from pkpdsim.pbpk import cite, sources


def _tables(node: object) -> list[dict[str, Any]]:
    """Get all tables of a parsed file which state a value."""
    found: list[dict[str, Any]] = []
    if isinstance(node, dict):
        if "value" in node:
            found.append(node)
        for child in node.values():
            found.extend(_tables(child))
    elif isinstance(node, list):
        for child in node:
            found.extend(_tables(child))
    return found


def test_every_source_has_citation() -> None:
    """Every source states its key and a citation."""
    assert {"ICRP89", "Jones2013", "Brown1997"} <= set(sources())
    for key, s in sources().items():
        assert s.key == key
        assert s.citation


def test_every_data_value_has_known_source() -> None:
    """Every number of every data file has a source of sources.toml."""
    known = sources()
    for path in resources.files("pkpdsim.pbpk").joinpath("data").iterdir():
        if not path.name.endswith(".toml") or path.name == "sources.toml":
            continue
        for table in _tables(tomllib.loads(path.read_text("utf-8"))):
            assert table.get("source") in known, f"{path.name}: {table}"


def test_cite_includes_doi() -> None:
    """The citation names the DOI."""
    assert "10.1038/psp.2013.41" in cite("Jones2013")


def test_cite_unknown_key_raises() -> None:
    """An unknown key raises a KeyError which names the key."""
    with pytest.raises(KeyError, match="nope"):
        cite("nope")
