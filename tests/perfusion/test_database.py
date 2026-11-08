"""Tests of the parameter database of a machine perfusion."""

from pkpdsim.database import tables

SECTIONS = ("oxygen", "enzymes", "temperature", "perfusate", "liver_rat")


def _entries() -> list[tuple[str, str]]:
    """Get the section and key of every number of the table."""
    data = tables.load("perfusion")
    return [(s, k) for s in SECTIONS for k in data[s]]


def test_the_table_has_its_sections() -> None:
    """The sections the library reads."""
    assert set(tables.load("perfusion")) >= set(SECTIONS)


def test_every_entry_has_a_source() -> None:
    """A source of `sources.toml`, which may be `ModelChoice` with a note."""
    for section, key in _entries():
        entry = tables.entry("perfusion", section, key)
        assert entry.source in tables.sources(), (section, key)
        if entry.source == "ModelChoice":
            assert entry.notes, (section, key)


def test_the_enzyme_activities_are_positive() -> None:
    """ALT, AST and LDH, in µkat per mmole of enzyme."""
    for key in ("ALT", "AST", "LDH"):
        entry = tables.entry("perfusion", "enzymes", key)
        assert isinstance(entry.value, float)
        assert entry.value > 0
        assert entry.unit == "ukat/mmol"


def test_the_oxygen_solubility() -> None:
    """The solubility has the unit mM/kPa and a plausible size at 37 °C."""
    entry = tables.entry("perfusion", "oxygen", "solubility")
    assert entry.unit == "mM/kPa"
    assert isinstance(entry.value, float)
    assert 0.008 < entry.value < 0.011
    coefficient = tables.entry("perfusion", "oxygen", "temperature_coefficient")
    assert isinstance(coefficient.value, float)
    assert coefficient.value < 0


def test_the_temperature_reference() -> None:
    """The reference is 37 °C."""
    assert tables.entry("perfusion", "temperature", "reference").value == 37.0
