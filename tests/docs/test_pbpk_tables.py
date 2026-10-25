"""The data tables of the PBPK pages are those of the data files of the library.

`scripts/pbpk_tables.py` writes the tables (reference physiology,
populations, impairment, arm vein presets, intestine, nephron) into the
marked regions of the pages of `docs/pbpk/` from the TOML files of
`pkpdsim/pbpk/data`. The test renders the pages again, so a change of a data
file or of the script fails until the pages are written again.
"""

import importlib.util
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).parent.parent.parent


def _script() -> ModuleType:
    """Import `scripts/pbpk_tables.py` from its path."""
    path = ROOT / "scripts" / "pbpk_tables.py"
    spec = importlib.util.spec_from_file_location("pbpk_tables", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_tables_are_current() -> None:
    """Every generated region of a page is what the script writes."""
    script = _script()
    stale = [
        page.name
        for page in sorted(script.PBPK_DIR.glob("*.md"))
        if script.render_page(page.read_text(encoding="utf-8"))
        != page.read_text(encoding="utf-8")
    ]
    assert not stale, (
        f"the tables of {stale} are out of date, run "
        "`uv run python scripts/pbpk_tables.py`"
    )


def test_every_table_is_on_a_page() -> None:
    """Every table of the script has a region on a page, and only one."""
    script = _script()
    text = "".join(
        p.read_text(encoding="utf-8") for p in sorted(script.PBPK_DIR.glob("*.md"))
    )
    counts = {name: text.count(script.begin_marker(name)) for name in script.TABLES}
    assert counts == dict.fromkeys(script.TABLES, 1)


def test_an_unknown_table_is_refused() -> None:
    """A region of a name the script does not know is an error, not left alone."""
    script = _script()
    page = f"{script.begin_marker('nothing')}\nold\n{script.END}\n"
    try:
        script.render_page(page)
    except KeyError as err:
        assert "nothing" in str(err)
    else:
        raise AssertionError("an unknown table was accepted")
