"""The table of the id scheme on the overview page is the id scheme of the library.

A row of the table `docs/pbpk/index.md` names the pattern of an id and the
function of `pkpdsim.pbpk.ids` which makes it:

    | `fu_<s>` | the unbound fraction in plasma | `ids.fu(s)` |

The test calls the function with example values of the placeholders and
compares its result with the pattern filled with the same values.
"""

import re
from pathlib import Path

from pkpdsim.pbpk import Organ, Part, Pool, ids

PAGE = Path(__file__).parent.parent.parent / "docs" / "pbpk" / "index.md"

#: a row of the id table
ROW = re.compile(
    r"^\| `(?P<pattern>[^`]+)` \| .+ \| `(?P<call>ids\.[^`]+)` \|$",
    re.M,
)

#: the values of the placeholders in the calls
NAMESPACE = {
    "ids": ids,
    "Organ": Organ,
    "Part": Part,
    "Pool": Pool,
    "s": "api",
    "o": Organ.LIVER,
    "pool": Pool.VENOUS,
    "part": Part.PLASMA,
    "protein": "CYP3A4",
    "seg": "du",
    "zone": "z1",
    "kind": "iv",
    "tseg": "pt",
}

#: the placeholders of the patterns and the same values
PLACEHOLDERS = {
    "<s>": "api",
    "<o>": "li",
    "<O>": "LI",
    "<pool>": "ve",
    "<part>": "pl",
    "<protein>": "cyp3a4",
    "<seg>": "du",
    "<p>": "API2M",
    "<zone>": "z1",
    "<kind>": "iv",
    "<tseg>": "pt",
}

#: the names of `ids.__all__` which make no id and have no row: the tokens of
#: the segments and the list of the species of a substance
NO_ROW = frozenset({"GUT_SEGMENTS", "TUBULE_SEGMENTS", "species_of"})


def rows() -> list[re.Match[str]]:
    """Get the rows of the id table."""
    return list(ROW.finditer(PAGE.read_text(encoding="utf-8")))


def test_the_table_is_found() -> None:
    """The overview page has the id table."""
    assert len(rows()) >= 20


def test_every_example_is_made_by_its_function() -> None:
    """The function of a row makes the id of its pattern."""
    problems = []
    for row in rows():
        made = eval(row["call"], dict(NAMESPACE))
        pattern = row["pattern"]
        for key, value in PLACEHOLDERS.items():
            pattern = pattern.replace(key, value)
        if made != pattern:
            problems.append(f"{row['call']} = {made}, the pattern gives {pattern}")
    assert not problems, "\n".join(problems)


def test_every_id_of_the_scheme_has_a_row() -> None:
    """Every name of `ids.__all__` is the call of a row, so a new id is documented."""
    documented = {row["call"].removeprefix("ids.").split("(")[0] for row in rows()}
    missing = sorted(set(ids.__all__) - documented - NO_ROW)
    assert not missing, f"ids without a row in docs/pbpk/index.md: {missing}"
    assert not NO_ROW & documented


def test_every_public_name_of_the_ids_is_exported() -> None:
    """The functions and constants of `pkpdsim.pbpk.ids` are its `__all__`."""
    public = {
        name
        for name, value in vars(ids).items()
        if not name.startswith("_")
        and getattr(value, "__module__", ids.__name__) == ids.__name__
        and (callable(value) or name.isupper())
    }
    assert public - {"annotations"} <= set(ids.__all__), sorted(
        public - set(ids.__all__)
    )
