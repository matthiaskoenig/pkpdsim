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
}


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
