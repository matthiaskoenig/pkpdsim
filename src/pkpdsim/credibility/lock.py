"""The lock of the evidence of a credibility declaration: splits fixed in advance.

A validation counts only if its data was set apart before the predictions were
seen. The kinds of the mappings (training or validation) are fixed by the
records of the adopted fits (`pkpdsim.credibility.records`); which validation
data answers which question is fixed by the declaration, which a modeller can
change after looking at a verdict. The lock records the mappings every evidence
and every limitation selects, `credibility.lock.json` next to the declaration
in the folder of the model, written on request (`run_credibility --lock`) and
reviewed in the history of the repository:

- `write_lock` writes the selection of a declaration on the table of the fits.
- `lock_differences` names every evidence whose mappings differ from the lock;
  `check_declaration` reports them and `assess` raises a `LockMismatchError`,
  so a changed split is renewed explicitly with `--lock`.

A model without a lock is assessed as before; the record says whether the
evidence matched a lock (`Provenance.lock`).
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Final

if TYPE_CHECKING:
    from pathlib import Path

    from pkpdsim.credibility.records import FitTable
    from pkpdsim.credibility.selection import EvidenceSelection, Selections
    from pkpdsim.model import PKPDModel

__all__ = [
    "LOCK_FILE",
    "LOCK_REMEDY",
    "LockMismatchError",
    "lock_differences",
    "lock_path",
    "read_lock",
    "selection_lock",
    "write_lock",
]

#: the file of the lock in the folder of the model
LOCK_FILE: Final = "credibility.lock.json"
#: the version of the format of the lock
_VERSION: Final = 1
#: the remedy of a selection which differs from the lock
LOCK_REMEDY: Final = (
    "Review the change of the split and renew the lock with "
    "`run_credibility <package> --lock`."
)

#: the mappings of the evidence by its key, `fit/experiment/mapping` each
type Lock = dict[str, list[str]]


class LockMismatchError(ValueError):
    """The evidence selects other mappings than the lock of the declaration.

    Attributes:
        problems: every evidence which differs, with its added and removed
            mappings.
    """

    def __init__(self, problems: list[str]) -> None:
        """Create the error.

        Args:
            problems: every evidence which differs.
        """
        self.problems = list(problems)
        super().__init__(
            "The evidence of the credibility declaration differs from its lock:\n"
            + "\n".join(f"- {problem}" for problem in self.problems)
            + f"\n{LOCK_REMEDY}"
        )


def lock_path(model: PKPDModel) -> Path:
    """Get the file of the lock of a model, `credibility.lock.json` in its folder."""
    return model.path / LOCK_FILE


def _key(selection: EvidenceSelection) -> str:
    """Get the key of an evidence: `context/name`, `/name` for one of the model."""
    return f"{selection.context}/{selection.evidence.name}"


def selection_lock(selections: Selections, fits: FitTable) -> Lock:
    """Get the mappings of every evidence and limitation of a selection.

    Args:
        selections: the selection of the declaration, see `select`.
        fits: the table of the fits the selection was made on.

    Returns:
        The sorted mappings, `fit/experiment/mapping`, by the key of the
        evidence, `context/name` (a limitation of the model `/name`).
    """
    table = fits.table
    lock: Lock = {}
    for selection in (*selections.evidence, *selections.limitation_evidence):
        rows = table.loc[selection.rows, ["fit", "experiment", "mapping"]].astype(str)
        mappings = {"/".join(row) for row in rows.itertuples(index=False)}
        lock.setdefault(_key(selection), [])
        lock[_key(selection)] = sorted({*lock[_key(selection)], *mappings})
    return dict(sorted(lock.items()))


def write_lock(model: PKPDModel, selections: Selections, fits: FitTable) -> Path:
    """Write the lock of the evidence of a model.

    Args:
        model: the model.
        selections: the selection of its declaration.
        fits: the table of the fits the selection was made on.

    Returns:
        The file of the lock.
    """
    path = lock_path(model)
    data = {"version": _VERSION, "evidence": selection_lock(selections, fits)}
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return path


def read_lock(model: PKPDModel) -> Lock | None:
    """Read the lock of a model, `None` without one.

    Raises:
        ValueError: if the file is no lock of this version.
    """
    path = lock_path(model)
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("version") != _VERSION:
        raise ValueError(f"The file '{path}' is no lock of version {_VERSION}.")
    return {str(k): [str(m) for m in v] for k, v in data["evidence"].items()}


def lock_differences(
    model: PKPDModel, selections: Selections, fits: FitTable
) -> list[str] | None:
    """Get every evidence whose mappings differ from the lock of the model.

    Args:
        model: the model.
        selections: the current selection of its declaration.
        fits: the table of the fits the selection was made on.

    Returns:
        A line per evidence which is new, gone or selects other mappings, with
        the added (`+`) and removed (`-`) mappings; empty if the selection
        matches the lock, `None` if the model has no lock.
    """
    locked = read_lock(model)
    if locked is None:
        return None
    current = selection_lock(selections, fits)
    problems: list[str] = []
    for key in sorted(locked.keys() | current.keys()):
        if key not in current:
            problems.append(f"the evidence '{key}' of the lock is gone")
            continue
        if key not in locked:
            problems.append(f"the evidence '{key}' is not in the lock")
            continue
        added = sorted(set(current[key]) - set(locked[key]))
        removed = sorted(set(locked[key]) - set(current[key]))
        if added or removed:
            changes = [*(f"+{m}" for m in added), *(f"-{m}" for m in removed)]
            problems.append(
                f"the evidence '{key}' selects other mappings: {', '.join(changes)}"
            )
    return problems
