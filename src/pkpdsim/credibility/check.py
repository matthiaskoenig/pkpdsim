"""The static check of the credibility declaration of a model, without a simulation."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pkpdsim.credibility.declaration import Credibility, DeclarationError
from pkpdsim.credibility.lock import LOCK_REMEDY, lock_differences
from pkpdsim.credibility.records import KIND_REMEDY, kind_differences, mapping_table
from pkpdsim.credibility.selection import SelectionError, select

if TYPE_CHECKING:
    from pkpdsim.model import PKPDModel

__all__ = ["check_declaration"]


def check_declaration(model: PKPDModel) -> list[str]:
    """Check the credibility declaration of a model against its adopted fits.

    The declaration must be valid, the kinds of the mappings of the current
    selection of the fit definitions must be the recorded kinds of the adopted
    fits, and the evidence must select valid data (see `select`), all from the
    fit definitions and the records, without a simulation, and select the
    mappings of its lock, if the model has one. A problem of the
    records, e.g. a record without the kinds of its mappings, two kinds of one
    mapping, no adopted fits or a subset adopted twice, is a problem as well,
    after which the kinds and the evidence are not checked.

    Args:
        model: the model with a declaration and adopted fits.

    Returns:
        Every problem, empty if the declaration holds.

    Raises:
        ValueError: if the model has no declaration.
    """
    problems: list[str] = []
    credibility: Credibility | None = None
    try:
        credibility = model.load_credibility()
    except DeclarationError as err:
        problems.extend(err.problems)
    try:
        fits = mapping_table(model)
    except (FileNotFoundError, ValueError) as err:
        return [*problems, str(err)]
    differences = kind_differences(fits.table, fits.recorded_kinds)
    if differences:
        problems.extend(differences)
        problems.append(
            "The kinds of the fit mappings differ from the kinds the adopted fits "
            f"recorded. {KIND_REMEDY}"
        )
    if credibility is not None:
        try:
            selections = select(credibility, fits)
        except SelectionError as err:
            problems.extend(err.problems)
        else:
            try:
                locked = lock_differences(model, selections, fits)
            except ValueError as err:
                problems.append(str(err))
            else:
                if locked:
                    problems.extend(locked)
                    problems.append(LOCK_REMEDY)
    return problems
