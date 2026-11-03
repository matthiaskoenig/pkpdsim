"""The static check of the credibility declaration of a model, without a simulation."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pkpdsim.credibility.declaration import Credibility, DeclarationError
from pkpdsim.credibility.records import kind_differences, mapping_table
from pkpdsim.credibility.selection import SelectionError, select

if TYPE_CHECKING:
    from pkpdsim.model import PKPDModel

__all__ = ["check_declaration"]


def check_declaration(model: PKPDModel) -> list[str]:
    """Check the credibility declaration of a model against its adopted fits.

    The declaration must be valid, the kinds of the mappings of the current
    selection of the fit definitions must be the recorded kinds of the adopted
    fits, and the evidence must select valid data (see `select`), all from the
    fit definitions and the records, without a simulation.

    Args:
        model: the model with a declaration and adopted fits.

    Returns:
        Every problem, empty if the declaration holds.

    Raises:
        ValueError: if the model has no declaration or no adopted fits.
    """
    problems: list[str] = []
    credibility: Credibility | None = None
    try:
        credibility = model.load_credibility()
    except DeclarationError as err:
        problems.extend(err.problems)
    try:
        fits = mapping_table(model)
    except FileNotFoundError as err:
        return [*problems, str(err)]
    problems.extend(kind_differences(fits.table, fits.recorded_kinds))
    if credibility is not None:
        try:
            select(credibility, fits)
        except SelectionError as err:
            problems.extend(err.problems)
    return problems
