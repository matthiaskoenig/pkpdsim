"""The data points of the evidence of the contexts of use of a model.

`select` applies the declaration (`Credibility`) to the table of the adopted fits
(`FitTable`): the points of every evidence, the limitation points and the
validation points which no evidence claims. Everything is keyed by the
experiment and the key of a mapping, never by the key alone, which two studies
can share. The rules:

- A mapping which is training data of one adopted fit is fitted data in every
  fit: its points in another fit, where it may be validation data, are no
  validation either.

- The quantities of a context of use restrict its evidence to their observables,
  compared without the brackets of a concentration (`[A]` is `A`).
- The limitations, of the model and of a context, are subtracted from every
  other evidence and counted; a limitation of fitted data is an error.
- A validation of fitted data is an error, its outliers are dropped and counted.
- A calibration keeps the fitted points of its selection.
- `Unclaimed` takes the validation points which no limitation and no other
  evidence of any context selects; two `Unclaimed` of different contexts may
  take the same points.
- An evidence which selects no point is an error, and so is an evidence of a
  fit the model has not adopted.
"""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pkpdsim.credibility.declaration import (
    ContextOfUse,
    Credibility,
    Evidence,
    EvidenceRole,
    Unclaimed,
)

if TYPE_CHECKING:
    from sbmlsim.fit.helpers import MappingFilter

    from pkpdsim.credibility.records import FitTable

__all__ = [
    "EvidenceSelection",
    "SelectionError",
    "Selections",
    "normalize_observable",
    "select",
]

type _Mask = np.ndarray[tuple[int], np.dtype[np.bool_]]


def normalize_observable(observable: str) -> str:
    """Get an observable without the brackets of a concentration.

    Args:
        observable: the observable, e.g. `[A]`, `A` or `pk.cmax`.

    Returns:
        The observable, e.g. `A` for `[A]`.
    """
    if observable.startswith("[") and observable.endswith("]"):
        return observable[1:-1]
    return observable


@dataclass(frozen=True, kw_only=True, eq=False)
class EvidenceSelection:
    """The data points of an evidence of a context of use.

    Attributes:
        context: the id of the context of use.
        evidence: the evidence.
        rows: the selected rows of the table of the fits (`FitTable.table`).
        outliers_dropped: the points of kind outlier which a validation
            selected and dropped.
        limitation_subtracted: the limitation points the evidence selected and
            which were subtracted.
    """

    context: str
    evidence: Evidence
    rows: pd.Index
    outliers_dropped: int
    limitation_subtracted: int


@dataclass(frozen=True, kw_only=True, eq=False)
class Selections:
    """The data points of the evidence of every context of use of a model.

    Attributes:
        evidence: the selection of every evidence, in the order of the contexts
            and their evidence.
        limitations: the rows of the limitations of the model and of every
            context.
        unclaimed: the validation points of no evidence and no limitation, with
            the columns `fit`, `experiment`, `mapping` and `n` (the points).
        limitation_evidence: the selection of every limitation, the ones of the
            model first, with the context `""`, then the ones of the contexts.
    """

    evidence: tuple[EvidenceSelection, ...]
    limitations: pd.Index
    unclaimed: pd.DataFrame
    limitation_evidence: tuple[EvidenceSelection, ...] = ()


class SelectionError(ValueError):
    """The evidence of a declaration selects invalid data.

    Attributes:
        problems: every problem of the selection.
    """

    def __init__(self, problems: list[str]) -> None:
        """Create the error.

        Args:
            problems: every problem of the selection.
        """
        self.problems = list(problems)
        super().__init__(
            "The evidence of the credibility declaration selects invalid data:\n"
            + "\n".join(f"- {problem}" for problem in self.problems)
        )


class _Selector:
    """The masks of the rows of a table of the fits, collecting the problems."""

    def __init__(self, fits: FitTable) -> None:
        """Index the table.

        Args:
            fits: the table of the fits.
        """
        table = fits.table
        self.fits = fits
        self.problems: list[str] = []
        self.size = len(table)
        self.fit = table["fit"].astype(str).to_numpy()
        self.subsets = sorted(set(self.fit))
        self.experiments = table["experiment"].astype(str).to_numpy()
        self.mappings = table["mapping"].astype(str).to_numpy()
        self.observables = (
            table["observable"].astype(str).map(normalize_observable).to_numpy()
        )
        kinds = table["kind"].astype(str).to_numpy()
        self.training: _Mask = kinds == "training"
        self.validation: _Mask = kinds == "validation"
        self.outlier: _Mask = kinds == "outlier"
        self.pairs = list(zip(self.experiments, self.mappings, strict=True))
        self.unique_pairs = list(dict.fromkeys(self.pairs))
        # the fits a mapping is training data of, it is fitted in every fit
        self.trained_in: dict[tuple[str, str], set[str]] = {}
        for pair, fit, training in zip(
            self.pairs, self.fit, self.training, strict=True
        ):
            if training:
                self.trained_in.setdefault(pair, set()).add(str(fit))
        self.fitted: _Mask = np.fromiter(
            (p in self.trained_in for p in self.pairs), bool, self.size
        )

    def none(self) -> _Mask:
        """Get the mask of no row."""
        return np.zeros(self.size, dtype=bool)

    def rows(self, mask: _Mask) -> pd.Index:
        """Get the rows of a mask as the index of the table."""
        return self.fits.table.index[mask]

    def describe(self, mask: _Mask) -> str:
        """Get the mappings of the rows of a mask, e.g. `CredA.fm_dose (DECAY)`."""
        names = {
            f"{experiment}.{mapping} ({fit})"
            for experiment, mapping, fit in zip(
                self.experiments[mask], self.mappings[mask], self.fit[mask], strict=True
            )
        }
        return ", ".join(sorted(names))

    def describe_fitted(self, mask: _Mask) -> str:
        """Get the fitted mappings of the rows of a mask with the fits they are trained in.

        E.g. `CredB.fm_dose (PD, training in PK)`.
        """
        names = {
            f"{experiment}.{mapping} ({fit}, training in "
            f"{', '.join(sorted(self.trained_in[(experiment, mapping)]))})"
            for experiment, mapping, fit in zip(
                self.experiments[mask], self.mappings[mask], self.fit[mask], strict=True
            )
        }
        return ", ".join(sorted(names))

    def accepted(self, selection: MappingFilter) -> _Mask:
        """Get the rows of the mappings a selection accepts, one call per mapping."""
        keep: set[tuple[str, str]] = {
            (experiment, mapping)
            for experiment, mapping in self.unique_pairs
            if selection(mapping, self.fits.lookup(experiment, mapping))
        }
        return np.fromiter((p in keep for p in self.pairs), bool, self.size)

    def of_fit(self, evidence: Evidence, who: str) -> _Mask | None:
        """Get the rows of the fit of an evidence, `None` for an unknown fit."""
        if evidence.fit is None:
            return np.ones(self.size, dtype=bool)
        if evidence.fit not in self.subsets:
            self.problems.append(
                f"{who}: the fit '{evidence.fit}' is no subset of the adopted "
                f"fits, the subsets are {self.subsets}"
            )
            return None
        return self.fit == evidence.fit

    def searched(self, evidence: Evidence) -> list[str]:
        """Get the subsets of the fits an evidence is searched in."""
        return self.subsets if evidence.fit is None else [evidence.fit]

    def quantities(self, cou: ContextOfUse) -> _Mask:
        """Get the rows of the quantities of a context, a problem per unknown one."""
        known = set(self.observables)
        ids = {normalize_observable(q.id) for q in cou.quantities}
        self.problems.extend(
            f"context of use '{cou.id}': the quantity '{q.id}' matches no "
            f"observable of the fits, the observables are {sorted(known)}"
            for q in cou.quantities
            if normalize_observable(q.id) not in known
        )
        return np.isin(self.observables, sorted(ids))

    def base(self, evidence: Evidence, who: str, quantities: _Mask) -> _Mask | None:
        """Get the rows of the mappings an evidence accepts.

        Restricted to the quantities and to the fit of the evidence; `None` for
        an evidence of an unknown fit.
        """
        fit = self.of_fit(evidence, who)
        if fit is None:
            return None
        if isinstance(evidence.selection, Unclaimed):
            raise TypeError(f"{who}: an unclaimed evidence has no selection")
        return self.accepted(evidence.selection) & quantities & fit

    def limitation(self, evidence: Evidence, who: str) -> _Mask:
        """Get the rows of a limitation, a problem for fitted data or none."""
        if isinstance(evidence.selection, Unclaimed):
            self.problems.append(
                f"{who}: a limitation selects mappings, it cannot be unclaimed"
            )
            return self.none()
        mask = self.base(evidence, who, np.ones(self.size, dtype=bool))
        if mask is None:
            return self.none()
        if not mask.any():
            self.problems.append(
                f"{who}: the limitation selects no point of the fits "
                f"{self.searched(evidence)}"
            )
        fitted = mask & self.fitted
        if fitted.any():
            self.problems.append(
                f"{who}: a limitation is data the model does not describe, but "
                f"it selects fitted data: {self.describe_fitted(fitted)}"
            )
        return mask

    def empty(self, mask: _Mask, cou: ContextOfUse, evidence: Evidence) -> None:
        """Add a problem if an evidence selects no point."""
        if not mask.any():
            self.problems.append(
                f"context of use '{cou.id}': the evidence '{evidence.name}' "
                f"selects no point of the fits {self.searched(evidence)}"
            )


def select(credibility: Credibility, fits: FitTable) -> Selections:
    """Select the data points of the evidence of every context of use.

    A mapping (its experiment and key) which is training data of any adopted fit
    counts as fitted in every fit: a validation which selects its points, in
    any fit, is an error which names the fits it is training data of, and
    neither `Unclaimed` nor the unclaimed validation points take them. Two
    `Unclaimed` of different contexts may take the same points. See the module
    for the other rules.

    Args:
        credibility: the declaration.
        fits: the table of the adopted fits, see `evaluate_fits` and
            `mapping_table`.

    Returns:
        The points of every evidence, the limitation points and the unclaimed
        validation points.

    Raises:
        SelectionError: with every problem of the selection.
    """
    selector = _Selector(fits)
    limitation_evidence: list[EvidenceSelection] = []

    def limitation(context: str, evidence: Evidence, who: str) -> _Mask:
        """Get the rows of a limitation and keep its selection."""
        mask = selector.limitation(evidence, who)
        limitation_evidence.append(
            EvidenceSelection(
                context=context,
                evidence=evidence,
                rows=selector.rows(mask),
                outliers_dropped=0,
                limitation_subtracted=0,
            )
        )
        return mask

    model_limitations = selector.none()
    for evidence in credibility.limitations:
        model_limitations |= limitation(
            "", evidence, f"model limitation '{evidence.name}'"
        )
    limitations: dict[str, _Mask] = {}
    quantities: dict[str, _Mask] = {}
    for cou in credibility.contexts:
        mask = model_limitations.copy()
        for evidence in cou.limitations:
            mask |= limitation(
                cou.id,
                evidence,
                f"context of use '{cou.id}': the limitation '{evidence.name}'",
            )
        limitations[cou.id] = mask
        quantities[cou.id] = selector.quantities(cou)
    all_limitations = model_limitations.copy()
    for mask in limitations.values():
        all_limitations |= mask

    # the evidence with a selection first, `Unclaimed` takes what it leaves
    selected: dict[tuple[int, int], EvidenceSelection] = {}
    claimed = selector.none()
    for i, cou in enumerate(credibility.contexts):
        for j, evidence in enumerate(cou.evidence):
            if isinstance(evidence.selection, Unclaimed):
                continue
            who = f"context of use '{cou.id}': the evidence '{evidence.name}'"
            if evidence.role is EvidenceRole.LIMITATION:
                selector.problems.append(
                    f"{who} has the role limitation, a limitation of a context "
                    f"is one of its limitations"
                )
                continue
            base = selector.base(evidence, who, quantities[cou.id])
            if base is None:
                continue
            subtracted = base & limitations[cou.id]
            rows = base & ~limitations[cou.id]
            outliers = selector.none()
            fitted = selector.none()
            if evidence.role is EvidenceRole.VALIDATION:
                fitted = rows & selector.fitted
                if fitted.any():
                    selector.problems.append(
                        f"{who} validates on fitted data: "
                        f"{selector.describe_fitted(fitted)}"
                    )
                outliers = rows & selector.outlier & ~selector.fitted
                rows = rows & selector.validation & ~selector.fitted
            else:
                rows = rows & selector.training
            if not fitted.any():
                # a validation of fitted data alone is one problem, not two
                selector.empty(rows, cou, evidence)
            claimed |= rows
            selected[(i, j)] = EvidenceSelection(
                context=cou.id,
                evidence=evidence,
                rows=selector.rows(rows),
                outliers_dropped=int(outliers.sum()),
                limitation_subtracted=int(subtracted.sum()),
            )

    # the validation points of a mapping which another fit fitted are fitted
    unclaimed = selector.validation & ~selector.fitted & ~all_limitations & ~claimed
    for i, cou in enumerate(credibility.contexts):
        for j, evidence in enumerate(cou.evidence):
            if not isinstance(evidence.selection, Unclaimed):
                continue
            who = f"context of use '{cou.id}': the evidence '{evidence.name}'"
            fit = selector.of_fit(evidence, who)
            if fit is None:
                continue
            base = quantities[cou.id] & selector.validation & ~selector.fitted & fit
            within = evidence.selection.within
            if within is not None:
                base &= selector.accepted(within)
            rows = base & ~all_limitations & ~claimed
            selector.empty(rows, cou, evidence)
            unclaimed &= ~rows
            selected[(i, j)] = EvidenceSelection(
                context=cou.id,
                evidence=evidence,
                rows=selector.rows(rows),
                outliers_dropped=0,
                limitation_subtracted=int((base & all_limitations).sum()),
            )

    if selector.problems:
        raise SelectionError(selector.problems)
    return Selections(
        evidence=tuple(selected[key] for key in sorted(selected)),
        limitations=selector.rows(all_limitations),
        unclaimed=_count(fits.table, unclaimed),
        limitation_evidence=tuple(limitation_evidence),
    )


def _count(table: pd.DataFrame, mask: _Mask) -> pd.DataFrame:
    """Get the points of the rows of a mask per fit and mapping."""
    keys: list[Hashable] = ["fit", "experiment", "mapping"]
    return (
        table.loc[mask]
        .groupby(keys, sort=True)
        .size()
        .rename("n")
        .reset_index()
        .astype({"n": int})
    )
