"""The ratios between the arms of a study and their limits of Guest et al. (2011).

A ratio quantity (`QuantityOfInterest(ratio=Ratio(...))`) asks how an arm
differs from its reference arm in the same study: an impaired over a healthy
group, a genotype over the wild type, a drug with its perpetrator over the drug
alone. A fold criterion on the absolute predictions is too lenient for such a
question near a ratio of 1: two predictions within 2-fold each may give a ratio
which is 4-fold off. The ratio is therefore compared directly:

- `ratio_pairs` forms the pairs of the decisive points of a context of use:
  per fit and study (the experiment) every mapping of the observable which the
  numerator selects with the one mapping of the study which the denominator
  selects; a study without a denominator gives no pair, a study with more than
  one is a problem of the declaration.
- The summary of a mapping is the same for the data and the predictions
  (`Summary`): the area under the points over their times by the trapezoid
  (`auc`, at least two points with a finite data value and prediction) or the
  value of its single point (`value`).
- `guest_limit` is the limit of Guest et al. (2011) for an observed ratio,
  `(delta + 2 (R - 1)) / R` with `R` the ratio or its inverse when below 1 and
  `delta` 1.25: a prediction is within when it lies between the observed ratio
  divided by and times the limit (1.25 at a ratio of 1, towards 2 for a large
  ratio).

Guest, E.J., Aarons, L., Houston, J.B., Rostami-Hodjegan, A., Galetin, A.
(2011). Critique of the two-fold measure of prediction success for ratios:
application for the assessment of drug-drug interactions. Drug Metab Dispos
39(2), 170-173.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Final, TypedDict

import numpy as np
import pandas as pd

from pkpdsim.credibility.declaration import GUEST_DELTA, Ratio, Summary
from pkpdsim.credibility.selection import normalize_observable

if TYPE_CHECKING:
    from pkpdsim.credibility.records import MappingLookup

__all__ = [
    "RATIO_COLUMNS",
    "Pairing",
    "RatioRecord",
    "guest_limit",
    "pairings",
    "ratio_pairs",
    "summarize",
]

#: the columns of the table of the ratio pairs, see `RatioRecord`
RATIO_COLUMNS: Final = (
    "quantity",
    "fit",
    "experiment",
    "numerator",
    "denominator",
    "observed",
    "predicted",
    "limit",
    "within",
)


class RatioRecord(TypedDict):
    """A pair of arms of a ratio quantity.

    Attributes:
        quantity: the id of the ratio quantity.
        fit: the subset of the fit of the points.
        experiment: the experiment, the study of the arms.
        numerator: the key of the mapping of the numerator.
        denominator: the key of the mapping of the denominator.
        observed: the ratio of the summaries of the data, `NaN` if one is not
            evaluable.
        predicted: the ratio of the summaries of the predictions.
        limit: the limit of Guest et al. for the observed ratio.
        within: whether the predicted ratio lies within the limits, `None` if
            the pair is not evaluable.
    """

    quantity: str
    fit: str
    experiment: str
    numerator: str
    denominator: str
    observed: float
    predicted: float
    limit: float
    within: bool | None


def guest_limit(observed: float, delta: float = GUEST_DELTA) -> float:
    """Get the limit of Guest et al. (2011) of an observed ratio.

    Args:
        observed: the observed ratio, positive.
        delta: the accepted variability at a ratio of 1.

    Returns:
        The limit, `NaN` for a ratio which is not positive and finite.
    """
    if not math.isfinite(observed) or observed <= 0:
        return math.nan
    ratio = observed if observed >= 1 else 1 / observed
    return (delta + 2 * (ratio - 1)) / ratio


def summarize(points: pd.DataFrame, summary: Summary) -> tuple[float, float]:
    """Get the summary of the data and of the predictions of the points of a mapping.

    Args:
        points: the points of one mapping with `DV`, `IPRED` and, for the area,
            the times `x`.
        summary: how the points are summarized.

    Returns:
        The summary of `DV` and of `IPRED`, `NaN` if not evaluable: the area
        needs two points with a finite data value and prediction, the value
        exactly one point.
    """
    dv = points["DV"].to_numpy(dtype=float)
    ipred = points["IPRED"].to_numpy(dtype=float)
    finite = np.isfinite(dv) & np.isfinite(ipred)
    if summary is Summary.VALUE:
        if len(points) != 1 or not finite.all():
            return math.nan, math.nan
        return float(dv[0]), float(ipred[0])
    if "x" not in points.columns:
        return math.nan, math.nan
    x = points["x"].to_numpy(dtype=float)
    finite &= np.isfinite(x)
    if finite.sum() < 2:
        return math.nan, math.nan
    order = np.argsort(x[finite], kind="stable")
    times = x[finite][order]
    return (
        float(np.trapezoid(dv[finite][order], times)),
        float(np.trapezoid(ipred[finite][order], times)),
    )


def _ratio(numerator: float, denominator: float) -> float:
    """Get a ratio of two summaries, `NaN` unless both are positive and finite."""
    if not (math.isfinite(numerator) and math.isfinite(denominator)):
        return math.nan
    if numerator <= 0 or denominator <= 0:
        return math.nan
    return numerator / denominator


@dataclass(frozen=True, kw_only=True)
class Pairing:
    """The arms of a ratio in a study: its numerators over its denominator.

    Attributes:
        fit: the subset of the fit of the points.
        experiment: the experiment, the study of the arms.
        numerators: the keys of the mappings of the numerator.
        denominator: the key of the mapping of the denominator.
    """

    fit: str
    experiment: str
    numerators: tuple[str, ...]
    denominator: str


def pairings(
    quantity: str, ratio: Ratio, points: pd.DataFrame, lookup: MappingLookup
) -> tuple[list[Pairing], list[str]]:
    """Get the arms of a ratio quantity in the points, per fit and study.

    Args:
        quantity: the id of the ratio quantity.
        ratio: its pairing.
        points: points or mappings with `fit`, `experiment`, `mapping` and
            `observable`.
        lookup: the fit mapping of an experiment and a key, which the filters
            of the numerator and the denominator read.

    Returns:
        The arms of every study with a numerator and a denominator, and the
        problems: a study with more than one mapping the denominator selects.
    """
    observable = normalize_observable(ratio.observable)
    rows = points[
        points["observable"].astype(str).map(normalize_observable) == observable
    ]
    found: list[Pairing] = []
    problems: list[str] = []
    for (fit, experiment), study in rows.groupby(["fit", "experiment"], sort=False):
        keys = list(dict.fromkeys(study["mapping"].astype(str)))
        denominators = [
            k for k in keys if ratio.denominator(k, lookup(str(experiment), k))
        ]
        numerators = tuple(
            k
            for k in keys
            if k not in denominators and ratio.numerator(k, lookup(str(experiment), k))
        )
        if not numerators or not denominators:
            continue
        if len(denominators) > 1:
            problems.append(
                f"the ratio '{quantity}' has more than one denominator in "
                f"{experiment} ({fit}): {', '.join(denominators)}"
            )
            continue
        found.append(
            Pairing(
                fit=str(fit),
                experiment=str(experiment),
                numerators=numerators,
                denominator=denominators[0],
            )
        )
    return found, problems


def ratio_pairs(
    quantity: str,
    ratio: Ratio,
    points: pd.DataFrame,
    lookup: MappingLookup,
    delta: float = GUEST_DELTA,
) -> tuple[list[RatioRecord], list[str]]:
    """Get the pairs of the arms of a ratio quantity in the decisive points.

    Args:
        quantity: the id of the ratio quantity.
        ratio: its pairing.
        points: the decisive points of the context of use with `fit`,
            `experiment`, `mapping`, `observable`, `DV`, `IPRED` and `x`.
        lookup: the fit mapping of an experiment and a key, which the filters
            of the numerator and the denominator read.
        delta: the accepted variability of the limit at a ratio of 1.

    Returns:
        The pairs, in the order of the points, and the problems of
        `pairings`.
    """
    found, problems = pairings(quantity, ratio, points, lookup)
    pairs: list[RatioRecord] = []
    for arms in found:
        study = points[
            (points["fit"].astype(str) == arms.fit)
            & (points["experiment"].astype(str) == arms.experiment)
        ]
        dv_den, ipred_den = summarize(
            study[study["mapping"] == arms.denominator], ratio.summary
        )
        for numerator in arms.numerators:
            dv_num, ipred_num = summarize(
                study[study["mapping"] == numerator], ratio.summary
            )
            observed = _ratio(dv_num, dv_den)
            predicted = _ratio(ipred_num, ipred_den)
            limit = guest_limit(observed, delta)
            evaluable = math.isfinite(limit) and math.isfinite(predicted)
            pairs.append(
                RatioRecord(
                    quantity=quantity,
                    fit=arms.fit,
                    experiment=arms.experiment,
                    numerator=numerator,
                    denominator=arms.denominator,
                    observed=observed,
                    predicted=predicted,
                    limit=limit,
                    within=(
                        bool(observed / limit <= predicted <= observed * limit)
                        if evaluable
                        else None
                    ),
                )
            )
    return pairs, problems


def ratio_table(pairs: list[RatioRecord]) -> pd.DataFrame:
    """Get the pairs as a table of `RATIO_COLUMNS`."""
    return pd.DataFrame(pairs, columns=list(RATIO_COLUMNS))
