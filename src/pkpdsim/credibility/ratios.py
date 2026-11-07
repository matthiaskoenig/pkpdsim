"""The ratios between the arms of a study and their limits of Guest et al. (2011).

A ratio quantity (`QuantityOfInterest(ratio=Ratio(...))` of
`sbmlsim.credibility`) asks how an arm differs from its reference arm in the
same study: an impaired over a healthy group, a genotype over the wild type, a
drug with its perpetrator over the drug alone. A fold criterion on the absolute
predictions is too lenient for such a question near a ratio of 1: two
predictions within 2-fold each may give a ratio which is 4-fold off. The ratio
is therefore compared directly:

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

The pairs are the evidence extension `RATIOS` (`RatiosExtension`, name
`ratios`) of the credibility of sbmlsim: it claims the `Ratio` pairings, its
rule (`GuestAcceptance`, the setting `ContextOfUse.extensions["ratios"]`, else
`DEFAULT_RATIO_ACCEPTANCE` of the risk) enters the verdict, the pairs cover the
ratio quantities, the factor relevance of the quantities states them, the
record of a context has the pairs in its section and the scorecard pools them
by a value of the metadata of their numerator (`RatiosExtension.pool`, the
qualification of a scaling across drugs).

Guest, E.J., Aarons, L., Houston, J.B., Rostami-Hodjegan, A., Galetin, A.
(2011). Critique of the two-fold measure of prediction success for ratios:
application for the assessment of drug-drug interactions. Drug Metab Dispos
39(2), 170-173.
"""

from __future__ import annotations

import math
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum, unique
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final, TypedDict, override

import numpy as np
import pandas as pd
from sbmlsim.credibility.assessment import applicable
from sbmlsim.credibility.declaration import (
    AcceptanceResult,
    Level,
    acceptance_result,
)
from sbmlsim.credibility.extensions import (
    ACCEPTANCE_EXTENSION,
    ContextBlock,
    ContextContribution,
    EvidenceExtension,
    FactorStatement,
)
from sbmlsim.credibility.selection import (
    SelectionError,
    in_domain,
    normalize_observable,
)

from pkpdsim.credibility.blocks import render_macro

if TYPE_CHECKING:
    from sbmlsim.credibility.declaration import ContextOfUse, Credibility
    from sbmlsim.credibility.extensions import Pairing as CorePairing
    from sbmlsim.credibility.records import FitTable, MappingLookup
    from sbmlsim.fit.helpers import MappingFilter
    from sbmlsim.project import Project

__all__ = [
    "DEFAULT_RATIO_ACCEPTANCE",
    "GUEST_DELTA",
    "QUALIFICATION_COLUMNS",
    "RATIOS",
    "RATIO_COLUMNS",
    "GuestAcceptance",
    "Pairing",
    "Ratio",
    "RatioRecord",
    "RatiosExtension",
    "Summary",
    "guest_limit",
    "pairings",
    "ratio_pairs",
    "ratio_rule",
    "section_pairs",
    "summarize",
]

#: the accepted variability of a ratio of 1 of the limits of Guest et al. (2011)
GUEST_DELTA: Final = 1.25

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

#: the columns of the qualification, see `RatiosExtension.pool`
QUALIFICATION_COLUMNS: Final = (
    "field",
    "value",
    "model",
    "n_pairs",
    "n_within",
    "within",
    "gmfe",
)

#: the field of the metadata the qualification pools by without an option
DEFAULT_FIELD: Final = "health"


@unique
class Summary(StrEnum):
    """How the points of an arm of a ratio are summarized, see `Ratio`."""

    AUC = "auc"
    VALUE = "value"


@dataclass(frozen=True, kw_only=True)
class Ratio:
    """The arms of a ratio quantity: a numerator over its reference in a study.

    The pairing of a ratio quantity (`sbmlsim.credibility.extensions.Pairing`),
    which `RatiosExtension` claims: per fit and study (the experiment) every
    decisive mapping of the observable which `numerator` selects is paired with
    the one which `denominator` selects, and the ratio of their summaries is
    compared for the data and the predictions.

    Attributes:
        observable: the observable of the arms, e.g. `Carm_hctz`; compared
            without the brackets of a concentration.
        numerator: selects the mappings of the arms, e.g. the renal impairment.
        denominator: selects the mapping of the reference arm, e.g. healthy.
        summary: how the points of an arm are summarized.
    """

    observable: str
    numerator: MappingFilter
    denominator: MappingFilter
    summary: Summary = Summary.AUC


@dataclass(frozen=True)
class GuestAcceptance:
    """The percentage of the ratio pairs within the limits of Guest et al. is at least `minimum`.

    The setting of the extension `ratios` of a context of use,
    `ContextOfUse(extensions={"ratios": GuestAcceptance(80)})`.

    Attributes:
        minimum: the minimum percentage (0 to 100).
        delta: the accepted variability at a ratio of 1, see `guest_limit`.
    """

    minimum: float
    delta: float = GUEST_DELTA

    @property
    def name(self) -> str:
        """Get the name of the rule."""
        return "Guest"

    @property
    def limit(self) -> str:
        """Get the limit, e.g. `>= 80`."""
        return f">= {self.minimum:g}"

    def evaluate(self, pairs: Collection[bool | None]) -> AcceptanceResult:
        """Evaluate the percentage of the evaluable pairs within the limits.

        Args:
            pairs: whether every pair is within, `None` for one which is not
                evaluable.

        Returns:
            The result, not evaluable without an evaluable pair.
        """
        evaluable = [within for within in pairs if within is not None]
        value = 100.0 * sum(evaluable) / len(evaluable) if evaluable else math.nan
        return acceptance_result(self.name, value, self.limit, value >= self.minimum)

    def problems(self) -> list[str]:
        """Get the problems of the rule: a minimum outside 0 to 100, a delta not above 1."""
        found: list[str] = []
        if not 0 <= self.minimum <= 100:
            found.append(
                f"the minimum of the ratio acceptance must be between 0 and 100, "
                f"not {self.minimum:g}"
            )
        if not self.delta > 1:
            found.append(
                f"the delta of the ratio acceptance must be above 1, not {self.delta:g}"
            )
        return found


#: the rule of the ratio pairs of a context of use which declares none, by risk
DEFAULT_RATIO_ACCEPTANCE: Final[Mapping[Level, GuestAcceptance]] = MappingProxyType(
    {
        Level.LOW: GuestAcceptance(50.0),
        Level.MEDIUM: GuestAcceptance(67.0),
        Level.HIGH: GuestAcceptance(80.0),
    }
)


def ratio_rule(cou: ContextOfUse) -> GuestAcceptance:
    """Get the rule of the ratio pairs of a context: the declared or the default of the risk.

    Args:
        cou: the context of use.

    Returns:
        The setting `extensions["ratios"]` of the context, else the default of
        its risk.

    Raises:
        TypeError: for a setting which is no `GuestAcceptance`.
    """
    rule = cou.extensions.get(RatiosExtension.name)
    if rule is None:
        return DEFAULT_RATIO_ACCEPTANCE[cou.risk]
    if not isinstance(rule, GuestAcceptance):
        raise TypeError(
            f"The ratio acceptance of the context of use '{cou.id}' is a "
            f"GuestAcceptance, not a {type(rule).__name__}."
        )
    return rule


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
        metadata: the values of the fields of the domain of the numerator, by
            name, e.g. `{"health": "RENAL_IMPAIRMENT"}`, which the
            qualification across models groups by.
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
    metadata: dict[str, str]


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
    fields: Sequence[str] | None = None,
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
        fields: the fields of the metadata of the domain the pairs record of
            their numerator, `pkpdsim.model.DOMAIN_FIELDS` without.

    Returns:
        The pairs, in the order of the points, and the problems of
        `pairings`.
    """
    if fields is None:
        from pkpdsim.model import DOMAIN_FIELDS

        fields = DOMAIN_FIELDS
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
            arm = study[study["mapping"] == numerator]
            dv_num, ipred_num = summarize(arm, ratio.summary)
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
                    metadata={
                        field: str(arm[field].iloc[0])
                        for field in fields
                        if field in arm.columns and pd.notna(arm[field].iloc[0])
                    },
                )
            )
    return pairs, problems


def _write(value: float) -> float | None:
    """Get a number for JSON, `None` for `NaN`."""
    return None if math.isnan(value) else value


def _read(value: float | None) -> float:
    """Get a number from JSON, `NaN` for `None`."""
    return math.nan if value is None else float(value)


def _to_dict(pair: RatioRecord) -> dict[str, Any]:
    """Get a pair as plain data, `NaN` as `None`."""
    return {
        **pair,
        "observed": _write(pair["observed"]),
        "predicted": _write(pair["predicted"]),
        "limit": _write(pair["limit"]),
        "metadata": dict(pair["metadata"]),
    }


def _from_dict(data: Mapping[str, Any]) -> RatioRecord:
    """Get a pair from plain data, `None` as `NaN`."""
    return RatioRecord(
        quantity=data["quantity"],
        fit=data["fit"],
        experiment=data["experiment"],
        numerator=data["numerator"],
        denominator=data["denominator"],
        observed=_read(data["observed"]),
        predicted=_read(data["predicted"]),
        limit=_read(data["limit"]),
        within=data["within"],
        metadata=dict(data.get("metadata", {})),
    )


def section_pairs(section: Mapping[str, Any] | None) -> list[RatioRecord]:
    """Get the ratio pairs of the section `ratios` of the record of a context.

    Args:
        section: the section, `{"ratios": [...]}`, also of a record of pkpdsim
            0.19; `None` without.

    Returns:
        The pairs, `NaN` for a missing number.
    """
    if section is None:
        return []
    return [_from_dict(r) for r in section.get("ratios", ())]


def _count(n: int, singular: str, plural: str | None = None) -> str:
    """Get a number of things in words, e.g. `1 point` and `2 points`."""
    return f"{n} {singular if n == 1 else plural or singular + 's'}"


class RatiosExtension(EvidenceExtension):
    """The ratio pairs of the ratio quantities, see the module.

    Prepared: the fields of the domain of the project, which the pairs record
    of their numerator.
    """

    name = "ratios"
    label = "ratio pairs"

    @override
    def claims(self, pairing: CorePairing) -> bool:
        """Claim the pairings of the arms of a study, `Ratio`."""
        return isinstance(pairing, Ratio)

    def _pairs(
        self,
        cou: ContextOfUse,
        decisive: pd.DataFrame,
        lookup: MappingLookup,
        fields: Sequence[str] | None,
    ) -> list[RatioRecord]:
        """Get the pairs of the assessed ratio quantities of a context.

        Raises:
            SelectionError: for a study with more than one denominator, which
                the selection refuses first.
        """
        rule = ratio_rule(cou)
        pairs: list[RatioRecord] = []
        for q in cou.ratios:
            if not isinstance(q.pairing, Ratio):
                continue
            found, problems = ratio_pairs(
                q.id, q.pairing, decisive, lookup, rule.delta, fields
            )
            if problems:
                raise SelectionError(
                    [f"context of use '{cou.id}': {p}" for p in problems]
                )
            pairs.extend(found)
        return pairs

    @override
    def selection_problems(
        self, cou: ContextOfUse, decisive: pd.DataFrame, lookup: MappingLookup
    ) -> list[str]:
        """Report every study of a ratio with more than one denominator."""
        problems: list[str] = []
        for q in cou.ratios:
            if isinstance(q.pairing, Ratio):
                problems += pairings(q.id, q.pairing, decisive, lookup)[1]
        return problems

    @override
    def prepare(
        self, project: Project, credibility: Credibility, fits: FitTable
    ) -> object:
        """Get the fields of the domain of the project."""
        return tuple(project.domain_fields())

    @override
    def assess_context(
        self,
        cou: ContextOfUse,
        decisive: pd.DataFrame,
        fits: FitTable,
        prepared: object,
    ) -> ContextContribution | None:
        """Assess the pairs of the ratio quantities of a context, none without."""
        ratios = [q for q in cou.ratios if isinstance(q.pairing, Ratio)]
        if not ratios:
            return None
        fields = [str(f) for f in prepared] if isinstance(prepared, tuple) else None
        rule = ratio_rule(cou)
        pairs = self._pairs(cou, decisive, fits.lookup, fields)
        result = rule.evaluate([p["within"] for p in pairs])
        inside = in_domain(cou, applicable(decisive))
        claimed = set(
            zip(
                inside["fit"].astype(str),
                inside["experiment"].astype(str),
                inside["mapping"].astype(str),
                strict=True,
            )
        )
        quantities = {
            q.id: sum(
                p["quantity"] == q.id
                and p["within"] is not None
                and (p["fit"], p["experiment"], p["numerator"]) in claimed
                for p in pairs
            )
            for q in ratios
        }
        evaluable = [p for p in pairs if p["within"] is not None]
        within = sum(p["within"] is True for p in evaluable)
        covered = sum(n > 0 for n in quantities.values())
        statement = (
            f"{within} of {_count(len(evaluable), 'ratio pair')} within the "
            f"limits of Guest et al. ({rule.limit} %), {covered} of "
            f"{_count(len(ratios), 'ratio quantity', 'ratio quantities')} covered."
        )
        return ContextContribution(
            section={"ratios": [_to_dict(p) for p in pairs]},
            acceptance=(result,),
            quantities=quantities,
            factors={"relevance_of_quantities": FactorStatement(statement)},
        )

    @override
    def acceptance(
        self, cou: ContextOfUse, decisive: pd.DataFrame, lookup: MappingLookup
    ) -> list[AcceptanceResult]:
        """Evaluate the rule of the ratio pairs on other predictions."""
        pairs = self._pairs(cou, decisive, lookup, None)
        return [ratio_rule(cou).evaluate([p["within"] for p in pairs])]

    @override
    def context_report_blocks(
        self, context: Mapping[str, Any], section: Mapping[str, Any] | None
    ) -> Sequence[ContextBlock]:
        """Show the ratio pairs after the acceptance of a context with the rule."""
        rows = context.get("acceptance", ())
        if not any(row.get(ACCEPTANCE_EXTENSION) == self.name for row in rows):
            return ()
        html = render_macro("ratio_pairs", section_pairs(section))
        return (ContextBlock(html, after="acceptance"),)

    @override
    def report_tables(self, record: Mapping[str, Any]) -> Mapping[str, pd.DataFrame]:
        """Write `ratios.tsv`: one row per context of use and pair of a ratio quantity."""
        rows = [
            {"context": c["id"], **pair}
            for c in record["contexts"]
            for pair in section_pairs(c.get("extensions", {}).get(self.name))
        ]
        return {"ratios": pd.DataFrame(rows, columns=["context", *RATIO_COLUMNS])}

    @override
    @override
    def pool_lines(
        self, table: pd.DataFrame, options: Mapping[str, Any]
    ) -> Sequence[str]:
        """Word the pooled rows (the model `all`) as the share within the limits of Guest."""
        pooled = table[table["model"] == "all"]
        return [
            f"  {row['value']}: {row['n_within']} of {row['n_pairs']} ratio pairs "
            f"within the limits of Guest ({row['within']:.0f} %), GMFE "
            f"{row['gmfe']:.3g}"
            for row in pooled.to_dict("records")
        ]

    def pool(
        self,
        contexts: Mapping[str, Sequence[Mapping[str, Any] | None]],
        options: Mapping[str, Any],
    ) -> pd.DataFrame | None:
        """Pool the ratio pairs of the models by a value of the metadata of their numerator.

        Every evaluable pair of every context of use counts once per study and
        arms: a pair which two contexts of a model share counts once; a pair
        whose numerator has the value `NR` (not reported) counts for no value.

        Args:
            contexts: the sections of the contexts of use, by model.
            options: `field`, the field of the metadata of the numerator,
                `health` without.

        Returns:
            The columns of `QUALIFICATION_COLUMNS`: per value of the field one
            row per model and one pooled row with the model `all`: the
            evaluable pairs, the pairs within the limits of Guest et al., their
            percentage and the geometric mean fold error of the predicted over
            the observed ratios; sorted by value with the pooled row first.
        """
        field = str(options.get("field", DEFAULT_FIELD))
        pairs: dict[tuple[str, str, str, str, str], tuple[str, float, float, bool]] = {}
        for model, sections in contexts.items():
            for section in sections:
                for r in section_pairs(section):
                    value = r["metadata"].get(field)
                    if r["within"] is None or value in (None, "NR"):
                        # a value which is not reported qualifies nothing
                        continue
                    key = (
                        model,
                        r["fit"],
                        r["experiment"],
                        r["numerator"],
                        r["denominator"],
                    )
                    pairs[key] = (
                        str(value),
                        r["observed"],
                        r["predicted"],
                        bool(r["within"]),
                    )
        table = pd.DataFrame(
            [
                {
                    "model": key[0],
                    "value": v,
                    "observed": o,
                    "predicted": p,
                    "within": w,
                }
                for key, (v, o, p, w) in pairs.items()
            ],
            columns=["model", "value", "observed", "predicted", "within"],
        )
        rows: list[dict[str, object]] = []
        for value, group in table.groupby("value", sort=True):
            for model, part in [("all", group), *group.groupby("model", sort=True)]:
                fold = np.abs(np.log10(part["predicted"] / part["observed"]))
                rows.append(
                    {
                        "field": field,
                        "value": value,
                        "model": model,
                        "n_pairs": len(part),
                        "n_within": int(part["within"].sum()),
                        "within": 100.0 * float(part["within"].mean()),
                        "gmfe": float(10 ** fold.mean()),
                    }
                )
        return pd.DataFrame(rows, columns=list(QUALIFICATION_COLUMNS))


#: the evidence extension of the ratio pairs, registered by `pkpdsim.model.PKPDModel`
RATIOS: Final = RatiosExtension()
