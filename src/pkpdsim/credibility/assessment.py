"""The assessment of the contexts of use of a model as plain data.

`assess` evaluates the adopted fits of a model (`evaluate_fits`) and assesses
its declaration on their data points, `assess_table` does the same on a given
table without a model. Per context of use:

- The metrics of every evidence: the points, mappings and studies,
  `criteria_metrics` and the values of the acceptance rules of the context.
- The acceptance: the rules on the union of the points of the decisive
  evidence (`ContextOfUse.decisive_role`, calibration for C1, else
  validation), each row of the table once.
- The coverage of every claimed value of the domain: `validated` with
  applicable points of validation evidence, `calibrated` with applicable points
  of calibration evidence only, `uncovered` without; a C1 context counts
  `calibrated` as covered. A point is applicable when the fold criteria can
  evaluate it, a positive and finite value and prediction (`applicable`): the
  points of a failed simulation or of zeros cover nothing. The values of the
  points outside the claimed domain are reported as unclaimed values, a point
  without the field has none.
- The coverage of every quantity of interest which is no ratio: `covered` with
  an applicable decisive point of its observable inside the claimed domain
  (every claimed field of the point has a claimed value), `uncovered` without.
- The minimum of evaluable points, `ContextOfUse.minimum_points`.
- The verification of the decisive predictions, when the assessment is
  verified (`pkpdsim.credibility.verification`): per check whether it passed,
  the largest relative difference of the alternative predictions and whether
  the acceptance changes with them, and the checks of the battery of the
  builds of the model (`pkpdsim.credibility.battery`), the same for every
  context.
- The verdict, see `verdict`, and the factors of V&V 40 (`FACTORS`).

The ratio quantities and the statement of the domain are recorded and shown,
not assessed: the ratios in step 3, the statement not at all.

The limitations are reported apart, with their metrics, and enter no
acceptance. The assessment writes and reads JSON (`to_json`, `from_json`),
with `NaN` as `null`.
"""

from __future__ import annotations

import dataclasses
import importlib.metadata
import json
import logging
import math
import subprocess
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime
from enum import Enum, StrEnum, unique
from itertools import chain
from typing import TYPE_CHECKING, Any, Final, Self, TypedDict

import numpy as np
import pandas as pd
from sbmlsim.fit.evaluation import Fold, criteria_metrics

from pkpdsim.credibility.declaration import (
    AcceptanceResult,
    AcceptanceRule,
    ContextOfUse,
    CouLevel,
    Credibility,
    EvidenceRole,
    Level,
    Selection,
    Unclaimed,
)
from pkpdsim.credibility.factors import FACTORS, FactorRecord, FactorStatus
from pkpdsim.credibility.ratios import RatioRecord, ratio_pairs
from pkpdsim.credibility.records import MappingLookup, evaluate_fits
from pkpdsim.credibility.selection import (
    EvidenceSelection,
    SelectionError,
    normalize_observable,
    select,
)
from pkpdsim.credibility.verification import (
    NUMERICAL_CHECKS,
    TIGHTENING,
    Check,
    CheckRecord,
    Verification,
    VerificationSummary,
    check_record_from_dict,
    check_record_to_dict,
    relative_differences,
    verify,
)
from pkpdsim.fit import MetadataFilter, NamedFilter

if TYPE_CHECKING:
    from pathlib import Path

    from pkpdsim.credibility.records import FitTable
    from pkpdsim.model import PKPDModel

__all__ = [
    "VERSIONED_PACKAGES",
    "AcceptanceRecord",
    "ContextResult",
    "Coverage",
    "CoverageEntry",
    "CredibilityAssessment",
    "EvidenceResult",
    "Provenance",
    "QuantityRecord",
    "UnclaimedRecord",
    "Verdict",
    "applicable",
    "assess",
    "assess_table",
    "describe_selection",
    "provenance",
    "verdict",
]

logger = logging.getLogger(__name__)

#: the packages whose versions the provenance records
VERSIONED_PACKAGES: Final = ("pkpdsim", "sbmlsim", "sbmlutils", "libroadrunner")

_VERDICT_LABELS: Final = {
    "supported": "supported",
    "supported_with_extrapolation": "supported with extrapolation",
    "not_supported": "not supported",
    "insufficient_evidence": "insufficient evidence",
    "not_verified": "not verified",
}


@unique
class Verdict(StrEnum):
    """The verdict of the assessment of a context of use."""

    SUPPORTED = "supported"
    SUPPORTED_WITH_EXTRAPOLATION = "supported_with_extrapolation"
    NOT_SUPPORTED = "not_supported"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"
    NOT_VERIFIED = "not_verified"

    @property
    def label(self) -> str:
        """Get the verdict in words, e.g. `supported with extrapolation`."""
        return _VERDICT_LABELS[self.value]


@unique
class Coverage(StrEnum):
    """The coverage of a claimed value of the domain of a context of use."""

    VALIDATED = "validated"
    CALIBRATED = "calibrated"
    UNCOVERED = "uncovered"


#: a quantity of interest with an applicable decisive point in the claimed domain
COVERED: Final = "covered"
#: a quantity of interest without one
UNCOVERED: Final = "uncovered"


class QuantityRecord(TypedDict):
    """A quantity of interest as plain data, see `QuantityOfInterest`.

    Attributes:
        id: the id of the quantity.
        description: what the quantity is.
        ratio: whether the quantity is a ratio between arms.
        coverage: `covered` with an applicable decisive point of the quantity
            inside the claimed domain, an assessed ratio with an evaluable
            pair whose numerator lies inside it, else `uncovered`; `None` for
            a ratio without a pairing and in a record of pkpdsim 0.17, which
            did not assess it.
        n_points: the applicable decisive points of the quantity inside the
            claimed domain, the evaluable pairs of an assessed ratio, 0 for a
            ratio without a pairing.
    """

    id: str
    description: str
    ratio: bool
    coverage: str | None
    n_points: int


class AcceptanceRecord(TypedDict):
    """The result of an acceptance rule as plain data, see `AcceptanceResult`."""

    name: str
    value: float
    limit: str
    met: bool | None


class UnclaimedRecord(TypedDict):
    """The unclaimed validation points of a fit mapping.

    Attributes:
        fit: the subset of the fit.
        experiment: the experiment.
        mapping: the key of the mapping.
        n: the points.
    """

    fit: str
    experiment: str
    mapping: str
    n: int


def _write(value: float) -> float | None:
    """Get a number for JSON, `None` for `NaN`."""
    return None if math.isnan(value) else value


def _read(value: float | None) -> float:
    """Get a number from JSON, `NaN` for `None`."""
    return math.nan if value is None else float(value)


@dataclass(frozen=True, kw_only=True)
class EvidenceResult:
    """The metrics of an evidence of a context of use or of a limitation.

    Attributes:
        context: the id of the context of use, `""` for a limitation of the
            model.
        name: the name of the evidence.
        role: the role of the evidence.
        description: what the evidence is.
        selection: what the evidence selects, see `describe_selection`.
        n_points: the data points.
        n_mappings: the fit mappings, by experiment and key.
        n_studies: the studies, the experiments of the mappings.
        outliers_dropped: the outlier points a validation dropped.
        limitation_subtracted: the limitation points subtracted.
        metrics: `criteria_metrics` and the values of the acceptance rules of
            the context, `NaN` if not evaluable.
    """

    context: str
    name: str
    role: EvidenceRole
    description: str
    selection: str
    n_points: int
    n_mappings: int
    n_studies: int
    outliers_dropped: int
    limitation_subtracted: int
    metrics: dict[str, float]

    def to_dict(self) -> dict[str, object]:
        """Get the result as plain data, `NaN` as `None`."""
        return {
            "context": self.context,
            "name": self.name,
            "role": self.role.value,
            "description": self.description,
            "selection": self.selection,
            "n_points": self.n_points,
            "n_mappings": self.n_mappings,
            "n_studies": self.n_studies,
            "outliers_dropped": self.outliers_dropped,
            "limitation_subtracted": self.limitation_subtracted,
            "metrics": {key: _write(value) for key, value in self.metrics.items()},
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the result from plain data, `None` as `NaN`."""
        return cls(
            context=data["context"],
            name=data["name"],
            role=EvidenceRole(data["role"]),
            description=data["description"],
            selection=data["selection"],
            n_points=data["n_points"],
            n_mappings=data["n_mappings"],
            n_studies=data["n_studies"],
            outliers_dropped=data["outliers_dropped"],
            limitation_subtracted=data["limitation_subtracted"],
            metrics={key: _read(value) for key, value in data["metrics"].items()},
        )


@dataclass(frozen=True, kw_only=True)
class CoverageEntry:
    """The coverage of a value of a field of the metadata by a context of use.

    Attributes:
        context: the id of the context of use.
        field: the field of the metadata, e.g. `health`.
        value: the name of the value, e.g. `RENAL_IMPAIRMENT`.
        coverage: the coverage of a claimed value, `None` for a value outside
            the claimed domain.
        claimed: whether the domain claims the value; the values of the points
            of the evidence outside the claim are not claimed.
        n_points: the points with the value: the applicable points of the
            validation evidence for a validated value, of the calibration
            evidence for a calibrated one, every point of the evidence for a
            value outside the claim.
    """

    context: str
    field: str
    value: str
    coverage: Coverage | None
    claimed: bool
    n_points: int

    def to_dict(self) -> dict[str, object]:
        """Get the entry as plain data."""
        return {
            "context": self.context,
            "field": self.field,
            "value": self.value,
            "coverage": None if self.coverage is None else self.coverage.value,
            "claimed": self.claimed,
            "n_points": self.n_points,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the entry from plain data."""
        coverage = data["coverage"]
        return cls(
            context=data["context"],
            field=data["field"],
            value=data["value"],
            coverage=None if coverage is None else Coverage(coverage),
            claimed=data["claimed"],
            n_points=data["n_points"],
        )


@dataclass(frozen=True, kw_only=True)
class ContextResult:
    """The assessment of a context of use.

    Attributes:
        id: the id of the context.
        title: the title.
        question: the question.
        level: the level on the ladder.
        influence: the influence of the model on the decision.
        consequence: the consequence of a wrong decision.
        risk: the model risk.
        rationale: why this level, influence and consequence.
        extrapolation: the argument for a domain without validation evidence.
        domain_statement: what the domain claims beyond the metadata.
        quantities: the quantities of interest (`id`, `description`, `ratio`,
            `coverage`, `n_points`, see `QuantityRecord`).
        acceptance: the results of the acceptance rules (`name`, `value`,
            `limit`, `met`, see `AcceptanceResult`).
        verdict: the verdict.
        factors: the factors of V&V 40 (`id`, `name`, `status`, `statement`,
            `step`), in the order of `FACTORS`.
        n_decisive: the points of the decisive evidence, each row once.
        n_evaluable: the decisive points which are applicable, see
            `applicable`.
        minimum_points: the evaluable points the verdict needs, see
            `ContextOfUse.minimum_points`; 0 in a record of pkpdsim 0.17.
        verification: the checks of the verification of the decisive
            predictions, in the order of `Check`; empty if the assessment was
            not verified.
        ratios: the pairs of the arms of the assessed ratio quantities, see
            `pkpdsim.credibility.ratios`; empty without one and in a record of
            pkpdsim 0.18 or earlier.
        ratio_acceptance: the result of the rule of the ratio pairs
            (`ContextOfUse.ratio_rule`), `None` without an assessed ratio.
    """

    id: str
    title: str
    question: str
    level: CouLevel
    influence: Level
    consequence: Level
    risk: Level
    rationale: str
    extrapolation: str
    domain_statement: str
    quantities: list[QuantityRecord]
    acceptance: list[AcceptanceRecord]
    verdict: Verdict
    factors: list[FactorRecord]
    n_decisive: int
    n_evaluable: int
    minimum_points: int
    verification: list[CheckRecord] = dataclasses.field(default_factory=list)
    ratios: list[RatioRecord] = dataclasses.field(default_factory=list)
    ratio_acceptance: AcceptanceRecord | None = None

    @property
    def verified(self) -> bool | None:
        """Get whether every check passed, `None` if not verified."""
        if not self.verification:
            return None
        return all(record["passed"] for record in self.verification)

    @property
    def decisive_role(self) -> EvidenceRole:
        """Get the role of the evidence which decides, see `CouLevel.decisive_role`."""
        return self.level.decisive_role

    @property
    def uncovered_quantities(self) -> list[str]:
        """Get the ids of the quantities of interest which are uncovered."""
        return [q["id"] for q in self.quantities if q["coverage"] == UNCOVERED]

    def to_dict(self) -> dict[str, object]:
        """Get the result as plain data, `NaN` as `None`."""
        return {
            "id": self.id,
            "title": self.title,
            "question": self.question,
            "level": self.level.value,
            "influence": self.influence.value,
            "consequence": self.consequence.value,
            "risk": self.risk.value,
            "rationale": self.rationale,
            "extrapolation": self.extrapolation,
            "domain_statement": self.domain_statement,
            "quantities": [dict(q) for q in self.quantities],
            "acceptance": [{**a, "value": _write(a["value"])} for a in self.acceptance],
            "verdict": self.verdict.value,
            "factors": [dict(f) for f in self.factors],
            "n_decisive": self.n_decisive,
            "n_evaluable": self.n_evaluable,
            "minimum_points": self.minimum_points,
            "verification": [check_record_to_dict(r) for r in self.verification],
            "ratios": [
                {
                    **r,
                    "observed": _write(r["observed"]),
                    "predicted": _write(r["predicted"]),
                    "limit": _write(r["limit"]),
                }
                for r in self.ratios
            ],
            "ratio_acceptance": (
                None
                if self.ratio_acceptance is None
                else {
                    **self.ratio_acceptance,
                    "value": _write(self.ratio_acceptance["value"]),
                }
            ),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the result from plain data, `None` as `NaN`.

        A record of pkpdsim 0.17 has no minimum of points, read as 0, and no
        coverage of the quantities, read as `None` with 0 points; a record of
        0.18 or earlier has no ratio pairs.
        """
        return cls(
            id=data["id"],
            title=data["title"],
            question=data["question"],
            level=CouLevel(data["level"]),
            influence=Level(data["influence"]),
            consequence=Level(data["consequence"]),
            risk=Level(data["risk"]),
            rationale=data["rationale"],
            extrapolation=data["extrapolation"],
            domain_statement=data["domain_statement"],
            quantities=[
                QuantityRecord(
                    id=q["id"],
                    description=q["description"],
                    ratio=q["ratio"],
                    coverage=q.get("coverage"),
                    n_points=q.get("n_points", 0),
                )
                for q in data["quantities"]
            ],
            acceptance=[
                AcceptanceRecord(
                    name=a["name"],
                    value=_read(a["value"]),
                    limit=a["limit"],
                    met=a["met"],
                )
                for a in data["acceptance"]
            ],
            verdict=Verdict(data["verdict"]),
            factors=[
                FactorRecord(
                    id=f["id"],
                    name=f["name"],
                    status=f["status"],
                    statement=f["statement"],
                    step=f["step"],
                )
                for f in data["factors"]
            ],
            n_decisive=data["n_decisive"],
            n_evaluable=data["n_evaluable"],
            minimum_points=data.get("minimum_points", 0),
            # a record of pkpdsim 0.17 was not verified
            verification=[
                check_record_from_dict(r) for r in data.get("verification", ())
            ],
            # a record of pkpdsim 0.18 or earlier assessed no ratio
            ratios=[
                RatioRecord(
                    quantity=r["quantity"],
                    fit=r["fit"],
                    experiment=r["experiment"],
                    numerator=r["numerator"],
                    denominator=r["denominator"],
                    observed=_read(r["observed"]),
                    predicted=_read(r["predicted"]),
                    limit=_read(r["limit"]),
                    within=r["within"],
                    metadata=dict(r.get("metadata", {})),
                )
                for r in data.get("ratios", ())
            ],
            ratio_acceptance=(
                None
                if data.get("ratio_acceptance") is None
                else AcceptanceRecord(
                    name=data["ratio_acceptance"]["name"],
                    value=_read(data["ratio_acceptance"]["value"]),
                    limit=data["ratio_acceptance"]["limit"],
                    met=data["ratio_acceptance"]["met"],
                )
            ),
        )


@dataclass(frozen=True, kw_only=True)
class Provenance:
    """What an assessment was computed from.

    Attributes:
        model: the name of the model.
        title: the title of the model.
        fits: the ids of the adopted fits.
        versions: the versions of the packages, see `VERSIONED_PACKAGES`.
        commit: the commit of the repository of the model, `None` outside one.
        dirty: whether the repository has changes, `None` outside one.
        date: the date of the assessment, ISO 8601 with the time zone.
        lock: whether the evidence matched the lock of the declaration
            (`matched`), the lock was written with the assessment (`written`),
            or `None` without a lock, see `pkpdsim.credibility.lock`.
    """

    model: str
    title: str
    fits: list[str]
    versions: dict[str, str]
    commit: str | None
    dirty: bool | None
    date: str
    lock: str | None = None

    def to_dict(self) -> dict[str, object]:
        """Get the provenance as plain data."""
        return {
            "model": self.model,
            "title": self.title,
            "fits": list(self.fits),
            "versions": dict(self.versions),
            "commit": self.commit,
            "dirty": self.dirty,
            "date": self.date,
            "lock": self.lock,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the provenance from plain data."""
        return cls(
            model=data["model"],
            title=data["title"],
            fits=list(data["fits"]),
            versions=dict(data["versions"]),
            commit=data["commit"],
            dirty=data["dirty"],
            date=data["date"],
            # a record of pkpdsim 0.18 or earlier knew no lock
            lock=data.get("lock"),
        )


@dataclass(frozen=True, kw_only=True)
class CredibilityAssessment:
    """The assessment of the contexts of use of a model.

    Attributes:
        provenance: what the assessment was computed from.
        contexts: the result of every context of use, in the declared order.
        evidence: the metrics of every evidence of every context.
        coverage: the coverage of the domain of every context.
        limitations: the metrics of the limitations, of the model (context
            `""`) and of the contexts; they enter no acceptance.
        unclaimed: the validation points which no evidence and no limitation
            selects, per fit and mapping (`fit`, `experiment`, `mapping`, `n`).
        notes: the notes of the declaration.
        verification: what the verification compared, `None` if the
            assessment was not verified.
    """

    provenance: Provenance
    contexts: list[ContextResult]
    evidence: list[EvidenceResult]
    coverage: list[CoverageEntry]
    limitations: list[EvidenceResult]
    unclaimed: list[UnclaimedRecord]
    notes: str
    verification: VerificationSummary | None = None

    def to_dict(self) -> dict[str, object]:
        """Get the assessment as plain data, `NaN` as `None`."""
        return {
            "provenance": self.provenance.to_dict(),
            "contexts": [c.to_dict() for c in self.contexts],
            "evidence": [e.to_dict() for e in self.evidence],
            "coverage": [c.to_dict() for c in self.coverage],
            "limitations": [e.to_dict() for e in self.limitations],
            "unclaimed": [dict(u) for u in self.unclaimed],
            "notes": self.notes,
            "verification": (
                None if self.verification is None else self.verification.to_dict()
            ),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the assessment from plain data, `None` as `NaN`."""
        return cls(
            provenance=Provenance.from_dict(data["provenance"]),
            contexts=[ContextResult.from_dict(c) for c in data["contexts"]],
            evidence=[EvidenceResult.from_dict(e) for e in data["evidence"]],
            coverage=[CoverageEntry.from_dict(c) for c in data["coverage"]],
            limitations=[EvidenceResult.from_dict(e) for e in data["limitations"]],
            unclaimed=[
                UnclaimedRecord(
                    fit=u["fit"],
                    experiment=u["experiment"],
                    mapping=u["mapping"],
                    n=u["n"],
                )
                for u in data["unclaimed"]
            ],
            notes=data["notes"],
            verification=(
                None
                if data.get("verification") is None
                else VerificationSummary.from_dict(data["verification"])
            ),
        )

    def to_json(self, path: Path) -> None:
        """Write the assessment as JSON.

        Args:
            path: the file.
        """
        text = json.dumps(self.to_dict(), indent=2, allow_nan=False)
        path.write_text(text + "\n", encoding="utf-8")

    @classmethod
    def from_json(cls, path: Path) -> CredibilityAssessment:
        """Read an assessment written with `to_json`.

        Args:
            path: the file.

        Returns:
            The assessment.
        """
        return cls.from_dict(json.loads(path.read_text(encoding="utf-8")))


def _names(values: Iterable[object]) -> str:
    """Show values sorted, the members of an enum by their name."""
    return ", ".join(sorted(v.name if isinstance(v, Enum) else str(v) for v in values))


def describe_selection(selection: Selection) -> str:
    """Describe what a selection of an evidence selects.

    Args:
        selection: the selection.

    Returns:
        The mapping keys by study for a `NamedFilter`, e.g.
        `CredA: fm_dose; CredB: fm_dose, fm_zero`; the clauses of a
        `MetadataFilter`, e.g. `route in [IV, PO] and health not in [HEALTHY]`;
        `unclaimed validation data` for `Unclaimed`, followed by `within` and
        the description of its filter; the qualified name otherwise.
    """
    if isinstance(selection, NamedFilter):
        return "; ".join(
            f"{study}: {_names(keys)}" for study, keys in sorted(selection.mappings)
        )
    if isinstance(selection, MetadataFilter):
        clauses = [f"{name} in [{_names(v)}]" for name, v in selection.allowed]
        clauses += [f"{name} not in [{_names(v)}]" for name, v in selection.excluded]
        return " and ".join(clauses)
    if isinstance(selection, Unclaimed):
        text = "unclaimed validation data"
        if selection.within is not None:
            text += f" within {describe_selection(selection.within)}"
        return text
    return getattr(selection, "__qualname__", type(selection).__qualname__)


def _points(
    table: pd.DataFrame, selections: Iterable[EvidenceSelection]
) -> pd.DataFrame:
    """Get the rows of the union of selections, each row once, in the order of the table."""
    rows = list(chain.from_iterable(s.rows for s in selections))
    return table.loc[table.index.isin(rows)]


def applicable(points: pd.DataFrame) -> pd.DataFrame:
    """Get the points the fold criteria can evaluate.

    The predicate of the ratio criteria of sbmlsim (`Fold`, GMFE, AFE): a
    positive and finite value (`DV`) and prediction (`IPRED`).

    Args:
        points: the data points.

    Returns:
        The applicable points.
    """
    accepted = Fold(2.0).accepted(points).to_numpy()
    return points.loc[np.isfinite(accepted)]


def _metrics(points: pd.DataFrame, rules: Sequence[AcceptanceRule]) -> dict[str, float]:
    """Get `criteria_metrics` and the values of the rules on points."""
    metrics = criteria_metrics(points)
    for rule in rules:
        metrics[rule.name] = rule.evaluate(points).value
    return metrics


def _counts(points: pd.DataFrame) -> tuple[int, int, int]:
    """Get the points, mappings (by experiment and key) and studies."""
    mappings = points[["experiment", "mapping"]].drop_duplicates()
    return len(points), len(mappings), int(points["experiment"].nunique())


def _evidence_result(
    selection: EvidenceSelection, table: pd.DataFrame, rules: Sequence[AcceptanceRule]
) -> EvidenceResult:
    """Get the metrics of the selection of an evidence."""
    points = _points(table, (selection,))
    n_points, n_mappings, n_studies = _counts(points)
    evidence = selection.evidence
    return EvidenceResult(
        context=selection.context,
        name=evidence.name,
        role=evidence.role,
        description=evidence.description,
        selection=describe_selection(evidence.selection),
        n_points=n_points,
        n_mappings=n_mappings,
        n_studies=n_studies,
        outliers_dropped=selection.outliers_dropped,
        limitation_subtracted=selection.limitation_subtracted,
        metrics=_metrics(points, rules),
    )


def _values(points: pd.DataFrame, field: str) -> pd.Series:
    """Get the values of a field of the metadata of points as text, empty if missing."""
    if field not in points.columns:
        return pd.Series("", index=points.index, dtype=str)
    return points[field].fillna("").astype(str)


def _coverage(
    cou: ContextOfUse, table: pd.DataFrame, selections: Sequence[EvidenceSelection]
) -> list[CoverageEntry]:
    """Get the coverage of the claimed values and the values outside the claim."""
    validation = applicable(
        _points(
            table,
            (s for s in selections if s.evidence.role is EvidenceRole.VALIDATION),
        )
    )
    calibration = applicable(
        _points(
            table,
            (s for s in selections if s.evidence.role is EvidenceRole.CALIBRATION),
        )
    )
    evidence = _points(table, selections)
    entries: list[CoverageEntry] = []
    for field, claimed in cou.domain.metadata.items():
        names = [value.name for value in claimed]
        validated = _values(validation, field).value_counts()
        calibrated = _values(calibration, field).value_counts()
        for name in names:
            if validated.get(name, 0):
                coverage, n = Coverage.VALIDATED, int(validated[name])
            elif calibrated.get(name, 0):
                coverage, n = Coverage.CALIBRATED, int(calibrated[name])
            else:
                coverage, n = Coverage.UNCOVERED, 0
            entries.append(
                CoverageEntry(
                    context=cou.id,
                    field=field,
                    value=name,
                    coverage=coverage,
                    claimed=True,
                    n_points=n,
                )
            )
        outside = _values(evidence, field).value_counts()
        entries.extend(
            CoverageEntry(
                context=cou.id,
                field=field,
                value=str(value),
                coverage=None,
                claimed=False,
                n_points=int(outside[value]),
            )
            for value in sorted(outside.index)
            # a point without the field reports no value
            if value and value not in names
        )
    return entries


def _in_domain(cou: ContextOfUse, points: pd.DataFrame) -> pd.DataFrame:
    """Get the points inside the claimed domain: a claimed value of every claimed field.

    A field without claimed values claims nothing and restricts no point.
    """
    inside = pd.Series(True, index=points.index)
    for field, claimed in cou.domain.metadata.items():
        if claimed:
            names = [value.name for value in claimed]
            inside &= _values(points, field).isin(names)
    return points.loc[inside]


def _quantities(
    cou: ContextOfUse, decisive: pd.DataFrame, pairs: Sequence[RatioRecord] = ()
) -> list[QuantityRecord]:
    """Get the quantities of interest with their coverage by the decisive points.

    An assessed ratio is covered by an evaluable pair whose numerator has a
    decisive point inside the claimed domain, its points are these pairs; a
    declared ratio without a pairing has no coverage.
    """
    inside = _in_domain(cou, applicable(decisive))
    observables = inside["observable"].astype(str).map(normalize_observable)
    claimed = set(
        zip(
            inside["fit"].astype(str),
            inside["experiment"].astype(str),
            inside["mapping"].astype(str),
            strict=True,
        )
    )
    records: list[QuantityRecord] = []
    for q in cou.quantities:
        if q.pairing is not None:
            n = sum(
                p["quantity"] == q.id
                and p["within"] is not None
                and (p["fit"], p["experiment"], p["numerator"]) in claimed
                for p in pairs
            )
            coverage = COVERED if n else UNCOVERED
        elif q.ratio:
            coverage, n = None, 0
        else:
            n = int((observables == normalize_observable(q.id)).sum())
            coverage = COVERED if n else UNCOVERED
        records.append(
            QuantityRecord(
                id=q.id,
                description=q.description,
                ratio=bool(q.ratio),
                coverage=coverage,
                n_points=n,
            )
        )
    return records


def _ratio_pairs(
    cou: ContextOfUse, decisive: pd.DataFrame, lookup: MappingLookup
) -> list[RatioRecord]:
    """Get the pairs of the assessed ratio quantities of a context."""
    pairs: list[RatioRecord] = []
    for q in cou.ratios:
        if q.pairing is None:
            continue
        found, problems = ratio_pairs(
            q.id, q.pairing, decisive, lookup, cou.ratio_rule.delta
        )
        if problems:
            # `select` refuses an ambiguous pairing first
            raise SelectionError([f"context of use '{cou.id}': {p}" for p in problems])
        pairs.extend(found)
    return pairs


def _covered(cou: ContextOfUse, entry: CoverageEntry) -> bool:
    """Check whether a claimed value is covered: validated, for C1 also calibrated."""
    if entry.coverage is Coverage.VALIDATED:
        return True
    return cou.level is CouLevel.C1 and entry.coverage is Coverage.CALIBRATED


def verdict(
    n_decisive: int,
    acceptance: Sequence[AcceptanceResult],
    covered: bool,
    extrapolation: str,
    *,
    level: CouLevel,
    n_evaluable: int,
    minimum: int,
    verified: bool | None = None,
) -> Verdict:
    """Get the verdict of a context of use.

    In this order: no decisive point is insufficient evidence; a failed check
    of the verification is not verified, since the comparison of the
    predictions with the data is not to be trusted; a rule which is not met is
    not supported, also on few points; a rule which is not evaluable
    is insufficient evidence; fewer evaluable points than the minimum is
    insufficient evidence; every claimed value and every quantity which is no
    ratio covered is supported, for C4 supported with extrapolation (a C4
    context extrapolates without data of its condition by definition, its
    extrapolation argument is required); an extrapolation argument is
    supported with extrapolation; else insufficient evidence.

    Args:
        n_decisive: the points of the decisive evidence.
        acceptance: the results of the acceptance rules.
        covered: whether every claimed value and every quantity which is no
            ratio is covered.
        extrapolation: the extrapolation argument of the context.
        level: the level of the context.
        n_evaluable: the decisive points which are applicable.
        minimum: the evaluable points the verdict needs.
        verified: whether every check of the verification passed, `None` if
            the assessment is not verified.

    Returns:
        The verdict.
    """
    if n_decisive == 0:
        return Verdict.INSUFFICIENT_EVIDENCE
    if verified is False:
        return Verdict.NOT_VERIFIED
    if any(result.met is False for result in acceptance):
        return Verdict.NOT_SUPPORTED
    if any(result.met is None for result in acceptance):
        return Verdict.INSUFFICIENT_EVIDENCE
    if n_evaluable < minimum:
        return Verdict.INSUFFICIENT_EVIDENCE
    if covered and level is not CouLevel.C4:
        return Verdict.SUPPORTED
    if extrapolation:
        return Verdict.SUPPORTED_WITH_EXTRAPOLATION
    return Verdict.INSUFFICIENT_EVIDENCE


def _count(n: int, singular: str, plural: str | None = None) -> str:
    """Get a number of things in words, e.g. `1 point` and `2 points`."""
    return f"{n} {singular if n == 1 else plural or singular + 's'}"


def _check_statement(record: CheckRecord, tightening: float) -> str:
    """Get the statement of the factor of a check of the verification."""
    result = "passed" if record["passed"] else "failed"
    if record["check"] not in {c.value for c in (*NUMERICAL_CHECKS, Check.SETUP)}:
        label = Check(record["check"]).label
        text = f"{label[0].upper()}{label[1:]} of {_count(record['n_points'], 'case')}"
        if not math.isnan(record["max_difference"]):
            text += f", largest deviation {record['max_difference']:.2g}"
            if record["at"]:
                text += f" at {record['at']}"
        problems = "".join(f"; {problem}" for problem in record["problems"])
        return f"{text}{problems}: {result}."
    if record["check"] == Check.SETUP.value:
        text = (
            "The builds are valid SBML with consistent units and the "
            f"{_count(record['n_points'], 'decisive prediction')} are finite"
            if not record["problems"]
            else "; ".join(record["problems"])
        )
        return f"{text}: {result}."
    compared = (
        f"the tolerances of the fits divided by {tightening:g}"
        if record["check"] == Check.TOLERANCES.value
        else "the BDF of CVODE of at most order 2 at the tightened tolerances"
    )
    if record["n_points"] == 0:
        difference = "no decisive prediction to compare"
    else:
        difference = (
            f"largest relative difference {record['max_difference']:.2g} of "
            f"{_count(record['n_points'], 'decisive prediction')} with {compared}, "
            f"at {record['at']}, acceptance "
            f"{'changed' if record['acceptance_changed'] else 'unchanged'}"
        )
    problems = "".join(f"; {problem}" for problem in record["problems"])
    return f"{difference[0].upper()}{difference[1:]}{problems}: {result}."


def _factors(
    cou: ContextOfUse,
    decisive: pd.DataFrame,
    acceptance: Sequence[AcceptanceResult],
    coverage: Sequence[CoverageEntry],
    quantities: Sequence[QuantityRecord],
    checks: Sequence[CheckRecord] = (),
    tightening: float = TIGHTENING,
    pairs: Sequence[RatioRecord] = (),
) -> list[FactorRecord]:
    """Get the factors of V&V 40 of a context, the assessed ones with its results.

    The factors of the checks of the verification are assessed when the
    assessment was verified, else they keep their status of `FACTORS`; the
    relevance of the quantities of interest is assessed when the context has
    an assessed ratio, by its pairs.
    """
    role = cou.decisive_role.value
    n_points, n_mappings, n_studies = _counts(decisive)
    claimed = [e for e in coverage if e.claimed]
    outside = [e for e in coverage if not e.claimed]
    met = sum(result.met is True for result in acceptance)
    evaluable = sum(result.met is not None for result in acceptance)
    covered = sum(_covered(cou, e) for e in claimed)
    absolute = [q for q in quantities if not q["ratio"]]
    covered_quantities = sum(q["coverage"] == COVERED for q in absolute)
    statements = {
        "test_samples": (
            f"{_count(n_points, 'point')} of {_count(n_mappings, 'mapping')} of "
            f"{_count(n_studies, 'study', 'studies')} of the {role} evidence "
            f"decide, {len(applicable(decisive))} of them evaluable "
            f"(minimum {cou.minimum_points})."
        ),
        "test_conditions": (
            f"{_count(len(claimed), 'claimed value')} of the metadata, "
            f"{_count(len(outside), 'value')} of the evidence outside the claim."
        ),
        "output_comparison": (
            f"{met} of {_count(len(acceptance), 'acceptance rule')} met on the "
            f"{role} evidence, {len(acceptance) - evaluable} not evaluable."
        ),
        "relevance_of_validation": (
            f"{covered} of {_count(len(claimed), 'claimed value')} and "
            f"{covered_quantities} of {_count(len(absolute), 'quantity', 'quantities')} "
            "of interest covered."
        ),
    }
    verified: dict[str, list[CheckRecord]] = {}
    for record in checks:
        verified.setdefault(record["factor"], []).append(record)
    factors: list[FactorRecord] = []
    for factor in FACTORS:
        data = factor.to_dict()
        if factor.status is FactorStatus.ASSESSED:
            data["statement"] = statements[factor.id]
        elif factor.id in verified:
            # a platform statement stays, the checks of the model follow it
            platform = (
                factor.statement if factor.status is FactorStatus.PLATFORM else ""
            )
            data["status"] = FactorStatus.ASSESSED.value
            data["statement"] = " ".join(
                [
                    *([platform] if platform else []),
                    *(_check_statement(r, tightening) for r in verified[factor.id]),
                ]
            )
            data["step"] = None
        elif factor.id == "relevance_of_quantities" and cou.ratios:
            evaluable = [p for p in pairs if p["within"] is not None]
            within = sum(p["within"] is True for p in evaluable)
            ratios = [q for q in quantities if q["id"] in {r.id for r in cou.ratios}]
            covered_ratios = sum(q["coverage"] == COVERED for q in ratios)
            data["status"] = FactorStatus.ASSESSED.value
            data["statement"] = (
                f"{within} of {_count(len(evaluable), 'ratio pair')} within the "
                f"limits of Guest et al. ({cou.ratio_rule.limit} %), "
                f"{covered_ratios} of "
                f"{_count(len(ratios), 'ratio quantity', 'ratio quantities')} "
                "covered."
            )
            data["step"] = None
        factors.append(data)
    return factors


def _at(point: pd.Series) -> str:
    """Name a data point, `experiment.mapping at x` (without `x` if it has none)."""
    name = f"{point['experiment']}.{point['mapping']}"
    x = point.get("x")
    return name if x is None or pd.isna(x) else f"{name} at {float(x):g}"


def _check_records(
    cou: ContextOfUse,
    table: pd.DataFrame,
    decisive: pd.DataFrame,
    acceptance: Sequence[AcceptanceResult],
    verification: Verification,
    lookup: MappingLookup | None = None,
) -> list[CheckRecord]:
    """Get the checks of the verification of the decisive predictions of a context.

    A numerical check passes when its evaluation of every fit of the decisive
    points succeeded, the largest relative difference of a decisive prediction
    is at most the limit of the verification and no rule of the acceptance
    changes; the setup passes when the builds have no issues and every decisive
    prediction is finite.
    """
    summary = verification.summary
    fits = set(decisive["fit"].astype(str)) if "fit" in decisive.columns else set()
    records: list[CheckRecord] = []
    for check in NUMERICAL_CHECKS:
        predictions = verification.predictions[check]
        differences = relative_differences(table, predictions, summary.scale_floor)
        compared = differences.loc[decisive.index].dropna()
        alternative = decisive.assign(IPRED=predictions.loc[decisive.index])
        alternatives = [rule.evaluate(alternative) for rule in cou.rules]
        if cou.ratios and lookup is not None:
            alternatives.append(
                cou.ratio_rule.evaluate(
                    [p["within"] for p in _ratio_pairs(cou, alternative, lookup)]
                )
            )
        changed = any(
            other.met != result.met
            for other, result in zip(alternatives, acceptance, strict=True)
        )
        problems = [
            message
            for fit, message in verification.failures[check].items()
            if fit in fits
        ]
        largest = float(compared.max()) if len(compared) else math.nan
        records.append(
            CheckRecord(
                check=check.value,
                factor=check.factor,
                n_points=len(compared),
                max_difference=largest,
                at=_at(decisive.loc[compared.idxmax()]) if len(compared) else "",
                acceptance_changed=changed,
                passed=not problems and not changed and not (largest > summary.limit),
                problems=problems,
            )
        )
    ipred = decisive["IPRED"].to_numpy(dtype=float)
    failed = decisive.loc[~np.isfinite(ipred)]
    problems = list(summary.problems)
    if len(failed):
        mappings = failed[["experiment", "mapping"]].drop_duplicates()
        problems.append(
            f"{len(failed)} of {_count(len(decisive), 'decisive prediction')} "
            f"{'is' if len(failed) == 1 else 'are'} not finite, of "
            + ", ".join(f"{e}.{m}" for e, m in mappings.itertuples(index=False))
        )
    records.append(
        CheckRecord(
            check=Check.SETUP.value,
            factor=Check.SETUP.factor,
            n_points=len(decisive),
            max_difference=math.nan,
            at="",
            acceptance_changed=None,
            passed=not problems,
            problems=problems,
        )
    )
    if summary.battery is not None:
        records.extend(
            CheckRecord(
                check=result.check.value,
                factor=result.check.factor,
                n_points=result.cases,
                max_difference=result.largest,
                at=result.at,
                acceptance_changed=None,
                passed=bool(result.passed),
                problems=list(result.problems),
            )
            for result in summary.battery.results
            if result.passed is not None
        )
    return records


def _context_result(
    cou: ContextOfUse,
    table: pd.DataFrame,
    selections: Sequence[EvidenceSelection],
    verification: Verification | None = None,
    lookup: MappingLookup | None = None,
) -> tuple[ContextResult, list[CoverageEntry]]:
    """Get the result and the coverage of a context of use.

    The lookup of the fit mappings pairs the arms of the assessed ratio
    quantities; a context without one needs none.
    """
    decisive = _points(
        table, (s for s in selections if s.evidence.role is cou.decisive_role)
    )
    acceptance = [rule.evaluate(decisive) for rule in cou.rules]
    pairs: list[RatioRecord] = []
    ratio_result: AcceptanceResult | None = None
    if cou.ratios:
        if lookup is None:
            raise ValueError(
                f"The context of use '{cou.id}' has a ratio, which needs the "
                "lookup of the fit mappings."
            )
        pairs = _ratio_pairs(cou, decisive, lookup)
        ratio_result = cou.ratio_rule.evaluate([p["within"] for p in pairs])
    coverage = _coverage(cou, table, selections)
    quantities = _quantities(cou, decisive, pairs)
    covered = all(_covered(cou, e) for e in coverage if e.claimed) and all(
        q["coverage"] == COVERED for q in quantities if q["coverage"] is not None
    )
    n_evaluable = len(applicable(decisive))
    checks = (
        []
        if verification is None
        else _check_records(
            cou,
            table,
            decisive,
            [*acceptance, *([] if ratio_result is None else [ratio_result])],
            verification,
            lookup,
        )
    )
    verified = None if not checks else all(record["passed"] for record in checks)
    result = ContextResult(
        id=cou.id,
        title=cou.title,
        question=cou.question,
        level=cou.level,
        influence=cou.influence,
        consequence=cou.consequence,
        risk=cou.risk,
        rationale=cou.rationale,
        extrapolation=cou.extrapolation,
        domain_statement=cou.domain.statement,
        quantities=quantities,
        acceptance=[
            AcceptanceRecord(name=a.name, value=a.value, limit=a.limit, met=a.met)
            for a in acceptance
        ],
        verdict=verdict(
            len(decisive),
            [*acceptance, *([] if ratio_result is None else [ratio_result])],
            covered,
            cou.extrapolation,
            level=cou.level,
            n_evaluable=n_evaluable,
            minimum=cou.minimum_points,
            verified=verified,
        ),
        factors=_factors(
            cou,
            decisive,
            [*acceptance, *([] if ratio_result is None else [ratio_result])],
            coverage,
            quantities,
            checks,
            TIGHTENING if verification is None else verification.summary.tightening,
            pairs,
        ),
        n_decisive=len(decisive),
        n_evaluable=n_evaluable,
        minimum_points=cou.minimum_points,
        verification=checks,
        ratios=pairs,
        ratio_acceptance=(
            None
            if ratio_result is None
            else AcceptanceRecord(
                name=ratio_result.name,
                value=ratio_result.value,
                limit=ratio_result.limit,
                met=ratio_result.met,
            )
        ),
    )
    return result, coverage


def assess_table(
    credibility: Credibility,
    fits: FitTable,
    provenance: Provenance,
    verification: Verification | None = None,
) -> CredibilityAssessment:
    """Assess the contexts of use of a declaration on the data points of the fits.

    Args:
        credibility: the declaration.
        fits: the data points of the adopted fits with `DV`, `IPRED` and
            `DV_SD`, see `evaluate_fits`.
        provenance: what the table was computed from.
        verification: the alternative predictions of the table and the issues
            of the builds, see `pkpdsim.credibility.verification.verify`; the
            assessment is not verified without.

    Returns:
        The assessment.

    Raises:
        SelectionError: if the evidence selects invalid data, see `select`.
    """
    selections = select(credibility, fits)
    table = fits.table
    evidence: list[EvidenceResult] = []
    contexts: list[ContextResult] = []
    coverage: list[CoverageEntry] = []
    for cou in credibility.contexts:
        own = [s for s in selections.evidence if s.context == cou.id]
        evidence.extend(_evidence_result(s, table, cou.rules) for s in own)
        result, entries = _context_result(cou, table, own, verification, fits.lookup)
        contexts.append(result)
        coverage.extend(entries)
    limitations = [
        _evidence_result(
            s, table, credibility.context(s.context).rules if s.context else ()
        )
        for s in selections.limitation_evidence
    ]
    unclaimed = [
        UnclaimedRecord(
            fit=str(row["fit"]),
            experiment=str(row["experiment"]),
            mapping=str(row["mapping"]),
            n=int(row["n"]),
        )
        for row in selections.unclaimed.to_dict("records")
    ]
    return CredibilityAssessment(
        provenance=provenance,
        contexts=contexts,
        evidence=evidence,
        coverage=coverage,
        limitations=limitations,
        unclaimed=unclaimed,
        notes=credibility.notes,
        verification=None if verification is None else verification.summary,
    )


def _git(model: PKPDModel, *args: str) -> str | None:
    """Run git in the directory of a model, `None` outside a repository."""
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=model.path,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None
    return result.stdout if result.returncode == 0 else None


def provenance(
    model: PKPDModel, source_commit: str | None = None, lock: str | None = None
) -> Provenance:
    """Get the provenance of an assessment of a model now.

    Args:
        model: the model.
        source_commit: the commit of the source repository of the model, given
            by a release, which assesses the model in its release repository:
            the commit is recorded and whether its tree was clean is not known.
        lock: the state of the lock of the evidence, see `Provenance.lock`.

    Returns:
        The versions of `VERSIONED_PACKAGES`, the commit and whether the folder
        of the model has changes (`git status --porcelain -- .` in `model.path`;
        `None` outside a repository, the dirty state `None` with a
        `source_commit`), the adopted fits, the date with the local time
        zone and the state of the lock.
    """
    if source_commit is None:
        commit = _git(model, "rev-parse", "HEAD")
        status = (
            None if commit is None else _git(model, "status", "--porcelain", "--", ".")
        )
        commit = None if commit is None else commit.strip()
        dirty = None if status is None else bool(status.strip())
    else:
        commit, dirty = source_commit, None
    return Provenance(
        model=model.name,
        title=model.title,
        fits=list(model.fits),
        versions={
            package: importlib.metadata.version(package)
            for package in VERSIONED_PACKAGES
        },
        commit=commit,
        dirty=dirty,
        date=datetime.now().astimezone().isoformat(timespec="seconds"),
        lock=lock,
    )


def assess(
    model: PKPDModel,
    source_commit: str | None = None,
    *,
    verification: bool = True,
    lock: bool = False,
) -> CredibilityAssessment:
    """Assess the contexts of use of a model on its adopted fits.

    The fits are evaluated fresh on the current model and data, see
    `evaluate_fits`, and, with `verification`, again with the alternatives of
    the checks of `pkpdsim.credibility.verification`.

    Args:
        model: the model with a declaration and adopted fits.
        source_commit: the commit of the source repository of a release, see
            `provenance`.
        verification: verify the decisive predictions; the factors of the
            verification are not assessed without.
        lock: write the lock of the evidence (`pkpdsim.credibility.lock`)
            instead of comparing the evidence with it.

    Returns:
        The assessment.

    Raises:
        ValueError: if the model has no declaration or no adopted fits.
        KindMismatchError: if a kind differs from the recorded one.
        SelectionError: if the evidence selects invalid data.
        LockMismatchError: if the evidence differs from the lock.
    """
    from pkpdsim.credibility.lock import (
        LockMismatchError,
        lock_differences,
        write_lock,
    )

    credibility = model.load_credibility()
    fits = evaluate_fits(model)
    selections = select(credibility, fits)
    state: str | None
    if lock:
        logger.info("Write the lock '%s'", write_lock(model, selections, fits))
        state = "written"
    else:
        differences = lock_differences(model, selections, fits)
        if differences:
            raise LockMismatchError(differences)
        state = None if differences is None else "matched"
    checks = verify(model, fits) if verification else None
    logger.info("Assess the contexts of use of '%s'", model.name)
    return assess_table(
        credibility, fits, provenance(model, source_commit, state), checks
    )
