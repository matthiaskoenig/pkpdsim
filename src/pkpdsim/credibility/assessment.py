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
- The verdict, see `verdict`, and the factors of V&V 40 (`FACTORS`).

The ratio quantities and the statement of the domain are recorded and shown,
not assessed: the ratios in step 3, the statement not at all.

The limitations are reported apart, with their metrics, and enter no
acceptance. The assessment writes and reads JSON (`to_json`, `from_json`),
with `NaN` as `null`.
"""

from __future__ import annotations

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
from pkpdsim.credibility.records import evaluate_fits
from pkpdsim.credibility.selection import EvidenceSelection, select
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
}


@unique
class Verdict(StrEnum):
    """The verdict of the assessment of a context of use."""

    SUPPORTED = "supported"
    SUPPORTED_WITH_EXTRAPOLATION = "supported_with_extrapolation"
    NOT_SUPPORTED = "not_supported"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"

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


class QuantityRecord(TypedDict):
    """A quantity of interest as plain data, see `QuantityOfInterest`."""

    id: str
    description: str
    ratio: bool


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
        quantities: the quantities of interest (`id`, `description`, `ratio`).
        acceptance: the results of the acceptance rules (`name`, `value`,
            `limit`, `met`, see `AcceptanceResult`).
        verdict: the verdict.
        factors: the factors of V&V 40 (`id`, `name`, `status`, `statement`,
            `step`), in the order of `FACTORS`.
        n_decisive: the points of the decisive evidence, each row once.
        n_evaluable: the decisive points which are applicable, see
            `applicable`.
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

    @property
    def decisive_role(self) -> EvidenceRole:
        """Get the role of the evidence which decides, see `CouLevel.decisive_role`."""
        return self.level.decisive_role

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
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the result from plain data, `None` as `NaN`."""
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
                    id=q["id"], description=q["description"], ratio=q["ratio"]
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
    """

    model: str
    title: str
    fits: list[str]
    versions: dict[str, str]
    commit: str | None
    dirty: bool | None
    date: str

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
    """

    provenance: Provenance
    contexts: list[ContextResult]
    evidence: list[EvidenceResult]
    coverage: list[CoverageEntry]
    limitations: list[EvidenceResult]
    unclaimed: list[UnclaimedRecord]
    notes: str

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
) -> Verdict:
    """Get the verdict of a context of use.

    In this order: no decisive point is insufficient evidence; a rule which is
    not met is not supported; a rule which is not evaluable is insufficient
    evidence; every claimed value covered is supported, for C4 supported with
    extrapolation (a C4 context extrapolates without data of its condition by
    definition, its extrapolation argument is required); an extrapolation
    argument is supported with extrapolation; else insufficient evidence.

    Args:
        n_decisive: the points of the decisive evidence.
        acceptance: the results of the acceptance rules.
        covered: whether every claimed value is covered.
        extrapolation: the extrapolation argument of the context.
        level: the level of the context.

    Returns:
        The verdict.
    """
    if n_decisive == 0:
        return Verdict.INSUFFICIENT_EVIDENCE
    if any(result.met is False for result in acceptance):
        return Verdict.NOT_SUPPORTED
    if any(result.met is None for result in acceptance):
        return Verdict.INSUFFICIENT_EVIDENCE
    if covered and level is not CouLevel.C4:
        return Verdict.SUPPORTED
    if extrapolation:
        return Verdict.SUPPORTED_WITH_EXTRAPOLATION
    return Verdict.INSUFFICIENT_EVIDENCE


def _count(n: int, singular: str, plural: str | None = None) -> str:
    """Get a number of things in words, e.g. `1 point` and `2 points`."""
    return f"{n} {singular if n == 1 else plural or singular + 's'}"


def _factors(
    cou: ContextOfUse,
    decisive: pd.DataFrame,
    acceptance: Sequence[AcceptanceResult],
    coverage: Sequence[CoverageEntry],
) -> list[FactorRecord]:
    """Get the factors of V&V 40 of a context, the assessed ones with its results."""
    role = cou.decisive_role.value
    n_points, n_mappings, n_studies = _counts(decisive)
    claimed = [e for e in coverage if e.claimed]
    outside = [e for e in coverage if not e.claimed]
    met = sum(result.met is True for result in acceptance)
    evaluable = sum(result.met is not None for result in acceptance)
    covered = sum(_covered(cou, e) for e in claimed)
    statements = {
        "test_samples": (
            f"{_count(n_points, 'point')} of {_count(n_mappings, 'mapping')} of "
            f"{_count(n_studies, 'study', 'studies')} of the {role} evidence "
            f"decide, {len(applicable(decisive))} of them evaluable."
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
            f"{covered} of {_count(len(claimed), 'claimed value')} covered."
        ),
    }
    factors: list[FactorRecord] = []
    for factor in FACTORS:
        data = factor.to_dict()
        if factor.status is FactorStatus.ASSESSED:
            data["statement"] = statements[factor.id]
        factors.append(data)
    return factors


def _context_result(
    cou: ContextOfUse,
    table: pd.DataFrame,
    selections: Sequence[EvidenceSelection],
) -> tuple[ContextResult, list[CoverageEntry]]:
    """Get the result and the coverage of a context of use."""
    decisive = _points(
        table, (s for s in selections if s.evidence.role is cou.decisive_role)
    )
    acceptance = [rule.evaluate(decisive) for rule in cou.rules]
    coverage = _coverage(cou, table, selections)
    covered = all(_covered(cou, e) for e in coverage if e.claimed)
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
        quantities=[
            QuantityRecord(id=q.id, description=q.description, ratio=q.ratio)
            for q in cou.quantities
        ],
        acceptance=[
            AcceptanceRecord(name=a.name, value=a.value, limit=a.limit, met=a.met)
            for a in acceptance
        ],
        verdict=verdict(
            len(decisive), acceptance, covered, cou.extrapolation, level=cou.level
        ),
        factors=_factors(cou, decisive, acceptance, coverage),
        n_decisive=len(decisive),
        n_evaluable=len(applicable(decisive)),
    )
    return result, coverage


def assess_table(
    credibility: Credibility, fits: FitTable, provenance: Provenance
) -> CredibilityAssessment:
    """Assess the contexts of use of a declaration on the data points of the fits.

    Args:
        credibility: the declaration.
        fits: the data points of the adopted fits with `DV`, `IPRED` and
            `DV_SD`, see `evaluate_fits`.
        provenance: what the table was computed from.

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
        result, entries = _context_result(cou, table, own)
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


def provenance(model: PKPDModel, source_commit: str | None = None) -> Provenance:
    """Get the provenance of an assessment of a model now.

    Args:
        model: the model.
        source_commit: the commit of the source repository of the model, given
            by a release, which assesses the model in its release repository:
            the commit is recorded and whether its tree was clean is not known.

    Returns:
        The versions of `VERSIONED_PACKAGES`, the commit and whether the folder
        of the model has changes (`git status --porcelain -- .` in `model.path`;
        `None` outside a repository, the dirty state `None` with a
        `source_commit`), the adopted fits and the date with the local time
        zone.
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
    )


def assess(model: PKPDModel, source_commit: str | None = None) -> CredibilityAssessment:
    """Assess the contexts of use of a model on its adopted fits.

    The fits are evaluated fresh on the current model and data, see
    `evaluate_fits`.

    Args:
        model: the model with a declaration and adopted fits.
        source_commit: the commit of the source repository of a release, see
            `provenance`.

    Returns:
        The assessment.

    Raises:
        ValueError: if the model has no declaration or no adopted fits.
        KindMismatchError: if a kind differs from the recorded one.
        SelectionError: if the evidence selects invalid data.
    """
    credibility = model.load_credibility()
    fits = evaluate_fits(model)
    logger.info("Assess the contexts of use of '%s'", model.name)
    return assess_table(credibility, fits, provenance(model, source_commit))
