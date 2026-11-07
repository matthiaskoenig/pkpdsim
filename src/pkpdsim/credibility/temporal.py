"""The temporal validation: the studies a model collected after its fit.

A study which the model collected after the day an adopted fit ran cannot have
shaped the fit, neither as training data nor by the choice of the model while
the modeller looked at it. Its decisive points test whether the acceptance of
a context of use holds on data which did not exist for the model when it was
fitted, the closest the existing models get to a prospective validation.

The day of the collection is the record of the study in the data of the model
(`pkpdsim.data.StudyRecord`, `.<study>_study.json`), which
`pkpdsim.data.write_study_tables` writes on the first collection and keeps
afterwards; the day of a fit is the date of its id (`sbmlsim.fit.result.fit_date`).
The release of a study in PK-DB is no such day, since a model can collect a
study from a checkout of pkdb_data before it is released. `chronology` reads
both for a model, and the evidence extension `TEMPORAL`
(`TemporalExtension`, name `temporal`) of the credibility of sbmlsim splits
the decisive points of every context of use above C1 by it (`TemporalRecord`,
the section `temporal` of the record of a context): the acceptance rules of
the context on the new points are recorded and shown, the verdict does not
change, and the factor output comparison states them; the scorecard has the
evaluable new points (`temporal_points`) and whether every rule is met on
them (`temporal_met`), a rule met on them before and no longer is a
regression.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Final, TypedDict, override

import pandas as pd
from sbmlsim.credibility.assessment import AcceptanceRecord, applicable
from sbmlsim.credibility.declaration import EvidenceRole
from sbmlsim.credibility.extensions import (
    ContextBlock,
    ContextContribution,
    EvidenceExtension,
    FactorStatement,
    ScorecardChange,
    scorecard_true,
)
from sbmlsim.fit.result import fit_date, fit_subset

from pkpdsim.credibility.blocks import render_macro
from pkpdsim.credibility.studies import study_folders
from pkpdsim.data import StudyRecord, read_study_record

if TYPE_CHECKING:
    from datetime import date

    from sbmlsim.credibility.declaration import ContextOfUse, Credibility
    from sbmlsim.credibility.records import FitTable
    from sbmlsim.project import Project

__all__ = [
    "TEMPORAL",
    "TEMPORAL_COLUMNS",
    "Chronology",
    "TemporalExtension",
    "TemporalRecord",
    "TemporalStudy",
    "chronology",
    "temporal_record",
    "temporal_statement",
]

#: the columns of the table of the new studies, see `TemporalStudy`
TEMPORAL_COLUMNS = (
    "study",
    "fit",
    "fit_date",
    "collected",
    "released",
    "pkdb_id",
    "n_points",
    "n_evaluable",
)


class TemporalStudy(TypedDict):
    """A study a model collected after the fit of its decisive points.

    Attributes:
        study: the study (the experiment).
        fit: the subset of the fit.
        fit_date: the day the fit ran, ISO format.
        collected: the day the model collected the study, ISO format.
        released: the day PK-DB released the study, `None` if not released.
        pkdb_id: the identifier of the release in PK-DB.
        n_points: the decisive points of the study of this fit.
        n_evaluable: those the fold criteria can evaluate.
    """

    study: str
    fit: str
    fit_date: str
    collected: str
    released: str | None
    pkdb_id: str | None
    n_points: int
    n_evaluable: int


class TemporalRecord(TypedDict):
    """The decisive points of a context of use which are new for their fit.

    Attributes:
        studies: the new studies, by study and fit.
        n_points: the decisive points of the new studies.
        n_evaluable: those the fold criteria can evaluate.
        acceptance: the acceptance rules of the context on the new points.
        n_known: the decisive studies which the model collected no later than
            the day of their fit.
        undated: the decisive studies without a day of collection.
    """

    studies: list[TemporalStudy]
    n_points: int
    n_evaluable: int
    acceptance: list[AcceptanceRecord]
    n_known: int
    undated: list[str]


@dataclass(frozen=True)
class Chronology:
    """When the adopted fits of a model ran and when it collected its studies.

    Attributes:
        fits: the day of every adopted fit, by its subset.
        studies: the record of every study of the data of the model, by study.
    """

    fits: Mapping[str, date]
    studies: Mapping[str, StudyRecord]

    def collected(self, study: str) -> date | None:
        """Get the day the model collected a study, `None` if not known."""
        record = self.studies.get(study)
        return None if record is None else record.collected

    def new(self, points: pd.DataFrame) -> pd.Series:
        """Get which points are of a study collected after the day of their fit.

        Args:
            points: data points with the columns `fit` and `experiment`.

        Returns:
            True for a point whose study the model collected after the day of
            its fit; False for a point without a day of collection or fit.
        """
        flags = [
            (day := self.collected(str(study))) is not None
            and str(fit) in self.fits
            and day > self.fits[str(fit)]
            for fit, study in zip(points["fit"], points["experiment"], strict=True)
        ]
        return pd.Series(flags, index=points.index, dtype=bool)


def _merged(first: StudyRecord, second: StudyRecord) -> StudyRecord:
    """Get one record of a study in two data folders: the earlier collection.

    A record without a day makes the day unknown, since the study may have
    been seen on any day.
    """
    if first.collected is None or second.collected is None:
        return first if first.collected is None else second
    return first if first.collected <= second.collected else second


def chronology(model: Project) -> Chronology:
    """Get when the adopted fits of a model ran and when it collected its studies.

    Args:
        model: the model.

    Returns:
        The days of the adopted fits of the model (`Project.adopted_fits`) and
        the records of the study folders of every data path of the model.

    Raises:
        ValueError: if a record is not readable, see
            `pkpdsim.data.read_study_record`.
    """
    studies: dict[str, StudyRecord] = {}
    for name, folders in study_folders(model).items():
        for folder in folders:
            record = read_study_record(folder)
            if record is None:
                continue
            earlier = studies.get(name)
            studies[name] = record if earlier is None else _merged(earlier, record)
    return Chronology(
        fits={
            fit_subset(fit_id): fit_date(fit_id)
            for fit_id in model.adopted_fits().values()
        },
        studies=studies,
    )


def _count(n: int, singular: str, plural: str | None = None) -> str:
    """Get a number of things in words, e.g. `1 point` and `2 points`."""
    return f"{n} {singular if n == 1 else plural or singular + 's'}"


def temporal_record(
    cou: ContextOfUse, decisive: pd.DataFrame, chronology: Chronology
) -> TemporalRecord | None:
    """Get the decisive points of a context which are new for their fit.

    Args:
        cou: the context of use.
        decisive: its decisive points with the columns `fit` and `experiment`.
        chronology: when the adopted fits ran and the model collected its
            studies.

    Returns:
        The record, `None` for a context decided on calibration data, which
        are never new for their fit.
    """
    if cou.decisive_role is not EvidenceRole.VALIDATION:
        return None
    new = chronology.new(decisive)
    points = decisive.loc[new]
    studies: list[TemporalStudy] = []
    for (fit, study), group in points.groupby(["fit", "experiment"], sort=True):
        record = chronology.studies[str(study)]
        studies.append(
            TemporalStudy(
                study=str(study),
                fit=str(fit),
                fit_date=chronology.fits[str(fit)].isoformat(),
                collected=str(record.collected),
                released=None if record.released is None else str(record.released),
                pkdb_id=record.pkdb_id,
                n_points=len(group),
                n_evaluable=len(applicable(group)),
            )
        )
    old = {str(study) for study in decisive.loc[~new, "experiment"]}
    undated = sorted(s for s in old if chronology.collected(s) is None)
    return TemporalRecord(
        studies=studies,
        n_points=len(points),
        n_evaluable=len(applicable(points)),
        acceptance=[
            AcceptanceRecord(name=a.name, value=a.value, limit=a.limit, met=a.met)
            for a in (rule.evaluate(points) for rule in cou.rules)
        ],
        n_known=len(old) - len(undated),
        undated=undated,
    )


def temporal_statement(temporal: TemporalRecord, minimum: int) -> str:
    """Describe the acceptance on the studies collected after their fit.

    Args:
        temporal: the record of a context.
        minimum: the evaluable points the verdict of the context needs.

    Returns:
        The statement the factor output comparison ends with.
    """
    before = f"{temporal['n_known']} collected before" + (
        f", {len(temporal['undated'])} without a date" if temporal["undated"] else ""
    )
    if not temporal["studies"]:
        return f"No decisive study was collected after its fit ({before})."
    rules = temporal["acceptance"]
    met = sum(a["met"] is True for a in rules)
    n_studies = len({t["study"] for t in temporal["studies"]})
    few = (
        f", fewer than the minimum of {minimum}"
        if temporal["n_evaluable"] < minimum
        else ""
    )
    return (
        f"Studies collected after their fit: "
        f"{_count(n_studies, 'study', 'studies')} ({before}), "
        f"{_count(temporal['n_evaluable'], 'evaluable point')}{few}, "
        f"{met} of {_count(len(rules), 'rule')} met."
    )


def _write(value: float) -> float | None:
    """Get a number for JSON, `None` for `NaN`."""
    return None if math.isnan(value) else value


def _read(value: float | None) -> float:
    """Get a number from JSON, `NaN` for `None`."""
    return math.nan if value is None else float(value)


def _to_dict(temporal: TemporalRecord) -> dict[str, Any]:
    """Get a temporal record as plain data, `NaN` as `None`."""
    return {
        **temporal,
        "studies": [dict(t) for t in temporal["studies"]],
        "acceptance": [
            {**a, "value": _write(a["value"])} for a in temporal["acceptance"]
        ],
        "undated": list(temporal["undated"]),
    }


def section_record(section: Mapping[str, Any] | None) -> TemporalRecord | None:
    """Get the temporal record of the section `temporal` of the record of a context.

    Args:
        section: the section, `{"temporal": {...}}`, also of a record of
            pkpdsim 0.19; `None` without.

    Returns:
        The record, `NaN` for a missing value of a rule; `None` without.
    """
    data = None if section is None else section.get("temporal")
    if data is None:
        return None
    return TemporalRecord(
        studies=[
            TemporalStudy(
                study=t["study"],
                fit=t["fit"],
                fit_date=t["fit_date"],
                collected=t["collected"],
                released=t["released"],
                pkdb_id=t["pkdb_id"],
                n_points=t["n_points"],
                n_evaluable=t["n_evaluable"],
            )
            for t in data["studies"]
        ],
        n_points=data["n_points"],
        n_evaluable=data["n_evaluable"],
        acceptance=[
            AcceptanceRecord(
                name=a["name"], value=_read(a["value"]), limit=a["limit"], met=a["met"]
            )
            for a in data["acceptance"]
        ],
        n_known=data["n_known"],
        undated=list(data["undated"]),
    )


def _temporal_met(temporal: TemporalRecord | None) -> bool | None:
    """Get whether every rule is met on the new points, `None` if none is evaluable."""
    if temporal is None:
        return None
    results = [a["met"] for a in temporal["acceptance"] if a["met"] is not None]
    return all(results) if results else None


class TemporalExtension(EvidenceExtension):
    """The temporal validation of the contexts of use, see the module.

    Prepared: the `Chronology` of the model.
    """

    name = "temporal"
    scorecard_names = ("temporal_points", "temporal_met")

    @override
    def prepare(
        self, project: Project, credibility: Credibility, fits: FitTable
    ) -> object:
        """Get the chronology of the model."""
        return chronology(project)

    @override
    def assess_context(
        self,
        cou: ContextOfUse,
        decisive: pd.DataFrame,
        fits: FitTable,
        prepared: object,
    ) -> ContextContribution | None:
        """Split the decisive points by the chronology, none without it."""
        if not isinstance(prepared, Chronology):
            return None
        temporal = temporal_record(cou, decisive, prepared)
        if temporal is None:
            return None
        return ContextContribution(
            section={"temporal": _to_dict(temporal)},
            factors={
                "output_comparison": FactorStatement(
                    temporal_statement(temporal, cou.minimum_points)
                )
            },
        )

    @override
    def context_report_blocks(
        self, context: Mapping[str, Any], section: Mapping[str, Any] | None
    ) -> Sequence[ContextBlock]:
        """Show the temporal validation after the acceptance of a context."""
        temporal = section_record(section)
        if temporal is None:
            return ()
        html = render_macro("temporal", temporal, context.get("minimum_points"))
        return (ContextBlock(html, after="acceptance"),)

    @override
    def report_tables(self, record: Mapping[str, Any]) -> Mapping[str, pd.DataFrame]:
        """Write `temporal.tsv`: one row per context of use and study collected after its fit."""
        rows = [
            {"context": c["id"], **study}
            for c in record["contexts"]
            if (t := section_record(c.get("extensions", {}).get(self.name))) is not None
            for study in t["studies"]
        ]
        return {"temporal": pd.DataFrame(rows, columns=["context", *TEMPORAL_COLUMNS])}

    @override
    def scorecard_columns(
        self, context: Mapping[str, Any] | None
    ) -> Mapping[str, object]:
        """Get the evaluable new points and whether every rule is met on them."""
        temporal = section_record(context)
        return {
            "temporal_points": None if temporal is None else temporal["n_evaluable"],
            "temporal_met": _temporal_met(temporal),
        }

    @override
    def regressions(
        self, baseline: Mapping[str, object], current: Mapping[str, object]
    ) -> list[ScorecardChange]:
        """Report the rules met on the new studies before which no longer are."""
        # a baseline of pkpdsim 0.18 or earlier has no temporal column
        if scorecard_true(baseline.get("temporal_met")) and not scorecard_true(
            current.get("temporal_met")
        ):
            return [
                ScorecardChange(
                    what="temporal_met",
                    baseline="True",
                    current=str(current.get("temporal_met")),
                )
            ]
        return []


#: the evidence extension of the temporal validation, registered by
#: `pkpdsim.model.PKPDModel`
TEMPORAL: Final = TemporalExtension()
