"""The temporal validation: the studies a model collected after its fit.

A study which the model collected after the day an adopted fit ran cannot have
shaped the fit, neither as training data nor by the choice of the model while
the modeller looked at it. Its decisive points test whether the acceptance of
a context of use holds on data which did not exist for the model when it was
fitted, the closest the existing models get to a prospective validation.

The day of the collection is the record of the study in the data of the model
(`pkpdsim.data.StudyRecord`, `.<study>_study.json`), which
`pkpdsim.data.write_study_tables` writes on the first collection and keeps
afterwards; the day of a fit is the date of its id (`pkpdsim.fit.fit_date`).
The release of a study in PK-DB is no such day, since a model can collect a
study from a checkout of pkdb_data before it is released. `chronology` reads
both for a model, and the assessment splits the decisive points of every
context of use above C1 by it (`TemporalRecord`).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, TypedDict

import pandas as pd

from pkpdsim.data import StudyRecord, read_study_record
from pkpdsim.fit import fit_date, fit_subset

if TYPE_CHECKING:
    from datetime import date

    from pkpdsim.credibility.assessment import AcceptanceRecord
    from pkpdsim.model import PKPDModel

__all__ = [
    "TEMPORAL_COLUMNS",
    "Chronology",
    "TemporalRecord",
    "TemporalStudy",
    "chronology",
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


def chronology(model: PKPDModel) -> Chronology:
    """Get when the adopted fits of a model ran and when it collected its studies.

    Args:
        model: the model.

    Returns:
        The days of the adopted fits of `model.fits` and the records of the
        study folders of every data path of the model.

    Raises:
        ValueError: if a record is not readable, see
            `pkpdsim.data.read_study_record`.
    """
    studies: dict[str, StudyRecord] = {}
    for data_path in model.data_paths:
        if not data_path.is_dir():
            continue
        for folder in sorted(p for p in data_path.iterdir() if p.is_dir()):
            record = read_study_record(folder)
            if record is None:
                continue
            earlier = studies.get(record.study)
            studies[record.study] = (
                record if earlier is None else _merged(earlier, record)
            )
    return Chronology(
        fits={fit_subset(fit_id): fit_date(fit_id) for fit_id in model.fits},
        studies=studies,
    )
