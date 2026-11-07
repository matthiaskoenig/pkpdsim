"""The provenance of the study data an assessment was computed on.

The simulations of a model run on the tables of its studies in its data
folders, which the collection writes from a checkout of pkdb_data together
with the record of the study (`pkpdsim.data.StudyRecord`): the days of the
first and the latest collection, the release in PK-DB, the source digest, the
commit of pkdb_data and the sha256 of every table. `study_data` gives, for
every study with data points in the adopted fits, its record and whether the
tables still are the ones the collection wrote, so an assessment names the
data it was computed on and a table changed by hand is visible.

The provenance is the evidence extension `STUDIES` (`StudiesExtension`, name
`studies`) of the credibility of sbmlsim: the section `studies` of the record
of an assessment, the table `studies.tsv` and the section Study data of the
report.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final, TypedDict, override

import pandas as pd
from sbmlsim.credibility.extensions import EvidenceExtension, ReportSection

from pkpdsim.credibility.blocks import render_macro
from pkpdsim.data import StudyRecord, read_study_record

if TYPE_CHECKING:
    from datetime import date

    from sbmlsim.credibility.declaration import Credibility
    from sbmlsim.credibility.records import FitTable
    from sbmlsim.project import Project

__all__ = [
    "STUDIES",
    "STUDY_COLUMNS",
    "StudiesExtension",
    "StudyData",
    "data_paths",
    "study_data",
    "study_folders",
]

#: the columns of the table of the study data, see `StudyData`
STUDY_COLUMNS = (
    "study",
    "folder",
    "fits",
    "n_points",
    "recorded",
    "collected",
    "dating",
    "updated",
    "released",
    "pkdb_id",
    "source_digest",
    "pkdb_data_commit",
    "pkdb_data_dirty",
    "n_tables",
    "changed",
)


class StudyData(TypedDict):
    """The data of a study the adopted fits of a model have points of.

    Attributes:
        study: the study (the experiment).
        folder: the folder of the study relative to the folder of the model,
            `None` if no data path of the model has it.
        fits: the subsets of the fits with points of the study.
        n_points: the data points of the study in the fits, a point of two
            fits twice.
        recorded: whether the folder has a record of its collection.
        collected: the day of the first collection, ISO format.
        dating: how the day of the first collection is known.
        updated: the day of the latest collection, ISO format.
        released: the day PK-DB released the study, ISO format.
        pkdb_id: the identifier of the release in PK-DB.
        source_digest: the digest of the source files in PK-DB.
        pkdb_data_commit: the commit of pkdb_data of the latest collection.
        pkdb_data_dirty: whether the study folder of that checkout had changes.
        n_tables: the tables of the record.
        changed: the tables which differ from the record (`-` removed, `+`
            added), `None` without the digests of the tables.
    """

    study: str
    folder: str | None
    fits: list[str]
    n_points: int
    recorded: bool
    collected: str | None
    dating: str | None
    updated: str | None
    released: str | None
    pkdb_id: str | None
    source_digest: str | None
    pkdb_data_commit: str | None
    pkdb_data_dirty: bool | None
    n_tables: int
    changed: list[str] | None


def data_paths(model: Project) -> tuple[Path, ...]:
    """Get the data paths of a model, `Project.data_path` as a tuple.

    Args:
        model: the model.

    Returns:
        The data paths, in their order; none without.
    """
    paths = model.data_path
    if paths is None:
        return ()
    if isinstance(paths, Path):
        return (paths,)
    return tuple(paths)


def study_folders(model: Project) -> dict[str, list[Path]]:
    """Get the folders of the studies in the data paths of a model.

    Args:
        model: the model.

    Returns:
        The folders of every study by its name, in the order of the data
        paths, whose first one the experiments read; a data path which is no
        folder is skipped.
    """
    folders: dict[str, list[Path]] = defaultdict(list)
    for data_path in data_paths(model):
        if not data_path.is_dir():
            continue
        for folder in sorted(p for p in data_path.iterdir() if p.is_dir()):
            folders[folder.name].append(folder)
    return dict(folders)


def _relative(folder: Path, model: Project) -> str:
    """Get a folder relative to the folder of the model where it lies below."""
    try:
        return folder.relative_to(model.root).as_posix()
    except ValueError:
        return folder.as_posix()


def study_data(model: Project, points: pd.DataFrame) -> list[StudyData]:
    """Get the provenance of the studies the data points are of.

    Args:
        model: the model.
        points: the data points of the adopted fits, with the columns `fit`
            and `experiment`, see `sbmlsim.credibility.records.evaluate_fits`.

    Returns:
        A record per study of the points, sorted by study.

    Raises:
        ValueError: if a record is not readable, see
            `pkpdsim.data.read_study_record`.
    """
    folders = study_folders(model)
    result: list[StudyData] = []
    for study, group in points.groupby("experiment", sort=True):
        name = str(study)
        folder = folders.get(name, [None])[0]
        record = None if folder is None else read_study_record(folder)
        result.append(
            StudyData(
                study=name,
                folder=None if folder is None else _relative(folder, model),
                fits=sorted({str(fit) for fit in group["fit"]}),
                n_points=len(group),
                **_recorded(record, folder),
            )
        )
    return result


class _Recorded(TypedDict):
    """The fields of `StudyData` from the record of a study."""

    recorded: bool
    collected: str | None
    dating: str | None
    updated: str | None
    released: str | None
    pkdb_id: str | None
    source_digest: str | None
    pkdb_data_commit: str | None
    pkdb_data_dirty: bool | None
    n_tables: int
    changed: list[str] | None


def _recorded(record: StudyRecord | None, folder: Path | None) -> _Recorded:
    """Get the fields of `StudyData` from the record of a study."""
    if record is None or folder is None:
        return _Recorded(
            recorded=False,
            collected=None,
            dating=None,
            updated=None,
            released=None,
            pkdb_id=None,
            source_digest=None,
            pkdb_data_commit=None,
            pkdb_data_dirty=None,
            n_tables=0,
            changed=None,
        )
    return _Recorded(
        recorded=True,
        collected=_iso(record.collected),
        dating=record.dating,
        updated=_iso(record.updated),
        released=_iso(record.released),
        pkdb_id=record.pkdb_id,
        source_digest=record.source_digest,
        pkdb_data_commit=record.pkdb_data_commit,
        pkdb_data_dirty=record.pkdb_data_dirty,
        n_tables=len(record.tables),
        changed=record.changed_tables(folder),
    )


def _iso(value: date | None) -> str | None:
    """Get a day in ISO format, `None` for a missing one."""
    return None if value is None else value.isoformat()


def section_studies(section: Mapping[str, Any] | None) -> list[StudyData]:
    """Get the study data of the section `studies` of the record of an assessment.

    Args:
        section: the section, `{"studies": [...]}`, also of a record of
            pkpdsim 0.19; `None` without.

    Returns:
        The study data, empty without.
    """
    if section is None:
        return []
    return [StudyData(**study) for study in section.get("studies", ())]


class StudiesExtension(EvidenceExtension):
    """The provenance of the study data of an assessment, see the module.

    Prepared: the `StudyData` of every study with points in the adopted fits.
    """

    name = "studies"

    @override
    def prepare(
        self, project: Project, credibility: Credibility, fits: FitTable
    ) -> object:
        """Get the study data of the points of the adopted fits."""
        return study_data(project, fits.table)

    @override
    def assessment_section(self, prepared: object) -> Mapping[str, Any] | None:
        """Get the study data as the section of the assessment, none without."""
        if not isinstance(prepared, list):
            return None
        return {"studies": [dict(study) for study in prepared]}

    @override
    def report_section(self, record: Mapping[str, Any]) -> ReportSection | None:
        """Show the section Study data, none without a study."""
        studies = section_studies(record.get("extensions", {}).get(self.name))
        if not studies:
            return None
        return ReportSection(
            id="studies", title="Study data", html=render_macro("studies", studies)
        )

    @override
    def report_tables(self, record: Mapping[str, Any]) -> Mapping[str, pd.DataFrame]:
        """Write `studies.tsv`: one row per study with points in the adopted fits.

        The fits and the changed tables are joined by `;`.
        """
        rows = [
            {
                **study,
                "fits": ";".join(study["fits"]),
                "changed": (
                    None if study["changed"] is None else ";".join(study["changed"])
                ),
            }
            for study in section_studies(record.get("extensions", {}).get(self.name))
        ]
        return {"studies": pd.DataFrame(rows, columns=list(STUDY_COLUMNS))}


#: the evidence extension of the study data, registered by `pkpdsim.model.PKPDModel`
STUDIES: Final = StudiesExtension()
