"""The temporal validation (`pkpdsim.credibility.temporal`).

The hand-built table of the test package has one fit, `DECAY`; the
chronologies date the studies `CredA` and `CredB` around the day of the fit.
"""

import json
from datetime import date, timedelta
from pathlib import Path

import pytest
from credibility_tables import fit_table
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim.fit.result import fit_date, fit_subset

from pkpdsim.credibility import (
    ContextResult,
    CredibilityAssessment,
    Provenance,
    assess,
    assess_table,
)
from pkpdsim.credibility.report import temporal_table, write_outputs
from pkpdsim.credibility.scorecard import compare, scorecard
from pkpdsim.credibility.temporal import Chronology, chronology
from pkpdsim.data import StudyRecord, study_record_path
from pkpdsim.model import PKPDModel

PROVENANCE = Provenance(
    model="credibility",
    title="Toy model with every kind of evidence",
    fits=["DECAY_20261001_120000__abcde"],
    versions={"pkpdsim": "0.0.0"},
    commit=None,
    dirty=None,
    date="2026-11-04T12:00:00+01:00",
)
FIT_DAY = date(2026, 10, 1)


def _chronology(a: date | None, b: date | None = FIT_DAY) -> Chronology:
    """Get a chronology with the days of collection of CredA and CredB."""
    return Chronology(
        fits={"DECAY": FIT_DAY},
        studies={
            "CredA": StudyRecord(
                study="CredA",
                collected=a,
                dating=None if a is None else "collection",
                released=date(2026, 9, 1),
                pkdb_id="PKDB09999",
            ),
            "CredB": StudyRecord(
                study="CredB", collected=b, dating=None if b is None else "collection"
            ),
        },
    )


def _context(assessment: CredibilityAssessment, id: str) -> ContextResult:
    return next(c for c in assessment.contexts if c.id == id)


def _statement(context: ContextResult) -> str:
    return next(
        f["statement"] for f in context.factors if f["id"] == "output_comparison"
    )


def test_the_date_of_a_fit() -> None:
    assert fit_date("PK_X_20261016_185018__3efe3") == date(2026, 10, 16)
    with pytest.raises(ValueError, match="is not the id of a fit"):
        fit_date("PK")
    with pytest.raises(ValueError, match="is no date"):
        fit_date("PK_20261399_185018__3efe3")


def test_the_new_points() -> None:
    table = fit_table().table
    new = _chronology(FIT_DAY + timedelta(days=1)).new(table)
    assert set(table.loc[new, "experiment"]) == {"CredA"}
    # the day of the fit is not after it, and an unknown day never is
    assert not _chronology(FIT_DAY).new(table).any()
    assert not _chronology(None, None).new(table).any()


def test_the_temporal_record_of_a_context() -> None:
    assessment = assess_table(
        CREDIBILITY,
        fit_table(),
        PROVENANCE,
        chronology=_chronology(FIT_DAY + timedelta(days=1)),
    )
    # a description is decided on calibration data, which are never new
    assert _context(assessment, "describe").temporal is None
    dose = _context(assessment, "dose")
    temporal = dose.temporal
    assert temporal is not None
    (study,) = temporal["studies"]
    assert study == {
        "study": "CredA",
        "fit": "DECAY",
        "fit_date": "2026-10-01",
        "collected": "2026-10-02",
        "released": "2026-09-01",
        "pkdb_id": "PKDB09999",
        "n_points": temporal["n_points"],
        "n_evaluable": temporal["n_evaluable"],
    }
    assert temporal["n_evaluable"] > 0
    assert [a["name"] for a in temporal["acceptance"]] == [
        a["name"] for a in dose.acceptance
    ]
    assert all(a["met"] for a in temporal["acceptance"])
    assert "Studies collected after their fit: 1 study" in _statement(dose)
    renal = _context(assessment, "renal")
    assert "fewer than the minimum of 20" in _statement(renal)
    # the verdict does not change
    plain = assess_table(CREDIBILITY, fit_table(), PROVENANCE)
    assert [c.verdict for c in assessment.contexts] == [
        c.verdict for c in plain.contexts
    ]
    assert all(c.temporal is None for c in plain.contexts)


def test_a_rule_which_fails_on_the_new_points() -> None:
    fits = fit_table(factors={("CredA", "fm_dose"): 4.0})
    temporal = _context(
        assess_table(
            CREDIBILITY,
            fits,
            PROVENANCE,
            chronology=_chronology(FIT_DAY + timedelta(days=1)),
        ),
        "dose",
    ).temporal
    assert temporal is not None
    assert any(a["met"] is False for a in temporal["acceptance"])


def test_no_new_study() -> None:
    assessment = assess_table(
        CREDIBILITY, fit_table(), PROVENANCE, chronology=_chronology(None, FIT_DAY)
    )
    dose = _context(assessment, "dose")
    assert dose.temporal is not None
    assert dose.temporal["studies"] == []
    assert dose.temporal["undated"] == ["CredA"]
    assert "No decisive study was collected after its fit" in _statement(dose)
    assert "1 without a date" in _statement(dose)


def test_the_record_and_the_outputs(tmp_path: Path) -> None:
    assessment = assess_table(
        CREDIBILITY,
        fit_table(),
        PROVENANCE,
        chronology=_chronology(FIT_DAY + timedelta(days=1)),
    )
    path = tmp_path / "credibility.json"
    assessment.to_json(path)
    read = CredibilityAssessment.from_json(path)
    assert read.to_dict() == assessment.to_dict()
    # a record of pkpdsim 0.18 has no temporal record
    data = json.loads(path.read_text(encoding="utf-8"))
    for context in data["contexts"]:
        del context["temporal"]
    assert all(
        c.temporal is None for c in CredibilityAssessment.from_dict(data).contexts
    )
    table = temporal_table(assessment)
    assert set(table["context"]) == {"dose", "renal"}
    assert set(table["study"]) == {"CredA"}
    output = write_outputs(assessment, tmp_path / "out")
    assert (
        (output / "temporal.tsv")
        .read_text(encoding="utf-8")
        .startswith("context\tstudy\tfit")
    )
    html = (output / "index.html").read_text(encoding="utf-8")
    assert "Temporal validation" in html
    assert "PKDB09999" in html


def test_the_scorecard() -> None:
    before = scorecard(
        {
            "toy": assess_table(
                CREDIBILITY,
                fit_table(),
                PROVENANCE,
                chronology=_chronology(FIT_DAY + timedelta(days=1)),
            )
        }
    )
    row = before.set_index("context").loc["dose"]
    assert row["temporal_met"] == True  # noqa: E712
    assert row["temporal_points"] > 0
    assert before.set_index("context").loc["describe", "temporal_met"] is None
    after = scorecard(
        {
            "toy": assess_table(
                CREDIBILITY,
                fit_table(factors={("CredA", "fm_dose"): 4.0}),
                PROVENANCE,
                chronology=_chronology(FIT_DAY + timedelta(days=1)),
            )
        }
    )
    regressions = [r for r in compare(before, after) if r.what == "temporal_met"]
    assert [r.context for r in regressions] == ["dose"]
    # a baseline without the temporal columns compares
    old = before.drop(columns=["temporal_points", "temporal_met"])
    assert not [r for r in compare(old, after) if r.what == "temporal_met"]


def _record(study: str, collected: date | None) -> StudyRecord:
    return StudyRecord(
        study=study,
        collected=collected,
        dating=None if collected is None else "collection",
    )


def _write_record(data: Path, record: StudyRecord) -> None:
    folder = data / record.study
    folder.mkdir(parents=True, exist_ok=True)
    study_record_path(folder).write_text(json.dumps(record.to_dict()))


def test_the_chronology_of_a_model(tmp_path: Path) -> None:
    first, second = tmp_path / "a", tmp_path / "b"
    _write_record(first, _record("CredA", date(2026, 10, 5)))
    _write_record(second, _record("CredA", date(2026, 10, 3)))
    _write_record(first, _record("CredB", date(2026, 10, 5)))
    _write_record(second, _record("CredB", None))
    (first / "Other").mkdir()
    model = MODEL.replace(
        adopted_ids=["DECAY_20261001_120000__abcde"],
        data_path=(first, second, tmp_path / "missing"),
    )
    found = chronology(model)
    assert found.fits == {"DECAY": FIT_DAY}
    # the earlier collection, and an unknown day stays unknown
    assert found.collected("CredA") == date(2026, 10, 3)
    assert found.collected("CredB") is None
    assert found.collected("Other") is None


def test_the_assessment_reads_the_chronology(
    credibility_model: PKPDModel, tmp_path: Path
) -> None:
    (fit_id,) = credibility_model.adopted_ids
    later = fit_date(fit_id) + timedelta(days=1)
    _write_record(tmp_path, _record("CredA", later))
    model = credibility_model.replace(data_path=(tmp_path,))
    assessment = assess(model, verification=False)
    dose = _context(assessment, "dose")
    assert dose.temporal is not None
    assert [(s["study"], s["fit"]) for s in dose.temporal["studies"]] == [
        ("CredA", fit_subset(fit_id))
    ]
    # the provenance of the study data names the record of CredA
    studies = {s["study"]: s for s in assessment.studies}
    assert set(studies) == {"CredA", "CredB"}
    assert studies["CredA"]["recorded"]
    assert studies["CredA"]["collected"] == later.isoformat()
    assert not studies["CredB"]["recorded"]
