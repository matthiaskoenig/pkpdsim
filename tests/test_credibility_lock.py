"""The lock of the evidence of a declaration (`pkpdsim.credibility.lock`).

The lock of a model lies in the folder of its package; the tests redirect it
into their temporary directory, so no test writes into the test package.
"""

import dataclasses
import json
from pathlib import Path

import pytest
from credibility_tables import fit_table
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim.fit.filters import named_filter

from pkpdsim.cli import credibility_command
from pkpdsim.credibility import (
    Credibility,
    CredibilityAssessment,
    Evidence,
    EvidenceRole,
    LockMismatchError,
    assess,
    check_declaration,
    select,
)
from pkpdsim.credibility import lock as lock_module
from pkpdsim.credibility.lock import (
    LOCK_FILE,
    LOCK_REMEDY,
    lock_differences,
    read_lock,
    selection_lock,
    write_lock,
)
from pkpdsim.model import PKPDModel


@pytest.fixture
def lock_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect the lock of every model into the temporary directory."""
    path = tmp_path / LOCK_FILE
    monkeypatch.setattr(lock_module, "lock_path", lambda model: path)
    return path


def _with_dose_subjects() -> Credibility:
    """Get the declaration whose context dose also validates on fm_subject."""
    contexts = []
    for cou in CREDIBILITY.contexts:
        if cou.id == "dose":
            extra = Evidence(
                "subjects",
                EvidenceRole.VALIDATION,
                named_filter({"CredA": {"fm_subject"}}),
            )
            cou = dataclasses.replace(cou, evidence=(cou.evidence[0], extra))
        contexts.append(cou)
    return dataclasses.replace(CREDIBILITY, contexts=tuple(contexts))


def test_the_lock_of_a_selection() -> None:
    fits = fit_table()
    lock = selection_lock(select(CREDIBILITY, fits), fits)
    assert lock["renal/renal impairment"] == ["DECAY/CredA/fm_renal"]
    assert lock["/scaled data"] == ["DECAY/CredA/fm_limit"]
    assert list(lock) == sorted(lock)


def test_write_and_compare(lock_file: Path) -> None:
    fits = fit_table()
    selections = select(CREDIBILITY, fits)
    assert lock_differences(MODEL, selections, fits) is None
    assert write_lock(MODEL, selections, fits) == lock_file
    data = json.loads(lock_file.read_text(encoding="utf-8"))
    assert data["version"] == 1
    assert read_lock(MODEL) == data["evidence"]
    assert lock_differences(MODEL, selections, fits) == []
    changed = select(_with_dose_subjects(), fits)
    differences = lock_differences(MODEL, changed, fits)
    assert differences is not None
    assert "the evidence 'dose/subjects' is not in the lock" in differences


def test_a_gone_and_a_changed_evidence(lock_file: Path) -> None:
    fits = fit_table()
    write_lock(MODEL, select(CREDIBILITY, fits), fits)
    data = json.loads(lock_file.read_text(encoding="utf-8"))
    data["evidence"]["renal/renal impairment"] = ["DECAY/CredA/fm_other"]
    data["evidence"]["renal/old"] = ["DECAY/CredA/fm_renal"]
    lock_file.write_text(json.dumps(data))
    differences = lock_differences(MODEL, select(CREDIBILITY, fits), fits)
    assert differences == [
        "the evidence 'renal/old' of the lock is gone",
        "the evidence 'renal/renal impairment' selects other mappings: "
        "+DECAY/CredA/fm_renal, -DECAY/CredA/fm_other",
    ]


def test_a_file_which_is_no_lock(lock_file: Path) -> None:
    lock_file.write_text(json.dumps({"version": 0}))
    with pytest.raises(ValueError, match="is no lock of version 1"):
        read_lock(MODEL)


def test_the_assessment_writes_and_matches_the_lock(
    credibility_model: PKPDModel, lock_file: Path
) -> None:
    written = assess(credibility_model, verification=False, lock=True)
    assert written.provenance.lock == "written"
    assert lock_file.is_file()
    matched = assess(credibility_model, verification=False)
    assert matched.provenance.lock == "matched"
    assert check_declaration(credibility_model) == []
    data = json.loads(lock_file.read_text(encoding="utf-8"))
    data["evidence"]["renal/renal impairment"] = []
    lock_file.write_text(json.dumps(data))
    with pytest.raises(LockMismatchError, match="renal/renal impairment"):
        assess(credibility_model, verification=False)
    problems = check_declaration(credibility_model)
    assert problems[-1] == LOCK_REMEDY
    assert "selects other mappings: +DECAY/CredA/fm_renal" in problems[0]


def test_without_a_lock(credibility_model: PKPDModel, lock_file: Path) -> None:
    assessment = assess(credibility_model, verification=False)
    assert assessment.provenance.lock is None
    record = assessment.provenance.to_dict()
    del record["lock"]
    assert type(assessment.provenance).from_dict(record).lock is None


def test_the_command_writes_the_lock(
    credibility_model: PKPDModel, lock_file: Path, tmp_path: Path
) -> None:
    output = credibility_command(
        credibility_model, ["-o", str(tmp_path / "out"), "--no-verification", "--lock"]
    )
    record = CredibilityAssessment.from_json(output / "credibility.json")
    assert record.provenance.lock == "written"
    assert "the lock of the evidence was written" in (output / "index.html").read_text(
        encoding="utf-8"
    )
    with pytest.raises(SystemExit) as exc:
        credibility_command(
            credibility_model,
            ["--from-json", str(output / "credibility.json"), "--lock"],
        )
    assert exc.value.code == 2
