"""The provenance of the study data of an assessment (`pkpdsim.credibility.studies`), the extension `studies`.

The hand-built table of the test package has points of the studies `CredA` and
`CredB`; the tests write a data folder with a record of `CredA` only.
"""

import json
from datetime import date
from pathlib import Path

from credibility_tables import assess, fit_table
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.credibility import CREDIBILITY
from sbmlsim.credibility import CredibilityAssessment
from sbmlsim.credibility.report import write_outputs

from pkpdsim.credibility import EXTENSIONS, PK_BATTERY, STUDIES
from pkpdsim.credibility.studies import section_studies, study_data, study_folders
from pkpdsim.data import StudyRecord, study_record_path, table_digests
from pkpdsim.model import PKPDModel


def _model(data: Path) -> PKPDModel:
    """Get the model of the test package with a data folder of the record of CredA."""
    folder = data / "CredA"
    folder.mkdir(parents=True)
    (folder / ".CredA_Fig1.tsv").write_text("time\tvalue\n0\t1\n", encoding="utf-8")
    record = StudyRecord(
        study="CredA",
        collected=date(2026, 9, 1),
        dating="git",
        dating_commit="0123456789abcdef",
        updated=date(2026, 9, 15),
        released=date(2026, 9, 20),
        pkdb_id="PKDB09999",
        source_digest="feedbeef" * 8,
        pkdb_data_commit="abcdef0123456789" * 2 + "01234567",
        pkdb_data_dirty=False,
        tables=table_digests(folder),
        versions={"pkpdsim": "0.0.0", "pkdb": "0.0.0"},
    )
    study_record_path(folder).write_text(json.dumps(record.to_dict()))
    # a data path which is no folder is skipped
    return MODEL.replace(data_path=(data, data / "missing"))


def test_the_study_data(tmp_path: Path) -> None:
    model = _model(tmp_path / "data")
    assert study_folders(model) == {"CredA": [tmp_path / "data" / "CredA"]}
    table = fit_table().table
    credA, credB = study_data(model, table)
    assert credA["study"] == "CredA"
    assert credA["folder"] == (tmp_path / "data" / "CredA").as_posix()
    assert credA["fits"] == ["DECAY"]
    assert credA["n_points"] == int((table["experiment"] == "CredA").sum())
    assert credA["recorded"]
    assert (credA["collected"], credA["dating"], credA["updated"]) == (
        "2026-09-01",
        "git",
        "2026-09-15",
    )
    assert (credA["n_tables"], credA["changed"]) == (1, [])
    assert credB["folder"] is None
    assert not credB["recorded"]
    assert credB["changed"] is None
    # a table changed by hand after the collection
    (tmp_path / "data" / "CredA" / ".CredA_Fig1.tsv").write_text("edited")
    assert study_data(model, table)[0]["changed"] == [".CredA_Fig1.tsv"]


def test_the_record_and_the_outputs(tmp_path: Path) -> None:
    model = _model(tmp_path / "data")
    assessment = assess(
        CREDIBILITY,
        fit_table(),
        prepared={"studies": study_data(model, fit_table().table)},
    )
    path = tmp_path / "credibility.json"
    assessment.to_json(path)
    read = CredibilityAssessment.from_json(path)
    assert read.to_dict() == assessment.to_dict()
    studies = section_studies(read.extensions.get("studies"))
    assert [s["study"] for s in studies] == ["CredA", "CredB"]
    # a record of pkpdsim 0.18 has no study data (nor the sections)
    data = json.loads(path.read_text(encoding="utf-8"))
    del data["extensions"]
    assert "studies" not in CredibilityAssessment.from_dict(data).extensions
    # a record of pkpdsim 0.19 has them as a field of its own
    data["studies"] = assessment.extensions["studies"]["studies"]
    legacy = CredibilityAssessment.from_dict(data).extensions["studies"]
    assert section_studies(legacy) == studies
    table = STUDIES.report_tables(assessment.to_dict())["studies"]
    assert list(table["study"]) == ["CredA", "CredB"]
    assert table.loc[0, "fits"] == "DECAY"
    output = write_outputs(assessment, tmp_path / "out", EXTENSIONS, PK_BATTERY)
    assert (
        (output / "studies.tsv")
        .read_text(encoding="utf-8")
        .startswith("study\tfolder\tfits")
    )
    html = (output / "index.html").read_text(encoding="utf-8")
    assert 'id="studies"' in html
    assert "2 studies, 1 with a record of their collection" in html
    assert "no data folder of the study" in html
    assert "abcdef012345" in html
