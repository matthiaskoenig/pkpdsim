"""The outputs, the report and the command of the credibility assessment."""

import dataclasses
import re
from pathlib import Path

import pandas as pd
import pytest
from credibility_tables import fit_table
from pkpdsim_credibility.credibility import CREDIBILITY

from pkpdsim.cli import credibility_command
from pkpdsim.credibility import (
    CredibilityAssessment,
    Provenance,
    Verdict,
    assess_table,
)
from pkpdsim.credibility.report import print_summary, render_html, write_outputs
from pkpdsim.model import PKPDModel

FILES = {"credibility.json", "contexts.tsv", "evidence.tsv", "domain.tsv", "index.html"}

PROVENANCE = Provenance(
    model="credibility",
    title="Toy model with every kind of evidence",
    fits=["DECAY_fit"],
    versions={"pkpdsim": "0.0.0"},
    commit="0123abc",
    dirty=False,
    date="2026-11-03T12:00:00+01:00",
)


@pytest.fixture(scope="module")
def assessment() -> CredibilityAssessment:
    """Get the assessment of the test package with a context of every verdict."""
    assessed = assess_table(CREDIBILITY, fit_table(), PROVENANCE)
    verdicts = {c.verdict for c in assessed.contexts}
    extra = [
        dataclasses.replace(assessed.contexts[1], id=f"other{index}", verdict=verdict)
        for index, verdict in enumerate(v for v in Verdict if v not in verdicts)
    ]
    return dataclasses.replace(assessed, contexts=[*assessed.contexts, *extra])


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def test_write_outputs(assessment: CredibilityAssessment, tmp_path: Path) -> None:
    output_dir = tmp_path / "credibility"
    assert write_outputs(assessment, output_dir) == output_dir
    assert {p.name for p in output_dir.iterdir()} == FILES
    assert (
        CredibilityAssessment.from_json(output_dir / "credibility.json").to_dict()
        == assessment.to_dict()
    )

    contexts = _tsv(output_dir / "contexts.tsv")
    assert list(contexts["id"]) == [c.id for c in assessment.contexts]
    assert list(contexts["verdict"]) == [c.verdict.value for c in assessment.contexts]
    assert list(contexts.columns[:9]) == [
        "id",
        "title",
        "level",
        "influence",
        "consequence",
        "risk",
        "verdict",
        "n_decisive",
        "n_evaluable",
    ]
    dose = contexts.set_index("id").loc["dose"]
    assert dose["within_2fold_limit"] == ">= 80"
    assert dose["GMFE_limit"] == "<= 1.5"
    assert float(dose["GMFE"]) == pytest.approx(1.0)

    evidence = _tsv(output_dir / "evidence.tsv")
    n_evidence = len(assessment.evidence) + len(assessment.limitations)
    assert len(evidence) == n_evidence
    assert {"context", "name", "role", "n", "GMFE", "AFE", "within_2fold"} <= set(
        evidence.columns
    )
    (limitation,) = evidence[evidence["role"] == "limitation"].itertuples()
    assert limitation.context == ""

    domain = _tsv(output_dir / "domain.tsv")
    assert list(domain.columns) == [
        "context",
        "field",
        "value",
        "claimed",
        "coverage",
        "n",
    ]
    assert len(domain) == len(assessment.coverage)

    html = (output_dir / "index.html").read_text(encoding="utf-8")
    for context in assessment.contexts:
        assert context.id in html
    for verdict in Verdict:
        assert f">{verdict.label}<" in html
    assert "CredB" in html and "fm_zero" in html  # the unclaimed validation data
    assert "0123abc" in html  # the provenance
    external = re.compile(r"<(?:script|link)\b[^>]*(?:src|href)=\"?https?://", re.I)
    assert not external.search(html)


def test_the_verdict_is_a_word_with_a_colour(assessment: CredibilityAssessment) -> None:
    html = render_html(assessment)
    for verdict in Verdict:
        assert f'class="verdict {verdict.value}">{verdict.label}<' in html


def test_write_outputs_replaces_its_own_directory(
    assessment: CredibilityAssessment, tmp_path: Path
) -> None:
    write_outputs(assessment, tmp_path)
    (tmp_path / "stale.txt").write_text("of the first run", encoding="utf-8")
    (tmp_path / "figures").mkdir()
    fewer = dataclasses.replace(assessment, contexts=assessment.contexts[:-1])
    write_outputs(fewer, tmp_path)
    assert {p.name for p in tmp_path.iterdir()} == FILES
    assert len(_tsv(tmp_path / "contexts.tsv")) == len(assessment.contexts) - 1


def test_write_outputs_refuses_a_foreign_directory(
    assessment: CredibilityAssessment, tmp_path: Path
) -> None:
    (tmp_path / "other.txt").write_text("not an assessment", encoding="utf-8")
    with pytest.raises(FileExistsError, match=r"holds no credibility\.json"):
        write_outputs(assessment, tmp_path)
    assert {p.name for p in tmp_path.iterdir()} == {"other.txt"}
    file = tmp_path / "other.txt"
    with pytest.raises(FileExistsError):
        write_outputs(assessment, file)


def test_write_outputs_into_an_empty_directory(
    assessment: CredibilityAssessment, tmp_path: Path
) -> None:
    write_outputs(assessment, tmp_path)
    assert {p.name for p in tmp_path.iterdir()} == FILES


def test_print_summary(
    assessment: CredibilityAssessment, capsys: pytest.CaptureFixture[str]
) -> None:
    print_summary(assessment)
    out = capsys.readouterr().out
    for context in assessment.contexts:
        assert context.id in out
        assert context.verdict.label in out


def test_the_command(credibility_model: PKPDModel, tmp_path: Path) -> None:
    output_dir = tmp_path / "assessed"
    assert credibility_command(credibility_model, ["-o", str(output_dir)]) == output_dir
    assert (output_dir / "credibility.json").exists()
    assert {p.name for p in output_dir.iterdir()} == FILES
    again = tmp_path / "again"
    assert (
        credibility_command(
            credibility_model,
            [
                "--from-json",
                str(output_dir / "credibility.json"),
                "--output_dir",
                str(again),
            ],
        )
        == again
    )
    assert (again / "contexts.tsv").read_text(encoding="utf-8") == (
        output_dir / "contexts.tsv"
    ).read_text(encoding="utf-8")


def test_the_command_writes_into_the_credibility_dir(
    credibility_model: PKPDModel, tmp_path: Path
) -> None:
    record = tmp_path / "credibility.json"
    write_outputs(assess_table(CREDIBILITY, fit_table(), PROVENANCE), tmp_path)
    model = dataclasses.replace(credibility_model, results_path=tmp_path / "results")
    assert credibility_command(model, ["--from-json", str(record)]) == (
        model.credibility_dir
    )
    assert (model.credibility_dir / "index.html").exists()


def test_the_command_reports_a_broken_declaration(
    credibility_model: PKPDModel, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    broken = dataclasses.replace(credibility_model, credibility=".credibility:BROKEN")
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(broken, ["-o", str(tmp_path / "broken")])
    assert exit_info.value.code == 1
    assert "CredB" in capsys.readouterr().out
    assert not (tmp_path / "broken").exists()


def test_the_command_reports_a_missing_record(
    credibility_model: PKPDModel, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(
            credibility_model,
            [
                "--from-json",
                str(tmp_path / "missing.json"),
                "-o",
                str(tmp_path / "out"),
            ],
        )
    assert exit_info.value.code == 1
    # the console folds a long path at its width
    assert "missing.json" in "".join(capsys.readouterr().out.split())
