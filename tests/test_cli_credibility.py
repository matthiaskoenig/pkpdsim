"""The command `run_credibility`, a wrapper of the kind `credibility` of sbmlsim.

The command keeps its options and output paths: the outputs of the kind
(`credibility/` of its output) are moved up into `--output_dir`, next to the
records of the run, `--from-json` is the input `record` and the lock of
`--lock` is moved into the folder of the model.
"""

import dataclasses
import shutil
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest
from sbmlsim.credibility import CredibilityAssessment
from sbmlsim.credibility import lock as locks
from sbmlsim.credibility.cross_validation import (
    CROSS_VALIDATION_FILE,
    read_cross_validations,
)
from sbmlsim.credibility.report import VERIFICATION_COLUMNS

from pkpdsim.cli import credibility_command, fit_command
from pkpdsim.model import PKPDModel
from pkpdsim.verify import verify

type AsModel = Callable[[PKPDModel], PKPDModel]


@pytest.fixture
def model(credibility_model: PKPDModel, as_model: AsModel) -> PKPDModel:
    """Get the credibility test model as the model object of its package."""
    return as_model(credibility_model)


@pytest.fixture
def lock_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect the lock of every model into the temporary directory."""
    path = tmp_path / locks.LOCK_FILE
    monkeypatch.setattr(locks, "lock_path", lambda project: path)
    return path


def test_the_command_writes_the_record_again(model: PKPDModel, tmp_path: Path) -> None:
    output_dir = credibility_command(
        model, ["-o", str(tmp_path / "assessed"), "--no-verification"]
    )
    again = tmp_path / "again"
    record = output_dir / "credibility.json"
    assert credibility_command(model, ["--from-json", str(record), "-o", str(again)])
    assert (again / "contexts.tsv").read_text(encoding="utf-8") == (
        output_dir / "contexts.tsv"
    ).read_text(encoding="utf-8")
    # the record is written again in its own directory
    assert credibility_command(model, ["--from-json", str(record), "-o", str(again)])
    assert credibility_command(
        model, ["--from-json", str(again / "credibility.json"), "-o", str(again)]
    )
    assert CredibilityAssessment.from_json(again / "credibility.json").contexts


def test_the_command_writes_into_the_credibility_dir(
    model: PKPDModel, as_model: AsModel, tmp_path: Path
) -> None:
    output_dir = credibility_command(
        model, ["-o", str(tmp_path / "assessed"), "--no-verification"]
    )
    into = as_model(model.replace(results_path=tmp_path / "results"))
    record = output_dir / "credibility.json"
    assert credibility_command(into, ["--from-json", str(record)]) == (
        into.credibility_dir
    )
    assert (into.credibility_dir / "index.html").exists()


def test_the_command_replaces_an_assessment_and_nothing_else(
    model: PKPDModel, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "out"
    credibility_command(model, ["-o", str(out), "--no-verification"])
    (out / "stale.tsv").write_text("x", encoding="utf-8")
    credibility_command(model, ["-o", str(out), "--no-verification"])
    assert not (out / "stale.tsv").exists()
    other = tmp_path / "other"
    other.mkdir()
    (other / "notes.txt").write_text("keep", encoding="utf-8")
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(model, ["-o", str(other), "--no-verification"])
    assert exit_info.value.code == 1
    assert "holdsnocredibility.json" in "".join(capsys.readouterr().out.split())
    assert [p.name for p in other.iterdir()] == ["notes.txt"]


def test_the_command_refuses_an_output_of_its_inputs(
    model: PKPDModel, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(model, ["-o", str(model.fits_dir)])
    assert exit_info.value.code == 2
    assert "is or contains the input" in capsys.readouterr().err


def test_the_command_reports_a_broken_declaration(
    model: PKPDModel,
    as_model: AsModel,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    broken = as_model(model.replace(credibility=".credibility:BROKEN"))
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(broken, ["-o", str(tmp_path / "broken")])
    assert exit_info.value.code == 1
    assert "CredB" in capsys.readouterr().out
    assert not (tmp_path / "broken").exists()
    # no staging is left behind
    assert list(tmp_path.iterdir()) == []


def test_the_command_reports_a_missing_record(
    model: PKPDModel, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    args = ["--from-json", str(tmp_path / "missing.json"), "-o", str(tmp_path / "out")]
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(model, args)
    assert exit_info.value.code == 1
    # the console folds a long path at its width
    assert "missing.json" in "".join(capsys.readouterr().out.split())


def test_the_command_reports_a_malformed_record(
    model: PKPDModel, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    record = tmp_path / "credibility.json"
    record.write_text('{"contexts": []}', encoding="utf-8")
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(
            model, ["--from-json", str(record), "-o", str(tmp_path / "o")]
        )
    assert exit_info.value.code == 1
    out = "".join(capsys.readouterr().out.split())
    assert "credibility.json" in out
    assert "isnoassessment" in out


def test_the_command_reports_a_declaration_of_another_type(
    model: PKPDModel,
    as_model: AsModel,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    other = as_model(model.replace(credibility=".fitting:DOSE"))
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(other, ["-o", str(tmp_path / "other")])
    assert exit_info.value.code == 1
    assert "not a Credibility" in "".join(capsys.readouterr().out.split("\n"))


def test_the_command_reports_a_malformed_record_of_a_fit(
    model: PKPDModel,
    as_model: AsModel,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    fits = tmp_path / "fits"
    shutil.copytree(model.fits_dir, fits)
    (fit_id,) = model.adopted_ids
    path = fits / fit_id / "metrics_mappings.tsv"
    pd.read_csv(path, sep="\t").drop(columns="kind").to_csv(path, sep="\t", index=False)
    malformed = as_model(model.replace(fits_dir=fits))
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(malformed, ["-o", str(tmp_path / "out")])
    assert exit_info.value.code == 1
    out = "".join(capsys.readouterr().out.split())
    assert "metrics_mappings.tsv" in out
    assert "adoptthefitagain" in out


@pytest.mark.parametrize(
    "args",
    [
        ["--source-commit", "0123abc"],
        ["--no-verification"],
        ["--lock"],
        ["--uncertainty"],
        ["--variability"],
    ],
)
def test_the_options_of_an_assessment_are_refused_with_a_record(
    model: PKPDModel,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    args: list[str],
) -> None:
    record = tmp_path / "credibility.json"
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(model, ["--from-json", str(record), *args])
    assert exit_info.value.code == 2
    assert args[0] in capsys.readouterr().err


@pytest.mark.parametrize(
    "args",
    [
        ["--draws", "8"],
        ["--uncertainty", "--draws", "1"],
        ["--subjects", "8"],
        ["--variability", "--subjects", "1"],
    ],
)
def test_the_options_of_the_uncertainty_and_the_variability(
    model: PKPDModel, tmp_path: Path, args: list[str]
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(model, ["-o", str(tmp_path / "out"), *args])
    assert exit_info.value.code == 2
    assert not (tmp_path / "out").exists()


def test_the_command_without_verification(model: PKPDModel, tmp_path: Path) -> None:
    output_dir = credibility_command(
        model, ["-o", str(tmp_path / "out"), "--no-verification"]
    )
    record = CredibilityAssessment.from_json(output_dir / "credibility.json")
    assert record.verification is None
    assert (output_dir / "verification.tsv").read_text(
        encoding="utf-8"
    ).strip() == "\t".join(VERIFICATION_COLUMNS)


def test_the_command_writes_the_lock(
    model: PKPDModel, lock_file: Path, tmp_path: Path
) -> None:
    output = credibility_command(
        model, ["-o", str(tmp_path / "out"), "--no-verification", "--lock"]
    )
    record = CredibilityAssessment.from_json(output / "credibility.json")
    assert record.provenance.lock == "written"
    # the lock is the one of the model, not an output
    assert lock_file.is_file()
    assert not (output / locks.LOCK_FILE).exists()
    assert "the lock of the evidence was written" in (output / "index.html").read_text(
        encoding="utf-8"
    )
    again = credibility_command(
        model, ["-o", str(tmp_path / "again"), "--no-verification"]
    )
    matched = CredibilityAssessment.from_json(again / "credibility.json")
    assert matched.provenance.lock == "matched"


def _assessment(model: str, renal_factor: float = 1.1) -> CredibilityAssessment:
    """Get an assessment of the test package with the ratio pairs as one of a model."""
    from credibility_tables import PROVENANCE, fit_table
    from sbmlsim.credibility import assess_table
    from test_credibility_ratios import _renal

    from pkpdsim.credibility import EXTENSIONS

    fits = fit_table(factors={("CredA", "fm_renal"): renal_factor}, times=True)
    provenance = dataclasses.replace(PROVENANCE, model=model)
    return assess_table(_renal(), fits, provenance, extensions=EXTENSIONS)


def _write(directory: Path, *assessments: CredibilityAssessment) -> Path:
    for assessment in assessments:
        folder = directory / assessment.provenance.model / "credibility"
        folder.mkdir(parents=True)
        assessment.to_json(folder / "credibility.json")
    return directory


def test_the_scorecard_command(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """`pkpdsim-verify scorecard` with the extensions and the qualification of pkpdsim."""
    root = _write(tmp_path / "records", _assessment("a"))
    out = tmp_path / "scorecard.tsv"
    qualified = tmp_path / "qualification.tsv"
    args = [
        "scorecard",
        str(root),
        "--out",
        str(out),
        "--qualification",
        str(qualified),
    ]
    assert verify(args) == 0
    printed = capsys.readouterr().out
    assert "1 models, 3 contexts of use" in printed
    assert "a/renal: supported" in printed
    assert "RENAL_IMPAIRMENT: 1 of 1 ratio pairs within the limits of Guest" in printed
    assert pd.read_csv(qualified, sep="\t")["n_pairs"].tolist() == [1, 1]
    table = pd.read_csv(out, sep="\t")
    assert {"temporal_points", "temporal_met", "inputs_differ", "Guest"} <= set(
        table.columns
    )
    # a worse renal context of a, and b which the baseline does not know
    worse = _write(tmp_path / "worse", _assessment("a", 2.0), _assessment("b"))
    assert verify(["scorecard", str(worse), "--baseline", str(out)]) == 1
    printed = capsys.readouterr().out
    assert "regressions against" in printed
    assert "a/renal: verdict supported -> not_supported" in printed
    assert verify(["scorecard", str(tmp_path / "nothing")]) == 2


def test_the_cross_validation_command(
    credibility_model: PKPDModel, tmp_path: Path
) -> None:
    """`run_fit --cross-validation` stores the result in the record of the adopted fit."""
    fits = tmp_path / "fits"
    shutil.copytree(credibility_model.fits_dir, fits)
    model = credibility_model.replace(fits_dir=fits)
    fit_id = model.adopted_ids[0]
    args = ["--cross-validation", fit_id, "--experiments", "CredB", "-o", str(tmp_path)]
    assert fit_command(model, args) == {}
    (result,) = read_cross_validations(model).values()
    assert [s.study for s in result.studies] == ["CredB"]
    assert (fits / fit_id / CROSS_VALIDATION_FILE).is_file()
    with pytest.raises(SystemExit) as exc:
        fit_command(model, ["--cross-validation", fit_id, "--experiments", "CredC"])
    assert exc.value.code == 1


def test_a_lock_which_is_not_moved_is_an_error(
    model: PKPDModel,
    lock_file: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def refuse(source: object, target: object) -> None:
        raise PermissionError("read-only")

    monkeypatch.setattr(shutil, "move", refuse)
    out = tmp_path / "out"
    with pytest.raises(SystemExit) as exit_info:
        credibility_command(model, ["-o", str(out), "--no-verification", "--lock"])
    assert exit_info.value.code == 1
    printed = "".join(capsys.readouterr().out.split())
    assert "isnotwritten(read-only)" in printed
    assert not lock_file.exists()
    # the assessment is written, its lock is left in the output
    assert (out / "credibility.json").is_file()
    assert (out / "credibility.lock.json").is_file()
