"""Tests of the report of an adopted fit, recreated from its record."""

import dataclasses
import re
import shutil
from pathlib import Path

import pytest

import pkpdsim_toy
from pkpdsim.cli import report_command, run_report_main
from pkpdsim.fit_report import fit_subset, report_adopted_fit
from pkpdsim_toy import MODEL


@pytest.mark.parametrize(
    ("fit_id", "subset"),
    [
        ("PK_20261016_185018__3efe3", "PK"),
        ("PD_LSQ_20261016_185440__b73fe", "PD_LSQ"),
    ],
)
def test_fit_subset(fit_id: str, subset: str) -> None:
    assert fit_subset(fit_id) == subset


@pytest.mark.parametrize(
    "fit_id",
    [
        "nonsense",
        "PK_20261016_185018",
        "_20261016_185018__3efe3",
        "../PK_20261016_185018__3efe3",
        "a/PK_20261016_185018__3efe3",
        "a\\PK_20261016_185018__3efe3",
        "..PK_20261016_185018__3efe3",
    ],
)
def test_fit_subset_of_no_fit_id_raises(fit_id: str) -> None:
    with pytest.raises(ValueError, match="is no id of a fit"):
        fit_subset(fit_id)


def headings(report: Path) -> set[str]:
    """Get the headings of the HTML report."""
    html = (report / "index.html").read_text(encoding="utf-8")
    return set(re.findall(r"<h[1-6][^>]*>([^<]+)</h[1-6]>", html))


def plots(report: Path) -> set[str]:
    """Get the names of the figures of a report."""
    return {path.name for path in (report / "plots").iterdir()}


def test_report_adopted_fit(fits_dir: Path, toy_fit: Path, tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, fits_path=fits_dir)
    report = report_adopted_fit(model, toy_fit.name, tmp_path)
    assert report == tmp_path / toy_fit.name
    # every figure and section of the report of the fit itself, which was made
    # without the figures of the fit mappings
    assert plots(toy_fit) <= plots(report)
    assert headings(toy_fit) <= headings(report)
    # the runs of the optimization, the figures of the fit mappings, the profiles
    assert {
        "traces.svg",
        "waterfall.svg",
        "fit_ToyDecay_fm_decay_A.svg",
        "profile_k.svg",
    } <= plots(report)
    # the identifiability section with the profiles and the Fisher information
    html = (report / "index.html").read_text(encoding="utf-8")
    assert 'id="identifiability"' in html
    assert ">Profile likelihood<" in html
    assert ">Fisher information<" in html


def test_report_fisher_only(fits_dir: Path, toy_fit: Path, tmp_path: Path) -> None:
    fits = tmp_path / "fits"
    shutil.copytree(fits_dir / toy_fit.name, fits / toy_fit.name)
    (fits / toy_fit.name / "identifiability.json").unlink()
    model = dataclasses.replace(MODEL, fits_path=fits)
    report = report_adopted_fit(model, toy_fit.name, tmp_path / "report")
    html = (report / "index.html").read_text(encoding="utf-8")
    assert 'id="identifiability"' in html
    assert ">Profile likelihood<" not in html
    assert ">Fisher information<" in html


def test_report_adopted_fit_of_no_record_raises(fits_dir: Path, tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, fits_path=fits_dir)
    with pytest.raises(FileNotFoundError, match="DECAY_20261016_185018__3efe3"):
        report_adopted_fit(model, "DECAY_20261016_185018__3efe3", tmp_path)


def test_report_adopted_fit_of_an_unknown_problem_raises(
    fits_dir: Path, toy_fit: Path, tmp_path: Path
) -> None:
    other = toy_fit.name.replace("DECAY", "PK", 1)
    shutil.copytree(fits_dir / toy_fit.name, tmp_path / "fits" / other)
    model = dataclasses.replace(MODEL, fits_path=tmp_path / "fits")
    with pytest.raises(ValueError, match="'PK'"):
        report_adopted_fit(model, other, tmp_path)


def test_run_report_adopted_cli(
    fits_dir: Path, toy_fit: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        pkpdsim_toy, "MODEL", dataclasses.replace(MODEL, fits_path=fits_dir)
    )
    run_report_main(["pkpdsim_toy", "--adopted", toy_fit.name, "-o", str(tmp_path)])
    assert (tmp_path / toy_fit.name / "index.html").exists()


def test_run_report_adopted_into_the_report_dir(
    fits_dir: Path, toy_fit: Path, tmp_path: Path
) -> None:
    model = dataclasses.replace(MODEL, fits_path=fits_dir, results_path=tmp_path)
    report = report_command(model, ["--adopted", toy_fit.name])
    assert report == model.report_dir / toy_fit.name
    assert (report / "index.html").exists()


@pytest.mark.parametrize(
    ("args", "message"),
    [
        (["--adopted", "nonsense"], "'nonsense' is no id of a fit"),
        (["--adopted", "DECAY_20261016_185018__3efe3"], "has no record"),
        (
            ["--adopted", "DECAY_20261016_185018__3efe3", "parameters.json"],
            "unrecognized arguments",
        ),
    ],
)
def test_run_report_adopted_usage_errors(
    fits_dir: Path,
    tmp_path: Path,
    args: list[str],
    message: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    model = dataclasses.replace(MODEL, fits_path=fits_dir)
    with pytest.raises(SystemExit) as exit_info:
        report_command(model, [*args, "-o", str(tmp_path)])
    assert exit_info.value.code == 2
    err = capsys.readouterr().err
    assert err.startswith(f"usage: run_report {MODEL.package}")
    assert message in err


def test_run_report_help_names_adopted(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exit_info:
        report_command(MODEL, ["--help"])
    assert exit_info.value.code == 0
    assert "--adopted FIT_ID" in capsys.readouterr().out
