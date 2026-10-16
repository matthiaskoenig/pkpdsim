"""The reports of the adopted fits of a model."""

import dataclasses
import json
from pathlib import Path

import pytest

from pkpdsim.cli import fit_command, main
from pkpdsim.model import PKPDModel
from pkpdsim.reports import adopted_reports
from pkpdsim.uncertainty import adopt_fit
from pkpdsim_toy import MODEL

FILES = ("index.html", "report.pdf", "report.typ", "references.bib", "references.csl")


@pytest.fixture(scope="module")
def adopted(tmp_path_factory: pytest.TempPathFactory) -> PKPDModel:
    """Get the toy model with two adopted fits of its subset."""
    base = tmp_path_factory.mktemp("adopted")
    fits = base / "fits"
    ids = []
    for index in range(2):
        out = base / f"fit{index}"
        runs = fit_command(
            MODEL,
            [
                "--runs",
                "1",
                "--identifiability",
                "--max-points",
                "5",
                "--no-fisher",
                "--no-mapping-figures",
                "--output_dir",
                str(out),
            ],
        )
        (opid,) = runs
        adopt_fit(out / opid, fits)
        ids.append(opid)
    return dataclasses.replace(
        MODEL, fits=ids, fits_path=fits, results_path=base / "results"
    )


def test_one_report_per_adopted_fit(adopted: PKPDModel) -> None:
    dirs = adopted_reports(adopted)
    assert [d.name for d in dirs] == list(adopted.fits)
    for directory in dirs:
        assert directory.parent == adopted.report_dir
        for name in FILES:
            assert (directory / name).exists(), name
        assert (directory / "report.pdf").read_bytes().startswith(b"%PDF")
        assert (directory / "plots").is_dir()
        typst = (directory / "report.typ").read_text(encoding="utf-8")
        assert adopted.title in typst
        assert "fit DECAY" in typst


def test_the_output_dir(adopted: PKPDModel, tmp_path: Path) -> None:
    dirs = adopted_reports(adopted, tmp_path)
    assert [d.parent for d in dirs] == [tmp_path, tmp_path]


def test_the_command_writes_the_reports(adopted: PKPDModel) -> None:
    from pkpdsim.cli import report_command

    assert report_command(adopted, ["--adopted"]) == adopted.report_dir
    for fit in adopted.fits:
        assert (adopted.report_dir / fit / "index.html").exists()


def test_the_console_script(
    adopted: PKPDModel, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr("pkpdsim.cli.resolve_model", lambda package: adopted)
    main("run_report", ["pkpdsim_toy", "--adopted", "--output_dir", str(tmp_path)])
    for fit in adopted.fits:
        assert (tmp_path / fit / "index.html").exists()


def test_a_record_without_provenance(adopted: PKPDModel, tmp_path: Path) -> None:
    fits = tmp_path / "fits"
    fits.mkdir()
    (fit,) = adopted.fits[:1]
    record = fits / fit
    record.mkdir()
    for source in (adopted.fits_dir / fit).iterdir():
        if source.is_dir():
            continue
        (record / source.name).write_bytes(source.read_bytes())
    path = record / "optimization_result.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data.pop("provenance") is not None
    path.write_text(json.dumps(data), encoding="utf-8")
    model = dataclasses.replace(adopted, fits=[fit], fits_path=fits)
    (directory,) = adopted_reports(model, tmp_path / "out")
    assert "not recorded" in (directory / "report.typ").read_text(encoding="utf-8")


def test_a_fit_of_an_unknown_subset(adopted: PKPDModel, tmp_path: Path) -> None:
    fits = tmp_path / "fits"
    fits.mkdir()
    fit = "OTHER_20260101_000000__abcde"
    (fits / fit).mkdir()
    model = dataclasses.replace(adopted, fits=[fit], fits_path=fits)
    with pytest.raises(ValueError, match=r"OTHER_20260101_000000__abcde.*DECAY"):
        adopted_reports(model, tmp_path / "out")


def test_a_model_without_adopted_fits(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="no adopted fits"):
        adopted_reports(MODEL, tmp_path)
