"""The command which reports the adopted fits of a model, `run_report --adopted`.

The reports themselves are those of the `report` kind of sbmlsim with `adopted`
(`sbmlsim.fit.cli.report_from_records`), tested in sbmlsim.
"""

from pathlib import Path

import pytest
from sbmlsim.fit.records import adopt_fit

from pkpdsim.cli import fit_command, main, report_command
from pkpdsim.model import PKPDModel
from pkpdsim_toy import MODEL


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
    return MODEL.replace(adopted_ids=ids, fits_dir=fits, results_path=base / "results")


def test_the_command_writes_the_reports(adopted: PKPDModel) -> None:
    assert report_command(adopted, ["--adopted"]) == adopted.report_dir
    for fit in adopted.adopted_ids:
        directory = adopted.report_dir / fit
        assert (directory / "index.html").exists()
        assert (directory / "report.pdf").read_bytes().startswith(b"%PDF")
        typst = (directory / "report.typ").read_text(encoding="utf-8")
        assert adopted.title in typst
        assert "fit DECAY" in typst


def test_the_console_script(
    adopted: PKPDModel, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr("pkpdsim.cli.resolve_model", lambda package: adopted)
    main("run_report", ["pkpdsim_toy", "--adopted", "--output_dir", str(tmp_path)])
    for fit in adopted.adopted_ids:
        assert (tmp_path / fit / "index.html").exists()


def test_a_model_without_adopted_fits(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit):
        report_command(MODEL, ["--adopted", "--output_dir", str(tmp_path)])
    assert "no adopted fits" in capsys.readouterr().out
