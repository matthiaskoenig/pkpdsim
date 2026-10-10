import dataclasses
from pathlib import Path

import pytest

from pkpdsim.cli import fit_command, fit_parser, report_command, with_output_dir
from pkpdsim_toy import MODEL


def test_with_output_dir() -> None:
    assert with_output_dir(["-x", "A"], Path("out")) == [
        "-x",
        "A",
        "--output_dir",
        "out",
    ]
    assert with_output_dir(["-o", "mine"], Path("out")) == ["-o", "mine"]
    assert with_output_dir(["--output_dir=mine"], Path("out")) == ["--output_dir=mine"]


def test_the_options_of_the_fit_command() -> None:
    text = fit_parser(MODEL, MODEL.load_fit_definitions()).format_help()
    for option in [
        "--subset",
        "--name",
        "--no-mapping-figures",
        "--cores",
        "--runs",
        "--seed",
        "--method",
        "--strategy",
        "--experiments",
        "--timeout",
        "--output_dir",
        "--identifiability",
        "--no-identifiability",
        "--parameter",
        "--alpha",
        "--df",
        "--max-points",
        "--no-reoptimize",
        "--no-fisher",
    ]:
        assert option in text, option
    assert "fit_toy" in text


def test_the_identifiability_default_comes_from_the_model() -> None:
    definitions = MODEL.load_fit_definitions()
    assert not fit_parser(MODEL, definitions).parse_args([]).identifiability
    model = dataclasses.replace(MODEL, fit_identifiability=True)
    assert fit_parser(model, definitions).parse_args([]).identifiability


def test_a_model_without_fits_has_no_fit_command() -> None:
    model = dataclasses.replace(MODEL, fit_definitions=None)
    with pytest.raises(ValueError, match="no fit definitions"):
        fit_command(model, [])


def test_fit_and_report(tmp_path: Path) -> None:
    runs = fit_command(
        MODEL,
        ["--runs", "1", "--output_dir", str(tmp_path / "fit"), "--no-mapping-figures"],
    )
    assert runs
    (parameters,) = (tmp_path / "fit").glob("**/parameters.json")
    report = report_command(
        MODEL,
        [
            "--subset",
            "DECAY",
            "--no-mapping-figures",
            "--output_dir",
            str(tmp_path / "report"),
            str(parameters),
        ],
    )
    assert report.is_relative_to(tmp_path / "report")


def test_fit_with_identifiability(tmp_path: Path) -> None:
    runs = fit_command(
        MODEL,
        [
            "--runs",
            "1",
            "--identifiability",
            "--max-points",
            "3",
            "--no-fisher",
            "--no-mapping-figures",
            "--output_dir",
            str(tmp_path),
        ],
    )
    assert runs


def test_the_fit_writes_into_the_results_of_the_model(tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, results_path=tmp_path)
    model.fit_main(["--runs", "1", "--no-mapping-figures"])
    assert list((tmp_path / "fit").glob("**/parameters.json"))
