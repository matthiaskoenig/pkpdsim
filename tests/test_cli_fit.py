import argparse
import dataclasses
import json
from pathlib import Path

import pytest

from pkpdsim.cli import (
    fit_command,
    fit_parser,
    identifiability_command,
    report_command,
    with_output_dir,
)
from pkpdsim_toy import MODEL


def test_with_output_dir() -> None:
    assert with_output_dir(["-x", "A"], Path("out")) == [
        "--output_dir",
        "out",
        "-x",
        "A",
    ]


@pytest.mark.parametrize("option", ["--out", "--output", "--output_dir", "-o"])
def test_an_output_dir_of_the_arguments_wins(option: str) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", default="unused")
    args = parser.parse_args(with_output_dir([option, "mine"], Path("out")))
    assert args.output_dir == "mine"
    assert parser.parse_args(with_output_dir([], Path("out"))).output_dir == "out"


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


def fitted_parameters(directory: Path, *extra: str) -> Path:
    """Fit the toy model once and get the stored parameters."""
    fit_command(
        MODEL,
        ["--runs", "1", "--no-mapping-figures", *extra, "--output_dir", str(directory)],
    )
    (parameters,) = directory.glob("**/parameters.json")
    return parameters


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
    (result,) = tmp_path.glob("**/identifiability.json")
    assert "k" in result.read_text(encoding="utf-8")
    assert list(tmp_path.glob("**/profile_k.svg"))


def test_a_fit_without_identifiability_has_no_analysis(tmp_path: Path) -> None:
    fitted_parameters(tmp_path)
    assert not list(tmp_path.glob("**/identifiability.json"))


def test_the_identifiability_of_stored_parameters(tmp_path: Path) -> None:
    parameters = fitted_parameters(tmp_path / "fit")
    report = identifiability_command(
        MODEL,
        [
            "--subset",
            "DECAY",
            "--max-points",
            "3",
            "--no-fisher",
            "--output_dir",
            str(tmp_path / "identifiability"),
            str(parameters),
        ],
    )
    assert report.is_relative_to(tmp_path / "identifiability")
    (result,) = (tmp_path / "identifiability").glob("**/identifiability.json")
    content = json.loads(result.read_text(encoding="utf-8"))
    assert content["opid"] == "DECAY"
    assert content["parameter_set"]["values"].keys() == {"k"}


def test_the_commands_of_the_model_object(tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, results_path=tmp_path)
    parameters = fitted_parameters(tmp_path / "fit")
    model.report_main(["--subset", "DECAY", "--no-mapping-figures", str(parameters)])
    assert list((tmp_path / "report").glob("**/index.html"))
    model.identifiability_main(
        ["--subset", "DECAY", "--max-points", "3", "--no-fisher", str(parameters)]
    )
    assert list((tmp_path / "identifiability").glob("**/identifiability.json"))


def test_a_model_with_an_empty_mapping_of_fit_definitions_fails() -> None:
    with pytest.raises(ValueError, match="'toy' has no fit definitions"):
        fit_parser(MODEL, {})


def test_the_fit_writes_into_the_results_of_the_model(tmp_path: Path) -> None:
    model = dataclasses.replace(MODEL, results_path=tmp_path)
    model.fit_main(["--runs", "1", "--no-mapping-figures"])
    assert list((tmp_path / "fit").glob("**/parameters.json"))


def test_the_help_states_the_default_of_the_identifiability() -> None:
    definitions = MODEL.load_fit_definitions()
    off = fit_parser(MODEL, definitions).format_help()
    off = " ".join(off.split())
    assert "default: False" in off
    assert "--identifiability analyses the identifiability" in off
    model = dataclasses.replace(MODEL, fit_identifiability=True)
    on = fit_parser(model, definitions).format_help()
    on = " ".join(on.split())
    assert "default: True" in on
    assert "--no-identifiability skips it" in on
