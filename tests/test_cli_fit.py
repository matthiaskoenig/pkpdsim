import argparse
import json
import os
from collections.abc import Callable
from pathlib import Path

import pytest
from sbmlsim.fit.cli import FitDefinition
from sbmlsim.fit.records import adopt_fit
from sbmlsim.project import SimulatorSettings
from sbmlsim.report.render import escape_typst

from pkpdsim.cli import (
    fit_command,
    fit_parser,
    identifiability_command,
    report_command,
    with_output_dir,
)
from pkpdsim.model import PKPDModel
from pkpdsim_toy import MODEL

type AsModel = Callable[[PKPDModel], PKPDModel]

#: the files of the records and the report of a fit, which the old command wrote
FIT_FILES = {
    "datapoints.tsv",
    "fisher.json",
    "fisher.tsv",
    "identifiability.json",
    "identifiability.tsv",
    "index.html",
    "metrics.tsv",
    "metrics_mappings.tsv",
    "metrics_observables.tsv",
    "optimization_result.json",
    "optimization_result.tsv",
    "parameters.json",
    "plots/bland_altman.png",
    "plots/bland_altman.svg",
    "plots/cost_bar.svg",
    "plots/goodness_of_fit.png",
    "plots/goodness_of_fit.svg",
    "plots/profile_k.svg",
    "plots/profiles.png",
    "plots/profiles.svg",
    "plots/residual_boxplot.svg",
    "plots/traces.png",
    "plots/traces.svg",
    "references.bib",
    "references.csl",
    "report.pdf",
    "report.tex",
    "report.txt",
    "report.typ",
    "report_standalone.tex",
}

#: the records of a run of a kind of the run layer
RUN_FILES = {"run.json", "spec.json", "summary.json", "progress.jsonl"}


def _files(directory: Path) -> set[str]:
    """Get the files below a directory, relative to it."""
    return {
        p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()
    }


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
        "--start",
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
    assert "usage: run_fit pkpdsim_toy " in text


def test_the_identifiability_default_comes_from_the_model() -> None:
    definitions = MODEL.load_fit_definitions()
    assert not fit_parser(MODEL, definitions).parse_args([]).identifiability
    model = MODEL.replace(fit_identifiability=True)
    assert fit_parser(model, definitions).parse_args([]).identifiability


def test_a_model_without_fits_has_no_fit_command() -> None:
    model = MODEL.replace(fits={})
    with pytest.raises(ValueError, match="no fit definitions"):
        fit_command(model, [])


def test_fit_and_report(tmp_path: Path) -> None:
    runs = fit_command(
        MODEL,
        ["--runs", "1", "--output_dir", str(tmp_path / "fit"), "--no-mapping-figures"],
    )
    ((opid, records),) = runs.items()
    assert records == tmp_path / "fit" / opid
    (parameters,) = (tmp_path / "fit").glob("**/parameters.json")
    assert parameters == records / "parameters.json"
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
    assert report == tmp_path / "report" / "report"
    assert (report / "index.html").is_file()
    assert _files(report) >= RUN_FILES


def test_the_fit_writes_the_files_of_the_old_command(tmp_path: Path) -> None:
    """The end-to-end run of `run_fit`: the records and the report of the fit.

    The kind `fit` writes them into the folder of the fit id of its output,
    which the command moves up into the output, `<output_dir>/<fit id>`.
    """
    runs = fit_command(
        MODEL,
        [
            "--runs",
            "1",
            "--identifiability",
            "--max-points",
            "3",
            "--no-mapping-figures",
            "--output_dir",
            str(tmp_path),
        ],
    )
    ((opid, records),) = runs.items()
    assert records == tmp_path / opid
    files = _files(records)
    assert files == FIT_FILES | RUN_FILES | {"fits.tsv", f"runs/{opid}_run_0.json"}
    spec = json.loads((records / "spec.json").read_text(encoding="utf-8"))
    assert spec["project"] == "pkpdsim_toy:MODEL"
    assert spec["opid"] == opid
    assert (spec["repeats"], spec["seed"], spec["algorithm"]) == (1, 1234, "LSQ")
    assert spec["record_adopted"] and spec["report"] and spec["identifiability"]


def test_the_identifiability_writes_the_report_of_the_old_command(
    tmp_path: Path,
) -> None:
    """The end-to-end run of `run_identifiability`: the report of the analysis."""
    parameters = fitted_parameters(tmp_path / "fit")
    report = identifiability_command(
        MODEL,
        [
            "--max-points",
            "3",
            "--no-mapping-figures",
            "--output_dir",
            str(tmp_path / "identifiability"),
            str(parameters),
        ],
    )
    assert report == tmp_path / "identifiability" / "identifiability"
    files = _files(report)
    assert files >= FIT_FILES
    # the runs of the kinds `identifiability` and `report`
    assert {f"report/{name}" for name in RUN_FILES} <= files
    assert {f"analysis/{name}" for name in RUN_FILES} <= files
    assert "analysis/identifiability.json" in files


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
    content = json.loads((report / "identifiability.json").read_text(encoding="utf-8"))
    assert content["opid"] == parameters.parent.name
    assert content["parameter_set"]["values"].keys() == {"k"}


def test_the_reports_write_into_the_results_of_the_model(
    tmp_path: Path, as_model: AsModel
) -> None:
    parameters = fitted_parameters(tmp_path / "fit")
    model = as_model(MODEL.replace(results_path=tmp_path))
    report_command(
        model, ["--subset", "DECAY", "--no-mapping-figures", str(parameters)]
    )
    assert list((tmp_path / "report").glob("**/index.html"))
    identifiability_command(
        model,
        ["--subset", "DECAY", "--max-points", "3", "--no-fisher", str(parameters)],
    )
    assert (tmp_path / "identifiability" / "identifiability" / "index.html").is_file()


def test_a_model_with_an_empty_mapping_of_fit_definitions_fails() -> None:
    with pytest.raises(ValueError, match="'toy' has no fit definitions"):
        fit_parser(MODEL, {})


def test_the_fit_writes_into_the_results_of_the_model(
    tmp_path: Path, as_model: AsModel
) -> None:
    model = as_model(MODEL.replace(results_path=tmp_path))
    fit_command(model, ["--runs", "1", "--no-mapping-figures"])
    assert list((tmp_path / "fit").glob("**/parameters.json"))


def test_the_help_states_the_default_of_the_identifiability() -> None:
    definitions = MODEL.load_fit_definitions()
    off = fit_parser(MODEL, definitions).format_help()
    off = " ".join(off.split())
    assert "default: False" in off
    assert "--identifiability analyses the identifiability" in off
    model = MODEL.replace(fit_identifiability=True)
    on = fit_parser(model, definitions).format_help()
    on = " ".join(on.split())
    assert "default: True" in on
    assert "--no-identifiability skips it" in on


def _report_dir(output: Path) -> Path:
    """Get the directory of the one report below an output directory."""
    (index,) = output.glob("**/index.html")
    return index.parent


def test_the_fit_report_has_the_model_context(tmp_path: Path) -> None:
    runs = fit_command(MODEL, ["--runs", "1", "--output_dir", str(tmp_path)])
    (opid,) = runs
    report = _report_dir(tmp_path)
    html = (report / "index.html").read_text(encoding="utf-8")
    assert "<math" in html and 'data-tab="report"' in html
    typst = (report / "report.typ").read_text(encoding="utf-8")
    assert MODEL.title in typst
    assert (report / "report.pdf").exists()
    (result,) = tmp_path.glob("**/optimization_result.json")
    arguments = json.loads(result.read_text(encoding="utf-8"))["provenance"]
    assert arguments["fixed_values"] == {}
    assert arguments["arguments"]["fixed_fits"] == []
    assert opid.startswith("DECAY_")
    # the command is the spec of the run of the kind `fit`
    spec = json.loads((report / "spec.json").read_text(encoding="utf-8"))
    assert spec["subset"] == "DECAY"


def _adopted_other(tmp_path: Path) -> Path:
    """Get a hand-made adopted fit of another subset, in `tmp_path / "fits"`.

    The record of a toy fit, copied under the id of a fit of the subset OTHER.
    """
    first = tmp_path / "first"
    fit_command(
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
            str(first),
        ],
    )
    (fit,) = [p for p in first.iterdir() if p.is_dir()]
    fits = tmp_path / "fits"
    adopt_fit(fit, fits)
    other = fits / "OTHER_20260101_000000__abcde"
    (fits / fit.name).rename(other)
    return other


def test_the_fit_names_the_fixed_adopted_fits(
    tmp_path: Path, as_model: AsModel
) -> None:
    other = _adopted_other(tmp_path)
    model = as_model(MODEL.replace(adopted_ids=[other.name], fits_dir=other.parent))
    output = tmp_path / "second"
    fit_command(
        model, ["--runs", "1", "--no-mapping-figures", "--output_dir", str(output)]
    )
    (result,) = output.glob("**/optimization_result.json")
    provenance = json.loads(result.read_text(encoding="utf-8"))["provenance"]
    assert provenance["arguments"]["fixed_fits"] == [other.name]
    # the fit optimizes `k` itself, the value of the adopted fit is not fixed
    assert provenance["fixed_values"] == {}
    typst = _report_dir(output).joinpath("report.typ").read_text(encoding="utf-8")
    # the Typst source escapes the underscores of the id
    assert (
        escape_typst(
            f"The parameters of the adopted fit {other.name} were fixed at their "
            f"fitted values."
        )
        in typst
    )


def _rename(value: object, old: str, new: str) -> object:
    """Rename a string and a key of a JSON document."""
    if isinstance(value, dict):
        return {
            (new if k == old else k): _rename(v, old, new) for k, v in value.items()
        }
    if isinstance(value, list):
        return [_rename(v, old, new) for v in value]
    return new if value == old else value


def test_the_fixed_values_leave_out_the_parameters_of_the_fit(
    tmp_path: Path, as_model: AsModel
) -> None:
    other = _adopted_other(tmp_path)
    # the adopted fit of the other subset fitted a parameter `k_other`
    for name in ("parameters.json", "optimization_result.json", "identifiability.json"):
        path = other / name
        data = json.loads(path.read_text(encoding="utf-8"))
        path.write_text(json.dumps(_rename(data, "k", "k_other")), encoding="utf-8")
    model = as_model(MODEL.replace(adopted_ids=[other.name], fits_dir=other.parent))
    value = model.fit_records().records[other.name].parameters["k_other"].value

    output = tmp_path / "second"
    fit_command(
        model, ["--runs", "1", "--no-mapping-figures", "--output_dir", str(output)]
    )
    (result,) = output.glob("**/optimization_result.json")
    provenance = json.loads(result.read_text(encoding="utf-8"))["provenance"]
    assert provenance["fixed_values"] == {"k_other": value}


def test_the_report_of_stored_parameters_has_the_context_of_its_record(
    tmp_path: Path, as_model: AsModel
) -> None:
    other = _adopted_other(tmp_path)
    model = as_model(MODEL.replace(adopted_ids=[other.name], fits_dir=other.parent))
    output = tmp_path / "second"
    fit_command(
        model, ["--runs", "1", "--no-mapping-figures", "--output_dir", str(output)]
    )
    parameters = _report_dir(output) / "parameters.json"
    # the adopted fits of the model changed since the fit ran
    later = as_model(MODEL.replace(adopted_ids=[], results_path=tmp_path / "results"))
    report = report_command(
        later, ["--subset", "DECAY", "--no-mapping-figures", str(parameters)]
    )
    typst = (report / "report.typ").read_text(encoding="utf-8")
    assert f"{MODEL.title}: fit DECAY" in typst
    assert escape_typst(f"adopted fit {other.name} were fixed") in typst


#: a fit of one run without the figures of the mappings
ONE_RUN = ("--runs", "1", "--no-mapping-figures")


def _optimization(directory: Path) -> dict:
    """Get the optimization result of the one fit in a directory."""
    (result,) = directory.glob("*/optimization_result.json")
    return json.loads(result.read_text(encoding="utf-8"))


def test_a_fit_starts_from_a_parameters_file(tmp_path: Path) -> None:
    """Least squares starts every run from the given point, not from a sample."""
    start = tmp_path / "start.json"
    start.write_text(json.dumps({"k": 0.3}), encoding="utf-8")
    output = tmp_path / "fit"
    fit_command(
        MODEL,
        [*ONE_RUN, "--start", str(start), "--output_dir", str(output)],
    )
    result = _optimization(output)
    assert result["provenance"]["start_values"] == {"k": 0.3}
    (run,) = result["fits"]
    assert run["x0"] == [0.3]
    assert result["provenance"]["arguments"]["start"] == str(start)
    # the definition of the model keeps its own start value
    (parameter,) = MODEL.load_fit_definitions()["DECAY"].parameters
    assert parameter.start_value == 0.5


def test_a_fit_starts_from_another_fit(toy_fit: Path) -> None:
    """A fit id is the best parameters of a finished fit in the output directory."""
    best = json.loads((toy_fit / "parameters.json").read_text(encoding="utf-8"))
    k = best["sets"][0]["values"]["k"]
    runs = fit_command(
        MODEL,
        [*ONE_RUN, "--start", toy_fit.name, "--output_dir", str(toy_fit.parent)],
    )
    (opid,) = runs
    result = json.loads(
        (toy_fit.parent / opid / "optimization_result.json").read_text("utf-8")
    )
    assert result["provenance"]["start_values"] == {"k": pytest.approx(k, rel=1e-12)}


def test_a_fit_starts_from_an_adopted_fit(
    fits_dir: Path, tmp_path: Path, as_model: AsModel
) -> None:
    """An adopted fit of the model is found by its id."""
    (record,) = fits_dir.iterdir()
    model = as_model(MODEL.replace(fits_dir=fits_dir))
    output = tmp_path / "fit"
    fit_command(
        model,
        [*ONE_RUN, "--start", record.name, "--output_dir", str(output)],
    )
    best = json.loads((record / "parameters.json").read_text(encoding="utf-8"))
    assert _optimization(output)["provenance"]["start_values"] == {
        "k": pytest.approx(best["sets"][0]["values"]["k"], rel=1e-12)
    }


@pytest.mark.parametrize(
    ("values", "extra", "message"),
    [
        ({"x": 0.3}, [], "not fitted"),
        ({"k": 2.0}, [], "outside of its bounds"),
        ({"k": "a"}, [], "not a number"),
        ({"k": 0.3}, ["--runs", "2"], "identical"),
        ({"k": 0.3}, ["--method", "DE"], "least squares"),
    ],
)
def test_a_start_which_cannot_be_used_is_refused(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    values: dict,
    extra: list[str],
    message: str,
) -> None:
    start = tmp_path / "start.json"
    start.write_text(json.dumps(values), encoding="utf-8")
    with pytest.raises(SystemExit) as exc:
        fit_command(MODEL, ["--runs", "1", *extra, "--start", str(start)])
    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_an_unknown_start_is_refused(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit):
        fit_command(
            MODEL, ["--runs", "1", "--start", "DECAY_nofit", "-o", str(tmp_path)]
        )
    assert "neither a json file nor a fit" in capsys.readouterr().err


def _level(definition: object) -> str:
    """Get the level of the problem of a fit definition."""
    assert isinstance(definition, FitDefinition)
    return str(definition.problem(opid="level").jit_optimization)


class _Stop(Exception):
    """Stops a fit once its problem is known."""


def test_the_fits_compile_at_the_level_of_the_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, as_model: AsModel
) -> None:
    """The problem of a fit has `jit_optimization` of the model, not the environment.

    The problem carries its level into the workers of a pool
    (`tests/test_fit.py::test_a_worker_loads_the_models_at_the_level_of_the_model`).
    A model which names no level is created with the level of pkpdsim, unless
    the variable of sbmlsim is set in the process which creates it: then the
    model names none and sbmlsim resolves the variable when it loads it.
    """
    import sbmlsim.fit.cli

    variable = "SBMLSIM_JIT_OPTIMIZATION"
    monkeypatch.setenv(variable, "O3")
    seen: list[str] = []

    def recording_run_fit(definition: object, **kwargs: object) -> dict:
        seen.append(_level(definition))
        assert os.environ[variable] == "O3"
        raise _Stop

    monkeypatch.setattr(sbmlsim.fit.cli, "run_fit", recording_run_fit)
    simulator = SimulatorSettings(
        absolute_tolerance=1e-10, relative_tolerance=1e-10, jit_optimization="O0"
    )
    model = as_model(MODEL.replace(simulator=simulator))
    with pytest.raises(_Stop):
        fit_command(model, ["--runs", "1", "--output_dir", str(tmp_path / "a")])
    assert seen == ["O0"]
    # a model created with the variable set names no level, so the level of the
    # environment is the one of its fits
    plain = as_model(MODEL.replace(simulator=SimulatorSettings()))
    assert plain.simulator.jit_optimization is None
    with pytest.raises(_Stop):
        fit_command(plain, ["--runs", "1", "--output_dir", str(tmp_path / "b")])
    assert seen == ["O0", "O3"]


def test_the_reports_compile_at_the_level_of_the_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`run_report` and `run_identifiability` build the problem at the level too."""
    import pkpdsim.cli

    seen: list[str] = []

    def recording(definitions: dict, **kwargs: object) -> Path:
        seen.extend(_level(d) for d in definitions.values())
        return tmp_path

    monkeypatch.setattr(pkpdsim.cli, "report_cli", recording)
    monkeypatch.setattr(pkpdsim.cli, "identifiability_cli", recording)
    simulator = SimulatorSettings(
        absolute_tolerance=1e-10, relative_tolerance=1e-10, jit_optimization="O1"
    )
    model = MODEL.replace(simulator=simulator)
    report_command(model, ["parameters.json"])
    identifiability_command(model, ["parameters.json"])
    assert seen == ["O1", "O1"]
