"""Tests of the pipeline of a model and of `run_sensitivity`."""

import json
import os
from pathlib import Path

import pandas as pd
import pytest

from pkpdsim.cli import pipeline_command, run_pipeline_main, run_sensitivity_main
from pkpdsim.pipeline import Pipeline, Step

TESTS = Path(__file__).parent


def test_a_step_is_a_word_with_known_placeholders() -> None:
    step = Step("fit", "run_fit", ("--cores", "{cores}", "-o", "{results}/fit"))
    assert step.arguments({"cores": 4, "results": Path("/r")}) == [
        "--cores",
        "4",
        "-o",
        "/r/fit",
    ]
    with pytest.raises(ValueError, match="no word"):
        Step("a fit", "run_fit")
    with pytest.raises(ValueError, match=r"placeholder '\{cpus\}'"):
        Step("fit", "run_fit", ("--cores", "{cpus}"))
    with pytest.raises(ValueError, match="a string"):
        Step("fit", "run_fit", "--cores")  # ty: ignore[invalid-argument-type]


def test_a_pipeline_selects_its_steps_in_their_order() -> None:
    pipeline = Pipeline(
        steps=tuple(Step(name, "run_simulation") for name in ("a", "b", "c", "d"))
    )
    assert [s.name for s in pipeline.select()] == ["a", "b", "c", "d"]
    assert [s.name for s in pipeline.select(only=["c", "a"])] == ["a", "c"]
    assert [s.name for s in pipeline.select(skip=["b"])] == ["a", "c", "d"]
    assert [s.name for s in pipeline.select(start="c")] == ["c", "d"]
    assert [s.name for s in pipeline.select(start="b", skip=["d"])] == ["b", "c"]
    with pytest.raises(ValueError, match="'x' is no step"):
        pipeline.select(only=["x"])
    with pytest.raises(ValueError, match="two steps"):
        Pipeline(steps=(Step("a", "run_fit"), Step("a", "run_fit")))
    with pytest.raises(ValueError, match="at least one step"):
        Pipeline(steps=())


def test_the_pipeline_of_a_model_is_resolved() -> None:
    from pkpdsim_sexes import MODEL
    from pkpdsim_sexes.pipeline import PIPELINE
    from pkpdsim_toy import MODEL as TOY

    assert MODEL.load_pipeline() is PIPELINE
    with pytest.raises(ValueError, match="declares no pipeline"):
        TOY.load_pipeline()
    with pytest.raises(TypeError, match="not a Pipeline"):
        MODEL.replace(pipeline=".pipeline:Step").load_pipeline()


def test_run_sensitivity_runs_the_analyses_into_the_results(tmp_path: Path) -> None:
    out = tmp_path / "sensitivity"
    run_sensitivity_main(["pkpdsim_sexes", "local", "--workers", "1", "-o", str(out)])
    for k, group in enumerate(("men", "women")):
        assert (out / "local" / f"local_d0.01_nvar3_{k:02d}_{group}.tsv").is_file()
    assert (out / "parameters.tsv").is_file()
    # the records of a run of the kind
    assert (out / "run.json").is_file()


def test_run_sensitivity_refuses_an_unknown_analysis(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        run_sensitivity_main(["pkpdsim_sexes", "sobol"])
    assert exit_info.value.code == 2
    assert "sobol is no analysis of 'decay'" in capsys.readouterr().err


def test_run_sensitivity_of_a_model_without_one_is_a_usage_error(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as exit_info:
        run_sensitivity_main(["pkpdsim_toy"])
    assert exit_info.value.code == 2
    assert "declares no sensitivity analysis" in capsys.readouterr().err


def test_the_pipeline_lists_its_steps(capsys: pytest.CaptureFixture[str]) -> None:
    run_pipeline_main(["pkpdsim_sexes", "--list", "-c", "3", "-r", "/results"])
    text = capsys.readouterr().out
    assert "experiments: run_simulation pkpdsim_sexes --action list_experiments" in text
    assert "local --workers 3 -o /results/sensitivity" in text


def test_the_pipeline_runs_its_steps_timed_and_logged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # the steps run in processes of their own, which import the test packages
    path = os.environ.get("PYTHONPATH")
    monkeypatch.setenv(
        "PYTHONPATH", os.pathsep.join([str(TESTS), *([path] if path else [])])
    )
    from pkpdsim_sexes import MODEL

    out = pipeline_command(MODEL, ["-c", "1", "-r", str(tmp_path)])
    assert out.parent == tmp_path / "pipeline"
    timings = pd.read_csv(out / "timings.tsv", sep="\t")
    assert list(timings["step"]) == ["experiments", "sensitivity"]
    assert list(timings["status"]) == ["ok", "ok"]
    assert (timings["wall_s"] > 0).all()
    assert (out / "logs" / "01_experiments.log").is_file()
    assert (tmp_path / "sensitivity" / "local").is_dir()
    record = json.loads((out / "pipeline.json").read_text(encoding="utf-8"))
    assert record["cores"] == 1
    assert record["versions"]["pkpdsim"]


def test_a_failed_step_skips_the_steps_after_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from pkpdsim_sexes import MODEL

    path = os.environ.get("PYTHONPATH")
    monkeypatch.setenv(
        "PYTHONPATH", os.pathsep.join([str(TESTS), *([path] if path else [])])
    )
    pipeline = Pipeline(
        steps=(
            Step("broken", "run_simulation", ("--action", "nothing")),
            Step("experiments", "run_simulation", ("--action", "list_experiments")),
        )
    )
    model = MODEL.replace(pipeline=pipeline)
    with pytest.raises(SystemExit) as exit_info:
        pipeline_command(model, ["-r", str(tmp_path)])
    assert exit_info.value.code == 1
    (out,) = (tmp_path / "pipeline").iterdir()
    timings = pd.read_csv(out / "timings.tsv", sep="\t")
    assert list(timings["status"]) == ["failed", "skipped"]
    assert timings["exit_code"].iloc[0] == 2
    assert "The steps broken failed" in capsys.readouterr().out


def test_run_sensitivity_replaces_the_results_of_a_sensitivity_script(
    tmp_path: Path,
) -> None:
    out = tmp_path / "sensitivity"
    (out / "local").mkdir(parents=True)
    (out / "parameters.tsv").write_text("uid\tvalue\n", encoding="utf-8")
    run_sensitivity_main(["pkpdsim_sexes", "local", "--workers", "1", "-o", str(out)])
    assert (out / "run.json").is_file()
    unknown = tmp_path / "other"
    unknown.mkdir()
    (unknown / "notes.txt").write_text("mine", encoding="utf-8")
    with pytest.raises(SystemExit) as exit_info:
        run_sensitivity_main(
            ["pkpdsim_sexes", "local", "--workers", "1", "-o", str(unknown)]
        )
    assert exit_info.value.code == 1
    assert (unknown / "notes.txt").is_file()


def test_the_progress_of_a_step_is_its_last_line_without_level_and_location(
    tmp_path: Path,
) -> None:
    from pkpdsim.pipeline import _last_line

    log = tmp_path / "step.log"
    log.write_text(
        "INFO     fit optimizing 20/48 (42%) elapsed 0:01:46  reporters.py:260\n\n",
        encoding="utf-8",
    )
    line, offset = _last_line(log, 0)
    assert line == "fit optimizing 20/48 (42%) elapsed 0:01:46"
    assert _last_line(log, offset) == ("", offset)
