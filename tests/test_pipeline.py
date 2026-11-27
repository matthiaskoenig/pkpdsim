"""Tests of the pipeline of a model and of `run_sensitivity`."""

import json
import math
import os
import sys
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from pkpdsim.cli import pipeline_command, run_pipeline_main, run_sensitivity_main
from pkpdsim.pipeline import Pipeline, Step, StepResult, run_pipeline

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
    results = Path("/results").resolve().as_posix()
    assert f"local --workers 3 -o {results}/sensitivity" in text
    assert "every core: the local sensitivity per sex" in text


def test_the_pipeline_runs_its_steps_timed_and_logged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # the steps run in processes of their own, which import the test packages
    path = os.environ.get("PYTHONPATH")
    monkeypatch.setenv(
        "PYTHONPATH", os.pathsep.join([str(TESTS), *([path] if path else [])])
    )
    from pkpdsim_sexes import MODEL

    out = pipeline_command(
        MODEL, ["-c", "1", "-j", "1", "--memory", "8", "-r", str(tmp_path)]
    )
    assert out.parent == tmp_path / "pipeline"
    timings = pd.read_csv(out / "timings.tsv", sep="\t")
    assert list(timings["step"]) == ["experiments", "sensitivity"]
    assert list(timings["status"]) == ["ok", "ok"]
    assert (timings["wall_s"] > 0).all()
    assert (out / "logs" / "01_experiments.log").is_file()
    assert (tmp_path / "sensitivity" / "local").is_dir()
    record = json.loads((out / "pipeline.json").read_text(encoding="utf-8"))
    assert record["cores"] == 1
    assert (record["jobs"], record["memory_gb"]) == (1, 8)
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


# ----------------------------------------------------------------------------
# the schedule of the steps, with the script `pipeline_step.py` as the program
# ----------------------------------------------------------------------------
#: the program of the steps of the tests of the schedule
PROGRAM = [sys.executable, str(TESTS / "pipeline_step.py")]


def fake(
    name: str,
    *args: str,
    after: tuple[str, ...] = (),
    cores: int | None = None,
    memory_gb: float = 0.0,
) -> Step:
    """Get a step of the script which records its run in `{results}/<name>.json`."""
    return Step(
        name,
        "fake",
        ("--record", f"{{results}}/{name}.json", "--cores", "{cores}", *args),
        after=after,
        cores=cores,
        memory_gb=memory_gb,
    )


def run(
    tmp_path: Path,
    steps: Sequence[Step],
    cores: int = 4,
    memory_gb: float = math.inf,
    jobs: int | None = None,
    keep_going: bool = False,
) -> tuple[list[StepResult], dict[str, dict[str, Any]]]:
    """Run steps of the script, get their results and records by name."""
    results = run_pipeline(
        "package",
        "Pipeline",
        steps,
        tmp_path,
        tmp_path / "run",
        cores=cores,
        memory_gb=memory_gb,
        jobs=jobs,
        keep_going=keep_going,
        program=PROGRAM,
    )
    records = {
        path.stem: json.loads(path.read_text(encoding="utf-8"))
        for path in tmp_path.glob("*.json")
    }
    return results, records


def overlap(a: dict[str, Any], b: dict[str, Any]) -> bool:
    """Whether two steps ran at the same time, by their records."""
    return a["start"] < b["end"] and b["start"] < a["end"]


def test_a_step_declares_its_steps_after_its_cores_and_its_memory() -> None:
    step = Step("fit", "run_fit", cores=8, memory_gb=2.5)
    assert step.given_cores(20) == 8
    assert step.given_cores(4) == 4
    assert Step("fit", "run_fit").given_cores(20) == 20
    with pytest.raises(ValueError, match="at least 1"):
        Step("fit", "run_fit", cores=0)
    with pytest.raises(ValueError, match="at least 0"):
        Step("fit", "run_fit", memory_gb=-1)
    with pytest.raises(ValueError, match="after of the step 'fit' are a string"):
        Step("fit", "run_fit", after="factory")  # ty: ignore[invalid-argument-type]
    with pytest.raises(ValueError, match="no step declared before it"):
        Pipeline(
            steps=(
                Step("fit", "run_fit", after=("factory",)),
                Step("factory", "run_simulation"),
            )
        )
    with pytest.raises(ValueError, match="no step declared before it"):
        Pipeline(steps=(Step("fit", "run_fit", after=("fit",)),))


def test_the_schedule_starts_the_steps_which_fit_in_their_order() -> None:
    from pkpdsim.pipeline import _Schedule

    def results(*steps: Step) -> list[StepResult]:
        return [
            StepResult(step, [], Path(f"{step.name}.log"), step.given_cores(4))
            for step in steps
        ]

    # a step which does not fit lets the steps after it start
    x, y, z = results(
        Step("x", "c", cores=3), Step("y", "c", cores=2), Step("z", "c", cores=1)
    )
    assert _Schedule([x, y, z], cores=4).ready() == [x, z]
    # a step of the whole budget runs alone, in its order
    w, v = results(Step("w", "c"), Step("v", "c", cores=1))
    assert _Schedule([w, v], cores=4).ready() == [w]
    w.status = "running"
    assert _Schedule([w, v], cores=4).ready() == []
    # the memory and the jobs bound the steps as well
    a, b = results(
        Step("a", "c", cores=1, memory_gb=3), Step("b", "c", cores=1, memory_gb=3)
    )
    assert _Schedule([a, b], cores=4, memory_gb=4).ready() == [a]
    assert _Schedule([a, b], cores=4, jobs=1).ready() == [a]
    # a step which exceeds the budget runs when no other step runs
    (big,) = results(Step("big", "c", cores=1, memory_gb=16))
    assert _Schedule([big], cores=4, memory_gb=4).ready() == [big]


def test_steps_without_cores_run_one_after_the_other_in_their_order(
    tmp_path: Path,
) -> None:
    steps = [fake(name, "--seconds", "0.3") for name in ("c", "a", "b")]
    results, records = run(tmp_path, steps)
    assert [r.status for r in results] == ["ok", "ok", "ok"]
    assert records["c"]["end"] <= records["a"]["start"]
    assert records["a"]["end"] <= records["b"]["start"]
    # the whole budget, as `{cores}` and as the CPUs of the processes
    assert {record["cores"] for record in records.values()} == {4}
    assert {record["cpus"] for record in records.values()} == {4}
    timings = pd.read_csv(tmp_path / "run" / "timings.tsv", sep="\t")
    assert list(timings["step"]) == ["c", "a", "b"]
    assert list(timings["log"]) == ["01_c.log", "02_a.log", "03_b.log"]
    assert list(timings["cores_given"]) == [4, 4, 4]
    assert (timings["start_s"].diff().iloc[1:] > 0.3).all()


def test_steps_with_cores_run_side_by_side_within_the_budget(
    tmp_path: Path,
) -> None:
    steps = [
        fake("a", "--seconds", "1.0", cores=2),
        fake("b", "--seconds", "1.0", cores=2),
        fake("c", "--seconds", "0.2", cores=3),
    ]
    results, records = run(tmp_path, steps)
    a, b, c = (records[name] for name in "abc")
    assert overlap(a, b)
    assert c["start"] >= min(a["end"], b["end"])
    # the cores of a step, as `{cores}` and as the CPUs of its processes
    assert [record["cores"] for record in (a, b, c)] == [2, 2, 3]
    assert [record["cpus"] for record in (a, b, c)] == [2, 2, 3]
    record = json.loads((tmp_path / "run" / "pipeline.json").read_text())
    assert record["steps_wall_s"] > record["wall_s"]
    assert [step["cores"] for step in record["steps"]] == [2, 2, 3]
    assert [r.step.name for r in results] == ["a", "b", "c"]


def test_a_step_runs_after_its_steps_and_alone_with_one_job(tmp_path: Path) -> None:
    steps = [
        fake("a", "--seconds", "0.6", cores=1),
        fake("b", "--seconds", "0.2", cores=1, after=("a",)),
        fake("c", "--seconds", "0.6", cores=1),
    ]
    _, records = run(tmp_path, steps)
    assert records["b"]["start"] >= records["a"]["end"]
    assert overlap(records["a"], records["c"])
    _, records = run(tmp_path, steps, jobs=1)
    assert records["a"]["end"] <= records["c"]["start"]
    assert records["c"]["cores"] == 1
    # the memory keeps two steps apart which fit into the cores
    steps = [fake(name, "--seconds", "0.4", cores=1, memory_gb=3) for name in "ab"]
    _, records = run(tmp_path, steps, memory_gb=4)
    assert not overlap(records["a"], records["b"])


def test_a_failed_step_skips_the_steps_which_did_not_start(tmp_path: Path) -> None:
    steps = [
        fake("broken", "--exit", "3", cores=1),
        fake("long", "--seconds", "1.0", cores=1),
        fake("after", cores=1, after=("broken",)),
        fake("other", cores=2),
    ]
    results, _ = run(tmp_path, steps, cores=2)
    assert [r.status for r in results] == ["failed", "ok", "skipped", "skipped"]
    assert results[0].exit_code == 3
    timings = pd.read_csv(tmp_path / "run" / "timings.tsv", sep="\t")
    assert list(timings["status"]) == ["failed", "ok", "skipped", "skipped"]
    assert timings["start_s"].isna().tolist() == [False, False, True, True]
    # with keep going the steps after it are skipped, the others run
    results, _ = run(tmp_path, steps, cores=2, keep_going=True)
    assert [r.status for r in results] == ["failed", "ok", "skipped", "ok"]


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="the CPU time of a step is measured where its processes can be read, "
    "not on Windows",
)
def test_the_cpu_time_of_steps_side_by_side_is_their_own(tmp_path: Path) -> None:
    steps = [
        fake("busy", "--burn", "1.0", cores=1),
        fake("idle", "--seconds", "1.0", cores=1),
    ]
    results, records = run(tmp_path, steps)
    assert overlap(records["busy"], records["idle"])
    busy, idle = results
    assert busy.cpu >= 0.9
    assert idle.cpu < 0.5


def test_an_interrupt_stops_the_processes_of_the_running_steps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import pkpdsim.pipeline

    waits = 0

    def wait(seconds: float) -> None:
        nonlocal waits
        waits += 1
        time.sleep(seconds)
        if waits == 5:
            raise KeyboardInterrupt

    monkeypatch.setattr(pkpdsim.pipeline, "_wait", wait)
    monkeypatch.setattr(pkpdsim.pipeline, "GRACE", 1.0)
    steps = [
        fake("ends", "--seconds", "60", "--child", cores=1),
        fake("stays", "--seconds", "60", "--child", "--ignore-interrupt", cores=1),
        fake("waits", cores=4),
    ]
    begin = time.monotonic()
    with pytest.raises(KeyboardInterrupt):
        run(tmp_path, steps)
    assert time.monotonic() - begin < 10
    timings = pd.read_csv(tmp_path / "run" / "timings.tsv", sep="\t")
    assert list(timings["status"]) == ["interrupted", "interrupted", "skipped"]
    for name in ("ends", "stays"):
        record = json.loads((tmp_path / f"{name}.json").read_text(encoding="utf-8"))
        assert "end" not in record
        for pid in (record["pid"], record["child"]):
            assert not alive(pid)


def alive(pid: int) -> bool:
    """Whether a process runs, after a second for its parent to reap it."""
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline:
        try:
            stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
        except OSError:
            return False
        if stat[stat.rfind(")") + 2] == "Z":
            return False
        time.sleep(0.1)
    return True
