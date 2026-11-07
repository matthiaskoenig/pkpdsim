"""The commands replace a result only when its run succeeded and never touch their inputs.

A kind of the run layer runs into a staging directory next to the result; the
earlier result is replaced only after the run succeeded, and an output which is
an input of the command, contains one or lies inside one is refused before
anything happens.
"""

import json
import math
from collections.abc import Callable
from pathlib import Path

import pytest
from sbmlsim.cgroup import available_cpus
from sbmlsim.experiment.runner import EXPERIMENTS_PER_WORKER
from sbmlsim.fit.records import adopt_fit
from sbmlsim.run import ConfigurationError

import pkpdsim.cli
from pkpdsim.cli import (
    fit_command,
    identifiability_command,
    report_command,
    run_command,
    run_workers,
)
from pkpdsim.model import ExtraAction, PKPDModel
from pkpdsim_toy import MODEL

type AsModel = Callable[[PKPDModel], PKPDModel]

ONE_RUN = ["--runs", "1", "--no-mapping-figures"]


def _snapshot(directory: Path) -> dict[str, bytes]:
    """Get the files below a directory with their contents."""
    return {
        p.relative_to(directory).as_posix(): p.read_bytes()
        for p in directory.rglob("*")
        if p.is_file()
    }


@pytest.fixture(scope="module")
def fit(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Get the records of a fit of the toy model with its identifiability."""
    out = tmp_path_factory.mktemp("fits")
    runs = fit_command(
        MODEL,
        [*ONE_RUN, "--identifiability", "--max-points", "3", "--output_dir", str(out)],
    )
    (records,) = runs.values()
    return records


@pytest.mark.parametrize("command", ["report", "identifiability"])
def test_an_output_in_its_own_input_is_refused(
    fit: Path, command: str, capsys: pytest.CaptureFixture[str]
) -> None:
    """`-o <fits> --name <fit id>` would replace the fit it reads."""
    before = _snapshot(fit)
    function = report_command if command == "report" else identifiability_command
    for name, output in [
        (fit.name, fit.parent),  # the fit itself
        ("plots", fit),  # inside the fit
        (fit.parent.name, fit.parent.parent),  # the folder of the fit
    ]:
        with pytest.raises(SystemExit) as info:
            function(
                MODEL,
                [
                    "--no-mapping-figures",
                    *("--name", name, "-o", str(output)),
                    str(fit / "parameters.json"),
                ],
            )
        assert info.value.code == 2
        assert "input" in capsys.readouterr().err
    assert _snapshot(fit) == before


def test_the_adopted_reports_are_not_written_into_the_adopted_fits(
    fit: Path, tmp_path: Path, as_model: AsModel, capsys: pytest.CaptureFixture[str]
) -> None:
    fits = tmp_path / "fits"
    adopt_fit(fit, fits)
    model = as_model(MODEL.replace(adopted_ids=[fit.name], fits_dir=fits))
    before = _snapshot(fits)
    with pytest.raises(SystemExit) as info:
        report_command(model, ["--adopted", "--no-mapping-figures", "-o", str(fits)])
    assert info.value.code == 2
    assert _snapshot(fits) == before


def test_a_failed_plan_keeps_the_earlier_result(fit: Path, tmp_path: Path) -> None:
    arguments = [
        "--no-mapping-figures",
        "--no-fisher",
        *("--max-points", "3", "-o", str(tmp_path)),
        str(fit / "parameters.json"),
    ]
    report = identifiability_command(MODEL, arguments)
    before = _snapshot(report)
    with pytest.raises(SystemExit) as info:
        identifiability_command(MODEL, ["-p", "unknown", *arguments])
    assert info.value.code == 2
    assert _snapshot(report) == before
    assert [p.name for p in tmp_path.iterdir()] == [report.name]


def test_a_failed_simulation_keeps_the_earlier_result(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    arguments = ["-a", "simulate", "-e", "studies", "-w", "1", "-r", str(tmp_path)]
    run_command(MODEL, arguments)
    out = tmp_path / "simulation" / "studies"
    before = _snapshot(out)

    def failing(*args: object, **kwargs: object) -> dict:
        raise RuntimeError("the run fails")

    monkeypatch.setattr(pkpdsim.cli, "execute", failing)
    with pytest.raises(RuntimeError, match="the run fails"):
        run_command(MODEL, arguments)
    assert _snapshot(out) == before
    assert [p.name for p in out.parent.iterdir()] == ["studies"]


def test_a_failed_move_is_a_message_and_keeps_the_earlier_result(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    target = tmp_path / "result"
    target.mkdir()
    (target / "earlier.txt").write_text("earlier", encoding="utf-8")

    def run(staging: Path) -> None:
        (staging / "inner").mkdir(parents=True)
        (staging / "inner" / "index.html").write_text("new", encoding="utf-8")
        # the records of the run of the same name as a file of the result
        (staging / "index.html").write_text("clash", encoding="utf-8")

    from pkpdsim.cli import _staged

    with pytest.raises(SystemExit) as info:
        _staged(target, run, inner="inner")
    assert info.value.code == 1
    assert "index.html" in capsys.readouterr().out
    assert [p.name for p in tmp_path.iterdir()] == ["result"]
    assert (target / "earlier.txt").read_text(encoding="utf-8") == "earlier"


@pytest.mark.parametrize(
    "fit_id", [".", "../fits", "DECAY_nofit", "DECAY_20260101_000000__00000"]
)
def test_adopt_checks_its_fit_first(
    fit: Path, fit_id: str, as_model: AsModel, tmp_path: Path
) -> None:
    model = as_model(MODEL.replace(fits_dir=tmp_path / "fits"))
    before = _snapshot(fit.parent)
    with pytest.raises(SystemExit) as info:
        fit_command(model, ["--adopt", fit_id, "--output_dir", str(fit.parent)])
    assert info.value.code in (1, 2)
    assert _snapshot(fit.parent) == before
    assert not (tmp_path / "fits").exists()


def test_an_adoption_writes_nothing_into_the_fit(
    fit: Path, as_model: AsModel, tmp_path: Path
) -> None:
    model = as_model(MODEL.replace(fits_dir=tmp_path / "fits"))
    before = _snapshot(fit)
    fit_command(model, ["--adopt", fit.name, "--output_dir", str(fit.parent)])
    assert (tmp_path / "fits" / fit.name / "parameters.json").is_file()
    assert _snapshot(fit) == before


def test_the_default_workers_are_the_old_default() -> None:
    """`model`: a worker per four experiments, at most one per CPU, serially below."""
    for n in (1, 4, 5, 9, 1000):
        expected = max(1, min(available_cpus(), math.ceil(n / EXPERIMENTS_PER_WORKER)))
        assert run_workers("model", n) == expected
    assert run_workers("model", 2) == 1
    assert run_workers("auto", 9) is None
    assert run_workers("3", 9) == 3


def test_auto_leaves_the_workers_to_sbmlsim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The kind resolves its workers, so `auto` runs the experiments of sbmlsim."""
    seen: list[object] = []
    original = pkpdsim.cli.Project.run_experiments

    def recording(project: object, *args: object, **kwargs: object) -> object:
        seen.append(kwargs.get("n_workers"))
        return original(project, *args, **kwargs)  # ty: ignore[invalid-argument-type]

    monkeypatch.setattr(pkpdsim.cli.Project, "run_experiments", recording)
    run_command(
        MODEL, ["-a", "simulate", "-e", "studies", "-w", "auto", "-r", str(tmp_path)]
    )
    assert seen == [None]
    out = tmp_path / "simulation" / "studies"
    assert (out / "index.html").is_file()
    assert list((out / "_figures").iterdir())


def test_an_error_of_a_run_which_is_no_usage_error_exits_with_1(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def failing(*args: object, **kwargs: object) -> dict:
        raise ConfigurationError("the records are not read")

    monkeypatch.setattr(pkpdsim.cli, "execute", failing)
    with pytest.raises(SystemExit) as info:
        run_command(
            MODEL, ["-a", "simulate", "-e", "studies", "-w", "1", "-r", str(tmp_path)]
        )
    assert info.value.code == 1
    assert "the records are not read" in capsys.readouterr().out


GREET = ExtraAction(name="greet", help="write a greeting", function=".extras:greet")
EXTRA = MODEL.replace(extra_actions=(GREET,), extra_arguments=".extras:add_arguments")


def test_an_extra_action_of_the_kind_checks_and_converts_its_arguments(
    tmp_path: Path,
) -> None:
    function = EXTRA.load_action("greet")
    function(EXTRA, tmp_path, workers=1, repeat="2")
    text = (tmp_path / "greeting.txt").read_text(encoding="utf-8")
    assert text == "hello toy hello toy"
    with pytest.raises(ValueError, match="greting"):
        function(EXTRA, tmp_path, workers=1, greting="hi")


def test_a_failed_fit_keeps_its_records_where_they_are(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A fit whose report fails has its records in the folder of its id."""
    import sbmlsim.fit.cli

    def failing(*args: object, **kwargs: object) -> object:
        raise RuntimeError("the report fails")

    monkeypatch.setattr(sbmlsim.fit.cli, "report_from_records", failing)
    with pytest.raises(RuntimeError, match="the report fails"):
        fit_command(MODEL, [*ONE_RUN, "--output_dir", str(tmp_path)])
    (out,) = tmp_path.iterdir()
    assert (out / "parameters.json").is_file()
    assert not (out / out.name).exists()
    run = json.loads((out / "run.json").read_text(encoding="utf-8"))
    assert run["status"] == "failed"


def test_the_message_of_a_replacement_shows_the_path(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    results = tmp_path / "a[b]"
    arguments = ["-a", "simulate", "-e", "studies", "-w", "1", "-r", str(results)]
    run_command(MODEL, arguments)
    capsys.readouterr()
    run_command(MODEL, arguments)
    assert "a[b]/simulation/studies" in "".join(capsys.readouterr().out.split())


def test_the_stagings_of_a_killed_run_are_removed(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    from pkpdsim.cli import _staged

    target = tmp_path / "result"
    stale = [tmp_path / f".result.{'a' * 32}", tmp_path / f".result.old.{'b' * 32}"]
    kept = [tmp_path / f".other.{'c' * 32}", tmp_path / ".result.keep"]
    for directory in [*stale, *kept]:
        directory.mkdir()
    _staged(target, lambda staging: staging.mkdir())
    assert target.is_dir()
    assert not any(directory.exists() for directory in stale)
    assert all(directory.exists() for directory in kept)
    assert sum("was killed" in record.message for record in caplog.records) == 2


def test_an_earlier_result_which_is_not_restored_is_named(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    from pkpdsim.cli import _staged

    target = tmp_path / "result"
    target.mkdir()
    (target / "earlier.txt").write_text("earlier", encoding="utf-8")
    original = Path.rename

    def rename(self: Path, other: Path) -> Path:
        if Path(other) == target:
            raise OSError("no rename onto the result")
        return original(self, other)

    monkeypatch.setattr(Path, "rename", rename)
    with pytest.raises(SystemExit) as info:
        _staged(target, lambda staging: staging.mkdir())
    assert info.value.code == 1
    out = "".join(capsys.readouterr().out.split())
    (earlier,) = tmp_path.glob(".result.old.*")
    assert (earlier / "earlier.txt").is_file()
    assert earlier.name in out


def test_an_error_of_the_run_itself_keeps_its_traceback(tmp_path: Path) -> None:
    from pkpdsim.cli import _staged

    def run(staging: Path) -> None:
        staging.mkdir()
        raise FileNotFoundError("a file of the run")

    with pytest.raises(FileNotFoundError, match="a file of the run"):
        _staged(tmp_path / "result", run)
    assert list(tmp_path.iterdir()) == []


def test_the_report_of_an_identifiability_names_its_analysis(
    fit: Path, tmp_path: Path
) -> None:
    report = identifiability_command(
        MODEL,
        [
            "--no-mapping-figures",
            "--no-fisher",
            *("--max-points", "3", "-o", str(tmp_path)),
            str(fit / "parameters.json"),
        ],
    )
    spec = json.loads((report / "report" / "spec.json").read_text(encoding="utf-8"))
    assert spec["fit"] == str((report / "analysis").resolve())
