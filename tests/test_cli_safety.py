"""The commands replace a result only when its run succeeded and never touch their inputs.

A kind of the run layer runs into a staging directory next to the result; the
earlier result is replaced only after the run succeeded, and an output which is
an input of the command, contains one or lies inside one is refused before
anything happens.
"""

import json
import math
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from textwrap import dedent

import pytest
from sbmlsim.cgroup import available_cpus
from sbmlsim.experiment.runner import EXPERIMENTS_PER_WORKER
from sbmlsim.fit.records import adopt_fit
from sbmlsim.run import ConfigurationError

import pkpdsim.cli
from pkpdsim.cli import (
    credibility_command,
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
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        f".{report.name}.lock",
        report.name,
    ]


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
    assert sorted(p.name for p in out.parent.iterdir()) == [".studies.lock", "studies"]


@pytest.mark.parametrize("dedicated_data", [False, True])
def test_simulation_protects_an_existing_study_directory(
    tmp_path: Path,
    as_model: AsModel,
    capsys: pytest.CaptureFixture[str],
    *,
    dedicated_data: bool,
) -> None:
    results = tmp_path / "data"
    model = as_model(MODEL.replace(data_path=results)) if dedicated_data else MODEL
    study = results / "simulation" / "all"
    study.mkdir(parents=True)
    (study / "study.tsv").write_text("study data", encoding="utf-8")
    # A website alone does not make a directory a simulation result.
    (study / "index.html").write_text("study page", encoding="utf-8")
    before = _snapshot(study)
    with pytest.raises(SystemExit) as info:
        run_command(
            model, ["-a", "simulate", "-e", "all", "-w", "1", "-r", str(results)]
        )
    out = capsys.readouterr()
    assert info.value.code == (2 if dedicated_data else 1)
    if dedicated_data:
        assert "input" in out.err
    else:
        assert "give another output" in out.out
    assert _snapshot(study) == before


def test_a_failed_experiment_preserves_the_published_report(
    tmp_path: Path, as_model: AsModel, capsys: pytest.CaptureFixture[str]
) -> None:
    arguments = ["-a", "simulate", "-e", "all", "-w", "1", "-r", str(tmp_path)]
    run_command(MODEL, arguments)
    output = tmp_path / "simulation" / "all"
    before = _snapshot(output)
    capsys.readouterr()
    broken = as_model(MODEL.replace(experiments=".broken:BROKEN"))
    with pytest.raises(SystemExit) as info:
        run_command(broken, arguments)
    assert info.value.code == 1
    out = capsys.readouterr().out
    assert "ToyCannotInitialize" in out
    assert "earlier results kept" in out
    assert _snapshot(output) == before
    assert sorted(path.name for path in output.parent.iterdir()) == [".all.lock", "all"]


def test_a_failed_retry_restores_the_only_earlier_result(tmp_path: Path) -> None:
    target = tmp_path / "simulation_uncertainty" / "all"
    target.mkdir(parents=True)
    (target / "index.html").write_text("prior results", encoding="utf-8")
    (target / "run.json").write_text('{"kind":"simulate"}', encoding="utf-8")
    before = _snapshot(target)
    # Terminate the writer immediately after it moved the earlier result aside.
    command = dedent("""
        import os
        import sys
        from pathlib import Path
        from pkpdsim.cli import _staged

        target = Path(sys.argv[1]).resolve()
        rename = Path.rename

        def interrupt(self, destination):
            renamed = rename(self, destination)
            if self == target:
                os._exit(19)
            return renamed

        Path.rename = interrupt
        _staged(target, Path.mkdir, kind="simulate")
        """)
    interrupted = subprocess.run(
        [sys.executable, "-c", command, str(target)],
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    assert interrupted.returncode == 19, interrupted.stderr
    assert not target.exists()
    (earlier,) = target.parent.glob(".all.old.*")
    assert _snapshot(earlier) == before
    with pytest.raises(SystemExit) as info:
        run_command(
            MODEL, ["-a", "simulate", "-e", "all", "-w", "1", "-r", str(tmp_path), "-u"]
        )
    assert info.value.code == 2
    assert _snapshot(target) == before
    assert not earlier.exists()


def test_a_second_process_keeps_the_live_staging_directory(tmp_path: Path) -> None:
    from pkpdsim.cli import _staged

    parent = tmp_path / "results"
    parent.mkdir()
    target = parent / "all"
    alias = tmp_path / "alias"
    try:
        alias.symlink_to(parent, target_is_directory=True)
    except OSError:
        # Windows can require administrator rights to create symbolic links.
        (parent / "nested").mkdir()
        alias = parent / "nested" / ".."

    def run(staging: Path) -> None:
        staging.mkdir()
        (staging / "index.html").write_text("first run", encoding="utf-8")
        attempt = subprocess.run(
            [
                sys.executable,
                "-c",
                "from pathlib import Path; from pkpdsim.cli import _staged; "
                "import sys; _staged(Path(sys.argv[1]), Path.mkdir, kind='simulate')",
                str(alias / target.name),
            ],
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
        assert attempt.returncode == 1
        assert "another run" in attempt.stdout
        assert (staging / "index.html").is_file()

    _staged(target, run, kind="simulate")
    assert (target / "index.html").read_text(encoding="utf-8") == "first run"


@pytest.mark.parametrize("function", [report_command, identifiability_command])
def test_a_report_of_plain_parameters_protects_the_data_directory(
    function: Callable[[PKPDModel, list[str]], Path],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    data = tmp_path / "data"
    model = MODEL.replace(data_path=data)
    with pytest.raises(SystemExit) as info:
        function(model, ["-o", str(data), str(tmp_path / "raw.json")])
    assert info.value.code == 2
    assert "input" in capsys.readouterr().err
    assert not data.exists()


def test_credibility_protects_a_separate_input_inside_its_output(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "assessment"
    source = out / "archive" / "credibility.json"
    source.parent.mkdir(parents=True)
    for record in [out / "credibility.json", source]:
        record.write_text("{}", encoding="utf-8")
    before = _snapshot(out)
    with pytest.raises(SystemExit) as info:
        credibility_command(MODEL, ["-o", str(out), "--from-json", str(source)])
    assert info.value.code == 2
    assert "input" in capsys.readouterr().err
    assert _snapshot(out) == before


@pytest.mark.parametrize("directory_link", [False, True])
def test_credibility_keeps_a_nested_record_aliased_as_the_direct_input(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    credibility_model: PKPDModel,
    as_model: AsModel,
    *,
    directory_link: bool,
) -> None:
    out = tmp_path / "assessment"
    direct = out / "credibility.json"
    archived = out / "archive" / "credibility.json"
    subdir = archived.parent / "subdir"
    model = as_model(credibility_model)
    credibility_command(model, ["--no-verification", "-o", str(out)])
    subdir.mkdir(parents=True)
    archived.write_bytes(direct.read_bytes())
    try:
        if directory_link:
            link = out / "link"
            link.symlink_to(subdir, target_is_directory=True)
            source = link / ".." / archived.name
        else:
            direct.unlink()
            direct.symlink_to(archived)
            source = direct
    except OSError:
        pytest.skip("symbolic links are unavailable")
    if source.resolve() != archived.resolve():
        pytest.skip("this platform resolves the alias to the direct record")
    before = _snapshot(out)
    with pytest.raises(SystemExit) as info:
        credibility_command(model, ["-o", str(out), "--from-json", str(source)])
    assert info.value.code == 2
    assert "input" in capsys.readouterr().err
    assert _snapshot(out) == before


def test_a_failed_move_is_a_message_and_keeps_the_earlier_result(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    target = tmp_path / "result"
    target.mkdir()
    (target / "earlier.txt").write_text("earlier", encoding="utf-8")
    (target / "run.json").write_text('{"kind":"simulate"}', encoding="utf-8")

    def run(staging: Path) -> None:
        (staging / "inner").mkdir(parents=True)
        (staging / "inner" / "index.html").write_text("new", encoding="utf-8")
        # the records of the run of the same name as a file of the result
        (staging / "index.html").write_text("clash", encoding="utf-8")

    from pkpdsim.cli import _staged

    with pytest.raises(SystemExit) as info:
        _staged(target, run, kind="simulate", inner="inner")
    assert info.value.code == 1
    assert "index.html" in capsys.readouterr().out
    assert sorted(p.name for p in tmp_path.iterdir()) == [".result.lock", "result"]
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
    args = ["-a", "simulate", "-e", "studies", "-r", str(tmp_path)]
    run_command(MODEL, [*args, "-w", "auto"])
    assert seen == [None]
    out = tmp_path / "simulation" / "studies"
    assert (out / "index.html").is_file()
    assert list((out / "_figures").iterdir())
    run_command(MODEL, [*args, "-w", "1"])
    assert (out / "run.json").is_file()


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
    # the console wraps a long path and writes the separators of the platform;
    # with markup the `[b]` of the path would vanish as a style
    out = "".join(capsys.readouterr().out.split()).replace("\\", "/")
    assert "earlierresultsin" in out
    assert "a[b]/simulation/studies'arereplaced" in out


def test_the_stagings_of_a_killed_run_are_removed(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    from pkpdsim.cli import _staged

    target = tmp_path / "result"
    target.mkdir()
    stale = [tmp_path / f".result.{'a' * 32}", tmp_path / f".result.old.{'b' * 32}"]
    kept = [tmp_path / f".other.{'c' * 32}", tmp_path / ".result.keep"]
    for directory in [*stale, *kept]:
        directory.mkdir()
    _staged(target, Path.mkdir, kind="simulate")
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
    (target / "run.json").write_text('{"kind":"simulate"}', encoding="utf-8")
    original = Path.rename

    def rename(self: Path, other: Path) -> Path:
        if Path(other) == target:
            raise OSError("no rename onto the result")
        return original(self, other)

    monkeypatch.setattr(Path, "rename", rename)
    with pytest.raises(SystemExit) as info:
        _staged(target, Path.mkdir, kind="simulate")
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
        _staged(tmp_path / "result", run, kind="simulate")
    assert [p.name for p in tmp_path.iterdir()] == [".result.lock"]


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
