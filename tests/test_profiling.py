"""The profile of a command with the sampling profiler of python 3.15.

The profiler is `sbmlsim.run.profiling`, tested in sbmlsim; the commands run
it on `python -m pkpdsim.cli` and write the profile into the results of the
model.
"""

from pathlib import Path

import pytest
from sbmlsim.run.profiling import profiler_available

from pkpdsim.cli import COMMANDS, command_parser, profile_dir, run_simulation_main
from pkpdsim_toy import MODEL

TESTS = Path(__file__).parent

needs_profiler = pytest.mark.skipif(
    not profiler_available(), reason="the sampling profiler of python 3.15"
)


def test_the_help_names_the_option() -> None:
    assert "--profile" in command_parser(COMMANDS["run_simulation"]).format_help()


def test_the_profile_is_in_the_results_of_the_model() -> None:
    folder = profile_dir(MODEL, "run_simulation")
    assert folder.parent == MODEL.results_dir / "profiling"
    assert folder.name.startswith("run_simulation_")


def test_without_the_profiler_the_option_is_an_error(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """With a python before 3.15 the command does not run."""
    monkeypatch.setattr("pkpdsim.cli.profiler_available", lambda: False)
    with pytest.raises(SystemExit) as error:
        run_simulation_main(["pkpdsim_toy", "--profile", "-a", "list_experiments"])
    assert error.value.code == 2
    assert "python 3.15" in capsys.readouterr().err


@needs_profiler
def test_a_profiled_command_writes_its_profile(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The flame graph, the heatmap, the statistics and the samples."""
    # the profiler runs the command in a python of its own, which imports the
    # toy model of the tests
    monkeypatch.setenv("PYTHONPATH", str(TESTS))
    folder = tmp_path / "profile"
    monkeypatch.setattr("pkpdsim.cli.profile_dir", lambda model, command: folder)
    results = tmp_path / "results"
    run_simulation_main(
        [
            *("pkpdsim_toy", "--profile", "-a", "simulate", "-e", "all", "-w", "1"),
            *("--results-dir", str(results)),
        ]
    )
    assert (results / "simulation" / "all" / "run.json").is_file()
    for name in ["flamegraph.html", "profile.pstats", "profile.bin", "summary.txt"]:
        assert (folder / name).is_file(), name


@needs_profiler
def test_the_exit_code_of_a_failing_command_is_kept(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PYTHONPATH", str(TESTS))
    monkeypatch.setattr(
        "pkpdsim.cli.profile_dir", lambda model, command: tmp_path / "profile"
    )
    with pytest.raises(SystemExit) as error:
        run_simulation_main(["pkpdsim_toy", "--profile", "-a", "simulate"])
    assert error.value.code == 2
