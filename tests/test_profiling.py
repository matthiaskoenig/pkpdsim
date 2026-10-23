"""The profile of a command with the sampling profiler of python 3.15."""

import sys
from pathlib import Path

import pytest

from pkpdsim import profiling
from pkpdsim.cli import COMMANDS, command_parser, run_simulation_main
from pkpdsim.profiling import profile_dir, profiler_available, run_profiled
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


def test_run_profiled_needs_the_profiler(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(profiling, "PROFILER_PYTHON", (sys.version_info[0] + 1, 0))
    with pytest.raises(RuntimeError, match="sampling profiler"):
        run_profiled(MODEL, "run_simulation", ["-a", "list_experiments"])


@needs_profiler
def test_a_profiled_command_writes_its_profile(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The flame graph, the heatmap, the statistics and the samples."""
    # the profiler runs the command in a python of its own, which imports the
    # toy model of the tests
    monkeypatch.setenv("PYTHONPATH", str(TESTS))
    folder = tmp_path / "profile"
    code = run_profiled(
        MODEL,
        "run_simulation",
        ["-a", "simulate", "-e", "all", "--results-dir", str(tmp_path / "results")],
        output=folder,
    )
    assert code == 0
    assert (tmp_path / "results").is_dir()
    for name in ["flamegraph.html", "profile.pstats", "profile.bin", "summary.txt"]:
        assert (folder / name).is_file(), name
    assert any((folder / "heatmap").iterdir())
    assert (
        (folder / "summary.txt")
        .read_text(encoding="utf-8")
        .startswith("Profile Stats:")
    )


@needs_profiler
def test_the_exit_code_of_a_failing_command_is_kept(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PYTHONPATH", str(TESTS))
    code = run_profiled(
        MODEL, "run_simulation", ["-a", "simulate"], output=tmp_path / "profile"
    )
    assert code != 0
