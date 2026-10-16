"""Tests of the build of a release repository."""

import dataclasses
from collections.abc import Sequence
from pathlib import Path

import pytest

from pkpdsim.release.build import build
from pkpdsim.release.source import ReleaseError
from pkpdsim_toy import MODEL

FIT_ID = "FIT_20260101_000000__abcde"
ARCHIVE = MODEL.path / "toy.omex"
MODEL_WITH_FIT = dataclasses.replace(MODEL, fits=[FIT_ID], archive_path=ARCHIVE)


class FakeRunner:
    """Records the commands and creates the files a real run would."""

    def __init__(
        self, fail_on: str | None = None, no_archive: bool = False, no_sim: bool = False
    ) -> None:
        self.no_archive = no_archive
        self.no_sim = no_sim
        self.commands: list[list[str]] = []
        self.streamed: list[list[str]] = []
        self.fail_on = fail_on

    def __call__(self, args: Sequence[str], cwd: Path, *, stream: bool = False) -> str:
        command = [str(a) for a in args]
        self.commands.append(command)
        if stream:
            self.streamed.append(command)
        joined = " ".join(command)
        if self.fail_on and self.fail_on in joined:
            raise ReleaseError(f"'{joined}' failed in {cwd}")
        if "--action factory" in joined:
            models = cwd / "src" / "pkpdsim_toy" / "models"
            models.mkdir(parents=True)
            (models / "toy.xml").write_text("<sbml/>")
            (models / "toy.md").write_text("# toy")
            (models / "other.txt").write_text("not copied")
            output = models.parent
            (output / "mit.txt").write_text("MIT")
            (output / "cc-by-sa-4.0.txt").write_text("CC")
            (output / "figures").mkdir()
            (output / "figures" / "model.svg").write_text("<svg/>")
            if not self.no_archive:
                (output / "toy.omex").write_text("omex")
        elif "--uncertainty" in command:
            self._write(cwd / "results/simulation_uncertainty/all/Exp/index.html")
        elif "--action simulate" in joined and not self.no_sim:
            self._write(cwd / "results/simulation/all/Exp/index.html")
        elif "run_report" in command:
            self._write(cwd / "results/fit_reports" / FIT_ID / "index.html")
        return ""

    @staticmethod
    def _write(path: Path) -> None:
        path.parent.mkdir(parents=True)
        path.write_text("<html/>")


def test_build_commands_in_order(tmp_path: Path) -> None:
    runner = FakeRunner()
    build(MODEL_WITH_FIT, tmp_path, runner=runner, workers="2")
    base = ["uv", "run", "run_simulation", "pkpdsim_toy"]
    simulate = [*base, "--action", "simulate", "--experiments", "all", "-r", "results"]
    simulate += ["--workers", "2"]
    report = ["uv", "run", "run_report", "pkpdsim_toy", "--adopted", FIT_ID]
    report += ["-o", "results/fit_reports"]
    assert runner.commands == [
        ["uv", "lock"],
        ["uv", "sync", "--locked"],
        [*base, "--action", "factory", "-r", "results"],
        simulate,
        [*simulate, "--uncertainty"],
        report,
    ]


def test_build_streams_and_prints_stages(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The steps stream their output and every stage prints a line."""
    runner = FakeRunner()
    build(MODEL_WITH_FIT, tmp_path, runner=runner)
    assert runner.streamed == runner.commands
    out = capsys.readouterr().out
    for line in (
        "Lock and sync of the environment",
        "Models",
        "Simulations",
        "Simulations with the confidence interval",
        f"Fit report {FIT_ID}",
    ):
        assert line in out


def test_build_copies_outputs(tmp_path: Path) -> None:
    written = build(MODEL_WITH_FIT, tmp_path, runner=FakeRunner())
    assert (tmp_path / "models/toy.xml").is_file()
    assert (tmp_path / "models/toy.md").is_file()
    assert not (tmp_path / "models/other.txt").exists()
    assert (tmp_path / "simulations/Exp/index.html").is_file()
    assert (tmp_path / "simulations_uncertainty/Exp/index.html").is_file()
    assert (tmp_path / "fit" / FIT_ID / "index.html").is_file()
    assert (tmp_path / "toy.omex").is_file()
    assert (tmp_path / "mit.txt").is_file()
    assert (tmp_path / "cc-by-sa-4.0.txt").is_file()
    assert (tmp_path / "figures/model.svg").is_file()
    assert Path("models/toy.xml") in written
    assert Path("uv.lock") in written
    assert Path("simulations") in written
    assert Path("fit") / FIT_ID in written


def test_build_replaces_stale_reports(tmp_path: Path) -> None:
    old = tmp_path / "fit" / "OLD_20200101_000000__aaaaa"
    old.mkdir(parents=True)
    (old / "index.html").write_text("old")
    build(MODEL_WITH_FIT, tmp_path, runner=FakeRunner())
    assert not old.exists()
    assert (tmp_path / "fit" / FIT_ID).is_dir()


def test_build_failure_names_step(tmp_path: Path) -> None:
    with pytest.raises(ReleaseError, match="--uncertainty"):
        build(MODEL_WITH_FIT, tmp_path, runner=FakeRunner(fail_on="--uncertainty"))


def test_build_missing_declared_archive_raises(tmp_path: Path) -> None:
    with pytest.raises(ReleaseError, match=r"toy\.omex"):
        build(MODEL_WITH_FIT, tmp_path, runner=FakeRunner(no_archive=True))


def test_build_removes_stale_models(tmp_path: Path) -> None:
    (tmp_path / "models").mkdir()
    (tmp_path / "models/old.xml").write_text("old")
    build(MODEL_WITH_FIT, tmp_path, runner=FakeRunner())
    assert not (tmp_path / "models/old.xml").exists()


def test_build_does_not_reuse_stale_results(tmp_path: Path) -> None:
    stale = tmp_path / "results/simulation/all/Exp"
    stale.mkdir(parents=True)
    (stale / "index.html").write_text("stale")
    with pytest.raises(ReleaseError, match="simulate"):
        build(MODEL_WITH_FIT, tmp_path, runner=FakeRunner(no_sim=True))
