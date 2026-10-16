"""Tests of the command `pkpd_release`, the complete release of a model repository.

The source and the release repository are real git repositories (`release_git`);
the build tools (`uv`, `run_simulation`, `run_report`), docker and `gh` are
replaced by a fake runner which writes what they would write, the publishing
steps by a fake publisher.
"""

import dataclasses
import hashlib
import json
import shutil
from collections.abc import Sequence
from pathlib import Path

import pytest
from release_git import git

from pkpdsim import __version__
from pkpdsim.model import PKPDModel
from pkpdsim.release.cli import main, release
from pkpdsim.release.manifest import MANIFEST_FILE
from pkpdsim.release.settings import ReleaseSettings
from pkpdsim.release.source import ReleaseError, Runner, run
from pkpdsim.uncertainty import adopt_fit
from pkpdsim_toy import MODEL, TOY_PATH

pytestmark = pytest.mark.usefixtures("git_env")

VERSION = "0.2.0"
PYPI = '{ registry = "https://pypi.org/simple" }'


def lock_text() -> str:
    """A lock of the release with the core libraries from PyPI."""
    packages = [("toy-model", VERSION, '{ editable = "." }')]
    packages += [(name, __version__, PYPI) for name in ("pkpdsim",)]
    packages += [(name, "1.0.0", PYPI) for name in ("sbmlsim", "sbmlutils")]
    blocks = [
        f'[[package]]\nname = "{name}"\nversion = "{version}"\nsource = {source}\n'
        for name, version, source in packages
    ]
    return 'version = 1\nrequires-python = ">=3.13"\n\n' + "\n".join(blocks)


class FakeRunner:
    """Runs git; fakes uv, the commands of the model, docker and gh."""

    def __init__(self, fail_on: str | None = None) -> None:
        self.commands: list[list[str]] = []
        self.fail_on = fail_on

    def __call__(self, args: Sequence[str], cwd: Path, *, stream: bool = False) -> str:
        command = [str(a) for a in args]
        self.commands.append(command)
        joined = " ".join(command)
        if self.fail_on and self.fail_on in joined:
            raise ReleaseError(f"'{joined}' failed in {cwd}")
        if command[0] == "git":
            return run(command, cwd)
        if command[:2] == ["uv", "lock"]:
            (cwd / "uv.lock").write_text(lock_text())
        elif "--action factory" in joined:
            output = cwd / "src" / "pkpdsim_toy"
            models = output / "models"
            models.mkdir(parents=True, exist_ok=True)
            shutil.copy(TOY_PATH / "models" / "toy.xml", models / "toy.xml")
            (models / "toy.md").write_text("# toy\n")
            (output / "toy.omex").write_text("omex\n")
            self._write(output / "figures" / "toy_model.png", "png")
        elif "--uncertainty" in command:
            self._write(cwd / "results/simulation_uncertainty/all/Decay/index.html")
        elif "--action simulate" in joined:
            self._write(cwd / "results/simulation/all/Decay/index.html")
        elif "run_report" in command:
            fit_id = command[command.index("--adopted") + 1]
            self._write(cwd / "results/fit_reports" / fit_id / "index.html")
        return ""

    @staticmethod
    def _write(path: Path, text: str = "<html/>") -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def ran(self, word: str) -> bool:
        """Tell whether a command with this word ran."""
        return any(word in command for command in self.commands)


class FakePublisher:
    """Records the calls of the publishing steps."""

    def __init__(self) -> None:
        self.calls: list[tuple[ReleaseSettings, Path, str, Path]] = []

    def __call__(
        self,
        settings: ReleaseSettings,
        release_path: Path,
        version: str,
        archive: Path,
        runner: Runner,
    ) -> list[str]:
        self.calls.append((settings, release_path, version, archive))
        return []


@pytest.fixture
def model(
    source_repo_model: PKPDModel, toy_fit: Path, release_settings: ReleaseSettings
) -> PKPDModel:
    """The toy model with an adopted fit and an archive, committed in the monorepo."""
    model = dataclasses.replace(
        source_repo_model,
        fits=[toy_fit.name],
        archive_path=source_repo_model.path / "toy.omex",
        release=release_settings,
    )
    adopt_fit(toy_fit, model.fits_dir)
    root = model.path.parent.parent
    git(root, "add", "-A")
    git(root, "commit", "-m", "adopt the fit")
    return model


@pytest.fixture
def clone(remote_and_clone: tuple[Path, Path]) -> Path:
    """The release repository with a file of an old layout."""
    _, clone = remote_and_clone
    (clone / "run_toy.py").write_text("print('old')\n")
    git(clone, "add", "-A")
    git(clone, "commit", "-m", "old layout")
    git(clone, "push", "origin", "main")
    return clone


def do_release(
    model: PKPDModel,
    clone: Path,
    *,
    publish: bool,
    runner: FakeRunner | None = None,
    publisher: FakePublisher | None = None,
) -> Path:
    return release(
        model,
        version=VERSION,
        repo=clone,
        publish=publish,
        pkpdsim_path=None,
        workers="2",
        runner=runner or FakeRunner(),
        publisher=publisher or FakePublisher(),
    )


def test_missing_release_settings(model: PKPDModel, clone: Path) -> None:
    runner = FakeRunner()
    with pytest.raises(ReleaseError, match=r"MODEL\.release is not set"):
        do_release(dataclasses.replace(model, release=None), clone, publish=False)
    assert model.release is not None
    incomplete = dataclasses.replace(model.release, doi="", funding="")
    with pytest.raises(ReleaseError, match="doi, funding"):
        do_release(
            dataclasses.replace(model, release=incomplete),
            clone,
            publish=False,
            runner=runner,
        )
    assert runner.commands == []


def test_pkpdsim_path_needs_no_publish(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    args = ["pkpdsim_toy", "--version", VERSION, "--pkpdsim-path", str(tmp_path)]
    with pytest.raises(SystemExit) as err:
        main(args)
    assert err.value.code == 2
    assert "--no-publish" in capsys.readouterr().err
    with pytest.raises(ReleaseError, match="--no-publish"):
        release(
            MODEL,
            version=VERSION,
            repo=None,
            publish=True,
            pkpdsim_path=tmp_path,
            workers="auto",
            runner=FakeRunner(),
        )


def test_no_publish_stops_before_commit(model: PKPDModel, clone: Path) -> None:
    head = git(clone, "rev-parse", "HEAD")
    runner, publisher = FakeRunner(), FakePublisher()
    path = do_release(model, clone, publish=False, runner=runner, publisher=publisher)
    assert path == clone
    assert publisher.calls == []
    assert not runner.ran("docker") and not runner.ran("gh")
    assert git(clone, "rev-parse", "HEAD") == head
    manifest = json.loads((clone / MANIFEST_FILE).read_text())
    assert manifest["version"] == VERSION
    assert "uv.lock" in manifest["files"]
    assert "docs/images/toy_model.png" in manifest["files"]
    assert (clone / "uv.lock").is_file()
    assert (clone / "docs" / "images" / "toy_model.png").is_file()
    assert (clone / "release-notes" / f"{VERSION}.md").is_file()
    assert (clone / "fit" / model.fits[0] / "index.html").is_file()
    assert not (clone / "run_toy.py").exists()


def test_stage_lines(
    model: PKPDModel, clone: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    do_release(model, clone, publish=False)
    out = capsys.readouterr().out
    for line in (
        "Export of pkpdsim_toy",
        "Project files",
        "Lock and sync of the environment",
        "Simulations with the confidence interval",
        "Release notes",
        "Docs and README",
        "Prune and manifest",
    ):
        assert line in out


def test_publish_after_review_skips_build(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    runner, publisher = FakeRunner(), FakePublisher()
    do_release(model, clone, publish=True, runner=runner, publisher=publisher)
    assert not runner.ran("run_simulation")
    assert runner.ran("docker") and runner.ran("gh")
    assert model.release is not None
    assert publisher.calls == [(model.release, clone, VERSION, clone / "toy.omex")]


def test_hand_edit_after_review_rebuilds(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    (clone / "docs" / "index.md").write_text("edited by hand\n")
    runner, publisher = FakeRunner(), FakePublisher()
    do_release(model, clone, publish=True, runner=runner, publisher=publisher)
    assert runner.ran("run_simulation")
    assert (clone / "docs" / "index.md").read_text() != "edited by hand\n"
    assert len(publisher.calls) == 1


def test_resume_existing_tag(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    git(clone, "add", "-A")
    git(clone, "commit", "-m", f"Release {VERSION}")
    git(clone, "tag", VERSION)
    runner, publisher = FakeRunner(), FakePublisher()
    do_release(model, clone, publish=False, runner=runner, publisher=publisher)
    assert not runner.ran("run_simulation")
    do_release(model, clone, publish=True, runner=runner, publisher=publisher)
    assert not runner.ran("run_simulation")
    assert len(publisher.calls) == 1


def test_dirty_unrelated_changes_refused(model: PKPDModel, clone: Path) -> None:
    (clone / "notes.txt").write_text("not from a release\n")
    runner = FakeRunner()
    with pytest.raises(ReleaseError, match="changes which no release made"):
        do_release(model, clone, publish=False, runner=runner)
    assert not runner.ran("run_simulation")
    assert not (clone / "src").exists()


def test_unrelated_change_after_review_refused(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    (clone / "notes.txt").write_text("not from a release\n")
    with pytest.raises(ReleaseError, match="changes which no release made"):
        do_release(model, clone, publish=True)


def test_missing_tools_fail_before_build(model: PKPDModel, clone: Path) -> None:
    runner = FakeRunner(fail_on="docker info")
    with pytest.raises(ReleaseError, match="docker"):
        do_release(model, clone, publish=True, runner=runner)
    assert not runner.ran("run_simulation")


def test_local_lock_after_review_rebuilds(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    lock = clone / "uv.lock"
    lock.write_text(lock.read_text().replace(PYPI, '{ directory = "../pkpdsim" }'))
    manifest = json.loads((clone / MANIFEST_FILE).read_text())
    manifest["files"]["uv.lock"] = hashlib.sha256(lock.read_bytes()).hexdigest()
    (clone / MANIFEST_FILE).write_text(json.dumps(manifest))
    runner, publisher = FakeRunner(), FakePublisher()
    do_release(model, clone, publish=True, runner=runner, publisher=publisher)
    assert runner.ran("run_simulation")
    assert len(publisher.calls) == 1


def test_fit_without_record_refused(model: PKPDModel, clone: Path) -> None:
    fit_id = model.fits[0].replace("__", "__f")
    runner = FakeRunner()
    with pytest.raises(ReleaseError, match=fit_id):
        do_release(
            dataclasses.replace(model, fits=[fit_id]),
            clone,
            publish=False,
            runner=runner,
        )
    assert not runner.ran("run_simulation")


def test_publish_needs_archive(model: PKPDModel, clone: Path) -> None:
    with pytest.raises(ReleaseError, match="archive"):
        do_release(dataclasses.replace(model, archive_path=None), clone, publish=True)


def test_main_reports_refusal(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def refuse(*args: object, **kwargs: object) -> Path:
        raise ReleaseError("the release repository has changes which no release made")

    monkeypatch.setattr("pkpdsim.release.cli.release", refuse)
    with pytest.raises(SystemExit) as err:
        main(["pkpdsim_toy", "--version", VERSION, "--no-publish"])
    assert err.value.code == 1
    assert "no release made" in capsys.readouterr().err


def test_main_no_publish(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    calls: list[dict[str, object]] = []

    def fake(model: PKPDModel, **kwargs: object) -> Path:
        calls.append(kwargs)
        return tmp_path

    monkeypatch.setattr("pkpdsim.release.cli.release", fake)
    main(["pkpdsim_toy", "--version", VERSION, "--no-publish", "--pkpdsim-path", "."])
    (kwargs,) = calls
    assert kwargs["publish"] is False
    assert kwargs["pkpdsim_path"] == Path.cwd()
    assert kwargs["workers"] == "auto"
    assert kwargs["repo"] is None
    out = capsys.readouterr().out
    assert str(tmp_path) in out
    assert "review, then run without --no-publish" in out


def commit_source(model: PKPDModel) -> None:
    """Move the source repository on by a commit which changes the package."""
    (model.path / "extra.py").write_text("EXTRA = 1\n")
    root = model.path.parent.parent
    git(root, "add", "-A")
    git(root, "commit", "-m", "more")


def commit_release(clone: Path, tag: bool) -> None:
    """Commit the build as the publisher does, with or without the tag."""
    git(clone, "add", "-A")
    git(clone, "commit", "-m", f"Release {VERSION}")
    if tag:
        git(clone, "tag", VERSION)


def test_failed_build_rerun(model: PKPDModel, clone: Path) -> None:
    with pytest.raises(ReleaseError, match="uncertainty"):
        do_release(
            model, clone, publish=False, runner=FakeRunner(fail_on="--uncertainty")
        )
    do_release(model, clone, publish=False)
    assert (clone / MANIFEST_FILE).is_file()
    assert not (clone / ".git" / "pkpd_release_build").exists()


def test_interrupted_rebuild_after_review(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    commit_source(model)  # the export writes a file which the manifest lacks
    with pytest.raises(ReleaseError, match="uncertainty"):
        do_release(
            model, clone, publish=True, runner=FakeRunner(fail_on="--uncertainty")
        )
    publisher = FakePublisher()
    do_release(model, clone, publish=True, publisher=publisher)
    assert len(publisher.calls) == 1
    assert (clone / "src" / "pkpdsim_toy" / "extra.py").is_file()


def test_results_survive_and_do_not_block(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    (clone / "results" / "x").write_text("mine\n")
    do_release(model, clone, publish=False)
    assert (clone / "results" / "x").is_file()
    assert (clone / "results/simulation/all/Decay/index.html").is_file()


def test_resume_tag_with_moved_source(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    commit_release(clone, tag=True)
    commit_source(model)
    runner, publisher = FakeRunner(), FakePublisher()
    do_release(model, clone, publish=True, runner=runner, publisher=publisher)
    assert not runner.ran("run_simulation")
    assert len(publisher.calls) == 1


def test_resume_release_commit_without_tag(model: PKPDModel, clone: Path) -> None:
    do_release(model, clone, publish=False)
    commit_release(clone, tag=False)
    head = git(clone, "rev-parse", "HEAD")
    commit_source(model)
    runner, publisher = FakeRunner(), FakePublisher()
    do_release(model, clone, publish=True, runner=runner, publisher=publisher)
    assert not runner.ran("run_simulation")
    assert len(publisher.calls) == 1
    assert git(clone, "rev-parse", "HEAD") == head
    assert not git(clone, "status", "--porcelain").strip()


def test_invalid_version_refused(model: PKPDModel, clone: Path) -> None:
    runner = FakeRunner()
    with pytest.raises(ReleaseError, match="not a version"):
        release(
            model,
            version="1.0",
            repo=clone,
            publish=True,
            pkpdsim_path=None,
            workers="auto",
            runner=runner,
        )
    assert runner.commands == []


@pytest.mark.parametrize("workers", ["0", "-2", "many"])
def test_workers_validated(workers: str, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as err:
        main(["pkpdsim_toy", "--version", VERSION, "--workers", workers])
    assert err.value.code == 2
    assert "--workers" in capsys.readouterr().err
