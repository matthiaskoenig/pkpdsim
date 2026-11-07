"""Tests of the command `pkpd_release`, the build and the publishing of sbmlsim.

The release is sbmlsim's (tested there); these tests check the wrapper: its
options and defaults, the output of the build next to the checkout, the build
skipped by a publishing after a review and the checkout fast-forwarded after
publishing. The source and the release repository are real git repositories
(`release_git`); uv, the analyses of the run layer, docker and `gh` are
replaced by a fake runner which writes what they would write, the publishing
steps by a fake publisher; `release` gets them as its `runner` and its
`publisher`.
"""

import json
import logging
import shutil
from collections.abc import Callable, Sequence
from pathlib import Path

import pytest
from release_git import git
from sbmlsim.fit.records import adopt_fit
from sbmlsim.release.manifest import MANIFEST_FILE
from sbmlsim.release.settings import ReleaseSettings
from sbmlsim.release.source import ReleaseError, Runner, run
from sbmlsim.run.context import RUN_FILE
from sbmlutils.console import console

from pkpdsim.model import PKPDModel
from pkpdsim.release.cli import default_out, main, release
from pkpdsim_toy import TOY_PATH

pytestmark = pytest.mark.usefixtures("git_env")

VERSION = "0.2.0"
PYPI = '{ registry = "https://pypi.org/simple" }'
#: the pkpdsim of the lock of the source repository (`tests/data/release/uv.lock`)
PKPDSIM = "0.8.0"


def lock_text() -> str:
    """A lock of the release with the core libraries from PyPI."""
    packages = [("toy-model", VERSION, '{ editable = "." }')]
    packages += [("pkpdsim", PKPDSIM, PYPI)]
    packages += [(name, "1.0.0", PYPI) for name in ("sbmlsim", "sbmlutils")]
    blocks = [
        f'[[package]]\nname = "{name}"\nversion = "{version}"\nsource = {source}\n'
        for name, version, source in packages
    ]
    return 'version = 1\nrequires-python = ">=3.14"\n\n' + "\n".join(blocks)


class FakeRunner:
    """Runs git; fakes uv, the analyses of the run layer, docker and gh."""

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
            (cwd / "uv.lock").write_text(lock_text(), encoding="utf-8")
        elif "sbmlsim.run" in command:
            self._analysis(command[5], cwd / command[8], cwd)
        return ""

    def _analysis(self, kind: str, out: Path, cwd: Path) -> None:
        if kind == "factory":
            models = out / "models"
            models.mkdir(parents=True, exist_ok=True)
            shutil.copy(TOY_PATH / "models" / "toy.xml", models / "toy.xml")
            (models / "toy.md").write_text("# toy\n", encoding="utf-8")
            (out / "toy.omex").write_text("omex\n", encoding="utf-8")
            self._write(out / "figures" / "toy_model.png", "png")
        elif kind == "simulate":
            self._write(out / "Decay" / "Decay.html")
            self._write(out / "index.html")
        elif kind == "report":
            for record in (cwd / "src").glob("**/parameters.json"):
                report = out / record.parent.name
                self._write(report / "index.html")
                self._write(report / "report.pdf", "%PDF")
        (out / "run.json").write_text("{}", encoding="utf-8")

    @staticmethod
    def _write(path: Path, text: str = "<html/>") -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def ran(self, word: str) -> bool:
        """Tell whether a command with this word ran."""
        return any(word in " ".join(command) for command in self.commands)


class FakePublisher:
    """Records the calls of the publishing, commits, tags and pushes the release."""

    def __init__(self, fail: bool = False) -> None:
        self.calls: list[tuple[ReleaseSettings, Path, str, Path]] = []
        self.fail = fail

    def __call__(
        self,
        settings: ReleaseSettings,
        release_path: Path,
        version: str,
        archive: Path,
        runner: Runner,
    ) -> list[str]:
        self.calls.append((settings, release_path, version, archive))
        # every step is skipped when it is done, as the publishing of sbmlsim
        if git(release_path, "status", "--porcelain").strip():
            git(release_path, "add", "-A")
            git(release_path, "commit", "-m", f"Release {version}")
        if self.fail:
            raise ReleaseError("the push failed")
        git(release_path, "tag", version)
        git(release_path, "push", "origin", "main", "--tags")
        return ["commit", "tag", "push"]


@pytest.fixture
def runner() -> FakeRunner:
    """The fake runner of the commands of the release, which `release` gets."""
    return FakeRunner()


@pytest.fixture
def model(
    source_repo_model: PKPDModel,
    toy_fit: Path,
    release_settings: ReleaseSettings,
    as_model: Callable[[PKPDModel], PKPDModel],
) -> PKPDModel:
    """The toy model with an adopted fit and an archive, committed in the monorepo.

    It is the model object of its package, which the kind `release` loads.
    """
    model = source_repo_model.replace(
        adopted_ids=[toy_fit.name],
        archive_path=source_repo_model.path / "toy.omex",
        release=release_settings,
    )
    adopt_fit(toy_fit, model.fits_dir)
    root = model.path.parent.parent
    git(root, "add", "-A")
    git(root, "commit", "-m", "adopt the fit")
    return as_model(model)


@pytest.fixture
def clone(remote_and_clone: tuple[Path, Path]) -> Path:
    """The checkout of the release repository, with a file of an old layout."""
    _, clone = remote_and_clone
    (clone / "run_toy.py").write_text("print('old')\n", encoding="utf-8")
    git(clone, "add", "-A")
    git(clone, "commit", "-m", "old layout")
    git(clone, "push", "origin", "main")
    return clone


def do_release(
    model: PKPDModel,
    clone: Path | None,
    *,
    publish: bool,
    runner: FakeRunner,
    publisher: FakePublisher | None = None,
    out: Path | None = None,
) -> Path:
    return release(
        model,
        version=VERSION,
        repo=clone,
        publish=publish,
        pkpdsim_path=None,
        workers="2",
        out=out,
        runner=runner,
        publisher=publisher or FakePublisher(),
    )


def test_no_publish_builds_next_to_the_checkout(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    head = git(clone, "rev-parse", "HEAD")
    publisher = FakePublisher()
    path = do_release(model, clone, publish=False, runner=runner, publisher=publisher)
    out = clone.parent / f"{clone.name}-release-{VERSION}"
    assert default_out(clone, VERSION) == out
    assert path == out / "repository"
    assert json.loads((out / RUN_FILE).read_text(encoding="utf-8"))["status"] == "done"
    assert publisher.calls == []
    assert not runner.ran("docker") and not runner.ran("gh")
    # the checkout is only read
    assert git(clone, "rev-parse", "HEAD") == head
    assert not (clone / MANIFEST_FILE).exists()
    manifest = json.loads((path / MANIFEST_FILE).read_text(encoding="utf-8"))
    assert manifest["version"] == VERSION
    assert manifest["framework"] == "pkpdsim"
    assert manifest["framework_version"] == PKPDSIM
    assert (path / "toy.omex").is_file()
    assert (path / "fit" / model.adopted_ids[0] / "index.html").is_file()
    assert not (path / "run_toy.py").exists()
    notes = (path / "release-notes" / f"{VERSION}.md").read_text(encoding="utf-8")
    assert f"with pkpdsim {PKPDSIM}" in notes
    simulate = next(c for c in runner.commands if "specs/simulate.json" in c)
    assert simulate[-2:] == ["--workers", "2"]


def test_a_second_build_builds_again_into_the_output(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    path = do_release(model, clone, publish=False, runner=runner)
    (path / "docs" / "index.md").write_text("edited by hand\n", encoding="utf-8")
    assert do_release(model, clone, publish=False, runner=runner) == path
    assert (path / "docs" / "index.md").read_text(encoding="utf-8") != (
        "edited by hand\n"
    )
    assert sum("specs/simulate.json" in c for c in runner.commands) == 2


def test_the_output_is_given(
    model: PKPDModel, clone: Path, runner: FakeRunner, tmp_path: Path
) -> None:
    out = tmp_path / "build"
    assert (
        do_release(model, clone, publish=False, runner=runner, out=out)
        == out / "repository"
    )
    assert (out / "repository" / MANIFEST_FILE).is_file()


def test_publish_after_review_skips_the_build(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    path = do_release(model, clone, publish=False, runner=runner)
    builds = sum("specs/simulate.json" in c for c in runner.commands)
    publisher = FakePublisher()
    assert (
        do_release(model, clone, publish=True, runner=runner, publisher=publisher)
        == path
    )
    assert sum("specs/simulate.json" in c for c in runner.commands) == builds
    assert runner.ran("docker info") and runner.ran("gh auth status")
    assert model.release is not None
    assert publisher.calls == [(model.release, path, VERSION, path / "toy.omex")]
    # the checkout is fast-forwarded to the release
    assert git(clone, "rev-parse", "HEAD") == git(path, "rev-parse", "HEAD")
    assert VERSION in git(clone, "tag", "--list").split()


def test_publish_builds_and_publishes(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    publisher = FakePublisher()
    path = do_release(model, clone, publish=True, runner=runner, publisher=publisher)
    assert runner.ran("specs/simulate.json")
    assert len(publisher.calls) == 1
    assert git(clone, "rev-parse", "HEAD") == git(path, "rev-parse", "HEAD")


def test_a_changed_checkout_is_not_fast_forwarded(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    (clone / "notes.txt").write_text("a note\n", encoding="utf-8")
    head = git(clone, "rev-parse", "HEAD")
    path = do_release(model, clone, publish=True, runner=runner)
    assert git(clone, "rev-parse", "HEAD") == head
    assert git(path, "rev-parse", "HEAD") != head


def test_a_failed_publishing_resumes_without_a_build(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    with pytest.raises(ReleaseError, match="the push failed"):
        do_release(
            model,
            clone,
            publish=True,
            runner=runner,
            publisher=FakePublisher(fail=True),
        )
    builds = sum("specs/simulate.json" in c for c in runner.commands)
    publisher = FakePublisher()
    path = do_release(model, clone, publish=True, runner=runner, publisher=publisher)
    assert sum("specs/simulate.json" in c for c in runner.commands) == builds
    assert [call[1] for call in publisher.calls] == [path]


def test_without_a_checkout_the_repository_is_cloned_from_github(
    model: PKPDModel,
    remote_and_clone: tuple[Path, Path],
    runner: FakeRunner,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    missing = tmp_path / "missing" / "toy-model"
    with caplog.at_level(logging.INFO, logger="pkpdsim.release.cli"):
        path = do_release(model, missing, publish=False, runner=runner)
    assert f"The checkout {missing} does not exist" in caplog.text
    assert "clones the release repository from GitHub" in caplog.text
    assert path == tmp_path / "missing" / f"toy-model-release-{VERSION}" / "repository"
    assert not missing.exists()
    remote, _ = remote_and_clone
    assert git(path, "remote", "get-url", "origin").strip() == f"file://{remote}"


def test_the_default_checkout_is_the_sibling_of_the_source(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    """The source is `<tmp>/monorepo`, the release repository `<tmp>/toy-model`."""
    assert model.release is not None
    assert clone.name == model.release.repo_name
    path = do_release(model, None, publish=False, runner=runner)
    assert path == clone.parent / f"toy-model-release-{VERSION}" / "repository"
    assert any(
        c[:3] == ["git", "clone", "--no-hardlinks"] and str(clone) in c
        for c in runner.commands
    )


def test_a_refused_release_builds_nothing(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    git(clone, "tag", "0.3.0")
    git(clone, "push", "origin", "--tags")
    with pytest.raises(ReleaseError, match="The release is refused"):
        do_release(model, clone, publish=False, runner=runner)
    assert not runner.ran("uv lock")


def test_missing_release_settings(
    model: PKPDModel,
    clone: Path,
    runner: FakeRunner,
    as_model: Callable[[PKPDModel], PKPDModel],
) -> None:
    with pytest.raises(ReleaseError, match="has no release settings"):
        do_release(
            as_model(model.replace(release=None)), clone, publish=False, runner=runner
        )
    assert runner.commands == []


def test_publish_needs_an_archive(
    model: PKPDModel,
    clone: Path,
    runner: FakeRunner,
    as_model: Callable[[PKPDModel], PKPDModel],
) -> None:
    with pytest.raises(ReleaseError, match="archive"):
        do_release(
            as_model(model.replace(archive_path=None)),
            clone,
            publish=True,
            runner=runner,
        )
    assert runner.commands == []


def test_missing_tools_fail_before_the_build(model: PKPDModel, clone: Path) -> None:
    runner = FakeRunner(fail_on="docker info")
    with pytest.raises(ReleaseError, match="docker"):
        do_release(model, clone, publish=True, runner=runner)
    assert not runner.ran("uv lock")


def test_a_copy_of_the_model_is_refused(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    with pytest.raises(ReleaseError, match="model object of the package"):
        do_release(model.replace(title="Another"), clone, publish=False, runner=runner)
    assert runner.commands == []


def test_pkpdsim_path_needs_no_publish(
    model: PKPDModel, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    args = ["pkpdsim_toy", "--version", VERSION, "--pkpdsim-path", str(tmp_path)]
    with pytest.raises(SystemExit) as err:
        main(args)
    assert err.value.code == 2
    assert "--no-publish" in capsys.readouterr().err
    with pytest.raises(ReleaseError, match="--no-publish"):
        release(
            model,
            version=VERSION,
            repo=None,
            publish=True,
            pkpdsim_path=tmp_path,
            workers="auto",
            runner=FakeRunner(),
        )


def test_invalid_version_refused(model: PKPDModel, runner: FakeRunner) -> None:
    with pytest.raises(ReleaseError, match=r"X\.Y\.Z"):
        release(
            model,
            version="0.2",
            repo=None,
            publish=False,
            pkpdsim_path=None,
            workers="auto",
            runner=runner,
        )
    assert runner.commands == []


@pytest.mark.parametrize("workers", ["0", "-1", "two"])
def test_workers_validated(workers: str, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as err:
        main(["pkpdsim_toy", "--version", VERSION, "--workers", workers])
    assert err.value.code == 2
    assert "auto or a positive int" in capsys.readouterr().err


def test_main_reports_refusal(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def refuse(*args: object, **kwargs: object) -> Path:
        raise ReleaseError("the release repository has changes which no release made")

    monkeypatch.setattr("pkpdsim.release.cli.release", refuse)
    with pytest.raises(SystemExit) as err:
        main(["pkpdsim_toy", "--version", VERSION, "--no-publish"])
    assert err.value.code == 1
    assert "pkpd_release: error: " in capsys.readouterr().err


def test_main_no_publish(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    calls: list[dict[str, object]] = []

    def fake(model: PKPDModel, **kwargs: object) -> Path:
        calls.append(kwargs)
        return tmp_path

    monkeypatch.setattr("pkpdsim.release.cli.release", fake)
    main(["pkpdsim_toy", "--version", VERSION, "--no-publish", "--pkpdsim-path", "."])
    options = ["--out", "build", "--repo", "repo", "--workers", "3"]
    main(["pkpdsim_toy", "--version", VERSION, "--no-publish", *options])
    first, second = calls
    assert first["publish"] is False
    assert first["pkpdsim_path"] == Path.cwd()
    assert first["workers"] == "auto"
    assert first["repo"] is None
    assert first["out"] is None
    assert second["out"] == Path.cwd() / "build"
    assert second["repo"] == Path.cwd() / "repo"
    assert second["workers"] == "3"
    out = capsys.readouterr().out
    assert str(tmp_path) in out
    assert "review, then run without --no-publish" in out


@pytest.mark.parametrize(
    ("args", "line"),
    [
        (["--no-publish"], "{path}"),
        ([], f"Released {VERSION} from {{path}}."),
    ],
)
def test_main_prints_the_release_path_verbatim(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
    args: list[str],
    line: str,
) -> None:
    """The path is one line on a narrow terminal, markup characters included."""
    path = tmp_path / ("long-folder-name-" * 4) / "[draft]"
    monkeypatch.setattr(console, "_width", 20)
    monkeypatch.setattr("pkpdsim.release.cli.release", lambda *a, **k: path)
    main(["pkpdsim_toy", "--version", VERSION, *args])
    assert line.format(path=path) in capsys.readouterr().out.splitlines()


def test_a_build_with_a_local_pkpdsim_has_an_output_of_its_own(
    model: PKPDModel, clone: Path, runner: FakeRunner, tmp_path: Path
) -> None:
    """A review with a local pkpdsim, then a release, builds twice, no refusal."""
    local = tmp_path / "pkpdsim"
    local.mkdir()
    path = release(
        model,
        version=VERSION,
        repo=clone,
        publish=False,
        pkpdsim_path=local,
        workers="2",
        runner=runner,
    )
    assert path == default_out(clone, VERSION, local=True) / "repository"
    assert path.parent.name == f"toy-model-release-{VERSION}-local"
    publisher = FakePublisher()
    released = do_release(
        model, clone, publish=True, runner=runner, publisher=publisher
    )
    assert released == default_out(clone, VERSION) / "repository"
    assert [call[1] for call in publisher.calls] == [released]
    assert sum("specs/simulate.json" in c for c in runner.commands) == 2


def test_an_output_of_another_build_is_refused(
    model: PKPDModel, clone: Path, runner: FakeRunner, tmp_path: Path
) -> None:
    out = tmp_path / "build"
    do_release(model, clone, publish=False, runner=runner, out=out)
    builds = sum(c[:2] == ["uv", "lock"] for c in runner.commands)
    other = tmp_path / "other" / "toy-model"
    with pytest.raises(ReleaseError) as err:
        do_release(model, other, publish=False, runner=runner, out=out)
    assert f"remove {out} or pass --out" in str(err.value)
    assert "spec file" not in str(err.value)
    assert sum(c[:2] == ["uv", "lock"] for c in runner.commands) == builds
    (tmp_path / "used").mkdir()
    (tmp_path / "used" / "notes.txt").write_text("mine\n", encoding="utf-8")
    with pytest.raises(ReleaseError, match="or pass --out"):
        do_release(model, clone, publish=False, runner=runner, out=tmp_path / "used")


def test_an_unpushed_commit_of_the_checkout_is_refused(
    model: PKPDModel, clone: Path, runner: FakeRunner
) -> None:
    (clone / "local.txt").write_text("not pushed\n", encoding="utf-8")
    git(clone, "add", "-A")
    git(clone, "commit", "-m", "not pushed")
    out = default_out(clone, VERSION)
    with pytest.raises(ReleaseError, match="ahead of origin/main") as err:
        do_release(model, clone, publish=False, runner=runner)
    assert f"push or reset in the checkout {clone}, or remove {out}" in str(err.value)
    assert not runner.ran("uv lock")
