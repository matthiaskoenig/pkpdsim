"""Test session configuration.

sbmlutils resolves the annotations of a model with the registry of
identifiers.org, which pymetadata downloads into its cache on first use. The
tests run in parallel, so on a machine without the cache every worker which
builds an annotated model, and every snippet of the documentation which does,
would download the registry and replace the cache file at the same time; on
Windows a replace of a file which another process reads fails
(matthiaskoenig/pymetadata#102). The registry is therefore loaded once in the
controller of the session, before the workers start, and the workers only read
the cache.

The SBML models of the toy model packages (`tests/pkpdsim_toy` and the toy
model per sex `tests/pkpdsim_sexes`) are created in the controller of the
session as well, before the workers start; the tests only read them.

A short fit of the toy model with its identifiability, adopted into a folder of
fits, is shared by the tests of the records (`toy_fit`, `fits_dir`); it is made
once per module which uses it.

A PBPK build writes its flat model through `pkpdsim.pbpk.builder.flat_model`,
which validates the comp model, flattens it and validates the flat model, most
of the time of a build. The flat model depends on the content of the comp model
alone, and many tests build the same model, so the flat models are cached by
the content of their comp model (without the dates of the model history) in a
folder of the session which every worker shares (`flat_model_cache`): a comp
model which another test has flattened is not validated and flattened again. A
test which checks what the validation or the flattening logs is marked
`uncached_build`.

The tests of the release share git repositories (`git_env`, `source_repo_model`,
`remote_and_clone`, `release_settings`, see `release_git`).

`--shard I/N` runs the share `I` of `N` of the tests, a job of the matrix of
the CI each (see `sharding`).
"""

import dataclasses
import hashlib
import os
import re
import shutil
import sys
from contextlib import suppress
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from pymetadata.webservices.registry import get_registry
from pymetadata.webservices.webservice import WebserviceError
from release_git import GIT_ENV, git, init_repo
from sharding import in_shard, parse_shard

from pkpdsim.cli import fit_command
from pkpdsim.model import PKPDModel
from pkpdsim.release import source as source_module
from pkpdsim.release.settings import ReleaseSettings
from pkpdsim.uncertainty import adopt_fit
from pkpdsim_sexes import MODEL_BASE_PATH as SEXES_BASE_PATH
from pkpdsim_sexes.factory import create_sexes_models
from pkpdsim_toy import MODEL, MODEL_BASE_PATH
from pkpdsim_toy.factory import create_toy_models

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

RELEASE_DATA = Path(__file__).parent / "data" / "release"


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add `--shard`, the share of the tests which runs."""
    parser.addoption(
        "--shard",
        default=None,
        metavar="I/N",
        help="run the share I of N of the tests (1 <= I <= N), by a hash of the test id",
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Deselect the tests outside of the share of `--shard`."""
    text = config.getoption("--shard")
    if text is None:
        return
    index, count = parse_shard(text)
    deselected = [item for item in items if not in_shard(item.nodeid, index, count)]
    if deselected:
        config.hook.pytest_deselected(items=deselected)
        items[:] = [item for item in items if in_shard(item.nodeid, index, count)]


def pytest_configure(config: pytest.Config) -> None:
    """Check `--shard`, load the registry of identifiers.org and create the toy models."""
    text = config.getoption("--shard")
    if text is not None:
        try:
            parse_shard(text)
        except ValueError as err:
            raise pytest.UsageError(str(err)) from None
    if hasattr(config, "workerinput"):
        # a worker of pytest-xdist, the controller has prepared the session
        return
    # offline without a cache the tests which resolve annotations fail with
    # this error themselves, the others still run
    with suppress(WebserviceError):
        get_registry()
    create_toy_models(MODEL_BASE_PATH)
    create_sexes_models(SEXES_BASE_PATH)


_HISTORY_DATE = re.compile(rb"<dcterms:W3CDTF>[^<]*</dcterms:W3CDTF>")
"""A date of the model history, which changes with every build."""


class _FlatModelCache:
    """The flat models of the session by the content of their comp model.

    Only a flat model which was validated is cached, so a comp model with
    errors fails in every test.
    """

    def __init__(
        self, folder: Path, flat_model: "Callable[[Path, Path], None]"
    ) -> None:
        self.folder = folder
        self.flat_model = flat_model
        #: off in a test marked `uncached_build`
        self.enabled = True

    def __call__(self, comp_path: Path, flat_path: Path) -> None:
        """Write the flat model of a comp model, from the cache if it is there."""
        if not self.enabled:
            self.flat_model(comp_path, flat_path)
            return
        content = _HISTORY_DATE.sub(b"", comp_path.read_bytes())
        entry = self.folder / f"{hashlib.sha256(content).hexdigest()}.xml"
        if entry.is_file():
            shutil.copyfile(entry, flat_path)
            return
        self.flat_model(comp_path, flat_path)
        # written in full before it is visible to the other workers
        partial = entry.with_suffix(f".{os.getpid()}.partial")
        shutil.copyfile(flat_path, partial)
        partial.replace(entry)


@pytest.fixture(scope="session", autouse=True)
def flat_model_cache(
    pytestconfig: pytest.Config, tmp_path_factory: pytest.TempPathFactory
) -> "Iterator[_FlatModelCache | None]":
    """Take the flat models of the PBPK builds from the cache of the session.

    Session-scoped, so the builds of the fixtures of a module are cached as
    well. The test modules which build have imported the builder when they
    were collected, before the first test; without them nothing is built.
    """
    builder = sys.modules.get("pkpdsim.pbpk.builder")
    if builder is None:
        yield None
        return
    base = tmp_path_factory.getbasetemp()
    if hasattr(pytestconfig, "workerinput"):
        # the folder of a worker of pytest-xdist lies in the one of the session
        base = base.parent
    folder = base / "flat_models"
    folder.mkdir(exist_ok=True)
    cache = _FlatModelCache(folder, builder.flat_model)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(builder, "flat_model", cache)
        yield cache


@pytest.fixture(autouse=True)
def uncached_build(
    request: pytest.FixtureRequest,
    monkeypatch: pytest.MonkeyPatch,
    flat_model_cache: "_FlatModelCache | None",
) -> None:
    """Validate and flatten every build of a test marked `uncached_build`."""
    if flat_model_cache is not None and request.node.get_closest_marker(
        "uncached_build"
    ):
        monkeypatch.setattr(flat_model_cache, "enabled", False)


@pytest.fixture
def no_pyplot(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make `matplotlib.pyplot.figure` raise, so a figure on pyplot fails the test.

    A pyplot figure belongs to the interactive backend (TkAgg with a display),
    whose objects abort the process when they are collected next to a pool of
    processes (`main thread is not in main loop`, `Tcl_AsyncDelete`).
    """
    import matplotlib.pyplot as plt

    def refuse(*args: object, **kwargs: object) -> None:
        raise AssertionError("A figure must be a matplotlib.figure.Figure, not pyplot.")

    monkeypatch.setattr(plt, "figure", refuse)


@pytest.fixture(scope="module")
def toy_fit(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Fit the toy model with its profile likelihood and Fisher information.

    Two runs, so that the report of the fit has the waterfall of the runs.
    """
    output = tmp_path_factory.mktemp("fit")
    runs = fit_command(
        MODEL,
        [
            "--runs",
            "2",
            "--identifiability",
            "--max-points",
            "5",
            "--no-mapping-figures",
            "--output_dir",
            str(output),
        ],
    )
    (opid,) = runs
    return output / opid


@pytest.fixture(scope="module")
def fits_dir(toy_fit: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Get a folder of fits with the adopted toy fit."""
    directory = tmp_path_factory.mktemp("fits")
    adopt_fit(toy_fit, directory)
    return directory


@pytest.fixture
def git_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """Run git without the configuration of the user, see `release_git.GIT_ENV`."""
    for key, value in GIT_ENV.items():
        monkeypatch.setenv(key, value)


@pytest.fixture
def source_repo_model(tmp_path: Path, git_env: None) -> PKPDModel:
    """Get the toy model in a committed copy of its package in a monorepo."""
    root = tmp_path / "monorepo"
    init_repo(root)
    (root / "pyproject.toml").write_text(
        '[project]\nname = "toy"\nrequires-python = ">=3.14"\n'
    )
    shutil.copy(RELEASE_DATA / "uv.lock", root / "uv.lock")
    package = root / "pkg" / "pkpdsim_toy"
    shutil.copytree(
        MODEL.path,
        package,
        ignore=shutil.ignore_patterns("__pycache__", "models", "results"),
    )
    (package / "TODO.md").write_text("todo\n")
    (package / "literature").mkdir()
    (package / "literature" / "x.txt").write_text("x\n")
    git(root, "add", "-A")
    git(root, "commit", "-m", "init")
    return dataclasses.replace(
        MODEL, path=package, model_path=package / "models" / "toy.xml"
    )


@pytest.fixture
def remote_and_clone(tmp_path: Path, git_env: None) -> tuple[Path, Path]:
    """Get a bare remote and a clone of the release repository, released 0.1.0."""
    remote = tmp_path / "remote.git"
    remote.mkdir()
    git(remote, "init", "--bare", "-b", "main")
    clone = tmp_path / "toy-model"
    git(tmp_path, "clone", f"file://{remote}", str(clone))
    git(clone, "checkout", "-B", "main")
    (clone / "README.md").write_text("readme\n")
    git(clone, "add", "-A")
    git(clone, "commit", "-m", "initial")
    git(clone, "tag", "0.1.0")
    git(clone, "push", "origin", "main", "--tags")
    return remote, clone


@pytest.fixture
def release_settings(
    monkeypatch: pytest.MonkeyPatch, remote_and_clone: tuple[Path, Path]
) -> ReleaseSettings:
    """Get the release settings of the toy model, its repository is the remote."""
    remote, _ = remote_and_clone
    monkeypatch.setattr(
        source_module, "clone_url", lambda repository: f"file://{remote}"
    )
    assert MODEL.release is not None
    return MODEL.release
