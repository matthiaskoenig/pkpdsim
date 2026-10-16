"""The command `pkpd_release`: the complete release of a model repository.

```bash
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
pkpd_release pkdb_models.models.hctz --version 0.8.0
pkpd_release pkdb_models.models.hctz --repo ../hctz-model --version 0.8.0
```

The command checks the source and the release repository, builds the release in
the release repository (export of the model package, project files, models,
simulations, fit reports, release notes, docs, manifest) and publishes it
(commit, tag, push, GitHub Pages, Docker image, GitHub release). `--no-publish`
stops before the commit, so the build can be reviewed; a later run without the
flag publishes the reviewed build without building it again. A run which failed
while publishing is finished by running the same command again. See `release`
for the order of the steps.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from pathlib import Path

from sbmlutils.console import console

from pkpdsim import __version__
from pkpdsim.cli import ModelNotFoundError, resolve_model
from pkpdsim.fit_report import fit_record_dir
from pkpdsim.model import PKPDModel
from pkpdsim.release.build import build
from pkpdsim.release.docs import write_docs
from pkpdsim.release.manifest import changed_by_build, matches, prune, write_manifest
from pkpdsim.release.notes import fits_prefix, write_release_notes
from pkpdsim.release.project import project_context, write_project_files
from pkpdsim.release.publish import check_lock, check_tools
from pkpdsim.release.publish import publish as publish_release
from pkpdsim.release.settings import ReleaseSettings
from pkpdsim.release.source import (
    ReleaseError,
    ReleaseRepo,
    Runner,
    SourceRepo,
    check_version,
    export_package,
    release_repo,
    run,
    source_repo,
)

Publisher = Callable[[ReleaseSettings, Path, str, Path, Runner], list[str]]
"""Publishes a built release, see `pkpdsim.release.publish.publish`."""


def _settings(model: PKPDModel) -> ReleaseSettings:
    """Get the complete release settings of a model."""
    settings = model.release
    if settings is None:
        raise ReleaseError(f"MODEL.release is not set for the model '{model.name}'.")
    missing = settings.missing()
    if missing:
        raise ReleaseError(
            f"MODEL.release of the model '{model.name}' lacks {', '.join(missing)}."
        )
    return settings


def _archive_name(model: PKPDModel) -> str:
    """Get the name of the COMBINE archive which the GitHub release attaches."""
    if model.archive_path is None:
        raise ReleaseError(
            f"The model '{model.name}' has no archive_path, the GitHub release "
            "attaches its COMBINE archive."
        )
    return model.archive_path.name


def _check_fits(model: PKPDModel) -> None:
    """Check that every adopted fit has its record in the model package."""
    if model.fits:
        fits_prefix(model)
    for fit_id in model.fits:
        try:
            fit_record_dir(model, fit_id)
        except (ValueError, FileNotFoundError) as err:
            raise ReleaseError(str(err)) from err


def _tagged_at_head(path: Path, version: str, runner: Runner) -> bool:
    """Tell whether the tag of the version exists and is the current commit."""
    if not runner(["git", "tag", "--list", version], path).strip():
        return False
    tagged = runner(["git", "rev-list", "-n", "1", version], path).strip()
    return tagged == runner(["git", "rev-parse", "HEAD"], path).strip()


def _lock_publishable(path: Path) -> bool:
    try:
        check_lock(path)
    except ReleaseError:
        return False
    return True


def _build(
    model: PKPDModel,
    settings: ReleaseSettings,
    source: SourceRepo,
    target: ReleaseRepo,
    version: str,
    pkpdsim_path: Path | None,
    workers: str,
    runner: Runner,
) -> None:
    """Build the release in the release repository and write its manifest."""
    path = target.path
    context = project_context(model, settings, source, version, pkpdsim_path)
    console.print(f"Export of {model.package} into {path}")
    written: list[Path] = export_package(model, settings, source, path, runner)
    written += write_project_files(context, path)
    console.print("Build of the models, simulations and fit reports")
    written += build(model, path, runner, workers)
    console.print("Release notes and docs")
    written.append(write_release_notes(context, target, model, source, runner))
    written += write_docs(context, path)
    prune(path, written, runner)
    write_manifest(path, version, source, written)


def release(
    model: PKPDModel,
    *,
    version: str,
    repo: Path | None,
    publish: bool,
    pkpdsim_path: Path | None,
    workers: str,
    runner: Runner = run,
    publisher: Publisher = publish_release,
) -> Path:
    """Make the complete release of a model repository.

    The steps, in order; the first failure stops the release:

    1. Checks, before anything is changed: the release settings of the model
       are complete, the source repository is clean (`source_repo`), the
       release repository is a checkout of the repository on `main` and up to
       date (`release_repo`, cloned when missing), the version is greater than
       the last release (`check_version`), every adopted fit has its record;
       when publishing, the model declares its COMBINE archive and docker and
       `gh` are ready (`check_tools`), so a missing tool fails before the long
       build.
    2. The build, skipped if a publishing run finds a build of this release
       (`matches`, i.e. a reviewed `--no-publish` build) with a publishable lock.
       Otherwise the release repository must be clean or changed only by an
       earlier build of this version (`changed_by_build`). The stages write, in
       this order: the export of the package, the project files, the build
       (`uv.lock`, models, simulations, fit reports), the release notes and the
       docs. `prune` then removes every file which no stage wrote, and the
       manifest `release.json` is written last, so it hashes the final tree
       which `matches` compares.
    3. `--no-publish` stops here, before the commit.
    4. The lock is checked (`check_lock`, it exists only after the build), then
       the release is published (`publisher`).

    If the tag of the version is at `HEAD` and the manifest matches, the release
    was committed by an earlier run which failed while publishing: the version
    check and the build are skipped and the publishing steps resume.

    Args:
        model: the model which is released, with its release settings.
        version: the version of the release, `X.Y.Z`.
        repo: the checkout of the release repository, default is a sibling of
            the source repository named like the release repository.
        publish: publish the release; `False` stops before the commit.
        pkpdsim_path: a local checkout of pkpdsim instead of the released
            pkpdsim, only without `publish`.
        workers: the number of workers of the simulations, `auto` or a number.
        runner: runs git, uv, the commands of the model, docker and gh.
        publisher: publishes the built release.

    Returns:
        The checkout of the release repository.

    Raises:
        ReleaseError: if a check or a step fails; the message is for the user.
    """
    if publish and pkpdsim_path is not None:
        raise ReleaseError("a local pkpdsim (pkpdsim_path) is only for --no-publish")
    settings = _settings(model)
    if publish:
        _archive_name(model)
    source = source_repo(model, runner)
    target = release_repo(settings, repo, source, runner)
    path = target.path
    resume = _tagged_at_head(path, version, runner) and matches(
        path, version, source.commit, __version__, runner
    )
    if not resume:
        check_version(version, target.last_tag)
    _check_fits(model)
    if publish:
        check_tools(path, runner)

    built = resume or (
        publish
        and matches(path, version, source.commit, __version__, runner)
        and _lock_publishable(path)
    )
    if not built:
        if not changed_by_build(path, version, runner):
            raise ReleaseError(
                "the release repository has changes which no release made, "
                f"discard them first (see git -C {path} status)"
            )
        _build(model, settings, source, target, version, pkpdsim_path, workers, runner)
    if not publish:
        return path

    check_lock(path)
    console.print(f"Publishing {settings.repository} {version}")
    publisher(settings, path, version, path / _archive_name(model), runner)
    return path


def parser() -> argparse.ArgumentParser:
    """Get the parser of the command `pkpd_release`."""
    p = argparse.ArgumentParser(
        prog="pkpd_release",
        description="Make the complete release of the repository of a model: "
        "build the models, simulations, fit reports, release notes and docs, "
        "commit, tag, push, push the Docker image and create the GitHub release.",
    )
    p.add_argument(
        "package",
        metavar="PACKAGE",
        help="the dotted name of the model package, e.g. pkdb_models.models.hctz",
    )
    p.add_argument("--version", required=True, help="the version of the release, X.Y.Z")
    p.add_argument(
        "--repo",
        type=Path,
        metavar="PATH",
        default=None,
        help="the checkout of the release repository; default is a sibling of "
        "the source repository named like the release repository, cloned when "
        "missing",
    )
    p.add_argument(
        "--no-publish",
        dest="publish",
        action="store_false",
        help="build the release but stop before the commit, for a review; a "
        "later run without the flag publishes the reviewed build",
    )
    p.add_argument(
        "--pkpdsim-path",
        type=Path,
        metavar="PATH",
        default=None,
        help="a local checkout of pkpdsim instead of the released one, only "
        "with --no-publish",
    )
    p.add_argument(
        "--workers",
        default="auto",
        help="the number of workers of the simulations (default: auto)",
    )
    return p


def main(argv: Sequence[str] | None = None) -> None:
    """The console script `pkpd_release`, see `release`.

    Args:
        argv: the command line arguments, `sys.argv[1:]` by default.

    Raises:
        SystemExit: with 2 for a usage error, with 1 if the release fails.
    """
    p = parser()
    options = p.parse_args(argv)
    if options.pkpdsim_path is not None and options.publish:
        p.error("--pkpdsim-path is only allowed with --no-publish")
    try:
        model = resolve_model(options.package)
    except ModelNotFoundError as err:
        p.error(str(err))
    try:
        path = release(
            model,
            version=options.version,
            repo=options.repo.resolve() if options.repo else None,
            publish=options.publish,
            pkpdsim_path=(
                options.pkpdsim_path.resolve() if options.pkpdsim_path else None
            ),
            workers=options.workers,
        )
    except ReleaseError as err:
        console.print(str(err), style="error", markup=False)
        raise SystemExit(1) from err
    if options.publish:
        console.print(f"Released {options.version} from {path}.")
    else:
        console.print(f"{path}\nreview, then run without --no-publish", markup=False)
