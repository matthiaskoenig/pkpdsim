"""The command `pkpd_release`: the complete release of a model repository.

```bash
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
pkpd_release pkdb_models.models.hctz --version 0.8.0
pkpd_release pkdb_models.models.hctz --repo ../hctz-model --version 0.8.0
```

The release is sbmlsim's (`sbmlsim.release`) with the release settings and the
release hook of pkpdsim (`pkpdsim.model.PKPDReleaseSettings`,
`pkpdsim.release.hook`); the command runs its two steps:

1. The build, the kind `release` of the run layer of sbmlsim
   (`sbmlsim.run.release_kinds`, `python -m sbmlsim.release build`), into an
   output (`--out`, by default `<repository>-release-<version>` next to the
   checkout of the release repository): the release repository is cloned into
   `<out>/repository`, from the checkout (`--repo`, by default the sibling of
   the source repository named like the release repository) when it exists,
   which is only read, else from GitHub, and the release is built there (export
   of the model package, project files, models, simulations, fit reports,
   credibility, release notes, docs, manifest).
2. The publishing (`python -m sbmlsim.release publish`): commit, tag, GitHub
   Pages, push, Docker image and GitHub release from `<out>/repository`;
   afterwards the checkout is fast-forwarded to the release when it is clean
   and on `main`.

`--no-publish` stops before the commit, so the build can be reviewed; a later
run without the flag publishes the reviewed build without building it again. A
run which failed while publishing is finished by running the same command
again. See `release` for the order of the steps.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Final

from sbmlsim.release.publish import check_tools
from sbmlsim.release.publish import publish as publish_steps
from sbmlsim.release.source import (
    ReleaseError,
    Runner,
    SourceRepo,
    parse_version,
    run,
    source_repo,
)
from sbmlsim.release.workflow import (
    REPOSITORY_DIR,
    Publisher,
    execute_release,
    publish_release,
    publishable,
    release_settings,
)
from sbmlsim.run.context import RUN_FILE, spec_hash
from sbmlsim.run.kinds import ConfigurationError
from sbmlsim.run.release_kinds import RELEASE, ReleaseSpec
from sbmlutils.console import console

from pkpdsim.cli import ModelNotFoundError, project_reference, resolve_model

if TYPE_CHECKING:
    from sbmlsim.release.settings import ReleaseSettings

    from pkpdsim.model import PKPDModel

__all__ = ["LOCAL_SUFFIX", "OUT_SUFFIX", "default_out", "main", "parser", "release"]

logger = logging.getLogger(__name__)

#: the suffix of the default output of a build, `<repository>-release-<version>`
OUT_SUFFIX: Final = "-release-"

#: the suffix of the default output of a build with a local pkpdsim, which is
#: never published, `<repository>-release-<version>-local`
LOCAL_SUFFIX: Final = "-local"

#: the refusals of the clone in the out which a change of the checkout or a new
#: out resolves, with what to do in the checkout
_CHECKOUT_REFUSALS: Final = {
    "behind origin/main": "pull in the checkout",
    "not 'main'": "pull in the checkout",
    "ahead of origin/main": "push or reset in the checkout",
}


def _checkout(settings: ReleaseSettings, repo: Path | None, src: SourceRepo) -> Path:
    """Get the checkout of the release repository, given or the default sibling."""
    return repo if repo is not None else src.root.parent / settings.repo_name


def default_out(checkout: Path, version: str, *, local: bool = False) -> Path:
    """Get the default output of the build of a release.

    Args:
        checkout: the checkout of the release repository, which need not exist.
        version: the version of the release.
        local: a build with a local pkpdsim (`--pkpdsim-path`), which is never
            published and gets an output of its own.

    Returns:
        `<checkout>-release-<version>` next to the checkout, with `-local` for
        a build with a local pkpdsim.
    """
    name = f"{checkout.name}{OUT_SUFFIX}{version}{LOCAL_SUFFIX if local else ''}"
    return checkout.parent / name


def _hint(err: ReleaseError, checkout: Path, out: Path) -> ReleaseError:
    """Add what to do to a refusal of the clone in the out.

    The clone of the checkout in the out is refused when it is behind
    `origin/main`, on another branch or ahead of it; the checkout is changed
    and the out removed, so the next build clones it again.
    """
    message = str(err)
    for refusal, action in _CHECKOUT_REFUSALS.items():
        if refusal in message and "remove" not in message:
            return ReleaseError(f"{message}; {action} {checkout}, or remove {out}")
    return err


def _check_out(out: Path, spec: dict[str, object]) -> None:
    """Refuse an out which holds another build than the one of the spec.

    An earlier build of the same spec is built again; the build of another
    version, another checkout or a local pkpdsim is no run the build resumes.

    Raises:
        ReleaseError: for an out which is not empty and holds no build of the
            spec, naming what to do.
    """
    if not out.is_dir() or not any(out.iterdir()):
        return
    try:
        record = json.loads((out / RUN_FILE).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        record = None
    given = spec_hash(ReleaseSpec.model_validate(spec), RELEASE.inputs)
    if not isinstance(record, dict) or (
        record.get("kind"),
        record.get("spec_hash"),
    ) != (RELEASE.name, given):
        raise ReleaseError(
            f"{out} holds another build than the one of this release (another "
            "version, --repo or --pkpdsim-path, or no build): remove "
            f"{out} or pass --out"
        )


def _build(
    reference: str,
    version: str,
    out: Path,
    checkout: Path,
    pkpdsim_path: Path | None,
    workers: str,
    runner: Runner,
) -> None:
    """Build the release into the output with the kind `release` of sbmlsim.

    The build resumes in an output of an earlier build of the same spec; an
    output of another build is refused (`_check_out`).

    Raises:
        ReleaseError: if the release is refused or a stage of the build fails.
    """
    repository = checkout.resolve() if checkout.exists() else None
    if repository is None and not (out / REPOSITORY_DIR).exists():
        logger.info(
            "The checkout %s does not exist, the build clones the release "
            "repository from GitHub.",
            checkout,
        )
    spec: dict[str, object] = {"project": reference, "version": version}
    if repository is not None:
        spec["repository"] = str(repository)
    if pkpdsim_path is not None:
        spec["framework_path"] = str(pkpdsim_path.resolve())
    _check_out(out, spec)
    try:
        execute_release(
            reference,
            version,
            out,
            repository=repository,
            framework_path=pkpdsim_path,
            workers=None if workers == "auto" else int(workers),
            resume=True,
            runner=runner,
        )
    except ConfigurationError as err:
        raise _hint(ReleaseError(str(err)), checkout, out) from err
    except ReleaseError as err:
        raise _hint(err, checkout, out) from err


def _fast_forward(checkout: Path, runner: Runner) -> None:
    """Fast-forward a clean checkout on `main` to the published release.

    A checkout with changes or on another branch is left as it is, with a
    warning.
    """
    if not (checkout / ".git").exists():
        return
    try:
        branch = runner(["git", "branch", "--show-current"], checkout).strip()
        clean = not runner(["git", "status", "--porcelain"], checkout).strip()
        if branch != "main" or not clean:
            logger.warning(
                "%s is not a clean checkout of main, pull the release into it.",
                checkout,
            )
            return
        runner(["git", "fetch", "--tags", "origin"], checkout)
        runner(["git", "merge", "--ff-only", "origin/main"], checkout)
    except ReleaseError as err:
        logger.warning("%s is not updated, pull the release into it: %s", checkout, err)


def release(
    model: PKPDModel,
    *,
    version: str,
    repo: Path | None,
    publish: bool,
    pkpdsim_path: Path | None,
    workers: str,
    out: Path | None = None,
    runner: Runner = run,
    publisher: Publisher = publish_steps,
) -> Path:
    """Make the complete release of a model repository.

    The steps, in order; the first failure stops the release:

    1. Checks, before anything is built: the version is `X.Y.Z`, a local
       pkpdsim only without publishing, the release settings of the model are
       complete, the model is the model object of its package (the kind loads
       it by its reference, `pkpdsim.cli.project_reference`), the source
       repository is clean (`source_repo`); when publishing, the model declares
       its COMBINE archive and docker, its buildx plugin and `gh` are ready
       (`check_tools`), so a missing tool fails before the long build.
    2. The build, the kind `release` of sbmlsim into `out`, which clones the
       release repository into `out/repository` (from the checkout `repo`
       when it exists, else from GitHub) and checks it, the version and the
       records of the adopted fits, see
       `sbmlsim.release.workflow.execute_release`. It is skipped if a
       publishing run finds a publishable build of this release in the output
       (`sbmlsim.release.workflow.publishable`: a reviewed `--no-publish`
       build, or the release commit of a publishing which failed half way);
       an earlier build of the same spec in the output is built again (a
       resume of the run layer).
    3. `--no-publish` stops here, before the commit.
    4. The publishing of the build of the output
       (`sbmlsim.release.workflow.publish_release`, `publisher`), then the
       checkout is fast-forwarded to the release when it is clean and on
       `main`.

    Args:
        model: the model which is released, with its release settings.
        version: the version of the release, `X.Y.Z`.
        repo: the checkout of the release repository, default is a sibling of
            the source repository named like the release repository; only
            read, and fast-forwarded after publishing.
        publish: publish the release; `False` stops before the commit.
        pkpdsim_path: a local checkout of pkpdsim instead of the released
            pkpdsim, only without `publish`.
        workers: the number of workers of the simulations, `auto` (the CPUs of
            the quota) or a number.
        out: the output of the build, `default_out` of the checkout by default
            (an output of its own for a build with a local pkpdsim).
        runner: runs git, uv, the analyses of the model, docker and gh.
        publisher: publishes the built release.

    Returns:
        The release repository of the build, `out/repository`.

    Raises:
        ReleaseError: if a check or a step fails; the message is for the user.
    """
    parse_version(version)
    if publish and pkpdsim_path is not None:
        raise ReleaseError("a local pkpdsim (pkpdsim_path) is only for --no-publish")
    settings = release_settings(model)
    if publish and settings.archive is None:
        raise ReleaseError(
            f"The model '{model.name}' has no archive_path, the GitHub release "
            "attaches its COMBINE archive."
        )
    try:
        reference = project_reference(model)
    except ValueError as err:
        raise ReleaseError(str(err)) from err
    src = source_repo(model, runner)
    checkout = _checkout(settings, repo, src)
    if out is None:
        out = default_out(checkout, version, local=pkpdsim_path is not None)
    path = out / REPOSITORY_DIR
    if publish:
        check_tools(model.root, runner)
    if not (publish and publishable(model, version, out, runner=runner) is None):
        _build(reference, version, out, checkout, pkpdsim_path, workers, runner)
    if not publish:
        return path
    try:
        publish_release(model, version, out, runner=runner, publisher=publisher)
    except ReleaseError as err:
        raise _hint(err, checkout, out) from err
    _fast_forward(checkout, runner)
    return path


def _workers(text: str) -> str:
    """Check the number of workers: `auto` or a positive int."""
    if text != "auto" and not (text.isdigit() and int(text) > 0):
        raise argparse.ArgumentTypeError(f"'{text}' is not auto or a positive int")
    return text


def parser() -> argparse.ArgumentParser:
    """Get the parser of the command `pkpd_release`."""
    p = argparse.ArgumentParser(
        prog="pkpd_release",
        description="Make the complete release of the repository of a model: "
        "build the models, simulations, fit reports, release notes and docs into "
        "a clone of the release repository, commit, tag, push, push the Docker "
        "image and create the GitHub release.",
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
        help="the checkout of the release repository, which the build clones and "
        "publishing fast-forwards; default is a sibling of the source repository "
        "named like the release repository, the repository on GitHub is cloned "
        "when it is missing",
    )
    p.add_argument(
        "--out",
        type=Path,
        metavar="DIR",
        default=None,
        help="the output of the build, which holds the release repository in "
        "DIR/repository; default is <repo>-release-<version> next to the checkout, "
        "<repo>-release-<version>-local with --pkpdsim-path",
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
        type=_workers,
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
            out=options.out.resolve() if options.out else None,
        )
    except ReleaseError as err:
        print(f"pkpd_release: error: {err}", file=sys.stderr)
        raise SystemExit(1) from err
    if options.publish:
        _print_verbatim(f"Released {options.version} from {path}.")
    else:
        _print_verbatim(f"{path}\nreview, then run without --no-publish")


def _print_verbatim(text: str) -> None:
    """Print text with paths, which are copied from the terminal.

    The lines are not wrapped at the terminal width, and their characters are no
    markup and not highlighted.
    """
    console.print(text, soft_wrap=True, markup=False, highlight=False)
