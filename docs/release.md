# Release of a model repository

Every model of pkdb_models is published as a repository of its own (`hctz-model`, `glimepiride-model`, ...), archived by Zenodo and cited by its DOI. `pkpd_release` makes the complete release of such a repository from the model package, the same way for every model:

```bash
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
pkpd_release pkdb_models.models.hctz --version 0.8.0
```

The first command builds the release into a clone of the release repository and stops before the commit, so the build can be looked at; the second one publishes it.
The release is the one of sbmlsim (`sbmlsim.release`): the build is the kind `release` of its run layer and the publishing the step `python -m sbmlsim.release publish`, see the page of the run layer of sbmlsim. `pkpd_release` runs both with the release settings and the release hook of pkpdsim, which add what is pharmacokinetic (see [What pkpdsim adds](#what-pkpdsim-adds)).

What a release repository holds and why nothing in it is edited by hand are described in [Release](https://matthiaskoenig.github.io/sbmlsim/release/) of sbmlsim.

## The release metadata of a model

What the monorepo does not know about a release is the `release` of the model object, a `PKPDReleaseSettings` (`pkpdsim.model`), the `ReleaseSettings` of sbmlsim (`sbmlsim.release.settings`) with the defaults of pkpdsim:

```python
from sbmlsim.release.settings import Author

from pkpdsim.model import PKPDReleaseSettings

release = PKPDReleaseSettings(
    title="Physiologically based pharmacokinetic model of hydrochlorothiazide",
    repository="matthiaskoenig/hctz-model",
    doi="10.5281/zenodo.15427418",
    authors=[
        Author(
            "König, Matthias",
            orcid="0000-0003-1725-179X",
            affiliation="Humboldt-University Berlin, Institute for Biology",
        ),
    ],
    keywords=["hydrochlorothiazide", "pharmacokinetics", "PBPK", "SBML"],
    docker_image="matthiaskoenig/hctz",
    funding="The funding of the work on the model.",
    description="A physiologically based pharmacokinetic model of hydrochlorothiazide.",
    first_year=2025,
)
print(release.site_url)
print(release.missing())
```

The fields (`title`, `repository`, `doi`, `authors`, `keywords`, `funding`, `description`, `docker_image`, `exclude`, `first_year`) and the checks of a release without them are those of [Release](https://matthiaskoenig.github.io/sbmlsim/release/) of sbmlsim.

The model object sets it as `PKPDModel(..., release=release)`. A model without `release`, with an empty field or a `first_year` which is no positive year cannot be released; the command names the missing fields. An empty `exclude` releases the whole model folder. The adopted fits of `MODEL.fits` are released with their records, the archive of `MODEL.archive_path` is attached to the GitHub release; it lies inside the model package.

The defaults of pkpdsim, which a model package does not write; a field it gives keeps its value (e.g. `framework="sbmlsim"`), and the release settings of sbmlsim themselves are refused by the model object, whose release is built with pkpdsim:

- `framework`: `pkpdsim`, the package the release is built with; the lock of the release must contain it, `release.json` records its version and the release notes say "with pkpdsim <version>".
- `core_libraries`: `pkpdsim`, `sbmlsim` and `sbmlutils`, whose release notes the release notes summarize; `library_owner`: `matthiaskoenig`, the owner of their repositories.
- `hook`: `pkpdsim.release.hook:HOOK`, the release hook of pkpdsim.
- `archive`: set by the model object, `archive_path` relative to the root of the model.

## One-time setup

The logins of docker and `gh`, the release repository and the Zenodo integration are the one-time setup of [Release](https://matthiaskoenig.github.io/sbmlsim/release/#one-time-setup) of sbmlsim. A release depends on released packages only: the version of pkpdsim which the monorepo locks must be on PyPI.

## The command

```bash
pkpd_release PACKAGE --version VERSION [--repo PATH] [--out DIR] [--no-publish] [--pkpdsim-path PATH] [--workers N]
```

- `PACKAGE`: the dotted name of the model package, e.g. `pkdb_models.models.hctz`.
- `--version`: the version of the release, `X.Y.Z`, greater than the last tag.
- `--repo`: the checkout of the release repository; by default a sibling of the monorepo named like the repository (`../hctz-model`). The build clones it and reads it only, a missing checkout is not created and the build clones the repository on GitHub instead (the command says so); commits of the checkout which the release repository lacks are refused (push them first), a clone in the output which is behind `origin/main` is fast-forwarded and one on another branch or with commits of its own is refused (remove the output or give another one); after publishing a clean checkout on `main` is fast-forwarded to the release.
- `--out`: the output of the build, which holds the release repository in `DIR/repository` and the records of the run of the kind (`run.json`, `spec.json`, `summary.json`, `progress.jsonl`); by default `<repo>-release-<version>` next to the checkout (`../hctz-model-release-0.8.0`), and `<repo>-release-<version>-local` for a build with `--pkpdsim-path`, which is never published. A run builds again into the output of an earlier run of the same release; an output of another build (another version, `--repo` or `--pkpdsim-path`) is refused: remove it or pass `--out`. The outputs of earlier releases stay next to the checkout until they are removed by hand; once a version is published its output is no longer needed.
- `--no-publish`: build the release and stop before the commit.
- `--pkpdsim-path`: a local checkout of pkpdsim instead of the released one, for a test of a pkpdsim which is not released yet; only with `--no-publish`.
- `--workers`: the number of workers of the simulations, `auto` (the CPUs of the quota) by default.

The command runs the checks, the build and the publishing of sbmlsim (the kind `release` and `python -m sbmlsim.release publish`), whose steps, the review of a build before it is published and the resumption of a publishing which failed ([Review, then publish](https://matthiaskoenig.github.io/sbmlsim/release/#review-then-publish)) are described in [Release](https://matthiaskoenig.github.io/sbmlsim/release/#the-build) of sbmlsim. The first failure stops the release and ends the command with exit code 1 and the reason. `pkpd_release` adds its own checks before anything is built: the version, the release metadata, a clean monorepo, and for a run which publishes the archive of the model, docker with its buildx plugin and `gh`.

## What pkpdsim adds

The release hook of pkpdsim (`pkpdsim.release.hook.PKReleaseHook`) adds the content of the pharmacokinetic models to the release of sbmlsim:

- The page `credibility` of the docs shows the evidence of the extensions of pkpdsim (see [Context of use and credibility](credibility.md)): after the acceptance table of a context of use, whose last row is the rule of the ratio pairs (`(ratio pairs)`), the table of the ratio pairs and the temporal validation; between the parameter uncertainty and the variability of a context the inputs of its decisive mappings; after the sections of the contexts the study data of the adopted fits with the records of their collection from pkdb_data.
- The release notes name a study of the data `<substance>/<study>`, a folder `data/<substance>/<study>` of the model package, as the releases of pkpdsim 0.19 named them, so the first release after the release moved into sbmlsim renames no study.
