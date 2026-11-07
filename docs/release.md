# Release of a model repository

Every model of pkdb_models is published as a repository of its own (`hctz-model`, `glimepiride-model`, ...), archived by Zenodo and cited by its DOI. `pkpd_release` makes the complete release of such a repository from the model package, the same way for every model:

```bash
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
pkpd_release pkdb_models.models.hctz --version 0.8.0
```

The first command builds the release into a clone of the release repository and stops before the commit, so the build can be looked at; the second one publishes it. Once a model has its release metadata and the one-time setup is done, a release needs no manual step before, during or after the command.

The release is the one of sbmlsim (`sbmlsim.release`): the build is the kind `release` of its run layer and the publishing the step `python -m sbmlsim.release publish`, see the page of the run layer of sbmlsim. `pkpd_release` runs both with the release settings and the release hook of pkpdsim, which add what is pharmacokinetic (see [What pkpdsim adds](#what-pkpdsim-adds)).

The release repository holds only what the build generates. The models, the simulations and the fit reports are created in the release repository from the exported package and its released dependencies, so a release proves that the model runs on its own. Nothing in a release repository is edited by hand: a change is made in the model package and released again.

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

- `title`: the title of the release: the Zenodo record, the citation, the heading of the README and the name of the docs site.
- `repository`: the GitHub repository as `owner/name`; the docs site is `https://<owner>.github.io/<name>/`.
- `doi`: the concept DOI of Zenodo, which cites every version.
- `authors`: the authors of the citation in order, each with the name as `"Family, Given"`, the ORCID and the affiliation.
- `keywords`, `funding`, `description`: the Zenodo record (`.zenodo.json`), the docs and the README.
- `docker_image`: the image which is pushed as `<docker_image>:<version>` and `:latest`. The image of the version is built and pushed while the registry lacks it, a pushed version is never rebuilt and `:latest` is then tagged on the registry from it.
- `exclude`: paths of the model folder which are not released, by default `TODO.md`, `literature/` and `thesis/`.
- `first_year`: the year of the first release, the start of the copyright.

The model object sets it as `PKPDModel(..., release=release)`. A model without `release`, with an empty field or a `first_year` which is no positive year cannot be released; the command names the missing fields. An empty `exclude` releases the whole model folder. The adopted fits of `MODEL.fits` are released with their records, the archive of `MODEL.archive_path` is attached to the GitHub release; it lies inside the model package.

The defaults of pkpdsim, which a model package does not write; a field it gives keeps its value (e.g. `framework="sbmlsim"`), and the release settings of sbmlsim themselves are refused by the model object, whose release is built with pkpdsim:

- `framework`: `pkpdsim`, the package the release is built with; the lock of the release must contain it, `release.json` records its version and the release notes say "with pkpdsim <version>".
- `core_libraries`: `pkpdsim`, `sbmlsim` and `sbmlutils`, whose release notes the release notes summarize; `library_owner`: `matthiaskoenig`, the owner of their repositories.
- `hook`: `pkpdsim.release.hook:HOOK`, the release hook of pkpdsim.
- `archive`: set by the model object, `archive_path` relative to the root of the model.

## One-time setup

Once per machine:

- `docker login` with an account which may push the `docker_image`; the buildx plugin of docker (part of Docker Desktop and Docker Engine) tags `:latest` on the registry.
- `gh auth login` with an account which may push to the repository and create its releases.

Once per release repository:

- The repository exists on GitHub (the command does not create it), it may be empty; the command clones it when there is no checkout.
- The GitHub integration of Zenodo is switched on for the repository (Zenodo, account settings, GitHub), so that Zenodo archives every GitHub release under the concept DOI.

A release depends on released packages only: the version of pkpdsim which the monorepo locks must be on PyPI.

## The command

```bash
pkpd_release PACKAGE --version VERSION [--repo PATH] [--out DIR] [--no-publish] [--pkpdsim-path PATH] [--workers N]
```

- `PACKAGE`: the dotted name of the model package, e.g. `pkdb_models.models.hctz`.
- `--version`: the version of the release, `X.Y.Z`, greater than the last tag.
- `--repo`: the checkout of the release repository; by default a sibling of the monorepo named like the repository (`../hctz-model`). The build clones it and reads it only, a missing checkout is not created and the build clones the repository on GitHub instead (the command says so); a clone which is behind `origin/main`, on another branch or ahead of it is refused with what to do in the checkout, or to remove the output; after publishing a clean checkout on `main` is fast-forwarded to the release.
- `--out`: the output of the build, which holds the release repository in `DIR/repository` and the records of the run of the kind (`run.json`, `spec.json`, `summary.json`, `progress.jsonl`); by default `<repo>-release-<version>` next to the checkout (`../hctz-model-release-0.8.0`), and `<repo>-release-<version>-local` for a build with `--pkpdsim-path`, which is never published. A run builds again into the output of an earlier run of the same release; an output of another build (another version, `--repo` or `--pkpdsim-path`) is refused: remove it or pass `--out`. The outputs of earlier releases stay next to the checkout until they are removed by hand; once a version is published its output is no longer needed.
- `--no-publish`: build the release and stop before the commit.
- `--pkpdsim-path`: a local checkout of pkpdsim instead of the released one, for a test of a pkpdsim which is not released yet; only with `--no-publish`.
- `--workers`: the number of workers of the simulations, `auto` (the CPUs of the quota) by default.

The steps, in order; the first failure stops the release and ends the command with exit code 1 and the reason:

1. **Checks** of `pkpd_release`, before anything is built: the version, the release metadata, the monorepo is clean, and for a run which publishes the archive of the model, docker with its buildx plugin and `gh`.
2. **Build**, the kind `release` of sbmlsim into the output: the release repository is cloned into `DIR/repository` and checked (on `main`, up to date, without local commits and without changes which no release made, the version greater than the last tag, every adopted fit with its record), then the export of the model package, the project files, the analyses of the run layer in the release repository (the factory, the simulations, with distributions of the fitted parameters the simulations with the confidence interval, the reports of the adopted fits, the assessment of the credibility of a model which declares its contexts of use), the release notes, the docs and the README, the prune and the manifest `release.json`. A build of the same spec in the output is built again, so a failed or interrupted build is finished by the next run.
3. `--no-publish` stops here.
4. **Publish** the build of the output, the step `python -m sbmlsim.release publish`: the lock must hold packages of PyPI only; then the Docker image, the commit, the tag, GitHub Pages, the push, the push of the Docker image and the GitHub release with the release notes and the COMBINE archive. Then the checkout of `--repo` is fast-forwarded to the release, a checkout with changes or on another branch is left as it is with a warning to pull it.

## Review, then publish

A run without `--no-publish` which finds the reviewed build of the same release in the output (same version, same commit of the monorepo, the same pkpdsim locked by the monorepo, every file as it was built) publishes it as it is, without running the build again. A file changed during the review makes the run build again, so the published release is always the output of the build.

A run which failed while publishing is finished by running the same command again: every publishing step checks whether it is done already (the commit, the tag, the push, the image, the GitHub release) and is skipped then. If the last commit of the release repository of the output is the release of the version (the tag is missing or at that commit), the build is skipped and the publishing steps resume, even if the monorepo moved on since.

## What pkpdsim adds

The release hook of pkpdsim (`pkpdsim.release.hook.PKReleaseHook`) adds the content of the pharmacokinetic models to the release of sbmlsim:

- The page `credibility` of the docs shows the evidence of the extensions of pkpdsim (see [Context of use and credibility](credibility.md)): after the acceptance table of a context of use, whose last row is the rule of the ratio pairs (`(ratio pairs)`), the table of the ratio pairs and the temporal validation; between the parameter uncertainty and the variability of a context the inputs of its decisive mappings; after the sections of the contexts the study data of the adopted fits with the records of their collection from pkdb_data.
- The release notes name a study of the data `<substance>/<study>`, a folder `data/<substance>/<study>` of the model package, as the releases of pkpdsim 0.19 named them, so the first release after the release moved into sbmlsim renames no study.
