# Release of a model repository

Every model of pkdb_models is published as a repository of its own (`hctz-model`, `glimepiride-model`, ...), archived by Zenodo and cited by its DOI. `pkpd_release` makes the complete release of such a repository from the model package, the same way for every model:

```bash
pkpd_release pkdb_models.models.hctz --version 0.8.0 --no-publish
pkpd_release pkdb_models.models.hctz --version 0.8.0
```

The first command builds the release in the release repository and stops before the commit, so the build can be looked at; the second one publishes it. Once a model has its release metadata and the one-time setup is done, a release needs no manual step before, during or after the command.

The release repository holds only what the command generates. The models, the simulations and the fit reports are created in the release repository from the exported package and its released dependencies, so a release proves that the model runs on its own. Nothing in a release repository is edited by hand: a change is made in the model package and released again.

## The release metadata of a model

What the monorepo does not know about a release is the `release` of the model object, a `ReleaseSettings` (`pkpdsim.release.settings`):

```python
from pkpdsim.release.settings import Author, ReleaseSettings

release = ReleaseSettings(
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

The model object sets it as `PKPDModel(..., release=release)`. A model without `release`, with an empty field or a `first_year` which is no positive year cannot be released; the command names the missing fields. An empty `exclude` releases the whole model folder. The adopted fits of `MODEL.fits` are released with their records, the archive of `MODEL.archive_path` is attached to the GitHub release.

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
pkpd_release PACKAGE --version VERSION [--repo PATH] [--no-publish] [--pkpdsim-path PATH] [--workers N]
```

- `PACKAGE`: the dotted name of the model package, e.g. `pkdb_models.models.hctz`.
- `--version`: the version of the release, `X.Y.Z`, greater than the last tag.
- `--repo`: the checkout of the release repository; by default a sibling of the monorepo named like the repository (`../hctz-model`), cloned when missing.
- `--no-publish`: build the release and stop before the commit.
- `--pkpdsim-path`: a local checkout of pkpdsim instead of the released one, for a test of a pkpdsim which is not released yet; only with `--no-publish`.
- `--workers`: the number of workers of the simulations, `auto` by default.

The steps, in order; the first failure stops the release. Every stage prints a line, and the long steps (`uv`, the simulations, the fit reports, `git push`, `docker build` and `docker push`) show their output as they run:

1. **Checks**, before anything is changed: the release metadata is complete, the monorepo is clean and committed (its commit is recorded in the release), the release repository is a checkout of the repository on `main`, up to date with the remote, without local commits (besides the release commit of a release which failed while publishing) and without changes which no release made, the version is greater than the last tag and every adopted fit has its record. A run which publishes also checks that the model has its archive, docker runs with its buildx plugin and `gh` is logged in.
2. **Export** of the files git tracks in the model folder into `src/<package as folders>/`, without the excluded paths.
3. **Project files** from the templates of pkpdsim: `pyproject.toml`, `Dockerfile`, the workflows, `zensical.toml`, `.zenodo.json`, `.gitignore`, `.python-version`. The Docker image and the workflows install the same uv (0.13.4, `pkpdsim.release.project.UV_VERSION`).
4. **Build** in the environment of the release repository (`uv lock`, locked to the library versions of the monorepo, and `uv sync --locked`): the models (`run_simulation --action factory`), the simulations, with a model which has adopted fits the simulations with the confidence interval and the reports of the adopted fits from their records (`run_report --adopted -o results/fit_reports`, one folder per fit with the HTML report and the methods report `report.pdf`, copied to `fit/<fit id>`).
5. **Release notes** `release-notes/<version>.md`, generated from the changes since the last release.
6. **Docs and README** of the release repository. Both start with the same badges: the DOI, the latest GitHub release, the documentation (the status of the docs workflow, linked to the site), the tests, the Docker image (an image of Docker Hub only), the Python version and the licenses of the code and of the models.
7. **Prune and manifest**: a written file which git ignores is refused, since the commit would lack it (the release repository ignores `results/`, `build/`, `dist/` and `site/` only at its top, a model package may track such folders); every file which the steps did not write is removed (besides `.git`, ignored files and the release notes of earlier versions), then `release.json` records the version, the commit of the monorepo, the versions of the core libraries and a hash of every generated file.
8. `--no-publish` stops here.
9. **Publish**: the lock must hold packages of PyPI only; then the build of the Docker image (for `linux/amd64`, before anything is published, so an image which does not build stops the release before the tag, which cannot be moved), the commit, the tag (the plain version), GitHub Pages (set up to deploy from the docs workflow, before the push which starts it), the push, the push of the Docker image (`:<version>` and `:latest`) and the GitHub release with the release notes and the COMBINE archive. Zenodo archives the release, the docs workflow deploys the site.

## Review, then publish

A run without `--no-publish` which finds the reviewed build of the same release (same version, same commit of the monorepo, the same pkpdsim locked by the monorepo, every file as it was built) commits it as it is, without running the simulations again. A file changed during the review makes the run build again, so the published release is always the output of the command.

A build which failed or was interrupted (a failing simulation, a stopped command) is built again by the next run; nothing has to be cleaned up by hand. The build records the files it changed before it stopped, so a file added by hand afterwards is still refused.

A run which failed while publishing is finished by running the same command again: every publishing step checks whether it is done already (the commit, the tag, the push, the image, the GitHub release) and is skipped then. If the last commit of the release repository is the release of the version (the tag is missing or at that commit), the build is skipped and the publishing steps resume, even if the monorepo moved on since.

The results of the build (`results/`), the environment (`.venv/`) and the settings of an editor (`.idea/`) are ignored by git and stay. A file added to the release repository by hand is refused; the release repository holds only what the command generates.
