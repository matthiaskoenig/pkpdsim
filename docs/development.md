# Development

Contributions are welcome. The repository is [matthiaskoenig/pkpdsim](https://github.com/matthiaskoenig/pkpdsim); development happens against the `develop` branch via pull requests.

## Branch model

Two branches are permanent:

- **`develop`** is the default branch and the branch everything is integrated into. The documentation on [matthiaskoenig.github.io/pkpdsim](https://matthiaskoenig.github.io/pkpdsim) is published from it.
- **`main`** tracks the latest published release. It is fast-forwarded to the released commit by the `sync-main` job of the `CI-CD` workflow after the package went to pypi, so `main` and the newest version on pypi always agree. Nothing is developed on `main` and nothing is merged into it by hand.

Work happens on short lived branches off `develop`, which GitHub deletes after the merge. Releases are tagged on `develop`, see [Release](#release).

## Pull requests

Neither branch accepts a direct push, every change goes through a pull request against `develop`. This includes the maintainer, there is no bypass.

A pull request can only be merged once the four required checks are green:

| check   | workflow      | content                                                              |
| ------- | ------------- | -------------------------------------------------------------------- |
| `tests` | `ci-cd.yml`   | the test matrix, python 3.14 on linux, macos and windows              |
| `ruff`  | `ruff.yml`    | `ruff check` and `ruff format --check`                                |
| `ty`    | `ty.yml`      | `tox r -e ty`                                                         |
| `docs`  | `docs.yml`    | the zensical build including the api reference and the agent files    |

`tests` aggregates the test matrix into a single job, so the name of the required check stays the same when the matrix changes.

Continuous integration is kept small: every workflow cancels its running build when a newer commit of the same branch or pull request arrives, uv caches the packages and the interpreters between runs, dependabot proposes its updates once a month, and the matrix tests only python 3.14. Python 3.13 is tested locally, see [Testing](#testing), before a pull request is opened.

Further rules of a pull request:

- conversations have to be resolved before the merge
- an approval is dismissed when new commits are pushed
- the history stays linear, i.e., a pull request is merged with squash or rebase; merge commits are disabled
- the maintainer is the code owner of the repository (`.github/CODEOWNERS`) and is requested for review on every pull request. A pull request of a contributor is therefore reviewed and merged by the maintainer, who has the only write access. The rulesets themselves do not require an approval: on a personal repository a ruleset cannot ask for an approval only from somebody else, and requiring one would block the pull requests of the maintainer, who cannot approve their own. Once a second person has write access, a ruleset requiring an approving review of a code owner can be added

[Auto-merge](https://docs.github.com/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/automatically-merging-a-pull-request) is enabled for the repository, so a pull request can be queued and is merged as soon as the checks pass and the required approval is there.

### Repository policies { #repository-policies }

The protection is implemented with [repository rulesets](https://docs.github.com/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/about-rulesets). They are part of the repository in `.github/rulesets/` instead of only living in the web interface, so a change to a policy is reviewed like any other change:

| ruleset                 | applies to | rules                                                                                                                                       |
| ----------------------- | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `develop.json`          | `develop`  | pull request required, the four checks above, resolved conversations, linear history, no force push, no deletion. **No bypass, for anybody.** |
| `main.json`             | `main`     | no force push, no deletion, no bypass. The fast-forward of the release workflow needs none, only a force push would be rejected |
| `tags.json`             | all tags   | a tag cannot be deleted or moved, so a release tag keeps pointing at what was released                                                       |

Changing a policy means changing the json and applying it:

```bash
.github/rulesets/apply.sh
```

The script is idempotent: it updates the rulesets which exist and creates the missing ones. It also sets the merge settings of the repository, i.e., auto-merge, delete branch on merge, and squash and rebase as the only merge methods. It needs the [github cli](https://cli.github.com) authenticated as a user with admin permission on the repository.

## Setup development environment

Development needs [uv](https://docs.astral.sh/uv/) and a checkout of the repository:

```bash
git clone https://github.com/matthiaskoenig/pkpdsim.git
cd pkpdsim
```

A single sync creates the virtual environment in `.venv`, installs `pkpdsim` into it in editable mode and adds the complete tooling:

```bash
uv sync --extra dev
```

The environment is resolved from `uv.lock`, which is committed, so local development, including `uv run ty check`, and the documentation build use the same versions. The tox environments, i.e. the test matrix and the `ty` check of continuous integration, resolve from `pyproject.toml` and do not use the lock; the lower bounds in `pyproject.toml` are what a user of the library installs against. After changing a dependency in `pyproject.toml` run `uv lock`: the `documentation` workflow syncs with `uv sync --locked`, which fails when `pyproject.toml` and the lock disagree. `uv lock --upgrade` moves the lock to the newest releases. Dependabot bumps the locked dependencies as one grouped monthly pull request, with `versioning-strategy: lockfile-only`, so the lower bounds in `pyproject.toml` are never raised by it and only ever changed by hand.

The `dev` extra contains everything used below, i.e., pytest, ruff, ty, tox, pre-commit, zensical and bump-my-version, so nothing has to be installed separately. The python version is taken from `.python-version` (currently 3.14); to work against the oldest supported version instead use `uv sync --extra dev --python 3.13`, which replaces the environment.

The tools are then run either with `uv run <command>`, which uses the environment without activating it, or from the activated environment:

```bash
source .venv/bin/activate        # Linux and macOS
.venv\Scripts\activate           # Windows
```

The commands in this document are written without the `uv run` prefix; prepend it if the environment is not activated.

The last step installs the git hook:

```bash
uv run pre-commit install          # install the hook, once per checkout
uv run pre-commit run --all-files  # check the current state of the repository
```

From now on every commit is checked with ruff (lint and format) and ty, i.e., the same checks that run in continuous integration. On a commit only the changed files are looked at, `--all-files` checks the whole repository and is what a newly added hook should be tried with.

### libroadrunner on Linux

The interpreters of uv keep `libpython`, which the extension of libroadrunner links against, outside the search path of the loader on Linux, see [Installation](installation.md#dependencies). Before running the tests:

```bash
export LD_LIBRARY_PATH="$(uv run python -c 'import sys; print(sys.base_prefix)')/lib"
```

The workflows do the same in `.github/actions/setup`.

## Testing

The tests are written with pytest, tox runs them against every supported python version.

The tox environments are named after the interpreter (`py3.13` and `py3.14`, see `envlist` in `tox.ini`), a single one is run with

```bash
tox r -e py3.14
```

and the complete matrix, including the `ty` environment, in parallel with

```bash
tox run-parallel
```

This needs the interpreters to be available, which uv installs with `uv python install 3.13 3.14`. The tox environments build the wheel and install it, so the tests run against the package as it is published. It is the complete test and is run before a pull request is opened: continuous integration runs only `py3.14`, on linux, macos and windows, as `uvx --with tox-uv tox -e py3.14 -- --durations=15`, which lists the slowest tests in the log.

To run the tests directly against the development environment use

```bash
pytest                                       # the full suite
pytest -n 0                                  # in one process, e.g. for --pdb
pytest tests/test_package.py                 # a single module
pytest tests/test_package.py::test_version   # a single test
```

The tests run in parallel, `addopts = "-n auto"` in `pyproject.toml` gives pytest-xdist one worker per core; `-n 0` on the command line runs everything in one process, which the debugger needs.

## Linting and formatting

Linting and formatting use [ruff](https://docs.astral.sh/ruff/):

```bash
ruff check     # lint
ruff format    # format
```

The `ruff` workflow installs the ruff version locked in `uv.lock`, so continuous integration and `uv run ruff` agree.

The docstring rules (`D`) are enforced for the package, not for `tests/` and `scripts/`, see `.ruff.toml`; the annotation rules (`ANN`) hold everywhere.

## Type checking

Type checking is performed with [ty](https://docs.astral.sh/ty/):

```bash
tox r -e ty
```

Or directly in the working tree:

```bash
uv run ty check
```

The configuration lives in `[tool.ty]` in `pyproject.toml`. Warnings are treated as errors, so the codebase is kept free of diagnostics. Suppress an unavoidable diagnostic with a rule specific `# ty: ignore[rule-name]` rather than a blanket comment.

## Documentation

The documentation is built with [Zensical](https://zensical.org/), the static site generator of the Material for MkDocs authors. The sources are markdown files in `docs/`, the site is configured in `zensical.toml` in the repository root. Nothing rendered is committed: the site is built by the `documentation` workflow on every push and published to [matthiaskoenig.github.io/pkpdsim](https://matthiaskoenig.github.io/pkpdsim) from the `develop` branch.

Build the site into `site/`:

```bash
uv run zensical build --clean --strict
```

The build is strict: a broken link or a missing page fails it, in continuous integration as well.

For writing, the preview rebuilds on save:

```bash
uv run zensical serve
```

The API reference is rendered from the docstrings by [mkdocstrings](https://mkdocstrings.github.io/); a page in `docs/api/` only contains the module directive:

```markdown
# pkpdsim

::: pkpdsim
```

Docstrings are therefore the place to document functions and classes, the markdown files provide the narrative around them. Adding a module to the reference means adding such a page, an entry to `nav` in `zensical.toml` and a row to `docs/api/index.md`.

### Snippets of the documentation

Every ` ```python ` block of the pages follows two rules:

- **It runs.** The first block of a page is self-contained (its imports, its data, the call and the output it prints) and runs from the root of the repository with warnings as errors. A later block of the same page may be a fragment, but then it names in a comment or in the sentence before it where every object it uses comes from. The output a snippet prints is shown below it as a `text` block and is pasted from a run, never written by hand.

    `tests/docs/test_snippets.py` keeps this honest: it runs the blocks of every page in the order they appear and in one namespace per page, in a temporary working directory, in a subprocess with `-W error`. A fragment which needs what a page cannot provide (a checkout of pkdb_data, a model of a reader) carries the comment `# not executed` as its first line and is skipped; every other block has to run.

- **It is formatted.** `ruff format` formats the code blocks of the markdown files as well, so `ruff format --check` covers the documentation and a snippet is written the way ruff would write it.

### Files for agents { #files-for-agents }

Agents and language models read markdown, not rendered html. `scripts/llms_txt.py` writes the files of the [llms.txt convention](https://llmstxt.org/) into the built site, i.e., [llms.txt](https://matthiaskoenig.github.io/pkpdsim/llms.txt) as an annotated index of all pages, [llms-full.txt](https://matthiaskoenig.github.io/pkpdsim/llms-full.txt) with the complete documentation in a single file, and the markdown of every page next to its html (`/installation.md` for `/installation/`). The markdown of the API reference is generated from the docstrings with `inspect`, since the pages themselves only contain the mkdocstrings directive.

```bash
uv run zensical build --clean --strict
uv run python scripts/llms_txt.py
```

The `documentation` workflow runs both steps, so the files are regenerated with every push. `docs/robots.txt` points crawlers at the sitemap and at these files. Zensical will provide agent context files itself at some point, then this script can go.

## Release

A release is made from `develop`. Since `develop` only accepts pull requests, the release is prepared on a branch and tagged once that pull request is merged:

1. branch off `develop`: `git switch -c release/x.y.z develop`
2. write the release notes for the version in `release-notes/x.y.z.md`
3. make sure everything passes: `tox run-parallel`, `ruff check`, `tox r -e ty`
4. check the version bump: `uvx bump-my-version bump [dev|major|minor|patch] --dry-run -vv`. A development version (`x.y.z.devN`) is finalized with `uvx bump-my-version bump dev`, which drops the `.devN` suffix (`0.2.0.dev0` becomes `0.2.0`); `major`, `minor` and `patch` start the next development cycle instead (`0.2.0.dev0` becomes `1.0.0.dev0`, `0.3.0.dev0`, `0.2.1.dev0`)
5. bump the version: `uvx bump-my-version bump [dev|major|minor|patch]`, which updates `src/pkpdsim/__init__.py`, `CITATION.cff`, `tests/test_package.py` and `.bumpversion.toml` (its own `current_version`) and commits. Use `dev` to release the current development version and `major|minor|patch` to open the next one. It does not create the tag; a squash or rebase merge would rewrite the commit and leave the tag behind on a commit which is not part of `develop`
6. push the branch, open the pull request against `develop` and merge it once the checks are green
7. tag the merged commit on `develop` and push the tag:

    ```bash
    git switch develop
    git pull
    git tag x.y.z
    git push origin x.y.z
    ```

    This starts the `CI-CD` workflow, which runs the test matrix, publishes to [pypi](https://pypi.org/project/pkpdsim/) through the trusted publisher of the GitHub environment `pypi`, creates the GitHub release from `release-notes/x.y.z.md` and fast-forwards `main` to the tagged commit. Zenodo archives the GitHub release. Check the version before pushing, a tag cannot be moved or deleted afterwards.

8. test the installation from pypi in a fresh environment:

    ```bash
    uv venv --python 3.14
    uv pip install pkpdsim
    ```

9. once Zenodo has archived the release, update the citation information, i.e., `date-released` and `doi` in `CITATION.cff` and the version, date and DOI of the release in the citation of `README.md` and `docs/index.md`. `bump-my-version` only updates the version, not the date and the DOI, which are only known after the release. These changes go in through a pull request like everything else
