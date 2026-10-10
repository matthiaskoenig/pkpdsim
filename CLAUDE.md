# CLAUDE.md

This file provides guidance when working with code in this repository.

## Project

`pkpdsim` is a python package with the building blocks shared by physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in SBML, i.e., the models of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models): what these models share moves here, so that the folder of a model imports nothing but itself and released packages ([pkdb_models#95](https://github.com/matthiaskoenig/pkdb_models/issues/95)). Pure library, no CLI entry points. Requires python >= 3.13 (tested on 3.13 and 3.14), packaged with hatchling (version is read from `src/pkpdsim/__init__.py`).

## Commands

```bash
# environment (uv based)
uv sync --extra dev
uv run pre-commit install
uv lock                                  # after every change of a dependency in pyproject.toml

# tests
pytest                                   # all tests, in parallel with pytest-xdist
pytest -n 0                              # in one process, e.g. for --pdb
pytest tests/test_package.py             # single file
pytest tests/test_package.py::test_version  # single test
tox r -e py3.14                          # single tox env (py3.13, py3.14 available)
tox run-parallel                         # full matrix + ty

# lint / format / types
ruff check
ruff format
tox -e ty                       # ty type check (config in [tool.ty] in pyproject.toml)
uv run ty check                 # same check, straight from the working tree

# docs
uv run zensical build --clean --strict
uv run python scripts/llms_txt.py
uv run zensical serve
```

`develop` is the default branch and takes every change through a pull request; direct pushes are rejected by the rulesets in `.github/rulesets/` (applied with `.github/rulesets/apply.sh`), which require the `tests`, `ruff`, `ty` and `docs` checks. `main` only tracks the latest release and is fast-forwarded by the `sync-main` job of the release workflow, never by hand.

Release steps are in `docs/development.md`: the release is prepared on a branch, `uvx bump-my-version bump [dev|major|minor|patch]` updates `src/pkpdsim/__init__.py`, `CITATION.cff` and `tests/test_package.py` and commits without tagging, and the tag (a plain version, `0.1.0`) is created on `develop` after the pull request was merged, which triggers the release job of `ci-cd.yml`: PyPI through the trusted publisher of the GitHub environment `pypi`, the GitHub release from `release-notes/<version>.md` (archived by Zenodo) and the fast-forward of `main`.

Documentation is [Zensical](https://zensical.org/): markdown sources in `docs/`, configured in `zensical.toml`, built into the gitignored `site/`. The API reference is rendered from the docstrings by mkdocstrings; a page in `docs/api/` is just `::: pkpdsim.<module>`. `scripts/llms_txt.py` runs after the build and writes `llms.txt`, `llms-full.txt` and the markdown of every page into `site/`.

## Architecture

`src/pkpdsim/__init__.py` holds the version and nothing else, so `import pkpdsim` stays free of the heavy dependencies; every topic is a module of its own which a model imports directly.

## Conventions

- Type checking is done with [ty](https://docs.astral.sh/ty/). `[tool.ty.terminal] error-on-warning = true` means warnings fail the check, so the tree must stay at zero diagnostics. Suppress a diagnostic with a rule-specific `# ty: ignore[rule-name]`, never a blanket `# type: ignore`. ty also runs as a pre-commit hook (`--extra dev`).
- Linting is ruff with the full rule set of the maintainer's packages including the annotations (`ANN`), see `.ruff.toml`; every module, class and function of the package carries full type annotations and a google style docstring (ruff `D`), `tests/` and `scripts/` are exempt from the docstring rules.
- Library code logs, it does not print; modules get their logger with `logging.getLogger(__name__)` and the package never configures logging.
- `tests/docs/test_snippets.py` executes every ` ```python ` fence of every page of `docs/` in the order it appears, one namespace per page, in a subprocess with warnings as errors and a temporary working directory. A fence whose first line is `# not executed` is a fragment and is skipped. An output a page pastes has to come from a real run of that snippet.
- `uv.lock` is committed and pins the local development environment (`uv sync`, `uv run`, so also `uv run ty check`), the `documentation` workflow, which runs `uv sync --locked`, and the ruff version of the `ruff` workflow; the tox environments, i.e. the test matrix and the `ty` check of continuous integration, resolve from `pyproject.toml` and do not use it. After changing a dependency in `pyproject.toml` run `uv lock`, otherwise the `documentation` workflow fails; `uv lock --upgrade` moves the lock to the newest releases. Dependabot bumps the lock weekly with `versioning-strategy: lockfile-only`, so the lower bounds in `pyproject.toml`, which users install against, are only changed by hand. The ruff and ty floors of the `dev` extra and the `rev` of their pre-commit hooks are the locked versions; update them together.
- libroadrunner links against `libpython`, which the interpreters of uv ship in their `lib` directory, outside the search path of the loader on Linux: export `LD_LIBRARY_PATH="$(uv run python -c 'import sys; print(sys.base_prefix)')/lib"` before running the tests; the workflows do it in `.github/actions/setup`.
- Markdown carries no hard line wraps: a paragraph, a list item or a table row is a single line. Code fences, headings and the rows of badges keep their line structure.
- Release notes go in `release-notes/` as part of a release commit.
