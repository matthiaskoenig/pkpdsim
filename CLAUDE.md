# CLAUDE.md

This file provides guidance when working with code in this repository.

## Project

`pkpdsim` is a python package with the building blocks shared by physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in SBML, i.e., the models of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models): what these models share moves here, so that the folder of a model imports nothing but itself and released packages ([pkdb_models#95](https://github.com/matthiaskoenig/pkdb_models/issues/95)). The five model commands take the dotted name of a model package first.

## Commands

```bash
# environment (uv based)
uv sync --extra dev
uv run pre-commit install
uv lock                                  # after every change of a dependency in pyproject.toml

# tests
tox r -e py3.14                          # single tox env (py3.14, py3.15 available)
tox run-parallel                         # full matrix + ty, from uv.lock (runner uv-venv-lock-runner of tox-uv, extra dev), run before opening a pull request

# lint / format / types
tox -e ty                       # ty type check (config in [tool.ty] in pyproject.toml)
uv run ty check                 # same check, straight from the working tree

# docs
uv run zensical build --clean --strict
uv run python scripts/llms_txt.py
uv run zensical serve
```

`develop` is the default branch and takes every change through a pull request; direct pushes are rejected by the rulesets in `.github/rulesets/`. Python 3.15 runs only locally. `main` only tracks the latest release and is fast-forwarded by the `sync-main` job of the release workflow, never by hand.

Release steps are in `docs/development.md`: the release is prepared on a branch and the tag (a plain version, `0.1.0`) is created on `develop` only after the pull request was merged, which triggers the release job of `ci-cd.yml`.

Documentation is [Zensical](https://zensical.org/) (`zensical.toml`, built into the gitignored `site/`); `scripts/llms_txt.py` runs after the build. A page in `docs/api/` is just `::: pkpdsim.<module>`. The pages of the PBPK library (`docs/pbpk/`) have generated tables, tested equations and inline SVG diagrams, see `pbpk/` below.

## Architecture

The module by module notes live next to this file; read the one that covers what you are touching.

- [.claude/architecture/core.md](.claude/architecture/core.md) - the shared modules: `__init__.py` and its lazy imports, `templates.py`, `annotations.py`, `data.py`, `initial_values.py`, `timecourse.py` and the licenses of a model.
- [.claude/architecture/models.md](.claude/architecture/models.md) - the model and run layer: `model.py`, the run of the experiments, `cli.py`, `fit.py`, `metadata.py`, `experiment.py`, `pk.py`, `sensitivity.py`, the uncertainty bands, the reports of the adopted fits, `credibility/`, `database/`, `impairment.py`, `verify.py` and `release/`.
- [.claude/architecture/pbpk-definition.md](.claude/architecture/pbpk-definition.md) - the PBPK definition: `pbpk/`, substances and proteins, chemistry, variants and routes, PD components and observers, slow pools, the generated steady state, the archive and the figures.
- [.claude/architecture/pbpk-physiology.md](.claude/architecture/pbpk-physiology.md) - physiology and populations: age and sex, species, circadian rhythms, pregnancy, adipose tissue, the build, the parts, the experiments of a library model, the kidney validity, the equilibration and the tested docs.
- [.claude/architecture/perfusion.md](.claude/architecture/perfusion.md) - the machine perfusion library `perfusion/`: the definition, the id scheme, the cellular core, the circuit, the builder, the observers, the figures, the archive and the experiments with medium changes.
- [.claude/architecture/tests.md](.claude/architecture/tests.md) - what the tests share: the flat model cache, the toy package and the `O0` compilation of the PBPK tests.

## Conventions

- Type checking is done with [ty](https://docs.astral.sh/ty/). The optional extra `amici` (the extra `amici` of sbmlsim, built from source, not in the dev environment) is imported by sbmlsim only, pkpdsim never imports it. `[tool.ty.terminal] error-on-warning = true` means warnings fail the check, so the tree must stay at zero diagnostics. Suppress a diagnostic with a rule-specific `# ty: ignore[rule-name]`, never a blanket `# type: ignore`.
- Modules get their logger with `logging.getLogger(__name__)` and the package never configures logging; the console output of `collect_tsv_files` is the exception, it reports the copies of a script, and its findings for the curation are `CollectionWarning`s (`warnings.warn`).
- A dependency is only added for an import of the package; the floors in `pyproject.toml` are the releases the package was tested with. The exception is a floor of a security fix of a transitive dependency: `certifi>=2026.10.19` (GHSA-6q48-6c7w-9gw2, the CA bundle of requests, which sbmlutils and pymetadata use).
- `tests/conftest.py` loads the registry of identifiers.org, which sbmlutils resolves annotations with, into the cache of pymetadata once in the controller of the session; without it the parallel workers download it at the same time and the replace of the cache file fails on Windows (matthiaskoenig/pymetadata#102).
- `tests/docs/test_snippets.py` executes every ` ```python ` fence of every page of `docs/` in the order it appears, one namespace per page, in a subprocess with warnings as errors and a temporary working directory. A fence whose first line is `# not executed` is a fragment and is skipped. An output a page pastes has to come from a real run of that snippet.
- `uv.lock` is committed and pins the local development environment (`uv sync`, `uv run`, so also `uv run ty check`), the `documentation` workflow, which runs `uv sync --locked`, and the ruff version of the `ruff` workflow; the tox environments, i.e. the test matrix and the `ty` check of continuous integration, too (the runner `uv-venv-lock-runner` of tox-uv, `uv sync --locked` with the extra `dev`), so a source of `[tool.uv.sources]` holds everywhere. After changing a dependency in `pyproject.toml` run `uv lock`, otherwise the `documentation` workflow fails; `uv lock --upgrade` moves the lock to the newest releases. Dependabot bumps the lock monthly with `versioning-strategy: lockfile-only`, so the lower bounds in `pyproject.toml`, which users install against, are only changed by hand. The ruff and ty floors of the `dev` extra and the `rev` of their pre-commit hooks are the locked versions; update them together.
- libroadrunner links against `libpython`, which the interpreters of uv ship in their `lib` directory, outside the search path of the loader on Linux: export `LD_LIBRARY_PATH="$(uv run python -c 'import sys; print(sys.base_prefix)')/lib"` before running the tests (a system interpreter with its own `libpython` does not need it); the workflows do it in `.github/actions/setup`. tox takes the interpreters of uv only (`uv_python_preference = only-managed`; the libroadrunner 2.11.0 wheel for 3.15 needs `__jit_debug_register_code` of a `libpython` built with the experimental JIT, which deadsnakes lacks, sys-bio/roadrunner#1406) and its inline plugin `toxfile.py` puts the `lib` of the interpreter of each environment on `LD_LIBRARY_PATH`, so `tox` needs no export (without it roadrunner loaded the system's `libpython` of the same version, a second runtime in the process).
- Markdown carries no hard line wraps: a paragraph, a list item or a table row is a single line. Code fences, headings and the rows of badges keep their line structure.
- Release notes go in `release-notes/` as part of a release commit.
