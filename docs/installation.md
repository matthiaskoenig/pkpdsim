# Installation

`pkpdsim` requires python >= 3.14 and is available from [pypi](https://pypi.org/project/pkpdsim/). It is tested on Linux, macOS and Windows with python 3.14 and on Linux with python 3.15. It simulates with sbmlsim 0.22 (`sbmlsim>=0.22.0`, installed with it), whose simulator [libroadrunner](https://libroadrunner.org) is not on PyPI: sbmlsim relies on fixes of the events of libroadrunner which are not released yet, so it needs a build of the `develop` branch, `libroadrunner>=2.11.1.dev655`, which the package index [roadrunner-wheels](https://github.com/matthiaskoenig/roadrunner-wheels) serves. Install the build from that index first and then `pkpdsim` from PyPI, whose requirement the installed build satisfies:

```bash
pip install --no-deps --index-url https://matthiaskoenig.github.io/roadrunner-wheels/simple/ "libroadrunner>=2.11.1.dev655"
pip install pkpdsim
```

The index serves the wheels of libroadrunner only, so the first command takes nothing else from it (`--no-deps`, the numpy of libroadrunner comes from PyPI with `pkpdsim`). `pip install pkpdsim` alone stops with `No matching distribution found for libroadrunner>=2.11.1.dev655`. Once libroadrunner is released with these fixes, the release on PyPI satisfies the requirement as well and the index is no longer needed.

## With uv

[uv](https://docs.astral.sh/uv/) is the recommended way to install the package. In a project the index of the builds of libroadrunner is declared in its `pyproject.toml`, which takes libroadrunner from it:

```toml
[tool.uv.sources]
libroadrunner = { index = "roadrunner-dev" }

[[tool.uv.index]]
name = "roadrunner-dev"
url = "https://matthiaskoenig.github.io/roadrunner-wheels/simple/"
explicit = true
```

and the package is added together with libroadrunner, which locks both with the rest of the environment:

```bash
uv add pkpdsim libroadrunner
```

uv applies a source to the direct dependencies of a project only, so `libroadrunner` is one of them; the index is `explicit`, so uv asks it for libroadrunner and takes every other package from PyPI. A project which already depends on `pkpdsim` adds the lines and runs `uv add libroadrunner`; `uv lock --upgrade-package libroadrunner` moves the lock to the latest build. The sources are not part of a wheel built from the project, so a package which depends on `pkpdsim` is installed with the index as well.

Into an existing virtual environment libroadrunner is installed from the index first and then the package, both through the pip interface of uv:

```bash
uv venv --python 3.14
uv pip install --no-deps --index-url https://matthiaskoenig.github.io/roadrunner-wheels/simple/ "libroadrunner>=2.11.1.dev655"
uv pip install pkpdsim
```

`--extra-index-url` in place of the two commands asks the index for every package of the install, and GitHub Pages, which serves it, answers a request for a package it does not have with a `503` now and then, which fails the install.

## With pip

The two commands at the top install the build of libroadrunner from the index and then the package from PyPI. A single

```bash
pip install --extra-index-url https://matthiaskoenig.github.io/roadrunner-wheels/simple/ pkpdsim
```

installs the build as well; it asks the index for every package of the install, which GitHub Pages, serving the index, answers with a `503` now and then for a package it does not have. pip skips such an answer, `uv pip install` fails on it.

## AMICI, the second simulator

The battery of the checks of a model compares its builds with [AMICI](https://github.com/AMICI-dev/AMICI), the third engine of sbmlsim, when it is installed, see [Battery of the builds](credibility.md#battery-of-the-builds). AMICI is the optional extra `amici`, the extra `amici` of sbmlsim; it is built from source and compiles every model, so it needs a C++ compiler, a BLAS library and the headers of HDF5 (see the installation of sbmlsim). Install it into the environment of `pkpdsim`, where the installed build of libroadrunner satisfies the requirement, so the index is needed only when the install upgrades libroadrunner:

```bash
uv pip install "pkpdsim[amici]"
```

## Development version

The current state of the `develop` branch is installed directly from GitHub:

```bash
uv add "pkpdsim @ git+https://github.com/matthiaskoenig/pkpdsim.git@develop" libroadrunner
```

or, with pip,

```bash
pip install --extra-index-url https://matthiaskoenig.github.io/roadrunner-wheels/simple/ git+https://github.com/matthiaskoenig/pkpdsim.git@develop
```

To work on the repository itself, with the test and documentation tooling, see [Development](development.md).

## Dependencies

| package | used for |
| --- | --- |
| [sbmlutils](https://github.com/matthiaskoenig/sbmlutils) | the factory of the models (creators, units, compartments), the qualifiers of the annotations and the console |
| [libroadrunner](https://libroadrunner.org) | the simulator, the initial values are selections of a roadrunner model |

libroadrunner ships binary wheels for Linux (x86_64), macOS and Windows. Its extension links against `libpython`: an interpreter which uv installs keeps `libpython` in its `lib` directory, where the loader of Linux does not look, and `import roadrunner` fails with `libpython3.14.so.1.0: cannot open shared object file`. Put the directory on the library search path, e.g. in a project managed by uv:

```bash
export LD_LIBRARY_PATH="$(uv run python -c 'import sys; print(sys.base_prefix)')/lib"
```

An interpreter of the system or of conda is not affected, nor are macOS and Windows. On python 3.15 the wheel of libroadrunner 2.11 takes the symbol `__jit_debug_register_code` from `libpython`, which only an interpreter built with the experimental JIT defines: the interpreters of uv are, others may not be (deadsnakes on Ubuntu is not), and `import roadrunner` fails with `undefined symbol: __jit_debug_register_code` ([roadrunner#1406](https://github.com/sys-bio/roadrunner/issues/1406)); use an interpreter of uv there.

## Logging

`pkpdsim` does not configure logging; it logs to the loggers below the `pkpdsim` logger and leaves handlers, levels and formatting to the application.
