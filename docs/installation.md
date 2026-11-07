# Installation

`pkpdsim` requires python >= 3.14 and is available from [pypi](https://pypi.org/project/pkpdsim/). It is tested on Linux, macOS and Windows with python 3.14 and on Linux with python 3.15. It simulates with sbmlsim 0.18 (`sbmlsim>=0.18.0,<0.19`, installed with it); pkpdsim 0.4 is the version for models on sbmlsim 0.8.

## With uv

[uv](https://docs.astral.sh/uv/) is the recommended way to install the package. In a project it is added as a dependency:

```bash
uv add pkpdsim
```

Into an existing virtual environment it is installed through the pip interface of uv:

```bash
uv venv --python 3.14
uv pip install pkpdsim
```

## With pip

```bash
pip install pkpdsim
```

## AMICI, the second simulator

The battery of the checks of a model compares its builds with [AMICI](https://github.com/AMICI-dev/AMICI), the third engine of sbmlsim, when it is installed, see [Battery of the builds](credibility.md#battery-of-the-builds). AMICI is the optional extra `amici`, the extra `amici` of sbmlsim; it is built from source and compiles every model, so it needs a C++ compiler, a BLAS library and the headers of HDF5 (see the installation of sbmlsim):

```bash
uv pip install "pkpdsim[amici]"
```

## Development version

The current state of the `develop` branch is installed directly from GitHub:

```bash
uv add "pkpdsim @ git+https://github.com/matthiaskoenig/pkpdsim.git@develop"
```

or, with pip,

```bash
pip install git+https://github.com/matthiaskoenig/pkpdsim.git@develop
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

An interpreter of the system or of conda is not affected, nor are macOS and Windows. On python 3.15 the wheel of libroadrunner 2.11.0 takes the symbol `__jit_debug_register_code` from `libpython`, which only an interpreter built with the experimental JIT defines: the interpreters of uv are, others may not be (deadsnakes on Ubuntu is not), and `import roadrunner` fails with `undefined symbol: __jit_debug_register_code` ([roadrunner#1406](https://github.com/sys-bio/roadrunner/issues/1406)); use an interpreter of uv there.

## Logging

`pkpdsim` does not configure logging; it logs to the loggers below the `pkpdsim` logger and leaves handlers, levels and formatting to the application.
