# Installation

`pkpdsim` requires python >= 3.13 and is available from [pypi](https://pypi.org/project/pkpdsim/). It is tested on Linux, macOS and Windows with python 3.13 and 3.14.

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
