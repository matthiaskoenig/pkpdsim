# pkpdsim: building blocks of PBPK/PD models
[![GitHub Actions CI/CD Status](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml)
[![Documentation](https://img.shields.io/badge/docs-pkpdsim-3f51b5.svg)](https://matthiaskoenig.github.io/pkpdsim)
[![Version](https://img.shields.io/pypi/v/pkpdsim.svg)](https://pypi.org/project/pkpdsim/)
[![Python Versions](https://img.shields.io/pypi/pyversions/pkpdsim.svg)](https://pypi.org/project/pkpdsim/)
[![MIT License](https://img.shields.io/pypi/l/pkpdsim.svg)](https://opensource.org/licenses/MIT)

`pkpdsim` is a python package with the building blocks shared by physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in [SBML](https://sbml.org). The models of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models) are written in python with the factory of [sbmlutils](https://github.com/matthiaskoenig/sbmlutils) and simulated with [sbmlsim](https://github.com/matthiaskoenig/sbmlsim); what they share lives in `pkpdsim`, so that the folder of a model imports nothing but itself and released packages.

The documentation is available at [https://matthiaskoenig.github.io/pkpdsim](https://matthiaskoenig.github.io/pkpdsim).

If you have any questions or issues please [open an issue](https://github.com/matthiaskoenig/pkpdsim/issues).

## How to cite

If you use `pkpdsim` please cite the archived software on [Zenodo](https://zenodo.org), whose DOI is added here with the first release.

## Installation

`pkpdsim` requires python >= 3.13 and is available from [pypi](https://pypi.org/project/pkpdsim/):

```bash
uv add pkpdsim
```

or with pip

```bash
pip install pkpdsim
```

See [Installation](https://matthiaskoenig.github.io/pkpdsim/installation/) for details and [Development](https://matthiaskoenig.github.io/pkpdsim/development/) for working on the repository.

## License

- Source Code: [MIT](https://opensource.org/license/MIT)
- Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

## Funding

Matthias König is supported by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151 "QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach)" by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA).

Matthias König was supported by the Federal Ministry of Education and Research (BMBF, Germany) within the research network Systems Medicine of the Liver (**LiSyM**, grant number 031L0054).

© 2026 Matthias König.
