# pkpdsim: building blocks of PBPK/PD models
[![GitHub Actions CI/CD Status](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml) [![Documentation](https://img.shields.io/badge/docs-pkpdsim-3f51b5.svg)](https://matthiaskoenig.github.io/pkpdsim) [![Version](https://img.shields.io/pypi/v/pkpdsim.svg)](https://pypi.org/project/pkpdsim/) [![Python Versions](https://img.shields.io/pypi/pyversions/pkpdsim.svg)](https://pypi.org/project/pkpdsim/) [![MIT License](https://img.shields.io/pypi/l/pkpdsim.svg)](https://opensource.org/licenses/MIT)

`pkpdsim` is a python package with the building blocks shared by physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in [SBML](https://sbml.org). The source code is available from [https://github.com/matthiaskoenig/pkpdsim](https://github.com/matthiaskoenig/pkpdsim).

## Background

The models of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models) describe the pharmacokinetics and pharmacodynamics of a substance in the whole body: one SBML model per tissue, composed into a whole-body model, written in python with the factory of [sbmlutils](https://github.com/matthiaskoenig/sbmlutils), simulated with [sbmlsim](https://github.com/matthiaskoenig/sbmlsim) against the studies curated in [PK-DB](https://pk-db.com). `pkpdsim` holds what these models share, so that the folder of a model imports nothing but itself and released packages and can live in a repository of its own.

## Quickstart

```python
import pkpdsim

print(pkpdsim.__version__)
```

```text
0.1.0
```

Continue with [Installation](installation.md).

## How to cite

If you use `pkpdsim` please cite the archived software on [Zenodo](https://zenodo.org), whose DOI is added here with the first release.

## License

- Source Code: [MIT](https://opensource.org/license/MIT)
- Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

## Funding

Matthias König is supported by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151 "QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach)" by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA).

Matthias König was supported by the Federal Ministry of Education and Research (BMBF, Germany) within the research network Systems Medicine of the Liver (**LiSyM**, grant number 031L0054).
