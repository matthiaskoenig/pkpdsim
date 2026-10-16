# pkpdsim: building blocks of PBPK/PD models
[![GitHub Actions CI/CD Status](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml)
[![Documentation](https://img.shields.io/badge/docs-pkpdsim-3f51b5.svg)](https://matthiaskoenig.github.io/pkpdsim)
[![Version](https://img.shields.io/pypi/v/pkpdsim.svg)](https://pypi.org/project/pkpdsim/)
[![Python Versions](https://img.shields.io/pypi/pyversions/pkpdsim.svg)](https://pypi.org/project/pkpdsim/)
[![MIT License](https://img.shields.io/pypi/l/pkpdsim.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23297325.svg)](https://doi.org/10.5281/zenodo.23297325)

`pkpdsim` is a python package with the building blocks shared by physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in [SBML](https://sbml.org). The models of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models) are written in python with the factory of [sbmlutils](https://github.com/matthiaskoenig/sbmlutils) and simulated with [sbmlsim](https://github.com/matthiaskoenig/sbmlsim); what they share lives in `pkpdsim`, so that the folder of a model imports nothing but itself and released packages.

Features include

- **templates and annotations** - the creators and the terms of use of a model, the annotations of the compartments which recur in the models (liver, kidneys, gut, blood, urine, ...)
- **study data** - the study folders of [PK-DB](https://pk-db.com) in a checkout of [pkdb_data](https://github.com/matthiaskoenig/pkdb_data), found through the environment variable `PKDB_DATA`, and the data of a model collected from them into its own folder
- **initial values** - the same initial concentration or amount of a substance in every tissue of a model, which removes the distribution kinetics
- **licenses** - the texts of CC BY-SA 4.0 and MIT for the COMBINE archive of a model, written byte for byte the same on every platform
- **models and their commands** - the model object and the commands of a model (`run`, `fit`, `report`, `identifiability`), the tests of a model as functions and the verification of a change, also across versions of sbmlsim (`pkpdsim-verify`)
- **experiments and their metadata** - the base class of the simulation experiments of a model (the SIds of its keys, the units of its data), the presets of the impaired organ function, the shared vocabulary of the metadata of a fit mapping and the filters which select the mappings of a fit by it
- **scans and sensitivity analyses** - the timecourse of a point of a scan as sbmlsim 0.8 gave it and the table of the parameters, the distributions and the runner of the sensitivity analyses of a script

The documentation is available at [https://matthiaskoenig.github.io/pkpdsim](https://matthiaskoenig.github.io/pkpdsim).

If you have any questions or issues please [open an issue](https://github.com/matthiaskoenig/pkpdsim/issues).

## Quickstart

```python
from pathlib import Path

from pkpdsim.annotations import compartments
from pkpdsim.licenses import License, write_license
from pkpdsim.templates import creators

print(creators[0].givenName, creators[0].familyName)
print([resource for _, resource in compartments["li"]])

results = Path("results")
results.mkdir()
print(write_license(results, License.MIT, copyright="2026 Matthias König").name)
```

```text
Matthias König
['fma/FMA:7197', 'bto/BTO:0000759', 'ncit/C12392']
mit.txt
```

The [user guide](https://matthiaskoenig.github.io/pkpdsim/templates/) of the documentation describes every module with a snippet which runs.

## How to cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23297325.svg)](https://doi.org/10.5281/zenodo.23297325)

If you use `pkpdsim` please cite the archived software on [Zenodo](https://doi.org/10.5281/zenodo.23297325):

> König, M. (2026). *pkpdsim: building blocks of physiologically based pharmacokinetic and pharmacodynamic models in SBML* (Version 0.5.4) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23357272

```bibtex
@software{konig_pkpdsim,
  author    = {König, Matthias},
  title     = {pkpdsim: building blocks of physiologically based pharmacokinetic and pharmacodynamic models in SBML},
  year      = {2026},
  month     = oct,
  version   = {0.5.4},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23357272},
  url       = {https://doi.org/10.5281/zenodo.23357272},
}
```

## Installation

`pkpdsim` requires python >= 3.13, simulates with sbmlsim 0.10 (`sbmlsim>=0.10.0,<0.11`, installed with it) and is available from [pypi](https://pypi.org/project/pkpdsim/):

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
