# pkpdsim: building blocks of PBPK/PD models
[![GitHub Actions CI/CD Status](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml) [![Documentation](https://img.shields.io/badge/docs-pkpdsim-3f51b5.svg)](https://matthiaskoenig.github.io/pkpdsim) [![Version](https://img.shields.io/pypi/v/pkpdsim.svg)](https://pypi.org/project/pkpdsim/) [![Python Versions](https://img.shields.io/pypi/pyversions/pkpdsim.svg)](https://pypi.org/project/pkpdsim/) [![MIT License](https://img.shields.io/pypi/l/pkpdsim.svg)](https://opensource.org/licenses/MIT)

`pkpdsim` is a python package with the building blocks shared by physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in [SBML](https://sbml.org). The source code is available from [https://github.com/matthiaskoenig/pkpdsim](https://github.com/matthiaskoenig/pkpdsim).

## Background

The models of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models) describe the pharmacokinetics and pharmacodynamics of a substance in the whole body: one SBML model per tissue, composed into a whole-body model, written in python with the factory of [sbmlutils](https://github.com/matthiaskoenig/sbmlutils), simulated with [sbmlsim](https://github.com/matthiaskoenig/sbmlsim) against the studies curated in [PK-DB](https://pk-db.com). `pkpdsim` holds what these models share, so that the folder of a model imports nothing but itself and released packages and can live in a repository of its own.

## Features

- **[Templates and annotations](templates.md)** - the creators and the terms of use of a model, the annotations of the compartments which recur in the models (liver, kidneys, gut, blood, urine, ...) with the terms of FMA, BTO, NCIT, UBERON and GO.
- **[Study data](data.md)** - the study folders of PK-DB in a checkout of pkdb_data, found through the environment variable `PKDB_DATA`, and the data of a model collected from them into its own folder.
- **[Initial values](initial_values.md)** - the same initial concentration or amount of a substance in every tissue of a model, which removes the distribution kinetics.
- **[Licenses](licenses.md)** - the texts of CC BY-SA 4.0 and MIT for the COMBINE archive of a model, written byte for byte the same on every platform.
- **[Models and their commands](models.md)** - the model object and the commands of a model (`run`, `fit`, `report`, `identifiability`), the tests of a model as functions and the verification of a change (`pkpdsim-verify`).
- **[Experiments and their metadata](experiments.md)** - the base class of the simulation experiments of a model, the presets of the impaired organ function, the shared vocabulary of the metadata of a fit mapping and the filters which select the mappings of a fit by it.

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

Continue with [Installation](installation.md) and the [Templates and annotations](templates.md) of a model.

## How to cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23297325.svg)](https://doi.org/10.5281/zenodo.23297325)

If you use `pkpdsim` please cite the archived software on [Zenodo](https://doi.org/10.5281/zenodo.23297325):

> König, M. (2026). *pkpdsim: building blocks of physiologically based pharmacokinetic and pharmacodynamic models in SBML* (Version 0.4.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23298294

## License

- Source Code: [MIT](https://opensource.org/license/MIT)
- Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

## Funding

Matthias König is supported by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151 "QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach)" by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA).

Matthias König was supported by the Federal Ministry of Education and Research (BMBF, Germany) within the research network Systems Medicine of the Liver (**LiSyM**, grant number 031L0054).
