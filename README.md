# pkpdsim: PBPK/PD models in SBML
[![GitHub Actions CI/CD Status](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml)
[![Documentation](https://img.shields.io/badge/docs-pkpdsim-3f51b5.svg)](https://matthiaskoenig.github.io/pkpdsim)
[![Version](https://img.shields.io/pypi/v/pkpdsim.svg)](https://pypi.org/project/pkpdsim/)
[![Python Versions](https://img.shields.io/pypi/pyversions/pkpdsim.svg)](https://pypi.org/project/pkpdsim/)
[![MIT License](https://img.shields.io/pypi/l/pkpdsim.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23297325.svg)](https://doi.org/10.5281/zenodo.23297325)

`pkpdsim` is a library for physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in [SBML](https://sbml.org): a standard physiology at its centre, the substances, chemistry, routes and observers of a model on top of it, and the simulation experiments, fits, scans, uncertainty and credibility of [sbmlsim](https://github.com/matthiaskoenig/sbmlsim) on top of those. A model package writes one definition of what is specific to its substances; the library adds the body and writes the SBML, a comp model with one submodel per organ and the flat model which is simulated. The models of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models) are being rebuilt from it, each in a repository of its own.

## Features

### The physiology

- **the standard physiology** - the organs with their volumes, blood flows and blood of ICRP 89, the topology of the circulation, and everything which follows from the covariates of an individual (`BW`, `HEIGHT`, `AGE`), see the [physiology](https://matthiaskoenig.github.io/pkpdsim/pbpk/physiology/)
- **species and populations** - the same definition on the physiology of the rat, the mouse or the pig, and the individuals a [virtual population](https://matthiaskoenig.github.io/pkpdsim/pbpk/populations/) draws
- **impaired organ function** - liver cirrhosis, chronic kidney disease and heart failure as parameter vectors on the physiology, by the clinical [classes](https://matthiaskoenig.github.io/pkpdsim/pbpk/classes/) the studies report
- **pregnancy, circadian rhythms and obesity** - the physiology over the gestational age with the fetoplacental unit, the clock time in a model, and the fat of an individual with its depots
- **organs and their variants** - the four-compartment [standard organ](https://matthiaskoenig.github.io/pkpdsim/pbpk/organ/) of every organ and the [variants](https://matthiaskoenig.github.io/pkpdsim/pbpk/variants/) within its boundary: the zonated liver, the nephron kidney, the intestines
- **tissue and perfusion** - the [zonated liver](https://matthiaskoenig.github.io/pkpdsim/zonation/) of `n` zones with flow and diffusion along the sinusoid, and the circuit of an ex vivo [machine perfusion](https://matthiaskoenig.github.io/pkpdsim/perfusion/) with one mounted organ

### The building blocks

- **substances and chemistry** - the mechanism flags of a [substance](https://matthiaskoenig.github.io/pkpdsim/pbpk/substances/) (binding, erythrocytes, distribution, excretion, inhibition, induction) and the [processes](https://matthiaskoenig.github.io/pkpdsim/pbpk/chemistry/) of an organ with the proteins which carry them
- **routes and the intestine** - the intravenous, oral tablet and oral solution routes, the lumen segments with dissolution and absorption, the bile and the enterohepatic circulation
- **observers** - what a model is read at: the arm vein, the hepatic vein, the free and the blood concentrations, the excretion rates and the [observers](https://matthiaskoenig.github.io/pkpdsim/pbpk/observers/) of a definition
- **annotation, templates and figures** - the SBO, PBPKO, UBERON, NCIT, UniProt and ChEBI terms of every created element, the creators and terms of use of a model, and the [figures](https://matthiaskoenig.github.io/pkpdsim/pbpk/figures/) generated from its definition
- **the parameter database** - the numbers the libraries are built from, each with the publication it is taken from, see the [parameter database](https://matthiaskoenig.github.io/pkpdsim/database/)

### The analysis

- **experiments and their metadata** - a published study as code (its data, its dosing protocol as simulations, its fit mappings) with the shared vocabulary which describes it and the filters which select the mappings of a fit by it, see [experiments](https://matthiaskoenig.github.io/pkpdsim/experiments/)
- **fits, scans and sensitivity analyses** - the fit settings, the PK observables and the library parameters of a [sensitivity analysis](https://matthiaskoenig.github.io/pkpdsim/scans/), which run on the runner of sbmlsim
- **uncertainty of the fitted parameters** - the adopted fits whose values the factory takes, their confidence intervals as `distrib` uncertainties of the SBML parameters and the bands of a run over their [draws](https://matthiaskoenig.github.io/pkpdsim/uncertainty/)
- **context of use and credibility** - what a model is credible for, which evidence supports each use and when the use is accepted, see [credibility](https://matthiaskoenig.github.io/pkpdsim/credibility/)
- **models and their commands** - the model object of a package and the [commands](https://matthiaskoenig.github.io/pkpdsim/commands/) which work with every model (`run_simulation`, `run_fit`, `run_report`, `run_identifiability`, `run_credibility`), the tests of a model as functions and the verification of a change (`pkpdsim-verify`)
- **release of a model repository** - `pkpd_release` makes the complete release of the repository of a model (models, simulations, fit reports, release notes, docs, Docker image, GitHub release) from its package, with `--no-publish` for a review before the commit

The documentation is available at [https://matthiaskoenig.github.io/pkpdsim](https://matthiaskoenig.github.io/pkpdsim), where the [PBPK library](https://matthiaskoenig.github.io/pkpdsim/pbpk/) links every component of the three layers.

If you have any questions or issues please [open an issue](https://github.com/matthiaskoenig/pkpdsim/issues).

## Quickstart

```python
from pkpdsim.annotations import compartments
from pkpdsim.templates import creators

print(creators[0].givenName, creators[0].familyName)
print([resource for _, resource in compartments["li"]])
```

```text
Matthias König
['fma/FMA:7197', 'bto/BTO:0000759', 'ncit/C12392']
```

The [tutorial](https://matthiaskoenig.github.io/pkpdsim/pbpk/tutorial/) builds a model and simulates an oral dose step by step; every page of the documentation describes its part with a snippet which runs.

## How to cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23297325.svg)](https://doi.org/10.5281/zenodo.23297325)

If you use `pkpdsim` please cite the archived software on [Zenodo](https://doi.org/10.5281/zenodo.23297325):

> König, M. (2026). *pkpdsim: building blocks of physiologically based pharmacokinetic and pharmacodynamic models in SBML* (Version 0.19.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23950939

```bibtex
@software{konig_pkpdsim,
  author    = {König, Matthias},
  title     = {pkpdsim: building blocks of physiologically based pharmacokinetic and pharmacodynamic models in SBML},
  year      = {2026},
  month     = nov,
  version   = {0.19.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23950939},
  url       = {https://doi.org/10.5281/zenodo.23950939},
}
```

## Installation

`pkpdsim` requires python >= 3.14, simulates with sbmlsim 0.18 (`sbmlsim>=0.18.0,<0.19`, installed with it) and is available from [pypi](https://pypi.org/project/pkpdsim/):

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
