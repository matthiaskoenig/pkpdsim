# pkpdsim: PBPK/PD models in SBML using python
[![GitHub Actions CI/CD Status](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml)
[![Documentation](https://img.shields.io/badge/docs-pkpdsim-3f51b5.svg)](https://matthiaskoenig.github.io/pkpdsim)
[![Version](https://img.shields.io/pypi/v/pkpdsim.svg)](https://pypi.org/project/pkpdsim/)
[![Python Versions](https://img.shields.io/pypi/pyversions/pkpdsim.svg)](https://pypi.org/project/pkpdsim/)
[![MIT License](https://img.shields.io/pypi/l/pkpdsim.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23297325.svg)](https://doi.org/10.5281/zenodo.23297325)

`pkpdsim` builds and analyses physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in [SBML](https://sbml.org). A whole-body model of a substance is composed from a short definition in python: its substances, the chemistry of its organs (metabolism, transport, custom reactions), its routes, organ variants and PD components. The library adds the physiology, the organs, the impaired organ function, the observers and the annotation of every element, writes the SBML with the factory of [sbmlutils](https://github.com/matthiaskoenig/sbmlutils), and simulates, fits and assesses the model with [sbmlsim](https://github.com/matthiaskoenig/sbmlsim).

[![The anatomical figure of a PBPK model of glimepiride built with pkpdsim](https://raw.githubusercontent.com/matthiaskoenig/pkpdsim/develop/docs/pbpk/img/figures/anatomy.svg)](https://matthiaskoenig.github.io/pkpdsim/pbpk/figures/)

*The anatomical figure of the glimepiride-like model of the [tutorial](https://matthiaskoenig.github.io/pkpdsim/pbpk/tutorial/), generated from its definition: the male and the female build with the organs sized by their volumes, the substances of every organ, the submodels of the liver and the kidney, the routes and the blood sample. The [interactive figure](https://matthiaskoenig.github.io/pkpdsim/pbpk/figures/) shows the details of every element on hover.*

## Features

- **[Model building](https://matthiaskoenig.github.io/pkpdsim/pbpk/)** - `build(definition, output_dir)` validates a definition and composes the whole-body model: a top model with the physiology, the blood pools, the substances, the sinks, the dose states and the observers, and one comp submodel per organ, written as the comp model `<sid>.xml` and the flat model which is simulated. Every port and every replacement is generated; a model package never writes a `ReplacedElement` and restates no equation of the library.
- **[Organs and their variants](https://matthiaskoenig.github.io/pkpdsim/pbpk/organ/)** - the 4-compartment standard organ of PK-Sim (plasma, blood cells, interstitial and intracellular space), of which only the parts a substance needs are built, and the [alternative components](https://matthiaskoenig.github.io/pkpdsim/pbpk/variants/) behind the same boundary: the zonated liver from periportal to pericentral, the nephron kidney with filtration, secretion and reabsorption, and the [segmented intestine](https://matthiaskoenig.github.io/pkpdsim/pbpk/intestine/) with dissolution, regional absorption and bile. Swapping a variant changes nothing outside the organ.
- **[Substances and chemistry](https://matthiaskoenig.github.io/pkpdsim/pbpk/substances/)** - a substance is a molecular weight and the mechanism flags which decide the equations built for it (binding to plasma proteins, erythrocytes, distribution, an endogenous baseline, a [slow pool](https://matthiaskoenig.github.io/pkpdsim/pbpk/substances/)); its [chemistry](https://matthiaskoenig.github.io/pkpdsim/pbpk/chemistry/) is written once per organ and placed into the parts of whatever variant the organ uses. The steady state of a declared process is generated, so an endogenous substance starts stationary.
- **[Physiology and individuals](https://matthiaskoenig.github.io/pkpdsim/pbpk/physiology/)** - the reference individual of ICRP 89 per sex with the covariates of an individual (`BW`, `HEIGHT`, `AGE`), from which the organ masses and volumes, the blood volume, the cardiac output, the blood flows and the GFR follow, from the newborn to 100 years, and [virtual populations](https://matthiaskoenig.github.io/pkpdsim/pbpk/populations/) drawn over age from the distributions of NHANES.
- **[Special populations and conditions](https://matthiaskoenig.github.io/pkpdsim/pbpk/impairment/)** - impaired organ function as parameter vectors on the physiology, never a structurally different model, with the [clinical classes](https://matthiaskoenig.github.io/pkpdsim/pbpk/classes/) the studies report (Child-Pugh, the GFR categories of CKD, NYHA), [pregnancy](https://matthiaskoenig.github.io/pkpdsim/pbpk/pregnancy/) over the gestational age with the fetoplacental unit, [obesity and adipose tissue](https://matthiaskoenig.github.io/pkpdsim/pbpk/adipose/), [circadian rhythms](https://matthiaskoenig.github.io/pkpdsim/pbpk/circadian/) for chronopharmacokinetics and the [species](https://matthiaskoenig.github.io/pkpdsim/pbpk/species/) of a laboratory animal (rat, mouse, pig).
- **[Observers](https://matthiaskoenig.github.io/pkpdsim/pbpk/observers/)** - the sampling sites of a study as assignment rules which conserve mass exactly: the arm vein of a blood sample, the amount in urine, the readouts of the intestine (the absorbed fraction and the rates of a dose), the changes against the start, and the observers a definition adds for its own readouts and PD components.
- **[Annotation](https://matthiaskoenig.github.io/pkpdsim/pbpk/annotation/)** - every compartment, species, parameter and reaction carries an SBO term and the ontology terms of what it is and where it is (FMA, BTO, NCIT, UBERON, GO, PBPKO), and the [creators and terms of use](https://matthiaskoenig.github.io/pkpdsim/templates/) of a model travel with it, with the [license texts](https://matthiaskoenig.github.io/pkpdsim/licenses/) of its COMBINE archive.
- **[Parameter database](https://matthiaskoenig.github.io/pkpdsim/database/)** - every number the libraries are built from with the publication it is taken from, its unit, its uncertainty and a remark, as tables next to the data of the studies they are integrated from; a number a library chooses where the literature gives none says so. The uncertainty is written into the SBML as a `distrib` uncertainty of the parameter.
- **[Figures](https://matthiaskoenig.github.io/pkpdsim/pbpk/figures/)** - the anatomical figure of the male and the female build, the schematic of the circulation, the one-page figure of a publication, the bodies of individuals and the pictograms of the impaired organs, all generated from the definition, pure python, SVG and deterministic, with an interactive page of each for the web.
- **[Simulation experiments](https://matthiaskoenig.github.io/pkpdsim/experiments/)** - the base class of the experiments of a model, a published study as code: the [study data](https://matthiaskoenig.github.io/pkpdsim/data/) of [PK-DB](https://pk-db.com) collected into the model, the dosing protocol as simulations, the data paired with the observables, the [initial values](https://matthiaskoenig.github.io/pkpdsim/initial_values/) which remove the distribution kinetics, the [timecourses in parts](https://matthiaskoenig.github.io/pkpdsim/timecourse/), the presets of the impaired organ function and the vocabulary which describes a fit mapping.
- **[Fits, identifiability and uncertainty](https://matthiaskoenig.github.io/pkpdsim/uncertainty/)** - the fit of a parameter subset with its identifiability (profile likelihood, Fisher information), the adopted fits whose values the factory takes and whose confidence intervals become the `distrib` uncertainties of the SBML parameters, and the run of every experiment over draws of the fitted parameters with uncertainty bands.
- **[Scans and sensitivity analyses](https://matthiaskoenig.github.io/pkpdsim/scans/)** - the PK observables of a non-compartmental analysis, the parameters of a model of the library and the local, sampling, FAST and Sobol analyses on the runner of sbmlsim, per sex and from the drug-free steady state of a condition.
- **[Context of use and credibility](https://matthiaskoenig.github.io/pkpdsim/credibility/)** - what a model is credible for, declared and assessed in the frame of ASME V&V 40: the contexts of use `C1` to `C5` with their quantities of interest, domains and evidence, the ratios of a special population, the temporal validation, the verification, the predictive check, the lock, the report and the scorecard.
- **[Machine perfusion](https://matthiaskoenig.github.io/pkpdsim/perfusion/)** - the SBML of an ex vivo machine perfusion circuit with one mounted organ from a definition, as the PBPK library builds a body: the reservoir, the tubing, the pump, the oxygenator, the temperature, the perfusate and the organ as a chain of zones with the cellular core, the liver with its bile and the kidney with its filtration and urine.
- **[Commands](https://matthiaskoenig.github.io/pkpdsim/commands/)** - the commands which work with every [model package](https://matthiaskoenig.github.io/pkpdsim/models/), the dotted name of the package as their first argument: `run_simulation`, `run_fit`, `run_report`, `run_identifiability` and `run_credibility`, the tests of a model as functions, the verification of a change, also across versions of sbmlsim (`pkpdsim-verify`), and the complete [release](https://matthiaskoenig.github.io/pkpdsim/release/) of the repository of a model (`pkpd_release`).

The documentation is available at [https://matthiaskoenig.github.io/pkpdsim](https://matthiaskoenig.github.io/pkpdsim); the [tutorial](https://matthiaskoenig.github.io/pkpdsim/pbpk/tutorial/) builds a PBPK model from its definition step by step.

If you have any questions or issues please [open an issue](https://github.com/matthiaskoenig/pkpdsim/issues).

## How to cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23297325.svg)](https://doi.org/10.5281/zenodo.23297325)

If you use `pkpdsim` please cite the archived software on [Zenodo](https://doi.org/10.5281/zenodo.23297325):

> König, M. (2026). *pkpdsim: building blocks of physiologically based pharmacokinetic and pharmacodynamic models in SBML* (Version 0.21.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23297325

```bibtex
@software{konig_pkpdsim,
  author    = {König, Matthias},
  title     = {pkpdsim: building blocks of physiologically based pharmacokinetic and pharmacodynamic models in SBML},
  year      = {2026},
  month     = nov,
  version   = {0.21.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23297325},
  url       = {https://doi.org/10.5281/zenodo.23297325},
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
