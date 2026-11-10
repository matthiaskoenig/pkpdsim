# pkpdsim: PBPK/PD models in SBML using python
[![GitHub Actions CI/CD Status](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/matthiaskoenig/pkpdsim/actions/workflows/ci-cd.yml) [![Documentation](https://img.shields.io/badge/docs-pkpdsim-3f51b5.svg)](https://matthiaskoenig.github.io/pkpdsim) [![Version](https://img.shields.io/pypi/v/pkpdsim.svg)](https://pypi.org/project/pkpdsim/) [![Python Versions](https://img.shields.io/pypi/pyversions/pkpdsim.svg)](https://pypi.org/project/pkpdsim/) [![MIT License](https://img.shields.io/pypi/l/pkpdsim.svg)](https://opensource.org/licenses/MIT)

`pkpdsim` is a library for physiologically based pharmacokinetic and pharmacodynamic (PBPK/PD) models in [SBML](https://sbml.org): a standard physiology at its centre, the substances, chemistry, routes and observers of a model on top of it, and the simulation experiments, fits, scans, uncertainty and credibility of [sbmlsim](https://github.com/matthiaskoenig/sbmlsim) on top of those. The source code is available from [https://github.com/matthiaskoenig/pkpdsim](https://github.com/matthiaskoenig/pkpdsim).

A model package writes one definition of what is specific to its substances: the substances with their mechanism flags, the chemistry of its organs, its choices and its PD components. The library adds the body and writes the SBML with the factory of [sbmlutils](https://github.com/matthiaskoenig/sbmlutils), a comp model with one submodel per organ and the flat model which is simulated. The [PBPK library](pbpk/index.md) links every component of the three layers below, and its [tutorial](pbpk/tutorial.md) builds such a model step by step.

<iframe class="pkpd-figure" src="pbpk/img/figures/anatomy.html" title="The interactive anatomical figure of the glimepiride-like model" style="aspect-ratio:1460 / 962" loading="lazy"></iframe>

The anatomical figure of the glimepiride-like model of the [tutorial](pbpk/tutorial.md), generated from its definition: the male and the female build with the organs sized by their volumes, the substances of every organ, the submodels of the liver and the kidney, the routes and the blood sample. Hover an organ, a label, a substance or an icon for its details; [Figures](pbpk/figures.md) describes all figures of a model.

## Background

The models of [pkdb_models](https://github.com/matthiaskoenig/pkdb_models) describe the pharmacokinetics and pharmacodynamics of a substance in the whole body, against the studies curated in [PK-DB](https://pk-db.com). They were written one by one, each with its own body, liver, kidney and intestine model in python with the factory of [sbmlutils](https://github.com/matthiaskoenig/sbmlutils); they are being rebuilt from this library, and each lives in a repository of its own, which imports nothing but itself and released packages (see [Migration](pbpk/migration.md) and [Release of a model repository](release.md)).

## The physiology

Every model is built on one physiology, scaled to the individual of a study and changed by its condition.

- **[Physiology](pbpk/physiology.md)** - the organs with their volumes, blood flows and blood (ICRP 89), the topology of the circulation, and everything which follows from the covariates `BW`, `HEIGHT` and `AGE`; [species](pbpk/species.md) builds the same definition on the rat, the mouse or the pig.
- **[Virtual populations](pbpk/populations.md)** - how an individual is drawn (sex, age, body weight, height, organ function) and how a population is simulated.
- **[Impairment](pbpk/impairment.md)** - liver cirrhosis, chronic kidney disease and heart failure as parameter vectors on the physiology, by the clinical [classes](pbpk/classes.md) the studies report.
- **[Pregnancy](pbpk/pregnancy.md)**, **[circadian rhythms](pbpk/circadian.md)** and **[adipose tissue](pbpk/adipose.md)** - the physiology over the gestational age with the fetoplacental unit, the clock time in a model, and the fat of an individual with its depots.
- **[Standard organ](pbpk/organ.md)** - the four compartments every organ is built from, and the [variants](pbpk/variants.md) within its boundary: the zonated liver, the nephron kidney, the intestines.
- **[Zonated liver](zonation.md)** and **[machine perfusion](perfusion.md)** - the sinusoid of `n` zones with flow and diffusion, and the circuit of an ex vivo perfusion with one mounted organ.

## The building blocks

On top of the physiology a model puts what is specific to its substances, and the library writes the SBML for it.

- **[Substances](pbpk/substances.md)** and **[chemistry](pbpk/chemistry.md)** - the mechanism flags of a substance (binding, erythrocytes, distribution, excretion, inhibition, induction) and the processes of an organ with the proteins which carry them.
- **[Intestine and routes](pbpk/intestine.md)** - how a dose is applied, the lumen segments with their dissolution and absorption, the mucosa, the bile and the enterohepatic circulation.
- **[Observers](pbpk/observers.md)** - what a model is read at: the arm vein, the hepatic vein, the free and the blood concentrations, the excretion rates and the observers of a definition.
- **[Annotation](pbpk/annotation.md)**, **[templates and annotations](templates.md)** and **[figures](pbpk/figures.md)** - the SBO, PBPKO, UBERON, NCIT, UniProt and ChEBI terms of every created element, the creators and the terms of use of a model, and the figures generated from its definition: the anatomical figure, the schematic, the publication figure, the bodies of individuals and the pictograms of the impaired organs.
- **[Parameter database](database.md)** - the numbers the libraries are built from, each with the publication it is taken from; **[initial values](initial_values.md)** put a substance into every tissue at the start.

## The analysis

The models are simulated, fitted and assessed with sbmlsim; these pages are its capabilities in the PBPK context.

- **[Experiments and their metadata](experiments.md)** - a published study as code (its data, its dosing protocol as simulations, its fit mappings) with the shared vocabulary which describes it and the filters which select the mappings of a fit by it.
- **[Scans and sensitivity analyses](scans.md)** - the designs which run a model over many parameter values, and the sensitivity of its PK observables.
- **[Uncertainty of the fitted parameters](uncertainty.md)** - adopted fits whose values the factory takes, their confidence intervals as distrib uncertainties of the SBML parameters, and the run of the experiments over draws of the fitted parameters with bands (`run_simulation <package> --uncertainty`).
- **[Context of use and credibility](credibility.md)** - what a model is credible for, which evidence supports each use and when the use is accepted.
- **[Models and their commands](models.md)** - the model object of a package, the [commands](commands.md) which work with every model, the tests of a model as functions and the verification of a change (`pkpdsim-verify`).
- **[Study data](data.md)**, **[timecourses in parts](timecourse.md)**, **[licenses](licenses.md)** and the **[release of a model repository](release.md)** - the studies of PK-DB a model keeps in its own folder, a model simulated directly with roadrunner in parts, the license texts of a COMBINE archive and the complete release of the repository of a model.

## Quickstart

A definition of what is specific to the substance, and the library writes the whole-body model: the organs with their volumes and blood flows, the oral route, the kidney which clears the drug and the observers it is read at. The numbers are illustrative.

```python
from pathlib import Path
from tempfile import mkdtemp

from sbmlsim.model.model_roadrunner import load_roadrunner

from pkpdsim.pbpk import Binding, Excretion, Kp, PBPKDefinition, Route, Substance, build

# what is specific to the substance: an oral drug which the kidney clears
definition = PBPKDefinition(
    sid="drg_body",
    substances=[
        Substance(
            "drg",
            name="drug",
            mr=450.0,
            binding=Binding(fu=0.1, protein="albumin"),
            distribution=Kp(1.5, ftissue=5.0),
            excretion=Excretion(urine=True),
        )
    ],
    routes=[Route.oral("drg")],
)

# the library adds the physiology, the organs, the routes and the observers
model = build(definition, Path(mkdtemp()))

r = load_roadrunner(model.flat_path)
r.timeCourseSelections = ["time", "Carm_drg"]
r["PODOSE_drg"] = 10.0  # mg
s = r.simulate(0, 24 * 60, 1441)  # min

carm = 1000 * s["Carm_drg"]  # µM
print(f"BW {r['BW']:.0f} kg, Vli {r['Vli']:.2f} l, Qh {r['Qh']:.2f} l/min")
print(f"Cmax {carm.max():.3f} µM, C(24 h) {carm[-1]:.3f} µM")
```

```text
BW 73 kg, Vli 2.33 l, Qh 1.79 l/min
Cmax 0.445 µM, C(24 h) 0.169 µM
```

Continue with [Installation](installation.md) and the [tutorial](pbpk/tutorial.md), which builds such a model step by step and adds what this example leaves out: the hepatic chemistry with its metabolites, the individual of a study, the impairment of an organ and the figures of the model.

## How to cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23297325.svg)](https://doi.org/10.5281/zenodo.23297325)

If you use `pkpdsim` please cite the archived software on [Zenodo](https://doi.org/10.5281/zenodo.23297325):

> König, M. (2026). *pkpdsim: building blocks of physiologically based pharmacokinetic and pharmacodynamic models in SBML* (Version 0.22.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23297325

## License

- Source Code: [MIT](https://opensource.org/license/MIT)
- Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

## Funding

Matthias König is supported by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151 "QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach)" by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA).

Matthias König was supported by the Federal Ministry of Education and Research (BMBF, Germany) within the research network Systems Medicine of the Liver (**LiSyM**, grant number 031L0054).
