# PBPK library

`pkpdsim.pbpk` builds the whole-body PBPK model of a substance from a short definition. A model package writes only what is specific to its substances: the substances with their mechanism flags, the chemistry of its organs (metabolism, transport), its choices (organ variants, routes, population, sex) and its PD components. The library adds a standard human physiology (ICRP 89), a standard organ for every organ, the alternative organ variants, the routes and sinks, the impaired organ function, the observers and the annotation of every element, and writes the SBML: a comp model with one submodel per organ and the flat model which is simulated.

Start with the [tutorial](tutorial.md), which writes a definition and simulates an oral dose step by step, and read [the build](build.md) next: how the parts are composed, what the library generates, how a build per sex works, where every number comes from and which ids a model has.

The pages below are the components of the library, in the order in which it is layered: the physiology, the building blocks a model puts on top of it and the analysis which runs on the models. Some of them stand in the user guide of the navigation; the tier says which layer a page belongs to, not where it is filed.

## Physiology

Every model is built on one physiology, scaled to the individual of a study and changed by its condition.

- [Physiology](physiology.md) - the organs with their volumes, blood flows and blood (ICRP 89), the topology of the circulation, and what follows from the covariates `BW`, `HEIGHT` and `AGE`.
- [Species](species.md) - what the same definition becomes on the physiology of the rat, the mouse or the pig instead of the human.
- [Virtual populations](populations.md) - how an individual is drawn (sex, age, body weight, height, organ function) and how a population is simulated.
- [Impairment](impairment.md) - the parameter vectors with which liver cirrhosis, chronic kidney disease and heart failure change the healthy individual.
- [Impairment classes](classes.md) - which value of a model parameter represents a Child-Pugh class, a CKD category or an NYHA class, and what each condition does to a drug.
- [Pregnancy](pregnancy.md) - the physiology over the gestational age and the fetoplacental unit, which gives the fetal exposure.
- [Circadian rhythms](circadian.md) - the clock time in a model, so that a dose in the morning and the same dose in the evening differ.
- [Adipose tissue and obesity](adipose.md) - the fat of an individual as the remainder of its body weight, its visceral and subcutaneous depots and a body weight which changes.
- [Standard organ](organ.md) - the four compartments every organ is built from, which of them a substance needs, and the validity of the kidney outflow.
- [Organ variants](variants.md) - the alternative component of an organ within the boundary of the standard organ: the zonated liver, the nephron kidney, the intestines.
- [Zonated liver](../zonation.md) - the sinusoid of `n` zones the zonation engine writes, its flow, diffusion and profiles, and the figures of a zonated result.
- [Machine perfusion](../perfusion.md) - the circuit of an ex vivo perfusion with one mounted organ, built from a definition as a body is.

## Building blocks

On top of the physiology a model puts what is specific to its substances, and the library writes the SBML for it.

- [Substances](substances.md) - the mechanism flags of a substance (binding, erythrocytes, distribution, excretion, inhibition, induction) and the equations each of them builds.
- [Chemistry](chemistry.md) - the processes of an organ (metabolism, active transport, custom reactions), their kinetics and the proteins which carry them.
- [Intestine and routes](intestine.md) - how a dose is applied, the lumen segments with their dissolution and absorption, the mucosa, the bile and the enterohepatic circulation.
- [Observers](observers.md) - what a model is read at: the arm vein, the hepatic vein, the free and the blood concentrations, the excretion rates and the observers of a definition.
- [Annotation](annotation.md) - the SBO, PBPKO, UBERON, NCIT, UniProt and ChEBI terms which every created element carries.
- [Figures](figures.md) - the anatomical figure, the schematic, the publication figure and the interactive pages, all generated from the definition.
- [Migration](migration.md) - how a hand-written model of pkdb_models moves onto the library, with glimepiride as the example.
- [Templates and annotations](../templates.md) - the creators of a model, its terms of use and the annotations of the compartments which recur in the models.
- [Parameter database](../database.md) - the numbers the libraries are built from, each with the publication it is taken from.
- [Initial values](../initial_values.md) - the same initial concentration or amount of a substance in every tissue, which removes the distribution kinetics.

## Analysis

The models are simulated, fitted and assessed with sbmlsim; these pages are its capabilities in the PBPK context.

- [Experiments and their metadata](../experiments.md) - a published study as code: its data, its dosing protocol as simulations, its fit mappings and the vocabulary which describes it.
- [Scans and sensitivity analyses](../scans.md) - the designs which run a model over many parameter values, and the sensitivity of its PK observables.
- [Uncertainty of the fitted parameters](../uncertainty.md) - the adopted fits, their confidence intervals in the SBML and the bands of a run over their draws.
- [Context of use and credibility](../credibility.md) - what a model is credible for, which evidence supports each use and when the use is accepted.
- [Models and their commands](../models.md) - how a model package describes itself once as a `PKPDModel`, and what the commands take from it.
- [Commands](../commands.md) - the six commands which work with every model package, and their options.
- [Study data](../data.md) - the studies curated in PK-DB and the subset of them a model keeps in its own folder.
- [Timecourses in parts](../timecourse.md) - a model simulated directly with roadrunner in parts, one part after the other in one timecourse.
- [Release of a model repository](../release.md) - the complete release of the repository of a model with `pkpd_release`.
- [Licenses](../licenses.md) - the license texts of the COMBINE archive of a model.
