"""pkpdsim: physiologically based PK/PD models in SBML.

A library for physiologically based pharmacokinetic and pharmacodynamic
(PBPK/PD) models in SBML: a standard physiology at its centre, the
substances, chemistry, routes and observers of a model on top of it, and the
simulation experiments, fits, scans, uncertainty and credibility of sbmlsim
on top of those. The models of
[pkdb_models](https://github.com/matthiaskoenig/pkdb_models) are built from
it, each in a repository of its own which imports nothing but itself and
released packages.

The physiology and what is built on it are the libraries: `pkpdsim.pbpk`
(the whole-body PBPK model of a substance from a `PBPKDefinition`),
`pkpdsim.zonation` (the zonated liver, a chain of zones with a cell model
each) and `pkpdsim.perfusion` (the circuit of an ex vivo machine perfusion
with one mounted organ). They write their SBML with one kit (`pkpdsim.sbml`)
and take their numbers from one database (`pkpdsim.database`).

Everything a model needs besides them is a module of its own, which it
imports directly: `pkpdsim.model` (the model object of a package) and
`pkpdsim.cli` (its commands), `pkpdsim.experiment` and `pkpdsim.metadata` (a
published study as code and the vocabulary which describes it),
`pkpdsim.fit`, `pkpdsim.sensitivity`, `pkpdsim.credibility` and
`pkpdsim.release` (the analysis of a model and its release),
`pkpdsim.templates` (the creators and the terms of use of a model),
`pkpdsim.annotations` (the annotations of the recurring compartments),
`pkpdsim.data` (the study data of a model from the study folders of PK-DB)
and `pkpdsim.initial_values` (the initial values of a substance in every
tissue). `import pkpdsim` imports none of them. The license texts of the
COMBINE archive of a model are `sbmlsim.licenses`.
"""

__version__ = "0.22.0.dev0"

__all__ = ["__version__"]
