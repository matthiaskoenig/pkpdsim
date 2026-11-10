"""The SBML writing kit, the layer below the libraries of pkpdsim.

The PBPK library (`pkpdsim.pbpk`), the zonation engine (`pkpdsim.zonation`)
and the machine perfusion library (`pkpdsim.perfusion`) write SBML with the
same pieces: the units of a model (`pkpdsim.sbml.units`), the numbers with
their provenance (`pkpdsim.sbml.values`), the parameters and the observers
these numbers become (`pkpdsim.sbml.parameters`), the initial assignment of a
constant and the taxon of a build (`pkpdsim.sbml.annotation`) and the
container which collects the objects of a part of a model
(`pkpdsim.sbml.objects`).

The kit is the lower layer: it imports `pkpdsim.database` (the parameter
database, for the citation of a source) and no library of pkpdsim, so the
libraries share these pieces without importing each other. The layering is
pinned by `tests/zonation/test_limits.py`, which walks the four packages and
reads the imports of every module.

`import pkpdsim.sbml` imports none of the modules.
"""
