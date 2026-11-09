# sbml

The SBML writing kit, the layer below the libraries of pkpdsim: the units of a model, the numbers with their provenance, the parameters and the observers these numbers become, the initial assignment of a constant with the taxon of a build, and the container which collects the objects of a part of a model. The PBPK library, the zonation engine and the machine perfusion library write SBML with these pieces, so they share them without importing each other; the numbers themselves are the [parameter database](../database.md).

Two names of the kit are the names a library has always called the object by and keep working: `PU` of the machine perfusion library ([perfusion](perfusion.md)) is `pkpdsim.sbml.units.OrganUnits`, and `CoreObjects` of its core is `pkpdsim.sbml.objects.SbmlObjects`.

::: pkpdsim.sbml
    options:
      members: false

## units

::: pkpdsim.sbml.units

## values

::: pkpdsim.sbml.values

## parameters

::: pkpdsim.sbml.parameters

## objects

::: pkpdsim.sbml.objects

## annotation

::: pkpdsim.sbml.annotation
