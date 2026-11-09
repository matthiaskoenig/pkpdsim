# perfusion

The machine perfusion library, see [Machine perfusion](../perfusion.md): the definition of a circuit with its organ, the cellular core, the circuit and the observers, the build, the archive, the figures and the base of the simulation experiments.

::: pkpdsim.perfusion
    options:
      members: false

## definition

::: pkpdsim.perfusion.definition

## core

::: pkpdsim.perfusion.core

## circuit

::: pkpdsim.perfusion.circuit

## observers

::: pkpdsim.perfusion.observers

## builder

::: pkpdsim.perfusion.builder

## experiment

::: pkpdsim.perfusion.experiment

## archive

::: pkpdsim.perfusion.archive

## figure

::: pkpdsim.perfusion.figure

## anatomy

::: pkpdsim.perfusion.anatomy

## ids

::: pkpdsim.perfusion.ids

## units

`PU`, the units of every perfusion model, is the kit's `OrganUnits` under the name this library has always called it by, documented at [sbml](sbml.md#units); `pkpdsim.perfusion.units` re-exports it and keeps the import path of the library. `CoreObjects` of the core above is the kit's `SbmlObjects` the same way.
