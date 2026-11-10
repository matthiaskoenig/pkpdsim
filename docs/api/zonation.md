# zonation

The zonation engine, see [Zonated liver](../zonation.md): the settings of a zonated liver, the geometry of its sinusoid, the gradients along the lobule, the transport of a substance, the flow through the sinusoids, the cell model per zone, the observers, the build, the id scheme, the figures and the bridge from a simulation result of sbmlsim to them. The engine writes plain SBML with the [sbml](sbml.md) kit and imports neither the PBPK library nor the perfusion library (of pkpdsim it reads `pkpdsim.sbml`, `pkpdsim.database.tables`, `pkpdsim.species` and `pkpdsim.zones`): a caller (a zonated liver of a body model, a perfused organ) hands over its boundary and its cell model.

::: pkpdsim.zonation
    options:
      members: false

## settings

::: pkpdsim.zonation.settings

## geometry

::: pkpdsim.zonation.geometry

## profiles

::: pkpdsim.zonation.profiles

## transport

::: pkpdsim.zonation.transport

## flow

::: pkpdsim.zonation.flow

## composition

::: pkpdsim.zonation.composition

## observers

::: pkpdsim.zonation.observers

## builder

::: pkpdsim.zonation.builder

## ids

::: pkpdsim.zonation.ids

## figures

::: pkpdsim.zonation.figures

## results

::: pkpdsim.zonation.results

## render

::: pkpdsim.zonation.render
