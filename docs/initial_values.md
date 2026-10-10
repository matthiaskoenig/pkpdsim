# Initial values

A substance which is given before a study starts, or which the body produces itself, is present in every tissue at the start of a simulation. Setting its initial value in one compartment, e.g. the venous blood, leaves the model to distribute it first, and the simulation starts with distribution kinetics which the study never saw. `pkpdsim.initial_values` sets the same initial value in every tissue at once, so that the simulation starts in (or close to) the steady state of the distribution.

The species of a substance are found by the naming pattern of the ids of the models: `C<tissue>_<substance>` for a concentration and `A<tissue>_<substance>` for an amount, e.g. `Cve_len`, `Cli_plasma_len` and `Aurine_len` for the venous blood, the plasma of the liver and the urine of lenvatinib (`len`). The urine collects what the body excretes and is left alone.

## Usage

A model with a substance `x` in the venous blood, the liver and the urine, simulated with [libroadrunner](https://libroadrunner.org):

```python
from pathlib import Path

import roadrunner
from sbmlutils.factory import Compartment, Model, Species, create_model

from pkpdsim.initial_values import init_concentrations_changes

model = Model(
    sid="distribution",
    compartments=[
        Compartment("Vve", 5.0),
        Compartment("Vli", 1.5),
        Compartment("Vurine", 1.0),
    ],
    species=[
        Species("Cve_x", "Vve", initialConcentration=0.0),
        Species("Cli_plasma_x", "Vli", initialConcentration=0.0),
        Species("Aurine_x", "Vurine", initialAmount=0.0, hasOnlySubstanceUnits=True),
    ],
)
create_model(model, filepath=Path("distribution.xml"), validate=False)
r = roadrunner.RoadRunner("distribution.xml")

changes = init_concentrations_changes(r, "x", 2.0)
print(changes)
```

```text
{'[Cve_x]': 2.0, '[Cli_plasma_x]': 2.0}
```

The keys are selections of roadrunner: `[sid]` the concentration of a species and `sid` its amount (`init_amounts_changes`). The changes are applied to the initial values of the model before the simulation starts: a `Simulation` of [sbmlsim](https://github.com/matthiaskoenig/sbmlsim) takes them as its `preinit_changes`, here they are set in roadrunner directly:

```python
for key, value in changes.items():
    r.setValue(f"init({key})", value)
r.reset()
print(r.getValue("[Cve_x]"), r.getValue("[Cli_plasma_x]"), r.getValue("Aurine_x"))
```

```text
2.0 2.0 0.0
```

`get_species_keys(skey, species_ids)` is the pattern alone; it does not see the species of the submodels of a composed model, only the ones of the top model, which a flattened model has all.
