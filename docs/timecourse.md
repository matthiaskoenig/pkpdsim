# Timecourses in parts

A model of a tissue, e.g. a zonated liver coupled to its perfusion, is often simulated directly with [libroadrunner](https://libroadrunner.org) instead of the simulation experiments of sbmlsim, and in parts: a part sets values of the model at its start (a perfusion, a substrate, a switch) and runs for its duration, the next part continues from the state in which the part before ends. `pkpdsim.timecourse` holds the one runner of such simulations, which the tissue and perfusion models of pkdb_models share.

- `SimPart(tend, changes, steps)`: the duration of a part in the time unit of the model, the values it sets at its start by the selection of roadrunner (`"[S_ext]"` for a concentration, `"S_ext"` for an amount or a parameter) and the number of its steps (a part has `steps + 1` points).
- `Simulation(sid, title, parts)`: the parts which run one after the other, with an id and a title for the results and the figures.
- `run_simulation(r, sim)`: resets the model (`resetAll`), runs the parts and returns one `pandas.DataFrame` with the `timeCourseSelections` of `r` as columns, whose time continues over the parts; the last point of a part and the first point of the next have the same time (the state before and after the changes).

The script which runs the model loads it, sets its tolerances and selects what the timecourse holds.

## Usage

A substance `S` flows in at the rate `k_in` and is degraded; the inflow stops after 10 minutes:

```python
from pathlib import Path

import roadrunner
from sbmlutils.factory import (
    Compartment,
    Model,
    Parameter,
    Reaction,
    Species,
    create_model,
)

from pkpdsim.timecourse import SimPart, Simulation, run_simulation

model = Model(
    sid="inflow",
    compartments=[Compartment("V", 1.0)],
    species=[Species("S", "V", initialConcentration=0.0)],
    parameters=[Parameter("k_in", 1.0, constant=False), Parameter("k_out", 0.5)],
    reactions=[
        Reaction("IN", equation="-> S", formula="k_in"),
        Reaction("OUT", equation="S ->", formula="k_out * S"),
    ],
)
create_model(model, filepath=Path("inflow.xml"), validate=False)
r = roadrunner.RoadRunner("inflow.xml")
r.timeCourseSelections = ["time", "[S]"]

sim = Simulation(
    sid="inflow_off",
    title="inflow, then no inflow",
    parts=[
        SimPart(tend=10.0, changes={}, steps=10),
        SimPart(tend=5.0, changes={"k_in": 0.0}, steps=5),
    ],
)
df = run_simulation(r, sim)
print(df.iloc[[0, 10, 11, 16]])
```

```text
    time       [S]
0    0.0  0.000000
10  10.0  1.986522
0   10.0  1.986522
5   15.0  0.163064
```
