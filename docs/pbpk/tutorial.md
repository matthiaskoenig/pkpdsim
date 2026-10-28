# Tutorial

This tutorial builds a glimepiride-like model step by step: the substances, their chemistry, the routes, the build, and an oral dose simulated with the arm vein concentration read from the result. Glimepiride is an oral antidiabetic which is almost completely bound to albumin and metabolised by CYP2C9 in the liver into the hydroxy metabolite M1, which cytosolic enzymes oxidise further into the carboxy metabolite M2; both metabolites are excreted into the urine. The numbers are illustrative, not the fitted values of the glimepiride model of pkdb_models. Every block of this page runs, in order (`tests/docs/test_snippets.py`).

## 1. Substances

A substance is a short id, a name, the molecular weight and its mechanism flags. Glimepiride binds to albumin (fu 0.005) and enters the tissues with one Kp (the merged tissue with the whole-body exchange `ftissue`, the default). With fu 0.005 the exchange of the liver which `ftissue` gives would limit the hepatic metabolism, so the liver gets its own exchange (`ps`). The metabolites stay in the blood (no distribution) and are filtered into the urine.

```python
from pkpdsim.pbpk import Binding, Excretion, Kp, Organ, Substance

glimepiride = Substance(
    "gli",
    name="glimepiride",
    mr=490.6,
    binding=Binding(fu=0.005, protein="albumin"),
    # ps per litre of liver: 42.9 * 2.33 l = 100 l/min for the whole liver
    distribution=Kp(0.25, ftissue=50.0, ps={Organ.LIVER: 42.9}),
)
m1 = Substance(
    "gli_m1",
    name="hydroxy glimepiride (M1)",
    mr=506.6,
    excretion=Excretion(urine=True),
)
m2 = Substance(
    "gli_m2",
    name="carboxy glimepiride (M2)",
    mr=520.6,
    excretion=Excretion(urine=True),
)
print([s.sid for s in (glimepiride, m1, m2)])
```

```text
['gli', 'gli_m1', 'gli_m2']
```

The metabolites take the prefix of the parent (`gli_m1`), so their ids stay unique when models are coupled. See [Substances](substances.md) for every flag.

## 2. Chemistry

The chemistry is defined once per organ, the organ places it. CYP2C9 converts glimepiride into M1 with Michaelis-Menten kinetics on the unbound concentration in the liver tissue (Vmax per litre of liver), and M1 becomes M2 by first order, a clearance per litre of liver (1/min, the default of `FirstOrder`), so the clearance of the liver follows its volume as the per-litre constant of the hand-written model did.

```python
from pkpdsim.pbpk import MM, FirstOrder, Metabolism

chemistry = [
    Metabolism(
        "GLI2M1",
        Organ.LIVER,
        "gli",
        products={"gli_m1": 1.0},
        kinetics=MM(vmax=0.04, km=0.01),
        protein="CYP2C9",
    ),
    Metabolism(
        "M12M2",
        Organ.LIVER,
        "gli_m1",
        products={"gli_m2": 1.0},
        kinetics=FirstOrder(k=0.2),
    ),
]
```

The process ids are those of the hand-written model (`GLI2M1`), so the fitted Vmax keeps its flat id `LI__GLI2M1_Vmax`. The protein brings the activity factor `LI__f_cyp2c9`, which the experiments of the CYP2C9 genotypes set. See [Chemistry](chemistry.md).

## 3. Routes and the definition

Glimepiride is given as tablet, with the first order dissolution of the library default. The definition collects the substances, the chemistry and the routes; everything else is a default: the standard organ in every organ, the segmented intestine, the male ICRP 89 reference individual.

```python
from pkpdsim.pbpk import PBPKDefinition, Route

definition = PBPKDefinition(
    sid="glimepiride_body",
    substances=[glimepiride, m1, m2],
    chemistry=chemistry,
    routes=[Route.oral("gli")],
)
definition.validate()
print(definition.variants.organs, definition.sex, definition.population)
```

```text
{} male european_icrp
```

`validate` names every problem at once (an unknown substance in a process, a fraction outside its interval, an id which collides with the id scheme); `build` calls it as well.

## 4. Build

`build` composes the comp model (one submodel per organ), writes it and flattens it; the flat model is the one which is simulated. In a model package this is the factory, `models/factory.py`, which writes into `models/output/`.

```python
from pathlib import Path

from pkpdsim.pbpk import build

output = Path("models")
output.mkdir()
result = build(definition, output)
print(result.comp_path.name, result.flat_path.name)
```

```text
glimepiride_body.xml glimepiride_body_flat.xml
```

## 5. Simulate an oral dose

The flat model is loaded as sbmlsim loads it (`load_roadrunner`, the production path of the experiments). The dose is the dose state `PODOSE_gli` in mg; the arm vein plasma concentration `Carm_gli` (mM) is the observable which plasma data map to.

```python
import numpy as np
from sbmlsim.model.model_roadrunner import load_roadrunner

r = load_roadrunner(result.flat_path)
r.timeCourseSelections = ["time", "Carm_gli", "Aurine_gli_m1", "Aurine_gli_m2"]
r["PODOSE_gli"] = 4.0  # mg
s = r.simulate(0, 24 * 60, 1441)  # min

carm = s["Carm_gli"] * 1000  # µM
i = int(np.argmax(carm))
dose = 4.0 / 490.6  # mmole
urine = (s["Aurine_gli_m1"][-1] + s["Aurine_gli_m2"][-1]) / dose
print(f"Cmax {carm[i]:.3f} µM at {s['time'][i] / 60:.2f} h")
print(f"C(24 h) {carm[-1]:.4f} µM, M1 + M2 in urine {100 * urine:.0f} % of the dose")
```

```text
Cmax 0.310 µM at 1.65 h
C(24 h) 0.0292 µM, M1 + M2 in urine 86 % of the dose
```

The other observables are read the same way: `Cve_gli` (the venous plasma), `Chv_gli` (the hepatic vein), `Carm_free_gli` (unbound), `LI__Cti_gli` (the liver tissue), `GU__Clumen_du_gli` (the duodenal lumen), `Afeces_gli` (the unabsorbed and excreted drug).

## 6. Change the individual

The covariates are initial values: an experiment of sbmlsim changes them before the simulation, and every organ volume and flow follows. With plain roadrunner the initial value is set and the model reset:

```python
r.resetAll()
for bw in (73.0, 100.0):
    r.setValue("init(BW)", bw)
    r.reset()
    print(
        f"BW {r['BW']:.0f} kg: Vad {r['Vad']:.1f} l, Vli {r['Vli']:.2f} l, CO {r['CO']:.2f} l/min"
    )
```

```text
BW 73 kg: Vad 19.0 l, Vli 2.33 l, CO 6.50 l/min
BW 100 kg: Vad 36.2 l, Vli 2.77 l, CO 7.72 l/min
```

An individual of 100 kg and the reference height has more adipose, a larger liver and a higher cardiac output: the lean tissues and the cardiac output scale with the fat-free mass, adipose is the rest of the body weight (at 73 kg more than the ICRP adipose mass, since it takes the mass of the tissues which are no organ of the model). The severities of the impairment work the same way (`init(CP_severity)`, `init(f_renal_function)`), see [Impairment](impairment.md).

## 7. In a model package

The factory of the model package builds the definition with the fitted values from the adopted fits of the model (`FITS = MODEL.fit_records()`), so the uncertainties of the fits reach the SBML; the observers the experiments read which the library does not compute (the sum of the metabolites in urine) and the creators of the model are part of the definition, and so are the constants of the model which no field of the library holds (`parameters`, see [Substances](substances.md#parameters-of-the-definition)):

```python
# not executed: the model package and its adopted fits; m1, m2, chemistry of above
from pathlib import Path
from typing import Any

from pkpdsim.pbpk import (
    Binding,
    Kp,
    Observer,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    U,
    create_archive,
    fitted,
)
from pkpdsim.templates import creators

from pkdb_models.models.glimepiride import MODEL_BASE_PATH
from pkdb_models.models.glimepiride.models.fits import FITS

glimepiride = Substance(
    "gli",
    name="glimepiride",
    mr=490.6,
    binding=Binding(fu=0.005, protein="albumin"),
    distribution=Kp(
        fitted(FITS.parameter("Kp_gli", unit=U.dimensionless)),
        ftissue=fitted(FITS.parameter("ftissue_gli", unit=U.l_per_min)),
        # per litre of liver, 42.9 * 2.33 l = 100 l/min for the whole liver
        ps={Organ.LIVER: 42.9},
    ),
)
DEFINITION = PBPKDefinition(
    sid="glimepiride_body",
    substances=[glimepiride, m1, m2],
    chemistry=chemistry,
    routes=[Route.oral("gli")],
    observers=[
        Observer(
            "Aurine_gli_m1_m2",
            "Aurine_gli_m1 + Aurine_gli_m2",
            U.mmole,
            "M1 and M2 in urine [mmole]",
        )
    ],
    creators=creators,
)


def create_glimepiride_models(model_output_dir: Path) -> dict[str, dict[str, Any]]:
    """Build the glimepiride model and its COMBINE archive."""
    return create_archive(
        DEFINITION,
        model_output_dir,
        model="glimepiride",
        copyright="2024-2026 Matthias König",
        readme=MODEL_BASE_PATH.parent / "README.md",
    )
```

`create_archive` builds the definition (`build`), writes the differential equations of the flat model as markdown and its antimony, the texts of the licenses (CC BY-SA 4.0 and MIT with the copyright) into the parent of `model_output_dir` and the COMBINE archive `<model>_model.omex` there, with the README, the licenses, the comp model, the flat model as the master entry, the markdown, the antimony, the SBML files of `extra_models` (a submodel which is simulated on its own) and the files of `extra_files` (the ODEs of an extra model as markdown, a figure; a file below the parent of `model_output_dir`, at its path relative to it, in the format of its suffix: `.md`, `.ant`, `.txt`, `.png`, `.svg`, `.csv`, `.tsv`, `.pdf`). It returns the files with their entries by the keys `README`, `cc-by`, `mit`, `<sid>`, `<sid>_flat`, `<sid>_flat_md`, `<sid>_flat_ant`, the stem of every extra model and `<stem>_<suffix>` of every extra file (`fluid_md`).

The experiments then set the doses, covariates and activity factors and read the observables as before; the [migration guide](migration.md) lists the ids which change.
