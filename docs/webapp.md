# Web app of a model

Every model has a web app: one HTML file next to its SBML which shows the readouts of the model and simulates them in the browser of its reader, written by the factory with the models (`run_simulation <package> --action factory`, sbmlsim's `Project.create_models`, see the page Web app of the sbmlsim docs) and shown by the release of the model in its documentation. A model which declares its own panels and app (`panels.py`, `App`) shows them; every other model shows the default app of pkpdsim (`pkpdsim.webapp`, `PKPDModel.default_app`), which is built from what the model and its simulation experiments state, so it needs nothing written by hand.

## The pages

| page | the knobs | from |
| --- | --- | --- |
| Dose | a dose per route and substance the experiments with data dose (`PODOSE_<s>`, `IVDOSE_<s>`, `SOLDOSE_<s>`, ...), at most three, the most studied first; a range from 0 to the largest studied dose with the studied doses as its stops, the most studied dose of the first route as its default | the pre-initialization of the simulations of the experiments with fit mappings |
| Dosing protocol | the doses and a choice of their interval (a single dose, once, twice or three times daily) and the days of the dosing, six days | the same doses, repeated by a dosing protocol of the panel (`sbmlsim.panel.PanelDosing`): the model is not changed, the doses are changes of the simulation |
| Renal impairment, Liver cirrhosis, Heart failure | the classes of the disease as a choice and the first dose | the classes the experiments of the model simulate (`renal_map`, `cirrhosis_map`, `cardiac_map` of the base of its experiments) on the parameter of the model (`f_renal_function` or `KI__f_renal_function`, `CP_severity` or `f_cirrhosis`, `HF_severity` or `f_cardiac_function`); a page for a disease of the subjects of its studies (`Health` of the metadata of the fit mappings) or a scan of its parameter |
| one per other scan | the scanned parameter as a range over the scanned values and the first dose | the `scan_map` of the scan experiments, e.g. the body weight or the activity of a transporter |

Every page shows the readouts as small panels side by side under the knobs it shares: the outputs the experiments with data compare with their studies, the observables of their fit mappings (at most eight, the most compared first), a concentration in uM and an amount in umol. A model of two builds (male and female) names both, and a reader switches between them; the build of each sex has the classes of a renal impairment at its own healthy GFR (`ScaledRenalMap.of_sex`), the panels of its page through `PKPDModel.default_variant_panels`, the override of sbmlsim's `Project.default_variant_panels`. A model whose experiments dose nothing shows its readouts on one page, and one whose experiments compare nothing with data the time courses of its species.

## The body figure of an app

A model which declares its own app can make its interactive body figure the map of its panels: a reader clicks an organ of the drawing and the timecourses the model reports there appear beside the plots of the app. `pkpdsim.webapp.body_figure(svg, sbml, *, simulation, title, credit="", substances=None, knobs=(), readouts=None)` returns the panels by their key, the figure and the keys which select nothing, ready to be declared in a `panels.py`.

```python
# not executed
from pkpdsim.webapp import body_figure

FIGURE_PANELS, FIGURE, FIGURE_DROPPED = body_figure(
    OUTPUT / "glimepiride_body_body.svg",
    OUTPUT / "glimepiride_body_flat.xml",
    simulation=SIMULATION,
    title="Glimepiride PBPK model",
    substances={"hydroxy-glimepiride (M1)": "M1"},
    knobs=KNOBS,
)
```

It does three things a consumer would otherwise do by hand, and the reason all three live here is that pkpdsim writes both the figure and the model, so it can read the ids back instead of asking a package to keep a list of what its model lacks.

- **The ids are restricted to the build.** An id of the sidecar is a proposal: the drawing shows every organ of the library while a model has the organs of its substances, and the readouts of the library are the ids a model has to request (`PBPKDefinition.readouts`). The entities of the flat build are read with libsbml, which costs a fraction of a second and no compile, and a key whose readouts are all absent is returned in `FIGURE_DROPPED` and stays a decoration of the drawing. Declaring a panel which selects nothing is refused by the registry, so this step is not optional.
- **Every panel carries the knobs of the app.** A simulation whose derived targets read a knob can only be evaluated on a panel which has that knob, so a generated panel without them is refused when the factory writes the annotation.
- **The curves are labelled as the panels of an app are.** A readout label is a legend of its own (`glimepiride, lung interstitium unbound`, 38 characters); a panel of one substance labels the compartment alone (`plasma`), a panel of several names the substance beside it in the short name `substances` gives (`M1 plasma`), and an unbound compartment reads `free <x>`.

The page of the app shows the figure with `AppPage(figure=True, shown=(...))`, which keeps the declared plots on the page and makes the generated panels the ones the drawing opens, and `App(figure=FIGURE)` carries the drawing itself.
