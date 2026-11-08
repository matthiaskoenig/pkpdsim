# Circadian rhythms

The physiology changes over the day and the night: the GFR is highest in the afternoon, a meal leaves the stomach more slowly in the evening, the clearance of midazolam (CYP3A4) is highest in the evening. A model with circadian rhythms (`PBPKDefinition(circadian=Circadian())`, `pkpdsim.pbpk.circadian`, #131) knows the clock time, so a dose in the morning and a dose in the evening (chronopharmacokinetics) are simulated with the same model.

- `CLOCK` is the clock time (hours) at the start of a simulation, time 0: an initial value, 8:00 by default (`Circadian(clock=8.0)`), which an experiment sets before the simulation (`PBPKSimulationExperiment.clock_changes(20.0)`, a `preinit_changes` of the simulation, or `init(CLOCK)` in roadrunner); `TOD` is the clock time during the simulation, 0 to 24 hours.
- Every rhythm is a cosinor relative to its mean over 24 hours, the factor `fcirc_<key> = 1 + A_circ_<key>*cos(2*pi*(CLOCK + time - phi_circ_<key>)/24 h)` with the relative amplitude `A_circ_<key>` and the acrophase `phi_circ_<key>` (the clock time of the maximum), constants of the model which a fit or an experiment may change. The factor is an assignment rule of the time and multiplies its quantity as the factors of the impairment do: `fimp_GFR` (the GFR), the transit time of the stomach in the fed state `T_st_fed` (an input of the intestine), `fimp_albumin` (the albumin of the binding, an optional rhythm) and `fimp_<o>_<protein>` (the activity of a protein which a process of the model uses).
- Every quantity which reads a factor follows it during the simulation: it is an assignment rule, as a quantity downstream of a parameter of `PBPKDefinition.dynamic`. A model without `circadian` has none of these elements; `Circadian(rhythms=("GFR", "albumin"))` chooses the rhythms by key, `Circadian()` takes the default ones.

<!-- generated: circadian-rhythms (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/circadian_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| rhythm | factor | amplitude | acrophase | 08:00 | 20:00 | default | sources |
| --- | --- | ---: | --- | ---: | ---: | --- | --- |
| GFR | `fcirc_GFR` | 0.15 | 16:30 | 0.9087 | 1.091 | yes | `Koopman1996` |
| t_stomach_fed | `fcirc_t_stomach_fed` | 0.21 | 20:00 | 0.79 | 1.21 | yes | `Goo1987` |
| albumin | `fcirc_albumin` | 0.05 | 16:30 | 0.9696 | 1.03 | no | `Koopman1996` |
| CYP3A4 liver | `fcirc_li_cyp3a4` | 0.071 | 18:50 | 0.9323 | 1.068 | yes | `VanRongen2015` |

??? note "Sources"

    - `Goo1987`: Goo RH, Moore JG, Greenberg E, Alazraki NP. Circadian variation in gastric emptying of meals in humans. Gastroenterology 93(3):515-518, 1987. [doi:10.1016/0016-5085(87)90913-9](https://doi.org/10.1016/0016-5085(87)90913-9)
    - `Koopman1996`: Koopman MG, Koomen GC, van Acker BA, Arisz L. Circadian rhythm in glomerular transport of macromolecules through large pores and shunt pathway. Kidney Int 49(5):1242-1249, 1996. [doi:10.1038/ki.1996.178](https://doi.org/10.1038/ki.1996.178)
    - `VanRongen2015`: van Rongen A, Kervezee L, Brill M, van Meir H, den Hartigh J, Guchelaar HJ, Meijer JH, Burggraaf J, van Oosterhout F. Population pharmacokinetic model characterizing 24-hour variation in the pharmacokinetics of oral and intravenous midazolam in healthy volunteers. CPT Pharmacometrics Syst Pharmacol 4(8):454-464, 2015. [doi:10.1002/psp4.12007](https://doi.org/10.1002/psp4.12007)
<!-- end generated -->

A dose in the morning and in the evening, with a meal: the evening dose leaves the stomach later (`fcirc_t_stomach_fed`), so its maximum is lower and later; the exposure over the day is the same here, the clearance of the drug of the example is mostly its metabolism, whose rhythm averages out over 24 hours.

```python
from pathlib import Path

import numpy as np
from sbmlsim.model.model_roadrunner import load_roadrunner

from pkpdsim.pbpk import (
    MM,
    Binding,
    Circadian,
    Excretion,
    Kp,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    build,
)

drug = Substance(
    sid="drg",
    name="drug",
    mr=300.0,
    binding=Binding(fu=0.1, protein="albumin"),
    distribution=Kp(2.0, ftissue=1.0),
    excretion=Excretion(urine=True),
)
definition = PBPKDefinition(
    sid="chrono",
    substances=[drug, Substance(sid="drg_m", name="metabolite", mr=316.0)],
    chemistry=[
        Metabolism(
            sid="DRG2M",
            organ=Organ.LIVER,
            substrate="drg",
            products={"drg_m": 1.0},
            kinetics=MM(vmax=1e-3, km=0.01),
            protein="CYP3A4",
        )
    ],
    routes=[Route.oral("drg")],
    circadian=Circadian(),
)
output = Path("models")
output.mkdir()
result = build(definition, output)

r = load_roadrunner(result.flat_path)
for clock in (8.0, 20.0):
    r.setValue("init(CLOCK)", clock)
    r.setValue("init(FED)", 1.0)
    r.resetAll()
    r.timeCourseSelections = ["time", "TOD", "Carm_drg", "fcirc_GFR"]
    r["PODOSE_drg"] = 10.0  # mg, with a meal
    s = r.simulate(0, 24 * 60, 1441)
    c = s["Carm_drg"] * 1000  # µM
    i = int(np.argmax(c))
    auc = float(np.trapezoid(c, s["time"])) / 60  # µM*h
    print(
        f"dose at {clock:04.1f} h: Cmax {c[i]:.3f} µM at {s['time'][i] / 60:.2f} h, "
        f"AUC(24 h) {auc:.2f} µM*h, GFR factor {s['fcirc_GFR'].min():.2f}-{s['fcirc_GFR'].max():.2f}"
    )
```

```text
dose at 08.0 h: Cmax 1.791 µM at 1.70 h, AUC(24 h) 8.32 µM*h, GFR factor 0.85-1.15
dose at 20.0 h: Cmax 1.532 µM at 1.90 h, AUC(24 h) 8.32 µM*h, GFR factor 0.85-1.15
```

<!-- generated: figures-circadian-1 (scripts/visual/docs_figures.py), do not edit -->
The content of the model `chrono`: the anatomical figure (A) and the schematic (B), hover them for details ([figures](figures.md)).

<iframe class="pkpd-figure" src="img/figures/models/circadian-1-anatomy.html" title="The anatomical figure of the model chrono" style="aspect-ratio:1460 / 898" loading="lazy"></iframe>

<iframe class="pkpd-figure" src="img/figures/models/circadian-1-schematic.html" title="The schematic of the model chrono" style="aspect-ratio:1384 / 902" loading="lazy"></iframe>
<!-- end generated -->

## No rhythm

The search (PubMed, Europe PMC; reviews Bicker 2020, Ballesta 2017, Baraldo 2008) found no rhythm, or no amplitude which could be used, for:

- the small intestinal transit of dosage forms (08:00 against 20:00, scintigraphy, Kambayashi 2020; the small bowel motility after a meal the same day and night), the fasted emptying of liquids (Goo 1987), the colonic transit (it spans more than a day);
- the renal plasma flow (day 843 and night 822 ml/min, Koopman 1996) and the renal clearance of drugs (gentamicin, Dickson 1986; levofloxacin, Kervezee 2016), CYP1A2 (caffeine at 10:00 and 22:00, Hashiguchi 1992), the other enzymes and transporters (no data of humans in vivo), AGP (amplitude unknown), the hematocrit (posture);
- the liver blood flow: a maximum at 08:00 (Lemmer and Nold 1991, ICG), amplitude not read; the cardiac output: about 29 % lower at night than at day (Veerman 1995), mostly activity and sleep and the flow of the muscle, not an endogenous rhythm of every flow; the gastric acid and pH (no gastric pH in the library).

## Limitations

- The amplitudes of the gastric emptying (Goo 1987) and the GFR (Koopman 1989) are from the abstracts; the amplitude of the gastric emptying is derived from a morning and an evening value. The GFR rhythm is robust under bed rest, but the renal clearance of drugs shows little of it.
- A periodic model has no steady state: the generated starts of an endogenous substance and the equilibration of a condition (`equilibrated`) are those at the clock time of the start, with the rhythms frozen at it.
- The rhythms are those of adults with a day-active schedule; shift work, the meals and the sleep are not modelled, nor the rhythms of the elderly and the children, nor an interaction with the age and the sex factors (they multiply).
- The PBPK tools have no circadian option of their own (Simcyp, PK-Sim); published chrono-PBPK models set the parameters over the day by hand (Marok 2024, Polak 2013).
