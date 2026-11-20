# Interactive liver lobe examples

Explore portal vein ligation and hepatectomy on the lobes of the rat, the mouse and the human. This explorer is part of the [liver lobe documentation](lobes.md#the-app-of-a-model-with-lobes): a procedure is planned by picking what is done and which lobes it acts on, and the five days which follow are the course of the flows, the volumes and the exposure of every lobe. The plan comes first, the dashboard of the run below it.

<iframe class="pkpd-figure" src="examples/lobes.html" title="Interactive liver lobe examples: portal vein ligation and hepatectomy" style="height:1200px" loading="lazy"></iframe>

[Open the explorer in its own window](examples/lobes.html){ target="_blank" rel="noopener" } for a larger liver and larger plots. The page runs locally in the browser on the files of this documentation: the precomputed reference courses of `lobes-data.json` (`female/lobes-data.json` for the female human, loaded when a reader asks for it) and, for anything the dataset does not hold, the compiled lobed livers of `solver_lobes/<species>/` and `female/solver_lobes/human/` on the sbmlsim-core WebAssembly runtime of the [cardiovascular explorer](cardiovascular_examples.md). It holds no model of a liver of its own, so what it shows is what the SBML says.

## Plan the procedure

The planner is the procedure and the lobes it acts on: pick the portal vein ligation, the hepatic artery ligation or the resection, then click the lobes it is done to. The named procedures of the parameter database fill both at once, and the line below the two says what the plan costs, the share of the liver mass it acts on and the remnant it leaves. A click on a block of the liver cycles the same states, so a reader who wants one lobe in another state does not have to leave the drawing.

## Read the liver

The area of the whole liver is the mass which is there, measured against the intact liver of the same individual, so a hepatectomy shrinks the drawing and the regrowth fills it again; the dashed outline is the liver before the operation. Within it every lobe is a block whose area is the share of the mass it holds at the displayed time, so the atrophy of a deportalised lobe and the hypertrophy of its neighbours are a change of area. The lobes stand in columns which come from the lobe table and never from the shares of the moment: a lobe keeps its place and its neighbours and only grows or shrinks. A resected lobe has no mass in the liver any more and therefore no block: it is listed under the liver as a box of its own colour, shaded as resected, which puts it back. The eight blocks of the human are the Couinaud segments, the five of the rat and of the mouse its anatomical lobes; the labels are the tokens of the lobe tables, which the ids of the model are built from.

Beside the liver stands the animal of the run, the anatomical figure of the library ([physiology figures](figures.md)) with its liver picked out, so a reader sees which individual the lobes belong to; the drawings are each scaled to their own box, so the species are not on one scale.

The switch beside the species is the sex of that individual. The parameter database has a female reference physiology of the human and the male alone of the rat and of the mouse, so the page offers the sex where the data carries it and says on the disabled button why it does not: there is no female rat to show. The female human is a second individual and not a second drawing, a liver of 1.79 l against the 2.33 l of the male and a hepatic inflow of 1.55 l/min against 1.79 l/min, with its own model, its own precomputed runs and its own compiled liver in the browser.

A lobe has one colour, in the drawing, in the curves and in the bars alike, and it keeps it when another lobe is resected and leaves the diagram. What was done to a lobe is a shading of that colour and not a colour of its own: a portal vein ligation hatches the block one way, a hepatic artery ligation the other.

```python
from pkpdsim.pbpk.lobes import lobe_table
from pkpdsim.species import Species

for lobe in lobe_table(Species.HUMAN):
    print(f"{lobe.token}: {lobe.name}, {lobe.mass_fraction.value:.3f} of the mass")
```

```
sg1: segment I (caudate lobe), 0.040 of the mass
sg2: segment II (left lateral superior), 0.069 of the mass
sg3: segment III (left lateral inferior), 0.095 of the mass
sg4: segment IV (left medial), 0.136 of the mass
sg5: segment V (right anterior inferior), 0.130 of the mass
sg6: segment VI (right posterior inferior), 0.079 of the mass
sg7: segment VII (right posterior superior), 0.190 of the mass
sg8: segment VIII (right anterior superior), 0.261 of the mass
```

Three vessels belong to every lobe: the hepatic artery and the portal vein branch which run to it from the hilum below, and the hepatic vein which carries its outflow to the vessel above. An edge becomes wider as its flow increases, with the width growing with the square root of the flow so that the branches of a small lobe stay visible, and a closed branch keeps a thin dashed path to show the connection while its flow is zero. A branch the plan ligated is cut across, two strokes on the vessel itself, where a branch which merely carries nothing is only dashed.

A click on a block, or Enter on a focused one, cycles the four states a lobe can be in. They are the interventions of the library and nothing else: the state of the page is a set of parameter changes applied to the same model.

| state | what the page changes | what it means |
| --- | --- | --- |
| intact | nothing | both branches open, the lobe perfused as the lobe table says |
| portal vein ligated | `g_po_<lobe>` 0 | the portal branch closed, the lobe on its artery alone, hypoperfused, atrophying |
| hepatic artery ligated | `g_ha_<lobe>` 0 | the arterial branch closed, the portal inflow of the lobe untouched |
| resected | `g_ha_<lobe>` 0, `g_po_<lobe>` 0, `f_present_<lobe>` 0 | the lobe gone: no inflow, no metabolic capacity, no growth |

## The named procedures

The presets are the procedures of the parameter database, one set per species, and the explorer adds a portal vein ligation and a resection of the largest single lobe for a reader who wants the simplest case.

```python
from pkpdsim.pbpk.interventions import scenarios
from pkpdsim.species import Species

for species in (Species.RAT, Species.MOUSE, Species.HUMAN):
    print(f"{species.name.lower()}: {', '.join(scenarios(species))}")
```

```
rat: pvl70, pvl90, hepatectomy70, hepatectomy90
mouse: pvl, hepatectomy
human: right_pvl, right_hepatectomy, left_hepatectomy, extended_right_hepatectomy
```

A preset is the changes of its scenario, which is what the page applies when you choose it; a right portal vein ligation of the human closes the portal branches of the four right segments and touches nothing else.

```python
from pkpdsim.pbpk.interventions import scenario_changes
from pkpdsim.species import Species

for pid, value in scenario_changes(Species.HUMAN)["right_pvl"].items():
    print(f"{pid} = {value}")
```

```
g_po_sg5 = 0.0
g_po_sg6 = 0.0
g_po_sg7 = 0.0
g_po_sg8 = 0.0
```

## The controls

The explorer runs the model the sliders and the lobe states define. The model behind it is a lobed liver with the flow driven remodelling and the hepatic arterial buffer response on the pressure network, the definition of `scripts/lobes_examples.py`:

```python
from pkpdsim.pbpk import HABR, FlowDriven, LobedLiver, Organ, PressureFlow, Variants

variants = Variants(
    circulation=PressureFlow(),
    organs={Organ.LIVER: LobedLiver(remodelling=FlowDriven(), habr=HABR())},
)
liver = variants.organs[Organ.LIVER]
print(f"growth {liver.remodelling.k.value} 1/min, buffer gain {liver.habr.gain.value}")
```

```
growth 0.0005 1/min, buffer gain 0.255
```

**The buffer response gain** is `habr_gain`, the fraction of a lost portal inflow the hepatic artery of that lobe compensates. Its default is the 0.255 of Lautt 1985; at 0 the arterial branches do not react and a deportalised lobe is as ischaemic as its portal loss makes it, which is what a model without a buffer response does.

**The growth constant** is `k_growth`, the relative growth of a lobe at twice the reference perfusion per volume. Its default restores the liver in about a week; at 0 the anatomy is frozen and the explorer shows the redistribution of the flows alone.

**The time of the surgery** decides how the intervention enters the simulation. At day 0 the changes are applied before the initialization, which is how a study of the days after an operation sets it up. At a later day they are applied to the running model, which is the surgery at a time point of a timecourse: the model has no events, so a change at a time is a change between two phases.

**The hepatic arterial inflow** and **the portal inflow** scale the gates of every branch which the plan left open, `g_ha_<lobe>` and `g_po_<lobe>`, the same parameters a ligation sets to zero. At 1 the liver is perfused as the lobe table says; below it the inflow of the whole liver is throttled, which a reader uses to ask what a lobe sees when the supply of the organ falls rather than the supply of one lobe.

A preset at the default sliders is a precomputed reference run of the dataset and appears at once. Moving a slider, or asking for a run of your own, puts the request to the compiled model in a Web Worker: a spinner runs while it is solved and the time it took stands beside the button afterwards, a few tens of milliseconds for a rat. The page states for every result whether it is the precomputed run or the live one.

The days after the intervention are played with the button beside the time, which walks the course from the operation to the fifth day and stops there; the slider beside it scrubs to any day by hand.

## What the plots show

The dashboard reads one run in two ways and the switch above it changes nothing else: **volumes** is what the remodelling does, **pressures and flows** is what drives it. No simulation is run by the switch, because both readings are columns of the same result.

**The volume of the lobes** is the five days of `Vli_<lobe>` of the lobes which are there, the course the remodelling drives. A resected lobe is left out of the plots as it is left out of the liver: the model keeps its volume where it was on the day of the operation, which is no statement about a lobe which is gone. A 70 % hepatectomy of the rat starts at a remnant of 0.300 of the intact volume and reaches 0.981 after five days with the default growth constant; the caudate lobe of that remnant grows from 0.87 ml to 2.85 ml. A right portal vein ligation of the human does not change the volume of the liver at all, it redistributes it: segment VIII falls from 0.608 l to 0.269 l while segment III rises from 0.221 l to 0.461 l.

**The blood flow per lobe** is `Qha_<lobe>` and `Qpo_<lobe>` at the displayed time, which is where the interventions act. The flows redistribute at the moment of the surgery and then stay where they are: the reference resistance of a branch is a constant of the build, so the lobes which stay open keep the reference shares of the hepatic inflow while their volumes move ([Gaps](lobes.md#gaps)). What the days change is the perfusion per volume, the stimulus of the growth law.

**The liver which is there** is `f_remnant`, the fraction of the intact volume which is present and perfused: it starts at the remnant of the operation and rises as the lobes which are left grow.

**The exposure of the lobes** is the area under the concentration leaving each lobe, `Cout_li_<lobe>_drg` over the five days, which is what a study of a regional intervention reads: a deportalised lobe sees a different exposure than its neighbours because it is perfused differently, not because its enzymes differ.

In the pressure view the same run answers the other question. **The pressure of the hepatic inflow** is `CV_P_po` and `CV_P_li`, the portal pressure and the pressure of the sinusoids: an intact rat liver sits at 10.0 mmHg in the portal vein, a 70 % portal vein ligation at 19.1 and a 90 % one at 36.8, which is the portal hypertension the remaining bed is driven at. **The pressure per lobe** is `CV_P_li_<lobe>`, where the intervention is read off directly: after the 90 % ligation of the rat the caudate lobe stands at 23.1 mmHg and the four deportalised lobes at 5.8, the pressure of their artery alone. **The hepatic inflow of the liver** is `Qh`, and **the perfusion per volume** is `q_<lobe>`, the flow per volume of a lobe, which is the stimulus the growth law reads.

The pressures do not move over the five days: the reference resistance of a branch is a constant of the build, so what the remodelling changes is the volume and not the conductance ([Gaps](lobes.md#gaps)). The plots are flat in time and differ between plans, which is where the pressure network says something.

The download gives the displayed course as comma separated text, the time, the remnant and the volume, the two inflows and the outflow concentration of every lobe.

## What this is not

The explorer simulates the model of the [liver lobe documentation](lobes.md) and inherits every one of its gaps. The lobes are in parallel on one hepatic artery, one portal vein and one hepatic vein, with no collateral vessel and no portosystemic shunt of a single lobe; a lobe is supplied in proportion to its mass, because no publication of the parameter database gives the supply of a single lobe; and the growth law is a model choice to be fitted, not a measured constant. The days of a course are the days of the model and not the days of a patient: a reader comparing a remnant of this page with the volumetry of a study is comparing with a time constant which was chosen to restore the liver of the rat in about a week.

The test drug of the dataset is a toy substance with a first order hepatic metabolism, chosen so that the exposure of a lobe is a readable number; it is not a model of a compound.
