# The lobes of the liver

The liver of the [standard organ](organ.md) is one well-stirred chain, and the [zonated liver](variants.md#zonated-liver) is that chain cut into zones along the sinusoid. Neither can hold what a portal vein ligation or a hepatectomy does, because both act on a *region* of the liver and not on the organ: the lobes which lose their portal inflow atrophy while their neighbours grow, and the drug which the liver clears is cleared in the lobes which are still perfused. `LobedLiver` adds that second spatial axis.

The lobes are perfused in **parallel**, the zones of a lobe in **series**: a lobe is a chain of its own with its share of the liver volume, its own hepatic artery and portal vein inflow and its own drainage, and the organ drains the flow weighted mix of the lobes. The two axes compose, so a site of the liver carries both (`Cpl_lll_z3_<s>` is the plasma of the third zone of the left lateral lobe of the rat), and either axis of length one disappears: a liver of one lobe is the standard or the zonated liver line by line, exactly as a zonated liver of one zone is the standard liver.

## The lobe sets

A lobe set is the anatomy of a species. The rodents get their anatomical lobes and the human the eight Couinaud segments, which are the units a resection is planned in; the shares of the mass come from `pkpdsim.database` and a lobe without a measured share is marked `ModelChoice` (see [Gaps](#gaps)).

```python
from pkpdsim.pbpk.lobes import lobe_table
from pkpdsim.species import Species

for species in (Species.HUMAN, Species.RAT, Species.MOUSE):
    lobes = lobe_table(species)
    shares = ", ".join(f"{lobe.token} {lobe.mass_fraction.value:.3f}" for lobe in lobes)
    print(f"{species.value}: {shares}")
```

```text
human: sg1 0.040, sg2 0.069, sg3 0.095, sg4 0.136, sg5 0.130, sg6 0.079, sg7 0.190, sg8 0.261
rat: cl 0.100, lll 0.300, lml 0.200, rl 0.200, rml 0.200
mouse: cl 0.081, ll 0.344, ml 0.262, ri 0.147, rs 0.166
```

Every set sums to 1 and `check_lobes` refuses one which does not, together with a token which is no lower case identifier, a repeated token, a negative share and a set of supply fractions which does not sum to 1. A species without a set says so:

```python
from pkpdsim.pbpk.lobes import lobe_table
from pkpdsim.species import Species

try:
    lobe_table(Species.PIG)
except ValueError as error:
    print(error)
```

```text
Species 'pig' has no liver lobe set; the species with a lobe set are human, mouse, rat.
```

A model may give its own lobes instead, which is how a study with a different partition or a lobe of a known supply fraction is built: `LobedLiver(lobes=[Lobe("left", "left lobe", Value(0.6, U.dimensionless, "Study2020")), ...])`.

## The structure

```python
from pkpdsim.pbpk import Organ, Variants, ids
from pkpdsim.pbpk.organs import Part
from pkpdsim.pbpk.variants.liver_lobes import LobedLiver

variants = Variants(organs={Organ.LIVER: LobedLiver()})
liver = variants.organs[Organ.LIVER]
print(liver.name, liver.lobes, liver.zonation)
print(ids.site("lll", None), ids.site("lll", ids.zone(3)), ids.site(None, ids.zone(3)))
print(ids.compartment(Part.PLASMA, ids.site("lll", ids.zone(3))))
print(ids.species(Part.PLASMA, "apap", ids.site("lll", ids.zone(3))))
```

```text
lobes None None
lll lll_z3 z3
Vpl_lll_z3
Cpl_lll_z3_apap
```

Every lobe holds its share `f_vol_<lobe>` of the volumes and the passive exchanges of the liver, and of its capacity: a metabolism or a transport of the liver is placed in every site with the share of its lobe times the share of its zone, so the capacity of the liver is the capacity of the standard liver whatever the two axes are. The blood of a lobe flows from the hepatic artery and the portal vein into its periportal end, through its zones, and out of its pericentral end into the venous pool; no reaction crosses from one lobe into another. Each lobe drains at its own concentration `Cout_<lobe>_<s>` and the organ at the flow weighted mix

\[
C_{out,s} = \frac{\sum_l Q_{h,l}\,C_{out,l,s}}{\sum_l Q_{h,l}},
\]

under the unchanged id `Cout_<s>`, so the hepatic and the arm vein observers, the extraction, the steady state and the pressure driven circulation read what they read before.

That mix is also the equivalence of the variant: when every lobe is supplied in proportion to its mass, every lobe has the same perfusion per volume and extracts the same fraction `CL/(Q + CL)` of what reaches it, whatever its share, so the mix of the lobes is the outflow of the lumped liver. A lobed liver is then the liver without lobes, which `tests/pbpk/test_lobes_build.py` asserts on the venous plasma curve of a whole day.

### The ids of a lobe

| id | what it is |
|---|---|
| `Vli_<lobe>` | the volume of the lobe (l); the volumes of the lobes which are in place sum to `Vli_lobes` |
| `f_vol_ref_<lobe>` | the share of the intact liver which the lobe holds, the constant of the lobe table |
| `f_vol_<lobe>` | `Vli_<lobe>/Vli`, the share the lobe holds now; it follows the volume when the anatomy changes |
| `f_present_<lobe>` | 1 is the lobe in place, 0 the resected one, which has no capacity and no growth |
| `g_ha_<lobe>`, `g_po_<lobe>` | the patency of the hepatic artery and the portal vein branch, 1 open and 0 ligated |
| `f_ha_<lobe>`, `f_po_<lobe>` | the reference shares of the two inflows, the share of the mass where no source gives one |
| `Qha_<lobe>`, `Qpo_<lobe>`, `Qh_<lobe>` | the blood flows of the lobe (l/min), which sum to `Qha`, `Qpo` and `Qh` while a branch of that kind is open |
| `q_<lobe>` | `Qh_<lobe>/Vli_<lobe>`, the perfusion per volume: the stimulus of the remodelling |
| `dVli_<lobe>` | the growth of the lobe (l/min), which a driver of the remodelling sets |
| `Vli_ref`, `Vli_lobes`, `f_remnant` | the intact liver, the liver which is there and their ratio |
| `Cout_li_<lobe>_<s>` | the plasma concentration leaving the lobe; its area is the exposure of the lobe |

A lobe is supplied in proportion to its reference share among the lobes whose branch is open,

\[
Q_{po,l} = Q_{po}\,\frac{g_{po,l}\,f_{po,l}}{\sum_j g_{po,j}\,f_{po,j}},
\]

so a lobe whose portal branch is ligated passes its flow to the lobes which are open, which is what happens when the splanchnic bed keeps delivering its blood to a liver with one branch less. The denominator is guarded, so a liver whose every branch is closed has the flow 0 in every lobe and no undefined rate. An occlusion always acts on the conductance of a branch and never on its resistance, for the same reason.

## The remodelling

After an intervention the anatomy of the liver follows the simulation over days, so the volume of a lobe is a state of the model and not a constant of the individual. `LobedLiver(remodelling=...)` says what moves it, and every driver moves it through the one quantity `dVli_<lobe>`, so the liver always knows how fast a lobe is losing tissue.

`FlowDriven` is the predictive driver. Its stimulus is the perfusion per volume of the lobe against the reference perfusion of the liver,

\[
\frac{dV_{li,l}}{dt} = f_{present,l}\;k_{growth}\;V_{li,l}\left(\frac{q_l}{q_{target}} - 1\right),
\]

which is what the hyperperfusion of the remnant after a resection and the hypoperfusion of a deportalised lobe after a ligation are: the remnant carries the whole hepatic inflow in a fraction of the volume and grows, a deportalised lobe keeps its arterial inflow only and atrophies, and the liver comes to rest when every lobe is at the reference perfusion again, which is the intact liver. The law therefore needs no ceiling of its own; `total_cap=True` adds one for a study which wants a hard limit, at the price of a much slower approach.

`PrescribedVolume` and `PrescribedFlow` map a measured course onto the model instead, the volume of a lobe over time as an imaging study reports it or its growth rate. A prescribed volume is tracked with the time constant `tau_track` rather than assigned, because the growth rate has to be a quantity of the model and the rate of change of an assignment rule is not available to a rate law; one hour is fast against the days of remodelling.

```python
from pkpdsim.pbpk.remodelling import FlowDriven, PrescribedVolume
from pkpdsim.sbml.values import Formula

flow_driven = FlowDriven()
print(f"k_growth {flow_driven.k.value} 1/min, source {flow_driven.k.source}")
print(f"total_cap {flow_driven.total_cap}")
mapped = PrescribedVolume({"lll": Formula("Vli_ref * 0.2 dimensionless", "l")})
print(f"tau_track {mapped.tau.value} min, source {mapped.tau.source}")
```

```text
k_growth 0.0005 1/min, source ModelChoice
total_cap False
tau_track 60.0 min, source ModelChoice
```

### What atrophy carries away

The state of a species of the liver is its amount, so a lobe which grows spreads the same amount over the larger volume: growth dilutes the tissue by itself, which is what growing tissue does. A lobe which shrinks would concentrate it instead, which no tissue does. The remodelling therefore emits one flux per species and part of a shrinking lobe,

\[
J_{remodel} = V_{part}\,C_{part}\,\frac{\max(0,\,-dV_{li,l})}{V_{li,l}},
\]

into the sink `Aremodelled_<s>` of the body, where the amount which left is auditable. The flux is off for a growing lobe, which `tests/pbpk/test_lobes_remodelling.py` asserts.

## The interventions

The models of the library have no events, so an intervention is a set of parameter changes, as the [impairment](impairment.md) maps are: a study applies them before a simulation or between the phases of one, which is what a surgery at a time point is. The patency of the branches and the presence of every lobe are variable for that reason, so the lobe flows follow a change at the point it is made.

| intervention | what it changes | what it is |
|---|---|---|
| `PortalLigation(lobes)` | `g_po_<lobe>` 0 | portal vein ligation or embolisation: the lobes keep their artery, stay perfused below the reference and atrophy while the others hypertrophy |
| `ArterialLigation(lobes)` | `g_ha_<lobe>` 0 | hepatic artery ligation; both together devascularise a lobe |
| `Resection(lobes)` | both branches and `f_present_<lobe>` 0 | a hepatectomy: the lobes are gone, hold no capacity and leave the remnant fraction |

The named procedures of the literature come with the species:

```python
from pkpdsim.pbpk.interventions import scenarios
from pkpdsim.species import Species

for species in (Species.HUMAN, Species.RAT, Species.MOUSE):
    for name, intervention in scenarios(species).items():
        kind = type(intervention).__name__
        print(f"{species.value:5s} {name:26s} {kind:15s} {list(intervention.lobes)}")
```

```text
human right_pvl                  PortalLigation  ['sg5', 'sg6', 'sg7', 'sg8']
human right_hepatectomy          Resection       ['sg5', 'sg6', 'sg7', 'sg8']
human left_hepatectomy           Resection       ['sg2', 'sg3', 'sg4']
human extended_right_hepatectomy Resection       ['sg4', 'sg5', 'sg6', 'sg7', 'sg8']
rat   pvl70                      PortalLigation  ['lll', 'lml', 'rml']
rat   pvl90                      PortalLigation  ['lll', 'lml', 'rml', 'rl']
rat   hepatectomy70              Resection       ['lll', 'lml', 'rml']
rat   hepatectomy90              Resection       ['lll', 'lml', 'rml', 'rl']
mouse pvl                        PortalLigation  ['ll', 'ml']
mouse hepatectomy                Resection       ['ll', 'ml']
```

`scenario_changes(species)` gives the same procedures as the parameter maps which `PBPKSimulationExperiment.lobe_map` takes, and asking a species for the procedure of another one says so:

```python
from pkpdsim.pbpk.interventions import scenarios
from pkpdsim.species import Species

rat = ("cl", "lll", "lml", "rl", "rml")
print(sorted(scenarios(Species.RAT)["hepatectomy70"].changes(rat)))
try:
    scenarios(Species.RAT)["right_hepatectomy"]
except ValueError as error:
    print(error)
```

```text
['f_present_lll', 'f_present_lml', 'f_present_rml', 'g_ha_lll', 'g_ha_lml', 'g_ha_rml', 'g_po_lll', 'g_po_lml', 'g_po_rml']
'right_hepatectomy' is no procedure of the liver of the rat; its procedures are ['hepatectomy70', 'hepatectomy90', 'pvl70', 'pvl90'].
```

## The lobes on the pressure network

With the [pressure driven circulation](cardiovascular.md) the flows of the liver are no longer shares of a prescribed cardiac output but the solution of a resistance network, and a lobed liver gets a node of its own per lobe: the sinusoidal pressure `CV_P_li_<lobe>`, with the hepatic artery and the portal vein flowing into it and the hepatic vein out of it. The balance of the three branches solves the pressure of the lobe,

\[
Q_{ha,l} + Q_{po,l} - Q_{h,l} = 0,
\qquad
Q_{ha,l} = \frac{(P_{ar} - P_{li,l})\,g_{ha,l}\,f_{habr,l}}{R_{ha,l}},
\quad
Q_{po,l} = \frac{(P_{po} - P_{li,l})\,g_{po,l}}{R_{po,l}},
\quad
Q_{h,l} = \frac{P_{li,l} - P_{ve}}{R_{h,l}},
\]

and each resistance `CV_R_ha_<lobe>`, `CV_R_po_<lobe>` and `CV_R_h_<lobe>` is calibrated from the reference pressure drop over its branch and the reference flow through it (`Qha_ref_<lobe>`, `Qpo_ref_<lobe>`, `Qh_ref_<lobe>`, the reference share of the lobe of the flows of the individual), as the resistance of an organ of the circulation is, so at baseline every lobe carries its share and the network agrees with the supply fraction mode above. The whole liver multipliers `CV_fR_hepatic`, `CV_fR_portal` and `CV_fR_li` apply to every lobe, which is what lets a cirrhosis or a vasodilator act on the liver as a whole while a ligation acts on one lobe. `Qha`, `Qpo` and `Qh` are the sums of the lobes, and `CV_portal_balance` stays on the sum, so the splanchnic bed delivers its blood into a liver with one branch fewer and its portal pressure rises, which is what a ligation does and what a prescribed flow cannot show.

The hepatic outflow of a lobe is deliberately written without the `max(..., 0 mmHg)` of the whole liver branch: the balance is then linear in the pressure of the lobe and has the one solution `P_li_lobe = P_ve + (Qha_lobe + Qpo_lobe) * R_h_lobe`, so a lobe whose both branches are closed sits at the venous pressure and carries no flow instead of leaving the solver an under determined node. The inflows are never negative, so the outflow never is either.

One occlusion has no solution on the network: closing the portal branch of *every* lobe leaves the splanchnic bed nowhere to drain, because the portal junction balance sends the flow of the portal organs into the liver, so the solver fails at the initial conditions. The supply fraction mode takes the same occlusion without failing, and what it then says is not a model of anything: the portal pool has one outlet in the body, so with every branch closed the splanchnic organs fill a dead end, the whole dose ends up in the portal vein and the systemic exposure goes to 0 (`test_a_total_portal_occlusion_strands_the_drug_in_the_portal_vein`). A total portal occlusion is a liver perfused by collaterals, so it needs the portosystemic shunt of the gaps below, and neither mode has it. Occlude the branches of some of the lobes, which is what a portal vein ligation is.

### The hepatic arterial buffer response

A liver which loses its portal inflow raises its arterial inflow, because the adenosine which the portal flow washes out of the space of Mall accumulates and dilates the arteriole (Lautt 1985). The response is of the territory and not of the organ, which is why it belongs to the lobes: `HABR(gain, maximum)` multiplies the arterial conductance of every lobe by

\[
f_{habr,l} = \min\left(f_{max},\;
1 + \mathrm{gain}\,\frac{\max(0,\,Q_{po,l}^{ref} - Q_{po,l})}{Q_{ha,l}^{ref}}\right),
\]

so a lobe which loses a part of its portal flow gains `gain` times that flow in its artery, up to the ceiling. The default gain is the measured 25.5 % of Lautt 1985 (`liver_remodelling.toml`, `habr_gain`), and the response is off unless a model asks for it. It acts on the arterial conductance of a lobe, which only the pressure driven circulation has: a lobed liver with a `HABR` and no `PressureFlow` is refused at the definition instead of building a model which ignores it.

```python
from pkpdsim.pbpk.lobes import HABR

habr = HABR()
print(f"gain {habr.gain.value} ({habr.gain.source}), maximum {habr.maximum}")
```

```text
gain 0.255 (Lautt1985), maximum 2.0
```

## The app of a model with lobes

The app which a model generates ([web app of a model](../webapp.md)) gets a page of the lobes when the model has them, one choice per intervention per lobe: the portal vein and the hepatic artery branch of the lobe open or ligated, and the lobe in place or resected. Intact is the default of every one, so the app opens on the healthy liver and a reader ligates or resects from there and sees the readouts of the model change. The [interactive explorer](lobes_examples.md) goes further: a procedure is planned there by picking what is done and the lobes it acts on, the named procedures of the database fill the plan, and the days of the remodelling are a course which is played or scrubbed.

## Gaps

What the lobe library takes as a model choice, every one of them marked `ModelChoice` in its table and listed here so a study knows what it is standing on:

- **The supply fractions.** No publication of the parameter database gives the share of the hepatic artery or the portal vein inflow of a single lobe in any species, so a lobe is supplied in proportion to its mass. A Couinaud segment is defined as a portal territory and the arterial branch runs with the portal branch in the portal triad, so the perfusion per volume is the same in every lobe, which is the assumption the equivalence with the lumped liver rests on.
- **The human segment volumes.** Only the share of segment VIII is sourced (Mise 2014); the others keep the published descending order VIII > VII > IV > V > III > VI > II > I and the bound of 10 % on the segments II and III, and are model choices within it.
- **The split of the lobes of the rat.** Madrahimov 2006 fixes the lobes of the rat around its hepatectomies and gives the sum of the left lateral, the median and the right lobe (the 90 % hepatectomy), not the split of that sum into the three lobes (0.30, 0.40, 0.20, the conventional split) nor the split of the median lobe into its left and its right half (taken as even): all four are `ModelChoice` rows of the table.
- **The growth constant.** No publication states a rate constant of the flow driven growth law, so the default is the one which restores the liver in about a week after a 70 % resection, the time course the hepatectomy literature reports. A measured course is mapped with `PrescribedVolume` instead of fitted through the law.
- **The tracking time of a prescribed volume**, one hour, is a numerical choice and not a physiological one.
- **The split of the hepatic resistances over the lobes.** The reference resistance of a branch of a lobe is the resistance of the whole liver divided by the reference share of the lobe, which is the one split which reproduces the baseline flows; no publication of the database gives the resistance of a single branch.
- **The absence of a portosystemic shunt of a lobe.** The shunt `f_shunt` of the liver is of the whole organ and the pressure network has no collateral vessel, so a liver whose every portal branch is closed has no solution on the network (above) and a portal pressure which rises is not relieved by collaterals, as it is in a cirrhotic liver. The functional share `f_functional` of a cirrhotic liver is of the whole organ too, so a cirrhosis is uniform over the lobes.
- **The conductance of a lobe does not follow its volume.** The reference resistance of a branch is a constant of the build, so the split of the hepatic inflow over the lobes which stay open keeps the reference shares while the volumes move: a lobe which has doubled does not draw a larger share of the portal inflow than it did on the day of the operation. What the remodelling changes is the perfusion per volume, which is the stimulus of the growth law and the exposure of the lobe, and what the interventions change is the split; a conductance which grows with the volume would couple the two and is not implemented.
- **The ceiling of the buffer response**, a doubling of the arterial inflow of a lobe, is a model choice: Lautt 1985 measured the gain and not a saturation.
- **The procedures of the mouse and the human** are the anatomical standard (the Brisbane 2000 terminology for the human), which no publication of the database states. The two thirds hepatectomy of the mouse is named `hepatectomy` and not `hepatectomy70`, because the lobes it removes are 0.61 of the mass with the shares of Brennan 2024 and not the two thirds the procedure is named after in the literature.
