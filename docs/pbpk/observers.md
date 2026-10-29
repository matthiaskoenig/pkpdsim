# Observers

The observers are assignment rules, no states: they read the pool species and the outflow concentrations of the organs and change nothing, so they conserve mass exactly. A definition adds observers of its own (below). The forearm compartment of some hand-written models lost mass and was removed from three of them; the arm vein of the library is an observer.

| observable | id | what |
|---|---|---|
| plasma of a pool | `Car_<s>`, `Cve_<s>`, `Cpo_<s>` | the species of the arterial and venous pool and the portal vein |
| unbound plasma of a pool | `C<pool>_free_<s>` | `fu*C` |
| blood of a pool | `C<pool>_blood_<s>` | `(1 - HCT)*C + HCT*C_bc`: `BP*C` unless the blood cells are kinetic |
| hepatic vein | `Chv_<s>`, `Chv_free_<s>`, `Chv_blood_<s>` | the outflow of the liver mixed with the portosystemic shunt |
| arm vein | `Carm_<s>`, `Carm_free_<s>`, `Carm_blood_<s>` | the peripheral venous sampling site, the default fit target |
| tissue | `<O>__Cti_<s>`, `<O>__Cor_<s>`, `<O>__Cic_<s>`, ... | the species of the organs (concentrations) |
| unbound tissue | `<O>__Cu_ti_<s>`, `<O>__Cu_ic_<s>`, ... | the unbound concentrations of the tissue parts, see [Standard organ](organ.md#kp-and-fu) |
| outflow of an organ | `Cout_<o>_<s>` | the plasma leaving an organ |
| element of an organ read by a PD component or an observer | `C<o>_<part>_<s>`, `C<o>_u_<part>_<s>`, `A<o>_<x>` | a species or an unbound concentration of an organ as a parameter of the top model, see below |
| urine, feces | `Aurine_<s>`, `Afeces_<s>` | the amounts in the sinks |
| urinary excretion rate | `Rurine_<s>` | the net flux of the kidney into the urine (mmole/min), an output of the kidney, see [Standard organ](organ.md#kidney) |
| lumen of the intestine | `Alumen_<s>` | the dissolved amount in the lumen, the sum over the segments (mmole), a readout, see below |
| absorption, efflux, feces | `Rabs_<s>`, `Reff_<s>`, `Rfeces_<s>` | the summed fluxes of the intestine (mmole/min), readouts, see below |
| bile | `Cbile_<s>` | the concentration in the bile |

The equations of the top model, for a substance without erythrocytes and for one with kinetic erythrocytes (`prm`), with the preset `pksim_default`:

```equations title="pools and hepatic vein"
Cve_free_<s> = fu_<s> * Cve_<s>
Cve_blood_<s> = BP_<s> * Cve_<s>
Cve_blood_prm = (1 - HCT) * Cve_prm + HCT * Cve_bc_prm
Chv_<s> = ((1 - f_shunt) * Qh * Cout_li_<s> + f_shunt * (Qha * Car_<s> + Qpo * Cpo_<s>)) / Qh
Chv_free_<s> = fu_<s> * Chv_<s>
```

```equations title="arm vein (preset pksim_default)"
Carm_<s> = w_arm_sk * Cout_sk_<s> + w_arm_mu * Cout_mu_<s>
Carm_free_<s> = fu_<s> * Carm_<s>
Carm_blood_<s> = BP_<s> * Carm_<s>
Carm_blood_prm = (1 - HCT) * Carm_prm + HCT * (w_arm_sk * Cout_sk_bc_prm + w_arm_mu * Cout_mu_bc_prm)
```

## Arm vein

Plasma and blood are typically sampled from an antecubital or forearm vein, whose blood is the venous outflow of the tissues of the arm, not the mixed venous blood of the venous pool. The arm vein observer is the flow-weighted mix of the venous plasma outflows of organs (and the arterial plasma) with the weights `w_arm_<site>` of a preset of `data/arm_vein.toml`, selected by `PBPKDefinition.arm_vein`. The default `pksim_default` is the mix of skin and muscle of the PK-Sim documentation ("about 70 %" skin and "30 %" muscle), taken as 0.7 and 0.3 as documented and not from the database of PK-Sim; the presets of Levitt 2004 and Musther 2015 add the arteriovenous shunt of the hand, and `arterialised` (a heated hand vein) is the arterial plasma.

<!-- generated: arm-vein (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/arm_vein.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

| preset | adipose | arterial | muscle | skin | source |
| --- | ---: | ---: | ---: | ---: | --- |
| `pksim_default` |  |  | 0.3 | 0.7 | `OSPDocumentation` |
| `levitt2004` | 0.075 |  | 0.05 | 0.875 | `Levitt2004` |
| `musther2015` | 0.075 | 0.1 | 0.05 | 0.775 | `Musther2015` |
| `arterialised` |  | 1 |  |  | `Liu1992` |

??? note "Sources"

    - `Levitt2004`: Levitt DG. Physiologically based pharmacokinetic modeling of arterial - antecubital vein concentration difference. BMC Clin Pharmacol 4:2, 2004. [doi:10.1186/1472-6904-4-2](https://doi.org/10.1186/1472-6904-4-2)
    - `Liu1992`: Liu D, Moberg E, Kollind M, Lins PE, Adamson U, Macdonald IA. Arterial, arterialized venous, venous and capillary blood glucose measurements in normal man during hyperinsulinaemic euglycaemia and hypoglycaemia. Diabetologia 35(3):287-290, 1992. [doi:10.1007/BF00400932](https://doi.org/10.1007/BF00400932)
    - `Musther2015`: Musther H, Gill KL, Chetty M, Rostami-Hodjegan A, Rowland M, Jamei M. Are physiologically based pharmacokinetic models reporting the right Cmax? Central venous versus peripheral sampling site. AAPS J 17(5):1268-1279, 2015. [doi:10.1208/s12248-015-9796-7](https://doi.org/10.1208/s12248-015-9796-7)
    - `OSPDocumentation`: Open Systems Pharmacology. PK-Sim documentation: Simulations, section on the venous blood outputs (Peripheral Venous Blood). https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-documentation/pk-sim-simulations (accessed 2026-10-25).
<!-- end generated -->

## Sampling sites

The site of a blood sample is part of the metadata of a fit mapping: `PKPDMappingMetaData.sampling_site`, the vocabulary `SamplingSite` of `pkpdsim.metadata`. Plasma, serum and blood data map to the arm vein by default (decided by Matthias König: this is how the data are typically sampled); a venous, arterial, arterialised or hepatic venous observable only where the study states that site. Arterial blood and arterial plasma name their site: they are `ARTERIAL` without a statement, and any other stated site on them is a `ValueError`. The arm vein observer is the default fit target of the models built on the library. A mapping of urine, feces, feces and urine, bile, saliva or no tissue has no blood sample: its site is `NOT_APPLICABLE` without a statement (the field holds `None`, `resolved_sampling_site()` the site, which the filters, `to_dict` and `check_sampling_sites` use), and a blood site there (or `NOT_APPLICABLE` on plasma, serum or blood) is a `ValueError`. `pkpdsim.fit.check_sampling_sites(fit_mappings, sites)` checks the site of every mapping against a table `{observable id: SamplingSite}` and returns one problem per mismatch or unknown observable; a model test asserts it is empty.

| `SamplingSite` | value | observable |
|---|---|---|
| `ARM_VEIN` (default) | `arm vein` | `Carm_<s>`, `Carm_free_<s>`, `Carm_blood_<s>` |
| `VENOUS` | `venous` | `Cve_<s>`, `Cve_free_<s>`, `Cve_blood_<s>` |
| `ARTERIAL` | `arterial` | `Car_<s>`, `Car_free_<s>`, `Car_blood_<s>` |
| `ARTERIALISED` | `arterialised` | `Car_<s>`; the arm vein of the preset `arterialised` is the same |
| `HEPATIC_VEIN` | `hepatic vein` | `Chv_<s>` |
| `CAPILLARY` | `capillary` | no observer of its own |
| `NOT_APPLICABLE` | `not applicable` | urine and feces, `Aurine_<s>`, `Afeces_<s>` |

```python
from pkpdsim.metadata import (
    Dosing,
    Health,
    PKPDMappingMetaData,
    Route,
    SamplingSite,
    Tissue,
)

metadata = PKPDMappingMetaData(
    tissue=Tissue.PLASMA, route=Route.PO, dosing=Dosing.SINGLE, health=Health.HEALTHY
)
print(metadata.resolved_sampling_site())
```

```text
arm vein
```

The urine and the feces are the amounts of the sinks; the bile `Cbile_<s>` drains into the duodenum, so a biliary excretion is read from the feces (or from the bile when a study samples it). The urinary excretion rate `Rurine_<s>` integrates to the amount in the urine: a study which reports the excretion per collection interval maps to the amount `Aurine_<s>` (reset at the start of the interval) or to the rate.

## Observers of a definition

A definition adds observers of its own, `Observer(sid, formula, unit, name, notes)` in `PBPKDefinition.observers`: parameters of the top model, appended after the observers of the library, for what a model package reads which the library does not compute (the sum of two metabolites, a ratio, the amount in urine of a parent and its metabolites, a tissue concentration in the unit of a study). The formula is SBML L3 infix over the ids of the top model: the pools, the sinks, the observers of the library and of the definition, the covariates, the [parameters of the definition](substances.md#parameters-of-the-definition), the outputs of the organs (`Cout_<o>_<s>`, `Rurine_<s>`) and the readouts; and a species or an unbound concentration of a built organ by its flat id (`MU__Cti_<s>`, `LI__Cu_ic_<s>`), which the organ exposes (see below). `validate` checks every observer against the model composed without them, with the parser of libsbml: an id which is an id of the body or of another observer, a symbol which is no value of the top model (a submodel `LI` or its model definition is none), a call of a function which SBML L3 infix does not define (the models have no function definitions), observers which read each other in a cycle.

```python
from pkpdsim.pbpk import (
    Binding,
    Excretion,
    Kp,
    Metabolism,
    MM,
    Observer,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    U,
    compose,
)
from pkpdsim.pbpk.units import ML_PER_MIN_1P73M2

definition = PBPKDefinition(
    sid="observed",
    substances=[
        Substance(
            "drg",
            name="drug",
            mr=300.0,
            binding=Binding(fu=0.1),
            distribution=Kp(2.0, ftissue=1.0),
            excretion=Excretion(urine=True),
        ),
        Substance(
            "drg_m", name="metabolite", mr=316.0, excretion=Excretion(urine=True)
        ),
    ],
    chemistry=[
        Metabolism("DRG2M", Organ.LIVER, "drg", {"drg_m": 1.0}, MM(vmax=1e-3, km=0.01))
    ],
    routes=[Route.iv("drg")],
    observers=[
        Observer(
            "Aurine_drg_total",
            "Aurine_drg + Aurine_drg_m",
            U.mmole,
            "drug and metabolite in urine [mmole]",
        ),
        Observer(
            "Cmu_drg_unbound",
            "MU__Cu_ti_drg",
            U.mM,
            "unbound drug in the muscle tissue [mM]",
        ),
        Observer(
            "GFRbsa",
            "GFR / BSA * 1730 ml_per_min_1p73m2 * 1 m2 / 1 l_per_min",
            ML_PER_MIN_1P73M2,
            "GFR per 1.73 m2 body surface area [ml/min/1.73 m2]",
        ),
    ],
)
definition.validate()
model = compose(definition)
observed = {"Cmu_drg_unbound", "GFRbsa"}
print("rules:", sorted(observed & {str(r.variable) for r in model.rules}))
print(
    "initial assignments:",
    sorted(observed & {str(a.symbol) for a in model.assignments}),
)
```

```text
rules: ['Cmu_drg_unbound']
initial assignments: ['GFRbsa']
```

An observer which reads a state, directly or through other observers (a species, an observer of the library, an output of an organ or a PD component, a readout, a dynamic parameter (`PBPKDefinition.dynamic`) or `time`), is an assignment rule. An observer of constants only (`GFRbsa` above, the GFR per body surface area) is a constant with an initial assignment: its value is fixed before the simulation and follows the changes of an experiment before it (sbmlsim evaluates the initial assignments again), and an assignment rule of constants would slow the compile of the model (the structural guard of the load time, see [the overview](index.md)). An observer of a constant observer is one too.

The unit of an observer is a unit definition of the model. The units of `U` are those of every model; a unit outside `U` is added to the model whose observers use it, so a model without it keeps its unit definitions. `pkpdsim.pbpk.units` has the unit per body surface area of the clearances and the GFR, `ML_PER_MIN_1P73M2` (`ml_per_min_1p73m2`, ml/min per 1.73 m2, `1e-3/1.73 l/min/m2`); a quantity in l/min divided by `BSA` (m2) is `1730` times its number in this unit, in consistent units `GFR / BSA * 1730 ml_per_min_1p73m2 * 1 m2 / 1 l_per_min` (the build checks the units of every formula).

## An element of an organ for a PD component or an observer

A PD component reads the elements of the top model through its input ports (`PDComponent.inputs`), an observer through the ids of its formula; an element of an organ is local to its submodel and no element of the top model. An input or an observer which names an element of an organ by its flat id reads it all the same, by one rule: `PDComponent(model, inputs={"C_ged": "LI__Cic_ged"}, outputs={})` reads the concentration in the cells of the liver, as an imaging signal of the liver does, and so does `Observer("Cli_ged", "LI__Cic_ged", U.mM, "liver cells")`. The organ exposes the element as an output, a parameter of the organ with the assignment rule of the element, which replaces the parameter of the top model `ids.organ_species` gives (the organ token after the first letter of the local id), as the outflow `Cout_<o>_<s>` of an organ does; the component and the observer read that parameter, one output for both. The elements an organ exposes:

| element of the organ | local id | id of the top model | unit |
|---|---|---|---|
| a concentration species | `C<part>_<s>` (`Cic_ged`) | `C<o>_<part>_<s>` (`Cli_ic_ged`) | mM |
| an unbound concentration | `Cu_<part>_<s>` (`Cu_ic_ged`) | `C<o>_u_<part>_<s>` (`Cli_u_ic_ged`) | mM |
| an amount species (`hasOnlySubstanceUnits`, a pool of a custom reaction) | `A<x>` (`Apool_crtn`) | `A<o>_<x>` (`Amu_pool_crtn`) | of the species (mmole) |

The element keeps its id in the flat model (`LI__Cic_ged`), the flat model has one assignment rule more, in the liver:

```equations title="a species of the liver read by a PD component (local ids of the liver)"
C<o>_ic_<s> = Cic_<s>
```

In the flat model it is `Cli_ic_<s> = LI__Cic_<s>`. The unbound concentration is that of the organ, `fu/K*C` of the part (see [Standard organ](organ.md#kp-and-fu)), the input of a PD effect on the unbound drug in an organ.

The build refuses an input or an observer which names no species and no unbound concentration of a built organ (a compartment, a parameter, an unknown species), a species which the organ reads from the body (a pool species, read it in the top model) and a species whose id is not of the id scheme (a concentration which is no `C...`, an amount which is no `A...`), naming the reader and the reason; an id with the prefix of an organ which is not built is an unknown id of the body. A component or an observer which reads only the top model adds nothing to the organs.

## Readouts of the intestine

The lumen and the fluxes of the intestine are spread over its segments (`ids.GUT_SEGMENTS`, the twelve tokens of the segmented intestine in the order of the transit) or its single lumen. A model which reads their sums names them in `PBPKDefinition.readouts`, by id; the intestine adds each as a parameter with its assignment rule and an output port, which replaces the parameter of its id in the top model, as the kidney does with `Rurine_<s>`. The sums are taken over the species and the reactions of the intestine, so the readouts have the same ids and meaning for the segmented intestine and the single lumen:

| readout | id | the sum of |
|---|---|---|
| lumen | `Alumen_<s>` (`ids.lumen_amount`) | the dissolved amount of the lumen species, `sum Vlumen_<seg>*Clumen_<seg>_<s>` (mmole); the undissolved drug of a tablet (`PODOSE_<s>`, `Asolid_<seg>_<s>`) is not in it |
| absorption | `Rabs_<s>` (`ids.absorption_rate`) | the fluxes from the lumen into the wall of the substance (the mucosa, else the plasma of the wall): the absorption of every segment and an apical influx (mmole/min) |
| efflux | `Reff_<s>` (`ids.efflux_rate`) | the fluxes from the wall into the lumen, the apical efflux of every segment (mmole/min) |
| feces | `Rfeces_<s>` (`ids.feces_rate`) | the net flux of the intestine into the feces (mmole/min): the transit out of the last segment, the unabsorbed fraction `1 - fabs_<s>`, the undissolved drug, the bile of a substance without lumen species; it integrates to `Afeces_<s>` |

```python
from pkpdsim.pbpk import ids

oral = PBPKDefinition(
    sid="oral",
    substances=definition.substances,
    routes=[Route.oral("drg")],
    readouts=(ids.lumen_amount("drg"), ids.absorption_rate("drg")),
)
model = compose(oral)
print([str(p.sid) for p in model.parameters if str(p.sid) in oral.readouts])
```

```text
['Alumen_drg', 'Rabs_drg']
```

`validate` refuses a readout which is no readout of a substance of the definition and a readout given twice; the build refuses the readouts of an intestine without a lumen (the standard organs), `Alumen_<s>`, `Rabs_<s>` and `Reff_<s>` of a substance without lumen species and a readout without a flux (an efflux without an apical efflux of the small intestine), naming the readout. A model without readouts has none of them (its SBML is unchanged).
