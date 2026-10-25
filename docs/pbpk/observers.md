# Observers

The observers are assignment rules, no states: they read the pool species and the outflow concentrations of the organs and change nothing, so they conserve mass exactly. The forearm compartment of some hand-written models lost mass and was removed from three of them; the arm vein of the library is an observer.

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
| urine, feces | `Aurine_<s>`, `Afeces_<s>` | the amounts in the sinks |
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

The site of a blood sample is part of the metadata of a fit mapping: `PKPDMappingMetaData.sampling_site`, the vocabulary `SamplingSite` of `pkpdsim.metadata`. Plasma and blood data map to the arm vein by default (decided by Matthias König: this is how the data are typically sampled); a venous, arterial, arterialised or hepatic venous observable only where the study states that site. The arm vein observer is the default fit target of the models built on the library.

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
print(metadata.sampling_site)
```

```text
arm vein
```

The urine and the feces are the amounts of the sinks; the bile `Cbile_<s>` drains into the duodenum, so a biliary excretion is read from the feces (or from the bile when a study samples it).
