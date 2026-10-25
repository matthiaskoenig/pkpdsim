# Standard organ

Every organ is a comp submodel built from the standard organ, unless the definition selects a [variant](variants.md). It is the 4-compartment organ of PK-Sim (Willmann 2005, Schmitt 2008, Kuepfer 2016): plasma, blood cells, interstitial and intracellular space, of which only the parts a substance needs are built. The organ owns its blood flows and its chemistry, and reads everything of the body (volume, flow, hematocrit, fractions, the parameters of the substances, the pool and sink species) through ports.

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/organ.svg"
</div>

## Parts

| part | compartment | volume | species | built for a substance |
|---|---|---|---|---|
| plasma | `Vpl` | `V<o>*f_vas*(1 - HCT)` | `Cpl_<s>` | always, but not in a flow-limited organ |
| blood cells | `Vbc` | `V<o>*f_vas*HCT` | `Cbc_<s>` | with kinetic erythrocytes |
| interstitial | `Vis` | `V<o>*f_int` | `Cis_<s>` | in the full form |
| intracellular | `Vic` | `V<o>*(1 - f_vas - f_int)` | `Cic_<s>` | in the full form |
| tissue | `Vti` | `V<o>*(1 - f_vas)` | `Cti_<s>` | in the merged tissue (interstitial and intracellular) |
| organ | `Vor` | `V<o>` | `Cor_<s>` | in the flow-limited organ (vascular and tissue) |

Every species is a concentration (mmole/l), every reaction an amount per time (mmole/min). The ids are local to the organ and prefixed by its submodel in the flat model (`MU__Cpl_api`); the inputs of the body keep the ids of the body (`Qmu`, `HCT`, `fu_api`, `Car_api`). The equations below are those of the muscle `mu` (`<o>`) for a substance `<s>`, copied from the generator; `tests/docs/test_pbpk_equations.py` builds toy models and checks every equation of this page against them.

```equations title="volumes of the parts"
Vpl = V<o> * f_vas_<o> * (1 - HCT)
Vbc = V<o> * f_vas_<o> * HCT
Vis = V<o> * f_int_<o>
Vic = V<o> * (1 - f_vas_<o> - f_int_<o>)
Vti = V<o> * (1 - f_vas_<o>)
Vor = V<o>
```

## Blood flows

The organ takes the plasma of the arterial pool (the venous pool for the lung, the hepatic artery and the portal vein for the liver) and returns the plasma leaving it, `Cout_<s>`, into the venous pool (the portal vein for the portal organs, the arterial pool for the lung). The outflow concentration is an output of the organ which replaces `Cout_<o>_<s>` of the top model, where the observers read it.

```equations title="plasma"
J_in_ar_pl_<s>: Car_<s> -> Cpl_<s> = Q<o> * (1 - HCT) * Car_<s>
J_out_pl_<s>: Cpl_<s> -> Cve_<s> = Q<o> * (1 - HCT) * Cout_<s>
Cout_<s> = Cpl_<s>
```

With kinetic erythrocytes the blood cells flow alongside and exchange with the plasma; with erythrocytes in equilibrium the plasma species carries the blood, `Q*BP*C_p`, in the capacity `V_vas*BP`, and no blood cell species is built.

```equations title="blood cells (kinetic erythrocytes)"
J_in_ar_bc_<s>: Car_bc_<s> -> Cbc_<s> = Q<o> * HCT * Car_bc_<s>
J_out_bc_<s>: Cbc_<s> -> Cve_bc_<s> = Q<o> * HCT * Cout_bc_<s>
J_bc_<s>: Cpl_<s> -> Cbc_<s> = kex_bc_<s> * Vbc * (Cpl_<s> - Cbc_<s> / Kbc_<s>)
Cout_bc_<s> = Cbc_<s>
```

## Collapsing

The flags of a substance and the chemistry of the organ decide which parts are built for the substance in the organ (`parts_for`), so a simple drug has about as many species as in the hand-written models:

| situation | parts built for the substance in the organ |
|---|---|
| no erythrocyte distribution | no blood cell species; the flows carry plasma only, `BP = 1 - HCT` |
| erythrocytes in equilibrium | no blood cell species; the flows carry `Q*BP*C_p`, the capacity is `V_vas*BP` |
| erythrocytes kinetic | blood cell species with exchange |
| no tissue distribution, or an organ the substance does not enter | plasma (and blood cells) only |
| `Mode.PERFUSION` | one species for the whole organ, flow-limited, the outflow `C_organ/K_organ` |
| `Mode.TISSUE` (the default with a Kp) | plasma and one tissue species with the exchange `PS*fu*(C_p - C_t/Kp)` |
| `Mode.PERMEABILITY`, or an active transport of the substance in the organ | interstitial and intracellular species |

A transport of a substance in an organ (basolateral, or an apical process) builds the full form there, also for a substance without distribution. The product of a metabolism is placed into the equivalent tissue part of the product, else into the plasma.

The collapsed forms are the limiting cases of the full form, not other models: `tests/pbpk/test_equivalence.py` shows that the merged tissue with a large exchange is the flow-limited organ, the full form with a large `PS_vas` the merged tissue, and fast kinetic erythrocytes the erythrocytes in equilibrium.

### Merged tissue

The default with a Kp: the interstitial and the intracellular space are one tissue species with the whole-tissue Kp. The exchange of the organ is its share of the whole-body exchange `ftissue_<s>` by its tissue volume, relative to the reference tissue volume of the organs the substance enters (a `ps` of the organ replaces it).

```equations title="merged tissue (Mode.TISSUE)"
J_ti_<s>: Cpl_<s> -> Cti_<s> = PS_<s> * fu_<s> * (Cpl_<s> - Cti_<s> / Kp_<s>)
PS_<s> = ftissue_<s> * V<o> * (1 - f_vas_<o>) / Vti_total_ref_<s>
Cu_ti_<s> = fu_<s> / Kp_<s> * Cti_<s>
```

### Flow-limited organ

One species for the vascular space and the tissue together, in equilibrium with the blood leaving the organ: the outflow is the organ concentration over the organ to plasma coefficient `K_or_<s>`, exactly flow-limited, without a fast rate constant (no stiffness). The coefficient of adipose is its override `Kp_ad_<s>` if the substance gives one.

```equations title="flow-limited organ (Mode.PERFUSION, erythrocytes in equilibrium)"
J_in_ar_pl_<s>: Car_<s> -> Cor_<s> = Q<o> * BP_<s> * Car_<s>
J_out_pl_<s>: Cor_<s> -> Cve_<s> = Q<o> * BP_<s> * Cout_<s>
K_or_<s> = f_vas_<o> * BP_<s> + (1 - f_vas_<o>) * Kp_<s>
Cout_<s> = Cor_<s> / K_or_<s>
Cu_or_<s> = fu_<s> / K_or_<s> * Cor_<s>
```

### Full form

Interstitial and intracellular species (PK-Sim 4-compartment structure, in the form of the OSP report, section 2.1). The exchange between plasma and interstitium crosses the endothelium, `PS_vas = k_vas*V<o>` with `k_vas` 1000 1/min, a choice of the library: the endothelium is not limiting for small molecules, as in the default of PK-Sim. The exchange between interstitium and cells is the share of `ftissue_<s>` (or the `ps` of the organ). The interstitial to plasma coefficient follows from fu (Schmitt 2008), and the intracellular one is derived so that both together are the Kp of the substance, `(f_int*K_int_pls + f_cell*K_cell_pls)/(f_int + f_cell) = Kp`; its floor `1e-6` keeps it positive should a fit or another individual reach the interstitial share.

```equations title="full form (Mode.PERMEABILITY)"
J_is_<s>: Cpl_<s> -> Cis_<s> = PS_vas * fu_<s> * (Cpl_<s> - Cis_<s> / K_int_pls_<s>)
J_ic_<s>: Cis_<s> -> Cic_<s> = PS_<s> * (fu_<s> / K_int_pls_<s> * Cis_<s> - fu_<s> / K_cell_pls_<s> * Cic_<s>)
PS_vas = k_vas * V<o>
K_int_pls_<s> = (f_water_int + r_protein_int_pls * (1 / fu_<s> - f_water_pls)) * fu_<s>
K_cell_pls_<s> = piecewise((Kp_<s> * (f_int_<o> + (1 - f_vas_<o> - f_int_<o>)) - f_int_<o> * K_int_pls_<s>) / (1 - f_vas_<o> - f_int_<o>), (Kp_<s> * (f_int_<o> + (1 - f_vas_<o> - f_int_<o>)) - f_int_<o> * K_int_pls_<s>) / (1 - f_vas_<o> - f_int_<o>) > 1e-6, 1e-6)
Cu_is_<s> = fu_<s> / K_int_pls_<s> * Cis_<s>
Cu_ic_<s> = fu_<s> / K_cell_pls_<s> * Cic_<s>
```

The water fractions of interstitium and plasma and the protein ratio of the interstitial fluid (`f_water_int` 0.935, `f_water_pls` 0.926, `r_protein_int_pls` 0.37) are those of Schmitt 2008, the values of the plan which are not yet checked against the table of the paper (`data/distribution.toml`).

## Kp and fu

The Kp of a substance (or the override of an organ, or the value of a Kp method) is the whole-tissue (interstitial and intracellular) to plasma coefficient in every form, so a fitted Kp means the same in every mode. Every rate acts on unbound concentrations: `fu*C_p` in plasma, `fu/Kp*C_t` in the merged tissue, `fu/K_or*C_or` in the flow-limited organ, `fu/K_int_pls*C_is` and `fu/K_cell_pls*C_ic` in the full form. These unbound concentrations are the observers `Cu_<part>_<s>` of the organ (`LI__Cu_ti_api` in the flat model).

## Liver

The liver is fed by the hepatic artery (`Qha`) and the portal vein (`Qpo`) and drains into the venous pool. Two parameters of the top model, healthy 0 and 1, carry the cirrhosis of the [impairment](impairment.md): the portosystemic shunt `f_shunt`, the fraction of the arterial and portal inflow which flows past the liver into the venous pool, and the functional fraction `f_functional`, which scales the tissue volumes of the liver and the capacity of its processes. An experiment may set both directly. The bile is a sink (`Cbile_<s>`) which drains into the duodenum of the intestine.

```equations title="liver, merged tissue"
J_in_ar_pl_<s>: Car_<s> -> Cpl_<s> = (1 - f_shunt) * Qha * (1 - HCT) * Car_<s>
J_in_po_pl_<s>: Cpo_<s> -> Cpl_<s> = (1 - f_shunt) * Qpo * (1 - HCT) * Cpo_<s>
J_shunt_ar_pl_<s>: Car_<s> -> Cve_<s> = f_shunt * Qha * (1 - HCT) * Car_<s>
J_shunt_po_pl_<s>: Cpo_<s> -> Cve_<s> = f_shunt * Qpo * (1 - HCT) * Cpo_<s>
J_out_pl_<s>: Cpl_<s> -> Cve_<s> = (1 - f_shunt) * Qh * (1 - HCT) * Cout_<s>
Vti = Vli * (1 - f_vas_li) * f_functional
PS_<s> = ftissue_<s> * Vli * (1 - f_vas_li) * f_functional / Vti_total_ref_<s>
```

The hepatic venous concentration is an observer of the top model, the outflow of the liver mixed with the shunt (see [Observers](observers.md)).

## Kidney

The kidney filters the unbound plasma concentration of a substance with urinary excretion into the urine, with the GFR of the physiology (scaled with the individual and the renal function):

```equations title="kidney"
J_gfr_<s>: Cpl_<s> -> Aurine_<s> = fu_<s> * GFR * Cpl_<s> * (1 - f_reabsorbed_<s>)
Rurine_<s> = J_gfr_<s>
```

The urinary excretion rate `Rurine_<s>` (mmole/min) is the sum of the fluxes of the kidney into the urine, the filtration and every apical efflux (`J_gfr_<s> + SEC` with a secretion `SEC`), for every substance with urinary excretion or an apical transport of the kidney; it is an output of the kidney which replaces `Rurine_<s>` of the top model, so experiments and observers read it there.

Tubular secretion is an apical efflux of the kidney into the urine (see [Chemistry](chemistry.md)); the standard kidney refuses an uptake from the urine, which is cumulative. Every metabolism and transport of the kidney runs with the tubular capacity `1 - f_tubular + f_tubular*f_renal_function`, which follows the renal function (the intact nephron hypothesis, see [Impairment](impairment.md#chronic-kidney-disease)). The [nephron kidney](variants.md#nephron-kidney) has the tubule with its segments, for secretion into and reabsorption from the tubular fluid.

## Boundary

The organ declares its boundary (`Boundary`): every input (`Qmu`, `HCT`, `fu_api`, `Car_api`, ...) with the element of the body which replaces it, and its outputs (`Cout_<s>` replaces `Cout_<o>_<s>` of the top model, the kidney's `Rurine_<s>` replaces `Rurine_<s>`). The builder generates the ports and the replacements from it. A variant has the same boundary, so the circulation sees nothing else and swapping a variant changes nothing outside the organ (see [Organ variants](variants.md)).
