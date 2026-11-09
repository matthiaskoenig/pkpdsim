# Organ variants

An organ variant is an alternative component for an organ with the boundary of the standard organ: blood in and blood out (plasma and blood cells, per substance), the flow, the organ volume, the hematocrit and the parameters of the substances and proteins it reads. The circulation sees nothing else, so swapping a variant changes nothing outside the organ: not the other organs, not the observers, not the experiments. A variant places the chemistry of its organ into its parts, so the chemistry is written once per organ whatever variant the model chooses.

## Catalogue

| organ | name | component | what it is |
|---|---|---|---|
| every organ | `standard` | `StandardOrgan()` | the [standard organ](organ.md), the default |
| liver | `zonated` | `ZonatedLiver(n_zones=10, profiles={}, geometry=Geometry.WELL_STIRRED, axial_diffusion=True, means=False)` | the sinusoid from periportal to pericentral as zones in series, the chemistry placed by volume or by a zonal profile of its protein; `Geometry.SINUSOID` adds the axial diffusion and the fenestrated endothelium of the zonation engine, `axial_diffusion=False` leaves the axial diffusion out, `means=True` adds the periportal and pericentral means |
| kidney | `nephron` | `NephronKidney()` | filtration into the tubule, the tubular segments, secretion into and reabsorption from the tubular fluid |
| small intestine | `segmented` | `SegmentedIntestine()` | the default intestine: the lumen segments of Thelen 2011, the mucosa, dissolution, absorption, bile; replaces small and large intestine with `GU` |
| small intestine | `single_lumen` | `SingleLumen()` | one lumen with first order absorption, the absorption of most hand-written models; replaces small and large intestine with `GU` |

A definition selects a variant by name or gives a component with its options:

```python
from pkpdsim.pbpk import Organ, Variants, ZonatedLiver

variants = Variants(
    organs={
        Organ.LIVER: ZonatedLiver(n_zones=4, profiles={"CYP3A4": (1.0, 1.0, 2.0, 4.0)}),
        Organ.KIDNEY: "nephron",
        Organ.SMALL_INTESTINE: "single_lumen",
    }
)
print({organ.name: getattr(v, "name", v) for organ, v in variants.organs.items()})
```

```text
{'LIVER': 'zonated', 'KIDNEY': 'nephron', 'SMALL_INTESTINE': 'single_lumen'}
```

The intestine variants are chosen on `Organ.SMALL_INTESTINE`; the large intestine has no variant of its own while an intestine replaces it. Coupled definitions (`PBPKDefinition.combine`) take the variant which one of them chooses for an organ the others leave at the default, and refuse two different choices.

## The equivalence contract

Every variant ships an equivalence test: in its limiting case it reproduces the standard organ, the outflow concentrations and the extraction at steady state.

| variant | limiting case | test |
|---|---|---|
| zonated liver | one zone is the standard liver, for every kinetics; uniform zones agree with the standard liver at a low extraction, and at a moderate extraction they extract as `n` well-stirred tanks in series; the sinusoid geometry without axial diffusion and with an open fenestral barrier is the well-stirred liver | `tests/pbpk/test_liver_zonated.py` |
| nephron kidney | without secretion and reabsorption the urinary excretion is `fu*GFR*C_p` of the standard kidney | `tests/pbpk/test_kidney_nephron.py` |
| intestines | the boundary of the standard organs of small and large intestine on the blood side; mass balance of an oral dose; the fraction absorbed of the compartmental absorption and transit model | `tests/pbpk/test_intestine.py` |

Every variant also has the same boundary as the standard organ (a test compares the inputs and outputs) and conserves mass in a closed system.

## Zonated liver

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/zonated.svg"
</div>

`ZonatedLiver` builds the liver as `n` well-stirred zones in series along the sinusoid, from the periportal zone `z1` to the pericentral zone `z<n>`, the pattern of the sinusoidal unit of the zonation models of pkdb_models (`misc/zonation`). Every zone has the parts of the standard liver which a substance needs, with the share `f_zone = 1/n` of their volumes and passive exchanges and the coefficients of the liver. The blood flows from the hepatic artery and the portal vein into the first zone, from zone to zone and from the last one into the venous pool; the shunt `f_shunt` bypasses all zones, the functional fraction `f_functional` scales the tissue and the capacity of every zone. The species and processes of a zone carry the zone in their id (`Cpl_z1_<s>`, `DRG2M_z1`), the outflow `Cout_<s>` is that of the last zone, so the hepatic and the arm vein observers are those of the standard liver.

Every metabolism and transport of the liver is placed into every zone with the share `w` of the capacity of the liver: by volume (`f_zone`), or by the profile of its protein, the relative abundance `p_<protein>_<zone>` per zone (a pericentral CYP3A4), `w = p*V_zone/sum(p*V_zone) = p/sum(p)`, so that the total capacity `Vmax*V` is that of the standard liver. A metabolism with a `profile` takes the profile of that name, any other process the profile of its protein. A profile is a value per zone, or a shape of the lobular position which the engine evaluates at the zone count (`pkpdsim.zonation.profiles`, so `Linear(1.0, 4.0)` is the four-fold pericentral gradient at any number of zones). A custom reaction is placed in every zone as well, with its rate law over the ids of the zone (`CUSTOM_z1`), and the constants of its rate law and the species it brings itself belong to the liver. The equations of a liver of three zones with a CYP3A4 profile `(1, 2, 4)`:

```equations title="zonated liver, three zones"
J_in_ar_pl_<s>: Car_<s> -> Cpl_z1_<s> = (1 - f_shunt) * Qha * (1 - HCT) * Car_<s>
J_sin_z1_pl_<s>: Cpl_z1_<s> -> Cpl_z2_<s> = (1 - f_shunt) * Qh * (1 - HCT) * Cpl_z1_<s>
J_sin_z2_pl_<s>: Cpl_z2_<s> -> Cpl_z3_<s> = (1 - f_shunt) * Qh * (1 - HCT) * Cpl_z2_<s>
J_out_pl_<s>: Cpl_z3_<s> -> Cve_<s> = (1 - f_shunt) * Qh * (1 - HCT) * Cout_<s>
Cout_<s> = Cpl_z3_<s>
Vti_z1 = Vli * (1 - f_vas_li) * f_functional * f_zone
J_ti_z1_<s>: Cpl_z1_<s> -> Cti_z1_<s> = PS_<s> * f_zone * fu_<s> * (Cpl_z1_<s> - Cti_z1_<s> / Kp_<s>)
w_cyp3a4_z1 = p_cyp3a4_z1 / (p_cyp3a4_z1 + p_cyp3a4_z2 + p_cyp3a4_z3)
DRG2M_z1: Cti_z1_drg -> Cpl_z1_drg_m = f_cyp3a4 * fimp_cyp3a4 * f_functional * w_cyp3a4_z1 * DRG2M_Vmax * Vli * (fu_drg / Kp_drg * Cti_z1_drg) / (DRG2M_Km + fu_drg / Kp_drg * Cti_z1_drg)
```

The zones of the default geometry are well-stirred tanks (`Geometry.WELL_STIRRED`). The sinusoid geometry (`Geometry.SINUSOID` of `pkpdsim.zonation.settings`) takes the two transport steps the well-stirred zone leaves out from the zonation engine (`pkpdsim.zonation`), with the geometry of the liver of the species of the build: the diffusion along the sinusoid between the zones (`D_z1_pl_<s>`, in the plasma and in the interstitium, `D_<s>*V*n^2/L_sin^2` over the distance `L_sin/n` of two zones) and the fenestrated endothelium, whose conductance `PS_vas_<zone>_<s> = P_fen_<s>*eps_fen_<zone>*A_pd_<zone>` replaces the generic endothelium `PS_vas` of the standard organ (the fenestral porosity rises from the periportal to the pericentral end, 0.06 to 0.08 for the human of the table of the engine, Braet and Wisse 2002; the lateral area `2*V_lumen/r_sin` follows the lumen of the zone, and the permeability of an open fenestra is the Renkin 1954 hindrance of the molecule). The layers of the engine are the parts of the liver: the plasma of the lumen is the plasma, the erythrocytes are the blood cells, the space of Disse is the interstitium and the hepatocyte is the intracellular space, so a substance whose form has no interstitium has no barrier either and the endothelial cell and the bile canaliculus are no part of the PBPK liver.

The two steps are of very different weight, and the axial diffusion is the one which changes the simulation. For a small molecule its conductance between two zones is of the order of the plasma flow through the chain (about 1.5 against 1 l/min for a drug of 0.3 kDa in ten zones, the Peclet number of about 6 against the `2n = 20` of a substance which does not diffuse), so it back-mixes the zones and partly undoes the zonation: such a drug deviates in the venous plasma from the well-stirred liver of the same zone count by about 7e-03 with the diffusion on and by 2e-05 with it off, and the gradient along the sinusoid is about a third flatter. A modeller who chooses a zone count, or who wants the gradient of the lobule and not its mean, needs that number; `axial_diffusion=False` leaves the term out. The fenestral path limits nothing for such a molecule: `PS_vas_<zone>_<s>` is hundreds of times the plasma flow, and the basolateral `PS_<s>` is the step which limits the transfer into the tissue. The lumen of a zone of this model is the whole vascular space of the liver (`f_vas`), where the engine takes the sinusoidal lumen only, so the lateral area and with it the conductance are about 2.1 times the engine's for the same liver, and the velocity of the same volume correspondingly lower; neither number limits anything, and the notes of the model state the convention.

The conductance of the fenestral path of a zone is the named parameter `PS_vas_<zone>_<s>` (l/min), a constant with an initial assignment as `PS_vas` of the standard organ is, so it follows the liver volume of the individual and a sensitivity analysis scans it and a fit changes it. A molecule whose radius reaches the radius of a fenestra is refused: the PBPK liver has no endothelial compartment and therefore no transcellular path through the endothelium, and the well-stirred geometry holds such a molecule. The chain of zones is upwind, so its convection adds the numerical dispersion `v*L_sin/(2n)` to the axial diffusion; the volume and the flow of the simulated individual are parameters of the model and no numbers at build time, so the notes of the model give the velocity of the blood, the dispersion and the Peclet number of the chain at one individual, the reference male of the species (the liver mass of the geometry table of the engine is his, so the vascular fraction and the hepatic blood flow are his too, whatever the sex of the build), and the notes of `D_<s>` the Peclet number with that coefficient.

`means=True` adds the two means of the zonation engine (`pkpdsim.zonation.observers`) to every zonated part and substance: `C<part>_pp_<s>` and `C<part>_pc_<s>` (mM), the mean concentration of the first (periportal) and the last (pericentral) third of the sinusoid, every zone weighted by its overlap with the third (`pkpdsim.zones.third_weights`, the weighting every zonated model of pkpdsim shares), so a zonated model of either library reports the same two numbers. They are assignment rules which read the zones and which nothing reads, and they are on request, as the observers of the engine are added on request and the readouts of the intestine are.

The zonated liver has two limits in which it is a model of the literature again, and a test of each. With one zone it is the standard liver, the well-stirred liver, for every kinetics. With the fenestral barrier open and the axial diffusion off (`axial_diffusion=False`, which emits no diffusion term at all: a `Substance` has no coefficient of 0, so the limit belongs to the geometry) the sinusoid geometry is the well-stirred geometry of the same number of zones. With `n` zones and linear kinetics the availability of a flow-limited liver at steady state is `(1 + x/n)^-n` with `x = fu*CLint/(Q*BP)` instead of the well-stirred `1/(1 + x)`: with the same intrinsic clearance the zonated liver extracts more, and for large `n` it approaches the parallel tube, `exp(-x)`. A protein which a substance inactivates or induces has one turnover state for the liver, driven by the unbound concentration in the cells of the zones weighted by the placement of the protein (exact for the same concentration in every zone, the mean of the zones along a gradient).

## Nephron kidney

`NephronKidney` builds the kidney with the tubule of the nephron. The unbound afferent (arterial) plasma concentration is filtered into the lumen of the proximal tubule, as in the standard kidney, and the tubular fluid flows through the lumens of the segments in series, proximal tubule `pt`, loop of Henle `lh`, distal tubule `dt` and collecting duct `cd`, into the urine. Every segment reabsorbs water: the flow leaving a segment is the fraction `f_tub_<seg>` of the GFR (Scotcher 2016), so the tubular flows follow the GFR of the individual and the renal function, and the drug concentrates as the water is reabsorbed. The lumens are well-stirred compartments of constant volume (Huang 2018, both kidneys), lumen and not tissue: the blood and tissue of the kidney keep the volumes of the standard kidney.

```equations title="nephron kidney"
J_gfr_<s>: Cpl_<s> -> Ctub_pt_<s> = fu_<s> * GFR * Car_<s>
J_tub_pt_<s>: Ctub_pt_<s> -> Ctub_lh_<s> = GFR * f_tub_pt * Ctub_pt_<s>
J_tub_lh_<s>: Ctub_lh_<s> -> Ctub_dt_<s> = GFR * f_tub_lh * Ctub_lh_<s>
J_tub_dt_<s>: Ctub_dt_<s> -> Ctub_cd_<s> = GFR * f_tub_dt * Ctub_dt_<s>
J_tub_cd_<s>: Ctub_cd_<s> -> Aurine_<s> = GFR * f_tub_cd * Ctub_cd_<s> * (1 - f_reabsorbed_<s>)
J_reab_<s>: Ctub_cd_<s> -> Cpl_<s> = GFR * f_tub_cd * Ctub_cd_<s> * f_reabsorbed_<s>
J_out_pl_<s>: Cpl_<s> -> Cve_<s> = Qki * (1 - HCT) * Cout_<s> - fu_<s> * GFR * Cout_<s> * (1 - f_reabsorbed_<s>)
Rurine_<s> = J_tub_cd_<s>
```

Without secretion and reabsorption every filtered molecule reaches the urine: the urinary excretion at steady state is the filtration of the standard kidney, whatever the volumes and flows of the segments; the tubule delays the excretion by minutes. The plasma leaving the kidney lacks the net filtrate of the standard kidney: at steady state the reabsorption returns the fraction reabsorbed of the filtrate into the plasma of the kidney, so the kidney plasma and the clearance of the arterial plasma are those of the standard kidney. The fraction reabsorbed (`Excretion.fraction_reabsorbed`) is an empirical passive reabsorption at the end of the collecting duct, which acts on secreted drug as well, the lumped `(fu*GFR + CL_sec)*(1 - F_reabs)` of the literature. A [threshold reabsorption](organ.md#threshold-reabsorption) (`Excretion.reabsorption=ThresholdReabsorption(...)`) returns the fraction of the standard kidney at the afferent plasma, `f(Car) = 1/(1 + (fu*|Car|/reab_threshold)^reab_n)^(1/reab_n)`, of the flow out of the collecting duct up to the filtered load `fu*GFR*Car` (`Car_<s>` a modifier of both reactions out of the collecting duct): it never exceeds the filtration, also with a secretion into the tubule, whose drug leaves into the urine, and at steady state without a secretion the circulation sees the standard kidney.

```equations title="nephron kidney, threshold reabsorption"
J_tub_cd_<s>: Ctub_cd_<s> -> Aurine_<s> = GFR * f_tub_cd * Ctub_cd_<s> - min(GFR * f_tub_cd * Ctub_cd_<s>, fu_<s> * GFR * Car_<s>) * 1 / (1 + (fu_<s> * abs(Car_<s>) / reab_threshold_<s>)^reab_n_<s>)^(1 / reab_n_<s>)
J_reab_<s>: Ctub_cd_<s> -> Cpl_<s> = min(GFR * f_tub_cd * Ctub_cd_<s>, fu_<s> * GFR * Car_<s>) * 1 / (1 + (fu_<s> * abs(Car_<s>) / reab_threshold_<s>)^reab_n_<s>)^(1 / reab_n_<s>)
```
 An apical efflux of the kidney secretes into the proximal tubule, an apical influx reabsorbs from it, with the inhibitors of its protein acting from the lumen. The urinary excretion rate `Rurine_<s>` is the outflow of the collecting duct, the output of the standard kidney with the same id.

<!-- generated: nephron-segments (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/kidney_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| segment | flow leaving, fraction of GFR [-] | volume [l] |
| --- | ---: | ---: |
| proximal tubule | 0.36 | 0.0915 |
| loop of Henle | 0.2 | 0.0054 |
| distal tubule | 0.09667 | 0.0194 |
| collecting duct | 0.008333 | 0.005 |

Sources of the columns: flow: `Scotcher2016`; volume: `Huang2018`, `ModelChoice`.

??? note "Sources"

    - `Huang2018`: Huang W, Isoherranen N. Development of a dynamic physiologically based mechanistic kidney model to predict renal clearance. CPT Pharmacometrics Syst Pharmacol 7(9):593-602, 2018. [doi:10.1002/psp4.12321](https://doi.org/10.1002/psp4.12321)
    - `ModelChoice`: Choice of a library of pkpdsim (PBPK or machine perfusion), not a value of the literature; the notes of the value give the reasoning.
    - `Scotcher2016`: Scotcher D, Jones C, Rostami-Hodjegan A, Galetin A. Novel minimal physiologically-based model for the prediction of passive tubular reabsorption and renal excretion clearance. Eur J Pharm Sci 94:59-71, 2016. [doi:10.1016/j.ejps.2016.03.018](https://doi.org/10.1016/j.ejps.2016.03.018)
<!-- end generated -->

## Intestine variants

The segmented intestine and the single lumen replace small and large intestine with one submodel `GU`, with the boundary of the two standard organs on the blood side; they own the intestinal lumen, apply the oral doses and take the bile. They are described with the routes on [Intestine and routes](intestine.md). Selecting the standard organs for the intestine (`Variants(organs={Organ.SMALL_INTESTINE: "standard"})`) builds no lumen: a model with an oral route or an apical transport of the intestine then fails with an error.

## Writing a variant

A variant is an `OrganComponent`: an object with a `name` and a method `build(ctx: OrganContext) -> tuple[ModelDefinition, Boundary]`. The context gives the organ, the substances with their parts in the organ (`parts_for`), the chemistry, the inhibitions and inductions and the routes; the boundary names, for every input of the model definition, the element of the body which replaces it, and the outputs which replace parameters of the body (`Cout_<s>` -> `Cout_<o>_<s>`). The builder generates the ports and replacements from the boundary; a variant writes no `ReplacedElement`.

An intestine which replaces several organs with one submodel is a `GutComponent`: in addition it names the organs it replaces (`replaces`, the small intestine first) and gives the parameters of the top model it reads (`body_objects(definition)`, e.g. `ka_<s>` and `FED`), as the segmented intestine and the single lumen do.

The base classes do the bookkeeping: `pkpdsim.pbpk.organ_builder.OrganBuilder` collects the inputs with their ports, the compartments, species, rules, initial assignments, constants with their provenance and the activities of the proteins, and `pkpdsim.pbpk.organ.StandardBuilder` builds the parts, flows and exchanges of the standard organ (with zones, as the zonated liver uses it); `pkpdsim.pbpk.placement` places the chemistry. A new variant derives from them, keeps the boundary of the standard organ, places the chemistry by volume fraction and ships its equivalence test.

A model package can write its own variant against the same interface first and pass the instance, `Variants(organs={Organ.LIVER: MyLiver()})`; a variant which several models need moves into the library and is registered by name (`pkpdsim.pbpk.variants.register(organ, name, component)`). The registry stores component instances, which are immutable dataclasses.

## Limitations

- The variants refuse a `CustomReaction` of their organ (with an error naming the process): an arbitrary formula cannot be split over zones or segments. A model with custom liver chemistry uses the standard liver until a hook for it exists.
- A zonal profile of zeros gives no placement weights: the build refuses an all-zero profile, but a profile changed to zeros during a simulation gives NaN.
- The volume of the collecting duct of the nephron is a choice of the library (0.005 l, no source could be verified), and the tubular lumens do not scale with the individual.
