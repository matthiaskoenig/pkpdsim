# PBPK library

`pkpdsim.pbpk` builds the whole-body PBPK model of a substance from a short definition. A model package writes only what is specific to its substances: the substances with their mechanism flags, the chemistry of its organs (metabolism, transport), its choices (organ variants, routes, population, sex) and its PD components. The library adds a standard human physiology (ICRP 89), a standard organ for every organ, the alternative organ variants, the routes and sinks, the impaired organ function, the observers and the annotation of every element, and writes the SBML: a comp model with one submodel per organ and the flat model which is simulated.

The pages of this section are for modellers first: this page and the [tutorial](tutorial.md) show how a model is built and how the parts fit together; the pages of the parts say what each part does, which equations it generates and where its numbers come from.

## How a model is built

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/composition.svg"
</div>

1. The model package describes its substances (`Substance`) and their chemistry (`Metabolism`, `Transport`, `CustomReaction`), and makes its choices in one `PBPKDefinition`: the organ variants (`Variants`, the standard organ everywhere else), the routes (`Route`), the PD components (`PDComponent`), the sex and the population.
2. `build(definition, output_dir)` in the factory of the model validates the definition and composes the comp model: a top model with the physiology, the blood pools, the parameters of the substances, the sinks, the dose states and the observers (of the library and of the definition), the creators of the definition as model history, and one submodel per organ (`LI`, `KI`, `GU`, ...), the standard organ or a variant, which owns everything inside its boundary. Every port and every replacement is generated; the model package never writes a `ReplacedElement`.
3. `build` writes the comp model `<sid>.xml` and flattens it into `<sid>_flat.xml`, the model which `MODEL_PATH` points at and which every simulation runs. Simulation never resolves comp.
4. The experiments of the model set the covariates (`BW`, `HEIGHT`), the doses (`PODOSE_<s>`), the impairment (`CP_severity`, `f_renal_function`) and the activities of the proteins (`LI__f_cyp2c9`) as they do today, and read the observers (`Carm_<s>`, `Aurine_<s>`).

## A definition

A simple substance needs a short specification and its chemistry: an oral drug, bound to albumin, with one Kp, renal filtration and a hepatic CYP3A4 metabolism into a metabolite.

```python
from pkpdsim.pbpk import (
    MM,
    Binding,
    Excretion,
    Kp,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
)

definition = PBPKDefinition(
    sid="api_body",
    substances=[
        Substance(
            "api",
            name="apixaban",
            mr=459.5,
            binding=Binding(fu=0.13, protein="albumin"),
            distribution=Kp(1.2, ftissue=5.0),
            excretion=Excretion(urine=True),
        ),
        Substance("api_m2", name="O-desmethyl apixaban", mr=445.5),
    ],
    chemistry=[
        Metabolism(
            "API2M2",
            Organ.LIVER,
            "api",
            products={"api_m2": 1.0},
            kinetics=MM(vmax=0.01, km=0.05),
            protein="CYP3A4",
        ),
    ],
    routes=[Route.iv("api"), Route.oral("api")],
)
definition.validate()
print([s.sid for s in definition.substances], definition.sex, definition.population)
```

```text
['api', 'api_m2'] male european_icrp
```

Everything else is a default: the standard organ in every organ (the segmented intestine for small and large intestine), the male ICRP 89 reference individual, the arm vein preset of PK-Sim. The [tutorial](tutorial.md) builds and simulates such a model step by step.

## How the parts fit

| part | what it is | page |
|---|---|---|
| physiology | the organs, their volumes, blood flows and blood (ICRP 89), the topology of the circulation, the scaling to the individual from `BW` and `HEIGHT`, the populations | [Physiology](physiology.md) |
| standard organ | plasma, blood cells, interstitial and intracellular space of every organ; only the parts a substance needs are built (collapsing) | [Standard organ](organ.md) |
| substances | the mechanism flags of a substance: protein binding, erythrocytes, tissue distribution, excretion, inhibition and induction | [Substances](substances.md) |
| chemistry | the processes of an organ, defined once per organ: metabolism, transport, custom reactions; the proteins with their activity factors | [Chemistry](chemistry.md) |
| organ variants | alternatives of an organ with the boundary of the standard organ: the zonated liver, the nephron kidney, the intestines | [Organ variants](variants.md) |
| intestine and routes | the lumen segments, dissolution, absorption, the mucosa, the bile and the enterohepatic circulation; the routes and the sinks | [Intestine and routes](intestine.md) |
| impairment | liver cirrhosis, chronic kidney disease and heart failure as parameter vectors on the physiology | [Impairment](impairment.md) |
| observers | the plasma, free and blood concentrations, the hepatic and the arm vein, the urinary excretion rate, the observers of a definition, the sampling sites of the fit mappings | [Observers](observers.md) |
| annotation | the SBO, PBPKO, UBERON, NCIT, UniProt and ChEBI terms of every element | [Annotation](annotation.md) |

The [migration guide](migration.md) explains the move of a hand-written model of pkdb_models onto the library, with the glimepiride model as example.

## What is generated

Every quantity which depends on parameters only (the physiology, the impairment, the compartments, `fu`, `BP`, the coefficients and exchanges of the organs) is a constant with an initial assignment: an experiment which changes `BW` before the simulation changes every organ volume and flow with it, because sbmlsim evaluates the initial assignments again after the changes of an experiment. A change during a simulation is not followed, unless the definition names the parameter in `PBPKDefinition.dynamic` (an impairment changed in the middle of a study, a PD output): then every quantity downstream of it is an assignment rule. The start of the species of an endogenous substance (`Substance.baseline`, and the starts of a definition at the steady state of a production, `PBPKDefinition.starts`) is an initial assignment as well, which stays one under a dynamic parameter (see [Substances](substances.md#endogenous-baseline)). Only what reads a species is always an assignment rule: the outflow concentrations of the organs, the urinary excretion rate, the readouts and the observers; an observer of a definition which reads constants only is a constant with an initial assignment (see [Observers](observers.md#observers-of-a-definition)). A model loaded with plain roadrunner sets a covariate with `r.setValue("init(BW)", 100.0)` and resets, so that the initial assignments are evaluated again.

Sex is a choice of the build (`PBPKDefinition.sex`): it selects the male or female reference values when the model is built, and a study of both sexes uses one model per sex. A sex switch inside the SBML made the models load about ten times slower.

## Values and their sources

Every number which the library writes into a model is a `Value`: the number, its unit and its source, a key of `data/sources.toml` (citation and DOI), and the generated SBML parameter carries the citation in its notes. The data files of the library (`pkpdsim/pbpk/data/*.toml`) give every number with its source; a number which no publication gives is a choice of the library with the source `ModelChoice` and a note which gives the reasoning. The tables of these pages are generated from the data files.

A number of a definition is a plain float (the unit of the model is implied), a `Value`, or a `Formula` over ids of the top model, which ties it to other parameters (a metabolite with the Kp of its parent, see [Substances](substances.md#values-fits-and-ids)); a fitted value enters with its uncertainty from the record of its fit, `fitted(FITS.parameter("Kp_api", unit=U.dimensionless))`, and keeps the id which `FITS.parameter` and `check_model` resolve. The units are those of the models: time in min, amount in mmole, volume in litre, area in m2.

## The id scheme

One id scheme for every model, built by the functions of `pkpdsim.pbpk.ids`, so that organs, observers, experiments and fits agree on the ids. The ids of the hand-written models are kept where they fit the scheme (`Cve_<s>`, `Aurine_<s>`, `BW`, `f_renal_function`, the submodel prefixes `LI__`, `KI__`, `GU__`). A local id of an organ is prefixed by its submodel in the flat model (`LI__Cpl_api`), an id of the top model is not. A substance id must be unique over all models (a metabolite takes the prefix of its parent, `api_m2`), must not contain `__` and must not start with a token of the scheme (`pl_`, `bc_`, `free_`, `z1_`, ...).

<div class="pbpk-ids" markdown>

| id | what | made by |
|---|---|---|
| `<O>__Cpl_<s>` | a local id of an organ in the flat model | `ids.flat(o, ids.species(part, s))` |
| `Cpl_<s>` | the plasma concentration in an organ (local) | `ids.species(part, s)` |
| `Cti_<s>` | the merged tissue (interstitial and intracellular) of an organ | `ids.species(Part.TISSUE, s)` |
| `Cu_ti_<s>` | the unbound tissue concentration of an organ | `ids.unbound(Part.TISSUE, s)` |
| `C<pool>_<s>` | the plasma concentration of a blood pool (`ar`, `ve`, `po`) | `ids.pool_species(pool, s)` |
| `C<pool>_bc_<s>` | the blood cells of a pool (kinetic erythrocytes) | `ids.pool_species(pool, s, Part.BLOOD_CELLS)` |
| `C<pool>_free_<s>` | the unbound plasma concentration of a pool | `ids.free(pool, s)` |
| `C<pool>_blood_<s>` | the blood concentration of a pool | `ids.blood(pool, s)` |
| `Carm_<s>` | the arm vein plasma concentration, the default fit target | `ids.arm(s)` |
| `Carm_free_<s>` | the unbound arm vein concentration | `ids.arm_free(s)` |
| `Chv_<s>` | the hepatic vein plasma concentration | `ids.hepatic_vein(s)` |
| `Cout_<o>_<s>` | the plasma leaving an organ | `ids.organ_outflow(o, s)` |
| `C<o>_<part>_<s>` | a species of an organ which a PD component or an observer reads (`Cli_ic_<s>`), an output of the organ | `ids.organ_species(o, ids.species(part, s))` |
| `Cu_<o>_<part>_<s>` | an unbound concentration of an organ part which a PD component or an observer reads (`Cu_li_ic_<s>`), an output of the organ | `ids.organ_species(o, ids.unbound(part, s))` |
| `Aurine_<s>` | the amount in urine | `ids.urine(s)` |
| `Rurine_<s>` | the urinary excretion rate (mmole/min), an output of the kidney | `ids.urine_rate(s)` |
| `Afeces_<s>` | the amount in feces | `ids.feces(s)` |
| `Alumen_<s>` | the dissolved amount in the intestinal lumen, the sum over the segments (mmole), a readout of the intestine | `ids.lumen_amount(s)` |
| `Rabs_<s>` | the absorption rate, the fluxes from the lumen into the intestinal wall (mmole/min), a readout of the intestine | `ids.absorption_rate(s)` |
| `Reff_<s>` | the efflux rate, the fluxes from the intestinal wall into the lumen (mmole/min), a readout of the intestine | `ids.efflux_rate(s)` |
| `Rfeces_<s>` | the fecal excretion rate, the net flux of the intestine into the feces (mmole/min), a readout of the intestine | `ids.feces_rate(s)` |
| `Cbile_<s>` | the concentration in bile | `ids.bile(s)` |
| `Mr_<s>` | the molecular weight | `ids.mr(s)` |
| `fu_ref_<s>` | the unbound fraction of the healthy individual (input) | `ids.fu_ref(s)` |
| `fu_<s>` | the unbound fraction, `fu_ref_<s>` changed by the impairment | `ids.fu(s)` |
| `Kp_<s>` | the tissue to plasma coefficient of the substance | `ids.kp(s)` |
| `Kp_<o>_<s>` | the coefficient of one organ (an override, a Kp method or the interstitial share) | `ids.kp_organ(o, s)` |
| `fKp_<s>` | the factor of the coefficients of a Kp method or of the interstitial shares | `ids.fkp(s)` |
| `ftissue_<s>` | the whole-body exchange between plasma and tissue | `ids.ftissue(s)` |
| `PS_<s>` | the exchange of an organ (l/min, local id) | `ids.ps(s)` |
| `kps_<s>` | the exchange of an organ per litre of organ (`ps`, 1/min, local id) | `ids.kps(s)` |
| `Cslow_<s>` | the slow (deep) pool of the tissue of an organ (`SlowPool`, local id) | `ids.slow_species(s)` |
| `kslow_<s>` | the exchange of the slow pool per litre of tissue (1/min) | `ids.kslow(s)` |
| `Kpslow_<s>` | the slow pool to tissue coefficient at equilibrium | `ids.kp_slow(s)` |
| `PS_ref_<s>` | the exchange of the whole healthy liver (`ps` with `ps_per_volume=False`, l/min, local id) | `ids.ps_ref(s)` |
| `BP_<s>` | the blood to plasma ratio | `ids.bp(s)` |
| `Kbc_<s>` | the blood cell to plasma coefficient | `ids.kbc(s)` |
| `f_reabsorbed_<s>` | the constant fraction of the filtrate which is reabsorbed | `ids.fraction_reabsorbed(s)` |
| `reab_threshold_<s>` | the renal threshold of a saturable reabsorption, a concentration of the filtrate (mM) | `ids.reabsorption_threshold(s)` |
| `reab_n_<s>` | the sharpness of the threshold of a saturable reabsorption | `ids.reabsorption_n(s)` |
| `baseline_<s>` | the plasma concentration of an endogenous substance at the start (mM) | `ids.baseline(s)` |
| `IVDOSE_<s>` | the intravenous dose state (mg) | `ids.iv_dose(s)` |
| `Ri_<s>` | the intravenous infusion (mg/min) | `ids.infusion_rate(s)` |
| `PODOSE_<s>` | the oral dose state of a tablet (mg) | `ids.po_dose(s)` |
| `SOLDOSE_<s>` | the oral dose state of a solution (mg) | `ids.sol_dose(s)` |
| `Ka_dis_<s>` | the first order dissolution of a tablet (1/hr) | `ids.dissolution(s)` |
| `ka_<s>` | the absorption rate constant of the small intestine (1/hr) | `ids.absorption(s)` |
| `Clumen_<seg>_<s>` | the dissolved substance in a lumen segment | `ids.lumen(s, seg)` |
| `Cmu_<seg>_<s>` | the mucosa (enterocytes) of a segment | `ids.mucosa(s, seg)` |
| `V<o>` | the volume of an organ (tissue and blood) | `ids.organ_volume(o)` |
| `Q<o>` | the blood flow of an organ (the hepatic artery for the liver) | `ids.organ_flow(o)` |
| `M<o>` | the tissue mass of an organ without blood | `ids.organ_mass(o)` |
| `f_vas_<o>` | the vascular fraction of the organ volume | `ids.vascular_fraction(o)` |
| `f_<protein>` | the activity of a protein in an organ (`LI__f_cyp3a4`) | `ids.activity(protein)` |
| `E_<protein>` | the turnover state of a protein in an organ | `ids.enzyme(protein)` |
| `Ki_<protein>_<s>` | the inhibition constant of a substance on a protein | `ids.ki(protein, s)` |
| `fimp_<o>_<protein>` | the factor of the impairment on a protein of an organ | `ids.impairment_activity(protein, o)` |
| `BW` | the body weight (kg) | `ids.BW` |
| `HEIGHT` | the height (cm) | `ids.HEIGHT` |
| `HCT` | the hematocrit | `ids.HCT` |
| `CO` | the cardiac output | `ids.CO` |
| `GFR` | the glomerular filtration rate | `ids.GFR` |
| `GFR_base` | the glomerular filtration rate at renal function 1 | `ids.GFR_BASE` |
| `f_renal_function` | the renal function, the GFR relative to healthy | `ids.F_RENAL_FUNCTION` |
| `f_cardiac_function` | the cardiac function, a factor on the cardiac output | `ids.F_CARDIAC_FUNCTION` |
| `CP_severity` | the severity of liver cirrhosis (0, Child-Pugh A 1, B 2, C 3) | `ids.CP_SEVERITY` |
| `HF_severity` | the severity of heart failure (0, NYHA I to IV 1 to 4) | `ids.HF_SEVERITY` |
| `f_shunt` | the portosystemic shunt of the liver | `ids.F_SHUNT` |
| `f_functional` | the functional fraction of the liver tissue | `ids.F_FUNCTIONAL` |
| `f_tubular` | the share of the tubular capacity of the kidney which follows the renal function (`KI__f_tubular`) | `ids.F_TUBULAR` |

</div>

Every other id of the scheme, the ids a model rarely reads: the physiology, the intestine and the nephron, the dissolution and the doses, the protein turnover, the starts and deltas of the observers and the reactions and local ids of the organs (`<zone>` a zone token `z1`, `<tseg>` a tubule segment `pt`, `<kind>` a dose `iv`, `po` or `sol`):

<div class="pbpk-ids" markdown>

| id | what | made by |
|---|---|---|
| `AGE` | the age (yr), a covariate | `ids.AGE` |
| `BMI` | the body mass index | `ids.BMI` |
| `BSA` | the body surface area (m2) | `ids.BSA` |
| `FFM` | the fat-free mass (kg) | `ids.FFM` |
| `Mlean` | the mass of the lean organ tissues without their blood (kg) | `ids.MLEAN` |
| `SV_tis` | the scaling factor of the lean tissues, the blood and the cardiac output | `ids.SV_TIS` |
| `Vblood` | the total blood volume (l) | `ids.VBLOOD` |
| `CO_base` | the cardiac output without impairment, of which the organ flows are fractions | `ids.CO_BASE` |
| `Qh` | the total liver blood flow | `ids.QH` |
| `Qha` | the blood flow of the hepatic artery | `ids.QHA` |
| `Qpo` | the blood flow of the portal vein | `ids.QPO` |
| `PS_vas` | the exchange of plasma and interstitium of an organ (local id) | `ids.PS_VAS` |
| `eGFR` | the GFR per 1.73 m2 body surface area of the renal function, which gives the CKD stage | `ids.EGFR` |
| `CKD_hepatic` | the switch (0 off, 1 on) of the hepatic changes of chronic kidney disease | `ids.CKD_HEPATIC` |
| `FED` | the prandial state, 0 fasted (default) and 1 fed, which slows the gastric emptying | `ids.FED` |
| `f_zone` | the volume fraction of a zone of a zonated organ, 1/n of n zones (local id) | `ids.F_ZONE` |
| `GU` | the submodel of the intestine (the small intestine, or the intestine which replaces small and large intestine) | `ids.GUT` |
| `Qgu` | the blood flow of small and large intestine (local id) | `ids.GUT_FLOW` |
| `Qgu_wall` | the blood flow of the intestinal wall past the mucosa (local id) | `ids.GUT_WALL_FLOW` |
| `Vgu` | the volume of the intestine which replaces small and large intestine (local id) | `ids.GUT_VOLUME` |
| `f_vas_gu` | the vascular fraction of the intestine (local id) | `ids.GUT_VASCULAR_FRACTION` |
| `f_int_gu` | the interstitial fraction of the intestine (local id) | `ids.GUT_INTERSTITIAL_FRACTION` |
| `f_wall` | the fraction of the intestinal tissue without the mucosa (local id) | `ids.WALL_FRACTION` |
| `f_mucosa_flow` | the fraction of the intestinal blood flow which perfuses the mucosa (local id) | `ids.F_MUCOSA_FLOW` |
| `f_mucosa_volume` | the fraction of the intestinal tissue which is mucosa (local id) | `ids.F_MUCOSA_VOLUME` |
| `Qmu_total` | the blood flow of the mucosa of all segments (local id) | `ids.MUCOSA_FLOW_TOTAL` |
| `Vmu_total` | the volume of the mucosa of all segments (local id) | `ids.MUCOSA_VOLUME_TOTAL` |
| `k_bile` | the rate of the transit of the bile into the intestinal lumen | `ids.K_BILE` |
| `t_lumen` | the transit time of the single intestinal lumen into the feces (local id) | `ids.T_LUMEN` |
| `T_st_fed` | the transit time of the stomach in the fed state (local id) | `ids.T_STOMACH_FED` |
| `Vbile` | the volume of the bile, the compartment of `Cbile_<s>` | `ids.VBILE` |
| `Vfeces` | the compartment of the feces amounts `Afeces_<s>` | `ids.VFECES` |
| `Vlumen` | the volume of the single intestinal lumen, the compartment of `Clumen_<s>` (local id) | `ids.VLUMEN` |
| `Vurine` | the compartment of the urine amounts `Aurine_<s>` | `ids.VURINE` |
| `<O>` | the submodel of an organ | `ids.submodel(o)` |
| `V<part>` | a compartment of an organ (local id) | `ids.compartment(part)` |
| `V<part>_<zone>` | a compartment of a zone of a zonated organ (local id) | `ids.compartment(part, zone)` |
| `z1` | the token of the first zone of a zonated organ | `ids.zone(1)` |
| `Vti_wall` | a tissue compartment of the intestinal wall without the mucosa (local id) | `ids.wall_compartment(Part.TISSUE)` |
| `V<pool>_<part>` | a compartment of a blood pool | `ids.pool_compartment(pool, part)` |
| `V<pool>` | the volume of a blood pool | `ids.pool_volume(pool)` |
| `fV<pool>` | the volume of a blood pool as fraction of the blood volume | `ids.pool_fraction(pool)` |
| `V<pool>_pl_<s>` | the compartment of a pool species with erythrocytes in equilibrium | `ids.pool_capacity(pool, s)` |
| `Vpl_<s>` | the compartment of the plasma species of an organ with erythrocytes in equilibrium (local id) | `ids.plasma_capacity(s)` |
| `fQ<o>` | the blood flow of an organ as fraction of the cardiac output | `ids.flow_fraction(o)` |
| `f_int_<o>` | the interstitial fraction of the organ volume | `ids.interstitial_fraction(o)` |
| `GFR_ref` | the reference value of a parameter (`GFR_ref`, `Mli_ref`) | `ids.ref(ids.GFR)` |
| `fimp_GFR` | the factor of all impairments on a quantity | `ids.impairment_factor(ids.GFR)` |
| `fcp_albumin` | the factor of one condition (Child-Pugh, CKD) on a quantity | `ids.condition_factor('cp', 'albumin')` |
| `Vti_total_ref_<s>` | the reference tissue volume of the organs the substance enters | `ids.tissue_volume_total_ref(s)` |
| `K_int_pls_<s>` | the interstitial to plasma coefficient (local id) | `ids.k_int_pls(s)` |
| `K_cell_pls_<s>` | the intracellular to plasma coefficient (local id) | `ids.k_cell_pls(s)` |
| `K_or_<s>` | the organ to plasma coefficient of the flow-limited organ (local id) | `ids.k_organ(s)` |
| `kex_bc_<s>` | the exchange rate between plasma and blood cells (kinetic erythrocytes) | `ids.kbc_exchange(s)` |
| `Carm_blood_<s>` | the blood concentration of the arm vein | `ids.arm_blood(s)` |
| `w_arm_sk` | the weight of an organ outflow (or a pool) in the arm vein (`w_arm_sk`) | `ids.arm_weight(Organ.SKIN)` |
| `Chv_free_<s>` | the unbound concentration of the hepatic vein | `ids.hepatic_vein_free(s)` |
| `Chv_blood_<s>` | the blood concentration of the hepatic vein | `ids.hepatic_vein_blood(s)` |
| `Cout_<s>` | the concentration leaving an organ (local id) | `ids.outflow_conc(s)` |
| `Cout_bc_<s>` | the blood cell concentration leaving an organ (local id) | `ids.outflow_conc_bc(s)` |
| `Cout_<o>_bc_<s>` | the blood cell concentration leaving an organ | `ids.organ_outflow_bc(o, s)` |
| `Kiv_<s>` | the rate of an intravenous injection | `ids.injection_rate(s)` |
| `ti_<s>` | the half-life of an intravenous injection (s) | `ids.injection_time(s)` |
| `Ka_sol_<s>` | the rate of the application of a solution | `ids.solution(s)` |
| `fabs_<s>` | the fraction of the absorbed flux which enters the body | `ids.fraction_absorbed(s)` |
| `f_ka_colon_<s>` | the factor of the absorption in the large intestine (segmented intestine) | `ids.colon_absorption_factor(s)` |
| `f_ka_<seg>_<s>` | the factor of the absorption in a lumen segment, the regional absorption profile (segmented intestine) | `ids.regional_absorption_factor(s, seg)` |
| `kconv_<s>_<m>` | the rate constant (1/hr) of the conversion of `<s>` into the substance `<m>` in the lumen (segmented intestine) | `ids.lumen_conversion_rate(s, m)` |
| `J_conv_<seg>_<s>_<m>` | the conversion of `<s>` into the substance `<m>` in a lumen segment (local id) | `ids.lumen_conversion_reaction(s, m, seg)` |
| `J_dis_<seg>_<s>` | the dissolution in a lumen segment (local id) | `ids.dissolution_reaction(s, seg)` |
| `kdis_<s>` | the dissolution rate of a tablet (1/min, local id) | `ids.dissolution_rate(s)` |
| `f_dis_fed_<s>` | the factor of the dissolution rate in the fed state | `ids.dissolution_fed_factor(s)` |
| `f_dis_release_<s>` | the factor of the dissolution rate of a modified-release formulation | `ids.dissolution_release_factor(s)` |
| `lag_dis_<s>` | the lag time of a Weibull dissolution (min) | `ids.dissolution_lag(s)` |
| `t50_dis_<s>` | the time of half dissolution of a Weibull dissolution (min) | `ids.dissolution_t50(s)` |
| `shape_dis_<s>` | the shape of a Weibull dissolution | `ids.dissolution_shape(s)` |
| `D50_dis_<s>` | the undissolved amount of half the dissolution rate (mg) | `ids.dissolution_dose50(s)` |
| `age_dis_<s>` | the mean age of the undissolved dose (min, local id) | `ids.dissolution_age(s)` |
| `W_dis_<s>` | the age-weighted undissolved amount (mmole*min, local id) | `ids.dissolution_clock(s)` |
| `Asolid_<seg>_<s>` | the undissolved amount in a lumen segment (local id) | `ids.solid(s, seg)` |
| `Vlumen_<seg>` | the volume of a lumen segment (local id) | `ids.lumen_volume(seg)` |
| `T_<seg>` | the transit time of a lumen segment (local id) | `ids.transit_time(seg)` |
| `f_area_<seg>` | the absorption scaling of a lumen segment (local id) | `ids.area_factor(seg)` |
| `Vmu_<seg>` | the volume of the mucosa of a segment (local id) | `ids.mucosa_volume(seg)` |
| `Qmu_<seg>` | the blood flow of the mucosa of a segment (local id) | `ids.mucosa_flow(seg)` |
| `K_mu_<s>` | the mucosa to plasma coefficient (local id) | `ids.mucosa_coefficient(s)` |
| `Cu_mu_<seg>_<s>` | the unbound concentration in the mucosa of a segment (local id) | `ids.mucosa_unbound(s, seg)` |
| `Ctub_<tseg>_<s>` | the concentration in the lumen of a tubule segment of the nephron kidney (local id) | `ids.tubule(s, tseg)` |
| `Vtub_<tseg>` | the volume of the lumen of a tubule segment (local id) | `ids.tubule_volume(tseg)` |
| `f_tub_<tseg>` | the tubular flow leaving a segment as fraction of the GFR (local id) | `ids.tubule_flow_fraction(tseg)` |
| `kdeg_<protein>` | the degradation rate of a protein in an organ (turnover) | `ids.kdeg(protein)` |
| `Kiu_<protein>_<s>` | the uncompetitive constant of a mixed inhibition | `ids.kiu(protein, s)` |
| `kinact_<protein>_<s>` | the maximal inactivation rate of a mechanism-based inhibition | `ids.kinact(protein, s)` |
| `Emax_<protein>_<s>` | the maximal induction of a protein | `ids.emax(protein, s)` |
| `EC50_<protein>_<s>` | the half-maximal inducing concentration | `ids.ec50(protein, s)` |
| `p_cyp3a4_<zone>` | the relative abundance of a profile in a zone (local id) | `ids.placement_profile('cyp3a4', zone)` |
| `w_cyp3a4_<seg>` | the share of a process in a segment or a zone (local id) | `ids.placement_weight('cyp3a4', seg)` |
| `API2M_<seg>` | a process of the chemistry in a segment or a zone (local id) | `ids.segment_process('API2M', seg)` |
| `<p>_Vmax_ratio` | the ratio of a constant of a process to the same constant of another process of the organ, `RateOf` (local id) | `ids.rate_ratio('API2M', 'Vmax')` |
| `f_<p>_<s>` | the fraction of a metabolism which forms a product, a value or the remainder `REST` (local id) | `ids.product_fraction('API2M', s)` |
| `<p>_to_<s>` | the formation of a product of a metabolism with product fractions (local id) | `ids.product_reaction('API2M', s)` |
| `Asi_solid_<seg>_<s>` | an amount species of an organ which a PD component or an observer reads (`A<x>` as `A<o>_<x>`), an output of the organ | `ids.organ_species(Organ.SMALL_INTESTINE, ids.solid(s, seg))` |
| `Carm_api_0` | the value of an observer at the start (`PBPKDefinition.deltas`) | `ids.start_value('Carm_api')` |
| `Carm_api_delta` | the change of an observer over its start (`PBPKDefinition.deltas`) | `ids.delta('Carm_api')` |
| `J_<kind>_<s>` | the application of a dose state (`iv`, `po`, `sol`) | `ids.dose_reaction(kind, s)` |
| `J_ti_<s>` | the passive exchange into a part of an organ (local id) | `ids.exchange_reaction(Part.TISSUE, s)` |
| `J_slow_<s>` | the exchange between the tissue and its slow pool (local id) | `ids.slow_reaction(s)` |
| `J_in_ar_<part>_<s>` | the inflow from a pool into an organ (local id) | `ids.inflow_reaction(Pool.ARTERIAL, part, s)` |
| `J_out_<part>_<s>` | the outflow from an organ into its pool (local id) | `ids.outflow_reaction(part, s)` |
| `J_shunt_po_<part>_<s>` | the shunt from a pool past the liver (local id) | `ids.shunt_reaction(Pool.PORTAL, part, s)` |
| `J_sin_<zone>_<part>_<s>` | the flow from a zone into the next one (local id) | `ids.sinusoid_reaction(part, s, zone)` |
| `J_bc_<pool>_<s>` | the exchange between plasma and blood cells of a pool | `ids.pool_exchange_reaction(pool, s)` |
| `J_gfr_<s>` | the glomerular filtration into the urine (local id) | `ids.filtration_reaction(s)` |
| `J_bile_<s>` | the transit of the bile into the intestinal lumen | `ids.bile_reaction(s)` |
| `J_abs_<seg>_<s>` | the absorption from a lumen segment (local id) | `ids.absorption_reaction(s, seg)` |
| `J_feces_<seg>_<s>` | the unabsorbed fraction of the absorption into the feces (local id) | `ids.feces_reaction(s, seg)` |
| `J_tr_<seg>_<s>` | the transit out of a lumen segment (local id) | `ids.transit_reaction(s, seg)` |
| `J_mu_in_<seg>_<s>` | a flow of the mucosa of a segment (local id) | `ids.mucosa_reaction('in', s, seg)` |
| `J_tub_<tseg>_<s>` | the tubular flow out of a tubule segment (local id) | `ids.tubule_reaction(s, tseg)` |
| `J_reab_<s>` | the reabsorbed fraction out of the collecting duct (local id) | `ids.tubule_reabsorption_reaction(s)` |

</div>

A process of the chemistry keeps the id its definition gives (`API2M2`), its parameters are `<p>_Vmax`, `<p>_Km`, `<p>_k`, `<p>_n` in the organ (`LI__API2M2_Vmax` in the flat model), the ratio of a shared constant `<p>_Vmax_ratio`, the product fractions `f_<p>_<s>`, and the parameters of a custom reaction keep their ids in the organ (`MU__k_deep`). The tables are checked against `pkpdsim.pbpk.ids` by `tests/docs/test_pbpk_ids.py`: every row is made by its function, and every name of `ids.__all__` has a row (but the tuples of the segment tokens and `species_of`, which lists ids), so a new id of the scheme fails the test until it is documented.
