# Migration

The models of pkdb_models were built from one body model which was copied and edited about 40 times, plus a liver, a kidney and an intestine model per substance. A model built on the library keeps its experiments, its data and its fits; its factory calls `build`, and its body, liver, kidney and intestine modules are deleted. This guide takes the glimepiride model of pkdb_models (`pkdb_models/models/glimepiride/models/`) as the example. The three pilots of the library are hctz, apixaban and glimepiride; the other models follow in the next projects, with what the pilots teach.

## What moves where

| hand-written model | lines | on the library |
|---|---|---|
| `model_body.py`: the whole-body model, organs, blood flows, dosing, the scaling with `BW`, the ports and replacements | 1472 | the physiology, the standard organs, the routes, the sinks and the observers of the library; nothing in the package |
| `model_liver.py`: `GLIIM` (uptake of glimepiride), `GLI2M1` (CYP2C9), `M1EX` and `M2EX` (export of the metabolites), `M12M2`, `M1BIL` and `M2BIL` (biliary excretion) | 374 | the distribution of `gli` (Kp, the exchange of the liver `ps`), `Metabolism("GLI2M1", ...)` and `Metabolism("M12M2", ...)`; the products are placed into the plasma of the liver (no export reaction); a biliary excretion is `Transport(..., Direction.EFFLUX, Membrane.APICAL, ...)` of the liver |
| `model_kidney.py`: `M1EX`, `M2EX` (`GFR*M1EX_f*m1_ext`, a fraction of the GFR) | 273 | `Excretion(urine=True)` of the metabolites: `fu*GFR*(1 - f_reabsorbed)`, with `f_reabsorbed = 1 - M1EX_f` for a fraction below 1; a fraction above 1 is a secretion, an apical efflux of the kidney, which does not follow `f_renal_function` yet (see [Known gaps](#known-gaps)) |
| `model_intestine.py`: dissolution of the tablet, `GLIABS` with the fraction `F_gli_abs`, the excretion of the lumen into the feces | 350 | `Route.oral("gli")` and the segmented intestine (dissolution, transit, absorption `ka_gli`, feces); or the single lumen (`fabs_gli`) |
| `factory.py`: create the four models, flatten, the ODEs, the COMBINE archive | 145 | `build(DEFINITION, output_dir)` for the comp and the flat model; the ODEs, the license texts and the COMBINE archive stay in the factory of the package |
| `annotations.py`, `templates.py` | | the annotation of the library (see [Annotation](annotation.md)); the substance annotations go into `Substance.annotations` |

The definition of glimepiride is the one of the [tutorial](tutorial.md) with the fitted values of the adopted fits, `fitted(FITS.parameter(...))`.

## Ids

The library keeps the ids of the hand-written models where they fit its id scheme, so experiments change as little as possible. The flat ids of glimepiride:

| hand-written model | on the library | note |
|---|---|---|
| `Cve_gli` (venous plasma) | `Cve_gli` | the same; plasma data map to `Carm_gli` (below) |
| `Car_gli`, `Cpo_gli` | `Car_gli`, `Cpo_gli` | the same |
| `Chv_gli` (a species) | `Chv_gli` (an observer) | the outflow of the liver mixed with the shunt |
| `m1`, `m2` (substance ids) | `gli_m1`, `gli_m2` | a substance id is unique over all models: `Cve_m1` becomes `Cve_gli_m1`, `Aurine_m1` becomes `Aurine_gli_m1` |
| `Aurine_m1`, `Afeces_gli` | `Aurine_gli_m1`, `Afeces_gli` | the sinks keep their pattern |
| `Cve_m1_m2` (the sum of the metabolites) | none | the library has no observers of the definition yet (see [Known gaps](#known-gaps)) |
| `Cli_plasma_gli`, `LI__gli` | `LI__Cpl_gli`, `LI__Cti_gli` | the plasma and the tissue of the liver (local ids of the organ) |
| `Clu_gli`, `Cre_gli` ("rest") | `LU__Cti_gli`, the organs `HE__`, `BR__`, `MU__`, ... | every organ of ICRP 89, no rest compartment |
| `GU__gli_stomach`, `GU__gli_lumen` | `PODOSE_gli`, `GU__Clumen_<seg>_gli`, `GU__Asolid_<seg>_gli` | the undissolved tablet, the lumen segments |
| `PODOSE_gli`, `IVDOSE_gli`, `ti_gli`, `Ri_gli` | the same | the dose states and the injection |
| `Ki_gli` (injection rate) | `Kiv_gli` | `Ki_<protein>_<s>` is an inhibition constant now |
| `GU__Ka_dis_gli` | `Ka_dis_gli` | a parameter of the top model (1/hr) |
| `GU__GLIABS_k` (1/min), `GU__F_gli_abs` | `ka_gli` (1/hr), `fabs_gli` (single lumen) | the absorption rate constant of the small intestine |
| `LI__GLI2M1_Vmax` | `LI__GLI2M1_Vmax` | the same when the process keeps its id: a Vmax per litre of liver in both |
| `LI__M12M2_k` (1/min, per litre of liver: `M12M2_k*Vli*m1`) | `LI__M12M2_k` (l/min, the whole liver) | the same id, another meaning: `FirstOrder.k` is a clearance of the whole organ on the unbound concentration, so an old value maps as `k_new = k_old*Vli` (about 2.3 l for the male reference) and does not scale with the liver of another individual (see [Known gaps](#known-gaps)) |
| `LI__GLI2M1_Km_gli` | `LI__GLI2M1_Km` | the Km of a process is `<p>_Km` |
| `LI__fu_gli` | `fu_ref_gli` (input), `fu_gli` | fu is a parameter of the substance, fitted as `fu_ref_gli` |
| `LI__f_cyp2c9` | `LI__f_cyp2c9` | the same: the activity of the protein in the liver |
| `Kp_gli`, `ftissue_gli` | `Kp_gli`, `ftissue_gli` | the same ids, other meaning of `ftissue` (below) |
| `BW`, `HEIGHT` | `BW`, `HEIGHT` | the same covariates |
| `KI__f_renal_function` | `f_renal_function` | a parameter of the top model; its GFR per 1.73 m2 gives the CKD stage |
| `f_cirrhosis` | `CP_severity` (or `f_shunt`, `f_functional`) | the Child-Pugh severity with the vectors of [Impairment](impairment.md); `f_cirrhosis` keeps its meaning in the hand-written models |
| `f_cardiac_function` | `f_cardiac_function` | the same; heart failure by class is `HF_severity` |
| `COBW`, `HR`, `FQ<o>`, `FV<o>`, `Fblood`, `HCT` 0.51 | `CO_ref`, `fQ<o>`, `M<o>`, `HCT` | the reference values of ICRP 89 with their sources |

## The arm vein observer

The hand-written models compared plasma data with the venous plasma `Cve_<s>`. Plasma is sampled from an arm vein, so plasma and blood data map to the arm vein observer `Carm_<s>` by default (decided by Matthias König): the flow-weighted venous outflow of skin and muscle, which differs from the mixed venous blood while the tissues of the arm take up or release the drug (Levitt 2004, Musther 2015). The fit mappings carry the site in their metadata, `PKPDMappingMetaData(sampling_site=SamplingSite.ARM_VEIN)` (the default); a study which states another site (venous catheter, arterial line, arterialised hand vein, hepatic vein) maps to that observer. See [Observers](observers.md).

## Why a refit is needed

The migration does not reproduce the simulations of the hand-written model, and is not meant to: many hand-written models have structural deficits which the library corrects, and every migrated model is refitted (decided by Matthias König). What changes the numbers:

- Blood and plasma: every vascular space has plasma and blood cells with the hematocrit of the sex (0.43 male), the blood flows carry plasma; the hand-written models had five conventions (whole-blood pools with plasma organ volumes, blood flows multiplied by plasma concentrations, hematocrits 0.51, 0.45, 0.40 or none).
- Tissue distribution: one convention, the whole-tissue Kp with the unbound plasma concentration, `PS*fu*(C_p - C_t/Kp)`. `ftissue_<s>` multiplies the unbound concentrations, so its fitted values change, and an organ has the share of its tissue volume of the reference individual. For a highly bound substance the exchange of the liver needs its own `ps` (glimepiride: fu 0.005).
- Physiology: the organs, volumes and flows of ICRP 89 with their sources, and the scaling of Huisinga 2012 (fat-free mass, adipose as the remainder) instead of `ftissue*BW/BW_ref`.
- First order processes: `FirstOrder.k` is a clearance of the whole organ, a per-litre constant of the hand-written model maps as `k*V` of the reference organ.
- Absorption: the segmented intestine has a transit through the segments into the feces, and the single lumen a physiological transit, `Fa = fabs*ka/(ka + 1/T)`; the fitted absorption constants of the hand-written models (without transit) do not carry over.
- Observables: the arm vein instead of the venous plasma.
- Impairment: Child-Pugh, CKD and heart failure as vectors with sources instead of one map `f_cirrhosis` (shunt = tissue loss), which over-predicted the exposure of low-extraction drugs.

The bootstrap of a pilot builds the model with the literature values and the current fitted values where the parameter still exists, then fits and adopts the fits (`run_fit <package> --adopt <fit id>`), so the factory takes the new values with their uncertainties.

## Known gaps

The library as it stands lacks what some hand-written models use. These gaps are known and are planned as the prerequisites of the pilots (plan 2, tasks P1 and P2 on this branch); the pages are updated when they land:

- A first order process (`FirstOrder.k`) is a clearance of the whole organ: it does not scale with the organ volume or the size of the individual, as the per-litre constants of the hand-written models did (`k*V`).
- A tubular secretion (an apical efflux of the kidney) does not follow the renal function `f_renal_function`, as the fraction of the GFR of the hand-written models did (`GFR*M1EX_f`); only the filtration does.
- There is no urinary excretion rate `Rurine_<s>`; the urine is read as the cumulative amount `Aurine_<s>`.
- A definition cannot declare observers of its own (a sum of metabolites such as `Cve_m1_m2`, a ratio).
- `build` writes no creators into the model; the factory of the package adds the model history it needs.

## The pilots are the validation gate

hctz, apixaban and glimepiride are rebuilt on the library in pkdb_models: their factories call `build`, their old body, liver, kidney and intestine modules are deleted, their experiments and fit mappings use the new observables, the models are refitted and the fits adopted, and the fast and slow tests of pkdb_models pass for the three models; `pkpdsim-verify` documents the change of the simulations. The pilots are the validation gate of the library (decided by Matthias König):

- the scientific tests `tests/models/test_<model>.py` keep their literature values and tolerances (clearance, half-life, bioavailability, excretion, impairment, genotypes); only ids and observables change;
- a refitted pilot must pass them and must not fit the data clearly worse than its current adopted fit, compared by the cost of the same fit mappings and noted in the pull request;
- where a pilot fails, the library is changed until it passes. During development the pilots use the branch of pkpdsim (`uv run --with-editable ../pkpdsim ...`); the maintainer releases pkpdsim before the pilot pull requests merge;
- the PD of the pilots (the fluid and electrolytes of hctz, the coagulation of apixaban) stays in the package as PD components with declared ports (`PDComponent`), coupled without `merge_models`.

The pilots are not released into their public model repositories in this project; the static model figure of their COMBINE archive no longer matches and is recorded in their `TODO.md`.

## Steps of a migration

1. Write the substances: short ids unique over all models (metabolites with the prefix of the parent), the molecular weights, the binding, the erythrocytes, the distribution and the excretion; the annotations of the old `annotations.py` into `Substance.annotations`.
2. Write the chemistry of every organ with the process ids of the hand-written model, so that the fitted parameters keep their flat ids (`LI__GLI2M1_Vmax`).
3. Write the definition with the routes, the variants (the zonated liver of a zonated model, the single lumen if the model keeps a first order absorption) and the PD components.
4. Replace the factory by `build` and delete the body, liver, kidney and intestine modules.
5. Change the experiments: the observables (`Carm_<s>` for plasma data, the metabolite ids with the prefix of the parent), the impairment parameters (`CP_severity`, `f_renal_function`, the presets `CHILD_PUGH`, `RENAL_FUNCTION`, `HEART_FAILURE` of `pkpdsim.impairment`), the sampling site of the fit mappings.
6. Bootstrap with the current fitted values where the parameters still exist, refit, adopt, and run the scientific tests of the model.
