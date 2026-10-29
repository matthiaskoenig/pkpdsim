# Migration

The models of pkdb_models were built from one body model which was copied and edited about 40 times, plus a liver, a kidney and an intestine model per substance. A model built on the library keeps its experiments, its data and its fits; its package writes one `PBPKDefinition` (`models/definition.py`), its factory is one call of `create_archive`, and its body, liver, kidney and intestine modules are deleted. This guide takes the glimepiride model of pkdb_models (`pkdb_models/models/glimepiride/models/`) as the example. Seven models are built on the library: the pilots hctz, apixaban and glimepiride and the first group torasemide, inulin, gdeobdtpa and creatinine; the other models follow with the procedure of this page and what the seven taught (the last section).

## What moves where

| hand-written model | lines | on the library |
|---|---|---|
| `model_body.py`: the whole-body model, organs, blood flows, dosing, the scaling with `BW`, the ports and replacements | 1472 | the physiology, the standard organs, the routes, the sinks and the observers of the library; nothing in the package |
| `model_liver.py`: `GLIIM` (uptake of glimepiride), `GLI2M1` (CYP2C9), `M1EX` and `M2EX` (export of the metabolites), `M12M2`, `M1BIL` and `M2BIL` (biliary excretion) | 374 | the distribution of `gli` (Kp, the exchange of the liver `ps`), `Metabolism("GLI2M1", ...)` and `Metabolism("M12M2", ...)`; the products are placed into the plasma of the liver (no export reaction); a biliary excretion is `Transport(..., Direction.EFFLUX, Membrane.APICAL, ...)` of the liver |
| `model_kidney.py`: `M1EX`, `M2EX` (`GFR*M1EX_f*m1_ext`, a fraction of the GFR) | 273 | `Excretion(urine=True)` of the metabolites: `fu*GFR*(1 - f_reabsorbed)`, with `f_reabsorbed = 1 - M1EX_f` for a fraction below 1; a fraction above 1 is a secretion, an apical efflux of the kidney, which follows `f_renal_function` with the tubular capacity of the kidney, as the fraction of the GFR did |
| `model_intestine.py`: dissolution of the tablet, `GLIABS` with the fraction `F_gli_abs`, the excretion of the lumen into the feces | 350 | `Route.oral("gli")` and the segmented intestine (dissolution, transit, absorption `ka_gli`, feces); or the single lumen (`fabs_gli`) |
| `factory.py`: create the four models, flatten, the ODEs, the COMBINE archive | 145 | `create_archive(DEFINITION, model_output_dir, model=..., copyright=..., readme=...)`: the comp and the flat model, the ODEs and the antimony of the flat model, the license texts and the COMBINE archive; a standalone model of the package is an `extra_models` entry, its ODEs an `extra_files` entry |
| `annotations.py`, `templates.py` | | the annotation of the library (see [Annotation](annotation.md)); the substance annotations go into `Substance.annotations`, the creators into `PBPKDefinition.creators` (the model history of the comp and the flat model) |
| the helpers of the experiments, the renal and hepatic maps | | `PBPKSimulationExperiment` of `pkpdsim.pbpk` (below) |

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
| the rate of a urinary excretion reaction (`KI__M1EX`, hctz `KI__HCTZEX`) | `Rurine_gli_m1` | the net urinary excretion rate: the fluxes of the kidney into the urine less a custom reaction from the urine (a reabsorption) |
| `Cve_m1_m2` (the sum of the metabolites) | `Cve_gli_m1_m2` | an observer of the definition, `Observer("Cve_gli_m1_m2", "Cve_gli_m1 + Cve_gli_m2", U.mM, ...)` (see [Observers](observers.md#observers-of-a-definition)) |
| `Cli_plasma_gli`, `LI__gli` | `LI__Cpl_gli`, `LI__Cti_gli` | the plasma and the tissue of the liver (local ids of the organ) |
| `Clu_gli`, `Cre_gli` ("rest") | `LU__Cti_gli`, the organs `HE__`, `BR__`, `MU__`, ... | every organ of ICRP 89, no rest compartment |
| `GU__gli_stomach`, `GU__gli_lumen` | `PODOSE_gli`, `GU__Clumen_<seg>_gli`, `GU__Asolid_<seg>_gli` | the undissolved tablet, the lumen segments |
| `PODOSE_gli`, `IVDOSE_gli`, `ti_gli`, `Ri_gli` | the same | the dose states and the injection |
| `Ki_gli` (injection rate) | `Kiv_gli` | `Ki_<protein>_<s>` is an inhibition constant now |
| `GU__Ka_dis_gli` | `Ka_dis_gli` | a parameter of the top model (1/hr); the fed state and a modified release are factors of it (`f_dis_fed_<s>`, `f_dis_release_<s>`) |
| `GU__GLIABS_k` (1/min), `GU__F_gli_abs` | `ka_gli` (1/hr), `fabs_gli` (single lumen) | the absorption rate constant of the small intestine, one value or one per substance (`PerSubstance`) |
| `LI__GLI2M1_Vmax` | `LI__GLI2M1_Vmax` | the same when the process keeps its id: a Vmax per litre of liver in both |
| `LI__M12M2_k` (1/min, per litre of liver: `M12M2_k*Vli*m1`) | `LI__M12M2_k` (1/min, per litre of liver) | the same id and value with `FirstOrder(k)` (per litre of organ, the default): the rate `M12M2_k*Vli*Cu` on the unbound concentration, which follows the liver of the individual (divided by `fu_ref` for a bound substance, whose old rate read the total concentration) |
| `LI__GLI2M1_Km_gli` | `LI__GLI2M1_Km` | the Km of a process is `<p>_Km` |
| `LI__fu_gli` | `fu_ref_gli` (input), `fu_gli` | fu is a parameter of the substance, fitted as `fu_ref_gli` |
| `LI__f_cyp2c9` | `LI__f_cyp2c9` | the same: the activity of the protein in the liver |
| `Kp_gli`, `ftissue_gli` | `Kp_gli`, `ftissue_gli` | the same ids, other meaning of `ftissue` (below) |
| `BW`, `HEIGHT` | `BW`, `HEIGHT` | the same covariates |
| `KI__f_renal_function` | `f_renal_function` | a parameter of the top model; its GFR per 1.73 m2 gives the CKD stage |
| `f_cirrhosis` | `CP_severity` (or `f_shunt`, `f_functional`) | the Child-Pugh severity with the vectors of [Impairment](impairment.md); `f_cirrhosis` keeps its meaning in the hand-written models |
| `f_cardiac_function` | `f_cardiac_function` | the same; heart failure by class is `HF_severity` |
| `COBW`, `HR`, `FQ<o>`, `FV<o>`, `Fblood`, `HCT` 0.51 | `CO_ref`, `fQ<o>`, `M<o>`, `HCT` | the reference values of ICRP 89 with their sources |

The tables of the [overview](index.md#the-id-scheme) list every id of the scheme.

## The arm vein observer and the sampling site

The hand-written models compared plasma data with the venous plasma `Cve_<s>`. Plasma is sampled from an arm vein, so plasma, serum and blood data map to the arm vein observer `Carm_<s>` by default (decided by Matthias König): the flow-weighted venous outflow of skin and muscle, which differs from the mixed venous blood while the tissues of the arm take up or release the drug (Levitt 2004, Musther 2015). The fit mappings carry the site in their metadata: `PKPDMappingMetaData.sampling_site` holds what a study states (`None` if nothing), and `resolved_sampling_site()` gives the site, which follows the tissue without a statement (plasma, serum and blood the arm vein, arterial blood and plasma `ARTERIAL`, urine, feces, bile, saliva and no tissue `NOT_APPLICABLE`; a contradiction is a `ValueError`). A study which states another site (venous catheter, arterial line, arterialised hand vein, hepatic vein) maps to that observer. A model test asserts that `pkpdsim.fit.check_sampling_sites(fit_mappings, SITES)` is empty, with `SITES` the site of every observable of the mappings (`{"Carm_gli": SamplingSite.ARM_VEIN, "Aurine_gli_m1": SamplingSite.NOT_APPLICABLE, ...}`). See [Observers](observers.md).

## Why a refit is needed

The migration does not reproduce the simulations of the hand-written model, and is not meant to: many hand-written models have structural deficits which the library corrects, and every migrated model is refitted (decided by Matthias König). What changes the numbers:

- Blood and plasma: every vascular space has plasma and blood cells with the hematocrit of the sex (0.43 male), the blood flows carry plasma; the hand-written models had five conventions (whole-blood pools with plasma organ volumes, blood flows multiplied by plasma concentrations, hematocrits 0.51, 0.45, 0.40 or none).
- Tissue distribution: one convention, the whole-tissue Kp with the unbound plasma concentration, `PS*fu*(C_p - C_t/Kp)`. `ftissue_<s>` multiplies the unbound concentrations, so its fitted values change (`ftissue_old*Kp_old/fu` where the old exchange read the total concentrations, `ftissue_old*Kp_old` where it carried fu), and an organ has the share of its tissue volume of the reference individual. For a highly bound substance the exchange of the liver needs its own `ps` (glimepiride: fu 0.005).
- Physiology: the organs, volumes and flows of ICRP 89 with their sources, and the scaling of Huisinga 2012 (fat-free mass, adipose as the remainder) instead of `ftissue*BW/BW_ref`.
- First order processes act on the unbound concentration: a per-litre constant of the hand-written model keeps its value with `FirstOrder(k)`, which is per litre of organ (divided by `fu_ref` where the old rate read the total concentration); a clearance of the whole organ (l/min) needs `FirstOrder(k, per_volume=False)`.
- Kidney: the glomerulus filters the afferent (arterial) plasma, `fu*GFR*(1 - f_reabsorbed)*Car`, and the plasma leaving the kidney lacks the filtrate, so the filtration clearance is `fu*GFR*(1 - f_reabsorbed)` and a secretion acts on the plasma flow past it (glomerulus and tubule in series). Every metabolism and transport of the kidney runs with the tubular capacity, which follows `f_renal_function` (the intact nephron hypothesis), as a secretion of a fraction of the GFR did.
- Absorption: the segmented intestine has a transit through the segments into the feces, and the single lumen a physiological transit, `Fa = fabs*ka/(ka + 1/T)`; the fitted absorption constants of the hand-written models (without transit) do not carry over.
- Observables: the arm vein instead of the venous plasma.
- Impairment: Child-Pugh, CKD and heart failure as vectors with sources instead of one map `f_cirrhosis` (shunt = tissue loss), which over-predicted the exposure of low-extraction drugs; the healthy GFR is that of the individual (`GFR_base`, it scales with the body surface area), not 101 ml/min.

## The package on the library

**The definition** (`models/definition.py`) holds the substances, the chemistry, the variants, the routes, the PD components (`PDComponent`, coupled by declared ports; an input may read a species of an organ by its flat id, `LI__Cic_<s>`), the observers the experiments read which the library does not compute (`PBPKDefinition.observers`), the creators, and:

- `parameters`: the constants of the model which no field of the library holds, `ModelParameter(sid, value, unit, name, notes)`, a number, a fitted `Value` or a `Formula`; the formulas, the observers and the custom reactions of the organs read them (a production rate which a baseline and a reaction share). A fitted constant of one custom reaction is a `ModelParameter` of the reaction, a local id of its organ (`MU__k_deep`, see [Custom reactions](chemistry.md#custom-reactions)); a capacity which two transports share is `RateOf` of the other transport ([Shared capacities](chemistry.md#shared-capacities)); fitted fractions of one metabolism are values of its products with the remainder `REST` ([Product fractions](chemistry.md#product-fractions)), not several metabolisms;
- `Substance.baseline`: the plasma concentration of an endogenous substance at the start (mM), every species of the substance at the equilibrium of its exchange with plasma;
- `starts`: the initial concentration of single species by their flat id, for the organs which a production or an elimination moves away from that equilibrium. The library owns the steady state of its equations: `pkpdsim.pbpk.steady_state` gives the formulas (`organ_plasma`, `organ_tissue`, `kidney_plasma` with a first order and a saturable secretion, `Saturable(vmax, km)`, at a given arterial plasma or at a given renal elimination, `kidney_arterial` (with the elimination it reads the kidney plasma once, which keeps the load time of the model), `liver_plasma`, `arm_vein`, ...), so the model writes no equation of the library (see [Substances](substances.md#starts-of-the-species)); a model with a production checks the starts with a stationarity test (every species over a day, several individuals);
- `deltas`: the observers whose change over the start a study reports (`Carm_<s>_delta`), e.g. the increase of an endogenous substance after a dose.

A Kp of an extracellular substance is the interstitial share of every organ, `Kp(interstitial=1.0, ftissue=...)` (no `Kp_<s>`; a fitted scale is `fitted(FITS.parameter("fKp_<s>"))`), and `Kp(overrides=...)` without a value covers every organ by overrides; a value of a hand-written model per kg or per litre of an organ is converted with the reference individual, `reference_individual(Sex.MALE)` (`bw` 73 kg, `volumes[Organ.LIVER]` 2.33 l), never with a copy of its numbers.

**The factory** (`models/factory.py`) is `create_archive(DEFINITION, model_output_dir, model=..., copyright=..., readme=..., extra_models=..., extra_files=...)`, which writes everything the factory of a hand-written model wrote. `models/output/models/` holds exactly what it writes; the old SBML is deleted.

**The experiments** derive from `PBPKSimulationExperiment` of `pkpdsim.pbpk` (`from pkpdsim.pbpk import PBPKSimulationExperiment, REFERENCE_GFR`) and bind the model object. It maps the classes of a study onto the library: `renal_map` (the classes of renal impairment keep their GFR in ml/min, relative to the healthy GFR of the model), `cirrhosis_map` (`CP_severity`), `cardiac_map` (`HF_severity`), and `renal_function(clearance, bw=..., height=...)` gives the renal function of a measured clearance, so the simulated GFR of the individual of the study is the measured one; `healthy_gfr` and `renal_function` return quantities for the changes of a simulation. The covariates and severities (`BW`, `HEIGHT`, `AGE`, `f_renal_function`, `CP_severity`, `HF_severity`, `FED`) go into `preinit_changes`, never into the changes of a simulation (their initial assignments are evaluated again before the start), unless the definition names the id `dynamic`. A model whose runs compile more than they simulate sets `SimulatorSettings(jit_optimization="O0")` (see [Models](../models.md)).

**The sensitivity analysis** names the parameters it varies: `pkpdsim.sensitivity.library_parameters(model_path, ids, fit_ids=..., scaling_ids=..., references=...)` with relative bounds, the severities of `SEVERITIES` from 0, and `exclude_derived` for a table of every constant without the initial assignments of the covariates (see [Scans](../scans.md)).

## The procedure of a migration

1. **Dump the old model.** In a checkout of the old model, build it and dump its experiments (`pkpdsim-verify dump --model pkdb_models.models.<m> --out base.pkl --atol model --rtol model`).
2. **The tables.** The id migration table (old id, new id, note; every id of the experiments, figures, PK modules, tests and sensitivity) and the fit parameter table (old parameter and adopted value, new parameter, unit, bootstrap value, bounds, with the conversions above: `ftissue`, per-litre rates, 1/min times 60 is 1/hr, bounds converted the same way).
3. **Failing tests first.** The model joins the library models of the tests of pkdb_models, its structural tests and its migrated scientific tests. The scientific tests keep every literature value and tolerance; only ids, observables, units of an id and model constants change. An assertion which is no literature value but a property of the old model (101 ml/min as healthy GFR, `f_cirrhosis` as a shunt, a fitted absorbed fraction, a cumulative dose state) is rewritten to the same claim on the library physiology before the refit, and named in the commit.
4. **The package.** The definition, the factory with `create_archive`, the base experiment on `PBPKSimulationExperiment`, the experiments, mapping collections, scans, PK module and sensitivity analysis on the library ids, the sampling sites of the mappings and their `check_sampling_sites` test. Delete the body, liver, kidney and intestine modules and their outputs. While bootstrapping the model has no adopted fits; a bootstrap value is a `Value` with a note naming the old fit.
5. **Compare.** Dump the bootstrap and compare it with the old model (`pkpdsim-verify compare base.pkl bootstrap.pkl --renames renames.json`, the renames from the id table): an unexpected difference is a migration error and is fixed, never listed as expected.
6. **Refit.** The fit parameters per the table, a smoke start, the fits (`run_fit <package> --subset <S> --method DE --runs 3 --seed 1234 --cores 4 --identifiability`, least squares with 24 starts where differential evolution does not finish within 6 hours), the adoption (`run_fit <package> --adopt <fit id>`), the old records deleted, the definition takes every fitted value as `fitted(FITS.parameter(pid, unit=...))`; a PD fit on the adopted PK.
7. **The gate.** The scientific tests pass; the training cost of the new fit is at most 1.25 times the cost of the old adopted fit re-evaluated on the current mappings (`run_report <package> --adopted` in the old checkout); the validation `within_2fold` not lower by more than 10 percentage points; no point of a profile likelihood below the cost of the fit (else refit from that point).
8. **Slow tests and figures.** `pytest -m slow --models <m> -n 4 --record-durations`, the experiments run, a look at the figures (legends, axes, units, curves on the data).
9. **TODO.md.** Delete the items the library resolved, a P1 item with the numbers (simulated against literature) for every strict xfail, the archive figure, women on the male build.

## When the gate fails

1. Reproduce with the smallest simulation, compare with the old model and name the mechanism.
2. A library deficit (an equation, a missing generic component, a wrong default) is fixed in the library, with a failing test first; a mechanism of one model is a `CustomReaction` or a PD component of the package.
3. A structural deficit the old model shares (its TODO.md): try the library option which addresses it (`Mode.PERMEABILITY`, per-organ Kp overrides, the interstitial shares, an exchange `ps`, a `CustomReaction` pool) and refit.
4. Only then a strict xfail (`pytest.mark.xfail(strict=True, raises=AssertionError, reason="structural limitation, TODO.md P1: <item>")`) with the TODO.md item; a literature value is never loosened.

What the seven models taught:

- **Measured clearances set the GFR.** A group whose GFR or inulin clearance was measured has that clearance as its GFR (`renal_function(clearance, bw=..., height=...)` on the body of the group); the population GFR of the library must not override the reference measurement (inulin, Prescott1990 and Prescott1991). A clearance per 70 kg or per 1.73 m2 is scaled to the group before.
- **Bounds from sources.** A parameter with a source is bounded by it (the slow muscle pool of inulin by Odeh1993, the hepatocyte transporters of gdeobdtpa by Weiss2021); a first order rate of a linear model bounds a saturable process at the operating concentration of its protocol, `Vmax = k*(Km + C_op)`. A fit which leaves a defining property of the substance (the hepatocyte accumulation of a hepatobiliary agent) is rejected even at a lower cost.
- **The interstitial share is the lower bound of a Kp.** An extracellular or highly bound substance distributes into the interstitial fluid of every organ, `Kp(interstitial=...)`; a single Kp below the shares (or many times the share of the brain, which an extracellular substance does not enter) is an artefact of the old rest compartment. The full form needs a Kp above the share of its organ.
- **Kept fits are checked.** A model whose adopted fit predates a change of a library equation is re-evaluated (`run_report --adopted`); a rise of the cost above 1 % triggers a refit, which is adopted only if it is better, not merely equivalent.
- **A window of a test is read from the publication before the model is evaluated against it.**
