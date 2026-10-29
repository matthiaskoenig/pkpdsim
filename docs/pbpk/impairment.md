# Impairment

Impaired organ function is a set of parameter vectors on the physiology, never a structurally different model: three severities of the top model set factors on the quantities of the healthy individual, the liver and the binding to plasma proteins. All severities 0 give the healthy model, every factor is 1. The numbers are those of `data/impairment.toml`, every number with its source; a number which no publication read gives is a choice of the library (`ModelChoice`) with a note.

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/impairment.svg"
</div>

| condition | severity | classes | acts on |
|---|---|---|---|
| liver cirrhosis | `CP_severity` | 0 healthy, 1, 2, 3 the Child-Pugh classes A, B, C, interpolated linearly | functional liver `f_functional`, shunt `f_shunt`, hepatic artery and kidney flow, cardiac index, GFR, albumin and AGP (so fu), the enzymes and transporters of the liver; hematocrit (factor 1, see [Limitations](#limitations)) |
| chronic kidney disease | `f_renal_function` | the GFR relative to healthy; its GFR per 1.73 m2 (`eGFR`) gives the stage | GFR, the binding capacity of the plasma proteins (so fu); hematocrit, kidney volume and flow, gastric emptying and small intestinal transit (factor 1, see [Limitations](#limitations)); with the switch `CKD_hepatic = 1` hepatic CYP2D6 and OATP1B |
| heart failure | `HF_severity` | 0 healthy, 1 to 4 the NYHA classes I to IV, interpolated linearly | the blood flows of liver and portal organs, kidney, skin, adipose and muscle; the cardiac output is their sum, the GFR follows the renal blood flow |

Each condition has its factors `fcp_<key>`, `fckd_<key>`, `fhf_<key>`; the factor of all conditions on a quantity is `fimp_<id>`, their product (`fimp_GFR`, `fimp_Qki`, `fimp_HCT`), which the [physiology](physiology.md#scaling-to-the-individual) multiplies. Every organ flow is its fraction of the cardiac output without impairment (`CO_base`) times `fimp_Q<o>`, and the cardiac output is `CO_base*fimp_CO`, so the flows sum to it. A rate of a protein in the liver runs with `f_<protein>*fimp_<protein>`, the base activity of the organ (genotype, interactions of an experiment) times the impairment, so both act multiplicatively.

The severities are initial values like a covariate: set them before the simulation, in the experiment. A severity which changes during a simulation (a study which changes the renal function between periods) is named in `PBPKDefinition.dynamic`, so its factors are assignment rules. `f_shunt` and `f_functional` stay parameters an experiment may set directly.

```equations title="factors of all conditions (male build)"
f_shunt = fcp_shunt
f_functional = fcp_functional_liver
fimp_Qli = fcp_hepatic_artery_flow * fhf_liver_flow
fimp_Qki = fcp_kidney_flow * fckd_kidney_flow * fhf_kidney_flow
fimp_Qmu = fcp_rest_flow * fhf_peripheral_flow
fimp_GFR = fcp_gfr * fhf_kidney_flow
fimp_HCT = fcp_hematocrit * fckd_hematocrit
fimp_li_cyp3a4 = fcp_cyp3a4 * (1 + CKD_hepatic * (fckd_cyp3a4 - 1))
fu_<s> = 1 / (1 + fimp_albumin / fimp_plasma_protein * (1 - fu_ref_<s>) / fu_ref_<s>)
```

The unbound fraction follows the binding (McNamara 2002): `fu = 1/(1 + f*(1 - fu_ref)/fu_ref)` with the binding factor `f`, the factor of the binding protein (`fimp_albumin`, `fimp_agp`, 1 for an unknown protein) over the CKD factor of the binding capacity `fimp_plasma_protein`. The fu of a strongly bound substance rises by the CKD factor, and fu stays below 1. A substance whose binding is measured unchanged in renal impairment (apixaban: 92-94 % bound in healthy and ESRD subjects, US label) opts out with `Binding(fu, protein, ckd=False)`: the formula then has no `fimp_plasma_protein`, `fu = 1/(1 + fimp_albumin * (1 - fu_ref)/fu_ref)`, so fu follows only the binding protein (cirrhosis still acts). With an unknown protein as well no factor is left and `fu = fu_ref`, the reference itself, which the notes of `fu_<s>` say (`unbound_fraction_notes`).

## Liver cirrhosis

`CP_severity` between the classes interpolates every factor linearly, so scans and fits use intermediate severities:

```equations title="the shunt of the Child-Pugh classes, interpolated"
fcp_shunt = piecewise(0, CP_severity <= 0, 0 + 0.399489795918367 * (CP_severity - 0), CP_severity < 1, 0.399489795918367 + 0.298469387755102 * (CP_severity - 1), CP_severity < 2, 0.697959183673469 + 0.114795918367347 * (CP_severity - 2), CP_severity < 3, 0.812755102040816)
fcp_rest_flow = 1 + (fcp_cardiac_index - 1 - (fQli * (fcp_hepatic_artery_flow - 1) + fQki * (fcp_kidney_flow - 1))) / (fQhe + fQmu + fQsk + fQad + fQbo + fQgo)
```

The organs without a factor of their own (heart, muscle, skin, adipose, bone, gonads) close the cardiac output to the cardiac index of the class (`fcp_rest_flow`). The brain (autoregulation of the cerebral blood flow, a choice of the library) and the portal organs (the shunt represents the portal blood which does not reach the liver) keep their flow.

<!-- generated: impairment-child-pugh (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/impairment.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

| factor | A | B | C | source |
| --- | ---: | ---: | ---: | --- |
| `functional_liver` | 0.69 | 0.55 | 0.28 | `ModelChoice` |
| `shunt` | 0.3995 | 0.698 | 0.8128 | `ModelChoice` |
| `hepatic_artery_flow` | 1.3 | 2.3 | 3.4 | `ModelChoice` |
| `kidney_flow` | 0.88 | 0.65 | 0.48 | `ModelChoice` |
| `cardiac_index` | 1.11 | 1.27 | 1.36 | `ModelChoice` |
| `gfr` | 1 | 0.7 | 0.36 | `ModelChoice` |
| `albumin` | 0.81 | 0.68 | 0.5 | `ModelChoice` |
| `agp` | 0.6 | 0.56 | 0.3 | `ModelChoice` |
| `hematocrit` | 1 | 1 | 1 | `ModelChoice` |

??? note "Sources"

    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
<!-- end generated -->

The physiological factors are those which Alasmari 2024 (Table 3) attribute to Edginton & Willmann 2008. Edginton 2008 was not read, so they are choices of the library, not verified against the paper, and are named in the notes as such. The shunt and the portal blood flow describe the same portal blood which does not reach the sinusoids; applied together they would count it twice. The library takes the shunt: the shunt fractions are the degrees of cirrhosis of the Child-Pugh classes of `pkpdsim.impairment.CIRRHOSIS` (the ICG models of the König lab, a choice of the library), the splanchnic flow of the portal organs is not reduced, and the fraction `f_shunt` of the liver inflow bypasses the liver, so the portal blood which perfuses the liver is 0.60, 0.30 and 0.19 of healthy (Edginton 2008: 0.40, 0.36, 0.04). The shunt acts on the hepatic arterial inflow as well.

The enzymes and transporters are the median abundances of El-Khateeb 2021 (the only source with all three classes); a class the figures do not annotate is the linear interpolation between healthy and the next annotated class, a choice of the library. Johnson 2010 (the abundance factors of the Simcyp cirrhosis populations) was not accessible and is not used; Prasad 2018 agrees with El-Khateeb 2021 for class C within the variability, but for CYP2D6, CYP2C8 and CYP2C9.

<!-- generated: impairment-child-pugh-proteins (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/impairment.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

| protein | A | B | C | source |
| --- | ---: | ---: | ---: | --- |
| `CYP1A2` | 0.69* | 0.38 | 0.13 | `ElKhateeb2021` |
| `CYP2A6` | 0.725* | 0.45 | 0.28 | `ElKhateeb2021` |
| `CYP2B6` | 1 | 1 | 1 | `ElKhateeb2021` |
| `CYP2C8` | 0.685* | 0.37 | 0.25 | `ElKhateeb2021` |
| `CYP2C9` | 0.73* | 0.46 | 0.23 | `ElKhateeb2021` |
| `CYP2C18` | 0.8* | 0.6* | 0.4 | `ElKhateeb2021` |
| `CYP2C19` | 0.51* | 0.02 | 0.05 | `ElKhateeb2021` |
| `CYP2D6` | 0.7* | 0.4 | 0.3 | `ElKhateeb2021` |
| `CYP2E1` | 0.695* | 0.39 | 0.14 | `ElKhateeb2021` |
| `CYP3A4` | 0.615* | 0.23 | 0.2 | `ElKhateeb2021` |
| `UGT1A1` | 1 | 1 | 1 | `ElKhateeb2021` |
| `UGT1A9` | 0.74* | 0.48 | 0.28 | `ElKhateeb2021` |
| `UGT2B7` | 0.65* | 0.3 | 0.13 | `ElKhateeb2021` |
| `SLCO1B1` | 0.46 | 0.42 | 0.3 | `ElKhateeb2021` |
| `SLCO1B3` | 1 | 1 | 1 | `ElKhateeb2021` |
| `ABCC2` | 0.46 | 0.34 | 0.18 | `ElKhateeb2021` |

A value marked `*` is a choice of the library, not a value of the source.

??? note "Sources"

    - `ElKhateeb2021`: El-Khateeb E, Achour B, Al-Majdoub ZM, Barber J, Rostami-Hodjegan A. Non-uniformity of changes in drug-metabolizing enzymes and transporters in liver cirrhosis: implications for drug dosage adjustment. Mol Pharm 18(9):3563-3577, 2021. [doi:10.1021/acs.molpharmaceut.1c00462](https://doi.org/10.1021/acs.molpharmaceut.1c00462)
<!-- end generated -->

## Chronic kidney disease

`f_renal_function` is the factor of the GFR, as in the hand-written models (`KI__f_renal_function` there, `f_renal_function` of the top model here): `GFR = GFR_base * f_renal_function`, with `GFR_base` the GFR at renal function 1, which carries the body surface area and the factor `fimp_GFR` of cirrhosis and heart failure (see [Physiology](physiology.md#scaling-to-the-individual)). Its GFR per 1.73 m2 of the reference individual of the build (`eGFR`) gives the stage: normal at least 90, mild 60-89, moderate 30-59, severe 15-29, end stage below 15 ml/min/1.73 m2. The factors are piecewise linear in the eGFR through anchor points at the stage midpoints: 1 from 90, the value of mild CKD at 75, of moderate CKD at 45, of severe CKD at 22.5 and of end stage CKD at 7.5 ml/min/1.73 m2 and below, linear in between (from 1 at 90 to the mild value at 75), as the factors of the Child-Pugh classes are between the classes. A scan or fit of `f_renal_function` therefore sees no step at a stage limit. The anchors are a choice of the library (`ModelChoice`: continuous covariate for scans and fits; stage midpoints), the values of the stages and their sources are unchanged.

The renal function scales the tubular capacity of the kidney as well, continuously and not by stage: every metabolism and transport of the kidney (a tubular secretion, a reabsorption by a transporter, a renal metabolism) runs with `1 - f_tubular + f_tubular*f_renal_function` of its healthy capacity (see [Chemistry](chemistry.md#kinetics)). The share `f_tubular` (a constant of the kidney, `KI__f_tubular`) is 1, the intact nephron hypothesis (Bricker 1960): the nephrons which remain in a diseased kidney work normally, so the tubular capacity falls with their number as the GFR does, as the secretion of a fraction of the GFR of the hand-written models did. The value is a choice of the library, the publication is not verified; the filtration `fu*GFR` follows the GFR itself.

```equations title="the stage and the binding capacity (male build)"
eGFR = 0.114467156025102 l_per_min * f_renal_function
fckd_plasma_protein = piecewise(1, eGFR >= 0.09, 1, eGFR >= 0.075, 1 + 0.45 * (0.075 - eGFR) / 0.03, eGFR >= 0.045, 1.45 + 0.75 * (0.045 - eGFR) / 0.0225, eGFR >= 0.0225, 2.2 + 0.3 * (0.0225 - eGFR) / 0.015, eGFR >= 0.0075, 2.5)
```

<!-- generated: impairment-ckd (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/impairment.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

Stage limits of the GFR per 1.73 m2 (ml/min/1.73 m2): 90, 60, 30, 15 (`Yoshida2016`).

Anchors of the interpolation, the eGFR (ml/min/1.73 m2) of the values of mild, moderate, severe and end stage: 75, 45, 22.5, 7.5 (`ModelChoice`).

| factor | mild | moderate | severe | end stage | source |
| --- | ---: | ---: | ---: | ---: | --- |
| `plasma_protein` | 1* | 1.45 | 2.2 | 2.5 | `Alasmari2023` |
| `hematocrit` | 1 | 1 | 1 | 1 | `ModelChoice` |
| `kidney_volume` | 1 | 1 | 1 | 1 | `ModelChoice` |
| `kidney_flow` | 1 | 1 | 1 | 1 | `ModelChoice` |
| `gastric_emptying_time` | 1 | 1 | 1 | 1 | `ModelChoice` |
| `si_transit_time` | 1 | 1 | 1 | 1 | `ModelChoice` |

A value marked `*` is a choice of the library, not a value of the source.

Share of the tubular capacity which follows the renal function (`f_tubular`, not per stage): 1 (`ModelChoice`).

??? note "Sources"

    - `Alasmari2023`: Alasmari F, Alasmari MS, Muwainea HM, Alomar HA, Alasmari AF, Alsanea S, Alshamsan A, Rasool MF, Alqahtani F. Physiologically-based pharmacokinetic modeling for single and multiple dosing regimens of ceftriaxone in healthy and chronic kidney disease populations: a tool for model-informed precision dosing. Front Pharmacol 14:1200828, 2023. [doi:10.3389/fphar.2023.1200828](https://doi.org/10.3389/fphar.2023.1200828)
    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
    - `Yoshida2016`: Yoshida K, Sun B, Zhang L, Zhao P, Abernethy DR, Nolin TD, Rostami-Hodjegan A, Zineh I, Huang SM. Systematic and quantitative assessment of the effect of chronic kidney disease on CYP2D6 and CYP3A4/5. Clin Pharmacol Ther 100(1):75-87, 2016. [doi:10.1002/cpt.337](https://doi.org/10.1002/cpt.337)
<!-- end generated -->

The binding effect is a factor on the binding capacity in the form of fu above, not a factor on fu: binding capacity is the mechanistic quantity and keeps fu below 1 (at `fu_ref` 0.5 the end stage gives fu 0.71). The hepatic proteins (switch `CKD_hepatic`, off by default) are the mean ratios of the unbound clearance of model drugs (Yoshida 2016, Tan 2018), taken as the activity of the pathway:

<!-- generated: impairment-ckd-proteins (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/impairment.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

| protein | mild | moderate | severe | end stage | source |
| --- | ---: | ---: | ---: | ---: | --- |
| `CYP2D6` | 1* | 0.53 | 0.41 | 0.5 | `Yoshida2016` |
| `CYP3A4` | 1 | 1 | 1 | 1 | `Yoshida2016` |
| `SLCO1B1` | 0.68 | 0.56 | 0.54 | 0.45 | `Tan2018` |
| `SLCO1B3` | 0.68 | 0.56 | 0.54 | 0.45 | `Tan2018` |

A value marked `*` is a choice of the library, not a value of the source.

??? note "Sources"

    - `Tan2018`: Tan ML, Yoshida K, Zhao P, Zhang L, Nolin TD, Piquette-Miller M, Galetin A, Huang SM. Effect of chronic kidney disease on nonrenal elimination pathways: a systematic assessment of CYP1A2, CYP2C8, CYP2C9, CYP2C19, and OATP. Clin Pharmacol Ther 103(5):854-867, 2018. [doi:10.1002/cpt.807](https://doi.org/10.1002/cpt.807)
    - `Yoshida2016`: Yoshida K, Sun B, Zhang L, Zhao P, Abernethy DR, Nolin TD, Rostami-Hodjegan A, Zineh I, Huang SM. Systematic and quantitative assessment of the effect of chronic kidney disease on CYP2D6 and CYP3A4/5. Clin Pharmacol Ther 100(1):75-87, 2016. [doi:10.1002/cpt.337](https://doi.org/10.1002/cpt.337)
<!-- end generated -->

## Heart failure

The flows of the organs are multiplied by the factors of the NYHA class and the cardiac output is their sum (Gu 2025, Table 1). The renal blood flow is not monotonic in the class (NYHA IV above NYHA III, as measured by Leithe 1984), and the GFR follows it.

<!-- generated: impairment-heart-failure (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/impairment.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

| factor | I | II | III | IV | source |
| --- | ---: | ---: | ---: | ---: | --- |
| `liver_flow` | 1* | 0.765 | 0.538 | 0.463 | `Rasool2015` |
| `kidney_flow` | 1* | 0.778 | 0.55 | 0.629 | `Rasool2015` |
| `peripheral_flow` | 1* | 0.57 | 0.44 | 0.28 | `Gu2025` |

A value marked `*` is a choice of the library, not a value of the source.

??? note "Sources"

    - `Gu2025`: Gu W, Shao Q, Jiang L. Predicting pharmacokinetics of drugs in patients with heart failure and optimizing their dosing strategies using a physiologically based pharmacokinetic model. Pharmaceutics 17(11):1394, 2025. [doi:10.3390/pharmaceutics17111394](https://doi.org/10.3390/pharmaceutics17111394)
    - `Rasool2015`: Rasool MF, Khalil F, Läer S. A physiologically based pharmacokinetic drug-disease model to predict carvedilol exposure in adult and paediatric heart failure patients by incorporating pathophysiological changes in hepatic and renal blood flows. Clin Pharmacokinet 54(9):943-962, 2015. [doi:10.1007/s40262-015-0253-7](https://doi.org/10.1007/s40262-015-0253-7)
<!-- end generated -->

## Classes in experiments

The experiments of a model built on the library map the classes of a study onto the severities with the presets of `pkpdsim.impairment`: `CHILD_PUGH` (colors `CIRRHOSIS_COLORS`) onto `CP_severity`, the `RENAL_FUNCTION` presets onto `f_renal_function`, `HEART_FAILURE` (colors `HEART_FAILURE_COLORS`) onto `HF_severity`. The base of their experiments is `PBPKSimulationExperiment` of [`pkpdsim.pbpk.experiment`](../api/pbpk_experiment.md), which has these maps: `renal_map` the classes of `RENAL_FUNCTION_GFR` (colors `RENAL_COLORS_END_STAGE`), `cirrhosis_map` `CHILD_PUGH` and `cardiac_map` `HEART_FAILURE`. The classes of renal impairment keep their GFR in ml/min: a preset is relative to a healthy GFR of 101 ml/min, the renal function of a class is its GFR over the healthy GFR `GFR_base` of the model (`ScaledRenalMap`, evaluated on the first access, so importing the experiments loads no model). This differs from the eGFR of the CKD factors above, which is relative to the GFR per 1.73 m2 of the reference individual. `healthy_gfr` reads `GFR_base` of an individual (body weight, height, `HF_severity`, `CP_severity`) from the flat model, a quantity in ml/min, and `renal_function(clearance, bw=..., height=...)` gives the renal function of a measured clearance in ml/min as a dimensionless quantity, so that the simulated GFR of the individual of the study is the measured one; both are quantities, which go into the changes of a simulation as they are. The renal function is limited by the validity of the kidney outflow (see [Standard organ](organ.md#validity-of-the-outflow)): the net filtered plasma flow of every filtered substance, `fu*GFR*(1 - f_reabsorbed)` (`fu*GFR` with a threshold reabsorption), must stay below the plasma flow of the kidney, `Qki*(1 - HCT)` (`Qki*BP` with erythrocytes in equilibrium). `kidney_outflow(model_path, changes)` of [`pkpdsim.pbpk.kidney_validity`](../api/pbpk_kidney_validity.md) gives both flows of an individual (`KidneyOutflow`, in ml/min, with the formula of its net filtered plasma flow, read as `healthy_gfr` reads the GFR (`individual_values`), from a simulation of one minute on a model loaded once per process, cached per model file and individual; only the changes of the parameters the flows depend on make the individual, a dose does not; a model without a filtered substance is not simulated). `renal_function` raises a `ValueError` which names the GFR and the plasma flow for a renal function beyond it, `renal_function_limits(bw=..., height=..., ...)` gives the range of an individual, from 0 to the renal function at which the net filtered plasma flow reaches the plasma flow (about 4 to 8 for an unbound substance), and a model narrows it with `renal_function_range` (`None` by default, e.g. `(0.05, 2.0)` for a model whose urine is not defined without filtration). The initialization of every experiment checks every individual it simulates (`simulated_individuals`), each task on the model it runs, at the start of every simulation and every point of a scan: its preinit changes, of a presimulation as well, and the values of the dimensions without a time. It does not check a change during a simulation (`Simulation.changes`, also at its start, of an id the definition names `dynamic`), a dimension with a time (`at`), a scan over models and the draws of the fitted parameters of a run over them; the margin is large (about 6 times at renal function 1). It fails with a message which names the experiment, the simulation, the GFR, the plasma flow and the changes of the individual. The names of `pkpdsim.pbpk.experiment` and `pkpdsim.pbpk.kidney_validity` are imported from `pkpdsim.pbpk` (`from pkpdsim.pbpk import REFERENCE_GFR, PBPKSimulationExperiment, kidney_outflow`), which exports them lazily. The 44 hand-written models keep their maps and ids unchanged: `f_cirrhosis` keeps its meaning there (shunt and tissue loss in one), and the library uses other ids.

## Limitations

- The physiological factors of the Child-Pugh classes are choices of the library from a review table (Alasmari 2024), not verified against Edginton & Willmann 2008.
- The anaemia of cirrhosis and CKD: the publications read state absolute hematocrits without the healthy value of their population, so the hematocrit factors are 1 until a source with its healthy reference is read.
- The kidney volume and blood flow and the gastric emptying and small intestinal transit of CKD are factor 1: the constants which reach a PBPK library only through the code of PK-Sim (GPL) are not shipped, and Malik 2020 was not read. A dropped effect is factor 1 with a note in the data file. The CKD physiology is thinner than that of PK-Sim until Malik 2020 is read.
- The binding factor of CKD acts on the binding capacity, PK-Sim multiplies fu by it: the two differ for weakly bound drugs (at `fu_ref` 0.5 and the factor 1.45 of moderate CKD, fu 0.59 against 0.725).
- A Kp of a Kp method is a number of the build: a change of fu by an impairment does not recompute it.
- The tubular capacity follows the renal function in proportion (`f_tubular` 1, the intact nephron hypothesis, Bricker 1960, not verified, a choice of the library); a transporter whose activity falls less or more than the GFR in CKD needs another share or an activity factor of the experiment.
