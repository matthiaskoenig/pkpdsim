# Physiology

The standard human physiology is a set of SBML parameters of the top model: the covariates of the individual, the reference values of the sex of the build (ICRP 89) and of the age, and everything which follows from them, the organ masses and volumes, the blood volume, the cardiac output, the blood flows and the GFR. An experiment changes the covariates (`BW`, `HEIGHT`, `AGE`) and the factors of organ function (`f_renal_function`, `f_cardiac_function`, the [impairment](impairment.md) severities) before the simulation, as it changes `BW` in the hand-written models, and the organs follow, from the newborn to 100 years ([physiology over age](#physiology-over-age)); [virtual populations](populations.md) combine these factors.

## Topology

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/circulation.svg"
</div>

The circulation of Jones & Rowland-Yeo 2013 and PK-Sim (Kuepfer 2016): an arterial and a venous blood pool, the lung in series (its flow is the cardiac output), the other organs in parallel between the arterial and the venous pool. Stomach, small intestine, large intestine, spleen and pancreas drain into the portal vein, which feeds the liver together with the hepatic artery; the liver drains into the venous pool. The organs (`pkpdsim.pbpk.Organ`, their two-letter tokens are the prefixes of their ids) are lung `lu`, heart `he`, brain `br`, liver `li`, kidney `ki`, stomach `st`, small intestine `si`, large intestine `la`, spleen `sp`, pancreas `pa`, muscle `mu`, skin `sk`, adipose `ad`, bone `bo` and gonads `go`; the blood pools (`Pool`) are arterial `ar`, venous `ve` and portal `po`.

Every organ owns its blood flows: it reads the plasma (and blood cell) species of the pool which supplies it and writes into the pool it drains into, through its ports (see [Standard organ](organ.md)). The pools have no reactions of their own besides the exchange with kinetic blood cells and the intravenous dose. The intestine variants replace small and large intestine with one submodel `GU` with the inflow of both.

## Reference values

The organ masses, densities, blood contents and blood flows are those of the ICRP 89 reference individuals (Valentin 2002), taken directly from its tables; the interstitial fractions are those of Schmitt 2008 as tabulated by httk 2.7.4 (`tissue.data`, MIT), where they were read. The data file `database/tables/physiology_human.toml` compares the male values with Jones & Rowland-Yeo 2013 in the notes of every value. The conventions of the data:

- `mass` is the mass of the organ without its blood (the blood is counted on its own); `density` is the specific gravity (1 kg/l where ICRP 89 states none); `blood` is the blood content of the organ as fraction of the blood volume.
- The organ volume of the model is the volume of its tissue and its blood, `V<o> = M<o>/density + blood*Vblood`; the vascular fraction is its blood, `f_vas_<o> = blood*Vblood/V<o>`, so that the vascular spaces of the organs and the pools are the blood volume for every individual.
- `f_int` is the interstitial fraction of the organ volume, the interstitial fraction of the tissue without blood of Schmitt 2008 (as tabulated by httk) times `1 - f_vas`; the intracellular fraction is `1 - f_vas - f_int`.
- `flow` is the blood flow as fraction of the cardiac output (ICRP 89 Table 2.40); the liver states its arterial flow (`Qli` is the hepatic artery, `Qha`), the portal organs feed it in addition. ICRP 89 gives about 8 % of the cardiac output to tissues which are no organ of the model (`[unmapped]` of the data): the flow fractions of the model, `fQ<o>`, are renormalised over the organs, so they sum to the cardiac output.
- The portal vein pool is the anatomical portal vein (Chau 2023, 4.4 ml); the splanchnic venous blood stays in the vascular volume of the portal organs. The blood of the unmapped tissues is part of the venous pool.
- Adipose closes the body weight: its mass is the remainder, `Mad = BW - Mlean - Vblood*rho_blood`; the ICRP adipose mass in the data is a reference value only.

<!-- generated: physiology-general (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/physiology_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| quantity | unit | male | female |
| --- | --- | ---: | ---: |
| body weight | kg | 73 | 60 |
| height | cm | 176 | 163 |
| body surface area | m2 | 1.9 | 1.66 |
| cardiac output | l/min | 6.5 | 5.9 |
| hematocrit | - | 0.43 | 0.38 |
| GFR | ml/min | 125 | 109 |
| hepatic artery share of the liver flow | - | 0.2549 | 0.2364 |
| blood volume | l | 5.3 | 3.9 |
| blood density | kg/l | 1.06 | 1.06 |
| arterial pool (fraction of blood) | - | 0.105 | 0.105 |
| venous pool (fraction of blood) | - | 0.2242 | 0.2239 |
| portal vein (fraction of blood) | - | 0.000828 | 0.00113 |

??? note "Sources"

    - `Chau2023`: Chau P, Yoon JS, Moses D, Pather N. A systematic review and meta-analysis of portal vein morphometry in pediatric and adult populations: drawing the line between normal and abnormal findings. Eur J Radiol 168:111016, 2023. [doi:10.1016/j.ejrad.2023.111016](https://doi.org/10.1016/j.ejrad.2023.111016)
    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
<!-- end generated -->

The numbers of the reference individuals are public, so that a model which states a value of a hand-written model per kg of body weight or per litre of an organ converts it with them and copies none: `reference_individual(Sex.MALE)` has the body weight `bw` (73 kg), the height `height` (176 cm), the healthy GFR `gfr` (0.125 l/min, 0.109 l/min of a woman: `GFR_base` of the default individual, the GFR at renal function 1 which `healthy_gfr` of the build of the sex gives in ml/min; the GFR of a component a model keeps outside the library, the fluid model of hctz, reads it), the organ volumes with their blood `volumes` (`V<o>`, the liver 2.33 l), the tissue volumes without blood `tissue_volumes`, the volumes of the cells `intracellular_volumes` (`V<o>*(1 - f_vas_<o> - f_int_<o>)`, the intracellular space `Vic` of the full form, the hepatocytes 1.3859 l) and the fluid volumes of the lumen segments `lumen_volumes` (`Vlumen_<seg>` by the tokens of `ids.GUT_SEGMENTS`, the same for both sexes), the values of the default individual of a male build (ICRP 89). A model writes no copy of these numbers (`73.0`, `2.33`, `1.3859`, `0.125`) into its definition or its fit parameters, it imports them.

```python
from pkpdsim.pbpk import Organ, Sex, reference_individual

man = reference_individual(Sex.MALE)
print(man.bw, man.height, round(man.volumes[Organ.LIVER], 3))
small = sum(man.lumen_volumes[seg] for seg in ("du", "jej1", "jej2", "ile1", "ile2"))
print(round(man.intracellular_volumes[Organ.LIVER], 4), round(small, 5))
print(man.gfr, round(reference_individual(Sex.FEMALE).gfr, 4))
```

```text
73.0 176.0 2.33
1.3859 0.105
0.125 0.109
```

### Male reference individual

<!-- generated: physiology-organs-male (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/physiology_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| organ | mass [kg] | density [kg/l] | blood [-] | f_vas [-] | f_int [-] | flow [-] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `lu` lung | 0.5 | 1 | 0.125 | 0.5699 | 0.08602 | 1 |
| `he` heart | 0.33 | 1 | 0.01 | 0.1384 | 0.2154 | 0.04 |
| `br` brain | 1.45 | 1.04 | 0.012 | 0.04363 | 0.004782 | 0.12 |
| `li` liver | 1.8 | 1 | 0.1 | 0.2275 | 0.1777 | 0.065 |
| `ki` kidney | 0.31 | 1.05 | 0.02 | 0.2642 | 0.1214 | 0.19 |
| `st` stomach | 0.15 | 1.045 | 0.01 | 0.2697 | 0.07011 | 0.01 |
| `si` small intestine | 0.65 | 1.045 | 0.038 | 0.2446 | 0.07252 | 0.1 |
| `la` large intestine | 0.37 | 1.045 | 0.022 | 0.2477 | 0.07222 | 0.04 |
| `sp` spleen | 0.15 | 1.06 | 0.014 | 0.344 | 0.1673 | 0.03 |
| `pa` pancreas | 0.14 | 1 | 0.006 | 0.1851 | 0.07823 | 0.01 |
| `mu` muscle | 29 | 1 | 0.14 | 0.02495 | 0.1463 | 0.17 |
| `sk` skin | 3.3 | 1.1 | 0.03 | 0.05033 | 0.5698 | 0.05 |
| `ad` adipose | 14.5 | 1 | 0.05 | 0.01795 | 0.1375 | 0.05 |
| `bo` bone | 10.5 | 1.3 | 0.07 | 0.04392 | 0.09561 | 0.05 |
| `go` gonads | 0.035 | 1.04 | 0.0004 | 0.05926 | 0.09031 | 0.0005 |

Sources of the columns: mass, density, blood, f_vas, flow: `ICRP89`; f_int: `httk`.

??? note "Sources"

    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `httk`: Pearce RG, Setzer RW, Strope CL, Wambaugh JF, Sipes NS. httk: R package for high-throughput toxicokinetics. J Stat Softw 79(4):1-26, 2017. R package httk version 2.7.4 (CRAN, 2025-12-08), MIT license, copyright 2025 John Wambaugh; the data set mecdt of HTTK-Pop (Ring 2017) holds the preprocessed NHANES 2013-2018 respondents. [doi:10.18637/jss.v079.i04](https://doi.org/10.18637/jss.v079.i04)
<!-- end generated -->

### Female reference individual

<!-- generated: physiology-organs-female (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/physiology_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| organ | mass [kg] | density [kg/l] | blood [-] | f_vas [-] | f_int [-] | flow [-] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `lu` lung | 0.42 | 1 | 0.125 | 0.5372 | 0.09256 | 1 |
| `he` heart | 0.25 | 1 | 0.01 | 0.1349 | 0.2163 | 0.05 |
| `br` brain | 1.3 | 1.04 | 0.012 | 0.03609 | 0.00482 | 0.12 |
| `li` liver | 1.4 | 1 | 0.1 | 0.2179 | 0.1799 | 0.065 |
| `ki` kidney | 0.275 | 1.05 | 0.02 | 0.2295 | 0.1271 | 0.17 |
| `st` stomach | 0.14 | 1.045 | 0.01 | 0.2255 | 0.07435 | 0.01 |
| `si` small intestine | 0.6 | 1.045 | 0.038 | 0.2052 | 0.0763 | 0.11 |
| `la` large intestine | 0.36 | 1.045 | 0.022 | 0.1994 | 0.07686 | 0.05 |
| `sp` spleen | 0.13 | 1.06 | 0.014 | 0.3081 | 0.1764 | 0.03 |
| `pa` pancreas | 0.12 | 1 | 0.006 | 0.1632 | 0.08033 | 0.01 |
| `mu` muscle | 17.5 | 1 | 0.105 | 0.02286 | 0.1466 | 0.12 |
| `sk` skin | 2.3 | 1.1 | 0.03 | 0.05299 | 0.5682 | 0.05 |
| `ad` adipose | 19 | 1 | 0.085 | 0.01715 | 0.1376 | 0.085 |
| `bo` bone | 7.8 | 1.3 | 0.07 | 0.04352 | 0.09565 | 0.05 |
| `go` gonads | 0.011 | 1.05 | 0.0002 | 0.0693 | 0.08935 | 0.0002 |

Sources of the columns: mass, density, blood, f_vas, flow: `ICRP89`; f_int: `httk`.

??? note "Sources"

    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `httk`: Pearce RG, Setzer RW, Strope CL, Wambaugh JF, Sipes NS. httk: R package for high-throughput toxicokinetics. J Stat Softw 79(4):1-26, 2017. R package httk version 2.7.4 (CRAN, 2025-12-08), MIT license, copyright 2025 John Wambaugh; the data set mecdt of HTTK-Pop (Ring 2017) holds the preprocessed NHANES 2013-2018 respondents. [doi:10.18637/jss.v079.i04](https://doi.org/10.18637/jss.v079.i04)
<!-- end generated -->

## Blood and plasma

Every vascular space has plasma `V_vas*(1 - HCT)` and blood cells `V_vas*HCT`, and the blood flows carry plasma, `Q*(1 - HCT)*C_p`, and blood cells, `Q*HCT*C_bc`. The hematocrit is that of the sex (0.43 male, 0.38 female, ICRP 89), not the uncited 0.51 of some hand-written models. The blood concentration and the blood to plasma ratio are observers, never a conversion written by hand; how a substance enters the erythrocytes is its `erythrocytes` flag (see [Substances](substances.md#erythrocytes)).

```equations title="blood pools (male build)"
HCT = HCT_ref * fimp_HCT
Vve = fVve * Vblood
Vve_pl = Vve * (1 - HCT)
Vve_bc = Vve * HCT
```

## Scaling to the individual

The covariates are the body weight `BW` (kg), the height `HEIGHT` (cm) and the age `AGE` (yr); their defaults are the reference individual of the sex of the population of the build, of ICRP 89 the adult of 35 years. The scaling is relative to the reference individual of the age (`BW_ref`, `HEIGHT_ref`, `M<o>_ref`, `Vblood_ref`, `CO_ref`, functions of `AGE`, see below), so the organs of a child of the size of its age are those of the reference child. The chain follows Huisinga 2012: body mass index, body surface area (Du Bois 1916) and fat-free mass (Janmahasatian 2005); brain is constant, skin scales with the body surface area, the other lean organs and the blood with the fat-free mass less brain and skin (`SV_tis`, Huisinga 2012 eqs. 9, 10 and 12), adipose is the remainder of the body weight. The cardiac output scales as the tissues (eq. 13), the organ flows are fractions of it, and the GFR scales with the body surface area. This replaces the `ftissue*BW/BW_ref` scaling of the hand-written models. The body composition which follows from it, the split of the adipose tissue into its depots, the physiology of overweight and obesity and a change of the body weight over time are on the page [Adipose tissue and obesity](adipose.md); `pkpdsim.pbpk.physiology.body_composition(sex, bw, height, age)` gives the quantities of the chain as numbers.

```equations title="covariates"
BMI = BW / 1 kg / (HEIGHT / 100 cm)^2 * 1 kg_per_m2
BSA = 0.007184 m2 * (BW / 1 kg)^0.425 * (HEIGHT / 1 cm)^0.725
FFM = BW * 9270 / (6680 + 216 * BMI / 1 kg_per_m2) * fage_FFM
BSA_ref = 0.007184 m2 * (BW_ref / 1 kg)^0.425 * (HEIGHT_ref / 1 cm)^0.725
SV_tis = (FFM - Mbr - Msk) / (FFM_ref - Mbr_ref - Msk_ref)
```

```equations title="organs (the liver and the lung as examples, male build)"
Mbr = Mbr_ref
Msk = Msk_ref * BSA / BSA_ref
M<o> = M<o>_ref * SV_tis
Vblood = Vblood_ref * SV_tis
Mad = BW - Mlean - Vblood * 1.06 kg_per_l
Vli = Mli / 1 kg_per_l + 0.1 * Vblood
f_vas_li = 0.1 * Vblood / Vli
f_int_lu = 0.2 * (1 - f_vas_lu)
```

```equations title="blood flows and GFR"
CO_base = CO_ref * SV_tis * f_cardiac_function
CO = CO_base * fimp_CO
Qlu = CO
Q<o> = fQ<o> * CO_base * fimp_Q<o>
Qha = Qli
Qpo = Qst + Qsi + Qla + Qsp + Qpa
Qh = Qha + Qpo
GFR_base = GFR_ref * BSA / 1.88918819606715 * fage_GFR * fimp_GFR
GFR = GFR_base * f_renal_function
```

`GFR_base` is the GFR of the individual at renal function 1, with the body surface area (relative to the adult reference, 1.889 m2 of the male build), the factor of the age (`fage_GFR`) and the other impairments: a model which needs the GFR of the kidney without chronic kidney disease (a PD model whose effect scales with the filtration relative to it, a measured clearance converted into a renal function) reads it, and it is defined at renal function 0, where the GFR is 0.

The factors `fimp_<id>` are those of the [impairment](impairment.md), 1 for the healthy individual; `f_cardiac_function` and `f_renal_function` are the factors of organ function of the hand-written models, with the same ids. All of these are initial assignments on constant parameters: a covariate changed before the simulation (a change of an experiment) changes every organ, a change during the simulation does not, unless the definition names it in `PBPKDefinition.dynamic`. The kidney mass carries the CKD factor of the kidney volume (`fimp_Mki`, 1 until its source is read, see [Impairment](impairment.md#limitations)).

## Physiology over age

The physiology over age is a choice of the build, as the sex and the species: `PBPKDefinition(age=True, ...)` makes every reference value a function of the covariate `AGE`, from the newborn to 100 years. Without it (the default) the model is the reference individual of 35 years, it has no `AGE` and every factor of the age is a constant 1, so a model of adults pays nothing for the curves; a build over age needs a roadrunner with the cache of the initial value resolvers (see [Limitations](#limitations)).

The reference values of the sex are functions of the covariate `AGE` (`pkpdsim.pbpk.age`), initial assignments as every other quantity of the physiology, so an experiment which sets `AGE` before the simulation gets the physiology of the age; at 35 years, the age of the adults of ICRP 89 (Section 1.4.2 para 14), every factor is 1 and the adult of ICRP 89 is the reference individual of the earlier versions. Every curve is written once (`pkpdsim.pbpk.curves`) and evaluated as a number and as the formula of the model; `tests/pbpk/test_age.py` compares both at 22 ages to 1e-12.

- **Growth.** The reference children of ICRP 89 at birth and at 1, 5, 10 and 15 years (Tables 2.8, 2.9, 2.39 and Section 7.4, `database/tables/children_human.toml`) and the adult at 20 years: the reference body weight and height, the tissue masses of the lean organs, the blood volume and the cardiac output, linear between the ages. An individual of an age is scaled to its body weight and height relative to the reference individual of its age (Huisinga 2012, who take the reference individual of the age class; the interpolation is a choice of the library); the fat-free mass of the children carries the factor of Al-Sallami 2015 (`fage_FFM`, 1 from 20 years).
- **Aging.** The adult values times a factor of the age at a fixed body size, 1 at 35 years (`database/tables/age_human.toml`): the liver mass (Imamura 2026), the kidney mass (Wang 2014), the cardiac output (Carlsson 2012, Patel 2021), the GFR (Delanaye 2022, Eriksen 2020, per sex), the hepatic and the renal blood flow (Zoli 1999, Davies and Shock 1950); a model which needs the body size of an older population takes the body weight and height of the population of its age (the [virtual populations](populations.md)).
- **Maturation.** The GFR matures with the postmenstrual age (Rhodin 2009, `fage_GFR`, 36 % at birth, 90 % at 1 year); the hematocrit of the neonates (Scholkmann 2019) and of every age of NHANES 2013-2018 (`fage_HCT`, `HCT_ref = HCT * fage_HCT` of the adult of the sex), the albumin (Johnson 2006, NHANES) and the alpha-1-acid glycoprotein (Maharaj 2018) of the binding (`fage_albumin`, `fage_agp`, factors on the binding protein of `fu` as those of the [impairment](impairment.md)).
- **Blood flows.** The fractions of the cardiac output of the adults times the factors of the age, which keep their sum 1: the brain takes its flow per mass of the age (Chiron 1992) of its mass, up to half of the cardiac output at 1 year (`fage_Qbr`), the other organs share the rest in the proportions of the adults (`fage_Qrest`), of these the liver with the portal organs and the kidney follow their own factors of the elderly (`fage_Qli`, `fage_Qki`) and the remaining organs close the sum (`fage_Qother`).
- **Ontogeny.** The enzymes and transporters of the liver, the kidney and the small intestine mature with the age (`fage_<o>_<protein>`, `database/tables/ontogeny.toml`): the Hill functions of the postmenstrual age of the PK-Sim ontogeny database for the CYPs and UGTs, of Hunt 2024 for the transporters of the kidney and of Prasad 2016 for OCT1, OATP1B3 and P-gp of the liver; a protein without data has no ontogeny. The factor multiplies the activity of the protein in the organ with the factors of the impairment (`fimp_<o>_<protein>`).

```equations title="hematocrit, GFR and flow fractions over age (male build)"
HCT_ref = 0.43 * fage_HCT
HCT = HCT_ref * fimp_HCT
fQbr = 0.12965964343598055 * fage_Qbr
```

The reference individual of an age is public as for the adult, `reference_individual(sex, age)`:

```python
from pkpdsim.pbpk import Organ, Sex, reference_individual

for age in (0.0, 5.0, 35.0, 80.0):
    r = reference_individual(Sex.MALE, age)
    print(age, r.bw, round(r.volumes[Organ.LIVER], 3), round(r.gfr * 1000, 1))
```

```text
0.0 3.5 0.157 5.0
5.0 19.0 0.71 49.7
35.0 73.0 2.33 125.0
80.0 73.0 2.033 89.2
```

### Male reference individual over age

<!-- generated: physiology-age-male (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/children_human.toml`, `physiology_human.toml` and `age_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| age [yr] | BW [kg] | HEIGHT [cm] | BSA [m2] | Mbr [kg] | Mli [kg] | Mki [kg] | Mmu [kg] | Vblood [l] | CO [l/min] | GFR [ml/min] | GFR [ml/min/1.73 m2] | HCT [-] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| newborn | 3.5 | 51 | 0.2116 | 0.38 | 0.13 | 0.025 | 0.8 | 0.27 | 0.6 | 4.966 | 40.6 | 0.4995 |
| 1 | 10 | 76 | 0.4415 | 0.95 | 0.33 | 0.07 | 1.9 | 0.5 | 1.2 | 26.4 | 103.5 | 0.3447 |
| 5 | 19 | 109 | 0.7533 | 1.31 | 0.57 | 0.11 | 5.6 | 1.4 | 3.4 | 49.75 | 114.2 | 0.3551 |
| 10 | 32 | 138 | 1.115 | 1.4 | 0.83 | 0.18 | 11 | 2.4 | 5 | 73.79 | 114.4 | 0.374 |
| 15 | 56 | 167 | 1.625 | 1.42 | 1.3 | 0.25 | 24 | 4.5 | 6.1 | 107.5 | 114.5 | 0.4163 |
| 20 | 73 | 176 | 1.889 | 1.45 | 1.8 | 0.31 | 29 | 5.3 | 6.792 | 125 | 114.5 | 0.4331 |
| 35 | 73 | 176 | 1.889 | 1.45 | 1.8 | 0.31 | 29 | 5.3 | 6.5 | 125 | 114.5 | 0.43 |
| 50 | 73 | 176 | 1.889 | 1.45 | 1.8 | 0.31 | 29 | 5.3 | 6.207 | 114.6 | 105 | 0.4261 |
| 65 | 73 | 176 | 1.889 | 1.45 | 1.651 | 0.2798 | 29 | 5.3 | 5.915 | 101.9 | 93.32 | 0.4202 |
| 80 | 73 | 176 | 1.889 | 1.45 | 1.503 | 0.2496 | 29 | 5.3 | 5.622 | 89.2 | 81.68 | 0.403 |
| 100 | 73 | 176 | 1.889 | 1.45 | 1.305 | 0.2093 | 29 | 5.3 | 5.232 | 72.25 | 66.16 | 0.3916 |

??? note "Sources"

    - `Carlsson2012`: Carlsson M, Andersson R, Bloch KM, Steding-Ehrenborg K, Mosén H, Stahlberg F, Ekmehag B, Arheden H. Cardiac output and cardiac index measured with cardiovascular magnetic resonance in healthy subjects, elite athletes and patients with congestive heart failure. J Cardiovasc Magn Reson 14(1):51, 2012. [doi:10.1186/1532-429X-14-51](https://doi.org/10.1186/1532-429X-14-51)
    - `Delanaye2022`: Delanaye P, Gaillard F, van der Weijden J, Mjøen G, Ferhman-Ekholm I, Dubourg L, Ebert N, Schaeffner E, Åkerfeldt T, Goffin K, et al. Age-adapted percentiles of measured glomerular filtration in healthy individuals: extrapolation to living kidney donors over 65 years. Clin Chem Lab Med 60(3):401-407, 2022. [doi:10.1515/cclm-2021-1011](https://doi.org/10.1515/cclm-2021-1011)
    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `Imamura2026`: Imamura H, Ito R, Ichida H, Yoshioka R, Mise Y, Nojiri S, Sano K, Yamataka A, Saiura A. Age-range-matched external validation of standard liver volume equations: methodological re-evaluation of 13 regression models. PLOS One 21:e0351983, 2026. [doi:10.1371/journal.pone.0351983](https://doi.org/10.1371/journal.pone.0351983)
    - `Jopling2009`: Jopling J, Henry E, Wiedmeier SE, Christensen RD. Reference ranges for hematocrit and blood hemoglobin concentration during the neonatal period: data from a multihospital health care system. Pediatrics 123(2):e333-e337, 2009. [doi:10.1542/peds.2008-2654](https://doi.org/10.1542/peds.2008-2654)
    - `NHANES2013_2018`: Centers for Disease Control and Prevention, National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES) 2013-2014, 2015-2016 and 2017-2018, demographics (DEMO), body measures (BMX), complete blood count (CBC) and standard biochemistry profile (BIOPRO) of the cycles H, I and J (public domain). https://wwwn.cdc.gov/nchs/nhanes/ (accessed 2026-11-06). Read from the public files by scripts/database/nhanes_populations.py, scripts/database/nhanes_blood.py and scripts/database/nhanes_waist.py.
    - `Rhodin2009`: Rhodin MM, Anderson BJ, Peters AM, Coulthard MG, Wilkins B, Cole M, Chatelut E, Grubb A, Veal GJ, Keir MJ, Holford NHG. Human renal function maturation: a quantitative description using weight and postmenstrual age. Pediatr Nephrol 24(1):67-76, 2009. [doi:10.1007/s00467-008-0997-5](https://doi.org/10.1007/s00467-008-0997-5)
    - `Scholkmann2019`: Scholkmann F, Ostojic D, Isler H, Bassler D, Wolf M, Karen T. Reference ranges for hemoglobin and hematocrit levels in neonates as a function of gestational age (22-42 weeks) and postnatal age (0-29 days): mathematical modeling. Children 6(3):38, 2019. [doi:10.3390/children6030038](https://doi.org/10.3390/children6030038)
    - `WangX2014`: Wang X, Vrtiska TJ, Avula RT, Walters LR, Chakkera HA, Kremers WK, Lerman LO, Rule AD. Age, kidney function, and risk factors associate differently with cortical and medullary volumes of the kidney. Kidney Int 85(3):677-685, 2014. [doi:10.1038/ki.2013.359](https://doi.org/10.1038/ki.2013.359)
<!-- end generated -->

### Female reference individual over age

<!-- generated: physiology-age-female (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/children_human.toml`, `physiology_human.toml` and `age_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| age [yr] | BW [kg] | HEIGHT [cm] | BSA [m2] | Mbr [kg] | Mli [kg] | Mki [kg] | Mmu [kg] | Vblood [l] | CO [l/min] | GFR [ml/min] | GFR [ml/min/1.73 m2] | HCT [-] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| newborn | 3.5 | 51 | 0.2116 | 0.38 | 0.13 | 0.025 | 0.8 | 0.27 | 0.6 | 4.976 | 40.68 | 0.5022 |
| 1 | 10 | 76 | 0.4415 | 0.95 | 0.33 | 0.07 | 1.9 | 0.5 | 1.2 | 26.46 | 103.7 | 0.3492 |
| 5 | 19 | 109 | 0.7533 | 1.18 | 0.57 | 0.11 | 5.6 | 1.4 | 3.4 | 49.85 | 114.5 | 0.3606 |
| 10 | 32 | 138 | 1.115 | 1.22 | 0.83 | 0.18 | 11 | 2.4 | 5 | 73.94 | 114.7 | 0.3758 |
| 15 | 53 | 161 | 1.546 | 1.3 | 1.3 | 0.24 | 17 | 3.3 | 6.1 | 102.5 | 114.7 | 0.3793 |
| 20 | 60 | 163 | 1.644 | 1.3 | 1.4 | 0.275 | 17.5 | 3.9 | 6.165 | 109 | 114.7 | 0.3829 |
| 35 | 60 | 163 | 1.644 | 1.3 | 1.4 | 0.275 | 17.5 | 3.9 | 5.9 | 109 | 114.7 | 0.38 |
| 50 | 60 | 163 | 1.644 | 1.3 | 1.4 | 0.275 | 17.5 | 3.9 | 5.635 | 99.94 | 105.2 | 0.3853 |
| 65 | 60 | 163 | 1.644 | 1.3 | 1.284 | 0.2482 | 17.5 | 3.9 | 5.369 | 85.79 | 90.27 | 0.3909 |
| 80 | 60 | 163 | 1.644 | 1.3 | 1.169 | 0.2214 | 17.5 | 3.9 | 5.104 | 71.63 | 75.38 | 0.3838 |
| 100 | 60 | 163 | 1.644 | 1.3 | 1.015 | 0.1856 | 17.5 | 3.9 | 4.75 | 52.76 | 55.51 | 0.3779 |

??? note "Sources"

    - `Carlsson2012`: Carlsson M, Andersson R, Bloch KM, Steding-Ehrenborg K, Mosén H, Stahlberg F, Ekmehag B, Arheden H. Cardiac output and cardiac index measured with cardiovascular magnetic resonance in healthy subjects, elite athletes and patients with congestive heart failure. J Cardiovasc Magn Reson 14(1):51, 2012. [doi:10.1186/1532-429X-14-51](https://doi.org/10.1186/1532-429X-14-51)
    - `Delanaye2022`: Delanaye P, Gaillard F, van der Weijden J, Mjøen G, Ferhman-Ekholm I, Dubourg L, Ebert N, Schaeffner E, Åkerfeldt T, Goffin K, et al. Age-adapted percentiles of measured glomerular filtration in healthy individuals: extrapolation to living kidney donors over 65 years. Clin Chem Lab Med 60(3):401-407, 2022. [doi:10.1515/cclm-2021-1011](https://doi.org/10.1515/cclm-2021-1011)
    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `Imamura2026`: Imamura H, Ito R, Ichida H, Yoshioka R, Mise Y, Nojiri S, Sano K, Yamataka A, Saiura A. Age-range-matched external validation of standard liver volume equations: methodological re-evaluation of 13 regression models. PLOS One 21:e0351983, 2026. [doi:10.1371/journal.pone.0351983](https://doi.org/10.1371/journal.pone.0351983)
    - `Jopling2009`: Jopling J, Henry E, Wiedmeier SE, Christensen RD. Reference ranges for hematocrit and blood hemoglobin concentration during the neonatal period: data from a multihospital health care system. Pediatrics 123(2):e333-e337, 2009. [doi:10.1542/peds.2008-2654](https://doi.org/10.1542/peds.2008-2654)
    - `NHANES2013_2018`: Centers for Disease Control and Prevention, National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES) 2013-2014, 2015-2016 and 2017-2018, demographics (DEMO), body measures (BMX), complete blood count (CBC) and standard biochemistry profile (BIOPRO) of the cycles H, I and J (public domain). https://wwwn.cdc.gov/nchs/nhanes/ (accessed 2026-11-06). Read from the public files by scripts/database/nhanes_populations.py, scripts/database/nhanes_blood.py and scripts/database/nhanes_waist.py.
    - `Rhodin2009`: Rhodin MM, Anderson BJ, Peters AM, Coulthard MG, Wilkins B, Cole M, Chatelut E, Grubb A, Veal GJ, Keir MJ, Holford NHG. Human renal function maturation: a quantitative description using weight and postmenstrual age. Pediatr Nephrol 24(1):67-76, 2009. [doi:10.1007/s00467-008-0997-5](https://doi.org/10.1007/s00467-008-0997-5)
    - `Scholkmann2019`: Scholkmann F, Ostojic D, Isler H, Bassler D, Wolf M, Karen T. Reference ranges for hemoglobin and hematocrit levels in neonates as a function of gestational age (22-42 weeks) and postnatal age (0-29 days): mathematical modeling. Children 6(3):38, 2019. [doi:10.3390/children6030038](https://doi.org/10.3390/children6030038)
    - `WangX2014`: Wang X, Vrtiska TJ, Avula RT, Walters LR, Chakkera HA, Kremers WK, Lerman LO, Rule AD. Age, kidney function, and risk factors associate differently with cortical and medullary volumes of the kidney. Kidney Int 85(3):677-685, 2014. [doi:10.1038/ki.2013.359](https://doi.org/10.1038/ki.2013.359)
<!-- end generated -->

### Factors of the age

The factors relative to 35 years: of the GFR (at a fixed body surface area), the hematocrit, the fat-free mass (at a fixed body weight and height), the binding proteins, and the brain share of the cardiac output with the factors of the flow fractions of the liver group (hepatic artery and portal organs), the kidney and the other organs.

<!-- generated: age-factors-male (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/age_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| age [yr] | fage_GFR | fage_HCT | fage_FFM | fage_albumin | fage_agp | brain share of CO | liver group | kidney | other organs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| newborn | 0.3547 | 1.162 | 0.8806 | 0.7695 | 0.2716 | 0.3608 | 0.7345 | 0.7029 | 0.7512 |
| 1 | 0.9038 | 0.8017 | 0.8806 | 0.9214 | 0.8937 | 0.5411 | 0.5272 | 0.5045 | 0.5392 |
| 5 | 0.9981 | 0.8258 | 0.8806 | 0.9628 | 0.9616 | 0.3117 | 0.7908 | 0.7567 | 0.8088 |
| 10 | 0.9998 | 0.8698 | 0.8835 | 0.9806 | 0.9793 | 0.2074 | 0.9106 | 0.8714 | 0.9313 |
| 15 | 0.9999 | 0.9682 | 0.9776 | 1.021 | 0.9874 | 0.1533 | 0.9729 | 0.931 | 0.9949 |
| 20 | 1 | 1.007 | 1 | 1.025 | 0.9922 | 0.1241 | 1.006 | 0.9631 | 1.029 |
| 35 | 1 | 1 | 1 | 1 | 1 | 0.1297 | 1 | 1 | 1 |
| 50 | 0.9169 | 0.991 | 1 | 0.9706 | 1.004 | 0.1358 | 0.993 | 0.915 | 1.034 |
| 65 | 0.8152 | 0.9771 | 1 | 0.9519 | 1.007 | 0.1425 | 0.904 | 0.7579 | 1.163 |
| 80 | 0.7136 | 0.9372 | 1 | 0.9343 | 1.008 | 0.1499 | 0.8156 | 0.5872 | 1.296 |
| 100 | 0.578 | 0.9106 | 1 | 0.9232 | 1.01 | 0.1611 | 0.6988 | 0.4682 | 1.413 |

??? note "Sources"

    - `AlSallami2015`: Al-Sallami HS, Goulding A, Grant A, Taylor R, Holford N, Duffull SB. Prediction of fat-free mass in children. Clin Pharmacokinet 54(11):1169-1178, 2015. [doi:10.1007/s40262-015-0277-z](https://doi.org/10.1007/s40262-015-0277-z)
    - `Carlsson2012`: Carlsson M, Andersson R, Bloch KM, Steding-Ehrenborg K, Mosén H, Stahlberg F, Ekmehag B, Arheden H. Cardiac output and cardiac index measured with cardiovascular magnetic resonance in healthy subjects, elite athletes and patients with congestive heart failure. J Cardiovasc Magn Reson 14(1):51, 2012. [doi:10.1186/1532-429X-14-51](https://doi.org/10.1186/1532-429X-14-51)
    - `Chiron1992`: Chiron C, Raynaud C, Maziere B, Zilbovicius M, Laflamme L, Masure MC, Dulac O, Bourguignon M, Syrota A. Changes in regional cerebral blood flow during brain maturation in children and adolescents. J Nucl Med 33(5):696-703, 1992. PMID 1569478.
    - `DaviesShock1950`: Davies DF, Shock NW. Age changes in glomerular filtration rate, effective renal plasma flow, and tubular excretory capacity in adult males. J Clin Invest 29(5):496-507, 1950. [doi:10.1172/JCI102286](https://doi.org/10.1172/JCI102286)
    - `Delanaye2022`: Delanaye P, Gaillard F, van der Weijden J, Mjøen G, Ferhman-Ekholm I, Dubourg L, Ebert N, Schaeffner E, Åkerfeldt T, Goffin K, et al. Age-adapted percentiles of measured glomerular filtration in healthy individuals: extrapolation to living kidney donors over 65 years. Clin Chem Lab Med 60(3):401-407, 2022. [doi:10.1515/cclm-2021-1011](https://doi.org/10.1515/cclm-2021-1011)
    - `Fomon1982`: Fomon SJ, Haschke F, Ziegler EE, Nelson SE. Body composition of reference children from birth to age 10 years. Am J Clin Nutr 35(5):1169-1175, 1982. [doi:10.1093/ajcn/35.5.1169](https://doi.org/10.1093/ajcn/35.5.1169)
    - `Johnson2006`: Johnson TN, Rostami-Hodjegan A, Tucker GT. Prediction of the clearance of eleven drugs and associated variability in neonates, infants and children. Clin Pharmacokinet 45(9):931-956, 2006. [doi:10.2165/00003088-200645090-00005](https://doi.org/10.2165/00003088-200645090-00005)
    - `Jopling2009`: Jopling J, Henry E, Wiedmeier SE, Christensen RD. Reference ranges for hematocrit and blood hemoglobin concentration during the neonatal period: data from a multihospital health care system. Pediatrics 123(2):e333-e337, 2009. [doi:10.1542/peds.2008-2654](https://doi.org/10.1542/peds.2008-2654)
    - `Maharaj2018`: Maharaj AR, Gonzalez D, Cohen-Wolkowiez M, Hornik CP, Edginton AN. Improving pediatric protein binding estimates: an evaluation of alpha1-acid glycoprotein maturation in healthy and infected subjects. Clin Pharmacokinet 57(5):577-589, 2018. [doi:10.1007/s40262-017-0576-7](https://doi.org/10.1007/s40262-017-0576-7)
    - `NHANES2013_2018`: Centers for Disease Control and Prevention, National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES) 2013-2014, 2015-2016 and 2017-2018, demographics (DEMO), body measures (BMX), complete blood count (CBC) and standard biochemistry profile (BIOPRO) of the cycles H, I and J (public domain). https://wwwn.cdc.gov/nchs/nhanes/ (accessed 2026-11-06). Read from the public files by scripts/database/nhanes_populations.py, scripts/database/nhanes_blood.py and scripts/database/nhanes_waist.py.
    - `Rhodin2009`: Rhodin MM, Anderson BJ, Peters AM, Coulthard MG, Wilkins B, Cole M, Chatelut E, Grubb A, Veal GJ, Keir MJ, Holford NHG. Human renal function maturation: a quantitative description using weight and postmenstrual age. Pediatr Nephrol 24(1):67-76, 2009. [doi:10.1007/s00467-008-0997-5](https://doi.org/10.1007/s00467-008-0997-5)
    - `Scholkmann2019`: Scholkmann F, Ostojic D, Isler H, Bassler D, Wolf M, Karen T. Reference ranges for hemoglobin and hematocrit levels in neonates as a function of gestational age (22-42 weeks) and postnatal age (0-29 days): mathematical modeling. Children 6(3):38, 2019. [doi:10.3390/children6030038](https://doi.org/10.3390/children6030038)
    - `Zoli1999`: Zoli M, Magalotti D, Bianchi G, Gueli C, Orlandini C, Grimaldi M, Marchesini G. Total and functional hepatic blood flow decrease in parallel with ageing. Age Ageing 28(1):29-33, 1999. [doi:10.1093/ageing/28.1.29](https://doi.org/10.1093/ageing/28.1.29)
<!-- end generated -->

<!-- generated: age-factors-female (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/age_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| age [yr] | fage_GFR | fage_HCT | fage_FFM | fage_albumin | fage_agp | brain share of CO | liver group | kidney | other organs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| newborn | 0.3547 | 1.322 | 1.081 | 0.8164 | 0.2716 | 0.3673 | 0.7275 | 0.6962 | 0.7425 |
| 1 | 0.9038 | 0.9188 | 1.07 | 0.9775 | 0.8937 | 0.551 | 0.5163 | 0.4941 | 0.527 |
| 5 | 0.9981 | 0.949 | 1.038 | 1.021 | 0.9616 | 0.2859 | 0.8212 | 0.7858 | 0.8381 |
| 10 | 0.9998 | 0.989 | 1.018 | 1.04 | 0.9793 | 0.1841 | 0.9383 | 0.8979 | 0.9576 |
| 15 | 0.9999 | 0.9982 | 1.007 | 1.053 | 0.9874 | 0.1613 | 0.9645 | 0.923 | 0.9844 |
| 20 | 1 | 1.008 | 1 | 1.013 | 0.9922 | 0.1569 | 0.9696 | 0.9278 | 0.9896 |
| 35 | 1 | 1 | 1 | 1 | 1 | 0.1639 | 0.9615 | 0.9615 | 0.9615 |
| 50 | 0.9169 | 1.014 | 1 | 0.9973 | 1.004 | 0.1661 | 0.9589 | 0.8836 | 0.9949 |
| 65 | 0.787 | 1.029 | 1 | 0.9976 | 1.007 | 0.1571 | 0.8894 | 0.7457 | 1.138 |
| 80 | 0.6572 | 1.01 | 1 | 0.9837 | 1.008 | 0.1652 | 0.8016 | 0.5771 | 1.266 |
| 100 | 0.484 | 0.9945 | 1 | 0.9756 | 1.01 | 0.1775 | 0.6857 | 0.4594 | 1.38 |

??? note "Sources"

    - `AlSallami2015`: Al-Sallami HS, Goulding A, Grant A, Taylor R, Holford N, Duffull SB. Prediction of fat-free mass in children. Clin Pharmacokinet 54(11):1169-1178, 2015. [doi:10.1007/s40262-015-0277-z](https://doi.org/10.1007/s40262-015-0277-z)
    - `Carlsson2012`: Carlsson M, Andersson R, Bloch KM, Steding-Ehrenborg K, Mosén H, Stahlberg F, Ekmehag B, Arheden H. Cardiac output and cardiac index measured with cardiovascular magnetic resonance in healthy subjects, elite athletes and patients with congestive heart failure. J Cardiovasc Magn Reson 14(1):51, 2012. [doi:10.1186/1532-429X-14-51](https://doi.org/10.1186/1532-429X-14-51)
    - `Chiron1992`: Chiron C, Raynaud C, Maziere B, Zilbovicius M, Laflamme L, Masure MC, Dulac O, Bourguignon M, Syrota A. Changes in regional cerebral blood flow during brain maturation in children and adolescents. J Nucl Med 33(5):696-703, 1992. PMID 1569478.
    - `DaviesShock1950`: Davies DF, Shock NW. Age changes in glomerular filtration rate, effective renal plasma flow, and tubular excretory capacity in adult males. J Clin Invest 29(5):496-507, 1950. [doi:10.1172/JCI102286](https://doi.org/10.1172/JCI102286)
    - `Delanaye2022`: Delanaye P, Gaillard F, van der Weijden J, Mjøen G, Ferhman-Ekholm I, Dubourg L, Ebert N, Schaeffner E, Åkerfeldt T, Goffin K, et al. Age-adapted percentiles of measured glomerular filtration in healthy individuals: extrapolation to living kidney donors over 65 years. Clin Chem Lab Med 60(3):401-407, 2022. [doi:10.1515/cclm-2021-1011](https://doi.org/10.1515/cclm-2021-1011)
    - `Fomon1982`: Fomon SJ, Haschke F, Ziegler EE, Nelson SE. Body composition of reference children from birth to age 10 years. Am J Clin Nutr 35(5):1169-1175, 1982. [doi:10.1093/ajcn/35.5.1169](https://doi.org/10.1093/ajcn/35.5.1169)
    - `Johnson2006`: Johnson TN, Rostami-Hodjegan A, Tucker GT. Prediction of the clearance of eleven drugs and associated variability in neonates, infants and children. Clin Pharmacokinet 45(9):931-956, 2006. [doi:10.2165/00003088-200645090-00005](https://doi.org/10.2165/00003088-200645090-00005)
    - `Jopling2009`: Jopling J, Henry E, Wiedmeier SE, Christensen RD. Reference ranges for hematocrit and blood hemoglobin concentration during the neonatal period: data from a multihospital health care system. Pediatrics 123(2):e333-e337, 2009. [doi:10.1542/peds.2008-2654](https://doi.org/10.1542/peds.2008-2654)
    - `Maharaj2018`: Maharaj AR, Gonzalez D, Cohen-Wolkowiez M, Hornik CP, Edginton AN. Improving pediatric protein binding estimates: an evaluation of alpha1-acid glycoprotein maturation in healthy and infected subjects. Clin Pharmacokinet 57(5):577-589, 2018. [doi:10.1007/s40262-017-0576-7](https://doi.org/10.1007/s40262-017-0576-7)
    - `NHANES2013_2018`: Centers for Disease Control and Prevention, National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES) 2013-2014, 2015-2016 and 2017-2018, demographics (DEMO), body measures (BMX), complete blood count (CBC) and standard biochemistry profile (BIOPRO) of the cycles H, I and J (public domain). https://wwwn.cdc.gov/nchs/nhanes/ (accessed 2026-11-06). Read from the public files by scripts/database/nhanes_populations.py, scripts/database/nhanes_blood.py and scripts/database/nhanes_waist.py.
    - `Rhodin2009`: Rhodin MM, Anderson BJ, Peters AM, Coulthard MG, Wilkins B, Cole M, Chatelut E, Grubb A, Veal GJ, Keir MJ, Holford NHG. Human renal function maturation: a quantitative description using weight and postmenstrual age. Pediatr Nephrol 24(1):67-76, 2009. [doi:10.1007/s00467-008-0997-5](https://doi.org/10.1007/s00467-008-0997-5)
    - `Scholkmann2019`: Scholkmann F, Ostojic D, Isler H, Bassler D, Wolf M, Karen T. Reference ranges for hemoglobin and hematocrit levels in neonates as a function of gestational age (22-42 weeks) and postnatal age (0-29 days): mathematical modeling. Children 6(3):38, 2019. [doi:10.3390/children6030038](https://doi.org/10.3390/children6030038)
    - `Zoli1999`: Zoli M, Magalotti D, Bianchi G, Gueli C, Orlandini C, Grimaldi M, Marchesini G. Total and functional hepatic blood flow decrease in parallel with ageing. Age Ageing 28(1):29-33, 1999. [doi:10.1093/ageing/28.1.29](https://doi.org/10.1093/ageing/28.1.29)
<!-- end generated -->

### Ontogeny of the enzymes and transporters

<!-- generated: ontogeny (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/ontogeny.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| organ | protein | newborn | 1 month | 6 months | 1 y | 2 y | 5 y | 18 y | sources |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `ki` kidney | ABCB1 | 0.2283 | 0.3215 | 0.7454 | 0.9304 | 0.9904 | 0.9997 | 1 | `Hunt2024` |
| `ki` kidney | SLC22A2 | 0.1962 | 0.3739 | 0.9499 | 0.9971 | 0.9999 | 1 | 1 | `Hunt2024` |
| `ki` kidney | SLC22A6 | 0.1171 | 0.2251 | 0.8577 | 0.987 | 0.9996 | 1 | 1 | `Hunt2024` |
| `ki` kidney | SLC22A8 | 0.1146 | 0.193 | 0.7196 | 0.9489 | 0.9963 | 1 | 1 | `Hunt2024` |
| `li` liver | ABCB1 | 0.2952 | 0.3424 | 0.4572 | 0.5383 | 0.6385 | 0.781 | 0.9442 | `Prasad2016` |
| `li` liver | CYP1A2 | 0.048 | 0.3813 | 1.151 | 1.394 | 1.516 | 1.467 | 1.04 | `Upreti2016`, `PKSimOntogeny` |
| `li` liver | CYP2C19 | 0.7567 | 0.8913 | 0.9971 | 0.9999 | 1 | 1 | 1 | `PKSimOntogeny` |
| `li` liver | CYP2C8 | 0.3598 | 0.5559 | 0.9651 | 0.9973 | 0.9999 | 1 | 1 | `PKSimOntogeny` |
| `li` liver | CYP2C9 | 0.6647 | 0.8211 | 0.9916 | 0.9994 | 1 | 1 | 1 | `PKSimOntogeny` |
| `li` liver | CYP2D6 | 0.5789 | 0.7057 | 0.9537 | 0.992 | 0.9993 | 1 | 1 | `PKSimOntogeny` |
| `li` liver | CYP3A4 | 0.207 | 0.9353 | 1.593 | 1.687 | 1.677 | 1.557 | 1.067 | `Upreti2016`, `PKSimOntogeny` |
| `li` liver | SLC22A1 | 0.1293 | 0.2794 | 0.5855 | 0.7211 | 0.8313 | 0.926 | 0.9865 | `Prasad2016` |
| `li` liver | SLCO1B3 | 0.3049 | 0.3049 | 0.532 | 0.9543 | 0.9983 | 1 | 1 | `Prasad2016` |
| `li` liver | UGT1A1 | 0.007234 | 0.05793 | 0.9958 | 1 | 1 | 1 | 1 | `PKSimOntogeny` |
| `li` liver | UGT1A9 | 0.1471 | 0.1803 | 0.3609 | 0.5535 | 0.7817 | 0.9535 | 0.9976 | `PKSimOntogeny` |
| `li` liver | UGT2B7 | 0.01995 | 0.03846 | 0.3523 | 0.8275 | 0.989 | 0.9999 | 1 | `PKSimOntogeny` |
| `si` small intestine | CYP3A4 | 0.3241 | 0.3529 | 0.4735 | 0.5778 | 0.7083 | 0.8656 | 0.9783 | `PKSimOntogeny` |

??? note "Sources"

    - `Hunt2024`: Hunt JP, Dubinsky S, McKnite AM, Cheung KWK, van Groen BD, Giacomini KM, de Wildt SN, Edginton AN, Watt KM. Maximum likelihood estimation of renal transporter ontogeny profiles for pediatric PBPK modeling. CPT Pharmacometrics Syst Pharmacol 13(4):576-588, 2024. [doi:10.1002/psp4.13102](https://doi.org/10.1002/psp4.13102)
    - `PKSimOntogeny`: Open Systems Pharmacology. PK-Sim Ontogeny Database. OSPSuite.Documentation, https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation, file PK-Sim Ontogeny Database.pdf (PDF of 2026-02-05, accessed 2026-11-06), GPLv2.
    - `Prasad2016`: Prasad B, Gaedigk A, Vrana M, Gaedigk R, Leeder JS, Salphati L, Chu X, Xiao G, Hop C, Evers R, Gan L, Unadkat JD. Ontogeny of hepatic drug transporters as quantified by LC-MS/MS proteomics. Clin Pharmacol Ther 100(4):362-370, 2016. [doi:10.1002/cpt.409](https://doi.org/10.1002/cpt.409)
    - `Upreti2016`: Upreti VV, Wahlstrom JL. Meta-analysis of hepatic cytochrome P450 ontogeny to underwrite the prediction of pediatric pharmacokinetics using physiologically based pharmacokinetic modeling. J Clin Pharmacol 56(3):266-283, 2016. [doi:10.1002/jcph.585](https://doi.org/10.1002/jcph.585)
<!-- end generated -->

### Limitations

A build over age (`PBPKDefinition(age=True)`) writes every reference value and every factor as a piecewise function of `AGE`. roadrunner inlines an initial assignment into each of its readers, so these curves are repeated through the chain of the physiology (the scaling of the organs, the flows, the volumes): the model grows from about 60,000 to 820,000 nodes and its compile from 7 to 66 seconds at LLJIT O0. The cache of the initial value resolvers of roadrunner (sys-bio/roadrunner#1435) generates the code of every initial assignment once and brings the build over age back to about 5 seconds; until a release carries it, a build over age is slow to load and the tests which build one are skipped (`over_age` of `tests/pbpk/conftest.py`). The paediatric fat-free mass is biased: against the reference children of Fomon 1982, whose composition comes from chemical analysis and balance data, the chain of Janmahasatian 2005 and `fage_FFM` (Al-Sallami 2015) gives 0.97 to 0.92 of their fat-free mass in boys and 0.94 to 0.91 in girls over 3 to 10 years, so a child's body fat comes out at 18 to 20 % (boys) and 23 to 27 % (girls) against their 13 to 18 % and 16 to 19 %, and the adipose compartment, which is the remainder of the body weight, inherits the difference. The maturation of the GFR is the parameterisation Rhodin 2009 estimated with a `(WT/70)^0.75` size term while the library scales the GFR with the body surface area; their Table 2 also gives the parameters of a body surface area size model (TM50 52.7 instead of 47.7 weeks of postmenstrual age, Hill 3.30 instead of 3.40, adult GFR 101 instead of 121 ml/min), which would be the self-consistent pairing and would make a term newborn 29 instead of 35 % of the adult; the library keeps their headline column, a decision of the maintainer, because they judge the body surface area scaling clearly worse and compute it by the Boyd formula from weight alone. The decline of the hepatic blood flow with age is applied to the hepatic artery and the portal organs together, while Zoli 1999 measure it in the portal limb alone (portal flow 1209 to 749 ml/min from below 45 to above 75 years, p = 0.0018; hepatic artery flow unchanged, p = 0.36). The interstitial fractions, the densities and the blood contents of the organs and the lumen of the intestine are those of the adults at every age (Fomon 1982 give the extracellular water falling from 42.5 % of the body weight at birth to 26.7 % at 10 years); the preterm neonate needs the gestational age, which the library does not model (`AGE` is the postnatal age of a birth at term). Beyond 85 years every curve is extrapolated from the adults of the data. The kidney and the liver of an infant get the share of the cardiac output which the brain leaves, which no organ data of infants confirm (the ICRP cardiac output of 1.2 l/min at 1 year is low against echocardiography, 1.5 l/min). The decline of the cardiac output with age is debated (1.0 to 0.64 at 80 years in the literature). The organs without data of their own over age (muscle, skin, adipose, bone, heart, lung, gonads) close the sum of the flows: their share of the cardiac output rises by 29 % to 80 years as the flows of the liver and the kidney fall faster than the cardiac output, 12 % in absolute flow, where the resting flow of the leg rather falls with age. The ontogeny of the hepatic CYP3A4 and CYP1A2 is the in vivo one of Upreti 2016, digitized from their Figure 1, and not the in vitro Hill of the PK-Sim database which the other proteins use: both pass the adult activity per amount of liver in the first year, peak at 1.69 and 1.53 times the adult at about 1 and 2 to 3 years and fall back to the adult value at 20 years. Upreti 2016 predict paediatric clearances within twofold in 100 % of their intravenous and 98 % of their oral cases with these functions against 51 % and 67 % with the in vitro ones, which is the reason for the choice; the in vitro functions disagree with them by an order of magnitude (half of the adult activity at 0.63 and 0.49 years in the PK-Sim database, 6.3 and 5.2 years in Upreti's own in vitro fits) and never exceed the adult. Their equations could not be used as printed and the curves were digitized instead, for the reasons in the notes of the two entries. The other ontogenies are in vitro, so a model whose drug is cleared by CYP3A4 and by another enzyme mixes the two kinds. The microsomal protein per gram of liver, which falls by a quarter from 35 to 70 years (Barter 2008), is not applied: the abundance per protein does not change with age and the smaller liver is the change of the elderly.

## Sex

The sex is a choice of the build: the female build has the reference values of the female adult and the female children of ICRP 89 (the children are alike up to 10 years but the brain and the gonads), the female formula of the fat-free mass and the curves of the women where they differ (the GFR of the elderly, the hematocrit and the albumin of NHANES relative to the women of 35 years). A study of both sexes uses a model per sex; a [virtual population](populations.md) assigns each individual to the build of its sex.

What differs beyond the reference anatomy between a woman and a man of the same age, body weight and height are the factors of the sex (`pkpdsim.pbpk.sex_differences`, `database/tables/sex_human.toml`): each is a function of `AGE`, an initial assignment `fsex_<key>` of the female build, as the factors of the age; the male build has none, its physiology is the reference of the parameters of a substance (fu, a clearance), which are measured in men or in groups of both sexes. The differences arise at puberty, so every factor is 1 in the children; some end at the menopause. Where the shape over age has no data of its own (the rise over puberty, the fall over the menopause) it is linear between the ages of the table, a choice.

- **Cardiac output** (`fsex_CO`, 0.914 of the adults). The cardiac index of women and men of the same size is the same: Williams 1994 (the source of the cardiac output of ICRP 89, Table 5, the same reference cardiac index of both sexes at every age; the 6.5 and 5.9 l/min of the adults are the different body sizes), Carlsson 2012 (CMR, 3.1 and 3.2 l/min/m2, n.s.), Patel 2021 (WASE, 18-40 years, equal), de Simone 1997 (no sex effect at the same body surface area) and Middlemiss 2019 (inert gas rebreathing, 0.94-0.96). The scaling of the cardiac output with the lean tissue of the sex gives a woman about 9 % more than a man of the same size (women have 20-25 % more cardiac output per fat-free mass and less fat-free mass at the same size); the factor removes it: computed from the physiology at every age (`same_cardiac_index`, the male over the female cardiac output of the scaling at the size between the two reference individuals of the age), 1 in the children, over the drift of the references of ICRP 89 through puberty to 0.914 from 20 years. The cardiac output of the sexes of the same size and age then agrees within 1.5 % up to a body mass index of about 30 (tested), and every organ flow of the women follows.
- **Brain flow** (`fsex_Qbr`, 1.257 of the brain fraction at 18-45 years, 1.096 from 65 years). The cerebral blood flow per mass of young women is higher: PET 15O-water +17 % at 20 years and +2 % at 65 years (Aanerud 2017), pCASL +16.5 % (Alisch 2021), Xe-133 +11 % (Rodriguez 1988), the total flow by phase-contrast MRI 21 % of the cardiac output in women and 15 % in men of 21-45 years (Xing 2017); no difference before puberty (Satterthwaite 2014, Wu 2016) and converged at about 65 years (Smith 2019). The factor gives the perfusion per mass of the female reference individual 1.17 times the male at 18-45 years and 1.02 from 65 years; the other organs share the rest of the cardiac output in their proportions.
- **Albumin** (`fsex_albumin`, 0.9398). The constant ratio of the women to the men of NHANES at 35 years: times the curve of the women relative to 35 years it gives the ratio of NHANES at every age, 0.93-0.95 at 18-49 years, 0.98-0.99 from 55 years and 0.997 in the children (Weaving 2016 alike, the oral contraceptives not the cause). It multiplies the binding factor of an albumin-bound substance, so fu of a highly bound substance of a woman of 35 years is about 6 % higher than of a man.
- **Gastric emptying of a meal** (`fsex_t_stomach_fed`, 1.2 from 15 to 45 years, 1 from 55 years). The emptying of solids and caloric liquids of premenopausal women is slower (Camilleri 2012, n = 319, 1.16-1.19; Vasavid 2014, n = 189, 1.17; the smaller studies 1.42-1.78, rounded by ICRP 89 to 1.33 and 1.40); the factor multiplies the transit time of the stomach in the fed state (`T_st_fed`), the fasted emptying of water does not differ (ICRP 89). It arises at puberty (Karrento 2025, 1.18; none before, Hauser 2016) and ends at the menopause (Hutson 1989, Gryback 2000).
- **Colonic transit** (`fsex_t_colon`, 1.333 from 20 years). ICRP 89: the colonic transit of the adult females is one third greater than of the males (48 and 36 h), Meier 1995 (n = 164) 1.40, Degen 1996 and Metcalf 1987 the men faster; on the transit times of every segment of the large intestine. The small intestine does not differ (ICRP 89, Degen 1996), nor the colon of children and adolescents (Velde 2013).
- **Enzymes** (`fsex_<o>_<protein>`, per amount of tissue, from 1 at 10 years to the adult value at 16 years): CYP3A4 of the liver 1.2 (the protein 2-fold by immunoblot, Wolbold 2003, and 1.18 by LC-MS/MS, Gaither 2025, the P-gp-independent marker 4beta-hydroxycholesterol +21 %; midazolam per kg +16 % oral, +20 % systemic, Hu and Zhao 2010, which 1.2 reproduces with the liver mass and the flows of the female build; no change at the menopause), CYP1A2 of the liver 0.9 (caffeine, Tantcheva-Poór 1999; the protein 0.74), and for a model which names them UGT2B17 0.4 (Bhatt 2018, rising at puberty, 9 to 15 years), UGT2B15 0.75 (Court 2004, Gaither 2025) and CYP2A6 1.2 until the menopause (nicotine, Benowitz 2006).

<!-- generated: sex-factors (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/sex_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| quantity | factor | 5 y | 12 y | 15 y | 18 y | 35 y | 50 y | 60 y | 80 y | sources |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| cardiac output | `fsex_CO` | 1 | 0.9684 | 0.9373 | 0.9199 | 0.9137 | 0.9137 | 0.9137 | 0.9137 | `Williams1994`, `Carlsson2012`, `Patel2021`, `Middlemiss2019`, `deSimone1997` |
| brain fraction of the cardiac output | `fsex_Qbr` | 1 | 1 | 1.128 | 1.257 | 1.257 | 1.217 | 1.136 | 1.096 | `Aanerud2017`, `Alisch2021`, `Rodriguez1988`, `Xing2017`, `Satterthwaite2014`, `Smith2019` |
| plasma albumin | `fsex_albumin` | 0.9398 | 0.9398 | 0.9398 | 0.9398 | 0.9398 | 0.9398 | 0.9398 | 0.9398 | `NHANES2013_2018`, `Weaving2016` |
| transit time of the stomach, fed | `fsex_t_stomach_fed` | 1 | 1.08 | 1.2 | 1.2 | 1.2 | 1.1 | 1 | 1 | `Camilleri2012`, `Vasavid2014`, `ICRP89`, `Degen1996`, `Hutson1989`, `Gryback2000`, `Karrento2025`, `Hauser2016` |
| transit times of the large intestine | `fsex_t_colon` | 1 | 1 | 1 | 1.2 | 1.333 | 1.333 | 1.333 | 1.333 | `ICRP89`, `Metcalf1987`, `Meier1995`, `Degen1996`, `Velde2013` |
| CYP1A2 liver | `fsex_li_cyp1a2` | 1 | 0.9667 | 0.9167 | 0.9 | 0.9 | 0.9 | 0.9 | 0.9 | `TantchevaPoor1999`, `Gaither2025` |
| CYP2A6 liver | `fsex_li_cyp2a6` | 1 | 1.067 | 1.167 | 1.2 | 1.2 | 1.1 | 1 | 1 | `Benowitz2006`, `Gaither2025` |
| CYP3A4 liver | `fsex_li_cyp3a4` | 1 | 1.067 | 1.167 | 1.2 | 1.2 | 1.2 | 1.2 | 1.2 | `HuZhao2010`, `Wolbold2003`, `Gaither2025` |
| UGT2B15 liver | `fsex_li_ugt2b15` | 1 | 0.9167 | 0.7917 | 0.75 | 0.75 | 0.75 | 0.75 | 0.75 | `Court2004`, `Gaither2025` |
| UGT2B17 liver | `fsex_li_ugt2b17` | 1 | 0.7 | 0.4 | 0.4 | 0.4 | 0.4 | 0.4 | 0.4 | `Bhatt2018`, `Gaither2025` |

??? note "Sources"

    - `Aanerud2017`: Aanerud J, Borghammer P, Rodell A, Jónsdottir KY, Gjedde A. Sex differences of human cortical blood flow and energy metabolism. J Cereb Blood Flow Metab 37(7):2433-2440, 2017. [doi:10.1177/0271678X16668536](https://doi.org/10.1177/0271678X16668536)
    - `Alisch2021`: Alisch JSR, Khattar N, Kim RW, Cortina LE, Rejimon AC, Qian W, Ferrucci L, Resnick SM, Spencer RG, Bouhrara M. Sex and age-related differences in cerebral blood flow investigated using pseudo-continuous arterial spin labeling magnetic resonance imaging. Aging (Albany NY) 13(4):4911-4925, 2021. [doi:10.18632/aging.202673](https://doi.org/10.18632/aging.202673)
    - `Benowitz2006`: Benowitz NL, Lessov-Schlaggar CN, Swan GE, Jacob P 3rd. Female sex and oral contraceptive use accelerate nicotine metabolism. Clin Pharmacol Ther 79(5):480-488, 2006. [doi:10.1016/j.clpt.2006.01.008](https://doi.org/10.1016/j.clpt.2006.01.008)
    - `Bhatt2018`: Bhatt DK, Basit A, Zhang H, Gaedigk A, Lee SB, Claw KG, Mehrotra A, Chaudhry AS, Pearce RE, Gaedigk R, Broeckel U, Thornton TA, Nickerson DA, Schuetz EG, Amory JK, Leeder JS, Prasad B. Hepatic abundance and activity of androgen- and drug-metabolizing enzyme UGT2B17 are associated with genotype, age, and sex. Drug Metab Dispos 46(6):888-896, 2018. [doi:10.1124/dmd.118.080952](https://doi.org/10.1124/dmd.118.080952)
    - `Camilleri2012`: Camilleri M, Iturrino J, Bharucha AE, Burton D, Shin A, Jeong ID, Zinsmeister AR. Performance characteristics of scintigraphic measurement of gastric emptying of solids in healthy participants. Neurogastroenterol Motil 24(12):1076-e562, 2012. [doi:10.1111/j.1365-2982.2012.01972.x](https://doi.org/10.1111/j.1365-2982.2012.01972.x)
    - `Carlsson2012`: Carlsson M, Andersson R, Bloch KM, Steding-Ehrenborg K, Mosén H, Stahlberg F, Ekmehag B, Arheden H. Cardiac output and cardiac index measured with cardiovascular magnetic resonance in healthy subjects, elite athletes and patients with congestive heart failure. J Cardiovasc Magn Reson 14(1):51, 2012. [doi:10.1186/1532-429X-14-51](https://doi.org/10.1186/1532-429X-14-51)
    - `Court2004`: Court MH, Hao Q, Krishnaswamy S, Bekaii-Saab T, Al-Rohaimi A, von Moltke LL, Greenblatt DJ. UDP-glucuronosyltransferase (UGT) 2B15 pharmacogenetics: UGT2B15 D85Y genotype and gender are major determinants of oxazepam glucuronidation by human liver. J Pharmacol Exp Ther 310(2):656-665, 2004. [doi:10.1124/jpet.104.067660](https://doi.org/10.1124/jpet.104.067660)
    - `Degen1996`: Degen LP, Phillips SF. Variability of gastrointestinal transit in healthy women and men. Gut 39(2):299-305, 1996. [doi:10.1136/gut.39.2.299](https://doi.org/10.1136/gut.39.2.299)
    - `Gaither2025`: Gaither KA, Yue G, Singh DK, Trudeau J, Ponraj K, Davydova NY, Lazarus P, Davydov DR, Prasad B. Effects of chronic alcohol intake on the composition of the ensemble of drug-metabolizing enzymes and transporters in the human liver. J Xenobiot 15(1):20, 2025. [doi:10.3390/jox15010020](https://doi.org/10.3390/jox15010020)
    - `Gryback2000`: Grybäck P, Hermansson G, Lyrenäs E, Beckman KW, Jacobsson H, Hellström PM. Nationwide standardisation and evaluation of scintigraphic gastric emptying: reference values and comparisons between subgroups in a multicentre trial. Eur J Nucl Med 27(6):647-655, 2000. [doi:10.1007/s002590050558](https://doi.org/10.1007/s002590050558)
    - `Hauser2016`: Hauser B, Roelants M, De Schepper J, Veereman G, Caveliers V, Devreker T, De Greef E, Vandenplas Y. Gastric emptying of solids in children: reference values for the 13C-octanoic acid breath test. Neurogastroenterol Motil 28(10):1480-1487, 2016. [doi:10.1111/nmo.12845](https://doi.org/10.1111/nmo.12845)
    - `HuZhao2010`: Hu ZY, Zhao YS. Sex-dependent differences in cytochrome P450 3A activity as assessed by midazolam disposition in humans: a meta-analysis. Drug Metab Dispos 38(5):817-823, 2010. [doi:10.1124/dmd.109.031328](https://doi.org/10.1124/dmd.109.031328)
    - `Hutson1989`: Hutson WR, Roehrkasse RL, Wald A. Influence of gender and menopause on gastric emptying and motility. Gastroenterology 96(1):11-17, 1989. [doi:10.1016/0016-5085(89)90758-0](https://doi.org/10.1016/0016-5085(89)90758-0)
    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `Karrento2025`: Karrento K, Chumpitazi BP, Zhang L, Simpson P, Shulman RJ. 13C-Spirulina gastric emptying breath test normative values in healthy children: a multicenter study with comparison to symptomatic children. Am J Gastroenterol 121(2):491-500, 2026 (epub 2025). [doi:10.14309/ajg.0000000000003545](https://doi.org/10.14309/ajg.0000000000003545)
    - `Meier1995`: Meier R, Beglinger C, Dederding JP, Meyer-Wyss B, Fumagalli M, Rowedder A, Turberg Y, Brignoli R. Influence of age, gender, hormonal status and smoking habits on colonic transit time. Neurogastroenterol Motil 7(4):235-238, 1995. [doi:10.1111/j.1365-2982.1995.tb00231.x](https://doi.org/10.1111/j.1365-2982.1995.tb00231.x)
    - `Metcalf1987`: Metcalf AM, Phillips SF, Zinsmeister AR, MacCarty RL, Beart RW, Wolff BG. Simplified assessment of segmental colonic transit. Gastroenterology 92(1):40-47, 1987. [doi:10.1016/0016-5085(87)90837-7](https://doi.org/10.1016/0016-5085(87)90837-7)
    - `Middlemiss2019`: Middlemiss JE, Cocks A, Paapstel K, Maki-Petaja KM, Sunita, Wilkinson IB, McEniery CM; ACCT Study Investigators. Evaluation of inert gas rebreathing for determination of cardiac output: influence of age, gender and body size. Hypertens Res 42(6):834-844, 2019. [doi:10.1038/s41440-018-0179-1](https://doi.org/10.1038/s41440-018-0179-1)
    - `NHANES2013_2018`: Centers for Disease Control and Prevention, National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES) 2013-2014, 2015-2016 and 2017-2018, demographics (DEMO), body measures (BMX), complete blood count (CBC) and standard biochemistry profile (BIOPRO) of the cycles H, I and J (public domain). https://wwwn.cdc.gov/nchs/nhanes/ (accessed 2026-11-06). Read from the public files by scripts/database/nhanes_populations.py, scripts/database/nhanes_blood.py and scripts/database/nhanes_waist.py.
    - `Patel2021`: Patel HN, Miyoshi T, Addetia K, Henry MP, Citro R, Daimon M, Gutierrez Fajardo P, Kasliwal RR, Kirkpatrick JN, Monaghan MJ, et al. Normal values of cardiac output and stroke volume according to measurement technique, age, sex, and ethnicity: results of the World Alliance of Societies of Echocardiography study. J Am Soc Echocardiogr 34(10):1077-1085.e1, 2021. [doi:10.1016/j.echo.2021.05.012](https://doi.org/10.1016/j.echo.2021.05.012)
    - `Rodriguez1988`: Rodriguez G, Warkentin S, Risberg J, Rosadini G. Sex differences in regional cerebral blood flow. J Cereb Blood Flow Metab 8(6):783-789, 1988. [doi:10.1038/jcbfm.1988.133](https://doi.org/10.1038/jcbfm.1988.133)
    - `Satterthwaite2014`: Satterthwaite TD, Shinohara RT, Wolf DH, Hopson RD, Elliott MA, Vandekar SN, Ruparel K, Calkins ME, Roalf DR, Gennatas ED, Jackson C, Erus G, Prabhakaran K, Davatzikos C, Detre JA, Hakonarson H, Gur RC, Gur RE. Impact of puberty on the evolution of cerebral perfusion during adolescence. Proc Natl Acad Sci U S A 111(23):8643-8648, 2014. [doi:10.1073/pnas.1400178111](https://doi.org/10.1073/pnas.1400178111)
    - `Smith2019`: Smith LA, Melbourne A, Owen D, Cardoso MJ, Sudre CH, Tillin T, Sokolska M, Atkinson D, Chaturvedi N, Ourselin S, Hughes AD, Barkhof F, Jäger HR. Cortical cerebral blood flow in ageing: effects of haematocrit, sex, ethnicity and diabetes. Eur Radiol 29(10):5549-5558, 2019. [doi:10.1007/s00330-019-06096-w](https://doi.org/10.1007/s00330-019-06096-w)
    - `TantchevaPoor1999`: Tantcheva-Poór I, Zaigler M, Rietbrock S, Fuhr U. Estimation of cytochrome P-450 CYP1A2 activity in 863 healthy Caucasians using a saliva-based caffeine test. Pharmacogenetics 9(2):131-144, 1999. PMID 10376760.
    - `Vasavid2014`: Vasavid P, Chaiwatanarat T, Pusuwan P, Sritara C, Roysri K, Namwongprom S, Kuanrakcharoen P, Premprabha T, Chunlertrith K, Thongsawat S, Sirinthornpunya S, Ovartlarnporn B, Kachintorn U, Leelakusolvong S, Kositchaiwat C, Chakkaphak S, Gonlachanvit S. Normal solid gastric emptying values measured by scintigraphy using Asian-style meal: a multicenter study in healthy volunteers. J Neurogastroenterol Motil 20(3):371-378, 2014. [doi:10.5056/jnm13114](https://doi.org/10.5056/jnm13114)
    - `Velde2013`: Van de Velde S, Notebaert A, Meersschaut V, Herregods N, Van Winckel M, Van Biervliet S. Colon transit time in healthy children and adolescents. Int J Colorectal Dis 28(12):1721-1724, 2013. [doi:10.1007/s00384-013-1750-5](https://doi.org/10.1007/s00384-013-1750-5)
    - `Weaving2016`: Weaving G, Batstone GF, Jones RG. Age and sex variation in serum albumin concentration: an observational study. Ann Clin Biochem 53(1):106-111, 2016. [doi:10.1177/0004563215593561](https://doi.org/10.1177/0004563215593561)
    - `Williams1994`: Williams LR. Reference values for total blood volume and cardiac output in humans. Report ORNL/TM-12814, Oak Ridge National Laboratory, Oak Ridge, TN, 1994. OSTI ID 10186900, https://www.osti.gov/biblio/10186900 (US government report, public).
    - `Wolbold2003`: Wolbold R, Klein K, Burk O, Nüssler AK, Neuhaus P, Eichelbaum M, Schwab M, Zanger UM. Sex is a major determinant of CYP3A4 expression in human liver. Hepatology 38(4):978-988, 2003. [doi:10.1053/jhep.2003.50393](https://doi.org/10.1053/jhep.2003.50393)
    - `Xing2017`: Xing CY, Tarumi T, Liu J, Zhang Y, Turner M, Riley J, Tinajero CD, Yuan LJ, Zhang R. Distribution of cardiac output to the brain across the adult lifespan. J Cereb Blood Flow Metab 37(8):2848-2856, 2017. [doi:10.1177/0271678X16676826](https://doi.org/10.1177/0271678X16676826)
    - `deSimone1997`: de Simone G, Devereux RB, Daniels SR, Mureddu G, Roman MJ, Kimball TR, Greco R, Witt S, Contaldo F. Stroke volume and cardiac output in normotensive children and adults: assessment of relations with body size and impact of overweight. Circulation 95(7):1837-1843, 1997. [doi:10.1161/01.CIR.95.7.1837](https://doi.org/10.1161/01.CIR.95.7.1837)
<!-- end generated -->

```equations title="the factors of the sex (female build)"
CO_base = CO_ref * SV_tis * f_cardiac_function * fsex_CO
fQbr = 0.130406433384047 * fage_Qbr * fsex_Qbr
fage_Qrest = (1 - 0.130406433384047 * fage_Qbr * fsex_Qbr) / 0.869593566615953
fimp_albumin = fcp_albumin * fage_albumin * fsex_albumin
fimp_li_cyp3a4 = fcp_cyp3a4 * (1 + CKD_hepatic * (fckd_cyp3a4 - 1)) * fage_li_cyp3a4 * fsex_li_cyp3a4
T_st_eff = (T_st + FED * (T_st_fed * fsex_t_stomach_fed - T_st)) * fimp_t_stomach
```

### No difference

The searches (PubMed, Europe PMC, Crossref; meta-analyses, large cohorts and measurements with a reference method first) found no difference of the sexes at the same size, or a difference too small or inconsistent for a factor, for:

- the fasted gastric emptying of water, the small intestinal transit, the fasted gastric and duodenal pH (median; Dressman 1990, Lindahl 1997), the gallbladder emptying; the fluid volumes of the small intestine and the colon were never compared by sex;
- the GFR per body surface area at 20-50 years (Pottel 2017, Fenton 2018, Daugirdas 2009; the faster decline of the women after 50 years is in the aging of the GFR), the renal plasma flow per body surface area (Berg 2006, Miller 1999), the hepatic plasma flow per body surface area (Martin 1975), the resting flow of the leg (Holland 1998);
- the blood volume at the same size (Williams 1994 and Nadler 1962 within 0-3 % of the scaling of ICRP 89), the hematocrit (NHANES per sex, already in);
- the alpha-1-acid glycoprotein (Ritchie 2000, Veering 1990; NHANES: no fall in girls over puberty; the 0.83-0.90 of small cohorts and of Simcyp are not taken), the measured plasma binding of propranolol, disopyramide, diazepam, flecainide, phenytoin and others;
- the proteins of the liver besides those above (CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP3A5, UGT1A1, UGT1A9, UGT2B7, SULT1A1, the transporters OATP1B1, OATP1B3, OATP2B1, OCT1, MATE1, P-gp, BCRP, MRP2, MRP3), every enzyme and transporter of the kidney (Oswald 2019, Prasad 2016, Cheung 2019) and of the intestine (CYP3A4 included, Paine 2005); the P-gp of the liver higher in men (Schuetz 1995) was not replicated (Wolbold 2003, Prasad 2014).

### Limitations

- The factor of the cardiac output removes the difference of the sexes at the same size; within a sex the cardiac output follows the lean tissue (Huisinga 2012), steeper than measured (an elasticity on the body surface area of 1.34 against 0.5-1.15), so the cardiac index of the female reference individual is 0.95 of the male one (measured 0.91-1.03).
- The factor of the brain rests on perfusion per mass (PET and ASL with small groups per decade) and the total flow by phase-contrast MRI; the convergence after the menopause is the least certain.
- The total body water of the women at the same size is 0.86 of the men in the library and 0.82 by Watson 1980: the water fractions of the tissues are of the rat (Rodgers and Rowland), with the same composition of a tissue in both sexes; the liver fat of the women at a normal body mass index is 0.83 of the men (Ulbrich 2015). The distribution has no factor of the sex.
- The colon segments of Metcalf 1987 are of a group of both sexes; the factor of ICRP 89 is applied to them. The fasted gastric volume of women is 6 ml smaller (Roelofs 2024, abstract), not applied.
- Not modelled: the menstrual cycle (no consistent effect on the GFR, CYP3A or the gastric emptying), pregnancy, and the oral contraceptives and the hormone replacement, which are co-medication and not a property of the sex (the oral contraceptives lower AGP to 0.80, CYP1A2 to 0.72 and accelerate the glucuronidation).

## Populations

A population is the reference individual of each sex, which becomes the default individual of the model (`BW`, `HEIGHT`, `AGE`), and the distribution of the height and the body weight of each sex over age (`Anthropometry`), from which the [virtual populations](populations.md) draw. The reference anatomy stays ICRP 89 (with its children) for every population: the organs of an individual of another population follow from the scaling at its body weight, height and age. The population is chosen when the model is built, `PBPKDefinition(population="white_american_nhanes", sex=Sex.FEMALE, ...)`, the European reference individuals of ICRP 89 are the default.

<!-- generated: populations (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/populations/*.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| population | sex | BW [kg] | HEIGHT [cm] | AGE [yr] | distribution over age |
| --- | --- | ---: | ---: | ---: | --- |
| European (ICRP 89) `european_icrp` | male | 73 | 176 | 35 | none |
| European (ICRP 89) `european_icrp` | female | 60 | 163 | 35 | none |
| Black American (NHANES 2013-2018) `black_american_nhanes` | male | 91.26 | 176.9 | 33.93 | 0-100 y, 37 ages; HEIGHT &plusmn; 7.138 cm, ln BW &plusmn; 0.2537, correlation 0.38 at 35 y |
| Black American (NHANES 2013-2018) `black_american_nhanes` | female | 86.21 | 163.5 | 34.91 | 0-100 y, 37 ages; HEIGHT &plusmn; 6.774 cm, ln BW &plusmn; 0.286, correlation 0.38 at 35 y |
| East Asian (Tanaka 1996) `east_asian_tanaka` | male | 60 | 170 | 35 | 20-49 y, 1 ages; HEIGHT &plusmn; 5.68 cm, ln BW &plusmn; 0.1383, independent |
| East Asian (Tanaka 1996) `east_asian_tanaka` | female | 51.28 | 155 | 35 | 20-49 y, 1 ages; HEIGHT &plusmn; 5.22 cm, ln BW &plusmn; 0.1409, independent |
| Mexican American (NHANES 2013-2018) `mexican_american_nhanes` | male | 89.04 | 171.5 | 33.99 | 0-100 y, 37 ages; HEIGHT &plusmn; 7.155 cm, ln BW &plusmn; 0.2148, correlation 0.51 at 35 y |
| Mexican American (NHANES 2013-2018) `mexican_american_nhanes` | female | 78.64 | 157.8 | 34.63 | 0-100 y, 37 ages; HEIGHT &plusmn; 6.7 cm, ln BW &plusmn; 0.2415, correlation 0.42 at 35 y |
| White American (NHANES 2013-2018) `white_american_nhanes` | male | 92.28 | 177.9 | 34.95 | 0-100 y, 37 ages; HEIGHT &plusmn; 7.2 cm, ln BW &plusmn; 0.2206, correlation 0.35 at 35 y |
| White American (NHANES 2013-2018) `white_american_nhanes` | female | 78.08 | 164.3 | 35.41 | 0-100 y, 37 ages; HEIGHT &plusmn; 6.419 cm, ln BW &plusmn; 0.2749, correlation 0.27 at 35 y |

??? note "Sources"

    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `NHANES2013_2018`: Centers for Disease Control and Prevention, National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES) 2013-2014, 2015-2016 and 2017-2018, demographics (DEMO), body measures (BMX), complete blood count (CBC) and standard biochemistry profile (BIOPRO) of the cycles H, I and J (public domain). https://wwwn.cdc.gov/nchs/nhanes/ (accessed 2026-11-06). Read from the public files by scripts/database/nhanes_populations.py, scripts/database/nhanes_blood.py and scripts/database/nhanes_waist.py.
    - `Tanaka1996`: Tanaka G, Kawamura H. Anatomical and physiological characteristics for Asian reference man: male and female of different ages: Tanaka model. NIRS-M-115, National Institute of Radiological Sciences, Hitachinaka, Japan, 1996. https://repo.qst.go.jp/records/73628 (accessed 2026-10-25).
<!-- end generated -->

The NHANES populations are generated from the public files of the examinations 2013-2018 (`scripts/database/nhanes_populations.py`, the women pregnant at the examination excluded): the reference individual the survey-weighted mean of the adults of 20 to 49 years, the distributions at 37 ages from 1.5 months to 80 years and over the survey-weighted kernel moments of the height and of the logarithm of the body weight and their correlation (a Gaussian kernel in age of 1.5 months for the infants, half a year for the children, 2.5 years for the adults; the row of 80 years holds up to 100 years). The East Asian population is the Asian Reference Man of Tanaka 1996 with the distributions of the Japanese of 20 to 49 years of the same tables (no correlation given, height and body weight are drawn independently). PK-Sim also offers a Japanese population of national statistics of 2015; its primary data were not read, so the library has none. ICRP 89 gives reference values, not distributions: the European population has no covariate distributions.

The covariate distributions are data of the library and never SBML `distrib` uncertainties: `run_simulation --uncertainty` draws every parameter with a `distrib` uncertainty and must draw only the uncertainty of the fits.
