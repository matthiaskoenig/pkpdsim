# Physiology

The standard human physiology is a set of SBML parameters of the top model: the covariates of the individual, the reference values of the sex of the build (ICRP 89), and everything which follows from them, the organ masses and volumes, the blood volume, the cardiac output, the blood flows and the GFR. An experiment changes the covariates (`BW`, `HEIGHT`) and the factors of organ function (`f_renal_function`, `f_cardiac_function`, the [impairment](impairment.md) severities) before the simulation, as it changes `BW` in the hand-written models, and the organs follow.

## Topology

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/circulation.svg"
</div>

The circulation of Jones & Rowland-Yeo 2013 and PK-Sim (Kuepfer 2016): an arterial and a venous blood pool, the lung in series (its flow is the cardiac output), the other organs in parallel between the arterial and the venous pool. Stomach, small intestine, large intestine, spleen and pancreas drain into the portal vein, which feeds the liver together with the hepatic artery; the liver drains into the venous pool. The organs (`pkpdsim.pbpk.Organ`, their two-letter tokens are the prefixes of their ids) are lung `lu`, heart `he`, brain `br`, liver `li`, kidney `ki`, stomach `st`, small intestine `si`, large intestine `la`, spleen `sp`, pancreas `pa`, muscle `mu`, skin `sk`, adipose `ad`, bone `bo` and gonads `go`; the blood pools (`Pool`) are arterial `ar`, venous `ve` and portal `po`.

Every organ owns its blood flows: it reads the plasma (and blood cell) species of the pool which supplies it and writes into the pool it drains into, through its ports (see [Standard organ](organ.md)). The pools have no reactions of their own besides the exchange with kinetic blood cells and the intravenous dose. The intestine variants replace small and large intestine with one submodel `GU` with the inflow of both.

## Reference values

The organ masses, densities, blood contents and blood flows are those of the ICRP 89 reference individuals (Valentin 2002), taken directly from its tables; the interstitial fractions are those of Schmitt 2008. The data file `data/physiology_human.toml` compares the male values with Jones & Rowland-Yeo 2013 in the notes of every value. The conventions of the data:

- `mass` is the mass of the organ without its blood (the blood is counted on its own); `density` is the specific gravity (1 kg/l where ICRP 89 states none); `blood` is the blood content of the organ as fraction of the blood volume.
- The organ volume of the model is the volume of its tissue and its blood, `V<o> = M<o>/density + blood*Vblood`; the vascular fraction is its blood, `f_vas_<o> = blood*Vblood/V<o>`, so that the vascular spaces of the organs and the pools are the blood volume for every individual.
- `f_int` is the interstitial fraction of the organ volume, the interstitial fraction of the tissue without blood of Schmitt 2008 times `1 - f_vas`; the intracellular fraction is `1 - f_vas - f_int`.
- `flow` is the blood flow as fraction of the cardiac output (ICRP 89 Table 2.40); the liver states its arterial flow (`Qli` is the hepatic artery, `Qha`), the portal organs feed it in addition. ICRP 89 gives about 8 % of the cardiac output to tissues which are no organ of the model (`[unmapped]` of the data): the flow fractions of the model, `fQ<o>`, are renormalised over the organs, so they sum to the cardiac output.
- The portal vein pool is the anatomical portal vein (Chau 2023, 4.4 ml); the splanchnic venous blood stays in the vascular volume of the portal organs. The blood of the unmapped tissues is part of the venous pool.
- Adipose closes the body weight: its mass is the remainder, `Mad = BW - Mlean - Vblood*rho_blood`; the ICRP adipose mass in the data is a reference value only.

<!-- generated: physiology-general (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/physiology_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

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

### Male reference individual

<!-- generated: physiology-organs-male (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/physiology_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

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

Sources of the columns: mass, density, blood, f_vas, flow: `ICRP89`; f_int: `Schmitt2008`.

??? note "Sources"

    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `Schmitt2008`: Schmitt W. General approach for the calculation of tissue to plasma partition coefficients. Toxicol In Vitro 22(2):457-467, 2008. [doi:10.1016/j.tiv.2007.09.010](https://doi.org/10.1016/j.tiv.2007.09.010)
<!-- end generated -->

### Female reference individual

<!-- generated: physiology-organs-female (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/physiology_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

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

Sources of the columns: mass, density, blood, f_vas, flow: `ICRP89`; f_int: `Schmitt2008`.

??? note "Sources"

    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `Schmitt2008`: Schmitt W. General approach for the calculation of tissue to plasma partition coefficients. Toxicol In Vitro 22(2):457-467, 2008. [doi:10.1016/j.tiv.2007.09.010](https://doi.org/10.1016/j.tiv.2007.09.010)
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

The covariates are the body weight `BW` (kg), the height `HEIGHT` (cm) and the age `AGE` (yr, no rule reads it yet); their defaults are the reference individual of the sex of the population of the build. The chain follows Huisinga 2012: body mass index, body surface area (Du Bois 1916) and fat-free mass (Janmahasatian 2005); brain is constant, skin scales with the body surface area, the other lean organs and the blood with the fat-free mass less brain and skin (`SV_tis`, Huisinga 2012 eqs. 9, 10 and 12), adipose is the remainder of the body weight. The cardiac output scales as the tissues (eq. 13), the organ flows are fractions of it, and the GFR scales with the body surface area. This replaces the `ftissue*BW/BW_ref` scaling of the hand-written models.

```equations title="covariates"
BMI = BW / 1 kg / (HEIGHT / 100 cm)^2 * 1 kg_per_m2
BSA = 0.007184 m2 * (BW / 1 kg)^0.425 * (HEIGHT / 1 cm)^0.725
FFM = BW * 9270 / (6680 + 216 * BMI / 1 kg_per_m2)
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
GFR = GFR_ref * BSA / BSA_ref * f_renal_function * fimp_GFR
```

The factors `fimp_<id>` are those of the [impairment](impairment.md), 1 for the healthy individual; `f_cardiac_function` and `f_renal_function` are the factors of organ function of the hand-written models, with the same ids. All of these are initial assignments on constant parameters: a covariate changed before the simulation (a change of an experiment) changes every organ, a change during the simulation does not, unless the definition names it in `PBPKDefinition.dynamic`. The kidney mass carries the CKD factor of the kidney volume (`fimp_Mki`, 1 until its source is read, see [Impairment](impairment.md#limitations)).

## Populations

A population is the reference individual of each sex, which becomes the default individual of the model (`BW`, `HEIGHT`, `AGE`), and the distributions of the covariates, data for a later sampler of virtual individuals. The reference anatomy stays ICRP 89 for every population: the organs of another population's reference individual follow from the scaling at its body weight and height. The population is chosen when the model is built, `PBPKDefinition(population="white_american_nhanes", sex=Sex.FEMALE, ...)`, the European reference individuals of ICRP 89 are the default.

<!-- generated: populations (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/populations/*.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

| population | sex | BW [kg] | HEIGHT [cm] | AGE [yr] | covariates |
| --- | --- | ---: | ---: | ---: | --- |
| European (ICRP 89) `european_icrp` | male | 73 | 176 | 40 | none |
| European (ICRP 89) `european_icrp` | female | 60 | 163 | 40 | none |
| Black American (NHANES 2013-2018) `black_american_nhanes` | male | 91.22 | 176.9 | 33.47 | BW lognormal 91.22 &plusmn; 25.11<br>HEIGHT normal 176.9 &plusmn; 6.93 |
| Black American (NHANES 2013-2018) `black_american_nhanes` | female | 86.2 | 163.5 | 34.05 | BW lognormal 86.2 &plusmn; 25.02<br>HEIGHT normal 163.5 &plusmn; 6.43 |
| East Asian (Tanaka 1996) `east_asian_tanaka` | male | 60 | 170 | 40 | BW lognormal 63.63 &plusmn; 8.84<br>HEIGHT normal 170 &plusmn; 5.68 |
| East Asian (Tanaka 1996) `east_asian_tanaka` | female | 51.28 | 155 | 40 | BW lognormal 52.34 &plusmn; 7.41<br>HEIGHT normal 155 &plusmn; 5.22 |
| Mexican American (NHANES 2013-2018) `mexican_american_nhanes` | male | 89.4 | 171.5 | 33.47 | BW lognormal 89.4 &plusmn; 20.37<br>HEIGHT normal 171.5 &plusmn; 7.17 |
| Mexican American (NHANES 2013-2018) `mexican_american_nhanes` | female | 78.76 | 157.8 | 33.95 | BW lognormal 78.76 &plusmn; 21.27<br>HEIGHT normal 157.8 &plusmn; 6.51 |
| White American (NHANES 2013-2018) `white_american_nhanes` | male | 92.32 | 178 | 34.56 | BW lognormal 92.32 &plusmn; 22.42<br>HEIGHT normal 178 &plusmn; 6.71 |
| White American (NHANES 2013-2018) `white_american_nhanes` | female | 78.25 | 164.3 | 34.74 | BW lognormal 78.25 &plusmn; 22.36<br>HEIGHT normal 164.3 &plusmn; 6.33 |

??? note "Sources"

    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `NHANES2013_2018`: Centers for Disease Control and Prevention, National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES) 2013-2014, 2015-2016 and 2017-2018, demographics, body measures and laboratory data (public domain). https://wwwn.cdc.gov/nchs/nhanes/ (accessed 2026-10-25). Read as the data set `mecdt` of the R package httk 2.7.4, see `httk`.
    - `Ring2017`: Ring CL, Pearce RG, Setzer RW, Wetmore BA, Wambaugh JF. Identifying populations sensitive to environmental chemicals by simulating toxicokinetic variability. Environ Int 106:105-118, 2017. [doi:10.1016/j.envint.2017.06.004](https://doi.org/10.1016/j.envint.2017.06.004)
    - `Tanaka1996`: Tanaka G, Kawamura H. Anatomical and physiological characteristics for Asian reference man: male and female of different ages: Tanaka model. NIRS-M-115, National Institute of Radiological Sciences, Hitachinaka, Japan, 1996. https://repo.qst.go.jp/records/73628 (accessed 2026-10-25).
<!-- end generated -->

The NHANES populations are the adults of 20 to 49 years of the examinations 2013-2018 as preprocessed by httk (HTTK-Pop, Ring 2017), survey-weighted means and standard deviations; the East Asian population is the Asian Reference Man of Tanaka 1996 with the distributions of the Japanese of the same tables. PK-Sim also offers a Japanese population of national statistics of 2015; its primary data were not read, so the library has none. ICRP 89 gives reference values, not distributions: the European population has no covariate distributions.

The covariate distributions are data of the library and never SBML `distrib` uncertainties: `run_simulation --uncertainty` draws every parameter with a `distrib` uncertainty and must draw only the uncertainty of the fits.
