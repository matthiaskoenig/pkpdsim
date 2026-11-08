# Adipose tissue and obesity

The adipose tissue of the library is the remainder of the body weight after the lean organs and the blood (`Mad = BW - Mlean - Vblood*rho_blood`, [physiology](physiology.md)): the fat-free mass is that of Janmahasatian 2005, every lean organ scales after Huisinga 2012, and what is left is adipose. The masses of an individual therefore always sum to its body weight and most of obesity is already in the model. `PBPKDefinition(adipose=Adipose())` (`pkpdsim.pbpk.adipose`, [#136](https://github.com/matthiaskoenig/pkpdsim/issues/136)) adds what the remainder does not give. It is a choice of the build; a build without it is unchanged.

- **The body composition**: the lean body weight `LBW = BW - Mad`, the body fat mass `FM = BW - FFM` with the body fat fraction `BFF = FM/BW` and the fat mass index `FMI = FM/height^2` (Kelly 2009), the waist circumference `WAIST` with the waist-to-height ratio `WHtR`, and the split of the adipose tissue into the visceral `VAT = f_vat*Mad` and the subcutaneous `SAT = Mad - VAT`. The classes of the body mass index are `pkpdsim.impairment.BMI_CLASSES`, which an experiment scans and draws like the classes of an impairment ([classes](classes.md)).
- **The physiology of overweight and obesity**: the factors `fobes_<key>` on the cardiac output, the flow fractions, the plasma proteins of the binding and the hepatic enzymes, each 1 up to a body mass index of 25 kg/m2. `Adipose(obesity=False)` leaves them out.
- **A change of the body weight over time**: `Adipose(effect=Formula(...))` makes `BW` a state of the model which approaches a set point the effect of the model sets, and the whole physiology follows it.

## Body composition

`FM` and `Mad` are two different quantities and the difference matters. `FM` is the fat of a body composition measurement (DXA, the four-compartment model); `Mad` is the adipose tissue compartment, which exceeds the fat mass by 2 to 4 kg over the whole range of the body mass index, because the neutral lipid is about 85 % of adipose tissue and because the compartment also holds the tissues of ICRP 89 which the library does not map (4.2 kg in the man, 4.3 in the woman). A model which needs the lipid pool uses the tissue composition of adipose ([organ](organ.md), Rodgers and Rowland); a model which compares with a DXA measurement uses `FM`.

The visceral share `f_vat` is the whole-body MRI regression of Shen 2009, the only cohort which measures the total and the visceral adipose tissue of the same subjects, which is what the remainder of the library is. It is a function of the age of the sex and not of the body mass index: the share roughly trebles between 20 and 60 years and is half as large in women.

<!-- generated: adipose-visceral (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/adipose_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| sex | 20 y | 30 y | 40 y | 50 y | 60 y | 70 y | sources |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| male | 0.04568 | 0.07967 | 0.1204 | 0.1612 | 0.1954 | 0.2166 | `Shen2009`, `Linge2018`, `Demerath2007` |
| female | 0.02402 | 0.03763 | 0.06114 | 0.08722 | 0.1086 | 0.1178 | `Shen2009`, `Linge2018`, `Demerath2007` |

??? note "Sources"

    - `Demerath2007`: Demerath EW, Shen W, Lee M, Choh AC, Czerwinski SA, Siervogel RM, Towne B. Approximation of total visceral adipose tissue with a single magnetic resonance image. Am J Clin Nutr 85(2):362-368, 2007. [doi:10.1093/ajcn/85.2.362](https://doi.org/10.1093/ajcn/85.2.362)
    - `Linge2018`: Linge J, Borga M, West J, Tuthill T, Miller MR, Dumitriu A, Thomas EL, Romu T, Tunon P, Bell JD, Dahlqvist Leinhard O. Body composition profiling in the UK Biobank imaging study. Obesity (Silver Spring) 26(11):1785-1795, 2018. [doi:10.1002/oby.22210](https://doi.org/10.1002/oby.22210)
    - `Shen2009`: Shen W, Punyanitya M, Silva AM, Chen J, Gallagher D, Sardinha LB, Allison DB, Heymsfield SB. Sexual dimorphism of adipose tissue distribution across the lifespan: a cross-sectional whole-body magnetic resonance imaging study. Nutr Metab (Lond) 6:17, 2009. [doi:10.1186/1743-7075-6-17](https://doi.org/10.1186/1743-7075-6-17)
<!-- end generated -->

Below 18 and above 70 years the values of 18 and 70 years are kept: the Shen cubic is not monotone in the children, whose visceral fat the regression does not resolve, and it turns over above about 70 years. A build which is not over age takes the share of 35 years. The visceral depot is a mass of the body composition and **not an organ of its own**: the adipose tissue of the library is one well-stirred organ which drains into the venous blood, while the visceral depot drains through the portal vein (about 3 to 10 % of the portal inflow at the perfusion of ICRP 89). Splitting it into two parallel sub-organs, the visceral one draining into the portal pool, is the natural extension and is not implemented.

### Waist circumference

The body mass index says how much adipose tissue an individual carries, the **waist circumference** says where it is, and it is the one covariate of the visceral depot a study is likely to report ([#154](https://github.com/matthiaskoenig/pkpdsim/issues/154)). Every build with the adipose option therefore carries it, measured or predicted:

```
WAIST_base = intercept + bw*BW + height*HEIGHT + age*AGE
WAIST      = f_waist * WAIST_base
f_vat      = f_vat(AGE, sex) * exp(k_vat_waist*(WAIST - WAIST_base))
```

`WAIST_base` is the waist of the population at the sex, body weight, height and age of the individual, and `f_waist` (default 1) is the covariate: the measured waist over that prediction, which `pkpdsim.pbpk.waist_factor(sex, waist, bw, height, age)` computes. The waist is the circumference at the site of the WHO protocol, the **midpoint between the lowest rib and the iliac crest**; the action points of the WHO are 94 and 102 cm in men and 80 and 88 cm in women, and the waist-to-height ratio `WHtR = WAIST/HEIGHT` is the index which needs no classification by sex, with the boundary of 0.5 every guideline shares.

```equations
WAIST = f_waist * WAIST_base
WHtR = WAIST / HEIGHT
```

<!-- generated: adipose-waist (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/adipose_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| sex | intercept [cm] | per kg | per cm | per year | WAIST of the reference individual [cm] | VAT per cm [%] | sources |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| male | 98.83 | 0.7657 | -0.4415 | 0.2056 | 84.21 | 2.958 | `NHANES2013_2018`, `Haueise2023` |
| female | 88.26 | 0.7904 | -0.4009 | 0.1482 | 75.53 | 2.715 | `NHANES2013_2018`, `Haueise2023` |

??? note "Sources"

    - `Haueise2023`: Haueise T, Schick F, Stefan N, Schlett CL, Weiss JB, Nattenmüller J, Göbel-Guéniot K, Norajitra T, et al. Analysis of volume and topography of adipose tissue in the trunk: results of MRI of 11,141 participants in the German National Cohort. Sci Adv 9(19):eadd0433, 2023. [doi:10.1126/sciadv.add0433](https://doi.org/10.1126/sciadv.add0433)
    - `NHANES2013_2018`: Centers for Disease Control and Prevention, National Center for Health Statistics. National Health and Nutrition Examination Survey (NHANES) 2013-2014, 2015-2016 and 2017-2018, demographics (DEMO), body measures (BMX), complete blood count (CBC) and standard biochemistry profile (BIOPRO) of the cycles H, I and J (public domain). https://wwwn.cdc.gov/nchs/nhanes/ (accessed 2026-11-06). Read from the public files by scripts/database/nhanes_populations.py, scripts/database/nhanes_blood.py and scripts/database/nhanes_waist.py.
<!-- end generated -->

The prediction is the survey-weighted least squares of NHANES 2013-2018 (`scripts/database/nhanes_waist.py`, 6985 men and 7291 women of 20 to 79 years): the covariates of a build already determine the waist to 4.6 cm in a man (R2 0.922) and 5.7 cm in a woman (R2 0.893), and the slopes do not depend on the population of the fit (the non-Hispanic white respondents alone give the same numbers to the third digit). Its **level** is not NHANES's: the intercept is shifted by -2.4 cm (men) and -6.2 cm (women) so that the prediction reproduces the mean waist of Haueise 2023, the cohort of the sensitivity below. NHANES measures the waist at the uppermost lateral border of the ilium and the German National Cohort at the WHO site, the two populations differ as well, and their means cannot separate the two causes; the shift is a model choice, and it puts `WAIST` on the site which the guidelines and the trials report.

The sensitivity `k_vat_waist` is the waist term of the multiple regression of Haueise 2023 (whole-trunk MRI of 11,141 participants of the German National Cohort, deep learning segmentation, 20 to 72 years, body mass index 15.6 to 54.6 kg/m2): at a fixed age, height and body weight the visceral depot rises by **3.0 % per cm of waist in men and 2.7 % in women** (their Table S1 model 2, standardized betas 0.734 and 0.576 with R2 0.75, over the standard deviations and means of their Table 2). The other three predictors of that regression are the covariates of the build, so only the waist is a residual over what the library already does, and the rise with age is not counted twice: the beta of the age falls from 0.301 to 0.152 when the waist enters, and the prediction `WAIST_base` carries that rise (0.21 cm per year in men).

An individual of average shape (`f_waist = 1`) therefore keeps the visceral share of its age and sex, and the build of a model which measures no waist is the one described above. A waist 10 cm above the prediction raises the visceral depot of the reference man by 34 % (1.86 to 2.50 kg) and lowers the subcutaneous depot by the same mass; the total adipose tissue, every volume and every flow are untouched, because the waist moves adipose tissue between the two depots and nothing else.

Against the cohort itself the library predicts **4.87 l of visceral adipose tissue for the men of Haueise 2023 against their measured 4.84 l**, and 2.99 against 2.51 l for the women, the same over-prediction in women which the check against Linge 2018 shows.

The body composition of the classes of the body mass index, at the height of the reference individual of each sex:

<!-- generated: adipose-classes (scripts/pbpk_tables.py), do not edit -->
*Generated from `pkpdsim.pbpk.physiology.body_composition` at the height of the reference individual of the sex (ICRP 89, 35 years) by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| class | sex | BMI | BW [kg] | FFM [kg] | FM [kg] | %fat | FMI | Mad [kg] | WAIST [cm] | WHtR | VAT [kg] | SAT [kg] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Underweight | male | 17 | 52.66 | 47.16 | 5.504 | 10.45 | 1.777 | 8.094 | 68.64 | 0.39 | 0.8061 | 7.288 |
| Normal weight | male | 23.57 | 73 | 57.49 | 15.51 | 21.24 | 5.006 | 18.7 | 84.21 | 0.4785 | 1.862 | 16.83 |
| Overweight | male | 27.5 | 85.18 | 62.57 | 22.61 | 26.55 | 7.3 | 26.1 | 93.54 | 0.5315 | 2.599 | 23.5 |
| Obesity class I | male | 32.5 | 100.7 | 68.12 | 32.55 | 32.34 | 10.51 | 36.36 | 105.4 | 0.5989 | 3.621 | 32.74 |
| Obesity class II | male | 37.5 | 116.2 | 72.86 | 43.3 | 37.28 | 13.98 | 47.38 | 117.3 | 0.6663 | 4.719 | 42.66 |
| Obesity class III | male | 45 | 139.4 | 78.79 | 60.6 | 43.48 | 19.56 | 65.02 | 135.1 | 0.7673 | 6.475 | 58.54 |
| Underweight | female | 17 | 45.17 | 32.39 | 12.78 | 28.3 | 4.81 | 14.57 | 63.81 | 0.3915 | 0.7085 | 13.87 |
| Normal weight | female | 22.58 | 60 | 38.92 | 21.08 | 35.13 | 7.933 | 23.26 | 75.53 | 0.4634 | 1.131 | 22.13 |
| Overweight | female | 27.5 | 73.06 | 43.73 | 29.34 | 40.15 | 11.04 | 31.81 | 85.86 | 0.5267 | 1.546 | 30.26 |
| Obesity class I | female | 32.5 | 86.35 | 47.9 | 38.45 | 44.52 | 14.47 | 41.16 | 96.36 | 0.5912 | 2.001 | 39.16 |
| Obesity class II | female | 37.5 | 99.63 | 51.51 | 48.12 | 48.3 | 18.11 | 51.05 | 106.9 | 0.6556 | 2.481 | 48.57 |
| Obesity class III | female | 45 | 119.6 | 56.09 | 63.47 | 53.09 | 23.89 | 66.67 | 122.6 | 0.7522 | 3.241 | 63.42 |
<!-- end generated -->

The WHO cut-offs (18.5, 25, 30, 35, 40 kg/m2) are those of WHO Technical Report Series 894. The WHO Expert Consultation 2004 retained them as the international classification and added the public health action points **23.0, 27.5, 32.5 and 37.5 kg/m2** for Asian populations, because the point of observed risk varies between a body mass index of 22 and 25 across them. The percent body fat of Gallagher 2000 reproduces those points: an Asian man of 22.7 and an Asian woman of 23.5 have the body fat of a white subject at 25, and 27.1 and 28.5 the body fat at 30. The library carries the action points here and not as a classification of its own.

## The physiology of obesity

Every factor is the **residual** of the obese repository of Berton 2022 over what the scaling of the library already does,

```
fobes(q, BMI) = [Berton(q, BMI)/library(q, BMI)] / [Berton(q, 25)/library(q, 25)]
```

at the reference height of the sex and at 35 years, so a factor of 1 means that the fat-free mass of Janmahasatian 2005 and the scaling of Huisinga 2012 reproduce the measured obese value and that a factor on top would count it twice. A factor on a flow fraction is the percentage of the cardiac output of Berton 2022 over its value at a body mass index of 25, which is a residual by construction (the fractions of the library are constant). The factors are 1 up to a body mass index of 25, so a normal weight individual of a build with `adipose` has the physiology of a build without it.

<!-- generated: adipose-factors (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/adipose_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| quantity | factor | on | BMI 25 | BMI 30 | BMI 35 | BMI 40 | BMI 50 | sources |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cardiac output | `fobes_CO` | every organ flow | 1 / 1 | 0.9475 / 0.9417 | 0.9235 / 0.913 | 0.9152 / 0.9006 | 0.9227 / 0.901 | `Berton2022`, `Collis2001`, `deDivitiis1981` |
| adipose flow | `fobes_Qad` | `Qad` | 1 / 1 | 1.077 / 1.062 | 1.161 / 1.128 | 1.25 / 1.198 | 1.451 / 1.352 | `Berton2022`, `Levitt2007`, `Blaak1995`, `Summers1996`, `Virtanen2002`, `Lesser1967` |
| muscle flow | `fobes_Qmu` | `Qmu` | 1 | 0.9418 | 0.8869 | 0.8353 | 0.7408 | `Berton2022` |
| brain flow | `fobes_Qbr` | `Qbr` | 1 | 0.9458 | 0.8917 | 0.8375 | 0.7292 | `Berton2022` |
| skin flow | `fobes_Qsk` | `Qsk` | 1 | 0.9648 | 0.9296 | 0.8944 | 0.824 | `Berton2022` |
| heart flow | `fobes_Qhe` | `Qhe` | 1 / 1 | 1.089 / 1.085 | 1.179 / 1.171 | 1.268 / 1.256 | 1.446 / 1.426 | `Berton2022` |
| gut flow | `fobes_Qst`, `fobes_Qsi`, `fobes_Qla` | `Qst`, `Qsi`, `Qla` | 1 / 1 | 0.9466 / 0.9406 | 0.8963 / 0.8853 | 0.8493 / 0.8341 | 0.7648 / 0.7441 | `Berton2022`, `Ghobadi2011` |
| pancreas flow | `fobes_Qpa` | `Qpa` | 1 / 1 | 0.9521 / 0.9566 | 0.9082 / 0.9168 | 0.8675 / 0.8797 | 0.7935 / 0.8122 | `Berton2022` |
| albumin | `fobes_albumin` | `fimp_albumin` | 1 | 0.9816 | 0.9633 | 0.9449 | 0.9082 | `Berton2022`, `Benedek1983`, `Benedek1984`, `Ghobadi2011` |
| agp | `fobes_agp` | `fimp_agp` | 1 | 1.077 | 1.154 | 1.231 | 1.385 | `Berton2022`, `Benedek1983`, `Benedek1984`, `Ghobadi2011` |
| CYP2E1 liver | `fobes_li_cyp2e1` | `fimp_li_cyp2e1` | 1 | 1.315 | 1.63 | 1.715 | 1.8 | `Ghobadi2011`, `Emery2003`, `OShea1994` |
| CYP3A4 liver | `fobes_li_cyp3a4` | `fimp_li_cyp3a4` | 1 | 0.9065 | 0.8344 | 0.7765 | 0.6886 | `Krogstad2021`, `Jamwal2018`, `Ulvestad2013`, `Ghobadi2011`, `Izat2026` |

??? note "Sources"

    - `Benedek1983`: Benedek IH, Fiske WD 3rd, Griffen WO, Bell RM, Blouin RA, McNamara PJ. Serum alpha 1-acid glycoprotein and the binding of drugs in obesity. Br J Clin Pharmacol 16(6):751-754, 1983. [doi:10.1111/j.1365-2125.1983.tb02258.x](https://doi.org/10.1111/j.1365-2125.1983.tb02258.x)
    - `Benedek1984`: Benedek IH, Blouin RA, McNamara PJ. Serum protein binding and the role of increased alpha 1-acid glycoprotein in moderately obese male subjects. Br J Clin Pharmacol 18(6):941-946, 1984. [doi:10.1111/j.1365-2125.1984.tb02567.x](https://doi.org/10.1111/j.1365-2125.1984.tb02567.x)
    - `Berton2022`: Berton M, Bettonte S, Stader F, Battegay M, Marzolini C. Repository describing the anatomical, physiological, and biological changes in an obese population to inform physiologically based pharmacokinetic models. Clin Pharmacokinet 61(9):1251-1270, 2022. [doi:10.1007/s40262-022-01132-3](https://doi.org/10.1007/s40262-022-01132-3)
    - `Blaak1995`: Blaak EE, van Baak MA, Kemerink GJ, Pakbiers MT, Heidendal GA, Saris WH. Beta-adrenergic stimulation and abdominal subcutaneous fat blood flow in lean, obese, and reduced-obese subjects. Metabolism 44(2):183-187, 1995. [doi:10.1016/0026-0495(95)90262-7](https://doi.org/10.1016/0026-0495(95)90262-7)
    - `Collis2001`: Collis T, Devereux RB, Roman MJ, de Simone G, Yeh J, Howard BV, Fabsitz RR, Welty TK. Relations of stroke volume and cardiac output to body composition: the Strong Heart Study. Circulation 103(6):820-825, 2001. [doi:10.1161/01.CIR.103.6.820](https://doi.org/10.1161/01.CIR.103.6.820)
    - `Emery2003`: Emery MG, Fisher JM, Chien JY, Kharasch ED, Dellinger EP, Kowdley KV, Thummel KE. CYP2E1 activity before and after weight loss in morbidly obese subjects with nonalcoholic fatty liver disease. Hepatology 38(2):428-435, 2003. [doi:10.1053/jhep.2003.50342](https://doi.org/10.1053/jhep.2003.50342)
    - `Ghobadi2011`: Ghobadi C, Johnson TN, Aarabi M, Almond LM, Allabi AC, Rowland-Yeo K, Jamei M, Rostami-Hodjegan A. Application of a systems approach to the bottom-up assessment of pharmacokinetics in obese patients: expected variations in clearance. Clin Pharmacokinet 50(12):809-822, 2011. [doi:10.2165/11594420-000000000-00000](https://doi.org/10.2165/11594420-000000000-00000)
    - `Izat2026`: Izat N, Kangne H, Jansson-Lofmark R, Hertel JK, Rostami-Hodjegan A, Robertsen I, Galetin A. Application of virtual twin PBPK models in individuals with obesity via CYP3A4 phenotyping using endogenous biomarker data. Clin Pharmacol Ther 120(3):762-773, 2026. [doi:10.1002/cpt.70369](https://doi.org/10.1002/cpt.70369)
    - `Jamwal2018`: Jamwal R, de la Monte SM, Ogasawara K, Adusumalli S, Barlock BB, Akhlaghi F. Nonalcoholic fatty liver disease and diabetes are associated with decreased CYP3A4 protein expression and activity in human liver. Mol Pharm 15(7):2621-2632, 2018. [doi:10.1021/acs.molpharmaceut.8b00159](https://doi.org/10.1021/acs.molpharmaceut.8b00159)
    - `Krogstad2021`: Krogstad V, Peric A, Robertsen I, Kringen MK, Vistnes M, Hjelmesaeth J, Sandbu R, Johnson LK, Angeles PC, Jansson-Lofmark R, Karlsson C, Andersson S, Asberg A, Andersson TB, Christensen H. Correlation of body weight and composition with hepatic activities of cytochrome P450 enzymes. J Pharm Sci 110(1):432-437, 2021. [doi:10.1016/j.xphs.2020.10.027](https://doi.org/10.1016/j.xphs.2020.10.027)
    - `Lesser1967`: Lesser GT, Deutsch S. Measurement of adipose tissue blood flow and perfusion in man by uptake of 85Kr. J Appl Physiol 23(5):621-630, 1967. [doi:10.1152/jappl.1967.23.5.621](https://doi.org/10.1152/jappl.1967.23.5.621)
    - `Levitt2007`: Levitt DG. Heterogeneity of human adipose blood flow. BMC Clin Pharmacol 7:1, 2007. [doi:10.1186/1472-6904-7-1](https://doi.org/10.1186/1472-6904-7-1)
    - `OShea1994`: O'Shea D, Davis SN, Kim RB, Wilkinson GR. Effect of fasting and obesity in humans on the 6-hydroxylation of chlorzoxazone: a putative probe of CYP2E1 activity. Clin Pharmacol Ther 56(4):359-367, 1994. [doi:10.1038/clpt.1994.150](https://doi.org/10.1038/clpt.1994.150)
    - `Summers1996`: Summers LK, Samra JS, Humphreys SM, Morris RJ, Frayn KN. Subcutaneous abdominal adipose tissue blood flow: variation within and between subjects and relationship to obesity. Clin Sci (Lond) 91(6):679-683, 1996. [doi:10.1042/cs0910679](https://doi.org/10.1042/cs0910679)
    - `Ulvestad2013`: Ulvestad M, Skottheim IB, Jakobsen GS, Bremer S, Molden E, Asberg A, Hjelmesaeth J, Andersson TB, Sandbu R, Christensen H. Impact of OATP1B1, MDR1, and CYP3A4 expression in liver and intestine on interpatient pharmacokinetic variability of atorvastatin in obese subjects. Clin Pharmacol Ther 93(3):275-282, 2013. [doi:10.1038/clpt.2012.261](https://doi.org/10.1038/clpt.2012.261)
    - `Virtanen2002`: Virtanen KA, Lonnroth P, Parkkola R, Peltoniemi P, Asola M, Viljanen T, Tolvanen T, Knuuti J, Ronnemaa T, Huupponen R, Nuutila P. Glucose uptake and perfusion in subcutaneous and visceral adipose tissue during insulin stimulation in nonobese and obese humans. J Clin Endocrinol Metab 87(8):3902-3910, 2002. [doi:10.1210/jcem.87.8.8761](https://doi.org/10.1210/jcem.87.8.8761)
    - `deDivitiis1981`: de Divitiis O, Fazio S, Petitto M, Maddalena G, Contaldo F, Mancini M. Obesity and cardiac function. Circulation 64(3):477-482, 1981. [doi:10.1161/01.CIR.64.3.477](https://doi.org/10.1161/01.CIR.64.3.477)
<!-- end generated -->

A cell with two numbers is the man and the woman; the others are alike in both sexes to within 0.0005. A protein factor exists only when a process of the model uses that protein in that organ.

The factor of the cardiac output carries **every** organ flow, as the factor of the pregnancy does, so the cardiac output of the model stays the sum of its flows; the factors of the individual organs then shape the distribution. The consequences for the man of 176 cm at a body mass index of 40, against the build without the factors:

| quantity | without | with | Berton 2022 |
|---|---:|---:|---:|
| cardiac output [l/min] | 8.55 | 7.37 | 6.63 |
| brain flow [l/min] | 1.109 | 0.850 | constant by construction |
| muscle flow [l/min] | 1.571 | 1.200 | nearly constant |
| adipose flow [l/min] | 0.462 | 0.529 | +16 % |
| adipose perfusion [ml/min/100 g] | 0.87 | 0.99 | 1.4 per kg of fat |
| total hepatic flow [l/min] | 2.36 | 1.95 | +22 % of 1.79, i.e. 2.19 |
| GFR [ml/min] | 156.5 | 156.5 | 156.3 |

The brain and the muscle keep their absolute flow (0.843 and 1.194 l/min at the reference individual), which is the physiologically correct behaviour and which the library does not have without the factors; the adipose flow rises while its perfusion per kg falls. The factors fix the **shape** of the obesity response and not the reference individual: the cardiac output of the reference man is 6.5 l/min in ICRP 89 and 5.5 in Berton 2022, whose reference is a regression of white adults of 20 to 50 years, so the library does not converge on their absolute values.

The GFR is the strongest result of the comparison and needs **no factor at all**: the `GFR_ref*BSA/BSA_ref` of the library reproduces the measured obese GFR to within 1 % at every body mass index (Berton 2022 from Tc-DTPA and iothalamate studies only), and the lean mass and not the fat drives it (Chew-Harris 2014: muscular men with a body mass index above 30 had the highest GFR, the obese did not; Chagnac 2003 is the upper bound, 145 ml/min at a body mass index of 48). The point to keep in mind when a study reports an eGFR: the **absolute** GFR rises about 30 % at a body mass index of 40 while the **eGFR per 1.73 m2** does not, so `f_renal_function` must not be lowered on the strength of a normal indexed eGFR.

## A change of the body weight

`Adipose(effect=Formula(...))` turns the body weight into a state:

```equations
BW0 = BW
BW_set = BW0 * (1 dimensionless + f_bw)
d(BW)/dt = k_bw * (BW_set - BW)
BW_change = BW - BW0
BW_ratio = BW / BW0
WAIST0 = WAIST
WAIST_change = WAIST - WAIST0
```

with `f_bw` the effect of the model. `BW0` is the body weight of the covariates, the weight at the start of the simulation; `f_bw` is the fractional change of the body weight at the set point (negative for a loss), an assignment rule which may read the species and observers of the top model, e.g. an Emax of the concentration of a GLP-1 agonist; `k_bw` is the rate at which the weight approaches the set point, a half time of 35 weeks by default, which a model fits. The model reports `BW_change` and `BW_ratio`, and `WAIST_change` beside them: the waist follows the weight through its prediction at the ratio `f_waist` of the individual, which is held, so the prediction falls by 0.77 cm per kg in a man and 0.79 in a woman. The trials of the GLP-1 receptor agonists report about 0.85 cm per kg (SURMOUNT-1 lost about 18 cm of waist at 20.9 % of 104.8 kg), so the model under-predicts the waist of a large weight loss by about a tenth: the weight these drugs remove is more central than the weight by which the individuals of a cross section differ.

Because `BW` is a dynamic parameter, every quantity downstream of it becomes an assignment rule: the fat-free mass, every organ mass, volume and flow, the blood volume, the cardiac output, the GFR, the body composition and the factors of obesity. **The physiology follows the weight**, which is the point: a weight loss lowers the cardiac output and the GFR, raises the hepatic CYP3A4 activity back towards its normal weight value and shrinks the adipose volume into which a lipophilic drug distributes.

The split of the lost weight into fat and lean mass is not a parameter: it falls out of the fat-free mass of Janmahasatian 2005. For an obese individual **0.68 to 0.73 of a change of the body weight is adipose tissue**, which is what the body composition substudies measure (0.74 in the SURMOUNT-1 DXA substudy of tirzepatide, 0.709 with an SD of 0.19 over 40 reports of GLP-1 based agents, 0.75 for diet-induced loss, 0.70 to 0.80 after bariatric surgery) and what the Forbes relation predicts at the body composition of those cohorts (`F/(F + 10.4 kg)` = 0.80 for an infinitesimal change at 41 kg of fat, 0.77 integrated over a 21 % loss). The lean limb matters: a model which clamps the lean mass gets the cardiac output, the blood volume and the GFR of a 20 % weight loss wrong by the whole lean limb.

A drug which lowers the body weight by 20 % at its set point, infused to a steady concentration, over a year:

```python
from pathlib import Path

from sbmlsim.model.model_roadrunner import load_roadrunner

from pkpdsim.pbpk import (
    MM,
    Adipose,
    Binding,
    Excretion,
    Formula,
    Kp,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
    U,
    build,
)

drug = Substance(
    sid="drg",
    name="drug",
    mr=300.0,
    binding=Binding(fu=0.1, protein="albumin"),
    distribution=Kp(2.0, ftissue=1.0),
    excretion=Excretion(urine=True),
)
definition = PBPKDefinition(
    sid="obese",
    substances=[drug, Substance(sid="drg_m", name="metabolite", mr=316.0)],
    chemistry=[
        Metabolism(
            sid="DRG2M",
            organ=Organ.LIVER,
            substrate="drg",
            products={"drg_m": 1.0},
            kinetics=MM(vmax=1e-3, km=0.01),
            protein="CYP3A4",
        )
    ],
    routes=[Route.iv("drg"), Route.oral("drg")],
    adipose=Adipose(
        effect=Formula(
            "-0.2 dimensionless * Cve_drg / (Cve_drg + 0.0002 mM)", U.dimensionless
        )
    ),
)
paths = build(definition, Path("."))

r = load_roadrunner(paths.flat_path)
r.setValue("init(BW)", 105.0)
r.setValue("init(Ri_drg)", 5.0e-3)
r.timeCourseSelections = ["time", "BW", "Mad", "BMI", "CO", "GFR", "fimp_li_cyp3a4"]
result = r.simulate(0.0, 52 * 7 * 24 * 60.0, 2000)
names = result.colnames
start, end = result[0], result[-1]


def at(row: "list[float]", name: str) -> float:
    return row[names.index(name)]


for label, row in (("week 0 ", start), ("week 52", end)):
    print(
        f"{label}: BW {at(row, 'BW'):.1f} kg, BMI {at(row, 'BMI'):.1f}, "
        f"Mad {at(row, 'Mad'):.1f} kg, CO {at(row, 'CO'):.2f} l/min, "
        f"GFR {1e3 * at(row, 'GFR'):.0f} ml/min"
    )
d_bw = at(end, "BW") - at(start, "BW")
d_fat = at(end, "Mad") - at(start, "Mad")
print(
    f"change:  {d_bw:.1f} kg ({100 * d_bw / at(start, 'BW'):.1f} %), of which "
    f"{d_fat / d_bw:.2f} adipose tissue and {1 - d_fat / d_bw:.2f} lean mass"
)
print(
    f"CYP3A4 of the liver: {at(start, 'fimp_li_cyp3a4'):.3f} -> "
    f"{at(end, 'fimp_li_cyp3a4'):.3f} of the normal weight activity"
)
```

<!-- generated: figures-adipose-1 (scripts/visual/docs_figures.py), do not edit -->
The content of the model `obese`: the anatomical figure (A) and the schematic (B), hover them for details ([figures](figures.md)).

<iframe class="pkpd-figure" src="img/figures/models/adipose-1-anatomy.html" title="The anatomical figure of the model obese" style="aspect-ratio:1460 / 898" loading="lazy"></iframe>

<iframe class="pkpd-figure" src="img/figures/models/adipose-1-schematic.html" title="The schematic of the model obese" style="aspect-ratio:1384 / 902" loading="lazy"></iframe>
<!-- end generated -->

```
week 0 : BW 105.0 kg, BMI 33.9, Mad 39.4 kg, CO 7.08 l/min, GFR 146 ml/min
week 52: BW 93.9 kg, BMI 30.3, Mad 31.7 kg, CO 6.92 l/min, GFR 139 ml/min
change:  -11.1 kg (-10.6 %), of which 0.68 adipose tissue and 0.32 lean mass
CYP3A4 of the liver: 0.849 -> 0.902 of the normal weight activity
```

The measured trajectories a model is fitted against: liraglutide 3.0 mg -7.9 % in 56 weeks (SCALE), semaglutide 2.4 mg -14.9 % in 68 weeks (STEP 1), tirzepatide 15 mg -20.9 % in 72 weeks (SURMOUNT-1), bariatric surgery -23.4 % in 2 years (SOS). The cleanest kinetic datum is the withdrawal limb: two thirds of the loss returns within 52 weeks of stopping the drug (STEP 1 extension, 11.6 of 17.3 percentage points; SURMOUNT-4, 14.0 of 20.9 %), a half time of 30 to 40 weeks, which is the default of `k_bw`. The curves are not exponential from day 0 because the dose is escalated over 16 to 20 weeks; a model which reproduces the escalation gets the shape from its own pharmacokinetics.

A weight change is slow against every drug half life (a half time of 35 weeks against hours to days), so a model whose weight is static can also scan over body weights with the equilibration of the library (`PBPKSimulationExperiment.equilibrated`), a sequence of static individuals, which is correct whenever the pharmacokinetics are fast compared with the weight change.

## No factor

The quantities for which the comparison with Berton 2022, with Simcyp (Ghobadi 2011) and with the primary literature found no factor the library needs, with the reason:

- **the GFR** (the scaling of the library reproduces the measured obese GFR, see above), the **haematocrit** (no change with the body mass index: more than 2500 NHANES III subjects in Ghobadi 2011, 764 subjects in Berton 2022, the sex the only covariate), the **heart rate** (unchanged; the cardiac output rise is stroke volume), the **kidney flow fraction** (a constant percentage of the cardiac output over 18 studies, so only the absolute flow rises) and the **absolute gastrointestinal blood flow** (unchanged, both platforms agree, which is what the gut fraction factor implements);
- the **kidney, heart and skin masses** (the library agrees with Berton 2022 to within 5 % to a body mass index of 40), the **muscle mass** (the library runs 7 % light at 40, within the noise to 35) and the **blood volume** (within 7 % to 40);
- the **liver mass**: the mass of the library is already a *lean* liver volume, 1.12 to 1.13 times the lean liver volume of Sinha 2020 at every body mass index, and the clearance must scale with the lean part. Berton 2022 add 0.01 kg of liver per kg of adipose tissue, which is the fat of the steatotic liver and raises the distribution volume by 10 % at a body mass index of 40; the library does not add it, because in the library the organ mass also scales the intrinsic clearance and the steatotic fat would raise it in the wrong direction. A model which needs the extra distribution volume adds the term itself;
- the **hepatic artery**: the two platforms disagree on the total hepatic blood flow. Simcyp holds the absolute flow constant (Alexander 1962), Berton 2022 hold the percentage of the cardiac output constant, so the absolute flow rises 22 % at a body mass index of 40 (11 studies, 336 subjects; Knibbe 2015 side with Berton). The library keeps the constant hepatic artery fraction, which with the falling gut fraction puts the total hepatic blood flow of the man of 176 cm at 1.95 l/min at a body mass index of 40 against 1.79 at the reference individual, i.e. **+9 %** (26.5 instead of 27.6 % of the cardiac output), between the two; without the factors it would rise 32 %;
- the **intestinal CYP3A4** (down, Ulvestad 2013 and Krogstad 2020, by an amount nobody has quantified), **OATP1B1** (down with the body weight, slope not published, Wegler 2021), **MRP3** (up in steatohepatitis, not quantified) and the **UGTs** (the in vivo clearance of lorazepam and oxazepam is raised, Abernethy 1983, with no abundance data; the lean-mass-driven liver of the library already raises their clearance roughly in proportion);
- the **tissue composition of adipose**: Martin 1994 give `f_lipid = 0.327 + 0.0124*(percent adiposity)` from six male cadavers, which implies that the conventional neutral lipid fraction of 0.853 already describes a 42 %-fat male. The fixed Rodgers and Rowland composition stays the default, which is what PK-Sim, Simcyp and the published partition coefficient validations use; a model which needs the dependence uses `Kp(...)` overrides;
- the **adipose exchange model**: adipose stays perfusion limited. Levitt 2010 give the quantitative criterion, `FA/F = 1 - exp(-3e6/Kow)`, so diffusion limitation sets in above a log Kow of about 6.2 in humans; the two-compartment adipose of Levitt 2007 (perfusions of 0.074 and 0.014 l/kg/min) maps onto the `SlowPool` the library already has, and a very lipophilic substance can use `KpMethod(Mode.PERMEABILITY, ps=...)` on adipose;
- the **gastric pH** (no correlation with the body mass index against young healthy volunteers, 6 studies), the **gastric emptying and small bowel transit** (the recent studies agree they are similar in the obese and the lean), the **bowel length** (not correlated with the body weight, Tacchino 2015) and **MPPGL/HPGL** (no data in obesity; both platforms keep the healthy per-gram values).

PK-Sim has **no obese population and no obesity disease state**: a high-BMI PK-Sim individual is built by the ordinary algorithm with the fat absorbing the residual weight and the fat blood flow proportional to the fat volume, which is what the library does without the factors. Ghobadi 2011 remains the only public description of Simcyp's `Sim-Obese` and `Sim-Morbidly-Obese`; the EMA qualification annex describes only the healthy North European Caucasian population and says nothing about obesity.

## Limitations

- **The range of validity is a body mass index of 18.5 to 40**, with a warning above it. Berton 2022 declare their own equations robust to 40 and usable with caution to 60, and the fat-free mass of Janmahasatian 2005 saturates (its body mass index term tends to 9270/216 kg/m) and was developed on a cohort reaching about 40: at a body mass index of 60 the fat-free mass of the library is 0.80 of Berton's and the liver, muscle and blood volume follow it down to 0.74 to 0.81. The table of the factors carries the values to 50 and 60 so that a scan does not fall off the end, not because they are supported there.
- `Mad` is a residual and is not `FM/0.853`: the implied neutral lipid fraction of the compartment runs from 0.77 at a body mass index of 20 to 0.95 at 40, against the 0.853 of the tissue composition the partition coefficients use, so the modelled lipid pool is about 18 % above the measured fat mass in a lean man and 5 % below it at a body mass index of 40. Making the adipose mass the fat mass over the lipid fraction plus the unmapped tissues would remove this and would change the SBML of every model above the reference size; it is not done here.
- The alpha-1-acid glycoprotein is the one numeric conflict left in the factors: the regression of Berton 2022 (12 studies, 455 subjects) raises it 1.23-fold at a body mass index of 40 while Simcyp doubles it (2.03-fold at 30 to 40 and 2.25 above 40, from two studies of 18 and 11 subjects). The library follows the larger set; a model whose drug binds mainly to the acid glycoprotein should scan the alternative.
- The visceral depot of a **woman** is over-predicted by 19 % against Haueise 2023 and by 35 % against Linge 2018, while the men agree to 1 % and 10 %; the cause is mostly that the adipose remainder of the library is about 3 l larger than the measured total adipose tissue of a woman (the fat-free mass of Janmahasatian 2005 over-predicts the fat of a lean woman), and partly that neither cohort defines the visceral depot as Shen 2009 do (the German National Cohort segments from the femoral heads to the cardiac apex, the UK Biobank differently again). The waist does not correct it: `f_waist` is 1 for an individual of average shape by construction, and what is wrong is the share, not the waist.
- Visceral fat falls **faster** than total fat during a weight loss (-40.1 % against -33.9 % for a -21.3 % weight change, SURMOUNT-1 DXA substudy). The waist reproduces part of it and only for an individual whose waist exceeds the prediction: the ratio `f_waist` is held while the predicted waist falls with the weight, so the deviation in cm shrinks and the visceral share falls with it, by about 4 % over a 20 % weight loss at `f_waist` = 1.1 against the 9 % the substudy measures. For an individual of average shape the two depots still fall in proportion. A faster visceral limb would need a longitudinal regression, which no open cohort publishes.
- The organs do not all follow the weight at the same rate. The CALERIE ancillary (whole body MRI, 27 adults, 25 % caloric restriction, 13 % weight loss over 12 months) gives elasticities against the body weight of 2.8 for adipose tissue, 0.57 for skeletal muscle, 0.48 for the liver, 0.42 for the heart and about 0 for the kidney and the brain, while the library applies one elasticity (`SV_tis`) to every lean organ. The liver is the clearest exception in the other direction: a low energy diet reduces the liver volume by 13 % within 2 to 4 weeks and then stops, while the body weight goes on falling (Nilsen 2025, 33 articles, 1868 patients), because it is losing steatosis, glycogen and water. A model which simulates the first weeks of a severe deficit or of bariatric surgery needs its own, faster liver term.
- No published PBPK model with a continuously time-varying body weight was found; every obesity and bariatric PBPK model in the literature is a static state. The circadian machinery of the library ([circadian rhythms](circadian.md)) is the architectural precedent and Marok 2024 its published use, at the wrong time scale.
- The factors are those of adults. Obesity in children has its own virtual population (Gerhart 2022) which the library does not carry, and the factors do not interact with the factors of the age, the sex, the pregnancy or the impairments otherwise than by multiplying.
