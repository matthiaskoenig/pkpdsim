# Species

The species is a choice of the build, as the sex: `PBPKDefinition(species=Species.RAT, ...)` builds the same definition (substances, chemistry, routes, variants) on the physiology of a laboratory animal instead of the human of ICRP 89 (`Species.HUMAN`, the default). The library has the rat, the mouse and the pig (#138, a first version). The structure of the model does not change: the same organs, blood pools, intestine and ids, with the reference animal of the species in place of the reference individual of ICRP 89. The model is annotated with the NCBI taxon of the species (`annotation.TAXA`).

A laboratory animal is its reference animal (`database/tables/physiology_<species>.toml`, the structure and conventions of `physiology_human.toml`) scaled to its body weight `BW`, the one covariate of an animal:

- the tissues of the lean organs and the blood in proportion to the body weight (`SV_tis = BW/BW_ref`), adipose the remainder of the body weight, as in the human;
- the cardiac output and the GFR allometrically, `CO = CO_ref*(BW/BW_ref)^b_CO` and `GFR = GFR_ref*(BW/BW_ref)^b_GFR`, with the exponents of the species (`ids.CO_EXPONENT`, `ids.GFR_EXPONENT`; 0.75 where the data give none);
- the flow fractions, the blood of the organs, the vascular and interstitial fractions, the hematocrit and the pools of the reference animal.

An animal has no height, no age, no factors of the age or the sex and no population: the human covariates (`HEIGHT`, `AGE`, `BSA`, `FFM`) are not in its model. The intestine is the segmented intestine with the lumen volumes, transit times and pH of the species (`database/tables/intestine_<species>.toml`; a species without a fed transit time of the stomach empties fed as fasted). The impairment factors are those of the human (1 at the default severities), and the arm vein observer is the venous outflow of skin and muscle of the animal. Not yet for an animal: the nephron kidney (its nephron is the human one, refused), the enzymes and transporters of the species (the activity of a protein is a parameter of the model, the expression profile of CYP3A4 along the intestine is human; an animal places a protein of the intestine by mucosa volume), the tissue composition of the partition coefficients (Rodgers and Rowland, of the rat for every species), and the experiments of the library which convert a measured clearance into a renal function of a human.

```python
from pkpdsim.pbpk import Organ, Sex, Species
from pkpdsim.pbpk.physiology import reference_tissue_volumes

for species in (Species.RAT, Species.MOUSE, Species.PIG):
    v = reference_tissue_volumes(Sex.MALE, species=species)
    print(
        species.value,
        round(v[Organ.LIVER] * 1000, 2),
        round(v[Organ.ADIPOSE] * 1000, 1),
    )
```

The tissue volumes of the liver and of adipose of the reference animals (ml, without their blood), adipose as the remainder of the body weight:

```text
rat 6.88 55.2
mouse 0.9 4.2
pig 528.0 3579.2
```

## Reference animals

The reference animals against the human of ICRP 89 (the human GFR scales with the body surface area); the adipose of the remainder of the body weight against the adipose tissue of the reference (in brackets): the remainder takes the tissues which are no organ of the model (the gut contents, the tissues not listed by the sources), as in the human.

<!-- generated: species-overview (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/physiology_<species>.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| species | sex | BW [kg] | CO [ml/min] | CO scales with | GFR [ml/min] | GFR scales with | blood [ml/kg] | HCT | adipose [% BW] (reference) |
| --- | --- | ---: | ---: | --- | ---: | --- | ---: | ---: | --- |
| human | male | 73 | 6500 | BSA | 125 | BSA | 72.6 | 0.43 | 25.61 (19.86) |
| human | female | 60 | 5900 | BSA | 109 | BSA | 65 | 0.38 | 38.77 (31.67) |
| rat | male | 0.25 | 83.09 | 0.75 | 2.475 | 1 | 63.08 | 0.46 | 20.23 (6.572) |
| mouse | male | 0.025 | 13.98 | 0.75 | 0.2367 | 0.75 | 84.72 | 0.45 | 15.25 (6.932) |
| pig | male | 30 | 4350 | 0.75 | 129.3 | 0.75 | 67.2 | 0.401 | 11.93 (14.93) |
| pig | female | 30 | 4350 | 0.75 | 129.3 | 0.75 | 67.2 | 0.401 | 13.22 (13) |

??? note "Sources"

    - `Brown1997`: Brown RP, Delp MD, Lindstedt SL, Rhomberg LR, Beliles RP. Physiological parameter values for physiologically based pharmacokinetic models. Toxicol Ind Health 13(4):407-484, 1997. [doi:10.1177/074823379701300401](https://doi.org/10.1177/074823379701300401)
    - `Carrara2016`: Carrara F, Azzollini N, Nattino G, Corna D, Villa S, Cerullo D, Zoja C, Abrante B, Luis-Lima S, Porrini E, Cannata A, Ferrari S, Fois M, Stucchi N, Gaspari F. Simplified method to measure glomerular filtration rate by iohexol plasma clearance in conscious rats. Nephron 133(1):62-70, 2016. [doi:10.1159/000445843](https://doi.org/10.1159/000445843)
    - `Davies1993`: Davies B, Morris T. Physiological parameters in laboratory animals and humans. Pharm Res 10(7):1093-1095, 1993. [doi:10.1023/a:1018943613122](https://doi.org/10.1023/a:1018943613122)
    - `Gasthuys2017`: Gasthuys E, Devreese M, Millecam J, Sys S, Vanderperren K, Delanghe J, Vande Walle J, Heyndrickx M, Croubels S. Postnatal maturation of the glomerular filtration rate in conventional growing piglets as potential juvenile animal model for preclinical pharmaceutical research. Front Pharmacol 8:431, 2017. [doi:10.3389/fphar.2017.00431](https://doi.org/10.3389/fphar.2017.00431)
    - `Hannon1985`: Hannon JP, Bossone CA, Rodkey WG. Splenic red cell sequestration and blood volume measurements in conscious pigs. Am J Physiol 248(3):R293-R301, 1985. [doi:10.1152/ajpregu.1985.248.3.r293](https://doi.org/10.1152/ajpregu.1985.248.3.r293)
    - `He2017`: He Q, Su G, Liu K, Zhang F, Jiang Y, Gao J, Liu L, Jiang Z, Jin M, Xie H. Sex-specific reference intervals of hematologic and biochemical analytes in Sprague-Dawley rats using the nonparametric rank percentile method. PLoS One 12(12):e0189837, 2017. [doi:10.1371/journal.pone.0189837](https://doi.org/10.1371/journal.pone.0189837)
    - `ICRP89`: Valentin J (ed). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89. Ann ICRP 32(3-4), 2002. [doi:10.1016/S0146-6453(03)00002-2](https://doi.org/10.1016/S0146-6453(03)00002-2)
    - `Lee1985`: Lee HB, Blaufox MD. Blood volume in the rat. J Nucl Med 26(1):72-76, 1985. PMID 3965655.
    - `Lin2020`: Lin Z, Li M, Wang YS, Tell LA, Baynes RE, Davis JL, Vickroy TW, Riviere JE. Physiological parameter values for physiologically based pharmacokinetic models in food-producing animals. Part I: Cattle and swine. J Vet Pharmacol Ther 43(5):385-420, 2020. [doi:10.1111/jvp.12861](https://doi.org/10.1111/jvp.12861)
    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
    - `Qi2004`: Qi Z, Whitt I, Mehta A, Jin J, Zhao M, Harris RC, Fogo AB, Breyer MD. Serial determination of glomerular filtration rate in conscious mice using FITC-inulin clearance. Am J Physiol Renal Physiol 286(3):F590-F596, 2004. [doi:10.1152/ajprenal.00324.2003](https://doi.org/10.1152/ajprenal.00324.2003)
    - `Riches1973`: Riches AC, Sharp JG, Thomas DB, Smith SV. Blood volume determination in the mouse. J Physiol 228(2):279-284, 1973. [doi:10.1113/jphysiol.1973.sp010086](https://doi.org/10.1113/jphysiol.1973.sp010086)
<!-- end generated -->

- **Rat** (male, 250 g, Sprague-Dawley or Wistar of 8-10 weeks): the organ weights, regional flows and residual blood of the compilation of Brown 1997 (Tables 5, 21, 22, 25, 30), the hematocrit of He 2017, the GFR of Carrara 2016 (iohexol, 0.99 ml/min/100 g; Davies and Morris 1993 give half of it), the blood volume of Lee 1985; the interstitial fractions of the rat of httk (Ruark 2014), the portal flow split onto the portal organs and the intestine (geometry, pH, transit of the small and large intestine) of the PK-Sim database (Thiel 2015).
- **Mouse** (male, 25 g): Brown 1997 (organ weights, cardiac output 13.98 ml/min, Table 22, flows, residual blood), the GFR of Qi 2004 (FITC-inulin), the blood volume of Riches 1973; adipose and bone flow of Davies and Morris 1993 (Brown 1997 has none), the interstitial fractions of the rat (httk has no mouse).
- **Pig** (domestic pig, Landrace x Yorkshire, growing, 30 kg, male and female): the meta-analysis of Lin 2020 (organ weights by sex, Tables 9-11; cardiac output 8.70 l/h/kg of 83 unanesthetized pigs, Table 17; flow fractions and perfusion, Tables 27-28; hematocrit), the blood volume of Hannon 1985 (67.2 ml/kg), the GFR of Gasthuys 2017 (iohexol, 4.31 ml/min/kg); the rest of the body of Lin 2020 (head, feet, tail, gut contents, connective tissue, 21 % of the body weight) is shared by muscle, skin and bone in proportion to their masses (else adipose would take a third of the body weight). The vascular fractions, densities and interstitial fractions of the organs are those of the human (no data of the pig). The Göttingen minipig of the pharmacology is not tabulated: its open data cover a few organs only (Ellegaard 2009); Suenderhauf and Parrott 2013 would complete it.

The organs of the reference rat:

<!-- generated: species-organs-rat (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/physiology_rat.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| organ | mass [kg] | density [kg/l] | blood [-] | f_vas [-] | f_int [-] | flow [-] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `lu` lung | 0.0008 | 1.05 | 0.02718 | 0.36 | 0.208 | 1 |
| `he` heart | 0.0006105 | 1.03 | 0.01321 | 0.26 | 0.1036 | 0.049 |
| `br` brain | 0.001382 | 1.035 | 0.002619 | 0.03 | 0.1571 | 0.02 |
| `li` liver | 0.007229 | 1.05 | 0.116 | 0.21 | 0.1003 | 0.024 |
| `ki` kidney | 0.001533 | 1.05 | 0.01763 | 0.16 | 0.1226 | 0.141 |
| `st` stomach | 0.001113 | 1.044 | 0.002235 | 0.032 | 0.09293 | 0.01749 |
| `si` small intestine | 0.003416 | 1.044 | 0.005102 | 0.024 | 0.0937 | 0.03858 |
| `la` large intestine | 0.00205 | 1.044 | 0.003061 | 0.024 | 0.0937 | 0.07718 |
| `sp` spleen | 0.00039 | 1.054 | 0.006618 | 0.22 | 0.1989 | 0.00977 |
| `pa` pancreas | 0.000656 | 1.044 | 0.008746 | 0.18 | 0.07872 | 0.007973 |
| `mu` muscle | 0.09703 | 1.041 | 0.2463 | 0.04 | 0.119 | 0.278 |
| `sk` skin | 0.04662 | 1.117 | 0.05402 | 0.02 | 0.3126 | 0.058 |
| `ad` adipose | 0.01643 | 0.916 | 0.01149 | 0.01 | 0.1386 | 0.07 |
| `bo` bone | 0.01752 | 1.49 | 0.03107 | 0.04 | 0.096 | 0.122 |
| `go` gonads | 0.00236 | 1 | 0.02436 | 0.14 | 0.08256 | 0.006341 |

Sources of the columns: mass, flow: `Brown1997`, `PKSimDB`; density: `Brown1997`, `httk`; blood, f_vas: `Brown1997`, `Kawai1994`, `PKSimDB`; f_int: `httk`.

??? note "Sources"

    - `Brown1997`: Brown RP, Delp MD, Lindstedt SL, Rhomberg LR, Beliles RP. Physiological parameter values for physiologically based pharmacokinetic models. Toxicol Ind Health 13(4):407-484, 1997. [doi:10.1177/074823379701300401](https://doi.org/10.1177/074823379701300401)
    - `Kawai1994`: Kawai R, Lemaire M, Steimer JL, Bruelisauer A, Niederberger W, Rowland M. Physiologically based pharmacokinetic study on a cyclosporin derivative, SDZ IMM 125. J Pharmacokinet Biopharm 22(5):327-365, 1994. [doi:10.1007/BF02353860](https://doi.org/10.1007/BF02353860)
    - `PKSimDB`: Open Systems Pharmacology. PK-Sim database src/Db/PKSimDB.sqlite, github.com/Open-Systems-Pharmacology/PK-Sim, branch develop, commit aaebb648cbc3e4bf71fb9768916932f881dde19d (downloaded 2026-11-06), GPL-2.0; the rat and mouse species of Thiel 2015.
    - `httk`: Pearce RG, Setzer RW, Strope CL, Wambaugh JF, Sipes NS. httk: R package for high-throughput toxicokinetics. J Stat Softw 79(4):1-26, 2017. R package httk version 2.7.4 (CRAN, 2025-12-08), MIT license, copyright 2025 John Wambaugh; the data set mecdt of HTTK-Pop (Ring 2017) holds the preprocessed NHANES 2013-2018 respondents. [doi:10.18637/jss.v079.i04](https://doi.org/10.18637/jss.v079.i04)
<!-- end generated -->

The organs of the reference mouse:

<!-- generated: species-organs-mouse (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/physiology_mouse.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| organ | mass [kg] | density [kg/l] | blood [-] | f_vas [-] | f_int [-] | flow [-] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `lu` lung | 9.125e-05 | 1.05 | 0.04104 | 0.5 | 0.1625 | 1 |
| `he` heart | 9.225e-05 | 1.03 | 0.01502 | 0.262 | 0.1033 | 0.066 |
| `br` brain | 0.0004001 | 1.035 | 0.005647 | 0.03 | 0.1571 | 0.033 |
| `li` liver | 0.000947 | 1.05 | 0.1914 | 0.31 | 0.08763 | 0.02 |
| `ki` kidney | 0.0003173 | 1.05 | 0.04507 | 0.24 | 0.111 | 0.091 |
| `st` stomach | 0.0001452 | 1.044 | 0.002171 | 0.032 | 0.09293 | 0.008855 |
| `si` small intestine | 0.0006173 | 1.044 | 0.006867 | 0.024 | 0.0937 | 0.08048 |
| `la` large intestine | 0.000266 | 1.044 | 0.002958 | 0.024 | 0.0937 | 0.04024 |
| `sp` spleen | 7.262e-05 | 1.054 | 0.006665 | 0.17 | 0.2117 | 0.007247 |
| `pa` pancreas | 0.0001268 | 1.044 | 0.01259 | 0.18 | 0.07872 | 0.004188 |
| `mu` muscle | 0.009216 | 1.041 | 0.1742 | 0.04 | 0.119 | 0.159 |
| `sk` skin | 0.004009 | 1.117 | 0.05242 | 0.03 | 0.3094 | 0.058 |
| `ad` adipose | 0.001733 | 0.916 | 0.009022 | 0.01 | 0.1386 | 0.027 |
| `bo` bone | 0.002387 | 1.49 | 0.09354 | 0.11 | 0.089 | 0.05828 |
| `go` gonads | 0.0002558 | 1 | 0.01967 | 0.14 | 0.08256 | 0.004085 |

Sources of the columns: mass: `Brown1997`, `PKSimDB`; density: `Brown1997`, `httk`; blood, f_vas: `Brown1997`, `Kawai1994`, `PKSimDB`; f_int: `httk`; flow: `Brown1997`, `Davies1993`, `PKSimDB`.

??? note "Sources"

    - `Brown1997`: Brown RP, Delp MD, Lindstedt SL, Rhomberg LR, Beliles RP. Physiological parameter values for physiologically based pharmacokinetic models. Toxicol Ind Health 13(4):407-484, 1997. [doi:10.1177/074823379701300401](https://doi.org/10.1177/074823379701300401)
    - `Davies1993`: Davies B, Morris T. Physiological parameters in laboratory animals and humans. Pharm Res 10(7):1093-1095, 1993. [doi:10.1023/a:1018943613122](https://doi.org/10.1023/a:1018943613122)
    - `Kawai1994`: Kawai R, Lemaire M, Steimer JL, Bruelisauer A, Niederberger W, Rowland M. Physiologically based pharmacokinetic study on a cyclosporin derivative, SDZ IMM 125. J Pharmacokinet Biopharm 22(5):327-365, 1994. [doi:10.1007/BF02353860](https://doi.org/10.1007/BF02353860)
    - `PKSimDB`: Open Systems Pharmacology. PK-Sim database src/Db/PKSimDB.sqlite, github.com/Open-Systems-Pharmacology/PK-Sim, branch develop, commit aaebb648cbc3e4bf71fb9768916932f881dde19d (downloaded 2026-11-06), GPL-2.0; the rat and mouse species of Thiel 2015.
    - `httk`: Pearce RG, Setzer RW, Strope CL, Wambaugh JF, Sipes NS. httk: R package for high-throughput toxicokinetics. J Stat Softw 79(4):1-26, 2017. R package httk version 2.7.4 (CRAN, 2025-12-08), MIT license, copyright 2025 John Wambaugh; the data set mecdt of HTTK-Pop (Ring 2017) holds the preprocessed NHANES 2013-2018 respondents. [doi:10.18637/jss.v079.i04](https://doi.org/10.18637/jss.v079.i04)
<!-- end generated -->

The organs of the reference pig (male):

<!-- generated: species-organs-pig (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/physiology_pig.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| organ | mass [kg] | density [kg/l] | blood [-] | f_vas [-] | f_int [-] | flow [-] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `lu` lung | 0.24 | 1 | 0.1577 | 0.5699 | 0.08602 | 1 |
| `he` heart | 0.102 | 1 | 0.008127 | 0.1384 | 0.2154 | 0.03 |
| `br` brain | 0.066 | 1.04 | 0.001436 | 0.04363 | 0.004782 | 0.015 |
| `li` liver | 0.528 | 1 | 0.07713 | 0.2275 | 0.1777 | 0.044 |
| `ki` kidney | 0.114 | 1.05 | 0.01934 | 0.2642 | 0.1214 | 0.114 |
| `st` stomach | 0.198 | 1.045 | 0.03471 | 0.2697 | 0.07011 | 0.021 |
| `si` small intestine | 0.657 | 1.045 | 0.101 | 0.2446 | 0.07252 | 0.153 |
| `la` large intestine | 0.441 | 1.045 | 0.06892 | 0.2477 | 0.07222 | 0.051 |
| `sp` spleen | 0.06 | 1.06 | 0.01472 | 0.344 | 0.1673 | 0.031 |
| `pa` pancreas | 0.045 | 1 | 0.00507 | 0.1851 | 0.07823 | 0.014 |
| `mu` muscle | 16.23 | 1 | 0.206 | 0.02495 | 0.1463 | 0.2679 |
| `sk` skin | 1.886 | 1.1 | 0.04506 | 0.05033 | 0.5698 | 0.035 |
| `ad` adipose | 4.479 | 1 | 0.04061 | 0.01795 | 0.1375 | 0.04098 |
| `bo` bone | 3.628 | 1.3 | 0.06359 | 0.04392 | 0.09561 | 0.01836 |
| `go` gonads | 0.08884 | 1.04 | 0.002669 | 0.05926 | 0.09031 | 0.001896 |

Sources of the columns: mass: `Ellegaard2009`, `Lin2020`, `ModelChoice`; density, blood, f_vas: `ModelChoice`; f_int: `httk`; flow: `Lin2020`, `ModelChoice`, `Moher2002`.

??? note "Sources"

    - `Ellegaard2009`: Ellegaard Göttingen Minipigs. Organ weight 2009: absolute (g) and relative (% of body weight) organ weight, group mean values, 8 weeks and 6 months. https://minipigs.dk/sites/default/files/2024-05/Organ%20Weights.pdf (accessed 2026-11-06).
    - `Lin2020`: Lin Z, Li M, Wang YS, Tell LA, Baynes RE, Davis JL, Vickroy TW, Riviere JE. Physiological parameter values for physiologically based pharmacokinetic models in food-producing animals. Part I: Cattle and swine. J Vet Pharmacol Ther 43(5):385-420, 2020. [doi:10.1111/jvp.12861](https://doi.org/10.1111/jvp.12861)
    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
    - `Moher2002`: Moher HE, Carey GB. Technical note: adipose tissue blood flow in miniature swine (Sus scrofa) using the 133xenon washout technique. J Anim Sci 80(5):1294-1298, 2002. [doi:10.2527/2002.8051294x](https://doi.org/10.2527/2002.8051294x)
    - `httk`: Pearce RG, Setzer RW, Strope CL, Wambaugh JF, Sipes NS. httk: R package for high-throughput toxicokinetics. J Stat Softw 79(4):1-26, 2017. R package httk version 2.7.4 (CRAN, 2025-12-08), MIT license, copyright 2025 John Wambaugh; the data set mecdt of HTTK-Pop (Ring 2017) holds the preprocessed NHANES 2013-2018 respondents. [doi:10.18637/jss.v079.i04](https://doi.org/10.18637/jss.v079.i04)
<!-- end generated -->

## Intestine

<!-- generated: species-intestine (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/intestine_<species>.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| species | stomach fasted [min] | stomach fed [min] | small intestine [min] | large intestine [h] | stomach [ml] | small intestine [ml] | large intestine [ml] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| human | 22.79 | 63.48 | 192 | 35.1 | 28 | 105 | 13 |
| rat | 30 |  | 88 | 13.5 | 1 | 1.2 | 1 |
| mouse | 20 |  | 40.9 | 4.7 | 0.2 | 0.2999 | 0.15 |
| pig | 121.2 | 121.2 | 210 | 30.48 | 81.07 | 139.6 | 229.7 |

??? note "Sources"

    - `Davies1993`: Davies B, Morris T. Physiological parameters in laboratory animals and humans. Pharm Res 10(7):1093-1095, 1993. [doi:10.1023/a:1018943613122](https://doi.org/10.1023/a:1018943613122)
    - `Davis2001`: Davis SS, Illum L, Hinchcliffe M. Gastrointestinal transit of dosage forms in the pig. J Pharm Pharmacol 53(1):33-39, 2001. [doi:10.1211/0022357011775163](https://doi.org/10.1211/0022357011775163)
    - `Merchant2011`: Merchant HA, McConnell EL, Liu F, Ramaswamy C, Kulkarni RP, Basit AW, Murdan S. Assessment of gastrointestinal pH, fluid and lymphoid tissue in the guinea pig, rabbit and pig, and implications for their use in drug development. Eur J Pharm Sci 42(1-2):3-10, 2011. [doi:10.1016/j.ejps.2010.09.019](https://doi.org/10.1016/j.ejps.2010.09.019)
    - `Metcalf1987`: Metcalf AM, Phillips SF, Zinsmeister AR, MacCarty RL, Beart RW, Wolff BG. Simplified assessment of segmental colonic transit. Gastroenterology 92(1):40-47, 1987. [doi:10.1016/0016-5085(87)90837-7](https://doi.org/10.1016/0016-5085(87)90837-7)
    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
    - `Mudie2010`: Mudie DM, Amidon GL, Amidon GE. Physiological parameters for oral delivery and in vitro testing. Mol Pharm 7(5):1388-1405, 2010. [doi:10.1021/mp100149j](https://doi.org/10.1021/mp100149j)
    - `PKSimDB`: Open Systems Pharmacology. PK-Sim database src/Db/PKSimDB.sqlite, github.com/Open-Systems-Pharmacology/PK-Sim, branch develop, commit aaebb648cbc3e4bf71fb9768916932f881dde19d (downloaded 2026-11-06), GPL-2.0; the rat and mouse species of Thiel 2015.
    - `Schiller2005`: Schiller C, Fröhlich CP, Giessmann T, Siegmund W, Mönnikes H, Hosten N, Weitschies W. Intestinal fluid volumes and transit of dosage forms as assessed by magnetic resonance imaging. Aliment Pharmacol Ther 22(10):971-979, 2005. [doi:10.1111/j.1365-2036.2005.02683.x](https://doi.org/10.1111/j.1365-2036.2005.02683.x)
<!-- end generated -->

## Limitations

- Many numbers of the rodents come from compilations (Brown 1997, Davies and Morris 1993 through httk, the PK-Sim database) and not from the primary studies; the notes of every number say which. The split of the portal flow of the rat (PK-Sim) gives the large intestine twice the flow of the small intestine; Delp 1991 was not read.
- The mouse has no measured flow of adipose, bone and the adrenals in Brown 1997: the tissues outside the organs take 34 % of its cardiac output, which the organs share in proportion.
- The arterial and venous pools of the rodents divide the blood outside the organs 25:75 (Brown 1997), their portal vein and the blood density are those of the human.
- The organs of an animal of another weight are in the proportions of the reference animal; the allometry of the organs within a species (the brain of a heavier rat) is not modelled, nor the growth of the young animal.
- Only the male rat and mouse: the organ weights, flows and residual blood of Brown 1997 are not given per sex (the female hematocrit, GFR and blood volume of the sources are listed in the notes).
