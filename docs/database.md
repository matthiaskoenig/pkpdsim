# Parameter database

The numbers which the PBPK library and the models are built from are the parameter database `pkpdsim.database`, apart from the code which uses them: the library reads the database and holds no numbers of its own. The database has three parts.

- **Sources** (`pkpdsim/database/sources.toml`): the publications every number is taken from, by key (e.g. `ICRP89`), with the citation and the DOI, and `ModelChoice` for a number the library chooses where the literature gives none.
- **Tables** (`pkpdsim/database/tables/*.toml`): the tables of the library, one file per topic: the physiology of the reference individuals, the populations, the distribution, the intestine, the kidney, the routes, the proteins and the impairment.
- **Studies** (`pkpdsim/database/studies/<source>/`): the data of the studies as their sources report them, the groups and the statistics of every measurement, from which the database integrates a consensus with its uncertainty.

## Tables

A number of a table is `{ value = 1.8, unit = "kg", source = "ICRP89", notes = "..." }`: the number, its unit, the key of its source and a remark. A number states its uncertainty where its source reports it: `sd` (the standard deviation of the individuals), `se` (the standard error of the number), `cv` (the coefficient of variation, a fraction), `n` (the number of individuals) and `ci` (`[lower, upper]`, the 95 % confidence interval). `pkpdsim.database.tables.entry` reads a number with its uncertainty, and the library writes the uncertainty into the SBML of the parameter as a distrib uncertainty (`pkpdsim.pbpk.data.value`): the standard deviation, the standard error, the coefficient of variation, the sample size and the confidence interval, no distribution, so the draws of `run_simulation --uncertainty` stay those of the fitted parameters.

```python
from pkpdsim.database import tables

e = tables.entry("routes", "iv", "ti")
print(e.value, e.unit, e.source)
```

## Studies

A study is a folder named by the key of its source. `study.toml` describes the tissue (`liver`, or `body` for the quantities of the subjects), the control group (left out for a study without, whose records give levels and no ratios), the groups (`condition` `control` or `cirrhosis`, the Child-Pugh class, the etiology, an `assumption` the source does not state), `via` (the key of a review whose table the numbers are read from, when the study itself was not read) and the measurements (`kind` `protein`, the abundance per amount of tissue, `mrna`, the expression, or `physiology`, a physiological or clinical quantity of the subjects; the unit, the method and where the values are in the source). Every measurement is a TSV file of records, one row per entity and group, with the statistics the source reports: `n`, `n_detected` (the samples above the limit of quantification), `mean`, `sd`, `se`, `cv_percent`, `median`, `q1` and `q3` (the quartiles), `min`, `max`, `unit` (of a record whose unit is not that of its measurement) and `ci_low`/`ci_high` (the 95 % confidence interval of the median). The entity of a protein or an mRNA is a gene symbol of HGNC (e.g. `CYP3A4` or `SLCO1B1`, with the name of the source), that of the kind `physiology` a quantity of `pkpdsim.database.studies.QUANTITIES` (e.g. `albumin`, `cardiac_index`, `meld`), whose records are all in the unit of the quantity: a source which prints another unit (albumin in g/dL) is converted when its numbers are entered, and the note of the record gives the printed numbers.

The records are the numbers of the source as printed, never derived ones. A statistic the source does not report is empty. A number which disagrees with the others of its row (a mean outside its range, a CV which is not the SD over the mean, a sum of the classes which is not their group) is kept as printed and has a note which says what the other numbers imply; `pkpdsim.database.studies.check` reports a disagreement without a note, and the tests keep it empty.

<!-- generated: database-studies (scripts/pbpk_tables.py), do not edit -->
*Generated from `pkpdsim/database/studies` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| study | measurement | unit | entities | records | groups | location |
| --- | --- | --- | ---: | ---: | ---: | --- |
| `Alberino2001` | a physiological or clinical quantity of the subjects |  | 1 | 3 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Annet2003` | a physiological or clinical quantity of the subjects |  | 3 | 12 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Boker2000` | a physiological or clinical quantity of the subjects |  | 1 | 2 | 2 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Chawla1998` | a physiological or clinical quantity of the subjects |  | 1 | 4 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Dincer2005` | a physiological or clinical quantity of the subjects |  | 1 | 4 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Drozdzik2020` | protein abundance per amount of tissue | fmol/mg tissue | 14 | 126 | 9 | Figures 1, 2, 5 and 6 (right panels), digitized |
| `Drozdzik2020` | mRNA expression | relative to the mean of the controls (RQ) | 14 | 126 | 9 | Figures 1, 2, 5 and 6 (left panels), digitized |
| `Drozdzik2021` | protein abundance per amount of tissue | fmol/mg tissue | 13 | 299 | 23 | Table S3 |
| `Drozdzik2021` | mRNA expression | relative to the mean of the controls (2^-ddCt) | 13 | 299 | 23 | Table S2 |
| `ElKhateeb2021` | protein abundance per amount of tissue | pmol/g liver tissue | 51 | 201 | 4 | Table S8 of the supporting information |
| `Froomes1999` | a physiological or clinical quantity of the subjects |  | 1 | 3 | 3 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Henriksen2001` | a physiological or clinical quantity of the subjects |  | 1 | 4 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Hong2025` | a physiological or clinical quantity of the subjects |  | 3 | 9 | 4 | Table 1 |
| `Iwao1997` | a physiological or clinical quantity of the subjects |  | 1 | 4 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Kayacetin2004` | a physiological or clinical quantity of the subjects |  | 2 | 8 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Kim2020` | a physiological or clinical quantity of the subjects |  | 3 | 9 | 4 | Table 1 |
| `Kwon2001` | a physiological or clinical quantity of the subjects |  | 3 | 9 | 3 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Lin1998` | a physiological or clinical quantity of the subjects |  | 1 | 4 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Liu2006` | a physiological or clinical quantity of the subjects |  | 1 | 2 | 3 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Matsui1996` | a physiological or clinical quantity of the subjects |  | 2 | 4 | 3 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Mesropyan2022` | a physiological or clinical quantity of the subjects |  | 7 | 28 | 4 | Table 1 |
| `Prasad2018` | protein abundance per amount of tissue | pmol/g tissue | 18 | 54 | 3 | Table 3 |
| `Sudhamshu2006` | a physiological or clinical quantity of the subjects |  | 2 | 8 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `VanBeers2001` | a physiological or clinical quantity of the subjects |  | 2 | 8 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Viani1992` | a physiological or clinical quantity of the subjects |  | 1 | 4 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Woitas2000` | a physiological or clinical quantity of the subjects |  | 1 | 3 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Wong2001` | a physiological or clinical quantity of the subjects |  | 1 | 3 | 3 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |
| `Xu2026` | a physiological or clinical quantity of the subjects |  | 3 | 9 | 4 | Table 1 |
| `Yang2020` | a physiological or clinical quantity of the subjects |  | 3 | 12 | 4 | Table 1 |
| `Zhu1999` | a physiological or clinical quantity of the subjects |  | 1 | 4 | 4 | Li 2015, Supplementary Materials, Physiological and biological values used in meta-analysis |

??? note "Sources"

    - `Alberino2001`: Alberino F, Gatta A, Amodio P, Merkel C, Di Pascoli L, Boffo G, Caregaro L. Nutrition and survival in patients with liver cirrhosis. Nutrition 17(6):445-450, 2001. [doi:10.1016/s0899-9007(01)00521-4](https://doi.org/10.1016/s0899-9007(01)00521-4)
    - `Annet2003`: Annet L, Materne R, Danse E, Jamart J, Horsmans Y, Van Beers BE. Hepatic Flow Parameters Measured with MR Imaging and Doppler US: Correlations with Degree of Cirrhosis and Portal Hypertension. Radiology 229(2):409-414, 2003. [doi:10.1148/radiol.2292021128](https://doi.org/10.1148/radiol.2292021128)
    - `Boker2000`: Böker KH, Franzen A, Wrobel M, Bahr MJ, Tietge U, Manns MP. Regulation of hepatic blood flow in patients with liver cirrhosis and after liver transplantation. Pathophysiology 6(4):251-260, 2000. [doi:10.1016/s0928-4680(99)00024-3](https://doi.org/10.1016/s0928-4680(99)00024-3)
    - `Chawla1998`: Chawla Y, Radha NS, Jang KD, Dilawari B. Portal Hemodynamics by Duplex Doppler Sonography in Different Grades of Cirrhosis. Dig Dis Sci 43(2):354-357, 1998. [doi:10.1023/a:1018814624307](https://doi.org/10.1023/a:1018814624307)
    - `Dincer2005`: Dincer D, Besisk F, Demirkol O, Demir K, Kaymakoglu S, Cakaloglu Y, Okten A. Relationships between hemodynamic alterations and Child-Pugh Score in patients with cirrhosis. Hepatogastroenterology 52(65):1521-1525, 2005. PMID 16201110.
    - `Drozdzik2020`: Drozdzik M, Szelag-Pieniek S, Post M, Zeair S, Wrzesinski M, Kurzawski M, Prieto J, Oswald S. Protein abundance of hepatic drug transporters in patients with different forms of liver damage. Clin Pharmacol Ther 107(5):1138-1148, 2020. [doi:10.1002/cpt.1717](https://doi.org/10.1002/cpt.1717)
    - `Drozdzik2021`: Drozdzik M, Lapczuk-Romanska J, Wenzel C, Szelag-Pieniek S, Post M, Skalski Ł, Kurzawski M, Oswald S. Gene expression and protein abundance of hepatic drug metabolizing enzymes in liver pathology. Pharmaceutics 13(9):1334, 2021. [doi:10.3390/pharmaceutics13091334](https://doi.org/10.3390/pharmaceutics13091334)
    - `ElKhateeb2021`: El-Khateeb E, Achour B, Al-Majdoub ZM, Barber J, Rostami-Hodjegan A. Non-uniformity of changes in drug-metabolizing enzymes and transporters in liver cirrhosis: implications for drug dosage adjustment. Mol Pharm 18(9):3563-3577, 2021. [doi:10.1021/acs.molpharmaceut.1c00462](https://doi.org/10.1021/acs.molpharmaceut.1c00462)
    - `Froomes1999`: Froomes PR, Morgan DJ, Smallwood RA, Angus PW. Comparative effects of oxygen supplementation on theophylline and acetaminophen clearance in human cirrhosis. Gastroenterology 116(4):915-920, 1999. [doi:10.1016/s0016-5085(99)70075-2](https://doi.org/10.1016/s0016-5085(99)70075-2)
    - `Henriksen2001`: Henriksen JH, Bendtsen F, Møller S. Normal red cell cardiac output in the hyperkinetic syndrome of alcoholic cirrhosis. J Hepatol 34(5):782-783, 2001. [doi:10.1016/s0168-8278(01)00019-8](https://doi.org/10.1016/s0168-8278(01)00019-8)
    - `Hong2025`: Hong SJ, Kim JE, Cho JM, Choi HC, Park MJ, Choi HY, Shin HS, Won JH, Yang W, Kim HO. Dual Energy CT-Derived Liver Extracellular Volume Fraction for Assessing Liver Functional Reserve in Patients with Liver Cirrhosis. Medicina 61(9):1561, 2025. [doi:10.3390/medicina61091561](https://doi.org/10.3390/medicina61091561)
    - `Iwao1997`: Iwao T, Oho K, Sakai T, Tayama C, Sato M, Nakano R, Yamawaki M, Toyonaga A, Tanikawa K. Splanchnic and extrasplanchnic arterial hemodynamics in patients with cirrhosis. J Hepatol 27(5):817-823, 1997. [doi:10.1016/s0168-8278(97)80318-2](https://doi.org/10.1016/s0168-8278(97)80318-2)
    - `Kayacetin2004`: Kayacetin E, Efe D, Doğan C. Portal and splenic hemodynamics in cirrhotic patients: relationship between esophageal variceal bleeding and the severity of hepatic failure. J Gastroenterol 39(7):661-667, 2004. [doi:10.1007/s00535-003-1362-x](https://doi.org/10.1007/s00535-003-1362-x)
    - `Kim2020`: Kim TH, Yun SG, Choi J, Goh HG, Lee HA, Yim SY, Choi SJ, Lee YS, Yoon EL, Jung YK, Seo YS, Kim JH, Yim HJ, Yeon JE, Byun KS, Um SH. Differential Impact of Serum 25-Hydroxyvitamin D3 Levels on the Prognosis of Patients with Liver Cirrhosis According to MELD and Child-Pugh Scores. J Korean Med Sci 35(19):e129, 2020. [doi:10.3346/jkms.2020.35.e129](https://doi.org/10.3346/jkms.2020.35.e129)
    - `Kwon2001`: Kwon HA, Matsui Y, Ha-Kawa SK, Kamiyama Y. Functional Hepatic Volume Measured by Technetium-99M-Galactosyl-Human Serum Albumin Liver Scintigraphy: Comparison Between Hepatocyte Volume and Liver Volume by Computed Tomography. Am J Gastroenterol 96(2):541-546, 2001. [doi:10.1111/j.1572-0241.2001.03556.x](https://doi.org/10.1111/j.1572-0241.2001.03556.x)
    - `Lin1998`: Lin XZ, Sun YN, Liu YH, Sheu BS, Cheng BN, Chen CY, Tsai HM, Shen CL. Liver volume in patients with or without chronic liver diseases. Hepatogastroenterology 45(22):1069-1074, 1998. PMID 9756008.
    - `Liu2006`: Liu T, Wong W, Hou M, Lin H, Chang F, Lee S. Hemorheology in patients with liver cirrhosis: Special emphasis on its relation to severity of esophageal variceal bleeding. J Gastroenterol Hepatol 21(5):908-913, 2006. [doi:10.1111/j.1440-1746.2006.04266.x](https://doi.org/10.1111/j.1440-1746.2006.04266.x)
    - `Matsui1996`: Matsui Y, Tu W, Kitade H, Nakagawa A, Kamiya T, Kwon A, Uetsuji S, Kamiyama Y. Hepatocyte volume as an indicator of hepatic functional reserve in cirrhotic patients with liver tumours. J Gastroenterol Hepatol 11(6):540-545, 1996. [doi:10.1111/j.1440-1746.1996.tb01699.x](https://doi.org/10.1111/j.1440-1746.1996.tb01699.x)
    - `Mesropyan2022`: Mesropyan N, Kupczyk PA, Dold L, Praktiknjo M, Chang J, Isaak A, Endler C, Kravchenko D, Bischoff LM, Sprinkart AM, Pieper CC, Kuetting D, Jansen C, Attenberger UI, Luetkens JA. Assessment of liver cirrhosis severity with extracellular volume fraction MRI. Sci Rep 12(1):9422, 2022. [doi:10.1038/s41598-022-13340-9](https://doi.org/10.1038/s41598-022-13340-9)
    - `Prasad2018`: Prasad B, Bhatt DK, Johnson K, Chapa R, Chu X, Salphati L, Xiao G, Lee C, Hop CE, Mathias A, Lai Y, Liao M, Humphreys WG, Kumer SC, Unadkat JD. Abundance of phase 1 and 2 drug-metabolizing enzymes in alcoholic and hepatitis C cirrhotic livers: a quantitative targeted proteomics study. Drug Metab Dispos 46(7):943-952, 2018. [doi:10.1124/dmd.118.080523](https://doi.org/10.1124/dmd.118.080523)
    - `Sudhamshu2006`: Sudhamshu KC, Matsutani S, Maruyama H, Akiike T, Saisho H. Doppler study of hepatic vein in cirrhotic patients: Correlation with liver dysfunction and hepatic hemodynamics. World J Gastroenterol 12(36):5853-5858, 2006. [doi:10.3748/wjg.v12.i36.5853](https://doi.org/10.3748/wjg.v12.i36.5853)
    - `VanBeers2001`: Van Beers BE, Leconte I, Materne R, Smith AM, Jamart J, Horsmans Y. Hepatic Perfusion Parameters in Chronic Liver Disease. AJR Am J Roentgenol 176(3):667-673, 2001. [doi:10.2214/ajr.176.3.1760667](https://doi.org/10.2214/ajr.176.3.1760667)
    - `Viani1992`: Viani A, Rizzo G, Carrai M, Pacifici GM. Interindividual variability in the concentrations of albumin and alpha-1-acid glycoprotein in patients with renal or liver disease, newborns and healthy subjects: implications for binding of drugs. Int J Clin Pharmacol Ther Toxicol 30(4):128-133, 1992. PMID 1572757.
    - `Woitas2000`: Woitas RP, Stoffel-Wagner B, Flommersfeld S, Poege U, Schiedermaier P, Klehr HU, Spengler U, Bidlingmaier F, Sauerbruch T. Correlation of Serum Concentrations of Cystatin C and Creatinine to Inulin Clearance in Liver Cirrhosis. Clin Chem 46(5):712-715, 2000. [doi:10.1093/clinchem/46.5.712](https://doi.org/10.1093/clinchem/46.5.712)
    - `Wong2001`: Wong F, Girgrah N, Graba J, Allidina Y, Liu P, Blendis L. The cardiac response to exercise in cirrhosis. Gut 49(2):268-275, 2001. [doi:10.1136/gut.49.2.268](https://doi.org/10.1136/gut.49.2.268)
    - `Xu2026`: Xu L, Jiang B, Ma Y, Pan J, Liang M, Du Q, Wu H, Cheng Z, Lu D. Comparison of Child-Pugh and MELD scores in predicting survival of hepatocellular carcinoma patients with splenomegaly undergoing TACE: a real-world study. Front Med 13:1804646, 2026. [doi:10.3389/fmed.2026.1804646](https://doi.org/10.3389/fmed.2026.1804646)
    - `Yang2020`: Yang M, Zhang Y, Zhao W, Cheng W, Wang H, Guo S. Evaluation of liver function using liver parenchyma, spleen and portal vein signal intensities during the hepatobiliary phase in Gd-EOB-D TPA-enhanced MRI. BMC Med Imaging 20(1):119, 2020. [doi:10.1186/s12880-020-00519-7](https://doi.org/10.1186/s12880-020-00519-7)
    - `Zhu1999`: Zhu JY, Leng XS, Dong N, Qi GY, Du RY. Measurement of liver volume and its clinical significance in cirrhotic portal hypertensive patients. World J Gastroenterol 5(6):525-526, 1999. [doi:10.3748/wjg.v5.i6.525](https://doi.org/10.3748/wjg.v5.i6.525)
<!-- end generated -->

Drozdzik 2021 gives the classes of every etiology and all diseases together, the mRNA as well as the protein; Drozdzik 2020 the transporters of the same cohorts, digitized from the box plots of its figures (see below); El-Khateeb 2021 the classes of a mixed cohort for 51 enzymes and transporters; Prasad 2018 end-stage alcoholic and HCV cirrhosis, which the authors assume to be Child-Pugh class C. Wang 2016 (the transporters of the livers of Prasad 2018, box plots only) is not yet in the database.

The physiology of cirrhosis comes from two groups of studies. The primary studies which Li 2015 tabulate for their meta-analysis of the physiology of the Child-Pugh classes (albumin, hematocrit, cardiac index, the flows of the portal vein and the hepatic artery, the perfusion and the volume of the liver) are entered with `via = "Li2015"`: the numbers as the supplement of Li 2015 gives them, the primary publication not read; the healthy groups which Li 2015 impute for a study without one are not taken. The open-access cohorts with the MELD score by Child-Pugh class (Mesropyan 2022, Hong 2025, Kim 2020, Yang 2020, Xu 2026) give the MELD, albumin, creatinine and, in Mesropyan 2022, the bilirubin, INR, hematocrit and the extracellular volume fraction of the liver.

Numbers which a source gives only in a figure are digitized from the image of the figure in the PDF: `scripts/database/digitize_boxplot.py` reads the box plots of R (median, quartiles, the smallest and largest value and the mean as `+`) of a configuration of the panels, every panel calibrated by two ticks of its y axis (a broken axis by the ticks of each segment), and the note of every record gives the resolution of a pixel. The digitized means of Drozdzik 2020 reproduce the changes the text states (`tests/test_database.py`).

## Integration into a consensus

`pkpdsim.database.integration.consensus(entity, condition)` integrates the studies which report a condition and its control into the ratio of the condition to the control:

1. The location of a group on the log scale, the log of its median (the geometric mean of log-normal data, so the means and the medians of different studies estimate the same ratio), with its standard error from the statistics its source reports, in this order: the median with its 95 % confidence interval, `ln(median)` with the standard error `(ln(upper) - ln(lower))/(2*1.96)`; the median with its quartiles, `sigma = (ln(q3) - ln(q1))/1.349` and the standard error of the log of a median `sqrt(pi/2)*sigma/sqrt(n)`; the median with the mean and SD, the log-normal SD of the log `sigma = sqrt(ln(1 + CV^2))` and the same standard error; the mean and SD alone, the geometric mean `ln(mean) - sigma^2/2` with the standard error `sigma/sqrt(n)`. A group without a spread takes the `sigma` of its study and kind of measurement (the median over its records), marked as imputed; a group whose median or mean is 0 (below the limit of quantification) has no log and is excluded with the reason. Data near the limit of quantification have no log-normal location: a comparison in which a group has a first quartile or a median of 0 or less than 75 % of its samples detected is the ratio of the arithmetic means (`ln(mean)` with the standard error `CV/sqrt(n)`), and one in which a group has a third quartile of 0 or less than 25 % of its samples detected is excluded.
2. The estimate of a study: its group of the condition without an etiology (or `mixed`) if it has one, else every group of the condition, strata of disjoint subjects pooled by their inverse variance; minus the control group of the study; the variance is the sum of both.
3. The consensus of the studies: the random-effects meta-analysis of DerSimonian and Laird, whose between-study variance `tau^2` widens the consensus by the heterogeneity of the studies (their populations, etiologies and methods); a single study has no `tau^2`.

The consensus is a log-normal factor `exp(mu)` with the 95 % confidence interval of the mean ratio and the uncertainty of the ratio of a new population, `sqrt(se^2 + tau^2)`. `integration.factor` takes the protein abundance and falls back to the mRNA expression only for an entity without a protein abundance of the condition.

```python
from pkpdsim.database import integration

c = integration.factor("CYP3A4", {"condition": "cirrhosis", "child_pugh": "C"})
print(c.text())
print(f"tau^2 {c.tau2:.3f}, prediction interval {c.prediction_interval}")
```

The levels of a quantity, the MELD score of a class, are integrated on the linear scale (`integration.level(entity, condition)`): the mean of a group with its standard error (the SE, or the SD over the root of `n`), or the median with the standard error of its quartiles; the groups of a study pooled by their inverse variance and the studies by the same random-effects meta-analysis.

```python
from pkpdsim.database import integration

meld = integration.level("meld", {"condition": "cirrhosis", "child_pugh": "B"})
print(meld.text())
```

## Liver cirrhosis

The factors of cirrhosis on the enzymes and transporters of the liver of the PBPK library (`CP_severity`, see [Impairment](pbpk/impairment.md)) are the consensus of each Child-Pugh class. The table gives the consensus of every entity of the studies:

<!-- generated: database-cirrhosis (scripts/pbpk_tables.py), do not edit -->
*Generated from the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| entity | A | B | C |
| --- | --- | --- | --- |
| `ABCB1` | 1.317 [0.6709-2.586] (1) | 1.174 [0.3323-4.148] (2) | 1.224 [0.2624-5.712] (2) |
| `ABCB11` | 0.6595 [0.3755-1.158] (2) | 0.3745 [0.2612-0.5371] (2) | 0.3364 [0.2152-0.5257] (2) |
| `ABCB4` | 0.7966 [0.1461-4.344] (1) | 0.4068 [0.1077-1.536] (1) | 0.339 [0.08424-1.364] (1) |
| `ABCC1` | 1.462 [0.9534-2.241] (1) (mRNA) | 1.54 [1.171-2.027] (1) (mRNA) | 1.708 [1.221-2.391] (1) (mRNA) |
| `ABCC2` | 0.2696 [0.08923-0.8143] (2) | 0.2959 [0.1799-0.4867] (2) | 0.1771 [0.122-0.2572] (2) |
| `ABCC3` | 0.474 [0.2742-0.8191] (2) | 0.4243 [0.3101-0.5807] (2) | 0.5212 [0.3782-0.7182] (2) |
| `ABCC4` | 0.3636 [0.03708-3.566] (1) | 0.6364 [0.0282-14.36] (1) | 0.2727 [0.02284-3.257] (1) |
| `ABCC6` | 0.6548 [0.2734-1.568] (1) | 0.2976 [0.145-0.6108] (1) | 0.2381 [0.1497-0.3786] (1) |
| `ABCG2` | 0.75 [0.1188-4.736] (1) | 1.25 [0.2579-6.059] (1) | 1 [0.2265-4.416] (1) |
| `ADH1A` |  |  | 0.3035 [0.2466-0.3735] (1) |
| `ADH1B` |  |  | 0.2841 [0.2303-0.3506] (1) |
| `ADH1C` |  |  | 0.2886 [0.2377-0.3505] (1) |
| `ALDH1A1` |  |  | 0.5949 [0.5415-0.6534] (1) |
| `ATP1A1` | 1.266 [0.5218-3.071] (1) | 0.5851 [0.3982-0.8598] (1) | 0.5266 [0.3473-0.7985] (1) |
| `CES1` | 0.5288 [0.2033-1.376] (1) | 0.2765 [0.1681-0.4549] (1) | 0.168 [0.05106-0.5526] (2) |
| `CES2` | 0.6257 [0.313-1.251] (1) | 0.4795 [0.2223-1.034] (1) | 0.4721 [0.1329-1.676] (2) |
| `CYP1A1` | 0.6149 [0.2785-1.357] (1) | 0.5711 [0.2827-1.154] (1) | 0.5042 [0.2264-1.123] (1) |
| `CYP1A2` | 0.481 [0.2853-0.811] (2) | 0.536 [0.336-0.8548] (2) | 0.1869 [0.14-0.2494] (3) |
| `CYP2A6` | 0.4408 [0.1969-0.9869] (1) | 0.4502 [0.206-0.9842] (1) | 0.3141 [0.2259-0.4368] (2) |
| `CYP2B6` | 0.7415 [0.5195-1.058] (2) | 0.6462 [0.4221-0.9892] (2) | 0.4604 [0.2971-0.7136] (2) |
| `CYP2C18` | 1.092 [0.3638-3.278] (1) | 0.4598 [0.1291-1.637] (1) | 0.4023 [0.1522-1.063] (1) |
| `CYP2C19` | 0.793 [0.5669-1.109] (1) | 0.1162 [0.002039-6.626] (2) | 0.2467 [0.01483-4.104] (2) |
| `CYP2C8` | 0.6561 [0.4722-0.9115] (2) | 0.3706 [0.2555-0.5375] (2) | 0.3428 [0.232-0.5065] (3) |
| `CYP2C9` | 0.9053 [0.6564-1.249] (2) | 0.7152 [0.3847-1.33] (2) | 0.415 [0.2853-0.6036] (3) |
| `CYP2D6` | 0.824 [0.5271-1.288] (2) | 0.7194 [0.2987-1.733] (2) | 0.5954 [0.415-0.8542] (3) |
| `CYP2E1` | 0.5138 [0.3773-0.6995] (2) | 0.4192 [0.306-0.5744] (2) | 0.2703 [0.1606-0.4548] (3) |
| `CYP2J2` | 0.4737 [0.1741-1.289] (1) | 0.2895 [0.1013-0.8274] (1) | 0.1 [0.03177-0.3148] (1) |
| `CYP3A4` | 0.5981 [0.3464-1.033] (2) | 0.3449 [0.1805-0.6591] (2) | 0.195 [0.101-0.3765] (3) |
| `CYP3A5` | 0.5281 [0.03754-7.427] (1) | 0.06068 [0.008073-0.4561] (1) | 0.5778 [0.04452-7.5] (1) |
| `CYP4F2` | 0.5732 [0.1958-1.678] (1) | 0.2878 [0.11-0.7532] (1) | 0.1683 [0.06503-0.4355] (1) |
| `EPHX1` | 0.5893 [0.3421-1.015] (1) | 0.3363 [0.2326-0.4861] (1) | 0.1467 [0.09512-0.2261] (1) |
| `FMO3` | 0.6261 [0.4236-0.9253] (1) | 0.4339 [0.2696-0.6983] (1) | 0.1911 [0.11-0.332] (1) |
| `FMO5` | 0.6944 [0.2549-1.891] (1) | 0.5036 [0.1679-1.51] (1) | 0.2187 [0.0808-0.5918] (1) |
| `MGST1` | 0.8407 [0.4791-1.475] (1) | 0.5641 [0.3731-0.8529] (1) | 0.3012 [0.1836-0.4941] (1) |
| `MGST3` | 0.4884 [0.2201-1.083] (1) | 0.4329 [0.2597-0.7219] (1) | 0.2791 [0.1642-0.4744] (1) |
| `POR` | 0.4808 [0.1252-1.846] (1) | 0.5265 [0.2358-1.176] (1) | 0.416 [0.3334-0.5191] (2) |
| `SLC10A1` | 0.626 [0.4217-0.9294] (2) | 0.4598 [0.2772-0.7626] (2) | 0.2351 [0.1527-0.3619] (2) |
| `SLC10A2` | 0.5 [0.08938-2.797] (1) | 0.75 [0.1524-3.69] (1) | 0.5 [0.09806-2.55] (1) |
| `SLC16A1` | 0.9906 [0.4375-2.243] (1) | 0.4688 [0.3101-0.7086] (1) | 0.3438 [0.2304-0.5129] (1) |
| `SLC22A1` | 0.4822 [0.3529-0.659] (2) | 0.3041 [0.1845-0.5009] (2) | 0.2795 [0.1533-0.5096] (2) |
| `SLC22A11` | 0.5 [0.2258-1.107] (1) | 0.2941 [0.1223-0.7073] (1) | 0.2353 [0.09607-0.5763] (1) |
| `SLC22A3` | 0.9468 [0.3491-2.568] (1) | 0.4894 [0.2864-0.8361] (1) | 0.5174 [0.3377-0.7927] (2) |
| `SLC22A7` | 0.9459 [0.4753-1.882] (2) | 0.5637 [0.1898-1.675] (2) | 0.3466 [0.1532-0.7842] (2) |
| `SLCO1A2` | 0.03333 [0.002647-0.4197] (1) | 0.1667 [0.0256-1.085] (1) | 0.06667 [0.004722-0.9412] (1) |
| `SLCO1B1` | 0.5187 [0.3636-0.7398] (2) | 0.4601 [0.3181-0.6653] (2) | 0.3746 [0.2395-0.5859] (2) |
| `SLCO1B3` | 0.9863 [0.3982-2.443] (1) | 0.3484 [0.122-0.9947] (2) | 0.2445 [0.03809-1.57] (2) |
| `SLCO2B1` | 0.6085 [0.3385-1.094] (2) | 0.3714 [0.242-0.57] (2) | 0.2638 [0.1657-0.4201] (2) |
| `UGT1A1` | 1.49 [1.046-2.123] (2) | 1.204 [0.8071-1.796] (2) | 0.9976 [0.6774-1.469] (2) |
| `UGT1A3` | 0.4279 [0.0584-3.135] (2) | 0.6337 [0.2801-1.434] (2) | 0.4518 [0.1047-1.949] (2) |
| `UGT1A4` | 0.5495 [0.4106-0.7355] (1) | 0.66 [0.2737-1.591] (1) | 0.2099 [0.1006-0.4381] (2) |
| `UGT1A6` | 0.7445 [0.4052-1.368] (1) | 0.4365 [0.2648-0.7197] (1) | 0.1174 [0.03049-0.4516] (2) |
| `UGT1A9` | 1.356 [0.4962-3.704] (1) | 0.4759 [0.2015-1.124] (1) | 0.2799 [0.1275-0.6143] (1) |
| `UGT2B15` | 0.6701 [0.2928-1.533] (2) | 0.6507 [0.464-0.9124] (2) | 0.4335 [0.2767-0.6793] (3) |
| `UGT2B4` | 0.6917 [0.3969-1.205] (1) | 0.404 [0.2666-0.6122] (1) | 0.2464 [0.1529-0.3973] (1) |
| `UGT2B7` | 0.7633 [0.5464-1.066] (2) | 0.4251 [0.2486-0.7271] (2) | 0.2091 [0.1223-0.3574] (3) |

The ratio of the class to the controls with its 95 % confidence interval and the number of studies, the protein abundance per amount of liver tissue, `(mRNA)` the mRNA expression of an entity without a protein abundance.

??? note "Sources"

    - `Drozdzik2020`: Drozdzik M, Szelag-Pieniek S, Post M, Zeair S, Wrzesinski M, Kurzawski M, Prieto J, Oswald S. Protein abundance of hepatic drug transporters in patients with different forms of liver damage. Clin Pharmacol Ther 107(5):1138-1148, 2020. [doi:10.1002/cpt.1717](https://doi.org/10.1002/cpt.1717)
    - `Drozdzik2021`: Drozdzik M, Lapczuk-Romanska J, Wenzel C, Szelag-Pieniek S, Post M, Skalski Ł, Kurzawski M, Oswald S. Gene expression and protein abundance of hepatic drug metabolizing enzymes in liver pathology. Pharmaceutics 13(9):1334, 2021. [doi:10.3390/pharmaceutics13091334](https://doi.org/10.3390/pharmaceutics13091334)
    - `ElKhateeb2021`: El-Khateeb E, Achour B, Al-Majdoub ZM, Barber J, Rostami-Hodjegan A. Non-uniformity of changes in drug-metabolizing enzymes and transporters in liver cirrhosis: implications for drug dosage adjustment. Mol Pharm 18(9):3563-3577, 2021. [doi:10.1021/acs.molpharmaceut.1c00462](https://doi.org/10.1021/acs.molpharmaceut.1c00462)
    - `Prasad2018`: Prasad B, Bhatt DK, Johnson K, Chapa R, Chu X, Salphati L, Xiao G, Lee C, Hop CE, Mathias A, Lai Y, Liao M, Humphreys WG, Kumer SC, Unadkat JD. Abundance of phase 1 and 2 drug-metabolizing enzymes in alcoholic and hepatitis C cirrhotic livers: a quantitative targeted proteomics study. Drug Metab Dispos 46(7):943-952, 2018. [doi:10.1124/dmd.118.080523](https://doi.org/10.1124/dmd.118.080523)
<!-- end generated -->

### Physiology and the MELD score

The factors of the Child-Pugh classes on the albumin, the hematocrit and the cardiac index of the library are the consensus of the ratios of the studies with a control group (`pkpdsim.pbpk.impairment.child_pugh_factor`); the MELD score of each class is the anchor of the mapping of a MELD score to the severity `CP_severity` (`cp_severity_of_meld`, see [Impairment](pbpk/impairment.md#the-meld-score)). The table gives the level of every quantity and its ratio to the controls:

<!-- generated: database-physiology (scripts/pbpk_tables.py), do not edit -->
*Generated from the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| quantity | controls | A | B | C |
| --- | --- | --- | --- | --- |
| `albumin` [g/L] | 46.34 [43.77-48.92] (5) | 40.04 [38.16-41.93] (9)<br>0.8645 [0.8077-0.9253] (4) | 31.73 [31.04-32.43] (10)<br>0.6798 [0.6327-0.7304] (5) | 24.83 [23.16-26.49] (10)<br>0.4892 [0.4303-0.5562] (5) |
| `hematocrit` [%] | 41.83 [40.74-42.93] (3) | 37.9 [37.19-38.62] (3)<br>0.9232 [0.892-0.9556] (2) | 30.7 [30.53-30.87] (1)<br>0.7411 [0.7124-0.771] (1) | 30.73 [24.46-37] (2)<br>0.7232 [0.6123-0.8541] (2) |
| `bilirubin` [mg/dL] | 0.78 [0.5669-0.9931] (1) | 1.02 [0.8352-1.205] (1)<br>1.408 [1.04-1.906] (1) | 1.89 [0.9753-2.805] (1)<br>1.685 [1.087-2.612] (1) | 2.97 [1.5-4.44] (1)<br>3.626 [2.187-6.013] (1) |
| `creatinine` [mg/dL] | 0.7107 [0.6767-0.7448] (2) | 0.8023 [0.7782-0.8264] (5)<br>1.192 [0.8475-1.676] (2) | 0.8524 [0.743-0.9618] (5)<br>1.375 [0.7053-2.681] (2) | 0.9902 [0.7604-1.22] (5)<br>1.667 [0.6626-4.195] (2) |
| `inr` [1] | 1.03 [0.9799-1.08] (1) | 1.12 [1.075-1.165] (1)<br>1.089 [1.022-1.159] (1) | 1.2 [1.131-1.269] (1)<br>1.157 [1.073-1.247] (1) | 1.54 [1.122-1.958] (1)<br>1.39 [1.066-1.812] (1) |
| `meld` [1] | 6.3 [6.007-6.593] (1) | 7.909 [7.554-8.265] (5)<br>1.233 [1.068-1.425] (2) | 12.78 [10.8-14.76] (5)<br>1.661 [1.452-1.9] (2) | 20.5 [17.01-23.99] (5)<br>2.555 [2.385-2.738] (2) |
| `cardiac_index` [L/min/m2] | 3.649 [2.868-4.43] (3) | 4.103 [3.227-4.98] (3)<br>1.099 [1.013-1.192] (3) | 4.477 [3.622-5.333] (3)<br>1.21 [1.103-1.327] (3) | 4.988 [3.974-6.001] (3)<br>1.36 [1.249-1.48] (3) |
| `portal_vein_flow` [mL/min] | 722.8 [147-1299] (3) | 725.3 [227.3-1223] (3)<br>0.9772 [0.863-1.107] (3) | 637.8 [355.4-920.3] (3)<br>0.9354 [0.5437-1.609] (3) | 531.4 [355.5-707.2] (3)<br>0.7875 [0.3781-1.64] (3) |
| `portal_perfusion` [mL/min/dL] | 70.05 [38.02-102.1] (2) | 36.67 [6.6-66.74] (2)<br>0.4768 [0.3353-0.678] (2) | 28.21 [7.049-49.38] (2)<br>0.3334 [0.2281-0.4874] (2) | 2.391 [-1.73-6.512] (2)<br>0.04108 [0.01085-0.1556] (2) |
| `splenic_vein_flow` [mL/min] | 447.3 [397.7-496.9] (1) | 535.8 [532.5-539.1] (1)<br>1.232 [1.104-1.374] (1) | 652 [558.1-745.9] (1)<br>1.449 [1.212-1.733] (1) | 689 [613.3-764.7] (1)<br>1.534 [1.316-1.789] (1) |
| `hepatic_artery_flow` [mL/min] | 99.4 [90.38-108.4] (1) | 153 [129.6-176.4] (1)<br>1.404 [1.185-1.662] (1) | 161 [138.3-183.7] (1)<br>1.558 [1.324-1.834] (1) | 195 [155.8-234.2] (1)<br>1.801 [1.461-2.22] (1) |
| `hepatic_arterial_perfusion` [mL/min/dL] | 10.82 [3.112-18.53] (2) | 14.19 [8.855-19.53] (2)<br>1.129 [0.7204-1.77] (2) | 22.14 [15.17-29.1] (2)<br>2.066 [1.389-3.073] (2) | 33.93 [22.98-44.88] (2)<br>3.262 [2.152-4.945] (2) |
| `liver_volume_per_bw` [cm3/kg] | 20.87 [19.78-21.95] (2) | 18.76 [17.81-19.71] (3)<br>0.8982 [0.8213-0.9823] (2) | 16.73 [13.31-20.16] (3)<br>0.7854 [0.5582-1.105] (2) | 14.6 [11.79-17.41] (1)<br>0.6655 [0.5446-0.8132] (1) |
| `liver_volume_per_height` [cm3/m] | 629.8 [581.9-677.7] (1) | <br>1.04 [0.8645-1.252] (1) | 533.1 [505.9-560.3] (1)<br>0.8615 [0.7867-0.9435] (1) | 393.8 [363.9-423.7] (1)<br>0.6345 [0.5703-0.7059] (1) |
| `hepatocyte_volume_per_bw` [cm3/kg] | 17.9 [16.14-19.66] (1) | 15.59 [14.57-16.61] (2)<br>0.8601 [0.7607-0.9723] (1) | 14.63 [12.13-17.13] (2)<br>0.8196 [0.661-1.016] (1) |  |
| `functional_liver_volume_per_bw` [cm3/kg] | 26.5 [24.36-28.64] (1) | 23.7 [22.14-25.26] (1)<br>0.8884 [0.8012-0.985] (1) | 24.1 [19.79-28.41] (1)<br>0.8897 [0.7335-1.079] (1) |  |
| `hepatic_ecv` [%] | 25.9 [24.02-27.78] (1) | 33.3 [31.04-35.56] (1)<br>1.284 [1.164-1.418] (1) | 39.6 [37.9-41.3] (1)<br>1.54 [1.416-1.675] (1) | 52.8 [52.02-53.58] (1)<br>2.069 [1.922-2.227] (1) |

Per cell the level, the mean or median of the groups pooled over the studies on the linear scale, and below it the ratio to the controls of the studies with a control group, pooled on the log scale; each with its 95 % confidence interval and the number of studies.

??? note "Sources"

    - `Alberino2001`: Alberino F, Gatta A, Amodio P, Merkel C, Di Pascoli L, Boffo G, Caregaro L. Nutrition and survival in patients with liver cirrhosis. Nutrition 17(6):445-450, 2001. [doi:10.1016/s0899-9007(01)00521-4](https://doi.org/10.1016/s0899-9007(01)00521-4)
    - `Annet2003`: Annet L, Materne R, Danse E, Jamart J, Horsmans Y, Van Beers BE. Hepatic Flow Parameters Measured with MR Imaging and Doppler US: Correlations with Degree of Cirrhosis and Portal Hypertension. Radiology 229(2):409-414, 2003. [doi:10.1148/radiol.2292021128](https://doi.org/10.1148/radiol.2292021128)
    - `Boker2000`: Böker KH, Franzen A, Wrobel M, Bahr MJ, Tietge U, Manns MP. Regulation of hepatic blood flow in patients with liver cirrhosis and after liver transplantation. Pathophysiology 6(4):251-260, 2000. [doi:10.1016/s0928-4680(99)00024-3](https://doi.org/10.1016/s0928-4680(99)00024-3)
    - `Chawla1998`: Chawla Y, Radha NS, Jang KD, Dilawari B. Portal Hemodynamics by Duplex Doppler Sonography in Different Grades of Cirrhosis. Dig Dis Sci 43(2):354-357, 1998. [doi:10.1023/a:1018814624307](https://doi.org/10.1023/a:1018814624307)
    - `Dincer2005`: Dincer D, Besisk F, Demirkol O, Demir K, Kaymakoglu S, Cakaloglu Y, Okten A. Relationships between hemodynamic alterations and Child-Pugh Score in patients with cirrhosis. Hepatogastroenterology 52(65):1521-1525, 2005. PMID 16201110.
    - `Froomes1999`: Froomes PR, Morgan DJ, Smallwood RA, Angus PW. Comparative effects of oxygen supplementation on theophylline and acetaminophen clearance in human cirrhosis. Gastroenterology 116(4):915-920, 1999. [doi:10.1016/s0016-5085(99)70075-2](https://doi.org/10.1016/s0016-5085(99)70075-2)
    - `Henriksen2001`: Henriksen JH, Bendtsen F, Møller S. Normal red cell cardiac output in the hyperkinetic syndrome of alcoholic cirrhosis. J Hepatol 34(5):782-783, 2001. [doi:10.1016/s0168-8278(01)00019-8](https://doi.org/10.1016/s0168-8278(01)00019-8)
    - `Hong2025`: Hong SJ, Kim JE, Cho JM, Choi HC, Park MJ, Choi HY, Shin HS, Won JH, Yang W, Kim HO. Dual Energy CT-Derived Liver Extracellular Volume Fraction for Assessing Liver Functional Reserve in Patients with Liver Cirrhosis. Medicina 61(9):1561, 2025. [doi:10.3390/medicina61091561](https://doi.org/10.3390/medicina61091561)
    - `Iwao1997`: Iwao T, Oho K, Sakai T, Tayama C, Sato M, Nakano R, Yamawaki M, Toyonaga A, Tanikawa K. Splanchnic and extrasplanchnic arterial hemodynamics in patients with cirrhosis. J Hepatol 27(5):817-823, 1997. [doi:10.1016/s0168-8278(97)80318-2](https://doi.org/10.1016/s0168-8278(97)80318-2)
    - `Kayacetin2004`: Kayacetin E, Efe D, Doğan C. Portal and splenic hemodynamics in cirrhotic patients: relationship between esophageal variceal bleeding and the severity of hepatic failure. J Gastroenterol 39(7):661-667, 2004. [doi:10.1007/s00535-003-1362-x](https://doi.org/10.1007/s00535-003-1362-x)
    - `Kim2020`: Kim TH, Yun SG, Choi J, Goh HG, Lee HA, Yim SY, Choi SJ, Lee YS, Yoon EL, Jung YK, Seo YS, Kim JH, Yim HJ, Yeon JE, Byun KS, Um SH. Differential Impact of Serum 25-Hydroxyvitamin D3 Levels on the Prognosis of Patients with Liver Cirrhosis According to MELD and Child-Pugh Scores. J Korean Med Sci 35(19):e129, 2020. [doi:10.3346/jkms.2020.35.e129](https://doi.org/10.3346/jkms.2020.35.e129)
    - `Kwon2001`: Kwon HA, Matsui Y, Ha-Kawa SK, Kamiyama Y. Functional Hepatic Volume Measured by Technetium-99M-Galactosyl-Human Serum Albumin Liver Scintigraphy: Comparison Between Hepatocyte Volume and Liver Volume by Computed Tomography. Am J Gastroenterol 96(2):541-546, 2001. [doi:10.1111/j.1572-0241.2001.03556.x](https://doi.org/10.1111/j.1572-0241.2001.03556.x)
    - `Lin1998`: Lin XZ, Sun YN, Liu YH, Sheu BS, Cheng BN, Chen CY, Tsai HM, Shen CL. Liver volume in patients with or without chronic liver diseases. Hepatogastroenterology 45(22):1069-1074, 1998. PMID 9756008.
    - `Liu2006`: Liu T, Wong W, Hou M, Lin H, Chang F, Lee S. Hemorheology in patients with liver cirrhosis: Special emphasis on its relation to severity of esophageal variceal bleeding. J Gastroenterol Hepatol 21(5):908-913, 2006. [doi:10.1111/j.1440-1746.2006.04266.x](https://doi.org/10.1111/j.1440-1746.2006.04266.x)
    - `Matsui1996`: Matsui Y, Tu W, Kitade H, Nakagawa A, Kamiya T, Kwon A, Uetsuji S, Kamiyama Y. Hepatocyte volume as an indicator of hepatic functional reserve in cirrhotic patients with liver tumours. J Gastroenterol Hepatol 11(6):540-545, 1996. [doi:10.1111/j.1440-1746.1996.tb01699.x](https://doi.org/10.1111/j.1440-1746.1996.tb01699.x)
    - `Mesropyan2022`: Mesropyan N, Kupczyk PA, Dold L, Praktiknjo M, Chang J, Isaak A, Endler C, Kravchenko D, Bischoff LM, Sprinkart AM, Pieper CC, Kuetting D, Jansen C, Attenberger UI, Luetkens JA. Assessment of liver cirrhosis severity with extracellular volume fraction MRI. Sci Rep 12(1):9422, 2022. [doi:10.1038/s41598-022-13340-9](https://doi.org/10.1038/s41598-022-13340-9)
    - `Sudhamshu2006`: Sudhamshu KC, Matsutani S, Maruyama H, Akiike T, Saisho H. Doppler study of hepatic vein in cirrhotic patients: Correlation with liver dysfunction and hepatic hemodynamics. World J Gastroenterol 12(36):5853-5858, 2006. [doi:10.3748/wjg.v12.i36.5853](https://doi.org/10.3748/wjg.v12.i36.5853)
    - `VanBeers2001`: Van Beers BE, Leconte I, Materne R, Smith AM, Jamart J, Horsmans Y. Hepatic Perfusion Parameters in Chronic Liver Disease. AJR Am J Roentgenol 176(3):667-673, 2001. [doi:10.2214/ajr.176.3.1760667](https://doi.org/10.2214/ajr.176.3.1760667)
    - `Viani1992`: Viani A, Rizzo G, Carrai M, Pacifici GM. Interindividual variability in the concentrations of albumin and alpha-1-acid glycoprotein in patients with renal or liver disease, newborns and healthy subjects: implications for binding of drugs. Int J Clin Pharmacol Ther Toxicol 30(4):128-133, 1992. PMID 1572757.
    - `Woitas2000`: Woitas RP, Stoffel-Wagner B, Flommersfeld S, Poege U, Schiedermaier P, Klehr HU, Spengler U, Bidlingmaier F, Sauerbruch T. Correlation of Serum Concentrations of Cystatin C and Creatinine to Inulin Clearance in Liver Cirrhosis. Clin Chem 46(5):712-715, 2000. [doi:10.1093/clinchem/46.5.712](https://doi.org/10.1093/clinchem/46.5.712)
    - `Wong2001`: Wong F, Girgrah N, Graba J, Allidina Y, Liu P, Blendis L. The cardiac response to exercise in cirrhosis. Gut 49(2):268-275, 2001. [doi:10.1136/gut.49.2.268](https://doi.org/10.1136/gut.49.2.268)
    - `Xu2026`: Xu L, Jiang B, Ma Y, Pan J, Liang M, Du Q, Wu H, Cheng Z, Lu D. Comparison of Child-Pugh and MELD scores in predicting survival of hepatocellular carcinoma patients with splenomegaly undergoing TACE: a real-world study. Front Med 13:1804646, 2026. [doi:10.3389/fmed.2026.1804646](https://doi.org/10.3389/fmed.2026.1804646)
    - `Yang2020`: Yang M, Zhang Y, Zhao W, Cheng W, Wang H, Guo S. Evaluation of liver function using liver parenchyma, spleen and portal vein signal intensities during the hepatobiliary phase in Gd-EOB-D TPA-enhanced MRI. BMC Med Imaging 20(1):119, 2020. [doi:10.1186/s12880-020-00519-7](https://doi.org/10.1186/s12880-020-00519-7)
    - `Zhu1999`: Zhu JY, Leng XS, Dong N, Qi GY, Du RY. Measurement of liver volume and its clinical significance in cirrhotic portal hypertensive patients. World J Gastroenterol 5(6):525-526, 1999. [doi:10.3748/wjg.v5.i6.525](https://doi.org/10.3748/wjg.v5.i6.525)
<!-- end generated -->

## Adding a study

A study is added as its folder: `study.toml` and one TSV file per measurement, the key of its source added to `sources.toml`. Extract the tables of the source as they are printed, with the location in the source; check the extraction by the statistics which constrain each other (the CV and the SD over the mean, the median and the mean within the range, the sum of the strata); write a note for every number which disagrees, and state an assumption of a group (e.g. its class) as its `assumption`. Then run the tests of `tests/test_database.py` and `uv run python scripts/pbpk_tables.py`, which writes the tables of this page and of the impairment of the library. A study changes the consensus, and with it the SBML of every model which uses an entity of the study; the release notes name such a change.
