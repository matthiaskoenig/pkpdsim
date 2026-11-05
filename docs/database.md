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

A study is a folder named by the key of its source. `study.toml` describes the tissue, the control group, the groups (`condition` `control` or `cirrhosis`, the Child-Pugh class, the etiology, an `assumption` the source does not state) and the measurements (`kind` `protein`, the abundance per amount of tissue, or `mrna`, the expression; the unit, the method and where the values are in the source). Every measurement is a TSV file of records, one row per entity (a gene symbol of HGNC, e.g. `CYP3A4` or `SLCO1B1`, with the name of the source) and group, with the statistics the source reports: `n`, `n_detected` (the samples above the limit of quantification), `mean`, `sd`, `se`, `cv_percent`, `median`, `min`, `max` and `ci_low`/`ci_high` (the 95 % confidence interval of the median).

The records are the numbers of the source as printed, never derived ones. A statistic the source does not report is empty. A number which disagrees with the others of its row (a mean outside its range, a CV which is not the SD over the mean, a sum of the classes which is not their group) is kept as printed and has a note which says what the other numbers imply; `pkpdsim.database.studies.check` reports a disagreement without a note, and the tests keep it empty.

<!-- generated: database-studies (scripts/pbpk_tables.py), do not edit -->
*Generated from `pkpdsim/database/studies` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| study | measurement | unit | entities | records | groups | location |
| --- | --- | --- | ---: | ---: | ---: | --- |
| `Drozdzik2021` | protein abundance per amount of tissue | fmol/mg tissue | 13 | 299 | 23 | Table S3 |
| `Drozdzik2021` | mRNA expression | relative to the mean of the controls (2^-ddCt) | 13 | 299 | 23 | Table S2 |
| `ElKhateeb2021` | protein abundance per amount of tissue | pmol/g liver tissue | 51 | 201 | 4 | Table S8 of the supporting information |
| `Prasad2018` | protein abundance per amount of tissue | pmol/g tissue | 18 | 54 | 3 | Table 3 |

??? note "Sources"

    - `Drozdzik2021`: Drozdzik M, Lapczuk-Romanska J, Wenzel C, Szelag-Pieniek S, Post M, Skalski Ł, Kurzawski M, Oswald S. Gene expression and protein abundance of hepatic drug metabolizing enzymes in liver pathology. Pharmaceutics 13(9):1334, 2021. [doi:10.3390/pharmaceutics13091334](https://doi.org/10.3390/pharmaceutics13091334)
    - `ElKhateeb2021`: El-Khateeb E, Achour B, Al-Majdoub ZM, Barber J, Rostami-Hodjegan A. Non-uniformity of changes in drug-metabolizing enzymes and transporters in liver cirrhosis: implications for drug dosage adjustment. Mol Pharm 18(9):3563-3577, 2021. [doi:10.1021/acs.molpharmaceut.1c00462](https://doi.org/10.1021/acs.molpharmaceut.1c00462)
    - `Prasad2018`: Prasad B, Bhatt DK, Johnson K, Chapa R, Chu X, Salphati L, Xiao G, Lee C, Hop CE, Mathias A, Lai Y, Liao M, Humphreys WG, Kumer SC, Unadkat JD. Abundance of phase 1 and 2 drug-metabolizing enzymes in alcoholic and hepatitis C cirrhotic livers: a quantitative targeted proteomics study. Drug Metab Dispos 46(7):943-952, 2018. [doi:10.1124/dmd.118.080523](https://doi.org/10.1124/dmd.118.080523)
<!-- end generated -->

Drozdzik 2021 gives the classes of every etiology and all diseases together, the mRNA as well as the protein; El-Khateeb 2021 the classes of a mixed cohort for 51 enzymes and transporters; Prasad 2018 end-stage alcoholic and HCV cirrhosis, which the authors assume to be Child-Pugh class C. Drozdzik 2020 (the transporters of the same cohorts as Drozdzik 2021) and Wang 2016 (the transporters of the livers of Prasad 2018, box plots only) are not yet in the database.

## Integration into a consensus

`pkpdsim.database.integration.consensus(entity, condition)` integrates the studies which report a condition and its control into the ratio of the condition to the control:

1. The location of a group on the log scale, the log of its median (the geometric mean of log-normal data, so the means and the medians of different studies estimate the same ratio), with its standard error from the statistics its source reports, in this order: the median with its 95 % confidence interval, `ln(median)` with the standard error `(ln(upper) - ln(lower))/(2*1.96)`; the median with the mean and SD, the log-normal SD of the log `sigma = sqrt(ln(1 + CV^2))` and the standard error of the log of a median `sqrt(pi/2)*sigma/sqrt(n)`; the mean and SD alone, the geometric mean `ln(mean) - sigma^2/2` with the standard error `sigma/sqrt(n)`. A group without a spread takes the `sigma` of its study and kind of measurement (the median over its records), marked as imputed; a group whose median or mean is 0 (below the limit of quantification) has no log and is excluded with the reason.
2. The estimate of a study: its group of the condition without an etiology (or `mixed`) if it has one, else every group of the condition, strata of disjoint subjects pooled by their inverse variance; minus the control group of the study; the variance is the sum of both.
3. The consensus of the studies: the random-effects meta-analysis of DerSimonian and Laird, whose between-study variance `tau^2` widens the consensus by the heterogeneity of the studies (their populations, etiologies and methods); a single study has no `tau^2`.

The consensus is a log-normal factor `exp(mu)` with the 95 % confidence interval of the mean ratio and the uncertainty of the ratio of a new population, `sqrt(se^2 + tau^2)`. `integration.factor` takes the protein abundance and falls back to the mRNA expression only for an entity without a protein abundance of the condition.

```python
from pkpdsim.database import integration

c = integration.factor("CYP3A4", {"condition": "cirrhosis", "child_pugh": "C"})
print(c.text())
print(f"tau^2 {c.tau2:.3f}, prediction interval {c.prediction_interval}")
```

## Liver cirrhosis

The factors of cirrhosis on the enzymes and transporters of the liver of the PBPK library (`CP_severity`, see [Impairment](pbpk/impairment.md)) are the consensus of each Child-Pugh class. The table gives the consensus of every entity of the studies:

<!-- generated: database-cirrhosis (scripts/pbpk_tables.py), do not edit -->
*Generated from the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| entity | A | B | C |
| --- | --- | --- | --- |
| `ABCB1` | 1.317 [0.6709-2.586] (1) | 0.6098 [0.3734-0.9957] (1) | 0.5488 [0.3085-0.9762] (1) |
| `ABCB11` | 0.5455 [0.2566-1.159] (1) | 0.3636 [0.2416-0.5473] (1) | 0.2727 [0.159-0.4678] (1) |
| `ABCB4` | 0.7966 [0.1461-4.344] (1) | 0.4068 [0.1077-1.536] (1) | 0.339 [0.08424-1.364] (1) |
| `ABCC2` | 0.4558 [0.2163-0.9604] (1) | 0.3401 [0.1942-0.5956] (1) | 0.1837 [0.1233-0.2736] (1) |
| `ABCC3` | 0.8111 [0.2722-2.417] (1) | 0.5556 [0.2972-1.038] (1) | 0.5556 [0.2843-1.086] (1) |
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
| `CYP1A1` | 0.6149 [0.2785-1.357] (1) | 0.5888 [0.298-1.164] (1) | 0.5433 [0.2574-1.147] (1) |
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
| `CYP3A4` | 0.5981 [0.3464-1.033] (2) | 0.3449 [0.1805-0.6591] (2) | 0.1816 [0.1048-0.3146] (3) |
| `CYP3A5` | 0.5281 [0.03754-7.427] (1) | 0.06068 [0.008073-0.4561] (1) | 0.5778 [0.04452-7.5] (1) |
| `CYP4F2` | 0.5732 [0.1958-1.678] (1) | 0.2878 [0.11-0.7532] (1) | 0.1683 [0.06503-0.4355] (1) |
| `EPHX1` | 0.5893 [0.3421-1.015] (1) | 0.3363 [0.2326-0.4861] (1) | 0.1467 [0.09512-0.2261] (1) |
| `FMO3` | 0.6261 [0.4236-0.9253] (1) | 0.4339 [0.2696-0.6983] (1) | 0.1911 [0.11-0.332] (1) |
| `FMO5` | 0.6944 [0.2549-1.891] (1) | 0.5036 [0.1679-1.51] (1) | 0.2187 [0.0808-0.5918] (1) |
| `MGST1` | 0.8407 [0.4791-1.475] (1) | 0.5641 [0.3731-0.8529] (1) | 0.3012 [0.1836-0.4941] (1) |
| `MGST3` | 0.4884 [0.2201-1.083] (1) | 0.4329 [0.2597-0.7219] (1) | 0.2791 [0.1642-0.4744] (1) |
| `POR` | 0.4808 [0.1252-1.846] (1) | 0.5265 [0.2358-1.176] (1) | 0.416 [0.3334-0.5191] (2) |
| `SLC10A1` | 0.5634 [0.3181-0.9977] (1) | 0.4028 [0.1826-0.8888] (1) | 0.2352 [0.1308-0.4229] (1) |
| `SLC10A2` | 0.5 [0.08938-2.797] (1) | 0.75 [0.1524-3.69] (1) | 0.5 [0.09806-2.55] (1) |
| `SLC16A1` | 0.9906 [0.4375-2.243] (1) | 0.4688 [0.3101-0.7086] (1) | 0.3438 [0.2304-0.5129] (1) |
| `SLC22A1` | 0.4498 [0.3165-0.6393] (1) | 0.2307 [0.1455-0.3657] (1) | 0.1845 [0.08515-0.3999] (1) |
| `SLC22A11` | 0.5 [0.2258-1.107] (1) | 0.2941 [0.1223-0.7073] (1) | 0.2353 [0.09607-0.5763] (1) |
| `SLC22A3` | 0.9468 [0.3491-2.568] (1) | 0.4894 [0.2864-0.8361] (1) | 0.3617 [0.1713-0.7637] (1) |
| `SLC22A7` | 0.6463 [0.2745-1.521] (1) | 0.3401 [0.1843-0.6277] (1) | 0.2483 [0.1367-0.4511] (1) |
| `SLCO1A2` | 0.03333 [0.002647-0.4197] (1) | 0.1667 [0.0256-1.085] (1) | 0.06667 [0.004722-0.9412] (1) |
| `SLCO1B1` | 0.4545 [0.2915-0.7087] (1) | 0.4182 [0.2635-0.6636] (1) | 0.297 [0.1337-0.6597] (1) |
| `SLCO1B3` |  | 0.2074 [0.07668-0.561] (1) | 0.08889 [0.03223-0.2451] (1) |
| `SLCO2B1` | 0.7805 [0.3432-1.775] (1) | 0.3512 [0.2143-0.5756] (1) | 0.2634 [0.1524-0.4553] (1) |
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

    - `Drozdzik2021`: Drozdzik M, Lapczuk-Romanska J, Wenzel C, Szelag-Pieniek S, Post M, Skalski Ł, Kurzawski M, Oswald S. Gene expression and protein abundance of hepatic drug metabolizing enzymes in liver pathology. Pharmaceutics 13(9):1334, 2021. [doi:10.3390/pharmaceutics13091334](https://doi.org/10.3390/pharmaceutics13091334)
    - `ElKhateeb2021`: El-Khateeb E, Achour B, Al-Majdoub ZM, Barber J, Rostami-Hodjegan A. Non-uniformity of changes in drug-metabolizing enzymes and transporters in liver cirrhosis: implications for drug dosage adjustment. Mol Pharm 18(9):3563-3577, 2021. [doi:10.1021/acs.molpharmaceut.1c00462](https://doi.org/10.1021/acs.molpharmaceut.1c00462)
    - `Prasad2018`: Prasad B, Bhatt DK, Johnson K, Chapa R, Chu X, Salphati L, Xiao G, Lee C, Hop CE, Mathias A, Lai Y, Liao M, Humphreys WG, Kumer SC, Unadkat JD. Abundance of phase 1 and 2 drug-metabolizing enzymes in alcoholic and hepatitis C cirrhotic livers: a quantitative targeted proteomics study. Drug Metab Dispos 46(7):943-952, 2018. [doi:10.1124/dmd.118.080523](https://doi.org/10.1124/dmd.118.080523)
<!-- end generated -->

## Adding a study

A study is added as its folder: `study.toml` and one TSV file per measurement, the key of its source added to `sources.toml`. Extract the tables of the source as they are printed, with the location in the source; check the extraction by the statistics which constrain each other (the CV and the SD over the mean, the median and the mean within the range, the sum of the strata); write a note for every number which disagrees, and state an assumption of a group (e.g. its class) as its `assumption`. Then run the tests of `tests/test_database.py` and `uv run python scripts/pbpk_tables.py`, which writes the tables of this page and of the impairment of the library. A study changes the consensus, and with it the SBML of every model which uses an entity of the study; the release notes name such a change.
