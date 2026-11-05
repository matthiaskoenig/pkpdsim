# Impairment

Impaired organ function is a set of parameter vectors on the physiology, never a structurally different model: three severities of the top model set factors on the quantities of the healthy individual, the liver and the binding to plasma proteins. All severities 0 give the healthy model, every factor is 1. The numbers are those of `database/tables/impairment.toml`, every number with its source; a number which no publication read gives is a choice of the library (`ModelChoice`) with a note.

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/impairment.svg"
</div>

| condition | severity | classes | acts on |
|---|---|---|---|
| liver cirrhosis | `CP_severity` | 0 healthy, 1, 2, 3 the Child-Pugh classes A, B, C, interpolated linearly | functional liver `f_functional`, shunt `f_shunt`, hepatic artery and kidney flow, cardiac index, GFR, albumin and AGP (so fu), the hematocrit, the enzymes and transporters of the liver; a MELD score gives the severity ([The MELD score](#the-meld-score)) |
| chronic kidney disease | `f_renal_function` | the GFR relative to healthy; its GFR per 1.73 m2 (`eGFR`) gives the stage | GFR, the binding to albumin and AGP by their CKD profiles and to an unknown protein by the binding capacity of the plasma proteins (so fu); the hematocrit and the renal blood flow; kidney volume, gastric emptying and small intestinal transit (factor 1, see [Limitations](#limitations)); with the switch `CKD_hepatic = 1` hepatic CYP2D6 and OATP1B |
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
fu_<s> = 1 / (1 + fimp_albumin * fimp_albumin_ckd * (1 - fu_ref_<s>) / fu_ref_<s>)
```

The unbound fraction follows the binding (McNamara 2002): `fu = 1/(1 + f*(1 - fu_ref)/fu_ref)` with the binding factor `f`, the factor of the binding protein in cirrhosis (`fimp_albumin`, `fimp_agp`, 1 for an unknown protein) times the CKD factor of the binding: for albumin and AGP the CKD profile of the protein (`fimp_albumin_ckd`, `fimp_agp_ckd`, see [Binding in chronic kidney disease](#binding-in-chronic-kidney-disease)), for an unknown protein one over the factor of the binding capacity `fimp_plasma_protein`. fu of an albumin-bound substance rises in CKD, fu of an AGP-bound substance falls, and fu stays below 1. A substance whose binding is measured unchanged in renal impairment (apixaban: 92-94 % bound in healthy and ESRD subjects, US label) opts out with `Binding(fu, protein, ckd=False)`: the formula then has no CKD factor, `fu = 1/(1 + fimp_albumin * (1 - fu_ref)/fu_ref)`, so fu follows only the binding protein (cirrhosis still acts). With an unknown protein as well no factor is left and `fu = fu_ref`, the reference itself; an unknown protein with the CKD factor has the CKD factor alone. The notes of `fu_<s>` say which case applies (`unbound_fraction_notes`).

## Liver cirrhosis

`CP_severity` between the classes interpolates every factor linearly, so scans and fits use intermediate severities:

```equations title="the shunt of the Child-Pugh classes, interpolated"
fcp_shunt = piecewise(0, CP_severity <= 0, 0 + 0.399489795918367 * (CP_severity - 0), CP_severity < 1, 0.399489795918367 + 0.298469387755102 * (CP_severity - 1), CP_severity < 2, 0.697959183673469 + 0.114795918367347 * (CP_severity - 2), CP_severity < 3, 0.812755102040816)
fcp_rest_flow = 1 + (fcp_cardiac_index - 1 - (fQli * (fcp_hepatic_artery_flow - 1) + fQki * (fcp_kidney_flow - 1))) / (fQhe + fQmu + fQsk + fQad + fQbo + fQgo)
```

The organs without a factor of their own (heart, muscle, skin, adipose, bone, gonads) close the cardiac output to the cardiac index of the class (`fcp_rest_flow`). The brain (autoregulation of the cerebral blood flow, a choice of the library) and the portal organs (the shunt represents the portal blood which does not reach the liver) keep their flow.

<!-- generated: impairment-child-pugh (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/impairment.toml` and the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| factor | A | B | C | source |
| --- | ---: | ---: | ---: | --- |
| `functional_liver` | 0.8982 | 0.7854 | 0.6655 | consensus of `liver_volume_per_bw`: `Kwon2001`, `Lin1998` |
| `shunt` | 0.3995 | 0.698 | 0.8128 | `ModelChoice` |
| `hepatic_artery_flow` | 1.3 | 2.3 | 3.4 | `Edginton2008` |
| `kidney_flow` | 0.88 | 0.65 | 0.48 | `Edginton2008` |
| `cardiac_index` | 1.099 | 1.21 | 1.36 | consensus of `cardiac_index`: `Dincer2005`, `Henriksen2001`, `Iwao1997` |
| `gfr` | 1 | 0.7 | 0.36 | `Edginton2008` |
| `albumin` | 0.8645 | 0.6798 | 0.4892 | consensus of `albumin`: `Chawla1998`, `Mesropyan2022`, `Viani1992`, `Yang2020`, `Froomes1999` |
| `agp` | 0.6 | 0.56 | 0.3 | `Edginton2008` |
| `hematocrit` | 0.9232 | 0.7411 | 0.7232 | consensus of `hematocrit`: `Mesropyan2022`, `Wong2001`, `Boker2000` |

A factor with a consensus is that of the studies of the database on the quantity, the ratio of the class to the controls of each study pooled over the studies, see [the database](../database.md); a value marked `*` is a choice of the library (for a consensus the interpolation of a class without a study).

??? note "Sources"

    - `Boker2000`: Böker KH, Franzen A, Wrobel M, Bahr MJ, Tietge U, Manns MP. Regulation of hepatic blood flow in patients with liver cirrhosis and after liver transplantation. Pathophysiology 6(4):251-260, 2000. [doi:10.1016/s0928-4680(99)00024-3](https://doi.org/10.1016/s0928-4680(99)00024-3)
    - `Chawla1998`: Chawla Y, Radha NS, Jang KD, Dilawari B. Portal Hemodynamics by Duplex Doppler Sonography in Different Grades of Cirrhosis. Dig Dis Sci 43(2):354-357, 1998. [doi:10.1023/a:1018814624307](https://doi.org/10.1023/a:1018814624307)
    - `Dincer2005`: Dincer D, Besisk F, Demirkol O, Demir K, Kaymakoglu S, Cakaloglu Y, Okten A. Relationships between hemodynamic alterations and Child-Pugh Score in patients with cirrhosis. Hepatogastroenterology 52(65):1521-1525, 2005. PMID 16201110.
    - `Edginton2008`: Edginton AN, Willmann S. Physiology-based simulations of a pathological condition: prediction of pharmacokinetics in patients with liver cirrhosis. Clin Pharmacokinet 47(11):743-752, 2008. [doi:10.2165/00003088-200847110-00005](https://doi.org/10.2165/00003088-200847110-00005)
    - `Froomes1999`: Froomes PR, Morgan DJ, Smallwood RA, Angus PW. Comparative effects of oxygen supplementation on theophylline and acetaminophen clearance in human cirrhosis. Gastroenterology 116(4):915-920, 1999. [doi:10.1016/s0016-5085(99)70075-2](https://doi.org/10.1016/s0016-5085(99)70075-2)
    - `Henriksen2001`: Henriksen JH, Bendtsen F, Møller S. Normal red cell cardiac output in the hyperkinetic syndrome of alcoholic cirrhosis. J Hepatol 34(5):782-783, 2001. [doi:10.1016/s0168-8278(01)00019-8](https://doi.org/10.1016/s0168-8278(01)00019-8)
    - `Iwao1997`: Iwao T, Oho K, Sakai T, Tayama C, Sato M, Nakano R, Yamawaki M, Toyonaga A, Tanikawa K. Splanchnic and extrasplanchnic arterial hemodynamics in patients with cirrhosis. J Hepatol 27(5):817-823, 1997. [doi:10.1016/s0168-8278(97)80318-2](https://doi.org/10.1016/s0168-8278(97)80318-2)
    - `Kwon2001`: Kwon HA, Matsui Y, Ha-Kawa SK, Kamiyama Y. Functional Hepatic Volume Measured by Technetium-99M-Galactosyl-Human Serum Albumin Liver Scintigraphy: Comparison Between Hepatocyte Volume and Liver Volume by Computed Tomography. Am J Gastroenterol 96(2):541-546, 2001. [doi:10.1111/j.1572-0241.2001.03556.x](https://doi.org/10.1111/j.1572-0241.2001.03556.x)
    - `Lin1998`: Lin XZ, Sun YN, Liu YH, Sheu BS, Cheng BN, Chen CY, Tsai HM, Shen CL. Liver volume in patients with or without chronic liver diseases. Hepatogastroenterology 45(22):1069-1074, 1998. PMID 9756008.
    - `Mesropyan2022`: Mesropyan N, Kupczyk PA, Dold L, Praktiknjo M, Chang J, Isaak A, Endler C, Kravchenko D, Bischoff LM, Sprinkart AM, Pieper CC, Kuetting D, Jansen C, Attenberger UI, Luetkens JA. Assessment of liver cirrhosis severity with extracellular volume fraction MRI. Sci Rep 12(1):9422, 2022. [doi:10.1038/s41598-022-13340-9](https://doi.org/10.1038/s41598-022-13340-9)
    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
    - `Viani1992`: Viani A, Rizzo G, Carrai M, Pacifici GM. Interindividual variability in the concentrations of albumin and alpha-1-acid glycoprotein in patients with renal or liver disease, newborns and healthy subjects: implications for binding of drugs. Int J Clin Pharmacol Ther Toxicol 30(4):128-133, 1992. PMID 1572757.
    - `Wong2001`: Wong F, Girgrah N, Graba J, Allidina Y, Liu P, Blendis L. The cardiac response to exercise in cirrhosis. Gut 49(2):268-275, 2001. [doi:10.1136/gut.49.2.268](https://doi.org/10.1136/gut.49.2.268)
    - `Yang2020`: Yang M, Zhang Y, Zhao W, Cheng W, Wang H, Guo S. Evaluation of liver function using liver parenchyma, spleen and portal vein signal intensities during the hepatobiliary phase in Gd-EOB-D TPA-enhanced MRI. BMC Med Imaging 20(1):119, 2020. [doi:10.1186/s12880-020-00519-7](https://doi.org/10.1186/s12880-020-00519-7)
<!-- end generated -->

The factors on the albumin, the hematocrit and the cardiac index are the consensus of the studies of the [parameter database](../database.md#physiology-and-the-meld-score) with a control group (`class_factors` of `pkpdsim.pbpk.impairment`; the entry of the table names the quantity, `consensus = "albumin"`): the ratio of each class to the controls of each study, pooled by the random-effects meta-analysis as for the proteins. They agree with the values of [Edginton & Willmann 2008](#ref-Edginton2008) (Table I) for the albumin (0.81, 0.68, 0.50) and the cardiac index (1.11, 1.27, 1.36); the hematocrit falls more than theirs (0.91, 0.86, 0.81). The liver tissue `f_functional`, which scales the tissue volume of the liver and the capacity of its processes, is the consensus of the liver volume per body weight ([Kwon 2001](#ref-Kwon2001), [Lin 1998](#ref-Lin1998)): the factors of the enzymes and transporters are abundances per amount of cirrhotic liver tissue, which holds fewer hepatocytes per gram, so the capacity of the liver is the abundance times the amount of tissue, and the functional liver mass of [Edginton & Willmann 2008](#ref-Edginton2008) (0.69, 0.55, 0.28, a count of the hepatocytes) would count their loss a second time. The other physiological factors are those of Table I of [Edginton & Willmann 2008](#ref-Edginton2008), each from the primary study its footnote names (the renal flow index of [Dincer 2005](#ref-Dincer2005), the GFR of inulin clearance, the AGP), and the notes of each factor name it; the factors of the hepatic arterial flow are the ratios of the arterial perfusion per volume of liver of [Annet 2003](#ref-Annet2003), larger than the absolute arterial flow of [Sudhamshu 2006](#ref-Sudhamshu2006) and the perfusion of the database times its liver volume (the notes give both). The shunt and the portal blood flow describe the same portal blood which does not reach the sinusoids; applied together they would count it twice. The library takes the shunt: the shunt fractions are the degrees of cirrhosis of the Child-Pugh classes of `pkpdsim.impairment.CIRRHOSIS` (the ICG models of the König lab, a choice of the library), the splanchnic flow of the portal organs is not reduced, and the fraction `f_shunt` of the liver inflow bypasses the liver, so the portal blood which perfuses the liver is 0.60, 0.30 and 0.19 of healthy ([Edginton 2008](#ref-Edginton2008): 0.40, 0.36, 0.04). The shunt acts on the hepatic arterial inflow as well.

The enzymes and transporters of the liver are the consensus of the studies of the [parameter database](../database.md): the ratio of the abundance per amount of liver tissue of each Child-Pugh class to the controls of each study on the log scale, pooled over the studies by a random-effects meta-analysis ([El-Khateeb 2021](#ref-ElKhateeb2021), [Drozdzik 2020](#ref-Drozdzik2020), [Drozdzik 2021](#ref-Drozdzik2021) and [Prasad 2018](#ref-Prasad2018), whose end-stage livers are assumed to be class C), the mRNA expression only for a protein without an abundance. The notes of every factor give the consensus of each class with its 95 % confidence interval and the ratio of every study. A class without a study is the interpolation between its neighbours (healthy 1) or the last class with a study, a choice of the library. [Johnson 2010](#ref-Johnson2010) (the abundance factors of the Simcyp cirrhosis populations) was not accessible and is not used.

<!-- generated: impairment-child-pugh-proteins (scripts/pbpk_tables.py), do not edit -->
*Generated from the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| protein | A | B | C | studies |
| --- | ---: | ---: | ---: | --- |
| `ABCB1` | 1.317 | 1.174 | 1.027 | `ElKhateeb2021`, `Drozdzik2020`, `Wang2016` |
| `ABCC2` | 0.2696 | 0.2959 | 0.265 | `Drozdzik2020`, `ElKhateeb2021`, `Wang2016` |
| `ABCC3` | 0.474 | 0.4243 | 0.7358 | `Drozdzik2020`, `ElKhateeb2021`, `Wang2016` |
| `ABCC4` | 0.3636 | 0.6364 | 0.2727 | `ElKhateeb2021` |
| `ABCG2` | 0.75 | 1.25 | 0.973 | `ElKhateeb2021`, `Wang2016` |
| `CYP1A2` | 0.481 | 0.536 | 0.1869 | `Drozdzik2021`, `ElKhateeb2021`, `Prasad2018` |
| `CYP2C19` | 0.793 | 0.1162 | 0.2467 | `Drozdzik2021`, `ElKhateeb2021` |
| `CYP2C8` | 0.6561 | 0.3706 | 0.3428 | `Drozdzik2021`, `ElKhateeb2021`, `Prasad2018` |
| `CYP2C9` | 0.9053 | 0.7152 | 0.415 | `Drozdzik2021`, `ElKhateeb2021`, `Prasad2018` |
| `CYP2D6` | 0.824 | 0.7194 | 0.5954 | `Drozdzik2021`, `ElKhateeb2021`, `Prasad2018` |
| `CYP3A4` | 0.5981 | 0.3449 | 0.195 | `Drozdzik2021`, `ElKhateeb2021`, `Prasad2018` |
| `SLC22A1` | 0.4822 | 0.3041 | 0.2935 | `Drozdzik2020`, `ElKhateeb2021`, `Wang2016` |
| `SLC47A1` | 1.037* | 1.075* | 1.112 | `Wang2016` |
| `SLCO1B1` | 0.5187 | 0.4601 | 0.4826 | `Drozdzik2020`, `ElKhateeb2021`, `Wang2016` |
| `SLCO1B3` | 0.9863 | 0.3484 | 0.2852 | `Drozdzik2020`, `ElKhateeb2021`, `Wang2016` |
| `UGT1A1` | 1.49 | 1.204 | 0.9976 | `Drozdzik2021`, `ElKhateeb2021` |
| `UGT1A9` | 1.356 | 0.4759 | 0.2799 | `ElKhateeb2021` |
| `UGT2B7` | 0.7633 | 0.4251 | 0.2091 | `Drozdzik2021`, `ElKhateeb2021`, `Prasad2018` |

The consensus of the studies of the database on the abundance per amount of liver tissue over the controls, see [the database](../database.md); a value marked `*` is the interpolation of a class without a study, a choice of the library.

??? note "Sources"

    - `Drozdzik2020`: Drozdzik M, Szelag-Pieniek S, Post M, Zeair S, Wrzesinski M, Kurzawski M, Prieto J, Oswald S. Protein abundance of hepatic drug transporters in patients with different forms of liver damage. Clin Pharmacol Ther 107(5):1138-1148, 2020. [doi:10.1002/cpt.1717](https://doi.org/10.1002/cpt.1717)
    - `Drozdzik2021`: Drozdzik M, Lapczuk-Romanska J, Wenzel C, Szelag-Pieniek S, Post M, Skalski Ł, Kurzawski M, Oswald S. Gene expression and protein abundance of hepatic drug metabolizing enzymes in liver pathology. Pharmaceutics 13(9):1334, 2021. [doi:10.3390/pharmaceutics13091334](https://doi.org/10.3390/pharmaceutics13091334)
    - `ElKhateeb2021`: El-Khateeb E, Achour B, Al-Majdoub ZM, Barber J, Rostami-Hodjegan A. Non-uniformity of changes in drug-metabolizing enzymes and transporters in liver cirrhosis: implications for drug dosage adjustment. Mol Pharm 18(9):3563-3577, 2021. [doi:10.1021/acs.molpharmaceut.1c00462](https://doi.org/10.1021/acs.molpharmaceut.1c00462)
    - `Prasad2018`: Prasad B, Bhatt DK, Johnson K, Chapa R, Chu X, Salphati L, Xiao G, Lee C, Hop CE, Mathias A, Lai Y, Liao M, Humphreys WG, Kumer SC, Unadkat JD. Abundance of phase 1 and 2 drug-metabolizing enzymes in alcoholic and hepatitis C cirrhotic livers: a quantitative targeted proteomics study. Drug Metab Dispos 46(7):943-952, 2018. [doi:10.1124/dmd.118.080523](https://doi.org/10.1124/dmd.118.080523)
    - `Wang2016`: Wang L, Collins C, Kelly EJ, Chu X, Ray AS, Salphati L, Xiao G, Lee C, Lai Y, Liao M, Mathias A, Evers R, Humphreys W, Hop CE, Kumer SC, Unadkat JD. Transporter Expression in Liver Tissue from Subjects with Alcoholic or Hepatitis C Cirrhosis Quantified by Targeted Quantitative Proteomics. Drug Metab Dispos 44(11):1752-1758, 2016. [doi:10.1124/dmd.116.071050](https://doi.org/10.1124/dmd.116.071050)
<!-- end generated -->

### Comparison with the Simcyp populations

The parameters of the cirrhosis populations of the Simcyp Simulator ([Johnson 2010](#ref-Johnson2010), Table III) come from a meta-analysis of the literature as well, with the activity of the enzymes from in vivo probes and in vitro studies, and the functional hepatocyte volume as the amount of liver. The table compares them with the factors of the library; the enzymes are compared as the capacity of the whole liver, the abundance per amount of tissue times the liver volume, since each source counts the loss of hepatocytes in one of its two terms.

<!-- generated: impairment-child-pugh-comparison (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/comparisons.toml` and the factors of the library by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| factor | kind | library A / B / C | Johnson 2010 A / B / C |
| --- | --- | --- | --- |
| `functional_liver` | physiology | 0.898 / 0.785 / 0.665 | 0.81 / 0.65 / 0.53 |
| `albumin` | physiology | 0.865 / 0.68 / 0.489 | 0.919 / 0.758 / 0.588 |
| `agp` | physiology | 0.6 / 0.56 / 0.3 | 0.712 / 0.65 / 0.575 |
| `hematocrit` | physiology | 0.923 / 0.741 / 0.723 | 0.895 / 0.804 / 0.78 |
| `cardiac_index` | physiology | 1.099 / 1.21 / 1.36 | 1.16 / 1.317 / 1.408 |
| `hepatic_artery_flow` | physiology | 1.3 / 2.3 / 3.4 | 1.407 / 1.623 / 1.915 |
| `gfr` | physiology | 1 / 0.7 / 0.36 | 0.698 / 0.583 / 0.554 |
| `CYP1A2` | capacity | 0.432 / 0.421 / 0.124 | 0.512 / 0.17 / 0.062 |
| `CYP2A6` | capacity | 0.396 / 0.354 / 0.209 | 0.717 / 0.4 / 0.17 |
| `CYP2B6` | capacity | 0.666 / 0.508 / 0.306 | 0.81 / 0.585 / 0.424 |
| `CYP2C8` | capacity | 0.589 / 0.291 / 0.228 | 0.56 / 0.339 / 0.174 |
| `CYP2C9` | capacity | 0.813 / 0.562 / 0.276 | 0.559 / 0.338 / 0.175 |
| `CYP2C19` | capacity | 0.712 / 0.091 / 0.164 | 0.26 / 0.167 / 0.064 |
| `CYP2D6` | capacity | 0.74 / 0.565 / 0.396 | 0.618 / 0.211 / 0.056 |
| `CYP2E1` | capacity | 0.461 / 0.329 / 0.18 | 0.599 / 0.312 / 0.058 |
| `CYP3A4` | capacity | 0.537 / 0.271 / 0.13 | 0.478 / 0.252 / 0.132 |

The physiological factors relative to healthy; for an enzyme the capacity of the liver relative to healthy, the abundance per amount of tissue times the liver volume in the library (`functional_liver`), the abundance per mg of microsomal protein times the functional hepatocyte volume in Johnson 2010 (Table III, the parameters of the Simcyp cirrhosis populations).

??? note "Sources"

    - `Johnson2010`: Johnson TN, Boussery K, Rowland-Yeo K, Tucker GT, Rostami-Hodjegan A. A semi-mechanistic model to predict the effects of liver cirrhosis on drug clearance. Clin Pharmacokinet 49(3):189-206, 2010. [doi:10.2165/11318160-000000000-00000](https://doi.org/10.2165/11318160-000000000-00000)
<!-- end generated -->

### The MELD score

A MELD score gives the severity of cirrhosis: `cp_severity_of_meld(meld)` of `pkpdsim.pbpk.impairment` is piecewise linear in the MELD through anchors, the MELD of healthy subjects at severity 0 and the MELD of each Child-Pugh class at its severity, and `meld_of_cp_severity(severity)` its inverse (`meld_anchors()` the anchors). Healthy is 6.43, the smallest score of the formula `3.78 ln(bilirubin) + 11.2 ln(INR) + 9.57 ln(creatinine) + 6.43`, whose laboratory values are bounded below by 1 (`MELD_HEALTHY`); the MELD of a class is the consensus of the cohorts of the database on the MELD of its patients (`integration.level`). Between the anchors the factors and the MELD are both linear in the severity, a choice of the library; a MELD below 6.43 is healthy and one above the MELD of class C is class C, since the library has no factors beyond it. A model sets `CP_severity` of a subject or a group of reported MELD from it:

```python
from pkpdsim.pbpk.impairment import cp_severity_of_meld, meld_of_cp_severity

print(f"MELD 15: CP_severity {cp_severity_of_meld(15):.3f}")
print(f"class B: MELD {meld_of_cp_severity(2):.3g}")
```

<!-- generated: impairment-meld (scripts/pbpk_tables.py), do not edit -->
*Generated from the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| `CP_severity` | class | MELD | 95 % CI | studies (MELD) |
| ---: | --- | ---: | --- | --- |
| 0 | healthy | 6.43 |  |  |
| 1 | A | 7.909 | 7.554-8.265 | `Hong2025` 8, `Kim2020` 8, `Mesropyan2022` 9.3, `Xu2026` 8, `Yang2020` 7.5 |
| 2 | B | 12.78 | 10.8-14.76 | `Hong2025` 14, `Kim2020` 13, `Mesropyan2022` 11.9, `Xu2026` 15, `Yang2020` 10.04 |
| 3 | C | 20.5 | 17.01-23.99 | `Hong2025` 22.1, `Kim2020` 22, `Mesropyan2022` 17.9, `Xu2026` 24, `Yang2020` 16.33 |

??? note "Sources"

    - `Hong2025`: Hong SJ, Kim JE, Cho JM, Choi HC, Park MJ, Choi HY, Shin HS, Won JH, Yang W, Kim HO. Dual Energy CT-Derived Liver Extracellular Volume Fraction for Assessing Liver Functional Reserve in Patients with Liver Cirrhosis. Medicina 61(9):1561, 2025. [doi:10.3390/medicina61091561](https://doi.org/10.3390/medicina61091561)
    - `Kim2020`: Kim TH, Yun SG, Choi J, Goh HG, Lee HA, Yim SY, Choi SJ, Lee YS, Yoon EL, Jung YK, Seo YS, Kim JH, Yim HJ, Yeon JE, Byun KS, Um SH. Differential Impact of Serum 25-Hydroxyvitamin D3 Levels on the Prognosis of Patients with Liver Cirrhosis According to MELD and Child-Pugh Scores. J Korean Med Sci 35(19):e129, 2020. [doi:10.3346/jkms.2020.35.e129](https://doi.org/10.3346/jkms.2020.35.e129)
    - `Mesropyan2022`: Mesropyan N, Kupczyk PA, Dold L, Praktiknjo M, Chang J, Isaak A, Endler C, Kravchenko D, Bischoff LM, Sprinkart AM, Pieper CC, Kuetting D, Jansen C, Attenberger UI, Luetkens JA. Assessment of liver cirrhosis severity with extracellular volume fraction MRI. Sci Rep 12(1):9422, 2022. [doi:10.1038/s41598-022-13340-9](https://doi.org/10.1038/s41598-022-13340-9)
    - `Xu2026`: Xu L, Jiang B, Ma Y, Pan J, Liang M, Du Q, Wu H, Cheng Z, Lu D. Comparison of Child-Pugh and MELD scores in predicting survival of hepatocellular carcinoma patients with splenomegaly undergoing TACE: a real-world study. Front Med 13:1804646, 2026. [doi:10.3389/fmed.2026.1804646](https://doi.org/10.3389/fmed.2026.1804646)
    - `Yang2020`: Yang M, Zhang Y, Zhao W, Cheng W, Wang H, Guo S. Evaluation of liver function using liver parenchyma, spleen and portal vein signal intensities during the hepatobiliary phase in Gd-EOB-D TPA-enhanced MRI. BMC Med Imaging 20(1):119, 2020. [doi:10.1186/s12880-020-00519-7](https://doi.org/10.1186/s12880-020-00519-7)
<!-- end generated -->

## Chronic kidney disease

`f_renal_function` is the factor of the GFR, as in the hand-written models (`KI__f_renal_function` there, `f_renal_function` of the top model here): `GFR = GFR_base * f_renal_function`, with `GFR_base` the GFR at renal function 1, which carries the body surface area and the factor `fimp_GFR` of cirrhosis and heart failure (see [Physiology](physiology.md#scaling-to-the-individual)). Its GFR per 1.73 m2 of the reference individual of the build (`eGFR`) gives the stage: normal at least 90, mild 60-89, moderate 30-59, severe 15-29, end stage below 15 ml/min/1.73 m2. The factors are piecewise linear in the eGFR through anchor points at the stage midpoints: 1 from 90, the value of mild CKD at 75, of moderate CKD at 45, of severe CKD at 22.5 and of end stage CKD at 7.5 ml/min/1.73 m2 and below, linear in between (from 1 at 90 to the mild value at 75), as the factors of the Child-Pugh classes are between the classes. A scan or fit of `f_renal_function` therefore sees no step at a stage limit. The anchors are a choice of the library (`ModelChoice`: continuous covariate for scans and fits; stage midpoints), the values of the stages and their sources are unchanged.

The renal function scales the tubular capacity of the kidney as well, continuously and not by stage: every metabolism and transport of the kidney (a tubular secretion, a reabsorption by a transporter, a renal metabolism) runs with `1 - f_tubular + f_tubular*f_renal_function` of its healthy capacity (see [Chemistry](chemistry.md#kinetics)). The share `f_tubular` (a constant of the kidney, `KI__f_tubular`) is 1, the intact nephron hypothesis ([Bricker 1960](#ref-Bricker1960)): the nephrons which remain in a diseased kidney work normally, so the tubular capacity falls with their number as the GFR does, as the secretion of a fraction of the GFR of the hand-written models did. The value is a choice of the library, the publication is not verified; the filtration `fu*GFR` follows the GFR itself.

```equations title="the stage, the binding capacity and the CKD profile of albumin (male build)"
eGFR = 0.114467156025102 l_per_min * f_renal_function
fimp_albumin_ckd = fckd_albumin * fckd_albumin_binding
fckd_plasma_protein = piecewise(1, eGFR >= 0.09, 1, eGFR >= 0.075, 1 + 0.45 * (0.075 - eGFR) / 0.03, eGFR >= 0.045, 1.45 + 0.75 * (0.045 - eGFR) / 0.0225, eGFR >= 0.0225, 2.2 + 0.3 * (0.0225 - eGFR) / 0.015, eGFR >= 0.0075, 2.5)
```

<!-- generated: impairment-ckd (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/impairment.toml` and the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

Stage limits of the GFR per 1.73 m2 (ml/min/1.73 m2): 90, 60, 30, 15 (`Yoshida2016`).

Anchors of the interpolation, the eGFR (ml/min/1.73 m2) of the values of mild, moderate, severe and end stage: 75, 45, 22.5, 7.5 (`ModelChoice`).

| factor | mild | moderate | severe | end stage | source |
| --- | ---: | ---: | ---: | ---: | --- |
| `plasma_protein` | 1* | 1.45 | 2.2 | 2.5 | `Alasmari2023` |
| `hematocrit` | 0.9349 | 0.8857 | 0.7909 | 0.7106 | consensus of `hemoglobin`: `Panjeta2017`, `Lin2017`, `Zakiyanov2013` |
| `kidney_volume` | 1 | 1 | 1 | 1 | `ModelChoice` |
| `kidney_flow` | 0.7712 | 0.5567 | 0.4183 | 0.3344 | consensus of `renal_cortical_perfusion`: `Brown2020`, `ZhangX2024`, `LiLP2017` |
| `gastric_emptying_time` | 1 | 1 | 1.6 | 1.6 | `Malik2020` |
| `si_transit_time` | 1 | 1 | 1.4 | 1.4 | `Malik2020` |

A factor with a consensus is that of the studies of the database on the quantity, the ratio of the class to the controls of each study pooled over the studies, see [the database](../database.md); a value marked `*` is a choice of the library (for a consensus the interpolation of a class without a study).

Share of the tubular capacity which follows the renal function (`f_tubular`, not per stage): 1 (`ModelChoice`).

??? note "Sources"

    - `Alasmari2023`: Alasmari F, Alasmari MS, Muwainea HM, Alomar HA, Alasmari AF, Alsanea S, Alshamsan A, Rasool MF, Alqahtani F. Physiologically-based pharmacokinetic modeling for single and multiple dosing regimens of ceftriaxone in healthy and chronic kidney disease populations: a tool for model-informed precision dosing. Front Pharmacol 14:1200828, 2023. [doi:10.3389/fphar.2023.1200828](https://doi.org/10.3389/fphar.2023.1200828)
    - `Brown2020`: Brown RS, Sun MRM, Stillman IE, Russell TL, Rosas SE, Wei JL. The utility of magnetic resonance imaging for noninvasive evaluation of diabetic nephropathy. Nephrol Dial Transplant 35(6):970-978, 2020. [doi:10.1093/ndt/gfz066](https://doi.org/10.1093/ndt/gfz066)
    - `LiLP2017`: Li LP, Tan H, Thacker JM, Li W, Zhou Y, Kohn O, Sprague SM, Prasad PV. Evaluation of Renal Blood Flow in Chronic Kidney Disease Using Arterial Spin Labeling Perfusion Magnetic Resonance Imaging. Kidney Int Rep 2(1):36-43, 2017. [doi:10.1016/j.ekir.2016.09.003](https://doi.org/10.1016/j.ekir.2016.09.003)
    - `Lin2017`: Lin HY, Lee YL, Lin KD, Chiu YW, Shin SJ, Hwang SJ, Chen HC, Hung CC. Association of Renal Elasticity and Renal Function Progression in Patients with Chronic Kidney Disease Evaluated by Real-Time Ultrasound Elastography. Sci Rep 7:43303, 2017. [doi:10.1038/srep43303](https://doi.org/10.1038/srep43303)
    - `Malik2020`: Malik PRV, Yeung CHT, Ismaeil S, Advani U, Djie S, Edginton AN. A physiological approach to pharmacokinetics in chronic kidney disease. J Clin Pharmacol 60(S1):S52-S62, 2020. [doi:10.1002/jcph.1713](https://doi.org/10.1002/jcph.1713)
    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
    - `Panjeta2017`: Panjeta M, Tahirović I, Sofić E, Ćorić J, Dervišević A. Interpretation of Erythropoietin and Haemoglobin Levels in Patients with Various Stages of Chronic Kidney Disease. J Med Biochem 36(2):145-152, 2017. [doi:10.1515/jomb-2017-0014](https://doi.org/10.1515/jomb-2017-0014)
    - `Yoshida2016`: Yoshida K, Sun B, Zhang L, Zhao P, Abernethy DR, Nolin TD, Rostami-Hodjegan A, Zineh I, Huang SM. Systematic and quantitative assessment of the effect of chronic kidney disease on CYP2D6 and CYP3A4/5. Clin Pharmacol Ther 100(1):75-87, 2016. [doi:10.1002/cpt.337](https://doi.org/10.1002/cpt.337)
    - `Zakiyanov2013`: Zakiyanov O, Kriha V, Vachek J, Zima T, Tesar V, Kalousova M. Placental growth factor, pregnancy-associated plasma protein-A, soluble receptor for advanced glycation end products, extracellular newly identified receptor for receptor for advanced glycation end products binding protein and high mobility group box 1 levels in patients with acute kidney injury: a cross sectional study. BMC Nephrol 14:245, 2013. [doi:10.1186/1471-2369-14-245](https://doi.org/10.1186/1471-2369-14-245)
    - `ZhangX2024`: Zhang X, Ye C, Lu F, Yang J, Xu Y, Wang C. Evaluation of renal oxygenation and perfusion in patients with chronic kidney disease: a preliminary prospective study based on functional magnetic resonance. Ren Fail 46(2):2428337, 2024. [doi:10.1080/0886022X.2024.2428337](https://doi.org/10.1080/0886022X.2024.2428337)
<!-- end generated -->

The factor `plasma_protein` of the binding capacity acts on the binding to an unknown protein only (the binding to albumin and AGP has its profile, below). It is a factor on the binding capacity in the form of fu above, not a factor on fu: binding capacity is the mechanistic quantity and keeps fu below 1 (at `fu_ref` 0.5 the end stage gives fu 0.71). Its source is Table 3 of [Alasmari 2023](#ref-Alasmari2023), the "correction factor for the effect of the CKD on the content of plasma proteins" of the PK-Sim plasma protein scale factor, taken from [Malik 2020](#ref-Malik2020) (doi:10.1002/jcph.1713, paywalled, not read); its direction (the capacity falls by it) is a choice of the library. [Malik 2020](#ref-Malik2020) parameterize CKD for drugs secreted by the organic cation system, which do not separate the proteins.

The hepatic proteins (switch `CKD_hepatic`, off by default) are the mean ratios of the unbound clearance of model drugs ([Yoshida 2016](#ref-Yoshida2016), [Tan 2018](#ref-Tan2018)), taken as the activity of the pathway:

<!-- generated: impairment-ckd-proteins (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/impairment.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

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

### Binding in chronic kidney disease

The binding to albumin and to alpha-1-acid glycoprotein has a CKD profile per protein, a factor on the binding of each stage, piecewise linear in the eGFR as the other factors; the factor of a protein (`fimp_albumin_ckd`, `fimp_agp_ckd`) is the product of its rows, and it multiplies the factor of cirrhosis in the binding factor of fu. Its parameters are created only when a substance binds to the protein with the CKD factor, so a model without such a substance keeps its SBML.

- **Albumin, concentration** (`albumin`): the consensus of the serum albumin of non-dialysis CKD against healthy controls ([Huang 2017](#ref-Huang2017), [Lin 2017](#ref-Lin2017), [Zakiyanov 2013](#ref-Zakiyanov2013)), 0.93, 0.91 and 0.86 of the controls in the moderate, severe and end stage; the mild stage has no study and is interpolated. The meta-analysis of [Schmulenson 2018](#ref-Schmulenson2018) (PAGE poster, not peer reviewed) gives 0.982, 0.953, 0.910 and 0.840 for mild CKD to end stage (with dialysis).
- **Albumin, binding function** (`albumin_binding`, [Klammt 2012](#ref-Klammt2012), Nephrol Dial Transplant, doi:10.1093/ndt/gfr616, abstract): the albumin-binding capacity of site II (ABiC: the plasma is diluted to an albumin concentration of 150 umol/l and incubated with the site II marker dansylsarcosine at 150 umol/l, and ABiC is the unbound marker of the reference, a pooled plasma, over the unbound marker of the sample, times 100, so a higher ABiC is more binding; the method as Klinkmann 2022, doi:10.3390/medicina58121780, describe it) of 104 CKD patients is 118, 111, 99 and 79 % in the stages 1/2, 3, 4 and 5, with the protein-bound uraemic toxin indoxyl sulphate rising from 3.9 to 56.1 umol/l. Relative to the stages 1/2 it is 0.941, 0.839 and 0.669: a capacity of the albumin molecule, measured at a fixed albumin concentration and so separate from the concentration; the two multiply. The ratio of the unbound marker is taken as the factor of the binding per molecule, which is the factor on fu of a strongly bound drug (approximate for a weakly bound one). The stages 1/2 are taken as healthy and site II for every albumin-bound drug, choices of the library. The other reference of the study, its pooled plasma (100 %), would give 1.18 for the stages 1/2 and 1.11, 0.99 and 0.79 for moderate CKD to end stage (a binding above healthy in early CKD), and at the end stage fu of a strongly bound albumin substance x1.40 instead of x1.65; the stages 1/2 are chosen as the patients of the study with the least loss of renal function, the pooled plasma is of unknown donors.
- **AGP, concentration** (`agp`): the consensus of the alpha-1-acid glycoprotein against healthy controls, 1.45 at the end stage (95 % CI 1.04-2.03; [Docci 1985](#ref-Docci1985) 1.23 by the serum creatinine, [Zakiyanov 2013](#ref-Zakiyanov2013) 1.72), the stages before it without a study and interpolated between healthy and the end stage (1.11, 1.23, 1.34). On dialysis the consensus is 1.69. The meta-analysis of [Schmulenson 2018](#ref-Schmulenson2018) (PAGE poster, not peer reviewed) gives 1.137, 1.328 and 2.062 for moderate, severe CKD and end stage renal disease, its end stage with dialysis patients; the abstracts of Vasson 1991 (doi:10.1177/039139889101400207) and Vasson 1993 (doi:10.1159/000187492) give 1.61 and 1.45 for non-dialysed uraemic patients.

The albumin concentration and the AGP of the end stage are of non-dialysis CKD 5, the library's end stage; the albumin binding function ([Klammt 2012](#ref-Klammt2012), stage 5) includes dialysis patients, whose binding differs between and on the days of dialysis. Mild CKD is measured by none of the sources: its value lies on the line from 1 at the eGFR 90 to the next stage with a study (albumin concentration, AGP), and the binding function of the stages 1/2 is 1, choices of the library. At the end stage the profile raises fu of a strongly bound albumin substance by 1/(0.855*0.669) = 1.75 (the binding capacity 2.5) and lowers fu of a strongly bound AGP substance to 0.69 of healthy (2.5 before). The measured rise of fu of highly bound acidic drugs in severe CKD reaches 2.71 (bumetanide, Tan 2022, doi:10.1002/cpt.2642, full text), while the concentration of albumin alone predicts at most 1.20 (Tan 2022); Vanholder 1988 (doi:10.1038/ki.1988.99, abstract) measured a falling binding of phenytoin, theophylline and methotrexate and a rising binding of the basic drugs propranolol, cimetidine and clonidine in dialysed patients.

<!-- generated: impairment-ckd-binding (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/impairment.toml` and the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| factor | mild | moderate | severe | end stage | source |
| --- | ---: | ---: | ---: | ---: | --- |
| `albumin` | 0.9654* | 0.9308 | 0.912 | 0.8553 | consensus of `albumin`: `Huang2017`, `Lin2017`, `Zakiyanov2013` |
| `albumin_binding` | 1* | 0.941 | 0.839 | 0.669 | `Klammt2012` |
| `agp` | 1.113* | 1.227* | 1.341* | 1.454 | consensus of `agp`: `Docci1985`, `Zakiyanov2013` |

A factor with a consensus is that of the studies of the database on the quantity, the ratio of the class to the controls of each study pooled over the studies, see [the database](../database.md); a value marked `*` is a choice of the library (for a consensus the interpolation of a class without a study).

The rows: `albumin` the plasma albumin concentration (protein `albumin`); `albumin_binding` the binding function of an albumin molecule (protein `albumin`); `agp` the plasma alpha-1-acid glycoprotein concentration (protein `agp`). The factor of a protein is the product of its rows.

??? note "Sources"

    - `Docci1985`: Docci D, Bilancioni R, Pistocchi E, Mosconi G, Turci F, Salvi G, Baldrati L, Orsi C. Serum alpha-1-acid glycoprotein in chronic renal failure. Nephron 39(3):160-163, 1985. [doi:10.1159/000183364](https://doi.org/10.1159/000183364)
    - `Huang2017`: Huang MJ, Wei RB, Wang Y, Su TY, Di P, Li QP, Yang X, Li P, Chen XM. Blood coagulation system in patients with chronic kidney disease: a prospective observational study. BMJ Open 7(5):e014294, 2017. [doi:10.1136/bmjopen-2016-014294](https://doi.org/10.1136/bmjopen-2016-014294)
    - `Klammt2012`: Klammt S, Wojak HJ, Mitzner A, Koball S, Rychly J, Reisinger EC, Mitzner S. Albumin-binding capacity (ABiC) is reduced in patients with chronic kidney disease along with an accumulation of protein-bound uraemic toxins. Nephrol Dial Transplant 27(6):2377-2383, 2012 (the abstract read). [doi:10.1093/ndt/gfr616](https://doi.org/10.1093/ndt/gfr616)
    - `Lin2017`: Lin HY, Lee YL, Lin KD, Chiu YW, Shin SJ, Hwang SJ, Chen HC, Hung CC. Association of Renal Elasticity and Renal Function Progression in Patients with Chronic Kidney Disease Evaluated by Real-Time Ultrasound Elastography. Sci Rep 7:43303, 2017. [doi:10.1038/srep43303](https://doi.org/10.1038/srep43303)
    - `Zakiyanov2013`: Zakiyanov O, Kriha V, Vachek J, Zima T, Tesar V, Kalousova M. Placental growth factor, pregnancy-associated plasma protein-A, soluble receptor for advanced glycation end products, extracellular newly identified receptor for receptor for advanced glycation end products binding protein and high mobility group box 1 levels in patients with acute kidney injury: a cross sectional study. BMC Nephrol 14:245, 2013. [doi:10.1186/1471-2369-14-245](https://doi.org/10.1186/1471-2369-14-245)
<!-- end generated -->

## Heart failure

The flows of the organs are multiplied by the factors of the NYHA class and the cardiac output is their sum ([Gu 2025](#ref-Gu2025), Table 1). The factors are the consensus of the database on the hepatic, renal and limb blood flow; the one study with the flows of single classes against controls is [Leithe 1984](#ref-Leithe1984) (16 normal subjects, 64 patients), whose subgroups are defined by the cardiac index (mild > 2.5, moderate 2.0-2.5, severe < 2.0 L/min/m2) and taken as NYHA II, III and IV, as [Rasool 2015](#ref-Rasool2015) and [Gu 2025](#ref-Gu2025) take them (the patients are NYHA II 12, III 28 and IV 24, the subgroups 18, 20 and 26). The means of its Figure 1 are digitized and give back the means of the whole group the text states. NYHA I has no study and is interpolated between healthy and NYHA II. The renal blood flow is not monotonic in the class (the severe subgroup above the moderate one, as measured), and the GFR follows it. The cardiac output of the model, the sum of the flows, falls to 0.89, 0.78, 0.65 and 0.62 of healthy (NYHA I to IV, the reference individual); the cardiac index of the subgroups of Leithe 1984 is 0.79, 0.63 and 0.50 of the normal subjects. The limb flow of Leithe 1984 (0.66, 0.57 and 0.45) replaces the factors of skin, adipose and muscle of [Gu 2025](#ref-Gu2025) (0.57, 0.44, 0.28), which are not its limb flow.

<!-- generated: impairment-heart-failure (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/impairment.toml` and the studies of `pkpdsim/database` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| factor | I | II | III | IV | source |
| --- | ---: | ---: | ---: | ---: | --- |
| `liver_flow` | 0.8728* | 0.7456 | 0.5267 | 0.4662 | consensus of `hepatic_blood_flow_per_bsa`: `Leithe1984` |
| `kidney_flow` | 0.8792* | 0.7584 | 0.5492 | 0.6298 | consensus of `renal_blood_flow_per_bsa`: `Leithe1984` |
| `peripheral_flow` | 0.8244* | 0.6488 | 0.5709 | 0.4654 | consensus of `forearm_blood_flow`: `Leithe1984` |

A factor with a consensus is that of the studies of the database on the quantity, the ratio of the class to the controls of each study pooled over the studies, see [the database](../database.md); a value marked `*` is a choice of the library (for a consensus the interpolation of a class without a study).

??? note "Sources"

    - `Leithe1984`: Leithe ME, Margorien RD, Hermiller JB, Unverferth DV, Leier CV. Relationship between central hemodynamics and regional blood flow in normal subjects and in patients with congestive heart failure. Circulation 69(1):57-64, 1984. [doi:10.1161/01.CIR.69.1.57](https://doi.org/10.1161/01.CIR.69.1.57)
<!-- end generated -->

## Classes in experiments

The experiments of a model built on the library map the classes of a study onto the severities with the presets of `pkpdsim.impairment`: `CHILD_PUGH` (colors `CIRRHOSIS_COLORS`) onto `CP_severity`, the `RENAL_FUNCTION` presets onto `f_renal_function`, `HEART_FAILURE` (colors `HEART_FAILURE_COLORS`) onto `HF_severity`. The base of their experiments is `PBPKSimulationExperiment` of [`pkpdsim.pbpk.experiment`](../api/pbpk_experiment.md), which has these maps: `renal_map` the classes of `RENAL_FUNCTION_GFR` (colors `RENAL_COLORS_END_STAGE`), `cirrhosis_map` `CHILD_PUGH` and `cardiac_map` `HEART_FAILURE`. The classes of renal impairment keep their GFR in ml/min: a preset is relative to a healthy GFR of 101 ml/min, the renal function of a class is its GFR over the healthy GFR `GFR_base` of the model (`ScaledRenalMap`, evaluated on the first access, so importing the experiments loads no model). This differs from the eGFR of the CKD factors above, which is relative to the GFR per 1.73 m2 of the reference individual. `healthy_gfr` reads `GFR_base` of an individual (body weight, height, `HF_severity`, `CP_severity`) from the flat model, a quantity in ml/min, and `renal_function(clearance, bw=..., height=...)` gives the renal function of a measured clearance in ml/min as a dimensionless quantity, so that the simulated GFR of the individual of the study is the measured one; both are quantities, which go into the changes of a simulation as they are. The renal function is limited by the validity of the kidney outflow (see [Standard organ](organ.md#validity-of-the-outflow)): the filtered plasma water, the GFR, and the net filtered plasma flow of every filtered substance, `fu*GFR*(1 - f_reabsorbed)` (`fu*GFR` with a threshold reabsorption), must stay below the plasma flow of the kidney, `Qki*(1 - HCT)` (`Qki*BP` with erythrocytes in equilibrium). `kidney_outflow(model_path, changes)` of [`pkpdsim.pbpk.kidney_validity`](../api/pbpk_kidney_validity.md) gives both flows of an individual, the water first (`KidneyOutflow`, in ml/min, with the formula of its net filtered plasma flow, `GFR` and the substance `None` for the water, read as `healthy_gfr` reads the GFR (`individual_values`), from a simulation of one minute on a model loaded once per process, cached per model file and individual; only the changes of the parameters the flows depend on make the individual, a dose does not; a model without a filtered substance is not simulated). `renal_function` raises a `ValueError` which names the GFR and the plasma flow for a renal function beyond it, `renal_function_limits(bw=..., height=..., ...)` gives the range of an individual, from 0 to the smaller of the renal functions at which the GFR and the net filtered plasma flow of a substance reach the plasma flow (the GFR, the tighter bound: 6.08 for the male and 6.20 for the female reference individual), and a model narrows it with `renal_function_range` (`None` by default, e.g. `(0.05, 2.0)` for a model whose urine is not defined without filtration). The initialization of every experiment checks every individual it simulates (`simulated_individuals`), each task on the model it runs, at the start of every simulation and every point of a scan: its preinit changes, of a presimulation as well, and the values of the dimensions without a time. It does not check a change during a simulation (`Simulation.changes`, also at its start, of an id the definition names `dynamic`), a dimension with a time (`at`), a scan over models and the draws of the fitted parameters of a run over them; the margin is large (about 6 times at renal function 1). It fails with a message which names the experiment, the simulation, the GFR, the plasma flow and the changes of the individual. All of this is of the build of the sex of the experiment (`sex`, see [Builds per sex](index.md#builds-per-sex)): the healthy GFR of a woman is `GFR_base` of the female build, so `healthy_gfr`, `renal_function`, `renal_function_limits` and `renal_map` of an experiment of women read the female build, each of them takes the sex of a simulation by its key with `simulation=` (`sex_of`, the sex of `simulation_sexes`: the group of women of a study of men, `renal_function(60.0, simulation="po_women")`, `renal_map_of(simulation="po_women")`) or a sex with `sex=`, and refuses a call with neither in an experiment of both sexes, and the kidney outflow of every task is checked on the build it runs. The names of `pkpdsim.pbpk.experiment` and `pkpdsim.pbpk.kidney_validity` are imported from `pkpdsim.pbpk` (`from pkpdsim.pbpk import REFERENCE_GFR, PBPKSimulationExperiment, kidney_outflow`), which exports them lazily. The 44 hand-written models keep their maps and ids unchanged: `f_cirrhosis` keeps its meaning there (shunt and tissue loss in one), and the library uses other ids.

## Equilibration of a condition

The state of a model at its start is the steady state of the healthy reference individual. A condition of a simulation, an impairment, a covariate or the parameters of a PD component with a homeostasis (the extracellular fluid and the blood pressure of a fluid model in renal impairment), moves a slow state away from it, which a study measures as the basal state of its subjects before the first dose. `PBPKSimulationExperiment.equilibrated(simulation, resets=...)` returns the simulation (or the scan) which starts from the drug-free steady state of its condition:

```python
# not executed
def simulations(self) -> dict[str, Simulation]:
    simulation = Simulation(
        time_unit="hr",
        end=24,
        preinit_changes={
            "f_renal_function": self.renal_function(29),
            "PODOSE_hctz": Q(25, "mg"),
        },
    )
    # the urine of the fluid PD component grows without a steady state
    resets = {"FLUID__Vurine": Q(1e-6, "l"), "FLUID__na_urine": Q(0, "mmole")}
    return {"hctz25": self.equilibrated(simulation, resets=resets)}
```

- The equilibration is the presimulation of sbmlsim (`SteadyState`): the model is initialized with the changes of the simulation before the initialization, the condition, and integrated until the rate of change of every state is within the tolerances, `|dx/dt| <= 1e-12 + 1e-10*|x|` per minute (`EQUILIBRATION_TOLERANCES`, tighter than the defaults of sbmlsim, whose relative tolerance of 1e-6 per minute stops a homeostasis of 7 days 5e-5 off its steady state); the simulation starts from that state, which is not in its output. The changes of a fit (the fitted parameters) and of a run over the draws of the fitted parameters apply during it, and every point of a scan equilibrates in its own condition.
- The doses before the initialization (the dose states and infusion rates of the library, `PODOSE_<s>`, `IVDOSE_<s>`, `SOLDOSE_<s>`, `Ri_<s>`) are given by a change at the start, after the equilibration; a dose at a later time stays. A dose before the initialization which a change at the start sets as well is refused, as is a dimension of a scan over a dose which the simulation does not give.
- The sinks have no steady state, their amounts grow: the urine and the feces of every substance (`Aurine_<s>`, `Afeces_<s>`), the readouts of the intestine which nothing drains, the amount absorbed and the oral dose applied (`Aabs_<s>`, `Aapp_<s>`, or `GU__Aabs_<s>`, `GU__Aapp_<s>` with the readout `Fa_<s>` alone; an endogenous substance which reaches the lumen at its baseline is absorbed again during the equilibration), and what a PD component accumulates (the urine of a fluid model, `resets`). Their rate falls within the relative tolerance as their amount grows (after about 1e10 minutes, within the longest equilibration of 1e14 minutes), and they are set back at the start: the library sinks and readouts to 0, the states of `resets` to their values (a species of `resets` named by its id is set as its quantity, `[sid]` for a concentration species, as the zeros below). A state which nothing drains (no rate reads it, no reaction consumes it, `accumulating_states`) and which grows during the equilibration is refused with its id unless it is a library sink or readout or in `resets` (`growing_states`, checked once per model file in a process, in the condition of the first equilibration); a model which reaches no steady state by `max_time` fails with sbmlsim's `SteadyStateError`.
- The species of the drugs start at exactly 0 (`drug_free_zeros`): the species of the substances without a baseline (the drugs and their metabolites, not an endogenous substance) which start at 0 and which no reaction feeds from a species with a start (a reactant which holds some, a product of a reversible reaction which does, a source which is no dose of the library; the direction is the `reversible` flag of the reaction, so a `CustomReaction` whose rate law can run backwards, an exchange `k*(A - B)`, is written reversible, `A <-> B`, else a species it feeds backwards is taken for 0) are 0 in exact arithmetic during a drug-free equilibration, and the presimulation leaves a residue of the rounding in them (1.5e-29 mM of hctz in its plasma), which a pre-dose data point with a prediction of zero would count in the fold metrics. They are set to 0 at the start, with the sinks, each as its quantity: `[sid]` for a species whose id is its concentration (`hasOnlySubstanceUnits` false), the id of an amount species, so that the start changes have a representation in PEtab, whose conditions set the concentration of such a species. A species which the condition or a dimension of a scan sets, what it feeds, and a species which a change of the simulation sets at its start keep their values; an endogenous substance (a baseline, creatinine) and what it feeds keep the equilibrated values.
- The base `<id>_0` of a delta observer (`ids.delta_base`, `PBPKDefinition.deltas`) is set to the value after the equilibration, so the delta is 0 at the start.

`equilibrated` has an overload for a simulation, one for a scan and one for the union `Simulation | Scan`, so a model which equilibrates a mapping of both calls it without a `match`. The equilibration lives in `pkpdsim.pbpk.equilibration` (`DOSE_PREFIXES`, `EQUILIBRATION_TOLERANCES`, `accumulating_states`, `drug_free_zeros`, `growing_states`; a model file is parsed once per process): `equilibrated_simulation(simulation, model_path, resets=...)`, also of `pkpdsim.pbpk` (exported lazily), reads the ids of the doses, the sinks, the readouts and the deltas from the model file, and `equilibrated` calls it with the file of the loaded model of the experiment. Outside an experiment the groups of a sensitivity analysis use it (`SensitivityGroup(prepare=Equilibration(resets=...))`, see [Scans and sensitivity](../scans.md)) instead of a pre-run of days from a negative start.

The presimulation was chosen over a hidden pre-run of a duration from a negative start (tested in `tests/pbpk/test_equilibration.py` on a homeostasis of a time constant of 7 days with and without a sink): the presimulation reaches the steady state to 1e-6, a pre-run of 28 days leaves the homeostasis 1.8 % off it, and its time before the dose would be in the output. On the hctz model with the fluid PD at the renal function 0.1 the equilibration takes a second and the extracellular fluid does not drift over 60 days after it.

## Limitations

- The physiological factors of the Child-Pugh classes but the liver volume, the albumin, the hematocrit and the cardiac index are those of [Edginton & Willmann 2008](#ref-Edginton2008), single studies each; the liver volume of class C rests on one study ([Lin 1998](#ref-Lin1998)).
- The anaemia of CKD is the consensus of the hemoglobin against healthy controls ([Panjeta 2017](#ref-Panjeta2017), [Lin 2017](#ref-Lin2017), [Zakiyanov 2013](#ref-Zakiyanov2013)), the hematocrit taken proportional to the hemoglobin; the renal blood flow of CKD is the consensus of the cortical perfusion by arterial spin labelling ([Zhang 2024](#ref-ZhangX2024), [Brown 2020](#ref-Brown2020), [Li 2017](#ref-LiLP2017)) . The kidney volume stays 1: the kidneys shrink in CKD ([Malik 2020](#ref-Malik2020), Figure 1), but the capacity of the processes of the kidney falls with the renal function already (the intact nephron hypothesis) and they are per litre of kidney, so a smaller kidney would count the loss of nephrons twice; the measured total renal flows of the database follow the ratio of the cortical perfusion.
- The factors of heart failure rest on one study, [Leithe 1984](#ref-Leithe1984), whose subgroups of the cardiac index are taken as NYHA classes; the other studies of the database give the flows of heart failure against controls only for groups of several classes (the notes of the factors compare them).
- The MELD score maps onto the Child-Pugh classes through the mean MELD of the patients of a class: a patient's MELD and class are correlated, not equivalent (the MELD of the classes overlaps), and a MELD beyond the mean of class C is class C.
- The gastric emptying and the small intestinal transit of CKD are those of [Malik 2020](#ref-Malik2020), the parameters of PK-Sim, one factor of the stages 4 and 5 (none for the stages 2 and 3); the kidney volume is 1 (see above).
- The binding factor of CKD of an unknown protein acts on the binding capacity, PK-Sim multiplies fu by it: the two differ for weakly bound drugs (at `fu_ref` 0.5 and the factor 1.45 of moderate CKD, fu 0.59 against 0.725). Its source ([Malik 2020](#ref-Malik2020)) was not read.
- The CKD profile of the binding is the mean of a population: the binding of a drug in CKD varies with its binding site and the uraemic toxins which displace it (Tan 2022: 1.0 to 2.7 times the fu of highly bound acidic drugs in severe CKD), so a drug whose fu is measured in renal impairment validates the profile with it. `Binding(ckd=False)` is only for a fraction unbound which is no binding measurement (a fitted renal factor) or whose binding is measured unchanged in renal impairment (apixaban); a drug with a measured change of fu in renal impairment keeps the CKD profile, and a measured change which differs from the profile is a finding of the model for its TODO, not a reason to switch the profile off. The AGP of the stages before the end stage is interpolated, no study with controls measures them; the binding function of albumin is that of site II.
- A Kp of a Kp method is a number of the build: a change of fu by an impairment does not recompute it. The interstitial to plasma coefficient of the full form does follow fu, so the interstitial share of a tissue rises with the fu of CKD while its Kp stays: where the share exceeds the Kp, the intracellular coefficient `K_cell_pls_<s>` is floored at `K_CELL_FLOOR` (1e-6), the cells hold almost nothing and the tissue holds its interstitial share, more than its Kp. The build warns (`KCellFloorWarning` of `pkpdsim.pbpk.distribution`, naming the substance, the organ, the stages with their fu and share and the sex of the build) when this happens at the fu of a CKD stage of its binding profile (`kp_floor_stages`, the fu of the stages from `ckd_unbound_fractions` of `pkpdsim.pbpk.impairment`); each build checks its own sex. It does not refuse: a model gives the Kp room for the fu of the end stage, or accepts the floor in an organ with its reason, `ckd_floor={organ: reason}` of the distribution (`Kp` or `KpMethod`). An accepted organ is not warned of, and the notes of its `K_cell_pls_<s>` carry the reason and the stages at which the floor holds in the build, so the archive records the decision; the acceptance is per substance and organ, every other organ is still checked, and a reason must not be empty (`validate`). The metabolite M1 of torasemide in the kidney reaches the floor at the end stage (Kp 0.08 against the share 0.0746 at the reference fu 0.14 and 0.0825 at the fu 0.221):

```python
# not executed
Kp(
    1.0,
    overrides={Organ.KIDNEY: 0.08},
    ckd_floor={
        Organ.KIDNEY: "the binding of M1 is not measured in renal impairment; "
        "the end stage keeps the interstitial share of the kidney"
    },
)
```
- The tubular capacity follows the renal function in proportion (`f_tubular` 1, the intact nephron hypothesis, [Bricker 1960](#ref-Bricker1960), not verified, a choice of the library); a transporter whose activity falls less or more than the GFR in CKD needs another share or an activity factor of the experiment.

## References

<!-- generated: references (scripts/pbpk_tables.py), do not edit -->
- <span id="ref-Alasmari2023"></span>Alasmari F, Alasmari MS, Muwainea HM, Alomar HA, Alasmari AF, Alsanea S, Alshamsan A, Rasool MF, Alqahtani F. Physiologically-based pharmacokinetic modeling for single and multiple dosing regimens of ceftriaxone in healthy and chronic kidney disease populations: a tool for model-informed precision dosing. Front Pharmacol 14:1200828, 2023. [doi:10.3389/fphar.2023.1200828](https://doi.org/10.3389/fphar.2023.1200828)
- <span id="ref-Annet2003"></span>Annet L, Materne R, Danse E, Jamart J, Horsmans Y, Van Beers BE. Hepatic Flow Parameters Measured with MR Imaging and Doppler US: Correlations with Degree of Cirrhosis and Portal Hypertension. Radiology 229(2):409-414, 2003. [doi:10.1148/radiol.2292021128](https://doi.org/10.1148/radiol.2292021128)
- <span id="ref-Bricker1960"></span>Bricker NS, Morrin PA, Kime SW Jr. The pathologic physiology of chronic Bright's disease. An exposition of the intact nephron hypothesis. Am J Med 28:77-98, 1960. [doi:10.1016/0002-9343(60)90225-4](https://doi.org/10.1016/0002-9343(60)90225-4)
- <span id="ref-Brown2020"></span>Brown RS, Sun MRM, Stillman IE, Russell TL, Rosas SE, Wei JL. The utility of magnetic resonance imaging for noninvasive evaluation of diabetic nephropathy. Nephrol Dial Transplant 35(6):970-978, 2020. [doi:10.1093/ndt/gfz066](https://doi.org/10.1093/ndt/gfz066)
- <span id="ref-Dincer2005"></span>Dincer D, Besisk F, Demirkol O, Demir K, Kaymakoglu S, Cakaloglu Y, Okten A. Relationships between hemodynamic alterations and Child-Pugh Score in patients with cirrhosis. Hepatogastroenterology 52(65):1521-1525, 2005. PMID 16201110.
- <span id="ref-Docci1985"></span>Docci D, Bilancioni R, Pistocchi E, Mosconi G, Turci F, Salvi G, Baldrati L, Orsi C. Serum alpha-1-acid glycoprotein in chronic renal failure. Nephron 39(3):160-163, 1985. [doi:10.1159/000183364](https://doi.org/10.1159/000183364)
- <span id="ref-Drozdzik2020"></span>Drozdzik M, Szelag-Pieniek S, Post M, Zeair S, Wrzesinski M, Kurzawski M, Prieto J, Oswald S. Protein abundance of hepatic drug transporters in patients with different forms of liver damage. Clin Pharmacol Ther 107(5):1138-1148, 2020. [doi:10.1002/cpt.1717](https://doi.org/10.1002/cpt.1717)
- <span id="ref-Drozdzik2021"></span>Drozdzik M, Lapczuk-Romanska J, Wenzel C, Szelag-Pieniek S, Post M, Skalski Ł, Kurzawski M, Oswald S. Gene expression and protein abundance of hepatic drug metabolizing enzymes in liver pathology. Pharmaceutics 13(9):1334, 2021. [doi:10.3390/pharmaceutics13091334](https://doi.org/10.3390/pharmaceutics13091334)
- <span id="ref-Edginton2008"></span>Edginton AN, Willmann S. Physiology-based simulations of a pathological condition: prediction of pharmacokinetics in patients with liver cirrhosis. Clin Pharmacokinet 47(11):743-752, 2008. [doi:10.2165/00003088-200847110-00005](https://doi.org/10.2165/00003088-200847110-00005)
- <span id="ref-ElKhateeb2021"></span>El-Khateeb E, Achour B, Al-Majdoub ZM, Barber J, Rostami-Hodjegan A. Non-uniformity of changes in drug-metabolizing enzymes and transporters in liver cirrhosis: implications for drug dosage adjustment. Mol Pharm 18(9):3563-3577, 2021. [doi:10.1021/acs.molpharmaceut.1c00462](https://doi.org/10.1021/acs.molpharmaceut.1c00462)
- <span id="ref-Gu2025"></span>Gu W, Shao Q, Jiang L. Predicting pharmacokinetics of drugs in patients with heart failure and optimizing their dosing strategies using a physiologically based pharmacokinetic model. Pharmaceutics 17(11):1394, 2025. [doi:10.3390/pharmaceutics17111394](https://doi.org/10.3390/pharmaceutics17111394)
- <span id="ref-Huang2017"></span>Huang MJ, Wei RB, Wang Y, Su TY, Di P, Li QP, Yang X, Li P, Chen XM. Blood coagulation system in patients with chronic kidney disease: a prospective observational study. BMJ Open 7(5):e014294, 2017. [doi:10.1136/bmjopen-2016-014294](https://doi.org/10.1136/bmjopen-2016-014294)
- <span id="ref-Johnson2010"></span>Johnson TN, Boussery K, Rowland-Yeo K, Tucker GT, Rostami-Hodjegan A. A semi-mechanistic model to predict the effects of liver cirrhosis on drug clearance. Clin Pharmacokinet 49(3):189-206, 2010. [doi:10.2165/11318160-000000000-00000](https://doi.org/10.2165/11318160-000000000-00000)
- <span id="ref-Klammt2012"></span>Klammt S, Wojak HJ, Mitzner A, Koball S, Rychly J, Reisinger EC, Mitzner S. Albumin-binding capacity (ABiC) is reduced in patients with chronic kidney disease along with an accumulation of protein-bound uraemic toxins. Nephrol Dial Transplant 27(6):2377-2383, 2012 (the abstract read). [doi:10.1093/ndt/gfr616](https://doi.org/10.1093/ndt/gfr616)
- <span id="ref-Kwon2001"></span>Kwon HA, Matsui Y, Ha-Kawa SK, Kamiyama Y. Functional Hepatic Volume Measured by Technetium-99M-Galactosyl-Human Serum Albumin Liver Scintigraphy: Comparison Between Hepatocyte Volume and Liver Volume by Computed Tomography. Am J Gastroenterol 96(2):541-546, 2001. [doi:10.1111/j.1572-0241.2001.03556.x](https://doi.org/10.1111/j.1572-0241.2001.03556.x)
- <span id="ref-Leithe1984"></span>Leithe ME, Margorien RD, Hermiller JB, Unverferth DV, Leier CV. Relationship between central hemodynamics and regional blood flow in normal subjects and in patients with congestive heart failure. Circulation 69(1):57-64, 1984. [doi:10.1161/01.CIR.69.1.57](https://doi.org/10.1161/01.CIR.69.1.57)
- <span id="ref-LiLP2017"></span>Li LP, Tan H, Thacker JM, Li W, Zhou Y, Kohn O, Sprague SM, Prasad PV. Evaluation of Renal Blood Flow in Chronic Kidney Disease Using Arterial Spin Labeling Perfusion Magnetic Resonance Imaging. Kidney Int Rep 2(1):36-43, 2017. [doi:10.1016/j.ekir.2016.09.003](https://doi.org/10.1016/j.ekir.2016.09.003)
- <span id="ref-Lin1998"></span>Lin XZ, Sun YN, Liu YH, Sheu BS, Cheng BN, Chen CY, Tsai HM, Shen CL. Liver volume in patients with or without chronic liver diseases. Hepatogastroenterology 45(22):1069-1074, 1998. PMID 9756008.
- <span id="ref-Lin2017"></span>Lin HY, Lee YL, Lin KD, Chiu YW, Shin SJ, Hwang SJ, Chen HC, Hung CC. Association of Renal Elasticity and Renal Function Progression in Patients with Chronic Kidney Disease Evaluated by Real-Time Ultrasound Elastography. Sci Rep 7:43303, 2017. [doi:10.1038/srep43303](https://doi.org/10.1038/srep43303)
- <span id="ref-Malik2020"></span>Malik PRV, Yeung CHT, Ismaeil S, Advani U, Djie S, Edginton AN. A physiological approach to pharmacokinetics in chronic kidney disease. J Clin Pharmacol 60(S1):S52-S62, 2020. [doi:10.1002/jcph.1713](https://doi.org/10.1002/jcph.1713)
- <span id="ref-Panjeta2017"></span>Panjeta M, Tahirović I, Sofić E, Ćorić J, Dervišević A. Interpretation of Erythropoietin and Haemoglobin Levels in Patients with Various Stages of Chronic Kidney Disease. J Med Biochem 36(2):145-152, 2017. [doi:10.1515/jomb-2017-0014](https://doi.org/10.1515/jomb-2017-0014)
- <span id="ref-Prasad2018"></span>Prasad B, Bhatt DK, Johnson K, Chapa R, Chu X, Salphati L, Xiao G, Lee C, Hop CE, Mathias A, Lai Y, Liao M, Humphreys WG, Kumer SC, Unadkat JD. Abundance of phase 1 and 2 drug-metabolizing enzymes in alcoholic and hepatitis C cirrhotic livers: a quantitative targeted proteomics study. Drug Metab Dispos 46(7):943-952, 2018. [doi:10.1124/dmd.118.080523](https://doi.org/10.1124/dmd.118.080523)
- <span id="ref-Rasool2015"></span>Rasool MF, Khalil F, Läer S. A physiologically based pharmacokinetic drug-disease model to predict carvedilol exposure in adult and paediatric heart failure patients by incorporating pathophysiological changes in hepatic and renal blood flows. Clin Pharmacokinet 54(9):943-962, 2015. [doi:10.1007/s40262-015-0253-7](https://doi.org/10.1007/s40262-015-0253-7)
- <span id="ref-Schmulenson2018"></span>Schmulenson E, Schlender JF, Frechen S, Jaehde U. A physiologically-based pharmacokinetic modeling approach to assess the impact of chronic kidney disease. PAGE 27 (2018), Abstr 8630, poster IV-61, https://www.page-meeting.org/?abstract=8630 (the poster read: https://www.page-meeting.org/pdf_assets/1160-Schmulenson_PAGE2018.pdf).
- <span id="ref-Sudhamshu2006"></span>Sudhamshu KC, Matsutani S, Maruyama H, Akiike T, Saisho H. Doppler study of hepatic vein in cirrhotic patients: Correlation with liver dysfunction and hepatic hemodynamics. World J Gastroenterol 12(36):5853-5858, 2006. [doi:10.3748/wjg.v12.i36.5853](https://doi.org/10.3748/wjg.v12.i36.5853)
- <span id="ref-Tan2018"></span>Tan ML, Yoshida K, Zhao P, Zhang L, Nolin TD, Piquette-Miller M, Galetin A, Huang SM. Effect of chronic kidney disease on nonrenal elimination pathways: a systematic assessment of CYP1A2, CYP2C8, CYP2C9, CYP2C19, and OATP. Clin Pharmacol Ther 103(5):854-867, 2018. [doi:10.1002/cpt.807](https://doi.org/10.1002/cpt.807)
- <span id="ref-Yoshida2016"></span>Yoshida K, Sun B, Zhang L, Zhao P, Abernethy DR, Nolin TD, Rostami-Hodjegan A, Zineh I, Huang SM. Systematic and quantitative assessment of the effect of chronic kidney disease on CYP2D6 and CYP3A4/5. Clin Pharmacol Ther 100(1):75-87, 2016. [doi:10.1002/cpt.337](https://doi.org/10.1002/cpt.337)
- <span id="ref-Zakiyanov2013"></span>Zakiyanov O, Kriha V, Vachek J, Zima T, Tesar V, Kalousova M. Placental growth factor, pregnancy-associated plasma protein-A, soluble receptor for advanced glycation end products, extracellular newly identified receptor for receptor for advanced glycation end products binding protein and high mobility group box 1 levels in patients with acute kidney injury: a cross sectional study. BMC Nephrol 14:245, 2013. [doi:10.1186/1471-2369-14-245](https://doi.org/10.1186/1471-2369-14-245)
- <span id="ref-ZhangX2024"></span>Zhang X, Ye C, Lu F, Yang J, Xu Y, Wang C. Evaluation of renal oxygenation and perfusion in patients with chronic kidney disease: a preliminary prospective study based on functional magnetic resonance. Ren Fail 46(2):2428337, 2024. [doi:10.1080/0886022X.2024.2428337](https://doi.org/10.1080/0886022X.2024.2428337)
<!-- end generated -->
