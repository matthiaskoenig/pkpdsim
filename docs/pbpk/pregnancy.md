# Pregnancy

The physiology of a woman changes over the pregnancy: the cardiac output and the blood volume rise, the blood is diluted, the GFR rises by a third, the albumin falls by a quarter and the activity of several enzymes changes. A model with a pregnancy (`PBPKDefinition(sex=Sex.FEMALE, pregnancy=Pregnancy())`, `pkpdsim.pbpk.pregnancy`, #134) knows the gestational age, and with the fetoplacental unit it gives the fetal exposure next to the maternal one.

- `GA` is the gestational age (weeks since the first day of the last menstrual period, 0 not pregnant, 40 at term): an initial value as the covariates, which an experiment sets before the simulation (`preinit_changes`, or `init(GA)` in roadrunner). `Pregnancy(ga=30.0)` is the default of the build.
- Every change of the maternal physiology is a factor `fpreg_<key>` of the published fits (`database/tables/pregnancy_human.toml`), 1 at GA 0 and linear between the gestational ages of the table. It multiplies its quantity as the factors of the impairment do: the cardiac output and with it every organ flow, the blood volume, the hematocrit, the GFR, the kidney flow, the plasma albumin and alpha-1-acid glycoprotein (so fu of a bound substance rises), the small intestinal transit and the activity of the enzymes and transporters of the liver and the kidney.
- The body weight `BW` stays the body weight before the pregnancy, from which the organs of the woman are scaled. The masses of the pregnancy (fetus, placenta, amniotic fluid and the gain of the uterus and of the breasts) are added to `BW_total`, the weight of the pregnant woman, and take no drug except the fetus and the placenta.
- The fits of Dallmann 2017 and Dallmann 2018 are of the fertilization age, two weeks after the last menstrual period, so every factor is 1 up to GA 2.

<!-- generated: pregnancy-factors (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/pregnancy_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

| quantity | factor | unit | 0 wk | 8 wk | 13 wk | 20 wk | 26 wk | 33 wk | 40 wk | sources |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| CO | `fpreg_CO` | - | 1 | 1.149 | 1.221 | 1.29 | 1.332 | 1.37 | 1.401 | `Dallmann2017` |
| upla | `fpreg_upla` | - | 0 | 0.00956 | 0.02356 | 0.05395 | 0.08355 | 0.1077 | 0.1186 | `Dallmann2017` |
| Vblood | `fpreg_Vblood` | - | 1 | 1.046 | 1.132 | 1.278 | 1.391 | 1.469 | 1.45 | `Dallmann2017` |
| HCT | `fpreg_HCT` | - | 1 | 0.9905 | 0.9671 | 0.9263 | 0.8976 | 0.8901 | 0.9289 | `Dallmann2017` |
| GFR | `fpreg_GFR` | - | 1 | 1.279 | 1.333 | 1.327 | 1.327 | 1.327 | 1.327 | `Dallmann2017` |
| Qki | `fpreg_Qki` | - | 1 | 1.299 | 1.277 | 1.175 | 1.132 | 1.098 | 0.8577 | `Dallmann2017` |
| albumin | `fpreg_albumin` | - | 1 | 0.9245 | 0.8757 | 0.8231 | 0.7897 | 0.7608 | 0.7396 | `Dallmann2017` |
| agp | `fpreg_agp` | - | 1 | 0.9549 | 0.9235 | 0.8884 | 0.8669 | 0.852 | 0.8476 | `Dallmann2017` |
| t_si | `fpreg_t_si` | - | 1 | 1 | 1.278 | 1.667 | 1.638 | 1.605 | 1.6 | `Lawson1985`, `Chiloiro2001`, `Szilagyi1996` |
| mass fetus | `fpreg_mass_fetus` | kg | 0 | 0.002353 | 0.03689 | 0.2969 | 0.9445 | 2.202 | 3.514 | `Dallmann2017` |
| mass placenta | `fpreg_mass_placenta` | kg | 0 | 0.02003 | 0.0762 | 0.2262 | 0.3917 | 0.5711 | 0.7077 | `Dallmann2017` |
| mass amniotic_fluid | `fpreg_mass_amniotic_fluid` | kg | 0 | 0.009208 | 0.08591 | 0.3278 | 0.5614 | 0.8194 | 0.6451 | `Dallmann2017` |
| mass uterus | `fpreg_mass_uterus` | kg | 0 | 0.05517 | 0.1565 | 0.4049 | 0.6466 | 0.8233 | 0.8948 | `Dallmann2017` |
| mass breasts | `fpreg_mass_breasts` | kg | 0 | 0.04818 | 0.09087 | 0.1546 | 0.2133 | 0.287 | 0.3664 | `Dallmann2017` |
| ABCB1 kidney | `fpreg_ki_abcb1` | - | 1 | 1.188 | 1.345 | 1.565 | 1.753 | 1.879 | 1.879 | `Hebert2008` |
| SLC22A2 kidney | `fpreg_ki_slc22a2` | - | 1 | 1.178 | 1.326 | 1.534 | 1.45 | 1.353 | 1.339 | `Eyal2010` |
| SLC22A6 kidney | `fpreg_ki_slc22a6` | - | 1 | 1.072 | 1.132 | 1.216 | 1.288 | 1.3 | 1.3 | `Thakur2026` |
| SLC22A8 kidney | `fpreg_ki_slc22a8` | - | 1 | 1.072 | 1.132 | 1.216 | 1.288 | 1.3 | 1.3 | `Thakur2026` |
| CYP1A2 liver | `fpreg_li_cyp1a2` | - | 1 | 0.8443 | 0.7308 | 0.5957 | 0.5028 | 0.4212 | 0.3676 | `Dallmann2018`, `Tracy2005` |
| CYP2A6 liver | `fpreg_li_cyp2a6` | - | 1 | 1.154 | 1.282 | 1.461 | 1.615 | 1.794 | 1.82 | `Dallmann2018`, `Dempsey2002` |
| CYP2B6 liver | `fpreg_li_cyp2b6` | - | 1 | 1.1 | 1.225 | 1.4 | 1.614 | 1.864 | 1.9 | `Ke2014` |
| CYP2C19 liver | `fpreg_li_cyp2c19` | - | 1 | 0.7933 | 0.6211 | 0.38 | 0.3543 | 0.3243 | 0.32 | `Ke2014` |
| CYP2C9 liver | `fpreg_li_cyp2c9` | - | 1 | 1.4 | 1.442 | 1.5 | 1.543 | 1.593 | 1.6 | `Ke2014` |
| CYP2D6 liver | `fpreg_li_cyp2d6` | - | 1 | 1.11 | 1.201 | 1.293 | 1.348 | 1.424 | 1.478 | `Tracy2005` |
| CYP3A4 liver | `fpreg_li_cyp3a4` | - | 1 | 1.3 | 1.55 | 1.6 | 1.6 | 1.6 | 1.6 | `Dallmann2018`, `Tracy2005`, `Hebert2008` |
| UGT1A1 liver | `fpreg_li_ugt1a1` | - | 1 | 1.24 | 1.407 | 1.457 | 1.5 | 1.55 | 1.6 | `Fischer2014` |
| UGT1A4 liver | `fpreg_li_ugt1a4` | - | 1 | 1.703 | 1.912 | 2.098 | 2.202 | 2.288 | 2.351 | `Karanam2025` |

??? note "Sources"

    - `Chiloiro2001`: Chiloiro M, Darconza G, Piccioli E, De Carne M, Clemente C, Riezzo G. Gastric emptying and orocecal transit time in pregnancy. J Gastroenterol 36(8):538-543, 2001. [doi:10.1007/s005350170056](https://doi.org/10.1007/s005350170056)
    - `Dallmann2017`: Dallmann A, Ince I, Meyer M, Willmann S, Eissing T, Hempel G. Gestation-specific changes in the anatomy and physiology of healthy pregnant women: an extended repository of model parameters for physiologically based pharmacokinetic modeling in pregnancy. Clin Pharmacokinet 56(11):1303-1330, 2017. [doi:10.1007/s40262-017-0539-z](https://doi.org/10.1007/s40262-017-0539-z)
    - `Dallmann2018`: Dallmann A, Ince I, Coboeken K, Eissing T, Hempel G. A physiologically based pharmacokinetic model for pregnant women to predict the pharmacokinetics of drugs metabolized via several enzymatic pathways. Clin Pharmacokinet 57(6):749-768, 2018. [doi:10.1007/s40262-017-0594-5](https://doi.org/10.1007/s40262-017-0594-5)
    - `Dempsey2002`: Dempsey D, Jacob P 3rd, Benowitz NL. Accelerated metabolism of nicotine and cotinine in pregnant smokers. J Pharmacol Exp Ther 301(2):594-598, 2002. [doi:10.1124/jpet.301.2.594](https://doi.org/10.1124/jpet.301.2.594)
    - `Eyal2010`: Eyal S, Easterling TR, Carr D, Umans JG, Miodovnik M, Hankins GD, Clark SM, Risler L, Wang J, Kelly EJ, Shen DD, Hebert MF. Pharmacokinetics of metformin during pregnancy. Drug Metab Dispos 38(5):833-840, 2010. [doi:10.1124/dmd.109.031245](https://doi.org/10.1124/dmd.109.031245)
    - `Fischer2014`: Fischer JH, Sarto GE, Hardman J, Endres L, Jenkins TM, Kilpatrick SJ, Jeong H, Geller S, Deyo K, Fischer PA, Rodvold KA. Influence of gestational age and body weight on the pharmacokinetics of labetalol in pregnancy. Clin Pharmacokinet 53(4):373-383, 2014. [doi:10.1007/s40262-013-0123-0](https://doi.org/10.1007/s40262-013-0123-0)
    - `Hebert2008`: Hebert MF, Easterling TR, Kirby B, Carr DB, Buchanan ML, Rutherford T, Thummel KE, Fishbein DP, Unadkat JD. Effects of pregnancy on CYP3A and P-glycoprotein activities as measured by disposition of midazolam and digoxin: a University of Washington Specialized Center of Research study. Clin Pharmacol Ther 84(2):248-253, 2008. [doi:10.1038/clpt.2008.1](https://doi.org/10.1038/clpt.2008.1)
    - `Karanam2025`: Karanam A, Pennell PB, Meador KJ, Long Y, Birnbaum AK. Characterization of lamotrigine disposition changes during and after pregnancy in women with epilepsy. Pharmacotherapy 45(1):33-42, 2025. [doi:10.1002/phar.4640](https://doi.org/10.1002/phar.4640)
    - `Ke2014`: Ke AB, Nallani SC, Zhao P, Rostami-Hodjegan A, Unadkat JD. Expansion of a PBPK model to predict disposition in pregnant women of drugs cleared via multiple CYP enzymes, including CYP2B6, CYP2C9 and CYP2C19. Br J Clin Pharmacol 77(3):554-570, 2014. [doi:10.1111/bcp.12207](https://doi.org/10.1111/bcp.12207)
    - `Lawson1985`: Lawson M, Kern F Jr, Everson GT. Gastrointestinal transit time in human pregnancy: prolongation in the second and third trimesters followed by postpartum normalization. Gastroenterology 89(5):996-999, 1985. [doi:10.1016/0016-5085(85)90199-4](https://doi.org/10.1016/0016-5085(85)90199-4)
    - `Szilagyi1996`: Szilagyi A, Salomon R, Smith BE, Martin M, Seidman E. Determinants of prolonged oral cecal transit time during late phase pregnancy. Clin Invest Med 19(1):20-27, 1996.
    - `Thakur2026`: Thakur A, Unadkat JD, Fay EE, Isoherranen N, Hebert MF. Effect of late second to early third trimester of pregnancy on the activity of renal organic anion transporters (OAT1 and OAT3): a biomarker study. Clin Pharmacol Ther 120(2):465-474, 2026. [doi:10.1002/cpt.70321](https://doi.org/10.1002/cpt.70321)
    - `Tracy2005`: Tracy TS, Venkataramanan R, Glover DD, Caritis SN. Temporal changes in drug metabolism (CYP1A2, CYP2D6 and CYP3A activity) during pregnancy. Am J Obstet Gynecol 192(2):633-639, 2005. [doi:10.1016/j.ajog.2004.08.030](https://doi.org/10.1016/j.ajog.2004.08.030)
<!-- end generated -->

## The fetoplacental unit

`Pregnancy(fetus=True)` (the default) adds the placenta and the fetus as two compartments of the top model; `Pregnancy(fetus=False)` gives the maternal physiology alone.

```equations title="the placenta and the fetus"
Qpla = CO_base * fpreg_Qpla
J_pla_in_<s>: Car_<s> -> Cpla_<s> = Qpla * (1 - HCT) * Car_<s>
J_pla_out_<s>: Cpla_<s> -> Cve_<s> = Qpla * (1 - HCT) * Cpla_<s> / Kp_pla_<s>
J_fet_<s>: Cpla_<s> -> Cfet_<s> = PS_pla_<s> * (fu_<s> * Cpla_<s> / Kp_pla_<s> - fu_<s> * Cfet_<s> / Kp_fet_<s>)
```

- The placenta `Vpla` (the mass of the table at 1 kg/l) is perfused with `Qpla`, the gain of the bilateral uterine arterial blood flow of the pregnancy (`fpreg_upla` of the cardiac output, 0.88 l/min at term): it takes the substance from the arterial blood and returns it to the venous blood, so the cardiac output is the flows of the organs and of the placenta. `Kp_pla_<s>` is the placenta to plasma coefficient, 1 by default.
- The fetus `Vfet` (3.5 l at term) exchanges with the placenta with the transfer clearance `PS_pla_<s>` on the unbound concentrations; `Kp_fet_<s>` is its tissue to plasma coefficient, 1 by default. `PS_pla_<s>` is the plasma flow of the placenta by default, so the transfer is limited by the flow and the fetal concentration follows the maternal one: a model fits it or sets it from a placenta perfusion study.
- `J_fet_<s>` is reversible (the flux is signed: the fetus gives the substance back when the maternal concentration falls). `Cfet_<s>` is the concentration of the fetus, the model of the cord blood, and the fetal to maternal ratio follows from it.

The same oral dose at four gestational ages: the maternal exposure falls (haemodilution, the induced CYP3A4 and the higher GFR), and the fetal concentration follows the maternal arterial one with the default flow-limited transfer.

```python
from pathlib import Path

import numpy as np
from sbmlsim.model.model_roadrunner import load_roadrunner

from pkpdsim.pbpk import (
    MM,
    Binding,
    Excretion,
    Kp,
    Metabolism,
    Organ,
    PBPKDefinition,
    Pregnancy,
    Route,
    Sex,
    Substance,
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
    sid="pregnant",
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
    routes=[Route.oral("drg")],
    sex=Sex.FEMALE,
    pregnancy=Pregnancy(),
)
output = Path("models")
output.mkdir()
result = build(definition, output)

r = load_roadrunner(result.flat_path)
for ga in (0.0, 13.0, 26.0, 40.0):
    r.setValue("init(GA)", ga)
    r.resetAll()
    r.timeCourseSelections = ["time", "Carm_drg", "[Cfet_drg]", "[Car_drg]"]
    r["PODOSE_drg"] = 10.0  # mg
    s = r.simulate(0, 48 * 60, 2881)
    carm = s["Carm_drg"] * 1000  # µM
    auc = float(np.trapezoid(carm, s["time"])) / 60  # µM*h
    ratio = float(np.trapezoid(s["[Cfet_drg]"], s["time"])) / float(
        np.trapezoid(s["[Car_drg]"], s["time"])
    )
    print(
        f"GA {ga:4.0f} weeks: Cmax {carm.max():.3f} µM, AUC(48 h) {auc:6.2f} µM*h, "
        f"fetal over maternal arterial AUC {ratio:.2f}"
    )
```

```text
GA    0 weeks: Cmax 2.151 µM, AUC(48 h)  13.53 µM*h, fetal over maternal arterial AUC 0.00
GA   13 weeks: Cmax 1.860 µM, AUC(48 h)  12.11 µM*h, fetal over maternal arterial AUC 1.00
GA   26 weeks: Cmax 1.533 µM, AUC(48 h)  11.77 µM*h, fetal over maternal arterial AUC 1.00
GA   40 weeks: Cmax 1.281 µM, AUC(48 h)  11.28 µM*h, fetal over maternal arterial AUC 0.99
```

<!-- generated: figures-pregnancy-1 (scripts/visual/docs_figures.py), do not edit -->
The content of the model `pregnant`: the anatomical figure (A) and the schematic (B), hover them for details ([figures](figures.md)).

<iframe class="pkpd-figure" src="img/figures/models/pregnancy-1-anatomy.html" title="The anatomical figure of the model pregnant" style="aspect-ratio:1460 / 898" loading="lazy"></iframe>

<iframe class="pkpd-figure" src="img/figures/models/pregnancy-1-schematic.html" title="The schematic of the model pregnant" style="aspect-ratio:1384 / 902" loading="lazy"></iframe>
<!-- end generated -->

## No change

The searches found no change of the pregnancy, or no number which could be used, for the gastric emptying (four scintigraphy and paracetamol studies), the gastric pH, UGT2B7 and the P-glycoprotein of the intestine; no data for the colonic transit and the other transporters. Not applied: the liver volume (+15 % by MRI, Bartlett 2021; the in vivo factors of the enzymes would count the growth twice), the kidney volume (the one study measures the dilated collecting system as well) and the body water of Abduljalil 2012 (it counts the fetus).

## Limitations

- The fetus neither metabolises nor excretes the substance: it is one well-stirred compartment which exchanges with the placenta alone, so its only way out is back through the placenta. The amniotic fluid takes no drug, and the placenta has no metabolism and no transporters of its own.
- `PS_pla_<s>` of the plasma flow of the placenta is a choice of the library, not a measurement: for a substance which crosses the placenta slowly (a large or a charged molecule, a P-glycoprotein substrate) a model has to set or fit it, otherwise the fetal exposure is over-predicted.
- The factors of the enzymes are the in vivo clearances of probe drugs against postpartum, so they carry the changes of the liver and of the binding of those studies; CYP2C9, CYP2C19, CYP2B6 and UGT1A1 rest on one study each (low confidence). The GFR of Dallmann 2017 (Eq. 27) is printed with two misprints, which the library reads as 0.0245, 113 and 150 ml/min, checked against the text and the tabulated curve of the PK-Sim database.
- The kidney keeps its measured renal blood flow while the other organs share the rest of the cardiac output, so the sum of the flows is about 2 % above the meta-analysed cardiac output at term (the cardiac output of the model is the sum of the flows).
- One pregnancy of one fetus, a healthy woman, and no labour: pre-eclampsia, twins, the postpartum period and the changes of the fetus itself over the trimesters (its binding proteins, its enzymes) are not modelled, nor is the milk of the lactation.
