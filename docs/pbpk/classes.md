# Impairment classes

The models describe impaired organ function by the clinical classes the studies report: the Child-Pugh class of liver cirrhosis, the GFR category of chronic kidney disease (CKD) and the NYHA class of heart failure. This page defines the classes, shows which value of a model parameter represents each class, and how a model scans and draws them. The factors which a class sets on the physiology are on the page [Impairment](impairment.md), the studies they come from in the [parameter database](../database.md).

Every class has one value of one parameter of the top model, set before the initialization like a covariate:

| condition | classification | the value of a class |
| --- | --- | --- |
| liver cirrhosis | Child-Pugh classes A, B, C (Pugh 1973) | `CP_severity`: healthy 0, A 1, B 2, C 3; a MELD score maps onto it ([The MELD score](impairment.md#the-meld-score)) |
| chronic kidney disease | KDIGO GFR categories G1 to G5 (KDIGO 2012) | `f_renal_function`: the eGFR of the category over the eGFR of the build, G1 the renal function 1 |
| heart failure | NYHA functional classes I to IV | `HF_severity`: healthy 0, I 1, II 2, III 3, IV 4 |

A value between two classes interpolates the factors linearly, so a scan or a fit can use any severity. The classifications are `CHILD_PUGH_CLASSES`, `KDIGO_GFR_CLASSES` and `NYHA_CLASSES` of `pkpdsim.impairment` (`CLASSIFICATIONS` by key): every class has its key, a label, a color for every figure, its definition and the value which represents it.

<!-- generated: classes-figure (scripts/pbpk_tables.py), do not edit -->
<div class="pbpk-figure pbpk-wide" markdown="0">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 362" role="img" aria-label="The clinical classes of the impairments and their scan values" font-family="system-ui, sans-serif" font-size="14" fill="currentColor">
<title>The clinical classes of the impairments and their scan values</title>
<text x="10" y="26" font-weight="600">Child-Pugh classes by the score 5-15, the CP_severity of each below</text>
<rect x="11" y="36" width="58" height="22" rx="3" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy: no liver disease</title></rect>
<text x="40" y="77" text-anchor="middle">healthy</text>
<rect x="71" y="36" width="114.4" height="22" rx="3" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A: well-compensated cirrhosis, score 5-6</title></rect>
<text x="128.2" y="77" text-anchor="middle">A</text>
<rect x="187.4" y="36" width="172.5" height="22" rx="3" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B: significant functional compromise, score 7-9</title></rect>
<text x="273.6" y="77" text-anchor="middle">B</text>
<rect x="361.9" y="36" width="347.1" height="22" rx="3" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C: decompensated cirrhosis, score 10-15</title></rect>
<text x="535.5" y="77" text-anchor="middle">C</text>
<line x1="40" x2="40" y1="33" y2="61" stroke="currentColor" stroke-width="2"/>
<text x="40" y="94" text-anchor="middle" font-size="12" opacity="0.8">0</text>
<line x1="128.2" x2="128.2" y1="33" y2="61" stroke="currentColor" stroke-width="2"/>
<text x="128.2" y="94" text-anchor="middle" font-size="12" opacity="0.8">1</text>
<line x1="273.6" x2="273.6" y1="33" y2="61" stroke="currentColor" stroke-width="2"/>
<text x="273.6" y="94" text-anchor="middle" font-size="12" opacity="0.8">2</text>
<line x1="535.5" x2="535.5" y1="33" y2="61" stroke="currentColor" stroke-width="2"/>
<text x="535.5" y="94" text-anchor="middle" font-size="12" opacity="0.8">3</text>
<text x="10" y="116" font-size="12" opacity="0.8">MELD of the classes (consensus of the database): healthy 6.43, A 7.91, B 12.8, C 20.5</text>
<text x="10" y="144" font-weight="600">KDIGO GFR categories by the eGFR (120 to 0 ml/min/1.73 m2), the eGFR of each below</text>
<rect x="11" y="154" width="173" height="22" rx="3" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal: normal or high, eGFR >= 90</title></rect>
<text x="97.5" y="195" text-anchor="middle">G1</text>
<rect x="186" y="154" width="173" height="22" rx="3" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild: mildly decreased, eGFR 60-89</title></rect>
<text x="272.5" y="195" text-anchor="middle">G2</text>
<rect x="361" y="154" width="85.5" height="22" rx="3" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate: mildly to moderately decreased, eGFR 45-59</title></rect>
<text x="403.8" y="195" text-anchor="middle">G3a</text>
<rect x="448.5" y="154" width="85.5" height="22" rx="3" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe: moderately to severely decreased, eGFR 30-44</title></rect>
<text x="491.2" y="195" text-anchor="middle">G3b</text>
<rect x="536" y="154" width="85.5" height="22" rx="3" fill="#006d2c" fill-opacity="0.9"><title>G4 severe: severely decreased, eGFR 15-29</title></rect>
<text x="578.8" y="195" text-anchor="middle">G4</text>
<rect x="623.5" y="154" width="85.5" height="22" rx="3" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure: kidney failure, eGFR < 15 (not on dialysis in the library)</title></rect>
<text x="666.2" y="195" text-anchor="middle">G5</text>
<line x1="272.5" x2="272.5" y1="151" y2="179" stroke="currentColor" stroke-width="2"/>
<text x="272.5" y="212" text-anchor="middle" font-size="12" opacity="0.8">75</text>
<line x1="406.7" x2="406.7" y1="151" y2="179" stroke="currentColor" stroke-width="2"/>
<text x="406.7" y="212" text-anchor="middle" font-size="12" opacity="0.8">52</text>
<line x1="494.2" x2="494.2" y1="151" y2="179" stroke="currentColor" stroke-width="2"/>
<text x="494.2" y="212" text-anchor="middle" font-size="12" opacity="0.8">37</text>
<line x1="578.8" x2="578.8" y1="151" y2="179" stroke="currentColor" stroke-width="2"/>
<text x="578.8" y="212" text-anchor="middle" font-size="12" opacity="0.8">22.5</text>
<line x1="666.2" x2="666.2" y1="151" y2="179" stroke="currentColor" stroke-width="2"/>
<text x="666.2" y="212" text-anchor="middle" font-size="12" opacity="0.8">7.5</text>
<text x="10" y="234" font-size="12" opacity="0.8">the eGFR of the scan of a category; G1 is the renal function 1 of the build (eGFR 114.5 male, 114.7 female)</text>
<text x="10" y="262" font-weight="600">NYHA functional classes, the HF_severity of each below</text>
<rect x="11" y="272" width="138" height="22" rx="3" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy: no heart failure</title></rect>
<text x="80" y="313" text-anchor="middle">healthy</text>
<rect x="151" y="272" width="138" height="22" rx="3" fill="#FDDBC7" fill-opacity="0.9"><title>NYHA I: no limitation of physical activity; ordinary activity causes no undue fatigue, palpitation or dyspnoea</title></rect>
<text x="220" y="313" text-anchor="middle">I</text>
<rect x="291" y="272" width="138" height="22" rx="3" fill="#F4A582" fill-opacity="0.9"><title>NYHA II: slight limitation; comfortable at rest, ordinary activity causes fatigue, palpitation or dyspnoea</title></rect>
<text x="360" y="313" text-anchor="middle">II</text>
<rect x="431" y="272" width="138" height="22" rx="3" fill="#D6604D" fill-opacity="0.9"><title>NYHA III: marked limitation; comfortable at rest, less than ordinary activity causes symptoms</title></rect>
<text x="500" y="313" text-anchor="middle">III</text>
<rect x="571" y="272" width="138" height="22" rx="3" fill="#67001F" fill-opacity="0.9"><title>NYHA IV: unable to carry on any physical activity without discomfort; symptoms at rest</title></rect>
<text x="640" y="313" text-anchor="middle">IV</text>
<line x1="80" x2="80" y1="269" y2="297" stroke="currentColor" stroke-width="2"/>
<text x="80" y="330" text-anchor="middle" font-size="12" opacity="0.8">0</text>
<line x1="220" x2="220" y1="269" y2="297" stroke="currentColor" stroke-width="2"/>
<text x="220" y="330" text-anchor="middle" font-size="12" opacity="0.8">1</text>
<line x1="360" x2="360" y1="269" y2="297" stroke="currentColor" stroke-width="2"/>
<text x="360" y="330" text-anchor="middle" font-size="12" opacity="0.8">2</text>
<line x1="500" x2="500" y1="269" y2="297" stroke="currentColor" stroke-width="2"/>
<text x="500" y="330" text-anchor="middle" font-size="12" opacity="0.8">3</text>
<line x1="640" x2="640" y1="269" y2="297" stroke="currentColor" stroke-width="2"/>
<text x="640" y="330" text-anchor="middle" font-size="12" opacity="0.8">4</text>
<text x="10" y="352" font-size="12" opacity="0.8">the classes grade the symptoms at physical activity</text>
</svg>
</div>
<!-- end generated -->

## Liver cirrhosis: Child-Pugh

The Child-Pugh score (Pugh 1973, a modification of Child and Turcotte 1964) adds 1 to 3 points for each of five criteria, two clinical and three laboratory values, from 5 to 15 points; the class is A for 5-6, B for 7-9 and C for 10-15 points. The FDA and EMA guidances on hepatic impairment grade the subjects of a pharmacokinetic study by it (mild A, moderate B, severe C). The score has no healthy class: a subject without liver disease has 5 points, the reference of the models (`CP_severity` 0).

<!-- generated: classes-child-pugh (scripts/pbpk_tables.py), do not edit -->
*Generated from `pkpdsim.impairment.CHILD_PUGH_CLASSES` by `scripts/pbpk_tables.py`.*

| criterion | 1 point | 2 points | 3 points |
| --- | --- | --- | --- |
| total bilirubin | < 2 mg/dL (< 34 umol/L) | 2-3 mg/dL (34-50 umol/L) | > 3 mg/dL (> 50 umol/L) |
| serum albumin | > 3.5 g/dL (> 35 g/L) | 2.8-3.5 g/dL (28-35 g/L) | < 2.8 g/dL (< 28 g/L) |
| prothrombin time: INR | < 1.7 | 1.7-2.3 | > 2.3 |
| (or its prolongation) | < 4 s | 4-6 s | > 6 s |
| ascites | none | mild (controlled by diuretics) | moderate to severe (refractory) |
| hepatic encephalopathy | none | grade I-II (or controlled by medication) | grade III-IV (refractory) |

| class | definition | score | CP_severity |
| --- | --- | --- | ---: |
| <span class="pbpk-swatch" style="background:black" title="black"></span> `healthy` Healthy | no liver disease |  | 0 |
| <span class="pbpk-swatch" style="background:#74a9cf" title="#74a9cf"></span> `A` Child-Pugh A | well-compensated cirrhosis, score 5-6 | 5-6 | 1 |
| <span class="pbpk-swatch" style="background:#2b8cbe" title="#2b8cbe"></span> `B` Child-Pugh B | significant functional compromise, score 7-9 | 7-9 | 2 |
| <span class="pbpk-swatch" style="background:#045a8d" title="#045a8d"></span> `C` Child-Pugh C | decompensated cirrhosis, score 10-15 | 10-15 | 3 |
<!-- end generated -->

The MELD score (model for end-stage liver disease, Kamath 2001) is continuous: `3.78 ln(bilirubin [mg/dL]) + 11.2 ln(INR) + 9.57 ln(creatinine [mg/dL]) + 6.43`, every value bounded below by 1 (UNOS also bounds the creatinine above by 4 mg/dL, 4 on dialysis, and multiplies by 10 and rounds the score of `0.957 ln(Cr) + 0.378 ln(bili) + 1.120 ln(INR) + 0.643`). The score of a healthy subject is 6.43. The library maps it onto `CP_severity` through the MELD of each class (`cp_severity_of_meld`, [The MELD score](impairment.md#the-meld-score)).

## Chronic kidney disease: KDIGO

CKD is an abnormality of the structure or the function of the kidneys present for more than 3 months (KDIGO 2012). It is graded by the cause, the GFR category and the albuminuria category. The GFR categories by the GFR per 1.73 m2 body surface area (the eGFR, estimated from creatinine or cystatin C, or measured) are G1 (>= 90, normal or high), G2 (60-89, mildly decreased), G3a (45-59, mildly to moderately decreased), G3b (30-44, moderately to severely decreased), G4 (15-29, severely decreased) and G5 (< 15, kidney failure, G5D on dialysis); G1 and G2 are CKD only with a marker of kidney damage. The albuminuria categories by the albumin excretion rate (mg/24 h) or the albumin-to-creatinine ratio (mg/g) are A1 (< 30, normal to mildly increased), A2 (30-300, moderately increased) and A3 (> 300, severely increased); the library has no albuminuria. The FDA guidance on renal impairment (2024) has the stages normal (>= 90), mild (60-89), moderate (30-59, G3a and G3b), severe (15-29) and kidney failure (< 15, not on dialysis, and on dialysis apart), by the eGFR or the creatinine clearance; the EMA guidance (2015) has the same limits.

The library has the stages of the FDA: a vector of a CKD factor gives mild, moderate, severe and end stage, and the factor is piecewise linear in the eGFR through the stage midpoints 75, 45, 22.5 and 7.5 ml/min/1.73 m2, 1 from 90 up ([Impairment](impairment.md#chronic-kidney-disease)). The eGFR of a model is `f_renal_function` times the eGFR of the normal renal function of its build (the GFR per 1.73 m2 of the reference individual, 114.5 ml/min/1.73 m2 male, 114.7 female). A GFR category is represented by the eGFR at its midpoint (G5 at 7.5, the anchor of the end stage), G1 by the renal function 1 of the build; the end stage of the library is not on dialysis.

<!-- generated: classes-kdigo (scripts/pbpk_tables.py), do not edit -->
*Generated from `pkpdsim.impairment.KDIGO_GFR_CLASSES` by `scripts/pbpk_tables.py`; the eGFR in ml/min/1.73 m2, `f_renal_function` of the male / female build.*

| class | definition | FDA stage | eGFR of the scan | f_renal_function |
| --- | --- | --- | --- | --- |
| <span class="pbpk-swatch" style="background:black" title="black"></span> `G1` G1 normal | normal or high, eGFR >= 90 | normal | reference | 1 / 1 |
| <span class="pbpk-swatch" style="background:#66c2a4" title="#66c2a4"></span> `G2` G2 mild | mildly decreased, eGFR 60-89 | mild | 75 | 0.6552 / 0.6539 |
| <span class="pbpk-swatch" style="background:#41ae76" title="#41ae76"></span> `G3a` G3a mild-moderate | mildly to moderately decreased, eGFR 45-59 | moderate | 52 | 0.4543 / 0.4534 |
| <span class="pbpk-swatch" style="background:#238b45" title="#238b45"></span> `G3b` G3b moderate-severe | moderately to severely decreased, eGFR 30-44 | moderate | 37 | 0.3232 / 0.3226 |
| <span class="pbpk-swatch" style="background:#006d2c" title="#006d2c"></span> `G4` G4 severe | severely decreased, eGFR 15-29 | severe | 22.5 | 0.1966 / 0.1962 |
| <span class="pbpk-swatch" style="background:#00441b" title="#00441b"></span> `G5` G5 kidney failure | kidney failure, eGFR < 15 (not on dialysis in the library) | end stage | 7.5 | 0.06552 / 0.06539 |
<!-- end generated -->

The renal function of a study which reports a measured or estimated clearance in ml/min is `renal_function(clearance, bw=..., height=...)` of the experiment, so the GFR of the simulated individual is the measured one; the classes of the earlier models keep their GFR in ml/min (`renal_map`, [Classes in experiments](impairment.md#classes-in-experiments)).

## Heart failure: NYHA

The NYHA functional classes (Criteria Committee of the New York Heart Association, 9th edition 1994) grade the symptoms of heart failure at physical activity. The ACC/AHA stages grade its development instead: A at risk without structural heart disease or symptoms, B structural heart disease without symptoms (pre-heart failure), C structural heart disease with current or previous symptoms and D advanced heart failure refractory to treatment; the NYHA classes grade the patients of the stages C and D, a patient of stage D is mostly NYHA IV. The ejection fraction separates heart failure with reduced (<= 40 %), mildly reduced (41-49 %) and preserved (>= 50 %) ejection fraction; the library has no such distinction.

<!-- generated: classes-nyha (scripts/pbpk_tables.py), do not edit -->
*Generated from `pkpdsim.impairment.NYHA_CLASSES` by `scripts/pbpk_tables.py`.*

| class | definition | HF_severity |
| --- | --- | ---: |
| <span class="pbpk-swatch" style="background:black" title="black"></span> `healthy` Healthy | no heart failure | 0 |
| <span class="pbpk-swatch" style="background:#FDDBC7" title="#FDDBC7"></span> `I` NYHA I | no limitation of physical activity; ordinary activity causes no undue fatigue, palpitation or dyspnoea | 1 |
| <span class="pbpk-swatch" style="background:#F4A582" title="#F4A582"></span> `II` NYHA II | slight limitation; comfortable at rest, ordinary activity causes fatigue, palpitation or dyspnoea | 2 |
| <span class="pbpk-swatch" style="background:#D6604D" title="#D6604D"></span> `III` NYHA III | marked limitation; comfortable at rest, less than ordinary activity causes symptoms | 3 |
| <span class="pbpk-swatch" style="background:#67001F" title="#67001F"></span> `IV` NYHA IV | unable to carry on any physical activity without discomfort; symptoms at rest | 4 |
<!-- end generated -->

## The factors over the classes

The physiological factors of each condition over its severity, a panel per factor: the line is the factor of the model, the points are the classes, the strip under a panel the classes of the axis in their colors (hover a point for its value). The factors of the enzymes and transporters of the liver are the tables of [Impairment](impairment.md).

<!-- generated: factors-child-pugh (scripts/pbpk_tables.py), do not edit -->
<div class="pbpk-figure pbpk-wide" markdown="0">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 730 364" role="img" aria-label="Liver cirrhosis: the factors over CP_severity (healthy, A, B, C)" font-family="system-ui, sans-serif" font-size="14" fill="currentColor">
<title>Liver cirrhosis: the factors over CP_severity (healthy, A, B, C)</title>
<text x="10" y="18" font-weight="600">Liver cirrhosis: the factors over CP_severity (healthy, A, B, C)</text>
<text x="10" y="45" font-size="13">functional_liver</text>
<rect x="10" y="125" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="44" y="125" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="78" y="125" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="112" y="125" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="10" x2="146" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="27,57.19 61,77.28 95,86.35 129,103.9" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="77.28" r="3.5" fill="currentColor"><title>A: 0.69</title></circle>
<circle cx="95" cy="86.35" r="3.5" fill="currentColor"><title>B: 0.55</title></circle>
<circle cx="129" cy="103.9" r="3.5" fill="currentColor"><title>C: 0.28</title></circle>
<text x="149" y="106.9" font-size="11" opacity="0.7">0.28</text>
<text x="190" y="45" font-size="13">shunt</text>
<rect x="190" y="125" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="224" y="125" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="258" y="125" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="292" y="125" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="190" x2="326" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="329" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="207,122 241,96.11 275,76.76 309,69.32" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="241" cy="96.11" r="3.5" fill="currentColor"><title>A: 0.3995</title></circle>
<circle cx="275" cy="76.76" r="3.5" fill="currentColor"><title>B: 0.698</title></circle>
<circle cx="309" cy="69.32" r="3.5" fill="currentColor"><title>C: 0.8128</title></circle>
<text x="329" y="125" font-size="11" opacity="0.7">0</text>
<text x="370" y="45" font-size="13">hepatic_artery_flow</text>
<rect x="370" y="125" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="404" y="125" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="438" y="125" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="472" y="125" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="370" x2="506" y1="102.9" y2="102.9" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="509" y="105.9" font-size="11" opacity="0.7">1</text>
<polyline points="387,102.9 421,97.22 455,78.15 489,57.19" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="421" cy="97.22" r="3.5" fill="currentColor"><title>A: 1.3</title></circle>
<circle cx="455" cy="78.15" r="3.5" fill="currentColor"><title>B: 2.3</title></circle>
<circle cx="489" cy="57.19" r="3.5" fill="currentColor"><title>C: 3.4</title></circle>
<text x="509" y="60.19" font-size="11" opacity="0.7">3.4</text>
<text x="550" y="45" font-size="13">kidney_flow</text>
<rect x="550" y="125" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="584" y="125" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="618" y="125" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="652" y="125" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="550" x2="686" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="689" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="567,57.19 601,64.96 635,79.87 669,90.89" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="601" cy="64.96" r="3.5" fill="currentColor"><title>A: 0.88</title></circle>
<circle cx="635" cy="79.87" r="3.5" fill="currentColor"><title>B: 0.65</title></circle>
<circle cx="669" cy="90.89" r="3.5" fill="currentColor"><title>C: 0.48</title></circle>
<text x="689" y="93.89" font-size="11" opacity="0.7">0.48</text>
<text x="10" y="155" font-size="13">cardiac_index</text>
<rect x="10" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="44" y="235" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="78" y="235" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="112" y="235" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="10" x2="146" y1="184.3" y2="184.3" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="187.3" font-size="11" opacity="0.7">1</text>
<polyline points="27,184.3 61,179.6 95,174.3 129,167.2" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="179.6" r="3.5" fill="currentColor"><title>A: 1.099</title></circle>
<circle cx="95" cy="174.3" r="3.5" fill="currentColor"><title>B: 1.21</title></circle>
<circle cx="129" cy="167.2" r="3.5" fill="currentColor"><title>C: 1.36</title></circle>
<text x="149" y="170.2" font-size="11" opacity="0.7">1.36</text>
<text x="190" y="155" font-size="13">gfr</text>
<rect x="190" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="224" y="235" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="258" y="235" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="292" y="235" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="190" x2="326" y1="167.2" y2="167.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="329" y="170.2" font-size="11" opacity="0.7">1</text>
<polyline points="207,167.2 241,167.2 275,186.6 309,208.7" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="241" cy="167.2" r="3.5" fill="currentColor"><title>A: 1</title></circle>
<circle cx="275" cy="186.6" r="3.5" fill="currentColor"><title>B: 0.7</title></circle>
<circle cx="309" cy="208.7" r="3.5" fill="currentColor"><title>C: 0.36</title></circle>
<text x="329" y="211.7" font-size="11" opacity="0.7">0.36</text>
<text x="370" y="155" font-size="13">albumin</text>
<rect x="370" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="404" y="235" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="438" y="235" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="472" y="235" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="370" x2="506" y1="167.2" y2="167.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="509" y="170.2" font-size="11" opacity="0.7">1</text>
<polyline points="387,167.2 421,176 455,187.9 489,200.3" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="421" cy="176" r="3.5" fill="currentColor"><title>A: 0.8645</title></circle>
<circle cx="455" cy="187.9" r="3.5" fill="currentColor"><title>B: 0.6798</title></circle>
<circle cx="489" cy="200.3" r="3.5" fill="currentColor"><title>C: 0.4892</title></circle>
<text x="509" y="203.3" font-size="11" opacity="0.7">0.4892</text>
<text x="550" y="155" font-size="13">agp</text>
<rect x="550" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="584" y="235" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="618" y="235" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="652" y="235" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="550" x2="686" y1="167.2" y2="167.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="689" y="170.2" font-size="11" opacity="0.7">1</text>
<polyline points="567,167.2 601,193.1 635,195.7 669,212.6" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="601" cy="193.1" r="3.5" fill="currentColor"><title>A: 0.6</title></circle>
<circle cx="635" cy="195.7" r="3.5" fill="currentColor"><title>B: 0.56</title></circle>
<circle cx="669" cy="212.6" r="3.5" fill="currentColor"><title>C: 0.3</title></circle>
<text x="689" y="215.6" font-size="11" opacity="0.7">0.3</text>
<text x="10" y="265" font-size="13">hematocrit</text>
<rect x="10" y="345" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="44" y="345" width="34" height="5" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A</title></rect>
<rect x="78" y="345" width="34" height="5" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B</title></rect>
<rect x="112" y="345" width="34" height="5" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C</title></rect>
<line x1="10" x2="146" y1="277.2" y2="277.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="280.2" font-size="11" opacity="0.7">1</text>
<polyline points="27,277.2 61,282.2 95,294 129,295.1" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="282.2" r="3.5" fill="currentColor"><title>A: 0.9232</title></circle>
<circle cx="95" cy="294" r="3.5" fill="currentColor"><title>B: 0.7411</title></circle>
<circle cx="129" cy="295.1" r="3.5" fill="currentColor"><title>C: 0.7232</title></circle>
<text x="149" y="298.1" font-size="11" opacity="0.7">0.7232</text>
</svg>
</div>
<!-- end generated -->

<!-- generated: factors-ckd (scripts/pbpk_tables.py), do not edit -->
<div class="pbpk-figure pbpk-wide" markdown="0">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 730 364" role="img" aria-label="Chronic kidney disease: the factors over the eGFR (120 to 0 ml/min/1.73 m2)" font-family="system-ui, sans-serif" font-size="14" fill="currentColor">
<title>Chronic kidney disease: the factors over the eGFR (120 to 0 ml/min/1.73 m2)</title>
<text x="10" y="18" font-weight="600">Chronic kidney disease: the factors over the eGFR (120 to 0 ml/min/1.73 m2)</text>
<text x="10" y="45" font-size="13">plasma_protein</text>
<rect x="10" y="125" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="44" y="125" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="78" y="125" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="95" y="125" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="112" y="125" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="129" y="125" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="10" x2="146" y1="96.07" y2="96.07" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="99.07" font-size="11" opacity="0.7">1</text>
<polyline points="10,96.07 44,96.07 61,96.07 95,84.41 120.5,64.96 137.5,57.19 146,57.19" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="96.07" r="3.5" fill="currentColor"><title>mild: 1</title></circle>
<circle cx="95" cy="84.41" r="3.5" fill="currentColor"><title>moderate: 1.45</title></circle>
<circle cx="120.5" cy="64.96" r="3.5" fill="currentColor"><title>severe: 2.2</title></circle>
<circle cx="137.5" cy="57.19" r="3.5" fill="currentColor"><title>end stage: 2.5</title></circle>
<text x="149" y="60.19" font-size="11" opacity="0.7">2.5</text>
<text x="190" y="45" font-size="13">hematocrit</text>
<rect x="190" y="125" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="224" y="125" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="258" y="125" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="275" y="125" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="292" y="125" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="309" y="125" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="190" x2="326" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="329" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="190,57.19 224,57.19 241,61.4 275,64.59 300.5,70.74 317.5,75.94 326,75.94" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="241" cy="61.4" r="3.5" fill="currentColor"><title>mild: 0.9349</title></circle>
<circle cx="275" cy="64.59" r="3.5" fill="currentColor"><title>moderate: 0.8857</title></circle>
<circle cx="300.5" cy="70.74" r="3.5" fill="currentColor"><title>severe: 0.7909</title></circle>
<circle cx="317.5" cy="75.94" r="3.5" fill="currentColor"><title>end stage: 0.7106</title></circle>
<text x="329" y="78.94" font-size="11" opacity="0.7">0.7106</text>
<text x="370" y="45" font-size="13">kidney_volume</text>
<rect x="370" y="125" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="404" y="125" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="438" y="125" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="455" y="125" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="472" y="125" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="489" y="125" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="370" x2="506" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="509" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="370,57.19 404,57.19 421,57.19 455,57.19 480.5,57.19 497.5,57.19 506,57.19" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="421" cy="57.19" r="3.5" fill="currentColor"><title>mild: 1</title></circle>
<circle cx="455" cy="57.19" r="3.5" fill="currentColor"><title>moderate: 1</title></circle>
<circle cx="480.5" cy="57.19" r="3.5" fill="currentColor"><title>severe: 1</title></circle>
<circle cx="497.5" cy="57.19" r="3.5" fill="currentColor"><title>end stage: 1</title></circle>
<text x="550" y="45" font-size="13">kidney_flow</text>
<rect x="550" y="125" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="584" y="125" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="618" y="125" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="635" y="125" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="652" y="125" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="669" y="125" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="550" x2="686" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="689" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="550,57.19 584,57.19 601,72.01 635,85.92 660.5,94.89 677.5,100.3 686,100.3" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="601" cy="72.01" r="3.5" fill="currentColor"><title>mild: 0.7712</title></circle>
<circle cx="635" cy="85.92" r="3.5" fill="currentColor"><title>moderate: 0.5567</title></circle>
<circle cx="660.5" cy="94.89" r="3.5" fill="currentColor"><title>severe: 0.4183</title></circle>
<circle cx="677.5" cy="100.3" r="3.5" fill="currentColor"><title>end stage: 0.3344</title></circle>
<text x="689" y="103.3" font-size="11" opacity="0.7">0.3344</text>
<text x="10" y="155" font-size="13">gastric_emptying_time</text>
<rect x="10" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="44" y="235" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="78" y="235" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="95" y="235" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="112" y="235" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="129" y="235" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="10" x2="146" y1="167.2" y2="167.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="170.2" font-size="11" opacity="0.7">1</text>
<polyline points="10,167.2 44,167.2 61,167.2 95,167.2 120.5,167.2 137.5,167.2 146,167.2" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="167.2" r="3.5" fill="currentColor"><title>mild: 1</title></circle>
<circle cx="95" cy="167.2" r="3.5" fill="currentColor"><title>moderate: 1</title></circle>
<circle cx="120.5" cy="167.2" r="3.5" fill="currentColor"><title>severe: 1</title></circle>
<circle cx="137.5" cy="167.2" r="3.5" fill="currentColor"><title>end stage: 1</title></circle>
<text x="190" y="155" font-size="13">si_transit_time</text>
<rect x="190" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="224" y="235" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="258" y="235" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="275" y="235" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="292" y="235" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="309" y="235" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="190" x2="326" y1="167.2" y2="167.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="329" y="170.2" font-size="11" opacity="0.7">1</text>
<polyline points="190,167.2 224,167.2 241,167.2 275,167.2 300.5,167.2 317.5,167.2 326,167.2" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="241" cy="167.2" r="3.5" fill="currentColor"><title>mild: 1</title></circle>
<circle cx="275" cy="167.2" r="3.5" fill="currentColor"><title>moderate: 1</title></circle>
<circle cx="300.5" cy="167.2" r="3.5" fill="currentColor"><title>severe: 1</title></circle>
<circle cx="317.5" cy="167.2" r="3.5" fill="currentColor"><title>end stage: 1</title></circle>
<text x="370" y="155" font-size="13">binding: albumin</text>
<rect x="370" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="404" y="235" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="438" y="235" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="455" y="235" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="472" y="235" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="489" y="235" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="370" x2="506" y1="167.2" y2="167.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="509" y="170.2" font-size="11" opacity="0.7">1</text>
<polyline points="370,167.2 404,167.2 421,169.4 455,171.7 480.5,172.9 497.5,176.5 506,176.5" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="421" cy="169.4" r="3.5" fill="currentColor"><title>mild: 0.9654</title></circle>
<circle cx="455" cy="171.7" r="3.5" fill="currentColor"><title>moderate: 0.9309</title></circle>
<circle cx="480.5" cy="172.9" r="3.5" fill="currentColor"><title>severe: 0.9124</title></circle>
<circle cx="497.5" cy="176.5" r="3.5" fill="currentColor"><title>end stage: 0.8557</title></circle>
<text x="550" y="155" font-size="13">binding: albumin_binding</text>
<rect x="550" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="584" y="235" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="618" y="235" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="635" y="235" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="652" y="235" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="669" y="235" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="550" x2="686" y1="167.2" y2="167.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="689" y="170.2" font-size="11" opacity="0.7">1</text>
<polyline points="550,167.2 584,167.2 601,167.2 635,171 660.5,177.6 677.5,188.6 686,188.6" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="601" cy="167.2" r="3.5" fill="currentColor"><title>mild: 1</title></circle>
<circle cx="635" cy="171" r="3.5" fill="currentColor"><title>moderate: 0.941</title></circle>
<circle cx="660.5" cy="177.6" r="3.5" fill="currentColor"><title>severe: 0.839</title></circle>
<circle cx="677.5" cy="188.6" r="3.5" fill="currentColor"><title>end stage: 0.669</title></circle>
<text x="689" y="191.6" font-size="11" opacity="0.7">0.669</text>
<text x="10" y="265" font-size="13">binding: agp</text>
<rect x="10" y="345" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="44" y="345" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="78" y="345" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="95" y="345" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="112" y="345" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="129" y="345" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="10" x2="146" y1="310.6" y2="310.6" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="313.6" font-size="11" opacity="0.7">1</text>
<polyline points="10,310.6 44,310.6 61,309.1 95,306.3 120.5,300.3 137.5,277.2 146,277.2" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="309.1" r="3.5" fill="currentColor"><title>mild: 1.046</title></circle>
<circle cx="95" cy="306.3" r="3.5" fill="currentColor"><title>moderate: 1.137</title></circle>
<circle cx="120.5" cy="300.3" r="3.5" fill="currentColor"><title>severe: 1.328</title></circle>
<circle cx="137.5" cy="277.2" r="3.5" fill="currentColor"><title>end stage: 2.062</title></circle>
<text x="149" y="280.2" font-size="11" opacity="0.7">2.062</text>
</svg>
</div>
<!-- end generated -->

<!-- generated: factors-heart-failure (scripts/pbpk_tables.py), do not edit -->
<div class="pbpk-figure pbpk-wide" markdown="0">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 730 144" role="img" aria-label="Heart failure: the factors over HF_severity (healthy, NYHA I to IV)" font-family="system-ui, sans-serif" font-size="14" fill="currentColor">
<title>Heart failure: the factors over HF_severity (healthy, NYHA I to IV)</title>
<text x="10" y="18" font-weight="600">Heart failure: the factors over HF_severity (healthy, NYHA I to IV)</text>
<text x="10" y="45" font-size="13">liver_flow</text>
<rect x="10" y="125" width="27.2" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="37.2" y="125" width="27.2" height="5" fill="#FDDBC7" fill-opacity="0.9"><title>NYHA I</title></rect>
<rect x="64.4" y="125" width="27.2" height="5" fill="#F4A582" fill-opacity="0.9"><title>NYHA II</title></rect>
<rect x="91.6" y="125" width="27.2" height="5" fill="#D6604D" fill-opacity="0.9"><title>NYHA III</title></rect>
<rect x="118.8" y="125" width="27.2" height="5" fill="#67001F" fill-opacity="0.9"><title>NYHA IV</title></rect>
<line x1="10" x2="146" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="23.6,57.19 50.8,57.19 78,72.42 105.2,87.13 132.4,91.99" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="50.8" cy="57.19" r="3.5" fill="currentColor"><title>I: 1</title></circle>
<circle cx="78" cy="72.42" r="3.5" fill="currentColor"><title>II: 0.765</title></circle>
<circle cx="105.2" cy="87.13" r="3.5" fill="currentColor"><title>III: 0.538</title></circle>
<circle cx="132.4" cy="91.99" r="3.5" fill="currentColor"><title>IV: 0.463</title></circle>
<text x="149" y="94.99" font-size="11" opacity="0.7">0.463</text>
<text x="190" y="45" font-size="13">kidney_flow</text>
<rect x="190" y="125" width="27.2" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="217.2" y="125" width="27.2" height="5" fill="#FDDBC7" fill-opacity="0.9"><title>NYHA I</title></rect>
<rect x="244.4" y="125" width="27.2" height="5" fill="#F4A582" fill-opacity="0.9"><title>NYHA II</title></rect>
<rect x="271.6" y="125" width="27.2" height="5" fill="#D6604D" fill-opacity="0.9"><title>NYHA III</title></rect>
<rect x="298.8" y="125" width="27.2" height="5" fill="#67001F" fill-opacity="0.9"><title>NYHA IV</title></rect>
<line x1="190" x2="326" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="329" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="203.6,57.19 230.8,57.19 258,71.57 285.2,86.35 312.4,81.23" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="230.8" cy="57.19" r="3.5" fill="currentColor"><title>I: 1</title></circle>
<circle cx="258" cy="71.57" r="3.5" fill="currentColor"><title>II: 0.778</title></circle>
<circle cx="285.2" cy="86.35" r="3.5" fill="currentColor"><title>III: 0.55</title></circle>
<circle cx="312.4" cy="81.23" r="3.5" fill="currentColor"><title>IV: 0.629</title></circle>
<text x="329" y="89.35" font-size="11" opacity="0.7">0.55</text>
<text x="370" y="45" font-size="13">peripheral_flow</text>
<rect x="370" y="125" width="27.2" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="397.2" y="125" width="27.2" height="5" fill="#FDDBC7" fill-opacity="0.9"><title>NYHA I</title></rect>
<rect x="424.4" y="125" width="27.2" height="5" fill="#F4A582" fill-opacity="0.9"><title>NYHA II</title></rect>
<rect x="451.6" y="125" width="27.2" height="5" fill="#D6604D" fill-opacity="0.9"><title>NYHA III</title></rect>
<rect x="478.8" y="125" width="27.2" height="5" fill="#67001F" fill-opacity="0.9"><title>NYHA IV</title></rect>
<line x1="370" x2="506" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="509" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="383.6,57.19 410.8,57.19 438,85.06 465.2,93.48 492.4,103.9" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="410.8" cy="57.19" r="3.5" fill="currentColor"><title>I: 1</title></circle>
<circle cx="438" cy="85.06" r="3.5" fill="currentColor"><title>II: 0.57</title></circle>
<circle cx="465.2" cy="93.48" r="3.5" fill="currentColor"><title>III: 0.44</title></circle>
<circle cx="492.4" cy="103.9" r="3.5" fill="currentColor"><title>IV: 0.28</title></circle>
<text x="509" y="106.9" font-size="11" opacity="0.7">0.28</text>
</svg>
</div>
<!-- end generated -->

## Scanning the classes in a model

An experiment derived from `PBPKSimulationExperiment` gets the changes of a class with `class_changes(classification, key)`, a simulation per class with `class_simulations(simulation, classification, prefix=...)` (the simulation of the healthy subjects with the changes of every class before the initialization, keys `<prefix>_<class key>`) and draws a selection of them with `add_classes(plot, classification, prefix=..., yid=...)`, a curve per class with its label and color:

```python
from pkpdsim.impairment import KDIGO_GFR_CLASSES
from pkpdsim.pbpk import PBPKSimulationExperiment

for c in KDIGO_GFR_CLASSES:
    changes = PBPKSimulationExperiment.class_changes(KDIGO_GFR_CLASSES, c.key)
    f = changes["f_renal_function"].m_as("dimensionless")
    print(f"{c.key:4} {c.label:20} f_renal_function {f:.4f}")
```

```python
# not executed
from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkpdsim.impairment import CLASSIFICATIONS


class ImpairmentScan(HctzSimulationExperiment):
    """A dose of 25 mg in every class of every condition."""

    def simulations(self) -> dict[str, Simulation]:
        dose = Simulation(end=48 * 60, preinit_changes={"PODOSE_hctz": Q(25, "mg")})
        sims = {}
        for key, classification in CLASSIFICATIONS.items():
            sims |= self.class_simulations(dose, classification, prefix=f"po_{key}")
        return sims

    def figures(self) -> dict[str, Figure]:
        fig = Figure(self, sid="Fig_classes", num_cols=len(CLASSIFICATIONS))
        plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)
        for plot, (key, classification) in zip(plots, CLASSIFICATIONS.items()):
            plot.set_title(classification.name)
            self.add_classes(plot, classification, prefix=f"po_{key}", yid="[Cve_hctz]")
        return self.figures_by_sid(fig)
```

A condition with a homeostasis of a PD component starts each class from its drug-free steady state, `self.equilibrated(simulation)` for each simulation of `class_simulations` ([Equilibration of a condition](impairment.md#equilibration-of-a-condition)). The colors of the classes are those of the maps of the experiments (`CIRRHOSIS_COLORS`, `HEART_FAILURE_COLORS`), so the figures of the classes and of the studies of a model match.
