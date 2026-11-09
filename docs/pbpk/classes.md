# Impairment classes

The models describe impaired organ function by the clinical classes the studies report: the Child-Pugh class of liver cirrhosis, the GFR category of chronic kidney disease (CKD) and the NYHA class of heart failure. This page explains what each condition does to the physiology and to the pharmacokinetics of a drug, defines the classes, shows which value of a model parameter represents each class, and how a model scans the conditions, by class or continuously and in combination. The references are at the end of the page. The factors which a class sets on the physiology are on the page [Impairment](impairment.md), the studies they come from in the [parameter database](../database.md).

Every class has one value of one parameter of the top model, set before the initialization like a covariate:

| condition | classification | the value of a class |
| --- | --- | --- |
| liver cirrhosis | Child-Pugh classes A, B, C ([Pugh 1973](#ref-Pugh1973)) | `CP_severity`: healthy 0, A 1, B 2, C 3; a MELD score maps onto it ([The MELD score](impairment.md#the-meld-score)) |
| chronic kidney disease | KDIGO GFR categories G1 to G5 ([KDIGO 2012](#ref-KDIGO2012)) | `f_renal_function`: the eGFR of the category over the eGFR of the build, G1 the renal function 1 |
| heart failure | NYHA functional classes I to IV ([NYHA 1994](#ref-NYHA1994)) | `HF_severity`: healthy 0, I 1, II 2, III 3, IV 4 |

A value between two classes interpolates the factors linearly, so a scan or a fit can use any severity. The classifications are `CHILD_PUGH_CLASSES`, `KDIGO_GFR_CLASSES` and `NYHA_CLASSES` of `pkpdsim.impairment` (`CLASSIFICATIONS` by key): every class has its key, a label, a color for every figure, its definition and the value which represents it.

A pictogram of the organ stands above every class, in the color of the class, its lost part pale (`impaired_organ` of `pkpdsim.pbpk.visual.icons`, the shapes of the anatomy resources): the liver by the portosystemic shunt of the class (`f_shunt`, the portal blood which bypasses the liver), a kidney by the eGFR of the class over that of the build, the heart by the cardiac output of the class relative to healthy (the sum of the flows which heart failure reduces). Hover a pictogram for its value.

<!-- generated: classes-figure (scripts/pbpk_tables.py), do not edit -->
<div class="pbpk-figure pbpk-wide" markdown="0">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 536" role="img" aria-label="The clinical classes of the impairments and their scan values" font-family="system-ui, sans-serif" font-size="14" fill="currentColor">
<title>The clinical classes of the impairments and their scan values</title>
<text x="10" y="26" font-weight="600">Child-Pugh classes by the score 5-15, the CP_severity of each below</text>
<g><title>Healthy: shunt 0.00, the portal blood which bypasses the liver</title><path d="M56.08,40.71 46.28,39.37 40.60,39.00 32.73,39.06 29.24,39.49 27.66,39.95 25.89,40.95 23.29,43.27 21.31,45.69 19.90,48.04 19.14,49.84 18.34,52.92 18.13,55.18 18.59,66.48 19.26,74.39 20.21,79.79 20.85,81.66 21.61,82.76 22.22,83.00 24.18,82.88 28.11,81.63 33.61,78.94 45.34,72.31 50.53,70.08 57.52,66.57 60.42,65.47 64.66,62.97 68.39,61.63 76.18,59.37 79.78,57.93 89.64,52.34 91.35,50.66 91.87,49.29 91.87,49.29 91.75,47.58 91.56,47.09 90.77,46.27 87.84,44.89 83.81,43.79 71.26,41.35 65.49,40.92Z" fill="#9e9e9e"/><path d="M56.08,40.71 46.28,39.37 40.60,39.00 32.73,39.06 29.24,39.49 27.66,39.95 25.89,40.95 23.29,43.27 21.31,45.69 19.90,48.04 19.14,49.84 18.34,52.92 18.13,55.18 18.59,66.48 19.26,74.39 20.21,79.79 20.85,81.66 21.61,82.76 22.22,83.00 24.18,82.88 28.11,81.63 33.61,78.94 45.34,72.31 50.53,70.08 57.52,66.57 60.42,65.47 64.66,62.97 68.39,61.63 76.18,59.37 79.78,57.93 89.64,52.34 91.35,50.66 91.87,49.29 91.75,47.58 91.56,47.09 90.77,46.27 87.84,44.89 83.81,43.79 71.26,41.35 65.49,40.92Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>Child-Pugh A: shunt 0.40, the portal blood which bypasses the liver</title><path d="M150.77,39.92 146.74,39.37 141.06,39.00 133.18,39.06 129.70,39.49 128.11,39.95 126.34,40.95 123.74,43.27 121.76,45.69 120.36,48.04 119.59,49.84 118.80,52.92 118.58,55.18 119.04,66.48 119.71,74.39 120.66,79.79 121.30,81.66 122.07,82.76 122.68,83.00 124.63,82.88 128.57,81.63 134.07,78.94 145.79,72.31 150.77,70.17Z" fill="#74a9cf"/><path d="M156.54,40.71 150.77,39.92 150.77,70.17 150.98,70.08 157.97,66.57 160.87,65.47 165.12,62.97 168.84,61.63 176.63,59.37 180.23,57.93 190.10,52.34 191.81,50.66 192.32,49.29 192.20,47.58 192.02,47.09 191.23,46.27 188.29,44.89 184.26,43.79 171.71,41.35 165.94,40.92Z" fill="#74a9cf" fill-opacity="0.15"/><path d="M156.54,40.71 146.74,39.37 141.06,39.00 133.18,39.06 129.70,39.49 128.11,39.95 126.34,40.95 123.74,43.27 121.76,45.69 120.36,48.04 119.59,49.84 118.80,52.92 118.58,55.18 119.04,66.48 119.71,74.39 120.66,79.79 121.30,81.66 122.07,82.76 122.68,83.00 124.63,82.88 128.57,81.63 134.07,78.94 145.79,72.31 150.98,70.08 157.97,66.57 160.87,65.47 165.12,62.97 168.84,61.63 176.63,59.37 180.23,57.93 190.10,52.34 191.81,50.66 192.32,49.29 192.20,47.58 192.02,47.09 191.23,46.27 188.29,44.89 184.26,43.79 171.71,41.35 165.94,40.92Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>Child-Pugh B: shunt 0.70, the portal blood which bypasses the liver</title><path d="M272.85,39.05 271.82,39.06 268.34,39.49 266.75,39.95 264.98,40.95 262.38,43.27 260.40,45.69 258.99,48.04 258.23,49.84 257.43,52.92 257.22,55.18 257.68,66.48 258.35,74.39 259.30,79.79 259.94,81.66 260.70,82.76 261.31,83.00 263.27,82.88 267.21,81.63 272.70,78.94 272.85,78.86Z" fill="#2b8cbe"/><path d="M295.17,40.71 285.37,39.37 279.69,39.00 272.85,39.05 272.85,78.86 284.43,72.31 289.62,70.08 296.61,66.57 299.51,65.47 303.76,62.97 307.48,61.63 315.27,59.37 318.87,57.93 328.73,52.34 330.44,50.66 330.96,49.29 330.84,47.58 330.66,47.09 329.86,46.27 326.93,44.89 322.90,43.79 310.35,41.35 304.58,40.92Z" fill="#2b8cbe" fill-opacity="0.15"/><path d="M295.17,40.71 285.37,39.37 279.69,39.00 271.82,39.06 268.34,39.49 266.75,39.95 264.98,40.95 262.38,43.27 260.40,45.69 258.99,48.04 258.23,49.84 257.43,52.92 257.22,55.18 257.68,66.48 258.35,74.39 259.30,79.79 259.94,81.66 260.70,82.76 261.31,83.00 263.27,82.88 267.21,81.63 272.70,78.94 284.43,72.31 289.62,70.08 296.61,66.57 299.51,65.47 303.76,62.97 307.48,61.63 315.27,59.37 318.87,57.93 328.73,52.34 330.44,50.66 330.96,49.29 330.84,47.58 330.66,47.09 329.86,46.27 326.93,44.89 322.90,43.79 310.35,41.35 304.58,40.92Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>Child-Pugh C: shunt 0.81, the portal blood which bypasses the liver</title><path d="M516.94,39.76 516.29,39.95 514.52,40.95 511.93,43.27 509.94,45.69 508.54,48.04 507.77,49.84 506.98,52.92 506.77,55.18 507.22,66.48 507.90,74.39 508.84,79.79 509.48,81.66 510.25,82.76 510.86,83.00 512.81,82.88 516.75,81.63 516.94,81.54Z" fill="#045a8d"/><path d="M544.72,40.71 534.92,39.37 529.24,39.00 521.36,39.06 517.88,39.49 516.94,39.76 516.94,81.54 522.25,78.94 533.97,72.31 539.16,70.08 546.16,66.57 549.06,65.47 553.30,62.97 557.03,61.63 564.81,59.37 568.41,57.93 578.28,52.34 579.99,50.66 580.51,49.29 580.38,47.58 580.20,47.09 579.41,46.27 576.48,44.89 572.45,43.79 559.90,41.35 554.12,40.92Z" fill="#045a8d" fill-opacity="0.15"/><path d="M544.72,40.71 534.92,39.37 529.24,39.00 521.36,39.06 517.88,39.49 516.29,39.95 514.52,40.95 511.93,43.27 509.94,45.69 508.54,48.04 507.77,49.84 506.98,52.92 506.77,55.18 507.22,66.48 507.90,74.39 508.84,79.79 509.48,81.66 510.25,82.76 510.86,83.00 512.81,82.88 516.75,81.63 522.25,78.94 533.97,72.31 539.16,70.08 546.16,66.57 549.06,65.47 553.30,62.97 557.03,61.63 564.81,59.37 568.41,57.93 578.28,52.34 579.99,50.66 580.51,49.29 580.38,47.58 580.20,47.09 579.41,46.27 576.48,44.89 572.45,43.79 559.90,41.35 554.12,40.92Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<rect x="11" y="94" width="88" height="22" rx="3" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy: no liver disease</title></rect>
<text x="55" y="135" text-anchor="middle">healthy</text>
<rect x="101" y="94" width="108.9" height="22" rx="3" fill="#74a9cf" fill-opacity="0.9"><title>Child-Pugh A: well-compensated cirrhosis, score 5-6</title></rect>
<text x="155.5" y="135" text-anchor="middle">A</text>
<rect x="211.9" y="94" width="164.4" height="22" rx="3" fill="#2b8cbe" fill-opacity="0.9"><title>Child-Pugh B: significant functional compromise, score 7-9</title></rect>
<text x="294.1" y="135" text-anchor="middle">B</text>
<rect x="378.3" y="94" width="330.7" height="22" rx="3" fill="#045a8d" fill-opacity="0.9"><title>Child-Pugh C: decompensated cirrhosis, score 10-15</title></rect>
<text x="543.6" y="135" text-anchor="middle">C</text>
<line x1="55" x2="55" y1="91" y2="119" stroke="currentColor" stroke-width="2"/>
<text x="55" y="152" text-anchor="middle" font-size="12" opacity="0.8">0</text>
<line x1="155.5" x2="155.5" y1="91" y2="119" stroke="currentColor" stroke-width="2"/>
<text x="155.5" y="152" text-anchor="middle" font-size="12" opacity="0.8">1</text>
<line x1="294.1" x2="294.1" y1="91" y2="119" stroke="currentColor" stroke-width="2"/>
<text x="294.1" y="152" text-anchor="middle" font-size="12" opacity="0.8">2</text>
<line x1="543.6" x2="543.6" y1="91" y2="119" stroke="currentColor" stroke-width="2"/>
<text x="543.6" y="152" text-anchor="middle" font-size="12" opacity="0.8">3</text>
<text x="10" y="174" font-size="12" opacity="0.8">MELD of the classes (consensus of the database): healthy 6.43, A 7.91, B 12.8, C 20.5</text>
<text x="10" y="202" font-weight="600">KDIGO GFR categories by the eGFR (120 to 0 ml/min/1.73 m2), the eGFR of each below</text>
<g><title>G1 normal: eGFR 114.5, 100 % of the renal function of the build</title><path d="M94.14,215.00 91.11,216.31 89.56,217.51 87.84,219.85 87.12,221.89 86.67,224.09 86.71,227.15 87.36,230.22 87.98,231.46 88.70,232.25 91.11,233.73 92.35,235.86 93.11,238.27 93.42,240.62 93.28,243.33 91.66,246.81 91.15,248.88 91.08,251.18 91.63,253.59 92.87,255.87 94.14,257.38 94.97,258.00 97.04,258.90 99.69,259.00 101.72,258.31 104.16,256.49 105.95,254.01 107.16,251.01 108.05,246.71 108.33,242.51 108.02,237.03 107.47,231.84 106.61,227.50 105.64,224.33 104.33,221.40 102.58,218.86 100.34,216.86 98.55,215.90 95.38,215.07 94.14,215.00Z" fill="#9e9e9e"/><path d="M94.14,215.00 91.11,216.31 89.56,217.51 87.84,219.85 87.12,221.89 86.67,224.09 86.71,227.15 87.36,230.22 87.98,231.46 88.70,232.25 91.11,233.73 92.35,235.86 93.11,238.27 93.42,240.62 93.28,243.33 91.66,246.81 91.15,248.88 91.08,251.18 91.63,253.59 92.87,255.87 94.14,257.38 94.97,258.00 97.04,258.90 99.69,259.00 101.72,258.31 104.16,256.49 105.95,254.01 107.16,251.01 108.05,246.71 108.33,242.51 108.02,237.03 107.47,231.84 106.61,227.50 105.64,224.33 104.33,221.40 102.58,218.86 100.34,216.86 98.55,215.90 95.38,215.07Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>G2 mild: eGFR 75, 66 % of the renal function of the build</title><path d="M262.24,229.66 262.36,230.22 262.98,231.46 263.70,232.25 266.11,233.73 267.35,235.86 268.11,238.27 268.42,240.62 268.28,243.33 266.66,246.81 266.15,248.88 266.08,251.18 266.63,253.59 267.87,255.87 269.14,257.38 269.97,258.00 272.04,258.90 274.69,259.00 276.72,258.31 279.16,256.49 280.95,254.01 282.16,251.01 283.05,246.71 283.33,242.51 283.02,237.03 282.47,231.84 282.03,229.66Z" fill="#66c2a4"/><path d="M269.14,215.00 266.11,216.31 264.56,217.51 262.84,219.85 262.12,221.89 261.67,224.09 261.71,227.15 262.24,229.66 282.03,229.66 281.61,227.50 280.64,224.33 279.33,221.40 277.58,218.86 275.34,216.86 273.55,215.90 270.38,215.07Z" fill="#66c2a4" fill-opacity="0.15"/><path d="M269.14,215.00 266.11,216.31 264.56,217.51 262.84,219.85 262.12,221.89 261.67,224.09 261.71,227.15 262.36,230.22 262.98,231.46 263.70,232.25 266.11,233.73 267.35,235.86 268.11,238.27 268.42,240.62 268.28,243.33 266.66,246.81 266.15,248.88 266.08,251.18 266.63,253.59 267.87,255.87 269.14,257.38 269.97,258.00 272.04,258.90 274.69,259.00 276.72,258.31 279.16,256.49 280.95,254.01 282.16,251.01 283.05,246.71 283.33,242.51 283.02,237.03 282.47,231.84 281.61,227.50 280.64,224.33 279.33,221.40 277.58,218.86 275.34,216.86 273.55,215.90 270.38,215.07Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>G3a mild-moderate: eGFR 52, 45 % of the renal function of the build</title><path d="M399.10,237.45 399.36,238.27 399.67,240.62 399.53,243.33 397.91,246.81 397.40,248.88 397.33,251.18 397.88,253.59 399.12,255.87 400.39,257.38 401.22,258.00 403.29,258.90 405.94,259.00 407.97,258.31 410.41,256.49 412.20,254.01 413.41,251.01 414.30,246.71 414.58,242.51 414.29,237.45Z" fill="#41ae76"/><path d="M400.39,215.00 397.36,216.31 395.81,217.51 394.09,219.85 393.37,221.89 392.92,224.09 392.96,227.15 393.61,230.22 394.23,231.46 394.95,232.25 397.36,233.73 398.60,235.86 399.10,237.45 414.29,237.45 414.27,237.03 413.72,231.84 412.86,227.50 411.89,224.33 410.58,221.40 408.83,218.86 406.59,216.86 404.80,215.90 401.63,215.07Z" fill="#41ae76" fill-opacity="0.15"/><path d="M400.39,215.00 397.36,216.31 395.81,217.51 394.09,219.85 393.37,221.89 392.92,224.09 392.96,227.15 393.61,230.22 394.23,231.46 394.95,232.25 397.36,233.73 398.60,235.86 399.36,238.27 399.67,240.62 399.53,243.33 397.91,246.81 397.40,248.88 397.33,251.18 397.88,253.59 399.12,255.87 400.39,257.38 401.22,258.00 403.29,258.90 405.94,259.00 407.97,258.31 410.41,256.49 412.20,254.01 413.41,251.01 414.30,246.71 414.58,242.51 414.27,237.03 413.72,231.84 412.86,227.50 411.89,224.33 410.58,221.40 408.83,218.86 406.59,216.86 404.80,215.90 401.63,215.07Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>G3b moderate-severe: eGFR 37, 32 % of the renal function of the build</title><path d="M487.02,243.36 485.41,246.81 484.90,248.88 484.83,251.18 485.38,253.59 486.62,255.87 487.89,257.38 488.72,258.00 490.79,258.90 493.44,259.00 495.47,258.31 497.91,256.49 499.70,254.01 500.91,251.01 501.80,246.71 502.02,243.36Z" fill="#238b45"/><path d="M487.89,215.00 484.86,216.31 483.31,217.51 481.59,219.85 480.87,221.89 480.42,224.09 480.46,227.15 481.11,230.22 481.73,231.46 482.45,232.25 484.86,233.73 486.10,235.86 486.86,238.27 487.17,240.62 487.03,243.33 487.02,243.36 502.02,243.36 502.08,242.51 501.77,237.03 501.22,231.84 500.36,227.50 499.39,224.33 498.08,221.40 496.33,218.86 494.09,216.86 492.30,215.90 489.13,215.07Z" fill="#238b45" fill-opacity="0.15"/><path d="M487.89,215.00 484.86,216.31 483.31,217.51 481.59,219.85 480.87,221.89 480.42,224.09 480.46,227.15 481.11,230.22 481.73,231.46 482.45,232.25 484.86,233.73 486.10,235.86 486.86,238.27 487.17,240.62 487.03,243.33 485.41,246.81 484.90,248.88 484.83,251.18 485.38,253.59 486.62,255.87 487.89,257.38 488.72,258.00 490.79,258.90 493.44,259.00 495.47,258.31 497.91,256.49 499.70,254.01 500.91,251.01 501.80,246.71 502.08,242.51 501.77,237.03 501.22,231.84 500.36,227.50 499.39,224.33 498.08,221.40 496.33,218.86 494.09,216.86 492.30,215.90 489.13,215.07Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>G4 severe: eGFR 22.5, 20 % of the renal function of the build</title><path d="M572.44,248.71 572.40,248.88 572.33,251.18 572.88,253.59 574.12,255.87 575.39,257.38 576.22,258.00 578.29,258.90 580.94,259.00 582.97,258.31 585.41,256.49 587.20,254.01 588.41,251.01 588.89,248.71Z" fill="#006d2c"/><path d="M575.39,215.00 572.36,216.31 570.81,217.51 569.09,219.85 568.37,221.89 567.92,224.09 567.96,227.15 568.61,230.22 569.23,231.46 569.95,232.25 572.36,233.73 573.60,235.86 574.36,238.27 574.67,240.62 574.53,243.33 572.91,246.81 572.44,248.71 588.89,248.71 589.30,246.71 589.58,242.51 589.27,237.03 588.72,231.84 587.86,227.50 586.89,224.33 585.58,221.40 583.83,218.86 581.59,216.86 579.80,215.90 576.63,215.07Z" fill="#006d2c" fill-opacity="0.15"/><path d="M575.39,215.00 572.36,216.31 570.81,217.51 569.09,219.85 568.37,221.89 567.92,224.09 567.96,227.15 568.61,230.22 569.23,231.46 569.95,232.25 572.36,233.73 573.60,235.86 574.36,238.27 574.67,240.62 574.53,243.33 572.91,246.81 572.40,248.88 572.33,251.18 572.88,253.59 574.12,255.87 575.39,257.38 576.22,258.00 578.29,258.90 580.94,259.00 582.97,258.31 585.41,256.49 587.20,254.01 588.41,251.01 589.30,246.71 589.58,242.51 589.27,237.03 588.72,231.84 587.86,227.50 586.89,224.33 585.58,221.40 583.83,218.86 581.59,216.86 579.80,215.90 576.63,215.07Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>G5 kidney failure: eGFR 7.5, 7 % of the renal function of the build</title><path d="M660.82,254.40 661.62,255.87 662.89,257.38 663.72,258.00 665.79,258.90 668.44,259.00 670.47,258.31 672.91,256.49 674.42,254.40Z" fill="#00441b"/><path d="M662.89,215.00 659.86,216.31 658.31,217.51 656.59,219.85 655.87,221.89 655.42,224.09 655.46,227.15 656.11,230.22 656.73,231.46 657.45,232.25 659.86,233.73 661.10,235.86 661.86,238.27 662.17,240.62 662.03,243.33 660.41,246.81 659.90,248.88 659.83,251.18 660.38,253.59 660.82,254.40 674.42,254.40 674.70,254.01 675.91,251.01 676.80,246.71 677.08,242.51 676.77,237.03 676.22,231.84 675.36,227.50 674.39,224.33 673.08,221.40 671.33,218.86 669.09,216.86 667.30,215.90 664.13,215.07Z" fill="#00441b" fill-opacity="0.15"/><path d="M662.89,215.00 659.86,216.31 658.31,217.51 656.59,219.85 655.87,221.89 655.42,224.09 655.46,227.15 656.11,230.22 656.73,231.46 657.45,232.25 659.86,233.73 661.10,235.86 661.86,238.27 662.17,240.62 662.03,243.33 660.41,246.81 659.90,248.88 659.83,251.18 660.38,253.59 661.62,255.87 662.89,257.38 663.72,258.00 665.79,258.90 668.44,259.00 670.47,258.31 672.91,256.49 674.70,254.01 675.91,251.01 676.80,246.71 677.08,242.51 676.77,237.03 676.22,231.84 675.36,227.50 674.39,224.33 673.08,221.40 671.33,218.86 669.09,216.86 667.30,215.90 664.13,215.07Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<rect x="11" y="270" width="173" height="22" rx="3" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal: normal or high, eGFR &gt;= 90</title></rect>
<text x="97.5" y="311" text-anchor="middle">G1</text>
<rect x="186" y="270" width="173" height="22" rx="3" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild: mildly decreased, eGFR 60-89</title></rect>
<text x="272.5" y="311" text-anchor="middle">G2</text>
<rect x="361" y="270" width="85.5" height="22" rx="3" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate: mildly to moderately decreased, eGFR 45-59</title></rect>
<text x="403.8" y="311" text-anchor="middle">G3a</text>
<rect x="448.5" y="270" width="85.5" height="22" rx="3" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe: moderately to severely decreased, eGFR 30-44</title></rect>
<text x="491.2" y="311" text-anchor="middle">G3b</text>
<rect x="536" y="270" width="85.5" height="22" rx="3" fill="#006d2c" fill-opacity="0.9"><title>G4 severe: severely decreased, eGFR 15-29</title></rect>
<text x="578.8" y="311" text-anchor="middle">G4</text>
<rect x="623.5" y="270" width="85.5" height="22" rx="3" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure: kidney failure, eGFR &lt; 15 (not on dialysis in the library)</title></rect>
<text x="666.2" y="311" text-anchor="middle">G5</text>
<line x1="272.5" x2="272.5" y1="267" y2="295" stroke="currentColor" stroke-width="2"/>
<text x="272.5" y="328" text-anchor="middle" font-size="12" opacity="0.8">75</text>
<line x1="406.7" x2="406.7" y1="267" y2="295" stroke="currentColor" stroke-width="2"/>
<text x="406.7" y="328" text-anchor="middle" font-size="12" opacity="0.8">52</text>
<line x1="494.2" x2="494.2" y1="267" y2="295" stroke="currentColor" stroke-width="2"/>
<text x="494.2" y="328" text-anchor="middle" font-size="12" opacity="0.8">37</text>
<line x1="578.8" x2="578.8" y1="267" y2="295" stroke="currentColor" stroke-width="2"/>
<text x="578.8" y="328" text-anchor="middle" font-size="12" opacity="0.8">22.5</text>
<line x1="666.2" x2="666.2" y1="267" y2="295" stroke="currentColor" stroke-width="2"/>
<text x="666.2" y="328" text-anchor="middle" font-size="12" opacity="0.8">7.5</text>
<text x="10" y="350" font-size="12" opacity="0.8">the eGFR of the scan of a category; G1 is the renal function 1 of the build (eGFR 114.5 male, 114.7 female)</text>
<text x="10" y="378" font-weight="600">NYHA functional classes, the HF_severity of each below</text>
<g><title>Healthy: cardiac output 1.00 of healthy (the sum of its flows)</title><path d="M72.90,393.39 71.31,391.80 70.29,392.06 69.76,393.01 70.82,394.90 69.27,396.72 68.89,399.07 66.69,399.03 65.71,399.48 64.65,400.39 64.27,401.07 64.65,402.25 62.41,402.62 62.07,403.42 62.87,403.87 65.52,404.10 65.36,404.86 63.21,405.28 63.36,406.45 64.42,406.79 66.08,406.49 66.65,414.70 67.49,418.94 67.56,423.34 67.30,427.20 67.90,428.15 68.85,429.06 69.91,429.55 72.29,429.70 72.71,429.21 73.47,424.97 76.50,427.43 79.87,428.71 80.28,433.98 80.93,434.51 82.10,434.96 83.77,435.00 85.24,434.55 86.04,434.02 86.49,429.89 90.96,429.93 92.36,429.36 95.20,427.01 97.13,423.91 97.85,421.10 97.93,419.51 97.21,414.93 96.26,412.05 94.48,408.04 91.19,402.55 95.09,402.40 95.70,401.75 94.94,401.15 92.78,400.69 93.95,400.47 94.37,399.56 93.80,399.03 91.53,400.24 86.91,400.05 91.95,398.19 91.98,397.10 91.53,395.81 90.96,395.24 83.73,396.79 81.61,394.71 78.81,393.23 79.34,391.98 78.31,391.42 77.60,391.87 77.03,392.86 76.01,393.01 76.46,391.95 75.97,391.00 75.97,391.00 75.21,391.38 74.19,393.12Z" fill="#9e9e9e"/><path d="M72.90,393.39 71.31,391.80 70.29,392.06 69.76,393.01 70.82,394.90 69.27,396.72 68.89,399.07 66.69,399.03 65.71,399.48 64.65,400.39 64.27,401.07 64.65,402.25 62.41,402.62 62.07,403.42 62.87,403.87 65.52,404.10 65.36,404.86 63.21,405.28 63.36,406.45 64.42,406.79 66.08,406.49 66.65,414.70 67.49,418.94 67.56,423.34 67.30,427.20 67.90,428.15 68.85,429.06 69.91,429.55 72.29,429.70 72.71,429.21 73.47,424.97 76.50,427.43 79.87,428.71 80.28,433.98 80.93,434.51 82.10,434.96 83.77,435.00 85.24,434.55 86.04,434.02 86.49,429.89 90.96,429.93 92.36,429.36 95.20,427.01 97.13,423.91 97.85,421.10 97.93,419.51 97.21,414.93 96.26,412.05 94.48,408.04 91.19,402.55 95.09,402.40 95.70,401.75 94.94,401.15 92.78,400.69 93.95,400.47 94.37,399.56 93.80,399.03 91.53,400.24 86.91,400.05 91.95,398.19 91.98,397.10 91.53,395.81 90.96,395.24 83.73,396.79 81.61,394.71 78.81,393.23 79.34,391.98 78.31,391.42 77.60,391.87 77.03,392.86 76.01,393.01 76.46,391.95 75.97,391.00 75.21,391.38 74.19,393.12Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>NYHA I: cardiac output 0.89 of healthy (the sum of its flows)</title><path d="M206.17,399.27 205.71,399.48 204.65,400.39 204.27,401.07 204.65,402.25 202.41,402.62 202.07,403.42 202.87,403.87 205.52,404.10 205.36,404.86 203.21,405.28 203.36,406.45 204.42,406.79 206.08,406.49 206.65,414.70 207.49,418.94 207.56,423.34 207.30,427.20 207.90,428.15 208.85,429.06 209.91,429.55 212.29,429.70 212.71,429.21 213.47,424.97 216.50,427.43 219.87,428.71 220.28,433.98 220.93,434.51 222.10,434.96 223.77,435.00 225.24,434.55 226.04,434.02 226.49,429.89 230.96,429.93 232.36,429.36 235.20,427.01 237.13,423.91 237.85,421.10 237.93,419.51 237.21,414.93 236.26,412.05 234.48,408.04 231.19,402.55 235.09,402.40 235.70,401.75 234.94,401.15 232.78,400.69 233.95,400.47 234.37,399.56 234.06,399.27 233.36,399.27 231.53,400.24 226.91,400.05 229.04,399.27Z" fill="#FDDBC7"/><path d="M212.90,393.39 211.31,391.80 210.29,392.06 209.76,393.01 210.82,394.90 209.27,396.72 208.89,399.07 206.69,399.03 206.17,399.27 234.06,399.27 233.80,399.03 233.36,399.27 229.04,399.27 231.95,398.19 231.98,397.10 231.53,395.81 230.96,395.24 223.73,396.79 221.61,394.71 218.81,393.23 219.34,391.98 218.31,391.42 217.60,391.87 217.03,392.86 216.01,393.01 216.46,391.95 215.97,391.00 215.21,391.38 214.19,393.12Z" fill="#FDDBC7" fill-opacity="0.15"/><path d="M212.90,393.39 211.31,391.80 210.29,392.06 209.76,393.01 210.82,394.90 209.27,396.72 208.89,399.07 206.69,399.03 205.71,399.48 204.65,400.39 204.27,401.07 204.65,402.25 202.41,402.62 202.07,403.42 202.87,403.87 205.52,404.10 205.36,404.86 203.21,405.28 203.36,406.45 204.42,406.79 206.08,406.49 206.65,414.70 207.49,418.94 207.56,423.34 207.30,427.20 207.90,428.15 208.85,429.06 209.91,429.55 212.29,429.70 212.71,429.21 213.47,424.97 216.50,427.43 219.87,428.71 220.28,433.98 220.93,434.51 222.10,434.96 223.77,435.00 225.24,434.55 226.04,434.02 226.49,429.89 230.96,429.93 232.36,429.36 235.20,427.01 237.13,423.91 237.85,421.10 237.93,419.51 237.21,414.93 236.26,412.05 234.48,408.04 231.19,402.55 235.09,402.40 235.70,401.75 234.94,401.15 232.78,400.69 233.95,400.47 234.37,399.56 233.80,399.03 231.53,400.24 226.91,400.05 231.95,398.19 231.98,397.10 231.53,395.81 230.96,395.24 223.73,396.79 221.61,394.71 218.81,393.23 219.34,391.98 218.31,391.42 217.60,391.87 217.03,392.86 216.01,393.01 216.46,391.95 215.97,391.00 215.21,391.38 214.19,393.12Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>NYHA II: cardiac output 0.78 of healthy (the sum of its flows)</title><path d="M342.17,403.18 342.07,403.42 342.87,403.87 345.52,404.10 345.36,404.86 343.21,405.28 343.36,406.45 344.42,406.79 346.08,406.49 346.65,414.70 347.49,418.94 347.56,423.34 347.30,427.20 347.90,428.15 348.85,429.06 349.91,429.55 352.29,429.70 352.71,429.21 353.47,424.97 356.50,427.43 359.87,428.71 360.28,433.98 360.93,434.51 362.10,434.96 363.77,435.00 365.24,434.55 366.04,434.02 366.49,429.89 370.96,429.93 372.36,429.36 375.20,427.01 377.13,423.91 377.85,421.10 377.93,419.51 377.21,414.93 376.26,412.05 374.48,408.04 371.57,403.18Z" fill="#F4A582"/><path d="M352.90,393.39 351.31,391.80 350.29,392.06 349.76,393.01 350.82,394.90 349.27,396.72 348.89,399.07 346.69,399.03 345.71,399.48 344.65,400.39 344.27,401.07 344.65,402.25 342.41,402.62 342.17,403.18 371.57,403.18 371.19,402.55 375.09,402.40 375.70,401.75 374.94,401.15 372.78,400.69 373.95,400.47 374.37,399.56 373.80,399.03 371.53,400.24 366.91,400.05 371.95,398.19 371.98,397.10 371.53,395.81 370.96,395.24 363.73,396.79 361.61,394.71 358.81,393.23 359.34,391.98 358.31,391.42 357.60,391.87 357.03,392.86 356.01,393.01 356.46,391.95 355.97,391.00 355.21,391.38 354.19,393.12Z" fill="#F4A582" fill-opacity="0.15"/><path d="M352.90,393.39 351.31,391.80 350.29,392.06 349.76,393.01 350.82,394.90 349.27,396.72 348.89,399.07 346.69,399.03 345.71,399.48 344.65,400.39 344.27,401.07 344.65,402.25 342.41,402.62 342.07,403.42 342.87,403.87 345.52,404.10 345.36,404.86 343.21,405.28 343.36,406.45 344.42,406.79 346.08,406.49 346.65,414.70 347.49,418.94 347.56,423.34 347.30,427.20 347.90,428.15 348.85,429.06 349.91,429.55 352.29,429.70 352.71,429.21 353.47,424.97 356.50,427.43 359.87,428.71 360.28,433.98 360.93,434.51 362.10,434.96 363.77,435.00 365.24,434.55 366.04,434.02 366.49,429.89 370.96,429.93 372.36,429.36 375.20,427.01 377.13,423.91 377.85,421.10 377.93,419.51 377.21,414.93 376.26,412.05 374.48,408.04 371.19,402.55 375.09,402.40 375.70,401.75 374.94,401.15 372.78,400.69 373.95,400.47 374.37,399.56 373.80,399.03 371.53,400.24 366.91,400.05 371.95,398.19 371.98,397.10 371.53,395.81 370.96,395.24 363.73,396.79 361.61,394.71 358.81,393.23 359.34,391.98 358.31,391.42 357.60,391.87 357.03,392.86 356.01,393.01 356.46,391.95 355.97,391.00 355.21,391.38 354.19,393.12Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>NYHA III: cardiac output 0.65 of healthy (the sum of its flows)</title><path d="M486.16,407.62 486.65,414.70 487.49,418.94 487.56,423.34 487.30,427.20 487.90,428.15 488.85,429.06 489.91,429.55 492.29,429.70 492.71,429.21 493.47,424.97 496.50,427.43 499.87,428.71 500.28,433.98 500.93,434.51 502.10,434.96 503.77,435.00 505.24,434.55 506.04,434.02 506.49,429.89 510.96,429.93 512.36,429.36 515.20,427.01 517.13,423.91 517.85,421.10 517.93,419.51 517.21,414.93 516.26,412.05 514.48,408.04 514.23,407.62Z" fill="#D6604D"/><path d="M492.90,393.39 491.31,391.80 490.29,392.06 489.76,393.01 490.82,394.90 489.27,396.72 488.89,399.07 486.69,399.03 485.71,399.48 484.65,400.39 484.27,401.07 484.65,402.25 482.41,402.62 482.07,403.42 482.87,403.87 485.52,404.10 485.36,404.86 483.21,405.28 483.36,406.45 484.42,406.79 486.08,406.49 486.16,407.62 514.23,407.62 511.19,402.55 515.09,402.40 515.70,401.75 514.94,401.15 512.78,400.69 513.95,400.47 514.37,399.56 513.80,399.03 511.53,400.24 506.91,400.05 511.95,398.19 511.98,397.10 511.53,395.81 510.96,395.24 503.73,396.79 501.61,394.71 498.81,393.23 499.34,391.98 498.31,391.42 497.60,391.87 497.03,392.86 496.01,393.01 496.46,391.95 495.97,391.00 495.21,391.38 494.19,393.12Z" fill="#D6604D" fill-opacity="0.15"/><path d="M492.90,393.39 491.31,391.80 490.29,392.06 489.76,393.01 490.82,394.90 489.27,396.72 488.89,399.07 486.69,399.03 485.71,399.48 484.65,400.39 484.27,401.07 484.65,402.25 482.41,402.62 482.07,403.42 482.87,403.87 485.52,404.10 485.36,404.86 483.21,405.28 483.36,406.45 484.42,406.79 486.08,406.49 486.65,414.70 487.49,418.94 487.56,423.34 487.30,427.20 487.90,428.15 488.85,429.06 489.91,429.55 492.29,429.70 492.71,429.21 493.47,424.97 496.50,427.43 499.87,428.71 500.28,433.98 500.93,434.51 502.10,434.96 503.77,435.00 505.24,434.55 506.04,434.02 506.49,429.89 510.96,429.93 512.36,429.36 515.20,427.01 517.13,423.91 517.85,421.10 517.93,419.51 517.21,414.93 516.26,412.05 514.48,408.04 511.19,402.55 515.09,402.40 515.70,401.75 514.94,401.15 512.78,400.69 513.95,400.47 514.37,399.56 513.80,399.03 511.53,400.24 506.91,400.05 511.95,398.19 511.98,397.10 511.53,395.81 510.96,395.24 503.73,396.79 501.61,394.71 498.81,393.23 499.34,391.98 498.31,391.42 497.60,391.87 497.03,392.86 496.01,393.01 496.46,391.95 495.97,391.00 495.21,391.38 494.19,393.12Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<g><title>NYHA IV: cardiac output 0.62 of healthy (the sum of its flows)</title><path d="M626.24,408.72 626.65,414.70 627.49,418.94 627.56,423.34 627.30,427.20 627.90,428.15 628.85,429.06 629.91,429.55 632.29,429.70 632.71,429.21 633.47,424.97 636.50,427.43 639.87,428.71 640.28,433.98 640.93,434.51 642.10,434.96 643.77,435.00 645.24,434.55 646.04,434.02 646.49,429.89 650.96,429.93 652.36,429.36 655.20,427.01 657.13,423.91 657.85,421.10 657.93,419.51 657.21,414.93 656.26,412.05 654.78,408.72Z" fill="#67001F"/><path d="M632.90,393.39 631.31,391.80 630.29,392.06 629.76,393.01 630.82,394.90 629.27,396.72 628.89,399.07 626.69,399.03 625.71,399.48 624.65,400.39 624.27,401.07 624.65,402.25 622.41,402.62 622.07,403.42 622.87,403.87 625.52,404.10 625.36,404.86 623.21,405.28 623.36,406.45 624.42,406.79 626.08,406.49 626.24,408.72 654.78,408.72 654.48,408.04 651.19,402.55 655.09,402.40 655.70,401.75 654.94,401.15 652.78,400.69 653.95,400.47 654.37,399.56 653.80,399.03 651.53,400.24 646.91,400.05 651.95,398.19 651.98,397.10 651.53,395.81 650.96,395.24 643.73,396.79 641.61,394.71 638.81,393.23 639.34,391.98 638.31,391.42 637.60,391.87 637.03,392.86 636.01,393.01 636.46,391.95 635.97,391.00 635.21,391.38 634.19,393.12Z" fill="#67001F" fill-opacity="0.15"/><path d="M632.90,393.39 631.31,391.80 630.29,392.06 629.76,393.01 630.82,394.90 629.27,396.72 628.89,399.07 626.69,399.03 625.71,399.48 624.65,400.39 624.27,401.07 624.65,402.25 622.41,402.62 622.07,403.42 622.87,403.87 625.52,404.10 625.36,404.86 623.21,405.28 623.36,406.45 624.42,406.79 626.08,406.49 626.65,414.70 627.49,418.94 627.56,423.34 627.30,427.20 627.90,428.15 628.85,429.06 629.91,429.55 632.29,429.70 632.71,429.21 633.47,424.97 636.50,427.43 639.87,428.71 640.28,433.98 640.93,434.51 642.10,434.96 643.77,435.00 645.24,434.55 646.04,434.02 646.49,429.89 650.96,429.93 652.36,429.36 655.20,427.01 657.13,423.91 657.85,421.10 657.93,419.51 657.21,414.93 656.26,412.05 654.48,408.04 651.19,402.55 655.09,402.40 655.70,401.75 654.94,401.15 652.78,400.69 653.95,400.47 654.37,399.56 653.80,399.03 651.53,400.24 646.91,400.05 651.95,398.19 651.98,397.10 651.53,395.81 650.96,395.24 643.73,396.79 641.61,394.71 638.81,393.23 639.34,391.98 638.31,391.42 637.60,391.87 637.03,392.86 636.01,393.01 636.46,391.95 635.97,391.00 635.21,391.38 634.19,393.12Z" fill="none" stroke="currentColor" stroke-opacity="0.7" stroke-width="1"/></g>
<rect x="11" y="446" width="138" height="22" rx="3" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy: no heart failure</title></rect>
<text x="80" y="487" text-anchor="middle">healthy</text>
<rect x="151" y="446" width="138" height="22" rx="3" fill="#FDDBC7" fill-opacity="0.9"><title>NYHA I: no limitation of physical activity; ordinary activity causes no undue fatigue, palpitation or dyspnoea</title></rect>
<text x="220" y="487" text-anchor="middle">I</text>
<rect x="291" y="446" width="138" height="22" rx="3" fill="#F4A582" fill-opacity="0.9"><title>NYHA II: slight limitation; comfortable at rest, ordinary activity causes fatigue, palpitation or dyspnoea</title></rect>
<text x="360" y="487" text-anchor="middle">II</text>
<rect x="431" y="446" width="138" height="22" rx="3" fill="#D6604D" fill-opacity="0.9"><title>NYHA III: marked limitation; comfortable at rest, less than ordinary activity causes symptoms</title></rect>
<text x="500" y="487" text-anchor="middle">III</text>
<rect x="571" y="446" width="138" height="22" rx="3" fill="#67001F" fill-opacity="0.9"><title>NYHA IV: unable to carry on any physical activity without discomfort; symptoms at rest</title></rect>
<text x="640" y="487" text-anchor="middle">IV</text>
<line x1="80" x2="80" y1="443" y2="471" stroke="currentColor" stroke-width="2"/>
<text x="80" y="504" text-anchor="middle" font-size="12" opacity="0.8">0</text>
<line x1="220" x2="220" y1="443" y2="471" stroke="currentColor" stroke-width="2"/>
<text x="220" y="504" text-anchor="middle" font-size="12" opacity="0.8">1</text>
<line x1="360" x2="360" y1="443" y2="471" stroke="currentColor" stroke-width="2"/>
<text x="360" y="504" text-anchor="middle" font-size="12" opacity="0.8">2</text>
<line x1="500" x2="500" y1="443" y2="471" stroke="currentColor" stroke-width="2"/>
<text x="500" y="504" text-anchor="middle" font-size="12" opacity="0.8">3</text>
<line x1="640" x2="640" y1="443" y2="471" stroke="currentColor" stroke-width="2"/>
<text x="640" y="504" text-anchor="middle" font-size="12" opacity="0.8">4</text>
<text x="10" y="526" font-size="12" opacity="0.8">the classes grade the symptoms at physical activity</text>
</svg>
</div>
<!-- end generated -->

## Liver cirrhosis

### Background

Cirrhosis is the end stage of chronic liver injury (alcohol, viral hepatitis, metabolic and autoimmune disease): fibrous septa and regenerative nodules replace the normal architecture of the liver. Three changes matter for a drug ([Verbeeck 2008](#ref-Verbeeck2008)):

- **Loss of functional hepatocytes.** The mass of working liver cells falls, and with it the amount of the enzymes and transporters; the decline differs between the proteins (CYP2C19 and CYP1A2 more than CYP2D6, glucuronidation relatively spared), which the factors of the proteins of the [parameter database](../database.md#liver-cirrhosis) take per protein.
- **Portosystemic shunting.** The fibrosis raises the resistance of the liver to the portal blood, the portal pressure rises (portal hypertension) and collaterals and intrahepatic shunts carry part of the portal blood around the hepatocytes ([Iwakiri 2014](#ref-Iwakiri2014)). The capillarization of the sinusoids (the loss of their fenestrae) hinders the exchange between blood and hepatocytes as well. The hepatic artery partly compensates the lost portal flow (the hepatic arterial buffer response).
- **Systemic consequences.** The synthesis of albumin and of the clotting factors falls (hypoalbuminaemia, a prolonged INR), the splanchnic vasodilation leads to a hyperdynamic circulation with a raised cardiac output, the renal perfusion and the GFR fall up to the hepatorenal syndrome ([Ginès 2018](#ref-Gines2018)), and anaemia and ascites are common.

Which of these dominates depends on the drug. The well-stirred model of the liver ([Wilkinson 1975](#ref-Wilkinson1975)), `CL_H = Q*fu*CL_int/(Q + fu*CL_int)`, separates the cases: a drug of high extraction (`fu*CL_int >> Q`) is cleared at the blood flow, and the shunt lets its oral dose bypass the first pass, so its bioavailability can rise several-fold; a drug of low extraction is cleared at `fu*CL_int`, which falls with the enzymes and rises with the fu of a lower albumin. The library has each of these mechanisms as a factor of the severity: the amount of liver tissue `f_functional` (the liver volume) times the abundance of every protein per amount of tissue, which counts the loss of hepatocytes once (the abundance per gram of cirrhotic tissue already holds it), the shunt `f_shunt`, the hepatic arterial and renal blood flow, the cardiac index, the GFR, the albumin and AGP of the binding and the hematocrit ([Impairment](impairment.md#liver-cirrhosis)). The physiology-based models of cirrhosis of PK-Sim and Simcyp are built in the same way ([Edginton 2008](#ref-Edginton2008), [Johnson 2010](#ref-Johnson2010)).

### Child-Pugh classification

The Child-Pugh score ([Pugh 1973](#ref-Pugh1973), a modification of [Child and Turcotte 1964](#ref-ChildTurcotte1964)) adds 1 to 3 points for each of five criteria, two clinical and three laboratory values, from 5 to 15 points; the class is A for 5-6, B for 7-9 and C for 10-15 points. The guidances on hepatic impairment grade the subjects of a pharmacokinetic study by it, mild A, moderate B and severe C ([FDA 2003](#ref-FDA2003Hepatic), [EMA 2005](#ref-EMA2005Hepatic)). The score was made to predict the outcome of surgery and correlates only loosely with the clearance of a drug: two of its criteria are clinical judgements, and the classes are wide. The score has no healthy class: a subject without liver disease has 5 points, the reference of the models (`CP_severity` 0).

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

### The MELD score

The MELD score (model for end-stage liver disease; [Malinchoc 2000](#ref-Malinchoc2000), [Kamath 2001](#ref-Kamath2001)) predicts the survival of a patient with cirrhosis from three laboratory values and since 2002 ranks the patients on the waiting list for a liver transplantation ([Wiesner 2003](#ref-Wiesner2003)); MELD-Na adds the serum sodium ([Kim 2008](#ref-Kim2008)). It is continuous: `3.78 ln(bilirubin [mg/dL]) + 11.2 ln(INR) + 9.57 ln(creatinine [mg/dL]) + 6.43`, every value bounded below by 1 (UNOS also bounds the creatinine above by 4 mg/dL, 4 on dialysis, and multiplies by 10 and rounds the score of `0.957 ln(Cr) + 0.378 ln(bili) + 1.120 ln(INR) + 0.643`). The score of a healthy subject is 6.43. The library maps it onto `CP_severity` through the MELD of each class (`cp_severity_of_meld`, [The MELD score](impairment.md#the-meld-score)).

## Chronic kidney disease

### Background

In chronic kidney disease the nephrons are lost over months to years (diabetes, hypertension, glomerulonephritis, polycystic and tubulointerstitial disease). The remaining nephrons adapt and keep filtering and secreting in proportion to each other, the intact nephron hypothesis ([Bricker 1960](#ref-Bricker1960)): the tubular secretion and reabsorption of a drug fall with the GFR, which is why the library scales every metabolism and transport of the kidney with the renal function. Beyond the renal clearance, CKD changes the body ([Verbeeck 2009](#ref-Verbeeck2009)):

- **Renal blood flow and kidney size.** The perfusion of the kidney falls with the loss of nephrons and the kidneys shrink in the late stages.
- **Anaemia.** The failing kidney produces less erythropoietin, the hemoglobin and the hematocrit fall from about stage 3 on.
- **Binding to plasma proteins.** Uraemic toxins (indoxyl sulfate, hippuric acid) occupy the binding sites of albumin and the albumin falls in inflammation and proteinuria, so the unbound fraction of an acidic drug rises; the acute-phase protein alpha-1-acid glycoprotein rises, so that of some basic drugs falls.
- **Non-renal clearance.** Uraemia inhibits hepatic enzymes and transporters (CYP2D6, CYP3A4, OATP1B) of drugs which the kidney does not clear ([Nolin 2008](#ref-Nolin2008)), measured as the unbound clearance of model drugs ([Yoshida 2016](#ref-Yoshida2016), [Tan 2018](#ref-Tan2018)).
- **Absorption.** The gastric emptying is delayed in late CKD.

The library has the renal function `f_renal_function` (the GFR relative to healthy) with its eGFR, and the anaemia, the renal blood flow, the binding and the hepatic proteins as factors of the eGFR ([Impairment](impairment.md#chronic-kidney-disease)).

### KDIGO classification

CKD is an abnormality of the structure or the function of the kidneys present for more than 3 months ([Levey 2005](#ref-Levey2005), [KDIGO 2012](#ref-KDIGO2012)). It is graded by the cause, the GFR category and the albuminuria category. The GFR categories by the GFR per 1.73 m2 body surface area (the eGFR, estimated from creatinine or cystatin C, or measured) are G1 (>= 90, normal or high), G2 (60-89, mildly decreased), G3a (45-59, mildly to moderately decreased), G3b (30-44, moderately to severely decreased), G4 (15-29, severely decreased) and G5 (< 15, kidney failure, G5D on dialysis); G1 and G2 are CKD only with a marker of kidney damage. The albuminuria categories by the albumin excretion rate (mg/24 h) or the albumin-to-creatinine ratio (mg/g) are A1 (< 30, normal to mildly increased), A2 (30-300, moderately increased) and A3 (> 300, severely increased); the library has no albuminuria. The FDA guidance on renal impairment ([FDA 2024](#ref-FDA2024Renal)) has the stages normal (>= 90), mild (60-89), moderate (30-59, G3a and G3b), severe (15-29) and kidney failure (< 15, not on dialysis, and on dialysis apart), by the eGFR or the creatinine clearance; the EMA guidance ([EMA 2015](#ref-EMA2015Renal)) has the same limits.

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

## Heart failure

### Background

In chronic heart failure the heart cannot pump the blood the body needs at normal filling pressures ([Bozkurt 2021](#ref-Bozkurt2021)). With a falling cardiac output the sympathetic nervous system and the renin-angiotensin-aldosterone system redistribute the blood: the flow to the heart and the brain is kept, the flow to the splanchnic organs, the kidneys, the skin and the muscle falls ([Leithe 1984](#ref-Leithe1984)). The congestion of the venous system (a raised central venous pressure) adds to it: the congested liver clears less, the congested gut may absorb less, and the venous pressure lowers the perfusion pressure of the kidney ([Shammas 1988](#ref-Shammas1988)). The kidney keeps its GFR longer than its blood flow by a higher filtration fraction; in advanced failure the GFR falls as well, the cardiorenal syndrome ([Ronco 2008](#ref-Ronco2008)). For a drug of high hepatic extraction the clearance follows the hepatic blood flow and falls with the severity (lidocaine is the classical example); the volume of distribution may change with the congestion.

The library has the blood flows of the liver and the portal organs, the kidney, and the skin, adipose and muscle as factors of the NYHA class; the cardiac output is their sum and the GFR follows the renal blood flow ([Impairment](impairment.md#heart-failure)).

### NYHA classification

The NYHA functional classes ([Criteria Committee of the New York Heart Association, 9th edition 1994](#ref-NYHA1994)) grade the symptoms of heart failure at physical activity. The ACC/AHA stages grade its development instead ([Heidenreich 2022](#ref-Heidenreich2022)): A at risk without structural heart disease or symptoms, B structural heart disease without symptoms (pre-heart failure), C structural heart disease with current or previous symptoms and D advanced heart failure refractory to treatment; the NYHA classes grade the patients of the stages C and D, a patient of stage D is mostly NYHA IV. The ejection fraction separates heart failure with reduced (<= 40 %), mildly reduced (41-49 %) and preserved (>= 50 %) ejection fraction ([Bozkurt 2021](#ref-Bozkurt2021)); the library has no such distinction.

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

## Age

### Background

The physiology is a continuous function of the age ([Physiology over age](physiology.md#physiology-over-age)); the classes of the age are those of the guidelines for the studies of special populations. ICH E11 classifies the paediatric patients by their development ([ICH E11](#ref-ICHE11)): preterm newborn infants, term newborn infants (0 to 27 days), infants and toddlers (28 days to 23 months), children (2 to 11 years) and adolescents (12 to 16-18 years, depending on the region; 12 to 17 years in the EMA guideline on paediatric pharmacokinetics), and notes that the age may be examined as a continuous covariate instead. ICH E7 defines the geriatric population as 65 years and older and asks for patients of 75 years and above; its questions and answers of 2010 name the subgroups 65-74, 75-84 and 85 years and older ([ICH E7](#ref-ICHE7)). The preterm newborn needs the gestational age, which the library does not model.

### Age classes

The value of a class is a representative age, the adult the adult of ICRP 89 (35 years); `class_changes` sets `AGE` and the body weight and height of the reference individual of the age and the sex (`reference_individual`), and a [virtual population](populations.md) draws the ages of a class uniformly with the body size of its population.

<!-- generated: classes-age (scripts/pbpk_tables.py), do not edit -->
*Generated from `pkpdsim.impairment.AGE_CLASSES` by `scripts/pbpk_tables.py`.*

| class | definition | age [yr] | AGE |
| --- | --- | --- | ---: |
| <span class="pbpk-swatch" style="background:#c6dbef" title="#c6dbef"></span> `newborn` Term newborn | term newborn infants, 0 to 27 days (ICH E11); preterm newborn infants need the gestational age, which the library does not model | 0-27 days | 0.03696 |
| <span class="pbpk-swatch" style="background:#6baed6" title="#6baed6"></span> `infant` Infant and toddler | infants and toddlers, 28 days to 23 months (ICH E11) | 0.07666-1 | 1 |
| <span class="pbpk-swatch" style="background:#2171b5" title="#2171b5"></span> `child` Child | children, 2 to 11 years (ICH E11) | 2-11 | 6.5 |
| <span class="pbpk-swatch" style="background:#08306b" title="#08306b"></span> `adolescent` Adolescent | adolescents, 12 to 17 years (ICH E11: 12 to 16-18 years by region; the EMA guideline on paediatric pharmacokinetics 12 to 17) | 12-17 | 15 |
| <span class="pbpk-swatch" style="background:black" title="black"></span> `adult` Adult | adults, 18 to 64 years; the value of the class is the adult of ICRP 89 (35 years) | 18-64 | 35 |
| <span class="pbpk-swatch" style="background:#fdae6b" title="#fdae6b"></span> `young_old` 65-74 years | geriatric patients of 65 to 74 years (ICH E7 Questions and Answers 2010; ICH E7: the geriatric population is 65 years or older) | 65-74 | 70 |
| <span class="pbpk-swatch" style="background:#e6550d" title="#e6550d"></span> `old_old` 75-84 years | geriatric patients of 75 to 84 years (ICH E7 Questions and Answers 2010) | 75-84 | 80 |
| <span class="pbpk-swatch" style="background:#7f2704" title="#7f2704"></span> `oldest_old` 85 years and older | geriatric patients of 85 years and older (ICH E7 Questions and Answers 2010); the physiology beyond 85 years is extrapolated | >= 85 | 90 |
<!-- end generated -->

## Body mass index

### Background

Overweight and obesity are graded by the body mass index, `BW/height^2`. WHO Technical Report Series 894 gives the international classification ([WHO 2000](#ref-WHO2000)): underweight below 18.5 (severe below 16, moderate 16 to 16.9, mild 17 to 18.4), normal weight 18.5 to 24.9, pre-obese (overweight) 25.0 to 29.9, obesity class I 30.0 to 34.9, class II 35.0 to 39.9 and class III (severe or morbid obesity) 40.0 and above. The WHO Expert Consultation of 2004 retained those cut-offs and added the public health action points 23.0, 27.5, 32.5 and 37.5 kg/m2 for Asian populations, whose point of observed risk lies between a body mass index of 22 and 25 ([WHO Asian 2004](#ref-WHOAsian2004)). What the physiology of an obese individual does beyond the larger adipose mass, and what the library changes with it, is on the page [Adipose tissue and obesity](adipose.md).

### Classes of the body mass index

The value of a class is the middle of its range (45 kg/m2 for class III, which is open); the class of normal weight has no value of its own and is the reference individual of the build. `class_changes` sets the height of the reference individual of the sex and the age and the body weight at the body mass index of the class, so a scan over the classes is a scan over body weights at a fixed height.

<!-- generated: classes-bmi (scripts/pbpk_tables.py), do not edit -->
*Generated from `pkpdsim.impairment.BMI_CLASSES` by `scripts/pbpk_tables.py`.*

| class | definition | BMI [kg/m2] | BMI | BW [kg] |
| --- | --- | --- | --- | --- |
| <span class="pbpk-swatch" style="background:black" title="black"></span> `normal` Normal weight | a body mass index of 18.5 to 24.9 kg/m2 (WHO Technical Report Series 894); the class has no value of its own, it is the reference individual of the build (a body mass index of 23.6 kg/m2 for the man and 22.6 for the woman of ICRP 89) | 18.5-24 | reference | male 73, female 60 |
| <span class="pbpk-swatch" style="background:#6baed6" title="#6baed6"></span> `underweight` Underweight | a body mass index below 18.5 kg/m2 (WHO: severe below 16, moderate 16 to 16.9, mild 17 to 18.4); the value of the class is the middle of the mild and moderate range | < 18.5 | 17 | male 52.66, female 45.17 |
| <span class="pbpk-swatch" style="background:#fdd0a2" title="#fdd0a2"></span> `overweight` Overweight | a body mass index of 25.0 to 29.9 kg/m2, pre-obese in WHO Technical Report Series 894; the value of the class is its middle | 25-29 | 27.5 | male 85.18, female 73.06 |
| <span class="pbpk-swatch" style="background:#fdae6b" title="#fdae6b"></span> `obesity_1` Obesity class I | a body mass index of 30.0 to 34.9 kg/m2; the value of the class is its middle | 30-34 | 32.5 | male 100.7, female 86.35 |
| <span class="pbpk-swatch" style="background:#e6550d" title="#e6550d"></span> `obesity_2` Obesity class II | a body mass index of 35.0 to 39.9 kg/m2; the value of the class is its middle | 35-39 | 37.5 | male 116.2, female 99.63 |
| <span class="pbpk-swatch" style="background:#7f2704" title="#7f2704"></span> `obesity_3` Obesity class III | a body mass index of 40.0 kg/m2 and above, severe or morbid obesity; the value of the class is 45 kg/m2, and the adipose model of the library is valid to about 40 (the equations of the fat-free mass and of the obese repository both saturate above it) | >= 40 | 45 | male 139.4, female 119.6 |
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
<polyline points="27,57.19 61,63.78 95,71.09 129,78.87" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="63.78" r="3.5" fill="currentColor"><title>A: 0.8982</title></circle>
<circle cx="95" cy="71.09" r="3.5" fill="currentColor"><title>B: 0.7854</title></circle>
<circle cx="129" cy="78.87" r="3.5" fill="currentColor"><title>C: 0.6655</title></circle>
<text x="149" y="81.87" font-size="11" opacity="0.7">0.6655</text>
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
<line x1="10" x2="146" y1="191.5" y2="191.5" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="194.5" font-size="11" opacity="0.7">1</text>
<polyline points="10,191.5 44,191.5 61,191.5 95,191.5 120.5,167.2 137.5,167.2 146,167.2" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="191.5" r="3.5" fill="currentColor"><title>mild: 1</title></circle>
<circle cx="95" cy="191.5" r="3.5" fill="currentColor"><title>moderate: 1</title></circle>
<circle cx="120.5" cy="167.2" r="3.5" fill="currentColor"><title>severe: 1.6</title></circle>
<circle cx="137.5" cy="167.2" r="3.5" fill="currentColor"><title>end stage: 1.6</title></circle>
<text x="149" y="170.2" font-size="11" opacity="0.7">1.6</text>
<text x="190" y="155" font-size="13">si_transit_time</text>
<rect x="190" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="224" y="235" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="258" y="235" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="275" y="235" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="292" y="235" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="309" y="235" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="190" x2="326" y1="185.7" y2="185.7" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="329" y="188.7" font-size="11" opacity="0.7">1</text>
<polyline points="190,185.7 224,185.7 241,185.7 275,185.7 300.5,167.2 317.5,167.2 326,167.2" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="241" cy="185.7" r="3.5" fill="currentColor"><title>mild: 1</title></circle>
<circle cx="275" cy="185.7" r="3.5" fill="currentColor"><title>moderate: 1</title></circle>
<circle cx="300.5" cy="167.2" r="3.5" fill="currentColor"><title>severe: 1.4</title></circle>
<circle cx="317.5" cy="167.2" r="3.5" fill="currentColor"><title>end stage: 1.4</title></circle>
<text x="329" y="170.2" font-size="11" opacity="0.7">1.4</text>
<text x="370" y="155" font-size="13">binding: albumin</text>
<rect x="370" y="235" width="34" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>G1 normal</title></rect>
<rect x="404" y="235" width="34" height="5" fill="#66c2a4" fill-opacity="0.9"><title>G2 mild</title></rect>
<rect x="438" y="235" width="17" height="5" fill="#41ae76" fill-opacity="0.9"><title>G3a mild-moderate</title></rect>
<rect x="455" y="235" width="17" height="5" fill="#238b45" fill-opacity="0.9"><title>G3b moderate-severe</title></rect>
<rect x="472" y="235" width="17" height="5" fill="#006d2c" fill-opacity="0.9"><title>G4 severe</title></rect>
<rect x="489" y="235" width="17" height="5" fill="#00441b" fill-opacity="0.9"><title>G5 kidney failure</title></rect>
<line x1="370" x2="506" y1="167.2" y2="167.2" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="509" y="170.2" font-size="11" opacity="0.7">1</text>
<polyline points="370,167.2 404,167.2 421,169.4 455,171.7 480.5,172.9 497.5,176.6 506,176.6" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="421" cy="169.4" r="3.5" fill="currentColor"><title>mild: 0.9654</title></circle>
<circle cx="455" cy="171.7" r="3.5" fill="currentColor"><title>moderate: 0.9308</title></circle>
<circle cx="480.5" cy="172.9" r="3.5" fill="currentColor"><title>severe: 0.912</title></circle>
<circle cx="497.5" cy="176.6" r="3.5" fill="currentColor"><title>end stage: 0.8553</title></circle>
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
<line x1="10" x2="146" y1="297.4" y2="297.4" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="149" y="300.4" font-size="11" opacity="0.7">1</text>
<polyline points="10,297.4 44,297.4 61,292.4 95,287.3 120.5,282.2 137.5,277.2 146,277.2" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="61" cy="292.4" r="3.5" fill="currentColor"><title>mild: 1.113</title></circle>
<circle cx="95" cy="287.3" r="3.5" fill="currentColor"><title>moderate: 1.227</title></circle>
<circle cx="120.5" cy="282.2" r="3.5" fill="currentColor"><title>severe: 1.341</title></circle>
<circle cx="137.5" cy="277.2" r="3.5" fill="currentColor"><title>end stage: 1.454</title></circle>
<text x="149" y="280.2" font-size="11" opacity="0.7">1.454</text>
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
<polyline points="23.6,57.19 50.8,65.43 78,73.67 105.2,87.86 132.4,91.78" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="50.8" cy="65.43" r="3.5" fill="currentColor"><title>I: 0.8728</title></circle>
<circle cx="78" cy="73.67" r="3.5" fill="currentColor"><title>II: 0.7456</title></circle>
<circle cx="105.2" cy="87.86" r="3.5" fill="currentColor"><title>III: 0.5267</title></circle>
<circle cx="132.4" cy="91.78" r="3.5" fill="currentColor"><title>IV: 0.4662</title></circle>
<text x="149" y="94.78" font-size="11" opacity="0.7">0.4662</text>
<text x="190" y="45" font-size="13">kidney_flow</text>
<rect x="190" y="125" width="27.2" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="217.2" y="125" width="27.2" height="5" fill="#FDDBC7" fill-opacity="0.9"><title>NYHA I</title></rect>
<rect x="244.4" y="125" width="27.2" height="5" fill="#F4A582" fill-opacity="0.9"><title>NYHA II</title></rect>
<rect x="271.6" y="125" width="27.2" height="5" fill="#D6604D" fill-opacity="0.9"><title>NYHA III</title></rect>
<rect x="298.8" y="125" width="27.2" height="5" fill="#67001F" fill-opacity="0.9"><title>NYHA IV</title></rect>
<line x1="190" x2="326" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="329" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="203.6,57.19 230.8,65.01 258,72.84 285.2,86.4 312.4,81.18" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="230.8" cy="65.01" r="3.5" fill="currentColor"><title>I: 0.8792</title></circle>
<circle cx="258" cy="72.84" r="3.5" fill="currentColor"><title>II: 0.7584</title></circle>
<circle cx="285.2" cy="86.4" r="3.5" fill="currentColor"><title>III: 0.5492</title></circle>
<circle cx="312.4" cy="81.18" r="3.5" fill="currentColor"><title>IV: 0.6298</title></circle>
<text x="329" y="89.4" font-size="11" opacity="0.7">0.5492</text>
<text x="370" y="45" font-size="13">peripheral_flow</text>
<rect x="370" y="125" width="27.2" height="5" fill="#9e9e9e" fill-opacity="0.35"><title>Healthy</title></rect>
<rect x="397.2" y="125" width="27.2" height="5" fill="#FDDBC7" fill-opacity="0.9"><title>NYHA I</title></rect>
<rect x="424.4" y="125" width="27.2" height="5" fill="#F4A582" fill-opacity="0.9"><title>NYHA II</title></rect>
<rect x="451.6" y="125" width="27.2" height="5" fill="#D6604D" fill-opacity="0.9"><title>NYHA III</title></rect>
<rect x="478.8" y="125" width="27.2" height="5" fill="#67001F" fill-opacity="0.9"><title>NYHA IV</title></rect>
<line x1="370" x2="506" y1="57.19" y2="57.19" stroke="currentColor" stroke-opacity="0.35" stroke-dasharray="3 3"/>
<text x="509" y="60.19" font-size="11" opacity="0.7">1</text>
<polyline points="383.6,57.19 410.8,68.57 438,79.95 465.2,85 492.4,91.84" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>
<circle cx="410.8" cy="68.57" r="3.5" fill="currentColor"><title>I: 0.8244</title></circle>
<circle cx="438" cy="79.95" r="3.5" fill="currentColor"><title>II: 0.6488</title></circle>
<circle cx="465.2" cy="85" r="3.5" fill="currentColor"><title>III: 0.5709</title></circle>
<circle cx="492.4" cy="91.84" r="3.5" fill="currentColor"><title>IV: 0.4654</title></circle>
<text x="509" y="94.84" font-size="11" opacity="0.7">0.4654</text>
</svg>
</div>
<!-- end generated -->

## Continuous severities and combinations

The classes are points of continuous severities: a model takes any value of `CP_severity` from 0 to 3, of `HF_severity` from 0 to 4 and of `f_renal_function` from 0 (its eGFR is `f_renal_function` times the eGFR of the normal function of the build). Every factor is piecewise linear between the classes (`CP_severity` 1.5 is half way between A and B in every factor) and between the stage anchors of the eGFR (75, 45, 22.5 and 7.5 ml/min/1.73 m2, 1 from 90 up), so a scan or a fit of a severity sees no step. The figures of the factors above are these functions.

The conditions occur together: cirrhosis lowers the GFR up to the hepatorenal syndrome ([Ginès 2018](#ref-Gines2018)), heart failure up to the cardiorenal syndrome ([Ronco 2008](#ref-Ronco2008)), and CKD and heart failure share their causes. The conditions combine without a rule of their own: each condition has its factors (`fcp_<key>`, `fckd_<key>`, `fhf_<key>`), and the factor of a quantity is the product of the factors of the conditions (`fimp_<id>`, see [Impairment](impairment.md)). Cirrhosis with a declining renal function (hepatorenal) is `CP_severity` with `f_renal_function`, heart failure with it (cardiorenal) `HF_severity` with `f_renal_function`, all three together as well:

- The GFR is `GFR_base * f_renal_function`, the GFR at the renal function 1 (`GFR_base`) carrying the functional decline of cirrhosis (`fcp_gfr`) and heart failure (`fhf_kidney_flow`), `f_renal_function` the intrinsic renal function. A measured GFR of a cirrhotic or heart failure patient is `renal_function(clearance, bw=..., height=..., cp_severity=..., hf_severity=...)`, which divides by the `GFR_base` of the impaired individual.
- No effect is counted twice. The CKD factors (anaemia, binding, renal blood flow) follow the eGFR of the intrinsic renal function, not the functional decline of cirrhosis or heart failure: the factors of the Child-Pugh and NYHA classes are measured in patients of these classes, their renal decline included, so the CKD factors on top of them would count it twice. A patient whose kidneys are damaged besides (a hepatorenal syndrome of long duration, CKD of another cause) has both, the class of the liver or the heart and the renal function of the intrinsic kidney disease.
- The flows multiply (the renal blood flow `fcp_kidney_flow * fckd_kidney_flow * fhf_kidney_flow`), the hematocrit (`fcp_hematocrit * fckd_hematocrit`) and the binding (the albumin of cirrhosis times the albumin and its binding function of CKD) as well.

<!-- generated: combinations-figure (scripts/pbpk_tables.py), do not edit -->
<div class="pbpk-figure pbpk-wide" markdown="0">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 760 310" role="img" aria-label="Combined impairments: GFR, renal blood flow and fu relative to healthy" font-family="system-ui, sans-serif" font-size="14" fill="currentColor">
<title>Combined impairments: GFR, renal blood flow and fu relative to healthy</title>
<text x="50" y="18" font-weight="600">GFR, hepatorenal</text>
<text transform="translate(22 140) rotate(-90)" text-anchor="middle" font-size="12">eGFR (ml/min/1.73 m2)</text>
<rect x="50" y="40" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 0, eGFR 117.5: 1.03</title></rect>
<rect x="58.33" y="40" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 0.13, eGFR 117.5: 1.03</title></rect>
<rect x="66.67" y="40" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 0.261, eGFR 117.5: 1.03</title></rect>
<rect x="75" y="40" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 0.391, eGFR 117.5: 1.03</title></rect>
<rect x="83.33" y="40" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 0.522, eGFR 117.5: 1.03</title></rect>
<rect x="91.67" y="40" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 0.652, eGFR 117.5: 1.03</title></rect>
<rect x="100" y="40" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 0.783, eGFR 117.5: 1.03</title></rect>
<rect x="108.3" y="40" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 0.913, eGFR 117.5: 1.03</title></rect>
<rect x="116.7" y="40" width="8.633" height="8.633" fill="#09336f"><title>CP_severity 1.043, eGFR 117.5: 1.01</title></rect>
<rect x="125" y="40" width="8.633" height="8.633" fill="#0d3e7b"><title>CP_severity 1.174, eGFR 117.5: 0.973</title></rect>
<rect x="133.3" y="40" width="8.633" height="8.633" fill="#114886"><title>CP_severity 1.304, eGFR 117.5: 0.933</title></rect>
<rect x="141.7" y="40" width="8.633" height="8.633" fill="#155292"><title>CP_severity 1.435, eGFR 117.5: 0.892</title></rect>
<rect x="150" y="40" width="8.633" height="8.633" fill="#195c9e"><title>CP_severity 1.565, eGFR 117.5: 0.852</title></rect>
<rect x="158.3" y="40" width="8.633" height="8.633" fill="#1d67a9"><title>CP_severity 1.696, eGFR 117.5: 0.812</title></rect>
<rect x="166.7" y="40" width="8.633" height="8.633" fill="#2171b5"><title>CP_severity 1.826, eGFR 117.5: 0.772</title></rect>
<rect x="175" y="40" width="8.633" height="8.633" fill="#2d7bba"><title>CP_severity 1.957, eGFR 117.5: 0.732</title></rect>
<rect x="183.3" y="40" width="8.633" height="8.633" fill="#3985c0"><title>CP_severity 2.087, eGFR 117.5: 0.688</title></rect>
<rect x="191.7" y="40" width="8.633" height="8.633" fill="#4790c6"><title>CP_severity 2.217, eGFR 117.5: 0.642</title></rect>
<rect x="200" y="40" width="8.633" height="8.633" fill="#549bcc"><title>CP_severity 2.348, eGFR 117.5: 0.597</title></rect>
<rect x="208.3" y="40" width="8.633" height="8.633" fill="#61a6d2"><title>CP_severity 2.478, eGFR 117.5: 0.551</title></rect>
<rect x="216.7" y="40" width="8.633" height="8.633" fill="#6fb0d7"><title>CP_severity 2.609, eGFR 117.5: 0.506</title></rect>
<rect x="225" y="40" width="8.633" height="8.633" fill="#7fb8dc"><title>CP_severity 2.739, eGFR 117.5: 0.46</title></rect>
<rect x="233.3" y="40" width="8.633" height="8.633" fill="#8fc0e0"><title>CP_severity 2.87, eGFR 117.5: 0.415</title></rect>
<rect x="241.7" y="40" width="8.633" height="8.633" fill="#a0c8e4"><title>CP_severity 3, eGFR 117.5: 0.369</title></rect>
<rect x="50" y="48.33" width="8.633" height="8.633" fill="#0c3b78"><title>CP_severity 0, eGFR 112.5: 0.983</title></rect>
<rect x="58.33" y="48.33" width="8.633" height="8.633" fill="#0c3b78"><title>CP_severity 0.13, eGFR 112.5: 0.983</title></rect>
<rect x="66.67" y="48.33" width="8.633" height="8.633" fill="#0c3b78"><title>CP_severity 0.261, eGFR 112.5: 0.983</title></rect>
<rect x="75" y="48.33" width="8.633" height="8.633" fill="#0c3b78"><title>CP_severity 0.391, eGFR 112.5: 0.983</title></rect>
<rect x="83.33" y="48.33" width="8.633" height="8.633" fill="#0c3b78"><title>CP_severity 0.522, eGFR 112.5: 0.983</title></rect>
<rect x="91.67" y="48.33" width="8.633" height="8.633" fill="#0c3b78"><title>CP_severity 0.652, eGFR 112.5: 0.983</title></rect>
<rect x="100" y="48.33" width="8.633" height="8.633" fill="#0c3b78"><title>CP_severity 0.783, eGFR 112.5: 0.983</title></rect>
<rect x="108.3" y="48.33" width="8.633" height="8.633" fill="#0c3b78"><title>CP_severity 0.913, eGFR 112.5: 0.983</title></rect>
<rect x="116.7" y="48.33" width="8.633" height="8.633" fill="#0e3e7b"><title>CP_severity 1.043, eGFR 112.5: 0.97</title></rect>
<rect x="125" y="48.33" width="8.633" height="8.633" fill="#114887"><title>CP_severity 1.174, eGFR 112.5: 0.931</title></rect>
<rect x="133.3" y="48.33" width="8.633" height="8.633" fill="#155292"><title>CP_severity 1.304, eGFR 112.5: 0.893</title></rect>
<rect x="141.7" y="48.33" width="8.633" height="8.633" fill="#195c9d"><title>CP_severity 1.435, eGFR 112.5: 0.854</title></rect>
<rect x="150" y="48.33" width="8.633" height="8.633" fill="#1d66a8"><title>CP_severity 1.565, eGFR 112.5: 0.816</title></rect>
<rect x="158.3" y="48.33" width="8.633" height="8.633" fill="#2070b3"><title>CP_severity 1.696, eGFR 112.5: 0.777</title></rect>
<rect x="166.7" y="48.33" width="8.633" height="8.633" fill="#2a79b9"><title>CP_severity 1.826, eGFR 112.5: 0.739</title></rect>
<rect x="175" y="48.33" width="8.633" height="8.633" fill="#3682be"><title>CP_severity 1.957, eGFR 112.5: 0.701</title></rect>
<rect x="183.3" y="48.33" width="8.633" height="8.633" fill="#428cc4"><title>CP_severity 2.087, eGFR 112.5: 0.659</title></rect>
<rect x="191.7" y="48.33" width="8.633" height="8.633" fill="#4e96c9"><title>CP_severity 2.217, eGFR 112.5: 0.615</title></rect>
<rect x="200" y="48.33" width="8.633" height="8.633" fill="#5ba1cf"><title>CP_severity 2.348, eGFR 112.5: 0.572</title></rect>
<rect x="208.3" y="48.33" width="8.633" height="8.633" fill="#68abd5"><title>CP_severity 2.478, eGFR 112.5: 0.528</title></rect>
<rect x="216.7" y="48.33" width="8.633" height="8.633" fill="#77b4d9"><title>CP_severity 2.609, eGFR 112.5: 0.484</title></rect>
<rect x="225" y="48.33" width="8.633" height="8.633" fill="#86bbdd"><title>CP_severity 2.739, eGFR 112.5: 0.441</title></rect>
<rect x="233.3" y="48.33" width="8.633" height="8.633" fill="#96c3e2"><title>CP_severity 2.87, eGFR 112.5: 0.397</title></rect>
<rect x="241.7" y="48.33" width="8.633" height="8.633" fill="#a5cbe6"><title>CP_severity 3, eGFR 112.5: 0.354</title></rect>
<rect x="50" y="56.67" width="8.633" height="8.633" fill="#114684"><title>CP_severity 0, eGFR 107.5: 0.939</title></rect>
<rect x="58.33" y="56.67" width="8.633" height="8.633" fill="#114684"><title>CP_severity 0.13, eGFR 107.5: 0.939</title></rect>
<rect x="66.67" y="56.67" width="8.633" height="8.633" fill="#114684"><title>CP_severity 0.261, eGFR 107.5: 0.939</title></rect>
<rect x="75" y="56.67" width="8.633" height="8.633" fill="#114684"><title>CP_severity 0.391, eGFR 107.5: 0.939</title></rect>
<rect x="83.33" y="56.67" width="8.633" height="8.633" fill="#114684"><title>CP_severity 0.522, eGFR 107.5: 0.939</title></rect>
<rect x="91.67" y="56.67" width="8.633" height="8.633" fill="#114684"><title>CP_severity 0.652, eGFR 107.5: 0.939</title></rect>
<rect x="100" y="56.67" width="8.633" height="8.633" fill="#114684"><title>CP_severity 0.783, eGFR 107.5: 0.939</title></rect>
<rect x="108.3" y="56.67" width="8.633" height="8.633" fill="#114684"><title>CP_severity 0.913, eGFR 107.5: 0.939</title></rect>
<rect x="116.7" y="56.67" width="8.633" height="8.633" fill="#124988"><title>CP_severity 1.043, eGFR 107.5: 0.927</title></rect>
<rect x="125" y="56.67" width="8.633" height="8.633" fill="#155393"><title>CP_severity 1.174, eGFR 107.5: 0.89</title></rect>
<rect x="133.3" y="56.67" width="8.633" height="8.633" fill="#195c9d"><title>CP_severity 1.304, eGFR 107.5: 0.853</title></rect>
<rect x="141.7" y="56.67" width="8.633" height="8.633" fill="#1d66a8"><title>CP_severity 1.435, eGFR 107.5: 0.816</title></rect>
<rect x="150" y="56.67" width="8.633" height="8.633" fill="#206fb3"><title>CP_severity 1.565, eGFR 107.5: 0.78</title></rect>
<rect x="158.3" y="56.67" width="8.633" height="8.633" fill="#2978b9"><title>CP_severity 1.696, eGFR 107.5: 0.743</title></rect>
<rect x="166.7" y="56.67" width="8.633" height="8.633" fill="#3481bd"><title>CP_severity 1.826, eGFR 107.5: 0.706</title></rect>
<rect x="175" y="56.67" width="8.633" height="8.633" fill="#3f89c2"><title>CP_severity 1.957, eGFR 107.5: 0.669</title></rect>
<rect x="183.3" y="56.67" width="8.633" height="8.633" fill="#4a93c7"><title>CP_severity 2.087, eGFR 107.5: 0.629</title></rect>
<rect x="191.7" y="56.67" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 2.217, eGFR 107.5: 0.588</title></rect>
<rect x="200" y="56.67" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 2.348, eGFR 107.5: 0.546</title></rect>
<rect x="208.3" y="56.67" width="8.633" height="8.633" fill="#6fb0d7"><title>CP_severity 2.478, eGFR 107.5: 0.505</title></rect>
<rect x="216.7" y="56.67" width="8.633" height="8.633" fill="#7eb8db"><title>CP_severity 2.609, eGFR 107.5: 0.463</title></rect>
<rect x="225" y="56.67" width="8.633" height="8.633" fill="#8dbfdf"><title>CP_severity 2.739, eGFR 107.5: 0.421</title></rect>
<rect x="233.3" y="56.67" width="8.633" height="8.633" fill="#9cc6e3"><title>CP_severity 2.87, eGFR 107.5: 0.38</title></rect>
<rect x="241.7" y="56.67" width="8.633" height="8.633" fill="#abcee8"><title>CP_severity 3, eGFR 107.5: 0.338</title></rect>
<rect x="50" y="65" width="8.633" height="8.633" fill="#155191"><title>CP_severity 0, eGFR 102.5: 0.895</title></rect>
<rect x="58.33" y="65" width="8.633" height="8.633" fill="#155191"><title>CP_severity 0.13, eGFR 102.5: 0.895</title></rect>
<rect x="66.67" y="65" width="8.633" height="8.633" fill="#155191"><title>CP_severity 0.261, eGFR 102.5: 0.895</title></rect>
<rect x="75" y="65" width="8.633" height="8.633" fill="#155191"><title>CP_severity 0.391, eGFR 102.5: 0.895</title></rect>
<rect x="83.33" y="65" width="8.633" height="8.633" fill="#155191"><title>CP_severity 0.522, eGFR 102.5: 0.895</title></rect>
<rect x="91.67" y="65" width="8.633" height="8.633" fill="#155191"><title>CP_severity 0.652, eGFR 102.5: 0.895</title></rect>
<rect x="100" y="65" width="8.633" height="8.633" fill="#155191"><title>CP_severity 0.783, eGFR 102.5: 0.895</title></rect>
<rect x="108.3" y="65" width="8.633" height="8.633" fill="#155191"><title>CP_severity 0.913, eGFR 102.5: 0.895</title></rect>
<rect x="116.7" y="65" width="8.633" height="8.633" fill="#165494"><title>CP_severity 1.043, eGFR 102.5: 0.884</title></rect>
<rect x="125" y="65" width="8.633" height="8.633" fill="#195d9f"><title>CP_severity 1.174, eGFR 102.5: 0.848</title></rect>
<rect x="133.3" y="65" width="8.633" height="8.633" fill="#1d66a9"><title>CP_severity 1.304, eGFR 102.5: 0.813</title></rect>
<rect x="141.7" y="65" width="8.633" height="8.633" fill="#206fb3"><title>CP_severity 1.435, eGFR 102.5: 0.778</title></rect>
<rect x="150" y="65" width="8.633" height="8.633" fill="#2978b9"><title>CP_severity 1.565, eGFR 102.5: 0.743</title></rect>
<rect x="158.3" y="65" width="8.633" height="8.633" fill="#3380bd"><title>CP_severity 1.696, eGFR 102.5: 0.708</title></rect>
<rect x="166.7" y="65" width="8.633" height="8.633" fill="#3e89c2"><title>CP_severity 1.826, eGFR 102.5: 0.673</title></rect>
<rect x="175" y="65" width="8.633" height="8.633" fill="#4891c6"><title>CP_severity 1.957, eGFR 102.5: 0.638</title></rect>
<rect x="183.3" y="65" width="8.633" height="8.633" fill="#539acb"><title>CP_severity 2.087, eGFR 102.5: 0.6</title></rect>
<rect x="191.7" y="65" width="8.633" height="8.633" fill="#5ea4d0"><title>CP_severity 2.217, eGFR 102.5: 0.56</title></rect>
<rect x="200" y="65" width="8.633" height="8.633" fill="#6aadd6"><title>CP_severity 2.348, eGFR 102.5: 0.521</title></rect>
<rect x="208.3" y="65" width="8.633" height="8.633" fill="#78b4da"><title>CP_severity 2.478, eGFR 102.5: 0.481</title></rect>
<rect x="216.7" y="65" width="8.633" height="8.633" fill="#86bbdd"><title>CP_severity 2.609, eGFR 102.5: 0.441</title></rect>
<rect x="225" y="65" width="8.633" height="8.633" fill="#94c2e1"><title>CP_severity 2.739, eGFR 102.5: 0.402</title></rect>
<rect x="233.3" y="65" width="8.633" height="8.633" fill="#a2c9e5"><title>CP_severity 2.87, eGFR 102.5: 0.362</title></rect>
<rect x="241.7" y="65" width="8.633" height="8.633" fill="#b1d0e9"><title>CP_severity 3, eGFR 102.5: 0.322</title></rect>
<rect x="50" y="73.33" width="8.633" height="8.633" fill="#195d9e"><title>CP_severity 0, eGFR 97.5: 0.852</title></rect>
<rect x="58.33" y="73.33" width="8.633" height="8.633" fill="#195d9e"><title>CP_severity 0.13, eGFR 97.5: 0.852</title></rect>
<rect x="66.67" y="73.33" width="8.633" height="8.633" fill="#195d9e"><title>CP_severity 0.261, eGFR 97.5: 0.852</title></rect>
<rect x="75" y="73.33" width="8.633" height="8.633" fill="#195d9e"><title>CP_severity 0.391, eGFR 97.5: 0.852</title></rect>
<rect x="83.33" y="73.33" width="8.633" height="8.633" fill="#195d9e"><title>CP_severity 0.522, eGFR 97.5: 0.852</title></rect>
<rect x="91.67" y="73.33" width="8.633" height="8.633" fill="#195d9e"><title>CP_severity 0.652, eGFR 97.5: 0.852</title></rect>
<rect x="100" y="73.33" width="8.633" height="8.633" fill="#195d9e"><title>CP_severity 0.783, eGFR 97.5: 0.852</title></rect>
<rect x="108.3" y="73.33" width="8.633" height="8.633" fill="#195d9e"><title>CP_severity 0.913, eGFR 97.5: 0.852</title></rect>
<rect x="116.7" y="73.33" width="8.633" height="8.633" fill="#1a5fa1"><title>CP_severity 1.043, eGFR 97.5: 0.84</title></rect>
<rect x="125" y="73.33" width="8.633" height="8.633" fill="#1e68ab"><title>CP_severity 1.174, eGFR 97.5: 0.807</title></rect>
<rect x="133.3" y="73.33" width="8.633" height="8.633" fill="#2170b4"><title>CP_severity 1.304, eGFR 97.5: 0.774</title></rect>
<rect x="141.7" y="73.33" width="8.633" height="8.633" fill="#2a78b9"><title>CP_severity 1.435, eGFR 97.5: 0.74</title></rect>
<rect x="150" y="73.33" width="8.633" height="8.633" fill="#3480bd"><title>CP_severity 1.565, eGFR 97.5: 0.707</title></rect>
<rect x="158.3" y="73.33" width="8.633" height="8.633" fill="#3d88c2"><title>CP_severity 1.696, eGFR 97.5: 0.674</title></rect>
<rect x="166.7" y="73.33" width="8.633" height="8.633" fill="#4790c6"><title>CP_severity 1.826, eGFR 97.5: 0.64</title></rect>
<rect x="175" y="73.33" width="8.633" height="8.633" fill="#5198ca"><title>CP_severity 1.957, eGFR 97.5: 0.607</title></rect>
<rect x="183.3" y="73.33" width="8.633" height="8.633" fill="#5ba1cf"><title>CP_severity 2.087, eGFR 97.5: 0.571</title></rect>
<rect x="191.7" y="73.33" width="8.633" height="8.633" fill="#66aad4"><title>CP_severity 2.217, eGFR 97.5: 0.533</title></rect>
<rect x="200" y="73.33" width="8.633" height="8.633" fill="#73b2d8"><title>CP_severity 2.348, eGFR 97.5: 0.495</title></rect>
<rect x="208.3" y="73.33" width="8.633" height="8.633" fill="#80b9dc"><title>CP_severity 2.478, eGFR 97.5: 0.458</title></rect>
<rect x="216.7" y="73.33" width="8.633" height="8.633" fill="#8ebfe0"><title>CP_severity 2.609, eGFR 97.5: 0.42</title></rect>
<rect x="225" y="73.33" width="8.633" height="8.633" fill="#9bc6e3"><title>CP_severity 2.739, eGFR 97.5: 0.382</title></rect>
<rect x="233.3" y="73.33" width="8.633" height="8.633" fill="#a9cde7"><title>CP_severity 2.87, eGFR 97.5: 0.344</title></rect>
<rect x="241.7" y="73.33" width="8.633" height="8.633" fill="#b6d3eb"><title>CP_severity 3, eGFR 97.5: 0.307</title></rect>
<rect x="50" y="81.67" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 0, eGFR 92.5: 0.808</title></rect>
<rect x="58.33" y="81.67" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 0.13, eGFR 92.5: 0.808</title></rect>
<rect x="66.67" y="81.67" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 0.261, eGFR 92.5: 0.808</title></rect>
<rect x="75" y="81.67" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 0.391, eGFR 92.5: 0.808</title></rect>
<rect x="83.33" y="81.67" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 0.522, eGFR 92.5: 0.808</title></rect>
<rect x="91.67" y="81.67" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 0.652, eGFR 92.5: 0.808</title></rect>
<rect x="100" y="81.67" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 0.783, eGFR 92.5: 0.808</title></rect>
<rect x="108.3" y="81.67" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 0.913, eGFR 92.5: 0.808</title></rect>
<rect x="116.7" y="81.67" width="8.633" height="8.633" fill="#1e6aae"><title>CP_severity 1.043, eGFR 92.5: 0.797</title></rect>
<rect x="125" y="81.67" width="8.633" height="8.633" fill="#2372b6"><title>CP_severity 1.174, eGFR 92.5: 0.766</title></rect>
<rect x="133.3" y="81.67" width="8.633" height="8.633" fill="#2c7aba"><title>CP_severity 1.304, eGFR 92.5: 0.734</title></rect>
<rect x="141.7" y="81.67" width="8.633" height="8.633" fill="#3582be"><title>CP_severity 1.435, eGFR 92.5: 0.702</title></rect>
<rect x="150" y="81.67" width="8.633" height="8.633" fill="#3e89c2"><title>CP_severity 1.565, eGFR 92.5: 0.671</title></rect>
<rect x="158.3" y="81.67" width="8.633" height="8.633" fill="#4791c6"><title>CP_severity 1.696, eGFR 92.5: 0.639</title></rect>
<rect x="166.7" y="81.67" width="8.633" height="8.633" fill="#5198ca"><title>CP_severity 1.826, eGFR 92.5: 0.608</title></rect>
<rect x="175" y="81.67" width="8.633" height="8.633" fill="#5aa0ce"><title>CP_severity 1.957, eGFR 92.5: 0.576</title></rect>
<rect x="183.3" y="81.67" width="8.633" height="8.633" fill="#64a8d3"><title>CP_severity 2.087, eGFR 92.5: 0.542</title></rect>
<rect x="191.7" y="81.67" width="8.633" height="8.633" fill="#6fb0d7"><title>CP_severity 2.217, eGFR 92.5: 0.506</title></rect>
<rect x="200" y="81.67" width="8.633" height="8.633" fill="#7cb6db"><title>CP_severity 2.348, eGFR 92.5: 0.47</title></rect>
<rect x="208.3" y="81.67" width="8.633" height="8.633" fill="#89bdde"><title>CP_severity 2.478, eGFR 92.5: 0.434</title></rect>
<rect x="216.7" y="81.67" width="8.633" height="8.633" fill="#95c3e2"><title>CP_severity 2.609, eGFR 92.5: 0.398</title></rect>
<rect x="225" y="81.67" width="8.633" height="8.633" fill="#a2c9e5"><title>CP_severity 2.739, eGFR 92.5: 0.362</title></rect>
<rect x="233.3" y="81.67" width="8.633" height="8.633" fill="#afd0e9"><title>CP_severity 2.87, eGFR 92.5: 0.327</title></rect>
<rect x="241.7" y="81.67" width="8.633" height="8.633" fill="#bcd6ec"><title>CP_severity 3, eGFR 92.5: 0.291</title></rect>
<rect x="50" y="90" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 0, eGFR 87.5: 0.764</title></rect>
<rect x="58.33" y="90" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 0.13, eGFR 87.5: 0.764</title></rect>
<rect x="66.67" y="90" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 0.261, eGFR 87.5: 0.764</title></rect>
<rect x="75" y="90" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 0.391, eGFR 87.5: 0.764</title></rect>
<rect x="83.33" y="90" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 0.522, eGFR 87.5: 0.764</title></rect>
<rect x="91.67" y="90" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 0.652, eGFR 87.5: 0.764</title></rect>
<rect x="100" y="90" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 0.783, eGFR 87.5: 0.764</title></rect>
<rect x="108.3" y="90" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 0.913, eGFR 87.5: 0.764</title></rect>
<rect x="116.7" y="90" width="8.633" height="8.633" fill="#2675b7"><title>CP_severity 1.043, eGFR 87.5: 0.754</title></rect>
<rect x="125" y="90" width="8.633" height="8.633" fill="#2f7cbb"><title>CP_severity 1.174, eGFR 87.5: 0.724</title></rect>
<rect x="133.3" y="90" width="8.633" height="8.633" fill="#3783bf"><title>CP_severity 1.304, eGFR 87.5: 0.694</title></rect>
<rect x="141.7" y="90" width="8.633" height="8.633" fill="#408bc3"><title>CP_severity 1.435, eGFR 87.5: 0.665</title></rect>
<rect x="150" y="90" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 1.565, eGFR 87.5: 0.635</title></rect>
<rect x="158.3" y="90" width="8.633" height="8.633" fill="#5299cb"><title>CP_severity 1.696, eGFR 87.5: 0.605</title></rect>
<rect x="166.7" y="90" width="8.633" height="8.633" fill="#5aa0cf"><title>CP_severity 1.826, eGFR 87.5: 0.575</title></rect>
<rect x="175" y="90" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 1.957, eGFR 87.5: 0.545</title></rect>
<rect x="183.3" y="90" width="8.633" height="8.633" fill="#6dafd6"><title>CP_severity 2.087, eGFR 87.5: 0.512</title></rect>
<rect x="191.7" y="90" width="8.633" height="8.633" fill="#79b5da"><title>CP_severity 2.217, eGFR 87.5: 0.478</title></rect>
<rect x="200" y="90" width="8.633" height="8.633" fill="#85bbdd"><title>CP_severity 2.348, eGFR 87.5: 0.445</title></rect>
<rect x="208.3" y="90" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.478, eGFR 87.5: 0.411</title></rect>
<rect x="216.7" y="90" width="8.633" height="8.633" fill="#9dc7e4"><title>CP_severity 2.609, eGFR 87.5: 0.377</title></rect>
<rect x="225" y="90" width="8.633" height="8.633" fill="#a9cde7"><title>CP_severity 2.739, eGFR 87.5: 0.343</title></rect>
<rect x="233.3" y="90" width="8.633" height="8.633" fill="#b5d3ea"><title>CP_severity 2.87, eGFR 87.5: 0.309</title></rect>
<rect x="241.7" y="90" width="8.633" height="8.633" fill="#c1d9ee"><title>CP_severity 3, eGFR 87.5: 0.275</title></rect>
<rect x="50" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>CP_severity 0, eGFR 82.5: 0.721</title></rect>
<rect x="58.33" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>CP_severity 0.13, eGFR 82.5: 0.721</title></rect>
<rect x="66.67" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>CP_severity 0.261, eGFR 82.5: 0.721</title></rect>
<rect x="75" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>CP_severity 0.391, eGFR 82.5: 0.721</title></rect>
<rect x="83.33" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>CP_severity 0.522, eGFR 82.5: 0.721</title></rect>
<rect x="91.67" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>CP_severity 0.652, eGFR 82.5: 0.721</title></rect>
<rect x="100" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>CP_severity 0.783, eGFR 82.5: 0.721</title></rect>
<rect x="108.3" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>CP_severity 0.913, eGFR 82.5: 0.721</title></rect>
<rect x="116.7" y="98.33" width="8.633" height="8.633" fill="#337fbd"><title>CP_severity 1.043, eGFR 82.5: 0.711</title></rect>
<rect x="125" y="98.33" width="8.633" height="8.633" fill="#3b86c0"><title>CP_severity 1.174, eGFR 82.5: 0.683</title></rect>
<rect x="133.3" y="98.33" width="8.633" height="8.633" fill="#438dc4"><title>CP_severity 1.304, eGFR 82.5: 0.655</title></rect>
<rect x="141.7" y="98.33" width="8.633" height="8.633" fill="#4b94c8"><title>CP_severity 1.435, eGFR 82.5: 0.627</title></rect>
<rect x="150" y="98.33" width="8.633" height="8.633" fill="#539bcb"><title>CP_severity 1.565, eGFR 82.5: 0.598</title></rect>
<rect x="158.3" y="98.33" width="8.633" height="8.633" fill="#5ca1cf"><title>CP_severity 1.696, eGFR 82.5: 0.57</title></rect>
<rect x="166.7" y="98.33" width="8.633" height="8.633" fill="#64a8d3"><title>CP_severity 1.826, eGFR 82.5: 0.542</title></rect>
<rect x="175" y="98.33" width="8.633" height="8.633" fill="#6cafd6"><title>CP_severity 1.957, eGFR 82.5: 0.514</title></rect>
<rect x="183.3" y="98.33" width="8.633" height="8.633" fill="#77b4d9"><title>CP_severity 2.087, eGFR 82.5: 0.483</title></rect>
<rect x="191.7" y="98.33" width="8.633" height="8.633" fill="#83badc"><title>CP_severity 2.217, eGFR 82.5: 0.451</title></rect>
<rect x="200" y="98.33" width="8.633" height="8.633" fill="#8ebfe0"><title>CP_severity 2.348, eGFR 82.5: 0.419</title></rect>
<rect x="208.3" y="98.33" width="8.633" height="8.633" fill="#99c5e3"><title>CP_severity 2.478, eGFR 82.5: 0.387</title></rect>
<rect x="216.7" y="98.33" width="8.633" height="8.633" fill="#a5cbe6"><title>CP_severity 2.609, eGFR 82.5: 0.355</title></rect>
<rect x="225" y="98.33" width="8.633" height="8.633" fill="#b0d0e9"><title>CP_severity 2.739, eGFR 82.5: 0.323</title></rect>
<rect x="233.3" y="98.33" width="8.633" height="8.633" fill="#bcd6ec"><title>CP_severity 2.87, eGFR 82.5: 0.291</title></rect>
<rect x="241.7" y="98.33" width="8.633" height="8.633" fill="#c7dbef"><title>CP_severity 3, eGFR 82.5: 0.259</title></rect>
<rect x="50" y="106.7" width="8.633" height="8.633" fill="#3d88c1"><title>CP_severity 0, eGFR 77.5: 0.677</title></rect>
<rect x="58.33" y="106.7" width="8.633" height="8.633" fill="#3d88c1"><title>CP_severity 0.13, eGFR 77.5: 0.677</title></rect>
<rect x="66.67" y="106.7" width="8.633" height="8.633" fill="#3d88c1"><title>CP_severity 0.261, eGFR 77.5: 0.677</title></rect>
<rect x="75" y="106.7" width="8.633" height="8.633" fill="#3d88c1"><title>CP_severity 0.391, eGFR 77.5: 0.677</title></rect>
<rect x="83.33" y="106.7" width="8.633" height="8.633" fill="#3d88c1"><title>CP_severity 0.522, eGFR 77.5: 0.677</title></rect>
<rect x="91.67" y="106.7" width="8.633" height="8.633" fill="#3d88c1"><title>CP_severity 0.652, eGFR 77.5: 0.677</title></rect>
<rect x="100" y="106.7" width="8.633" height="8.633" fill="#3d88c1"><title>CP_severity 0.783, eGFR 77.5: 0.677</title></rect>
<rect x="108.3" y="106.7" width="8.633" height="8.633" fill="#3d88c1"><title>CP_severity 0.913, eGFR 77.5: 0.677</title></rect>
<rect x="116.7" y="106.7" width="8.633" height="8.633" fill="#3f8ac2"><title>CP_severity 1.043, eGFR 77.5: 0.668</title></rect>
<rect x="125" y="106.7" width="8.633" height="8.633" fill="#4790c6"><title>CP_severity 1.174, eGFR 77.5: 0.642</title></rect>
<rect x="133.3" y="106.7" width="8.633" height="8.633" fill="#4f97c9"><title>CP_severity 1.304, eGFR 77.5: 0.615</title></rect>
<rect x="141.7" y="106.7" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 1.435, eGFR 77.5: 0.589</title></rect>
<rect x="150" y="106.7" width="8.633" height="8.633" fill="#5ea3d0"><title>CP_severity 1.565, eGFR 77.5: 0.562</title></rect>
<rect x="158.3" y="106.7" width="8.633" height="8.633" fill="#66aad4"><title>CP_severity 1.696, eGFR 77.5: 0.536</title></rect>
<rect x="166.7" y="106.7" width="8.633" height="8.633" fill="#6eafd7"><title>CP_severity 1.826, eGFR 77.5: 0.509</title></rect>
<rect x="175" y="106.7" width="8.633" height="8.633" fill="#77b4d9"><title>CP_severity 1.957, eGFR 77.5: 0.483</title></rect>
<rect x="183.3" y="106.7" width="8.633" height="8.633" fill="#82b9dc"><title>CP_severity 2.087, eGFR 77.5: 0.454</title></rect>
<rect x="191.7" y="106.7" width="8.633" height="8.633" fill="#8cbedf"><title>CP_severity 2.217, eGFR 77.5: 0.424</title></rect>
<rect x="200" y="106.7" width="8.633" height="8.633" fill="#97c4e2"><title>CP_severity 2.348, eGFR 77.5: 0.394</title></rect>
<rect x="208.3" y="106.7" width="8.633" height="8.633" fill="#a2c9e5"><title>CP_severity 2.478, eGFR 77.5: 0.364</title></rect>
<rect x="216.7" y="106.7" width="8.633" height="8.633" fill="#adcee8"><title>CP_severity 2.609, eGFR 77.5: 0.334</title></rect>
<rect x="225" y="106.7" width="8.633" height="8.633" fill="#b7d4eb"><title>CP_severity 2.739, eGFR 77.5: 0.304</title></rect>
<rect x="233.3" y="106.7" width="8.633" height="8.633" fill="#c2d9ee"><title>CP_severity 2.87, eGFR 77.5: 0.274</title></rect>
<rect x="241.7" y="106.7" width="8.633" height="8.633" fill="#caddf0"><title>CP_severity 3, eGFR 77.5: 0.244</title></rect>
<rect x="50" y="115" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 0, eGFR 72.5: 0.633</title></rect>
<rect x="58.33" y="115" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 0.13, eGFR 72.5: 0.633</title></rect>
<rect x="66.67" y="115" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 0.261, eGFR 72.5: 0.633</title></rect>
<rect x="75" y="115" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 0.391, eGFR 72.5: 0.633</title></rect>
<rect x="83.33" y="115" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 0.522, eGFR 72.5: 0.633</title></rect>
<rect x="91.67" y="115" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 0.652, eGFR 72.5: 0.633</title></rect>
<rect x="100" y="115" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 0.783, eGFR 72.5: 0.633</title></rect>
<rect x="108.3" y="115" width="8.633" height="8.633" fill="#4992c7"><title>CP_severity 0.913, eGFR 72.5: 0.633</title></rect>
<rect x="116.7" y="115" width="8.633" height="8.633" fill="#4c94c8"><title>CP_severity 1.043, eGFR 72.5: 0.625</title></rect>
<rect x="125" y="115" width="8.633" height="8.633" fill="#539acb"><title>CP_severity 1.174, eGFR 72.5: 0.6</title></rect>
<rect x="133.3" y="115" width="8.633" height="8.633" fill="#5aa0ce"><title>CP_severity 1.304, eGFR 72.5: 0.575</title></rect>
<rect x="141.7" y="115" width="8.633" height="8.633" fill="#61a6d2"><title>CP_severity 1.435, eGFR 72.5: 0.551</title></rect>
<rect x="150" y="115" width="8.633" height="8.633" fill="#68acd5"><title>CP_severity 1.565, eGFR 72.5: 0.526</title></rect>
<rect x="158.3" y="115" width="8.633" height="8.633" fill="#71b1d8"><title>CP_severity 1.696, eGFR 72.5: 0.501</title></rect>
<rect x="166.7" y="115" width="8.633" height="8.633" fill="#7ab5da"><title>CP_severity 1.826, eGFR 72.5: 0.476</title></rect>
<rect x="175" y="115" width="8.633" height="8.633" fill="#82badc"><title>CP_severity 1.957, eGFR 72.5: 0.451</title></rect>
<rect x="183.3" y="115" width="8.633" height="8.633" fill="#8cbedf"><title>CP_severity 2.087, eGFR 72.5: 0.425</title></rect>
<rect x="191.7" y="115" width="8.633" height="8.633" fill="#96c3e2"><title>CP_severity 2.217, eGFR 72.5: 0.396</title></rect>
<rect x="200" y="115" width="8.633" height="8.633" fill="#a0c8e5"><title>CP_severity 2.348, eGFR 72.5: 0.368</title></rect>
<rect x="208.3" y="115" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.478, eGFR 72.5: 0.34</title></rect>
<rect x="216.7" y="115" width="8.633" height="8.633" fill="#b4d2ea"><title>CP_severity 2.609, eGFR 72.5: 0.312</title></rect>
<rect x="225" y="115" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 2.739, eGFR 72.5: 0.284</title></rect>
<rect x="233.3" y="115" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 2.87, eGFR 72.5: 0.256</title></rect>
<rect x="241.7" y="115" width="8.633" height="8.633" fill="#cddff1"><title>CP_severity 3, eGFR 72.5: 0.228</title></rect>
<rect x="50" y="123.3" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 0, eGFR 67.5: 0.59</title></rect>
<rect x="58.33" y="123.3" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 0.13, eGFR 67.5: 0.59</title></rect>
<rect x="66.67" y="123.3" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 0.261, eGFR 67.5: 0.59</title></rect>
<rect x="75" y="123.3" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 0.391, eGFR 67.5: 0.59</title></rect>
<rect x="83.33" y="123.3" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 0.522, eGFR 67.5: 0.59</title></rect>
<rect x="91.67" y="123.3" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 0.652, eGFR 67.5: 0.59</title></rect>
<rect x="100" y="123.3" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 0.783, eGFR 67.5: 0.59</title></rect>
<rect x="108.3" y="123.3" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 0.913, eGFR 67.5: 0.59</title></rect>
<rect x="116.7" y="123.3" width="8.633" height="8.633" fill="#589ece"><title>CP_severity 1.043, eGFR 67.5: 0.582</title></rect>
<rect x="125" y="123.3" width="8.633" height="8.633" fill="#5fa4d1"><title>CP_severity 1.174, eGFR 67.5: 0.559</title></rect>
<rect x="133.3" y="123.3" width="8.633" height="8.633" fill="#66aad4"><title>CP_severity 1.304, eGFR 67.5: 0.536</title></rect>
<rect x="141.7" y="123.3" width="8.633" height="8.633" fill="#6dafd6"><title>CP_severity 1.435, eGFR 67.5: 0.513</title></rect>
<rect x="150" y="123.3" width="8.633" height="8.633" fill="#75b3d9"><title>CP_severity 1.565, eGFR 67.5: 0.49</title></rect>
<rect x="158.3" y="123.3" width="8.633" height="8.633" fill="#7db7db"><title>CP_severity 1.696, eGFR 67.5: 0.466</title></rect>
<rect x="166.7" y="123.3" width="8.633" height="8.633" fill="#85bbdd"><title>CP_severity 1.826, eGFR 67.5: 0.443</title></rect>
<rect x="175" y="123.3" width="8.633" height="8.633" fill="#8ebfdf"><title>CP_severity 1.957, eGFR 67.5: 0.42</title></rect>
<rect x="183.3" y="123.3" width="8.633" height="8.633" fill="#97c4e2"><title>CP_severity 2.087, eGFR 67.5: 0.395</title></rect>
<rect x="191.7" y="123.3" width="8.633" height="8.633" fill="#a0c8e5"><title>CP_severity 2.217, eGFR 67.5: 0.369</title></rect>
<rect x="200" y="123.3" width="8.633" height="8.633" fill="#a9cde7"><title>CP_severity 2.348, eGFR 67.5: 0.343</title></rect>
<rect x="208.3" y="123.3" width="8.633" height="8.633" fill="#b3d1ea"><title>CP_severity 2.478, eGFR 67.5: 0.317</title></rect>
<rect x="216.7" y="123.3" width="8.633" height="8.633" fill="#bcd6ec"><title>CP_severity 2.609, eGFR 67.5: 0.291</title></rect>
<rect x="225" y="123.3" width="8.633" height="8.633" fill="#c5dbef"><title>CP_severity 2.739, eGFR 67.5: 0.265</title></rect>
<rect x="233.3" y="123.3" width="8.633" height="8.633" fill="#cbdef1"><title>CP_severity 2.87, eGFR 67.5: 0.238</title></rect>
<rect x="241.7" y="123.3" width="8.633" height="8.633" fill="#d0e1f2"><title>CP_severity 3, eGFR 67.5: 0.212</title></rect>
<rect x="50" y="131.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 0, eGFR 62.5: 0.546</title></rect>
<rect x="58.33" y="131.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 0.13, eGFR 62.5: 0.546</title></rect>
<rect x="66.67" y="131.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 0.261, eGFR 62.5: 0.546</title></rect>
<rect x="75" y="131.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 0.391, eGFR 62.5: 0.546</title></rect>
<rect x="83.33" y="131.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 0.522, eGFR 62.5: 0.546</title></rect>
<rect x="91.67" y="131.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 0.652, eGFR 62.5: 0.546</title></rect>
<rect x="100" y="131.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 0.783, eGFR 62.5: 0.546</title></rect>
<rect x="108.3" y="131.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 0.913, eGFR 62.5: 0.546</title></rect>
<rect x="116.7" y="131.7" width="8.633" height="8.633" fill="#65a9d3"><title>CP_severity 1.043, eGFR 62.5: 0.539</title></rect>
<rect x="125" y="131.7" width="8.633" height="8.633" fill="#6baed6"><title>CP_severity 1.174, eGFR 62.5: 0.517</title></rect>
<rect x="133.3" y="131.7" width="8.633" height="8.633" fill="#73b2d8"><title>CP_severity 1.304, eGFR 62.5: 0.496</title></rect>
<rect x="141.7" y="131.7" width="8.633" height="8.633" fill="#7ab5da"><title>CP_severity 1.435, eGFR 62.5: 0.475</title></rect>
<rect x="150" y="131.7" width="8.633" height="8.633" fill="#82b9dc"><title>CP_severity 1.565, eGFR 62.5: 0.453</title></rect>
<rect x="158.3" y="131.7" width="8.633" height="8.633" fill="#89bdde"><title>CP_severity 1.696, eGFR 62.5: 0.432</title></rect>
<rect x="166.7" y="131.7" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 1.826, eGFR 62.5: 0.411</title></rect>
<rect x="175" y="131.7" width="8.633" height="8.633" fill="#99c5e3"><title>CP_severity 1.957, eGFR 62.5: 0.389</title></rect>
<rect x="183.3" y="131.7" width="8.633" height="8.633" fill="#a1c9e5"><title>CP_severity 2.087, eGFR 62.5: 0.366</title></rect>
<rect x="191.7" y="131.7" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.217, eGFR 62.5: 0.342</title></rect>
<rect x="200" y="131.7" width="8.633" height="8.633" fill="#b2d1ea"><title>CP_severity 2.348, eGFR 62.5: 0.318</title></rect>
<rect x="208.3" y="131.7" width="8.633" height="8.633" fill="#bbd6ec"><title>CP_severity 2.478, eGFR 62.5: 0.293</title></rect>
<rect x="216.7" y="131.7" width="8.633" height="8.633" fill="#c4daee"><title>CP_severity 2.609, eGFR 62.5: 0.269</title></rect>
<rect x="225" y="131.7" width="8.633" height="8.633" fill="#c9ddf0"><title>CP_severity 2.739, eGFR 62.5: 0.245</title></rect>
<rect x="233.3" y="131.7" width="8.633" height="8.633" fill="#cee0f2"><title>CP_severity 2.87, eGFR 62.5: 0.221</title></rect>
<rect x="241.7" y="131.7" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 3, eGFR 62.5: 0.197</title></rect>
<rect x="50" y="140" width="8.633" height="8.633" fill="#70b1d7"><title>CP_severity 0, eGFR 57.5: 0.502</title></rect>
<rect x="58.33" y="140" width="8.633" height="8.633" fill="#70b1d7"><title>CP_severity 0.13, eGFR 57.5: 0.502</title></rect>
<rect x="66.67" y="140" width="8.633" height="8.633" fill="#70b1d7"><title>CP_severity 0.261, eGFR 57.5: 0.502</title></rect>
<rect x="75" y="140" width="8.633" height="8.633" fill="#70b1d7"><title>CP_severity 0.391, eGFR 57.5: 0.502</title></rect>
<rect x="83.33" y="140" width="8.633" height="8.633" fill="#70b1d7"><title>CP_severity 0.522, eGFR 57.5: 0.502</title></rect>
<rect x="91.67" y="140" width="8.633" height="8.633" fill="#70b1d7"><title>CP_severity 0.652, eGFR 57.5: 0.502</title></rect>
<rect x="100" y="140" width="8.633" height="8.633" fill="#70b1d7"><title>CP_severity 0.783, eGFR 57.5: 0.502</title></rect>
<rect x="108.3" y="140" width="8.633" height="8.633" fill="#70b1d7"><title>CP_severity 0.913, eGFR 57.5: 0.502</title></rect>
<rect x="116.7" y="140" width="8.633" height="8.633" fill="#73b2d8"><title>CP_severity 1.043, eGFR 57.5: 0.496</title></rect>
<rect x="125" y="140" width="8.633" height="8.633" fill="#7ab5da"><title>CP_severity 1.174, eGFR 57.5: 0.476</title></rect>
<rect x="133.3" y="140" width="8.633" height="8.633" fill="#81b9dc"><title>CP_severity 1.304, eGFR 57.5: 0.456</title></rect>
<rect x="141.7" y="140" width="8.633" height="8.633" fill="#88bcde"><title>CP_severity 1.435, eGFR 57.5: 0.437</title></rect>
<rect x="150" y="140" width="8.633" height="8.633" fill="#8fc0e0"><title>CP_severity 1.565, eGFR 57.5: 0.417</title></rect>
<rect x="158.3" y="140" width="8.633" height="8.633" fill="#96c3e2"><title>CP_severity 1.696, eGFR 57.5: 0.397</title></rect>
<rect x="166.7" y="140" width="8.633" height="8.633" fill="#9dc7e4"><title>CP_severity 1.826, eGFR 57.5: 0.378</title></rect>
<rect x="175" y="140" width="8.633" height="8.633" fill="#a4cae6"><title>CP_severity 1.957, eGFR 57.5: 0.358</title></rect>
<rect x="183.3" y="140" width="8.633" height="8.633" fill="#abcee8"><title>CP_severity 2.087, eGFR 57.5: 0.337</title></rect>
<rect x="191.7" y="140" width="8.633" height="8.633" fill="#b3d2ea"><title>CP_severity 2.217, eGFR 57.5: 0.314</title></rect>
<rect x="200" y="140" width="8.633" height="8.633" fill="#bbd6ec"><title>CP_severity 2.348, eGFR 57.5: 0.292</title></rect>
<rect x="208.3" y="140" width="8.633" height="8.633" fill="#c3daee"><title>CP_severity 2.478, eGFR 57.5: 0.27</title></rect>
<rect x="216.7" y="140" width="8.633" height="8.633" fill="#c9ddf0"><title>CP_severity 2.609, eGFR 57.5: 0.248</title></rect>
<rect x="225" y="140" width="8.633" height="8.633" fill="#cde0f1"><title>CP_severity 2.739, eGFR 57.5: 0.225</title></rect>
<rect x="233.3" y="140" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 2.87, eGFR 57.5: 0.203</title></rect>
<rect x="241.7" y="140" width="8.633" height="8.633" fill="#d6e5f4"><title>CP_severity 3, eGFR 57.5: 0.181</title></rect>
<rect x="50" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 0, eGFR 52.5: 0.459</title></rect>
<rect x="58.33" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 0.13, eGFR 52.5: 0.459</title></rect>
<rect x="66.67" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 0.261, eGFR 52.5: 0.459</title></rect>
<rect x="75" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 0.391, eGFR 52.5: 0.459</title></rect>
<rect x="83.33" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 0.522, eGFR 52.5: 0.459</title></rect>
<rect x="91.67" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 0.652, eGFR 52.5: 0.459</title></rect>
<rect x="100" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 0.783, eGFR 52.5: 0.459</title></rect>
<rect x="108.3" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 0.913, eGFR 52.5: 0.459</title></rect>
<rect x="116.7" y="148.3" width="8.633" height="8.633" fill="#82b9dc"><title>CP_severity 1.043, eGFR 52.5: 0.453</title></rect>
<rect x="125" y="148.3" width="8.633" height="8.633" fill="#88bdde"><title>CP_severity 1.174, eGFR 52.5: 0.435</title></rect>
<rect x="133.3" y="148.3" width="8.633" height="8.633" fill="#8fc0e0"><title>CP_severity 1.304, eGFR 52.5: 0.417</title></rect>
<rect x="141.7" y="148.3" width="8.633" height="8.633" fill="#95c3e2"><title>CP_severity 1.435, eGFR 52.5: 0.399</title></rect>
<rect x="150" y="148.3" width="8.633" height="8.633" fill="#9cc6e3"><title>CP_severity 1.565, eGFR 52.5: 0.381</title></rect>
<rect x="158.3" y="148.3" width="8.633" height="8.633" fill="#a2c9e5"><title>CP_severity 1.696, eGFR 52.5: 0.363</title></rect>
<rect x="166.7" y="148.3" width="8.633" height="8.633" fill="#a9cce7"><title>CP_severity 1.826, eGFR 52.5: 0.345</title></rect>
<rect x="175" y="148.3" width="8.633" height="8.633" fill="#afd0e9"><title>CP_severity 1.957, eGFR 52.5: 0.327</title></rect>
<rect x="183.3" y="148.3" width="8.633" height="8.633" fill="#b6d3eb"><title>CP_severity 2.087, eGFR 52.5: 0.307</title></rect>
<rect x="191.7" y="148.3" width="8.633" height="8.633" fill="#bdd7ed"><title>CP_severity 2.217, eGFR 52.5: 0.287</title></rect>
<rect x="200" y="148.3" width="8.633" height="8.633" fill="#c4daef"><title>CP_severity 2.348, eGFR 52.5: 0.267</title></rect>
<rect x="208.3" y="148.3" width="8.633" height="8.633" fill="#c9ddf0"><title>CP_severity 2.478, eGFR 52.5: 0.246</title></rect>
<rect x="216.7" y="148.3" width="8.633" height="8.633" fill="#cde0f1"><title>CP_severity 2.609, eGFR 52.5: 0.226</title></rect>
<rect x="225" y="148.3" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 2.739, eGFR 52.5: 0.206</title></rect>
<rect x="233.3" y="148.3" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 2.87, eGFR 52.5: 0.185</title></rect>
<rect x="241.7" y="148.3" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 3, eGFR 52.5: 0.165</title></rect>
<rect x="50" y="156.7" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 0, eGFR 47.5: 0.415</title></rect>
<rect x="58.33" y="156.7" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 0.13, eGFR 47.5: 0.415</title></rect>
<rect x="66.67" y="156.7" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 0.261, eGFR 47.5: 0.415</title></rect>
<rect x="75" y="156.7" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 0.391, eGFR 47.5: 0.415</title></rect>
<rect x="83.33" y="156.7" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 0.522, eGFR 47.5: 0.415</title></rect>
<rect x="91.67" y="156.7" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 0.652, eGFR 47.5: 0.415</title></rect>
<rect x="100" y="156.7" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 0.783, eGFR 47.5: 0.415</title></rect>
<rect x="108.3" y="156.7" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 0.913, eGFR 47.5: 0.415</title></rect>
<rect x="116.7" y="156.7" width="8.633" height="8.633" fill="#91c1e1"><title>CP_severity 1.043, eGFR 47.5: 0.409</title></rect>
<rect x="125" y="156.7" width="8.633" height="8.633" fill="#97c4e2"><title>CP_severity 1.174, eGFR 47.5: 0.393</title></rect>
<rect x="133.3" y="156.7" width="8.633" height="8.633" fill="#9dc7e4"><title>CP_severity 1.304, eGFR 47.5: 0.377</title></rect>
<rect x="141.7" y="156.7" width="8.633" height="8.633" fill="#a3cae5"><title>CP_severity 1.435, eGFR 47.5: 0.361</title></rect>
<rect x="150" y="156.7" width="8.633" height="8.633" fill="#a9cce7"><title>CP_severity 1.565, eGFR 47.5: 0.345</title></rect>
<rect x="158.3" y="156.7" width="8.633" height="8.633" fill="#aecfe9"><title>CP_severity 1.696, eGFR 47.5: 0.328</title></rect>
<rect x="166.7" y="156.7" width="8.633" height="8.633" fill="#b4d2ea"><title>CP_severity 1.826, eGFR 47.5: 0.312</title></rect>
<rect x="175" y="156.7" width="8.633" height="8.633" fill="#bad5ec"><title>CP_severity 1.957, eGFR 47.5: 0.296</title></rect>
<rect x="183.3" y="156.7" width="8.633" height="8.633" fill="#c0d8ed"><title>CP_severity 2.087, eGFR 47.5: 0.278</title></rect>
<rect x="191.7" y="156.7" width="8.633" height="8.633" fill="#c7dbef"><title>CP_severity 2.217, eGFR 47.5: 0.26</title></rect>
<rect x="200" y="156.7" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 2.348, eGFR 47.5: 0.241</title></rect>
<rect x="208.3" y="156.7" width="8.633" height="8.633" fill="#cee0f1"><title>CP_severity 2.478, eGFR 47.5: 0.223</title></rect>
<rect x="216.7" y="156.7" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 2.609, eGFR 47.5: 0.205</title></rect>
<rect x="225" y="156.7" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 2.739, eGFR 47.5: 0.186</title></rect>
<rect x="233.3" y="156.7" width="8.633" height="8.633" fill="#d8e7f5"><title>CP_severity 2.87, eGFR 47.5: 0.168</title></rect>
<rect x="241.7" y="156.7" width="8.633" height="8.633" fill="#dce9f6"><title>CP_severity 3, eGFR 47.5: 0.149</title></rect>
<rect x="50" y="165" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0, eGFR 42.5: 0.371</title></rect>
<rect x="58.33" y="165" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.13, eGFR 42.5: 0.371</title></rect>
<rect x="66.67" y="165" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.261, eGFR 42.5: 0.371</title></rect>
<rect x="75" y="165" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.391, eGFR 42.5: 0.371</title></rect>
<rect x="83.33" y="165" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.522, eGFR 42.5: 0.371</title></rect>
<rect x="91.67" y="165" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.652, eGFR 42.5: 0.371</title></rect>
<rect x="100" y="165" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.783, eGFR 42.5: 0.371</title></rect>
<rect x="108.3" y="165" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.913, eGFR 42.5: 0.371</title></rect>
<rect x="116.7" y="165" width="8.633" height="8.633" fill="#a1c9e5"><title>CP_severity 1.043, eGFR 42.5: 0.366</title></rect>
<rect x="125" y="165" width="8.633" height="8.633" fill="#a6cbe6"><title>CP_severity 1.174, eGFR 42.5: 0.352</title></rect>
<rect x="133.3" y="165" width="8.633" height="8.633" fill="#abcee8"><title>CP_severity 1.304, eGFR 42.5: 0.337</title></rect>
<rect x="141.7" y="165" width="8.633" height="8.633" fill="#b0d0e9"><title>CP_severity 1.435, eGFR 42.5: 0.323</title></rect>
<rect x="150" y="165" width="8.633" height="8.633" fill="#b6d3eb"><title>CP_severity 1.565, eGFR 42.5: 0.308</title></rect>
<rect x="158.3" y="165" width="8.633" height="8.633" fill="#bbd5ec"><title>CP_severity 1.696, eGFR 42.5: 0.294</title></rect>
<rect x="166.7" y="165" width="8.633" height="8.633" fill="#c0d8ed"><title>CP_severity 1.826, eGFR 42.5: 0.279</title></rect>
<rect x="175" y="165" width="8.633" height="8.633" fill="#c5dbef"><title>CP_severity 1.957, eGFR 42.5: 0.265</title></rect>
<rect x="183.3" y="165" width="8.633" height="8.633" fill="#c9ddf0"><title>CP_severity 2.087, eGFR 42.5: 0.249</title></rect>
<rect x="191.7" y="165" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 2.217, eGFR 42.5: 0.232</title></rect>
<rect x="200" y="165" width="8.633" height="8.633" fill="#cfe1f2"><title>CP_severity 2.348, eGFR 42.5: 0.216</title></rect>
<rect x="208.3" y="165" width="8.633" height="8.633" fill="#d2e3f3"><title>CP_severity 2.478, eGFR 42.5: 0.199</title></rect>
<rect x="216.7" y="165" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 2.609, eGFR 42.5: 0.183</title></rect>
<rect x="225" y="165" width="8.633" height="8.633" fill="#d8e7f5"><title>CP_severity 2.739, eGFR 42.5: 0.167</title></rect>
<rect x="233.3" y="165" width="8.633" height="8.633" fill="#dce9f6"><title>CP_severity 2.87, eGFR 42.5: 0.15</title></rect>
<rect x="241.7" y="165" width="8.633" height="8.633" fill="#dfebf7"><title>CP_severity 3, eGFR 42.5: 0.134</title></rect>
<rect x="50" y="173.3" width="8.633" height="8.633" fill="#afcfe9"><title>CP_severity 0, eGFR 37.5: 0.328</title></rect>
<rect x="58.33" y="173.3" width="8.633" height="8.633" fill="#afcfe9"><title>CP_severity 0.13, eGFR 37.5: 0.328</title></rect>
<rect x="66.67" y="173.3" width="8.633" height="8.633" fill="#afcfe9"><title>CP_severity 0.261, eGFR 37.5: 0.328</title></rect>
<rect x="75" y="173.3" width="8.633" height="8.633" fill="#afcfe9"><title>CP_severity 0.391, eGFR 37.5: 0.328</title></rect>
<rect x="83.33" y="173.3" width="8.633" height="8.633" fill="#afcfe9"><title>CP_severity 0.522, eGFR 37.5: 0.328</title></rect>
<rect x="91.67" y="173.3" width="8.633" height="8.633" fill="#afcfe9"><title>CP_severity 0.652, eGFR 37.5: 0.328</title></rect>
<rect x="100" y="173.3" width="8.633" height="8.633" fill="#afcfe9"><title>CP_severity 0.783, eGFR 37.5: 0.328</title></rect>
<rect x="108.3" y="173.3" width="8.633" height="8.633" fill="#afcfe9"><title>CP_severity 0.913, eGFR 37.5: 0.328</title></rect>
<rect x="116.7" y="173.3" width="8.633" height="8.633" fill="#b0d0e9"><title>CP_severity 1.043, eGFR 37.5: 0.323</title></rect>
<rect x="125" y="173.3" width="8.633" height="8.633" fill="#b5d3ea"><title>CP_severity 1.174, eGFR 37.5: 0.31</title></rect>
<rect x="133.3" y="173.3" width="8.633" height="8.633" fill="#b9d5ec"><title>CP_severity 1.304, eGFR 37.5: 0.298</title></rect>
<rect x="141.7" y="173.3" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 1.435, eGFR 37.5: 0.285</title></rect>
<rect x="150" y="173.3" width="8.633" height="8.633" fill="#c3d9ee"><title>CP_severity 1.565, eGFR 37.5: 0.272</title></rect>
<rect x="158.3" y="173.3" width="8.633" height="8.633" fill="#c7dbef"><title>CP_severity 1.696, eGFR 37.5: 0.259</title></rect>
<rect x="166.7" y="173.3" width="8.633" height="8.633" fill="#c9ddf0"><title>CP_severity 1.826, eGFR 37.5: 0.246</title></rect>
<rect x="175" y="173.3" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 1.957, eGFR 37.5: 0.234</title></rect>
<rect x="183.3" y="173.3" width="8.633" height="8.633" fill="#cee0f2"><title>CP_severity 2.087, eGFR 37.5: 0.22</title></rect>
<rect x="191.7" y="173.3" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 2.217, eGFR 37.5: 0.205</title></rect>
<rect x="200" y="173.3" width="8.633" height="8.633" fill="#d4e4f4"><title>CP_severity 2.348, eGFR 37.5: 0.191</title></rect>
<rect x="208.3" y="173.3" width="8.633" height="8.633" fill="#d7e6f4"><title>CP_severity 2.478, eGFR 37.5: 0.176</title></rect>
<rect x="216.7" y="173.3" width="8.633" height="8.633" fill="#d9e8f5"><title>CP_severity 2.609, eGFR 37.5: 0.161</title></rect>
<rect x="225" y="173.3" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 2.739, eGFR 37.5: 0.147</title></rect>
<rect x="233.3" y="173.3" width="8.633" height="8.633" fill="#dfebf7"><title>CP_severity 2.87, eGFR 37.5: 0.132</title></rect>
<rect x="241.7" y="173.3" width="8.633" height="8.633" fill="#e2edf8"><title>CP_severity 3, eGFR 37.5: 0.118</title></rect>
<rect x="50" y="181.7" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 0, eGFR 32.5: 0.284</title></rect>
<rect x="58.33" y="181.7" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 0.13, eGFR 32.5: 0.284</title></rect>
<rect x="66.67" y="181.7" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 0.261, eGFR 32.5: 0.284</title></rect>
<rect x="75" y="181.7" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 0.391, eGFR 32.5: 0.284</title></rect>
<rect x="83.33" y="181.7" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 0.522, eGFR 32.5: 0.284</title></rect>
<rect x="91.67" y="181.7" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 0.652, eGFR 32.5: 0.284</title></rect>
<rect x="100" y="181.7" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 0.783, eGFR 32.5: 0.284</title></rect>
<rect x="108.3" y="181.7" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 0.913, eGFR 32.5: 0.284</title></rect>
<rect x="116.7" y="181.7" width="8.633" height="8.633" fill="#c0d8ed"><title>CP_severity 1.043, eGFR 32.5: 0.28</title></rect>
<rect x="125" y="181.7" width="8.633" height="8.633" fill="#c4daee"><title>CP_severity 1.174, eGFR 32.5: 0.269</title></rect>
<rect x="133.3" y="181.7" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 1.304, eGFR 32.5: 0.258</title></rect>
<rect x="141.7" y="181.7" width="8.633" height="8.633" fill="#c9ddf0"><title>CP_severity 1.435, eGFR 32.5: 0.247</title></rect>
<rect x="150" y="181.7" width="8.633" height="8.633" fill="#cbdef1"><title>CP_severity 1.565, eGFR 32.5: 0.236</title></rect>
<rect x="158.3" y="181.7" width="8.633" height="8.633" fill="#cde0f1"><title>CP_severity 1.696, eGFR 32.5: 0.225</title></rect>
<rect x="166.7" y="181.7" width="8.633" height="8.633" fill="#cfe1f2"><title>CP_severity 1.826, eGFR 32.5: 0.213</title></rect>
<rect x="175" y="181.7" width="8.633" height="8.633" fill="#d2e3f3"><title>CP_severity 1.957, eGFR 32.5: 0.202</title></rect>
<rect x="183.3" y="181.7" width="8.633" height="8.633" fill="#d4e4f4"><title>CP_severity 2.087, eGFR 32.5: 0.19</title></rect>
<rect x="191.7" y="181.7" width="8.633" height="8.633" fill="#d6e6f4"><title>CP_severity 2.217, eGFR 32.5: 0.178</title></rect>
<rect x="200" y="181.7" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 2.348, eGFR 32.5: 0.165</title></rect>
<rect x="208.3" y="181.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 2.478, eGFR 32.5: 0.153</title></rect>
<rect x="216.7" y="181.7" width="8.633" height="8.633" fill="#deeaf7"><title>CP_severity 2.609, eGFR 32.5: 0.14</title></rect>
<rect x="225" y="181.7" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 2.739, eGFR 32.5: 0.127</title></rect>
<rect x="233.3" y="181.7" width="8.633" height="8.633" fill="#e2eef8"><title>CP_severity 2.87, eGFR 32.5: 0.115</title></rect>
<rect x="241.7" y="181.7" width="8.633" height="8.633" fill="#e5eff9"><title>CP_severity 3, eGFR 32.5: 0.102</title></rect>
<rect x="50" y="190" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 0, eGFR 27.5: 0.24</title></rect>
<rect x="58.33" y="190" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 0.13, eGFR 27.5: 0.24</title></rect>
<rect x="66.67" y="190" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 0.261, eGFR 27.5: 0.24</title></rect>
<rect x="75" y="190" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 0.391, eGFR 27.5: 0.24</title></rect>
<rect x="83.33" y="190" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 0.522, eGFR 27.5: 0.24</title></rect>
<rect x="91.67" y="190" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 0.652, eGFR 27.5: 0.24</title></rect>
<rect x="100" y="190" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 0.783, eGFR 27.5: 0.24</title></rect>
<rect x="108.3" y="190" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 0.913, eGFR 27.5: 0.24</title></rect>
<rect x="116.7" y="190" width="8.633" height="8.633" fill="#cbdef1"><title>CP_severity 1.043, eGFR 27.5: 0.237</title></rect>
<rect x="125" y="190" width="8.633" height="8.633" fill="#cddff1"><title>CP_severity 1.174, eGFR 27.5: 0.228</title></rect>
<rect x="133.3" y="190" width="8.633" height="8.633" fill="#cfe1f2"><title>CP_severity 1.304, eGFR 27.5: 0.218</title></rect>
<rect x="141.7" y="190" width="8.633" height="8.633" fill="#d0e2f2"><title>CP_severity 1.435, eGFR 27.5: 0.209</title></rect>
<rect x="150" y="190" width="8.633" height="8.633" fill="#d2e3f3"><title>CP_severity 1.565, eGFR 27.5: 0.199</title></rect>
<rect x="158.3" y="190" width="8.633" height="8.633" fill="#d4e4f4"><title>CP_severity 1.696, eGFR 27.5: 0.19</title></rect>
<rect x="166.7" y="190" width="8.633" height="8.633" fill="#d6e5f4"><title>CP_severity 1.826, eGFR 27.5: 0.181</title></rect>
<rect x="175" y="190" width="8.633" height="8.633" fill="#d8e6f5"><title>CP_severity 1.957, eGFR 27.5: 0.171</title></rect>
<rect x="183.3" y="190" width="8.633" height="8.633" fill="#dae8f5"><title>CP_severity 2.087, eGFR 27.5: 0.161</title></rect>
<rect x="191.7" y="190" width="8.633" height="8.633" fill="#dce9f6"><title>CP_severity 2.217, eGFR 27.5: 0.15</title></rect>
<rect x="200" y="190" width="8.633" height="8.633" fill="#deeaf7"><title>CP_severity 2.348, eGFR 27.5: 0.14</title></rect>
<rect x="208.3" y="190" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 2.478, eGFR 27.5: 0.129</title></rect>
<rect x="216.7" y="190" width="8.633" height="8.633" fill="#e2edf8"><title>CP_severity 2.609, eGFR 27.5: 0.118</title></rect>
<rect x="225" y="190" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 2.739, eGFR 27.5: 0.108</title></rect>
<rect x="233.3" y="190" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 2.87, eGFR 27.5: 0.0971</title></rect>
<rect x="241.7" y="190" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 3, eGFR 27.5: 0.0865</title></rect>
<rect x="50" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0, eGFR 22.5: 0.197</title></rect>
<rect x="58.33" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0.13, eGFR 22.5: 0.197</title></rect>
<rect x="66.67" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0.261, eGFR 22.5: 0.197</title></rect>
<rect x="75" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0.391, eGFR 22.5: 0.197</title></rect>
<rect x="83.33" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0.522, eGFR 22.5: 0.197</title></rect>
<rect x="91.67" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0.652, eGFR 22.5: 0.197</title></rect>
<rect x="100" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0.783, eGFR 22.5: 0.197</title></rect>
<rect x="108.3" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0.913, eGFR 22.5: 0.197</title></rect>
<rect x="116.7" y="198.3" width="8.633" height="8.633" fill="#d3e4f3"><title>CP_severity 1.043, eGFR 22.5: 0.194</title></rect>
<rect x="125" y="198.3" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.174, eGFR 22.5: 0.186</title></rect>
<rect x="133.3" y="198.3" width="8.633" height="8.633" fill="#d6e6f4"><title>CP_severity 1.304, eGFR 22.5: 0.179</title></rect>
<rect x="141.7" y="198.3" width="8.633" height="8.633" fill="#d8e7f5"><title>CP_severity 1.435, eGFR 22.5: 0.171</title></rect>
<rect x="150" y="198.3" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.565, eGFR 22.5: 0.163</title></rect>
<rect x="158.3" y="198.3" width="8.633" height="8.633" fill="#dbe8f6"><title>CP_severity 1.696, eGFR 22.5: 0.155</title></rect>
<rect x="166.7" y="198.3" width="8.633" height="8.633" fill="#dce9f6"><title>CP_severity 1.826, eGFR 22.5: 0.148</title></rect>
<rect x="175" y="198.3" width="8.633" height="8.633" fill="#deeaf7"><title>CP_severity 1.957, eGFR 22.5: 0.14</title></rect>
<rect x="183.3" y="198.3" width="8.633" height="8.633" fill="#dfebf7"><title>CP_severity 2.087, eGFR 22.5: 0.132</title></rect>
<rect x="191.7" y="198.3" width="8.633" height="8.633" fill="#e1edf8"><title>CP_severity 2.217, eGFR 22.5: 0.123</title></rect>
<rect x="200" y="198.3" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 2.348, eGFR 22.5: 0.114</title></rect>
<rect x="208.3" y="198.3" width="8.633" height="8.633" fill="#e4eff9"><title>CP_severity 2.478, eGFR 22.5: 0.106</title></rect>
<rect x="216.7" y="198.3" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 2.609, eGFR 22.5: 0.0969</title></rect>
<rect x="225" y="198.3" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 2.739, eGFR 22.5: 0.0882</title></rect>
<rect x="233.3" y="198.3" width="8.633" height="8.633" fill="#e9f2fb"><title>CP_severity 2.87, eGFR 22.5: 0.0795</title></rect>
<rect x="241.7" y="198.3" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 3, eGFR 22.5: 0.0707</title></rect>
<rect x="50" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0, eGFR 17.5: 0.153</title></rect>
<rect x="58.33" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.13, eGFR 17.5: 0.153</title></rect>
<rect x="66.67" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.261, eGFR 17.5: 0.153</title></rect>
<rect x="75" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.391, eGFR 17.5: 0.153</title></rect>
<rect x="83.33" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.522, eGFR 17.5: 0.153</title></rect>
<rect x="91.67" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.652, eGFR 17.5: 0.153</title></rect>
<rect x="100" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.783, eGFR 17.5: 0.153</title></rect>
<rect x="108.3" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.913, eGFR 17.5: 0.153</title></rect>
<rect x="116.7" y="206.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 1.043, eGFR 17.5: 0.151</title></rect>
<rect x="125" y="206.7" width="8.633" height="8.633" fill="#ddeaf6"><title>CP_severity 1.174, eGFR 17.5: 0.145</title></rect>
<rect x="133.3" y="206.7" width="8.633" height="8.633" fill="#deebf7"><title>CP_severity 1.304, eGFR 17.5: 0.139</title></rect>
<rect x="141.7" y="206.7" width="8.633" height="8.633" fill="#dfebf7"><title>CP_severity 1.435, eGFR 17.5: 0.133</title></rect>
<rect x="150" y="206.7" width="8.633" height="8.633" fill="#e0ecf8"><title>CP_severity 1.565, eGFR 17.5: 0.127</title></rect>
<rect x="158.3" y="206.7" width="8.633" height="8.633" fill="#e1edf8"><title>CP_severity 1.696, eGFR 17.5: 0.121</title></rect>
<rect x="166.7" y="206.7" width="8.633" height="8.633" fill="#e2eef8"><title>CP_severity 1.826, eGFR 17.5: 0.115</title></rect>
<rect x="175" y="206.7" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 1.957, eGFR 17.5: 0.109</title></rect>
<rect x="183.3" y="206.7" width="8.633" height="8.633" fill="#e5eff9"><title>CP_severity 2.087, eGFR 17.5: 0.102</title></rect>
<rect x="191.7" y="206.7" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 2.217, eGFR 17.5: 0.0957</title></rect>
<rect x="200" y="206.7" width="8.633" height="8.633" fill="#e7f1fa"><title>CP_severity 2.348, eGFR 17.5: 0.0889</title></rect>
<rect x="208.3" y="206.7" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 2.478, eGFR 17.5: 0.0821</title></rect>
<rect x="216.7" y="206.7" width="8.633" height="8.633" fill="#eaf3fb"><title>CP_severity 2.609, eGFR 17.5: 0.0754</title></rect>
<rect x="225" y="206.7" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 2.739, eGFR 17.5: 0.0686</title></rect>
<rect x="233.3" y="206.7" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 2.87, eGFR 17.5: 0.0618</title></rect>
<rect x="241.7" y="206.7" width="8.633" height="8.633" fill="#eef5fc"><title>CP_severity 3, eGFR 17.5: 0.055</title></rect>
<rect x="50" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0, eGFR 12.5: 0.109</title></rect>
<rect x="58.33" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0.13, eGFR 12.5: 0.109</title></rect>
<rect x="66.67" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0.261, eGFR 12.5: 0.109</title></rect>
<rect x="75" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0.391, eGFR 12.5: 0.109</title></rect>
<rect x="83.33" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0.522, eGFR 12.5: 0.109</title></rect>
<rect x="91.67" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0.652, eGFR 12.5: 0.109</title></rect>
<rect x="100" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0.783, eGFR 12.5: 0.109</title></rect>
<rect x="108.3" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0.913, eGFR 12.5: 0.109</title></rect>
<rect x="116.7" y="215" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 1.043, eGFR 12.5: 0.108</title></rect>
<rect x="125" y="215" width="8.633" height="8.633" fill="#e5eff9"><title>CP_severity 1.174, eGFR 12.5: 0.103</title></rect>
<rect x="133.3" y="215" width="8.633" height="8.633" fill="#e5f0f9"><title>CP_severity 1.304, eGFR 12.5: 0.0992</title></rect>
<rect x="141.7" y="215" width="8.633" height="8.633" fill="#e6f0fa"><title>CP_severity 1.435, eGFR 12.5: 0.0949</title></rect>
<rect x="150" y="215" width="8.633" height="8.633" fill="#e7f1fa"><title>CP_severity 1.565, eGFR 12.5: 0.0907</title></rect>
<rect x="158.3" y="215" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 1.696, eGFR 12.5: 0.0864</title></rect>
<rect x="166.7" y="215" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 1.826, eGFR 12.5: 0.0821</title></rect>
<rect x="175" y="215" width="8.633" height="8.633" fill="#eaf2fb"><title>CP_severity 1.957, eGFR 12.5: 0.0778</title></rect>
<rect x="183.3" y="215" width="8.633" height="8.633" fill="#eaf3fb"><title>CP_severity 2.087, eGFR 12.5: 0.0732</title></rect>
<rect x="191.7" y="215" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 2.217, eGFR 12.5: 0.0684</title></rect>
<rect x="200" y="215" width="8.633" height="8.633" fill="#ecf4fc"><title>CP_severity 2.348, eGFR 12.5: 0.0635</title></rect>
<rect x="208.3" y="215" width="8.633" height="8.633" fill="#edf5fc"><title>CP_severity 2.478, eGFR 12.5: 0.0587</title></rect>
<rect x="216.7" y="215" width="8.633" height="8.633" fill="#eef5fc"><title>CP_severity 2.609, eGFR 12.5: 0.0538</title></rect>
<rect x="225" y="215" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 2.739, eGFR 12.5: 0.049</title></rect>
<rect x="233.3" y="215" width="8.633" height="8.633" fill="#f0f6fd"><title>CP_severity 2.87, eGFR 12.5: 0.0441</title></rect>
<rect x="241.7" y="215" width="8.633" height="8.633" fill="#f1f7fd"><title>CP_severity 3, eGFR 12.5: 0.0393</title></rect>
<rect x="50" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0, eGFR 7.5: 0.0655</title></rect>
<rect x="58.33" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.13, eGFR 7.5: 0.0655</title></rect>
<rect x="66.67" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.261, eGFR 7.5: 0.0655</title></rect>
<rect x="75" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.391, eGFR 7.5: 0.0655</title></rect>
<rect x="83.33" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.522, eGFR 7.5: 0.0655</title></rect>
<rect x="91.67" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.652, eGFR 7.5: 0.0655</title></rect>
<rect x="100" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.783, eGFR 7.5: 0.0655</title></rect>
<rect x="108.3" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.913, eGFR 7.5: 0.0655</title></rect>
<rect x="116.7" y="223.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 1.043, eGFR 7.5: 0.0646</title></rect>
<rect x="125" y="223.3" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 1.174, eGFR 7.5: 0.0621</title></rect>
<rect x="133.3" y="223.3" width="8.633" height="8.633" fill="#edf5fc"><title>CP_severity 1.304, eGFR 7.5: 0.0595</title></rect>
<rect x="141.7" y="223.3" width="8.633" height="8.633" fill="#eef5fc"><title>CP_severity 1.435, eGFR 7.5: 0.057</title></rect>
<rect x="150" y="223.3" width="8.633" height="8.633" fill="#eef5fc"><title>CP_severity 1.565, eGFR 7.5: 0.0544</title></rect>
<rect x="158.3" y="223.3" width="8.633" height="8.633" fill="#eff5fc"><title>CP_severity 1.696, eGFR 7.5: 0.0518</title></rect>
<rect x="166.7" y="223.3" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 1.826, eGFR 7.5: 0.0493</title></rect>
<rect x="175" y="223.3" width="8.633" height="8.633" fill="#f0f6fd"><title>CP_severity 1.957, eGFR 7.5: 0.0467</title></rect>
<rect x="183.3" y="223.3" width="8.633" height="8.633" fill="#f0f6fd"><title>CP_severity 2.087, eGFR 7.5: 0.0439</title></rect>
<rect x="191.7" y="223.3" width="8.633" height="8.633" fill="#f1f7fd"><title>CP_severity 2.217, eGFR 7.5: 0.041</title></rect>
<rect x="200" y="223.3" width="8.633" height="8.633" fill="#f1f7fd"><title>CP_severity 2.348, eGFR 7.5: 0.0381</title></rect>
<rect x="208.3" y="223.3" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 2.478, eGFR 7.5: 0.0352</title></rect>
<rect x="216.7" y="223.3" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 2.609, eGFR 7.5: 0.0323</title></rect>
<rect x="225" y="223.3" width="8.633" height="8.633" fill="#f3f8fe"><title>CP_severity 2.739, eGFR 7.5: 0.0294</title></rect>
<rect x="233.3" y="223.3" width="8.633" height="8.633" fill="#f3f9fe"><title>CP_severity 2.87, eGFR 7.5: 0.0265</title></rect>
<rect x="241.7" y="223.3" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 3, eGFR 7.5: 0.0236</title></rect>
<rect x="50" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0, eGFR 2.5: 0.0218</title></rect>
<rect x="58.33" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.13, eGFR 2.5: 0.0218</title></rect>
<rect x="66.67" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.261, eGFR 2.5: 0.0218</title></rect>
<rect x="75" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.391, eGFR 2.5: 0.0218</title></rect>
<rect x="83.33" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.522, eGFR 2.5: 0.0218</title></rect>
<rect x="91.67" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.652, eGFR 2.5: 0.0218</title></rect>
<rect x="100" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.783, eGFR 2.5: 0.0218</title></rect>
<rect x="108.3" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.913, eGFR 2.5: 0.0218</title></rect>
<rect x="116.7" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 1.043, eGFR 2.5: 0.0215</title></rect>
<rect x="125" y="231.7" width="8.633" height="8.633" fill="#f5f9fe"><title>CP_severity 1.174, eGFR 2.5: 0.0207</title></rect>
<rect x="133.3" y="231.7" width="8.633" height="8.633" fill="#f5f9fe"><title>CP_severity 1.304, eGFR 2.5: 0.0198</title></rect>
<rect x="141.7" y="231.7" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 1.435, eGFR 2.5: 0.019</title></rect>
<rect x="150" y="231.7" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 1.565, eGFR 2.5: 0.0181</title></rect>
<rect x="158.3" y="231.7" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 1.696, eGFR 2.5: 0.0173</title></rect>
<rect x="166.7" y="231.7" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 1.826, eGFR 2.5: 0.0164</title></rect>
<rect x="175" y="231.7" width="8.633" height="8.633" fill="#f6faff"><title>CP_severity 1.957, eGFR 2.5: 0.0156</title></rect>
<rect x="183.3" y="231.7" width="8.633" height="8.633" fill="#f6faff"><title>CP_severity 2.087, eGFR 2.5: 0.0146</title></rect>
<rect x="191.7" y="231.7" width="8.633" height="8.633" fill="#f6faff"><title>CP_severity 2.217, eGFR 2.5: 0.0137</title></rect>
<rect x="200" y="231.7" width="8.633" height="8.633" fill="#f6faff"><title>CP_severity 2.348, eGFR 2.5: 0.0127</title></rect>
<rect x="208.3" y="231.7" width="8.633" height="8.633" fill="#f6fbff"><title>CP_severity 2.478, eGFR 2.5: 0.0117</title></rect>
<rect x="216.7" y="231.7" width="8.633" height="8.633" fill="#f6fbff"><title>CP_severity 2.609, eGFR 2.5: 0.0108</title></rect>
<rect x="225" y="231.7" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 2.739, eGFR 2.5: 0.0098</title></rect>
<rect x="233.3" y="231.7" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 2.87, eGFR 2.5: 0.00883</title></rect>
<rect x="241.7" y="231.7" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 3, eGFR 2.5: 0.00786</title></rect>
<line x1="46" x2="50" y1="90" y2="90" stroke="currentColor"/>
<text x="44" y="94" text-anchor="end" font-size="11">90</text>
<line x1="46" x2="50" y1="140" y2="140" stroke="currentColor"/>
<text x="44" y="144" text-anchor="end" font-size="11">60</text>
<line x1="46" x2="50" y1="190" y2="190" stroke="currentColor"/>
<text x="44" y="194" text-anchor="end" font-size="11">30</text>
<line x1="46" x2="50" y1="215" y2="215" stroke="currentColor"/>
<text x="44" y="219" text-anchor="end" font-size="11">15</text>
<line x1="54.17" x2="54.17" y1="240" y2="244" stroke="currentColor"/>
<text x="54.17" y="257" text-anchor="middle" font-size="11">0</text>
<line x1="118.1" x2="118.1" y1="240" y2="244" stroke="currentColor"/>
<text x="118.1" y="257" text-anchor="middle" font-size="11">A</text>
<line x1="181.9" x2="181.9" y1="240" y2="244" stroke="currentColor"/>
<text x="181.9" y="257" text-anchor="middle" font-size="11">B</text>
<line x1="245.8" x2="245.8" y1="240" y2="244" stroke="currentColor"/>
<text x="245.8" y="257" text-anchor="middle" font-size="11">C</text>
<text x="150" y="273" text-anchor="middle" font-size="12">CP_severity</text>
<rect x="50" y="282" width="10.3" height="8" fill="#f7fbff"/>
<rect x="60" y="282" width="10.3" height="8" fill="#edf4fc"/>
<rect x="70" y="282" width="10.3" height="8" fill="#e2eef8"/>
<rect x="80" y="282" width="10.3" height="8" fill="#d8e7f5"/>
<rect x="90" y="282" width="10.3" height="8" fill="#cee0f2"/>
<rect x="100" y="282" width="10.3" height="8" fill="#c1d9ee"/>
<rect x="110" y="282" width="10.3" height="8" fill="#aecfe8"/>
<rect x="120" y="282" width="10.3" height="8" fill="#9bc6e3"/>
<rect x="130" y="282" width="10.3" height="8" fill="#88bcde"/>
<rect x="140" y="282" width="10.3" height="8" fill="#75b3d9"/>
<rect x="150" y="282" width="10.3" height="8" fill="#63a8d3"/>
<rect x="160" y="282" width="10.3" height="8" fill="#549bcc"/>
<rect x="170" y="282" width="10.3" height="8" fill="#448ec5"/>
<rect x="180" y="282" width="10.3" height="8" fill="#3481be"/>
<rect x="190" y="282" width="10.3" height="8" fill="#2574b7"/>
<rect x="200" y="282" width="10.3" height="8" fill="#1d67a9"/>
<rect x="210" y="282" width="10.3" height="8" fill="#18599a"/>
<rect x="220" y="282" width="10.3" height="8" fill="#134b8a"/>
<rect x="230" y="282" width="10.3" height="8" fill="#0d3e7b"/>
<rect x="240" y="282" width="10.3" height="8" fill="#08306b"/>
<text x="50" y="304" font-size="11">0.0079</text>
<text x="250" y="304" font-size="11" text-anchor="end">1</text>
<text x="300" y="18" font-weight="600">renal blood flow, cardiorenal</text>
<text transform="translate(272 140) rotate(-90)" text-anchor="middle" font-size="12">eGFR (ml/min/1.73 m2)</text>
<rect x="300" y="40" width="8.633" height="8.633" fill="#08306b"><title>HF_severity 0, eGFR 117.5: 1</title></rect>
<rect x="308.3" y="40" width="8.633" height="8.633" fill="#0b3773"><title>HF_severity 0.174, eGFR 117.5: 0.979</title></rect>
<rect x="316.7" y="40" width="8.633" height="8.633" fill="#0d3d7a"><title>HF_severity 0.348, eGFR 117.5: 0.958</title></rect>
<rect x="325" y="40" width="8.633" height="8.633" fill="#104482"><title>HF_severity 0.522, eGFR 117.5: 0.937</title></rect>
<rect x="333.3" y="40" width="8.633" height="8.633" fill="#124b8a"><title>HF_severity 0.696, eGFR 117.5: 0.916</title></rect>
<rect x="341.7" y="40" width="8.633" height="8.633" fill="#155291"><title>HF_severity 0.87, eGFR 117.5: 0.895</title></rect>
<rect x="350" y="40" width="8.633" height="8.633" fill="#175899"><title>HF_severity 1.043, eGFR 117.5: 0.874</title></rect>
<rect x="358.3" y="40" width="8.633" height="8.633" fill="#1a5fa1"><title>HF_severity 1.217, eGFR 117.5: 0.853</title></rect>
<rect x="366.7" y="40" width="8.633" height="8.633" fill="#1d66a8"><title>HF_severity 1.391, eGFR 117.5: 0.832</title></rect>
<rect x="375" y="40" width="8.633" height="8.633" fill="#1f6cb0"><title>HF_severity 1.565, eGFR 117.5: 0.811</title></rect>
<rect x="383.3" y="40" width="8.633" height="8.633" fill="#2373b6"><title>HF_severity 1.739, eGFR 117.5: 0.79</title></rect>
<rect x="391.7" y="40" width="8.633" height="8.633" fill="#2b79ba"><title>HF_severity 1.913, eGFR 117.5: 0.769</title></rect>
<rect x="400" y="40" width="8.633" height="8.633" fill="#3682be"><title>HF_severity 2.087, eGFR 117.5: 0.74</title></rect>
<rect x="408.3" y="40" width="8.633" height="8.633" fill="#438dc4"><title>HF_severity 2.261, eGFR 117.5: 0.704</title></rect>
<rect x="416.7" y="40" width="8.633" height="8.633" fill="#5098ca"><title>HF_severity 2.435, eGFR 117.5: 0.667</title></rect>
<rect x="425" y="40" width="8.633" height="8.633" fill="#5da3d0"><title>HF_severity 2.609, eGFR 117.5: 0.631</title></rect>
<rect x="433.3" y="40" width="8.633" height="8.633" fill="#6baed6"><title>HF_severity 2.783, eGFR 117.5: 0.595</title></rect>
<rect x="441.7" y="40" width="8.633" height="8.633" fill="#7bb6da"><title>HF_severity 2.957, eGFR 117.5: 0.558</title></rect>
<rect x="450" y="40" width="8.633" height="8.633" fill="#7ab5da"><title>HF_severity 3.13, eGFR 117.5: 0.56</title></rect>
<rect x="458.3" y="40" width="8.633" height="8.633" fill="#74b2d8"><title>HF_severity 3.304, eGFR 117.5: 0.574</title></rect>
<rect x="466.7" y="40" width="8.633" height="8.633" fill="#6eafd7"><title>HF_severity 3.478, eGFR 117.5: 0.588</title></rect>
<rect x="475" y="40" width="8.633" height="8.633" fill="#68abd5"><title>HF_severity 3.652, eGFR 117.5: 0.602</title></rect>
<rect x="483.3" y="40" width="8.633" height="8.633" fill="#63a7d2"><title>HF_severity 3.826, eGFR 117.5: 0.616</title></rect>
<rect x="491.7" y="40" width="8.633" height="8.633" fill="#5ea3d0"><title>HF_severity 4, eGFR 117.5: 0.63</title></rect>
<rect x="300" y="48.33" width="8.633" height="8.633" fill="#08306b"><title>HF_severity 0, eGFR 112.5: 1</title></rect>
<rect x="308.3" y="48.33" width="8.633" height="8.633" fill="#0b3773"><title>HF_severity 0.174, eGFR 112.5: 0.979</title></rect>
<rect x="316.7" y="48.33" width="8.633" height="8.633" fill="#0d3d7a"><title>HF_severity 0.348, eGFR 112.5: 0.958</title></rect>
<rect x="325" y="48.33" width="8.633" height="8.633" fill="#104482"><title>HF_severity 0.522, eGFR 112.5: 0.937</title></rect>
<rect x="333.3" y="48.33" width="8.633" height="8.633" fill="#124b8a"><title>HF_severity 0.696, eGFR 112.5: 0.916</title></rect>
<rect x="341.7" y="48.33" width="8.633" height="8.633" fill="#155291"><title>HF_severity 0.87, eGFR 112.5: 0.895</title></rect>
<rect x="350" y="48.33" width="8.633" height="8.633" fill="#175899"><title>HF_severity 1.043, eGFR 112.5: 0.874</title></rect>
<rect x="358.3" y="48.33" width="8.633" height="8.633" fill="#1a5fa1"><title>HF_severity 1.217, eGFR 112.5: 0.853</title></rect>
<rect x="366.7" y="48.33" width="8.633" height="8.633" fill="#1d66a8"><title>HF_severity 1.391, eGFR 112.5: 0.832</title></rect>
<rect x="375" y="48.33" width="8.633" height="8.633" fill="#1f6cb0"><title>HF_severity 1.565, eGFR 112.5: 0.811</title></rect>
<rect x="383.3" y="48.33" width="8.633" height="8.633" fill="#2373b6"><title>HF_severity 1.739, eGFR 112.5: 0.79</title></rect>
<rect x="391.7" y="48.33" width="8.633" height="8.633" fill="#2b79ba"><title>HF_severity 1.913, eGFR 112.5: 0.769</title></rect>
<rect x="400" y="48.33" width="8.633" height="8.633" fill="#3682be"><title>HF_severity 2.087, eGFR 112.5: 0.74</title></rect>
<rect x="408.3" y="48.33" width="8.633" height="8.633" fill="#438dc4"><title>HF_severity 2.261, eGFR 112.5: 0.704</title></rect>
<rect x="416.7" y="48.33" width="8.633" height="8.633" fill="#5098ca"><title>HF_severity 2.435, eGFR 112.5: 0.667</title></rect>
<rect x="425" y="48.33" width="8.633" height="8.633" fill="#5da3d0"><title>HF_severity 2.609, eGFR 112.5: 0.631</title></rect>
<rect x="433.3" y="48.33" width="8.633" height="8.633" fill="#6baed6"><title>HF_severity 2.783, eGFR 112.5: 0.595</title></rect>
<rect x="441.7" y="48.33" width="8.633" height="8.633" fill="#7bb6da"><title>HF_severity 2.957, eGFR 112.5: 0.558</title></rect>
<rect x="450" y="48.33" width="8.633" height="8.633" fill="#7ab5da"><title>HF_severity 3.13, eGFR 112.5: 0.56</title></rect>
<rect x="458.3" y="48.33" width="8.633" height="8.633" fill="#74b2d8"><title>HF_severity 3.304, eGFR 112.5: 0.574</title></rect>
<rect x="466.7" y="48.33" width="8.633" height="8.633" fill="#6eafd7"><title>HF_severity 3.478, eGFR 112.5: 0.588</title></rect>
<rect x="475" y="48.33" width="8.633" height="8.633" fill="#68abd5"><title>HF_severity 3.652, eGFR 112.5: 0.602</title></rect>
<rect x="483.3" y="48.33" width="8.633" height="8.633" fill="#63a7d2"><title>HF_severity 3.826, eGFR 112.5: 0.616</title></rect>
<rect x="491.7" y="48.33" width="8.633" height="8.633" fill="#5ea3d0"><title>HF_severity 4, eGFR 112.5: 0.63</title></rect>
<rect x="300" y="56.67" width="8.633" height="8.633" fill="#08306b"><title>HF_severity 0, eGFR 107.5: 1</title></rect>
<rect x="308.3" y="56.67" width="8.633" height="8.633" fill="#0b3773"><title>HF_severity 0.174, eGFR 107.5: 0.979</title></rect>
<rect x="316.7" y="56.67" width="8.633" height="8.633" fill="#0d3d7a"><title>HF_severity 0.348, eGFR 107.5: 0.958</title></rect>
<rect x="325" y="56.67" width="8.633" height="8.633" fill="#104482"><title>HF_severity 0.522, eGFR 107.5: 0.937</title></rect>
<rect x="333.3" y="56.67" width="8.633" height="8.633" fill="#124b8a"><title>HF_severity 0.696, eGFR 107.5: 0.916</title></rect>
<rect x="341.7" y="56.67" width="8.633" height="8.633" fill="#155291"><title>HF_severity 0.87, eGFR 107.5: 0.895</title></rect>
<rect x="350" y="56.67" width="8.633" height="8.633" fill="#175899"><title>HF_severity 1.043, eGFR 107.5: 0.874</title></rect>
<rect x="358.3" y="56.67" width="8.633" height="8.633" fill="#1a5fa1"><title>HF_severity 1.217, eGFR 107.5: 0.853</title></rect>
<rect x="366.7" y="56.67" width="8.633" height="8.633" fill="#1d66a8"><title>HF_severity 1.391, eGFR 107.5: 0.832</title></rect>
<rect x="375" y="56.67" width="8.633" height="8.633" fill="#1f6cb0"><title>HF_severity 1.565, eGFR 107.5: 0.811</title></rect>
<rect x="383.3" y="56.67" width="8.633" height="8.633" fill="#2373b6"><title>HF_severity 1.739, eGFR 107.5: 0.79</title></rect>
<rect x="391.7" y="56.67" width="8.633" height="8.633" fill="#2b79ba"><title>HF_severity 1.913, eGFR 107.5: 0.769</title></rect>
<rect x="400" y="56.67" width="8.633" height="8.633" fill="#3682be"><title>HF_severity 2.087, eGFR 107.5: 0.74</title></rect>
<rect x="408.3" y="56.67" width="8.633" height="8.633" fill="#438dc4"><title>HF_severity 2.261, eGFR 107.5: 0.704</title></rect>
<rect x="416.7" y="56.67" width="8.633" height="8.633" fill="#5098ca"><title>HF_severity 2.435, eGFR 107.5: 0.667</title></rect>
<rect x="425" y="56.67" width="8.633" height="8.633" fill="#5da3d0"><title>HF_severity 2.609, eGFR 107.5: 0.631</title></rect>
<rect x="433.3" y="56.67" width="8.633" height="8.633" fill="#6baed6"><title>HF_severity 2.783, eGFR 107.5: 0.595</title></rect>
<rect x="441.7" y="56.67" width="8.633" height="8.633" fill="#7bb6da"><title>HF_severity 2.957, eGFR 107.5: 0.558</title></rect>
<rect x="450" y="56.67" width="8.633" height="8.633" fill="#7ab5da"><title>HF_severity 3.13, eGFR 107.5: 0.56</title></rect>
<rect x="458.3" y="56.67" width="8.633" height="8.633" fill="#74b2d8"><title>HF_severity 3.304, eGFR 107.5: 0.574</title></rect>
<rect x="466.7" y="56.67" width="8.633" height="8.633" fill="#6eafd7"><title>HF_severity 3.478, eGFR 107.5: 0.588</title></rect>
<rect x="475" y="56.67" width="8.633" height="8.633" fill="#68abd5"><title>HF_severity 3.652, eGFR 107.5: 0.602</title></rect>
<rect x="483.3" y="56.67" width="8.633" height="8.633" fill="#63a7d2"><title>HF_severity 3.826, eGFR 107.5: 0.616</title></rect>
<rect x="491.7" y="56.67" width="8.633" height="8.633" fill="#5ea3d0"><title>HF_severity 4, eGFR 107.5: 0.63</title></rect>
<rect x="300" y="65" width="8.633" height="8.633" fill="#08306b"><title>HF_severity 0, eGFR 102.5: 1</title></rect>
<rect x="308.3" y="65" width="8.633" height="8.633" fill="#0b3773"><title>HF_severity 0.174, eGFR 102.5: 0.979</title></rect>
<rect x="316.7" y="65" width="8.633" height="8.633" fill="#0d3d7a"><title>HF_severity 0.348, eGFR 102.5: 0.958</title></rect>
<rect x="325" y="65" width="8.633" height="8.633" fill="#104482"><title>HF_severity 0.522, eGFR 102.5: 0.937</title></rect>
<rect x="333.3" y="65" width="8.633" height="8.633" fill="#124b8a"><title>HF_severity 0.696, eGFR 102.5: 0.916</title></rect>
<rect x="341.7" y="65" width="8.633" height="8.633" fill="#155291"><title>HF_severity 0.87, eGFR 102.5: 0.895</title></rect>
<rect x="350" y="65" width="8.633" height="8.633" fill="#175899"><title>HF_severity 1.043, eGFR 102.5: 0.874</title></rect>
<rect x="358.3" y="65" width="8.633" height="8.633" fill="#1a5fa1"><title>HF_severity 1.217, eGFR 102.5: 0.853</title></rect>
<rect x="366.7" y="65" width="8.633" height="8.633" fill="#1d66a8"><title>HF_severity 1.391, eGFR 102.5: 0.832</title></rect>
<rect x="375" y="65" width="8.633" height="8.633" fill="#1f6cb0"><title>HF_severity 1.565, eGFR 102.5: 0.811</title></rect>
<rect x="383.3" y="65" width="8.633" height="8.633" fill="#2373b6"><title>HF_severity 1.739, eGFR 102.5: 0.79</title></rect>
<rect x="391.7" y="65" width="8.633" height="8.633" fill="#2b79ba"><title>HF_severity 1.913, eGFR 102.5: 0.769</title></rect>
<rect x="400" y="65" width="8.633" height="8.633" fill="#3682be"><title>HF_severity 2.087, eGFR 102.5: 0.74</title></rect>
<rect x="408.3" y="65" width="8.633" height="8.633" fill="#438dc4"><title>HF_severity 2.261, eGFR 102.5: 0.704</title></rect>
<rect x="416.7" y="65" width="8.633" height="8.633" fill="#5098ca"><title>HF_severity 2.435, eGFR 102.5: 0.667</title></rect>
<rect x="425" y="65" width="8.633" height="8.633" fill="#5da3d0"><title>HF_severity 2.609, eGFR 102.5: 0.631</title></rect>
<rect x="433.3" y="65" width="8.633" height="8.633" fill="#6baed6"><title>HF_severity 2.783, eGFR 102.5: 0.595</title></rect>
<rect x="441.7" y="65" width="8.633" height="8.633" fill="#7bb6da"><title>HF_severity 2.957, eGFR 102.5: 0.558</title></rect>
<rect x="450" y="65" width="8.633" height="8.633" fill="#7ab5da"><title>HF_severity 3.13, eGFR 102.5: 0.56</title></rect>
<rect x="458.3" y="65" width="8.633" height="8.633" fill="#74b2d8"><title>HF_severity 3.304, eGFR 102.5: 0.574</title></rect>
<rect x="466.7" y="65" width="8.633" height="8.633" fill="#6eafd7"><title>HF_severity 3.478, eGFR 102.5: 0.588</title></rect>
<rect x="475" y="65" width="8.633" height="8.633" fill="#68abd5"><title>HF_severity 3.652, eGFR 102.5: 0.602</title></rect>
<rect x="483.3" y="65" width="8.633" height="8.633" fill="#63a7d2"><title>HF_severity 3.826, eGFR 102.5: 0.616</title></rect>
<rect x="491.7" y="65" width="8.633" height="8.633" fill="#5ea3d0"><title>HF_severity 4, eGFR 102.5: 0.63</title></rect>
<rect x="300" y="73.33" width="8.633" height="8.633" fill="#08306b"><title>HF_severity 0, eGFR 97.5: 1</title></rect>
<rect x="308.3" y="73.33" width="8.633" height="8.633" fill="#0b3773"><title>HF_severity 0.174, eGFR 97.5: 0.979</title></rect>
<rect x="316.7" y="73.33" width="8.633" height="8.633" fill="#0d3d7a"><title>HF_severity 0.348, eGFR 97.5: 0.958</title></rect>
<rect x="325" y="73.33" width="8.633" height="8.633" fill="#104482"><title>HF_severity 0.522, eGFR 97.5: 0.937</title></rect>
<rect x="333.3" y="73.33" width="8.633" height="8.633" fill="#124b8a"><title>HF_severity 0.696, eGFR 97.5: 0.916</title></rect>
<rect x="341.7" y="73.33" width="8.633" height="8.633" fill="#155291"><title>HF_severity 0.87, eGFR 97.5: 0.895</title></rect>
<rect x="350" y="73.33" width="8.633" height="8.633" fill="#175899"><title>HF_severity 1.043, eGFR 97.5: 0.874</title></rect>
<rect x="358.3" y="73.33" width="8.633" height="8.633" fill="#1a5fa1"><title>HF_severity 1.217, eGFR 97.5: 0.853</title></rect>
<rect x="366.7" y="73.33" width="8.633" height="8.633" fill="#1d66a8"><title>HF_severity 1.391, eGFR 97.5: 0.832</title></rect>
<rect x="375" y="73.33" width="8.633" height="8.633" fill="#1f6cb0"><title>HF_severity 1.565, eGFR 97.5: 0.811</title></rect>
<rect x="383.3" y="73.33" width="8.633" height="8.633" fill="#2373b6"><title>HF_severity 1.739, eGFR 97.5: 0.79</title></rect>
<rect x="391.7" y="73.33" width="8.633" height="8.633" fill="#2b79ba"><title>HF_severity 1.913, eGFR 97.5: 0.769</title></rect>
<rect x="400" y="73.33" width="8.633" height="8.633" fill="#3682be"><title>HF_severity 2.087, eGFR 97.5: 0.74</title></rect>
<rect x="408.3" y="73.33" width="8.633" height="8.633" fill="#438dc4"><title>HF_severity 2.261, eGFR 97.5: 0.704</title></rect>
<rect x="416.7" y="73.33" width="8.633" height="8.633" fill="#5098ca"><title>HF_severity 2.435, eGFR 97.5: 0.667</title></rect>
<rect x="425" y="73.33" width="8.633" height="8.633" fill="#5da3d0"><title>HF_severity 2.609, eGFR 97.5: 0.631</title></rect>
<rect x="433.3" y="73.33" width="8.633" height="8.633" fill="#6baed6"><title>HF_severity 2.783, eGFR 97.5: 0.595</title></rect>
<rect x="441.7" y="73.33" width="8.633" height="8.633" fill="#7bb6da"><title>HF_severity 2.957, eGFR 97.5: 0.558</title></rect>
<rect x="450" y="73.33" width="8.633" height="8.633" fill="#7ab5da"><title>HF_severity 3.13, eGFR 97.5: 0.56</title></rect>
<rect x="458.3" y="73.33" width="8.633" height="8.633" fill="#74b2d8"><title>HF_severity 3.304, eGFR 97.5: 0.574</title></rect>
<rect x="466.7" y="73.33" width="8.633" height="8.633" fill="#6eafd7"><title>HF_severity 3.478, eGFR 97.5: 0.588</title></rect>
<rect x="475" y="73.33" width="8.633" height="8.633" fill="#68abd5"><title>HF_severity 3.652, eGFR 97.5: 0.602</title></rect>
<rect x="483.3" y="73.33" width="8.633" height="8.633" fill="#63a7d2"><title>HF_severity 3.826, eGFR 97.5: 0.616</title></rect>
<rect x="491.7" y="73.33" width="8.633" height="8.633" fill="#5ea3d0"><title>HF_severity 4, eGFR 97.5: 0.63</title></rect>
<rect x="300" y="81.67" width="8.633" height="8.633" fill="#08306b"><title>HF_severity 0, eGFR 92.5: 1</title></rect>
<rect x="308.3" y="81.67" width="8.633" height="8.633" fill="#0b3773"><title>HF_severity 0.174, eGFR 92.5: 0.979</title></rect>
<rect x="316.7" y="81.67" width="8.633" height="8.633" fill="#0d3d7a"><title>HF_severity 0.348, eGFR 92.5: 0.958</title></rect>
<rect x="325" y="81.67" width="8.633" height="8.633" fill="#104482"><title>HF_severity 0.522, eGFR 92.5: 0.937</title></rect>
<rect x="333.3" y="81.67" width="8.633" height="8.633" fill="#124b8a"><title>HF_severity 0.696, eGFR 92.5: 0.916</title></rect>
<rect x="341.7" y="81.67" width="8.633" height="8.633" fill="#155291"><title>HF_severity 0.87, eGFR 92.5: 0.895</title></rect>
<rect x="350" y="81.67" width="8.633" height="8.633" fill="#175899"><title>HF_severity 1.043, eGFR 92.5: 0.874</title></rect>
<rect x="358.3" y="81.67" width="8.633" height="8.633" fill="#1a5fa1"><title>HF_severity 1.217, eGFR 92.5: 0.853</title></rect>
<rect x="366.7" y="81.67" width="8.633" height="8.633" fill="#1d66a8"><title>HF_severity 1.391, eGFR 92.5: 0.832</title></rect>
<rect x="375" y="81.67" width="8.633" height="8.633" fill="#1f6cb0"><title>HF_severity 1.565, eGFR 92.5: 0.811</title></rect>
<rect x="383.3" y="81.67" width="8.633" height="8.633" fill="#2373b6"><title>HF_severity 1.739, eGFR 92.5: 0.79</title></rect>
<rect x="391.7" y="81.67" width="8.633" height="8.633" fill="#2b79ba"><title>HF_severity 1.913, eGFR 92.5: 0.769</title></rect>
<rect x="400" y="81.67" width="8.633" height="8.633" fill="#3682be"><title>HF_severity 2.087, eGFR 92.5: 0.74</title></rect>
<rect x="408.3" y="81.67" width="8.633" height="8.633" fill="#438dc4"><title>HF_severity 2.261, eGFR 92.5: 0.704</title></rect>
<rect x="416.7" y="81.67" width="8.633" height="8.633" fill="#5098ca"><title>HF_severity 2.435, eGFR 92.5: 0.667</title></rect>
<rect x="425" y="81.67" width="8.633" height="8.633" fill="#5da3d0"><title>HF_severity 2.609, eGFR 92.5: 0.631</title></rect>
<rect x="433.3" y="81.67" width="8.633" height="8.633" fill="#6baed6"><title>HF_severity 2.783, eGFR 92.5: 0.595</title></rect>
<rect x="441.7" y="81.67" width="8.633" height="8.633" fill="#7bb6da"><title>HF_severity 2.957, eGFR 92.5: 0.558</title></rect>
<rect x="450" y="81.67" width="8.633" height="8.633" fill="#7ab5da"><title>HF_severity 3.13, eGFR 92.5: 0.56</title></rect>
<rect x="458.3" y="81.67" width="8.633" height="8.633" fill="#74b2d8"><title>HF_severity 3.304, eGFR 92.5: 0.574</title></rect>
<rect x="466.7" y="81.67" width="8.633" height="8.633" fill="#6eafd7"><title>HF_severity 3.478, eGFR 92.5: 0.588</title></rect>
<rect x="475" y="81.67" width="8.633" height="8.633" fill="#68abd5"><title>HF_severity 3.652, eGFR 92.5: 0.602</title></rect>
<rect x="483.3" y="81.67" width="8.633" height="8.633" fill="#63a7d2"><title>HF_severity 3.826, eGFR 92.5: 0.616</title></rect>
<rect x="491.7" y="81.67" width="8.633" height="8.633" fill="#5ea3d0"><title>HF_severity 4, eGFR 92.5: 0.63</title></rect>
<rect x="300" y="90" width="8.633" height="8.633" fill="#0d3c79"><title>HF_severity 0, eGFR 87.5: 0.962</title></rect>
<rect x="308.3" y="90" width="8.633" height="8.633" fill="#0f4380"><title>HF_severity 0.174, eGFR 87.5: 0.942</title></rect>
<rect x="316.7" y="90" width="8.633" height="8.633" fill="#124988"><title>HF_severity 0.348, eGFR 87.5: 0.921</title></rect>
<rect x="325" y="90" width="8.633" height="8.633" fill="#14508f"><title>HF_severity 0.522, eGFR 87.5: 0.901</title></rect>
<rect x="333.3" y="90" width="8.633" height="8.633" fill="#175696"><title>HF_severity 0.696, eGFR 87.5: 0.881</title></rect>
<rect x="341.7" y="90" width="8.633" height="8.633" fill="#195c9e"><title>HF_severity 0.87, eGFR 87.5: 0.861</title></rect>
<rect x="350" y="90" width="8.633" height="8.633" fill="#1c63a5"><title>HF_severity 1.043, eGFR 87.5: 0.841</title></rect>
<rect x="358.3" y="90" width="8.633" height="8.633" fill="#1e69ac"><title>HF_severity 1.217, eGFR 87.5: 0.82</title></rect>
<rect x="366.7" y="90" width="8.633" height="8.633" fill="#2170b4"><title>HF_severity 1.391, eGFR 87.5: 0.8</title></rect>
<rect x="375" y="90" width="8.633" height="8.633" fill="#2776b8"><title>HF_severity 1.565, eGFR 87.5: 0.78</title></rect>
<rect x="383.3" y="90" width="8.633" height="8.633" fill="#2e7cbb"><title>HF_severity 1.739, eGFR 87.5: 0.76</title></rect>
<rect x="391.7" y="90" width="8.633" height="8.633" fill="#3682be"><title>HF_severity 1.913, eGFR 87.5: 0.74</title></rect>
<rect x="400" y="90" width="8.633" height="8.633" fill="#408ac3"><title>HF_severity 2.087, eGFR 87.5: 0.712</title></rect>
<rect x="408.3" y="90" width="8.633" height="8.633" fill="#4d95c8"><title>HF_severity 2.261, eGFR 87.5: 0.677</title></rect>
<rect x="416.7" y="90" width="8.633" height="8.633" fill="#599fce"><title>HF_severity 2.435, eGFR 87.5: 0.642</title></rect>
<rect x="425" y="90" width="8.633" height="8.633" fill="#66aad4"><title>HF_severity 2.609, eGFR 87.5: 0.607</title></rect>
<rect x="433.3" y="90" width="8.633" height="8.633" fill="#75b3d9"><title>HF_severity 2.783, eGFR 87.5: 0.572</title></rect>
<rect x="441.7" y="90" width="8.633" height="8.633" fill="#84badd"><title>HF_severity 2.957, eGFR 87.5: 0.537</title></rect>
<rect x="450" y="90" width="8.633" height="8.633" fill="#84badd"><title>HF_severity 3.13, eGFR 87.5: 0.538</title></rect>
<rect x="458.3" y="90" width="8.633" height="8.633" fill="#7eb7db"><title>HF_severity 3.304, eGFR 87.5: 0.552</title></rect>
<rect x="466.7" y="90" width="8.633" height="8.633" fill="#78b4d9"><title>HF_severity 3.478, eGFR 87.5: 0.565</title></rect>
<rect x="475" y="90" width="8.633" height="8.633" fill="#72b1d8"><title>HF_severity 3.652, eGFR 87.5: 0.579</title></rect>
<rect x="483.3" y="90" width="8.633" height="8.633" fill="#6baed6"><title>HF_severity 3.826, eGFR 87.5: 0.592</title></rect>
<rect x="491.7" y="90" width="8.633" height="8.633" fill="#66aad4"><title>HF_severity 4, eGFR 87.5: 0.606</title></rect>
<rect x="300" y="98.33" width="8.633" height="8.633" fill="#165595"><title>HF_severity 0, eGFR 82.5: 0.886</title></rect>
<rect x="308.3" y="98.33" width="8.633" height="8.633" fill="#185b9b"><title>HF_severity 0.174, eGFR 82.5: 0.867</title></rect>
<rect x="316.7" y="98.33" width="8.633" height="8.633" fill="#1b60a2"><title>HF_severity 0.348, eGFR 82.5: 0.848</title></rect>
<rect x="325" y="98.33" width="8.633" height="8.633" fill="#1d66a9"><title>HF_severity 0.522, eGFR 82.5: 0.83</title></rect>
<rect x="333.3" y="98.33" width="8.633" height="8.633" fill="#1f6cb0"><title>HF_severity 0.696, eGFR 82.5: 0.811</title></rect>
<rect x="341.7" y="98.33" width="8.633" height="8.633" fill="#2272b6"><title>HF_severity 0.87, eGFR 82.5: 0.793</title></rect>
<rect x="350" y="98.33" width="8.633" height="8.633" fill="#2978b9"><title>HF_severity 1.043, eGFR 82.5: 0.774</title></rect>
<rect x="358.3" y="98.33" width="8.633" height="8.633" fill="#307dbc"><title>HF_severity 1.217, eGFR 82.5: 0.755</title></rect>
<rect x="366.7" y="98.33" width="8.633" height="8.633" fill="#3783bf"><title>HF_severity 1.391, eGFR 82.5: 0.737</title></rect>
<rect x="375" y="98.33" width="8.633" height="8.633" fill="#3e89c2"><title>HF_severity 1.565, eGFR 82.5: 0.718</title></rect>
<rect x="383.3" y="98.33" width="8.633" height="8.633" fill="#448ec5"><title>HF_severity 1.739, eGFR 82.5: 0.7</title></rect>
<rect x="391.7" y="98.33" width="8.633" height="8.633" fill="#4b94c8"><title>HF_severity 1.913, eGFR 82.5: 0.681</title></rect>
<rect x="400" y="98.33" width="8.633" height="8.633" fill="#549bcc"><title>HF_severity 2.087, eGFR 82.5: 0.656</title></rect>
<rect x="408.3" y="98.33" width="8.633" height="8.633" fill="#60a5d1"><title>HF_severity 2.261, eGFR 82.5: 0.623</title></rect>
<rect x="416.7" y="98.33" width="8.633" height="8.633" fill="#6caed6"><title>HF_severity 2.435, eGFR 82.5: 0.591</title></rect>
<rect x="425" y="98.33" width="8.633" height="8.633" fill="#7ab6da"><title>HF_severity 2.609, eGFR 82.5: 0.559</title></rect>
<rect x="433.3" y="98.33" width="8.633" height="8.633" fill="#89bdde"><title>HF_severity 2.783, eGFR 82.5: 0.527</title></rect>
<rect x="441.7" y="98.33" width="8.633" height="8.633" fill="#97c4e2"><title>HF_severity 2.957, eGFR 82.5: 0.494</title></rect>
<rect x="450" y="98.33" width="8.633" height="8.633" fill="#97c4e2"><title>HF_severity 3.13, eGFR 82.5: 0.496</title></rect>
<rect x="458.3" y="98.33" width="8.633" height="8.633" fill="#91c1e0"><title>HF_severity 3.304, eGFR 82.5: 0.508</title></rect>
<rect x="466.7" y="98.33" width="8.633" height="8.633" fill="#8cbedf"><title>HF_severity 3.478, eGFR 82.5: 0.521</title></rect>
<rect x="475" y="98.33" width="8.633" height="8.633" fill="#86bbdd"><title>HF_severity 3.652, eGFR 82.5: 0.533</title></rect>
<rect x="483.3" y="98.33" width="8.633" height="8.633" fill="#80b9dc"><title>HF_severity 3.826, eGFR 82.5: 0.545</title></rect>
<rect x="491.7" y="98.33" width="8.633" height="8.633" fill="#7bb6da"><title>HF_severity 4, eGFR 82.5: 0.558</title></rect>
<rect x="300" y="106.7" width="8.633" height="8.633" fill="#1f6db0"><title>HF_severity 0, eGFR 77.5: 0.809</title></rect>
<rect x="308.3" y="106.7" width="8.633" height="8.633" fill="#2372b6"><title>HF_severity 0.174, eGFR 77.5: 0.792</title></rect>
<rect x="316.7" y="106.7" width="8.633" height="8.633" fill="#2977b8"><title>HF_severity 0.348, eGFR 77.5: 0.775</title></rect>
<rect x="325" y="106.7" width="8.633" height="8.633" fill="#2f7dbb"><title>HF_severity 0.522, eGFR 77.5: 0.758</title></rect>
<rect x="333.3" y="106.7" width="8.633" height="8.633" fill="#3582be"><title>HF_severity 0.696, eGFR 77.5: 0.741</title></rect>
<rect x="341.7" y="106.7" width="8.633" height="8.633" fill="#3b87c1"><title>HF_severity 0.87, eGFR 77.5: 0.724</title></rect>
<rect x="350" y="106.7" width="8.633" height="8.633" fill="#428cc4"><title>HF_severity 1.043, eGFR 77.5: 0.707</title></rect>
<rect x="358.3" y="106.7" width="8.633" height="8.633" fill="#4891c6"><title>HF_severity 1.217, eGFR 77.5: 0.69</title></rect>
<rect x="366.7" y="106.7" width="8.633" height="8.633" fill="#4e96c9"><title>HF_severity 1.391, eGFR 77.5: 0.673</title></rect>
<rect x="375" y="106.7" width="8.633" height="8.633" fill="#549bcc"><title>HF_severity 1.565, eGFR 77.5: 0.656</title></rect>
<rect x="383.3" y="106.7" width="8.633" height="8.633" fill="#5aa0cf"><title>HF_severity 1.739, eGFR 77.5: 0.639</title></rect>
<rect x="391.7" y="106.7" width="8.633" height="8.633" fill="#60a5d1"><title>HF_severity 1.913, eGFR 77.5: 0.622</title></rect>
<rect x="400" y="106.7" width="8.633" height="8.633" fill="#69acd5"><title>HF_severity 2.087, eGFR 77.5: 0.599</title></rect>
<rect x="408.3" y="106.7" width="8.633" height="8.633" fill="#76b3d9"><title>HF_severity 2.261, eGFR 77.5: 0.57</title></rect>
<rect x="416.7" y="106.7" width="8.633" height="8.633" fill="#83badd"><title>HF_severity 2.435, eGFR 77.5: 0.54</title></rect>
<rect x="425" y="106.7" width="8.633" height="8.633" fill="#90c0e0"><title>HF_severity 2.609, eGFR 77.5: 0.511</title></rect>
<rect x="433.3" y="106.7" width="8.633" height="8.633" fill="#9dc7e4"><title>HF_severity 2.783, eGFR 77.5: 0.481</title></rect>
<rect x="441.7" y="106.7" width="8.633" height="8.633" fill="#aacde7"><title>HF_severity 2.957, eGFR 77.5: 0.452</title></rect>
<rect x="450" y="106.7" width="8.633" height="8.633" fill="#aacde7"><title>HF_severity 3.13, eGFR 77.5: 0.453</title></rect>
<rect x="458.3" y="106.7" width="8.633" height="8.633" fill="#a5cbe6"><title>HF_severity 3.304, eGFR 77.5: 0.464</title></rect>
<rect x="466.7" y="106.7" width="8.633" height="8.633" fill="#a0c8e4"><title>HF_severity 3.478, eGFR 77.5: 0.476</title></rect>
<rect x="475" y="106.7" width="8.633" height="8.633" fill="#9bc6e3"><title>HF_severity 3.652, eGFR 77.5: 0.487</title></rect>
<rect x="483.3" y="106.7" width="8.633" height="8.633" fill="#96c3e2"><title>HF_severity 3.826, eGFR 77.5: 0.498</title></rect>
<rect x="491.7" y="106.7" width="8.633" height="8.633" fill="#90c1e0"><title>HF_severity 4, eGFR 77.5: 0.51</title></rect>
<rect x="300" y="115" width="8.633" height="8.633" fill="#317ebc"><title>HF_severity 0, eGFR 72.5: 0.753</title></rect>
<rect x="308.3" y="115" width="8.633" height="8.633" fill="#3783bf"><title>HF_severity 0.174, eGFR 72.5: 0.737</title></rect>
<rect x="316.7" y="115" width="8.633" height="8.633" fill="#3c88c1"><title>HF_severity 0.348, eGFR 72.5: 0.722</title></rect>
<rect x="325" y="115" width="8.633" height="8.633" fill="#428cc4"><title>HF_severity 0.522, eGFR 72.5: 0.706</title></rect>
<rect x="333.3" y="115" width="8.633" height="8.633" fill="#4891c6"><title>HF_severity 0.696, eGFR 72.5: 0.69</title></rect>
<rect x="341.7" y="115" width="8.633" height="8.633" fill="#4e96c9"><title>HF_severity 0.87, eGFR 72.5: 0.674</title></rect>
<rect x="350" y="115" width="8.633" height="8.633" fill="#539acb"><title>HF_severity 1.043, eGFR 72.5: 0.658</title></rect>
<rect x="358.3" y="115" width="8.633" height="8.633" fill="#599fce"><title>HF_severity 1.217, eGFR 72.5: 0.643</title></rect>
<rect x="366.7" y="115" width="8.633" height="8.633" fill="#5fa4d1"><title>HF_severity 1.391, eGFR 72.5: 0.627</title></rect>
<rect x="375" y="115" width="8.633" height="8.633" fill="#65a9d3"><title>HF_severity 1.565, eGFR 72.5: 0.611</title></rect>
<rect x="383.3" y="115" width="8.633" height="8.633" fill="#6aadd6"><title>HF_severity 1.739, eGFR 72.5: 0.595</title></rect>
<rect x="391.7" y="115" width="8.633" height="8.633" fill="#71b1d8"><title>HF_severity 1.913, eGFR 72.5: 0.579</title></rect>
<rect x="400" y="115" width="8.633" height="8.633" fill="#7bb6da"><title>HF_severity 2.087, eGFR 72.5: 0.558</title></rect>
<rect x="408.3" y="115" width="8.633" height="8.633" fill="#87bcde"><title>HF_severity 2.261, eGFR 72.5: 0.53</title></rect>
<rect x="416.7" y="115" width="8.633" height="8.633" fill="#94c2e1"><title>HF_severity 2.435, eGFR 72.5: 0.503</title></rect>
<rect x="425" y="115" width="8.633" height="8.633" fill="#a0c8e5"><title>HF_severity 2.609, eGFR 72.5: 0.475</title></rect>
<rect x="433.3" y="115" width="8.633" height="8.633" fill="#accee8"><title>HF_severity 2.783, eGFR 72.5: 0.448</title></rect>
<rect x="441.7" y="115" width="8.633" height="8.633" fill="#b8d4eb"><title>HF_severity 2.957, eGFR 72.5: 0.421</title></rect>
<rect x="450" y="115" width="8.633" height="8.633" fill="#b8d4eb"><title>HF_severity 3.13, eGFR 72.5: 0.422</title></rect>
<rect x="458.3" y="115" width="8.633" height="8.633" fill="#b3d2ea"><title>HF_severity 3.304, eGFR 72.5: 0.432</title></rect>
<rect x="466.7" y="115" width="8.633" height="8.633" fill="#aecfe9"><title>HF_severity 3.478, eGFR 72.5: 0.443</title></rect>
<rect x="475" y="115" width="8.633" height="8.633" fill="#aacde7"><title>HF_severity 3.652, eGFR 72.5: 0.453</title></rect>
<rect x="483.3" y="115" width="8.633" height="8.633" fill="#a5cbe6"><title>HF_severity 3.826, eGFR 72.5: 0.464</title></rect>
<rect x="491.7" y="115" width="8.633" height="8.633" fill="#a0c8e5"><title>HF_severity 4, eGFR 72.5: 0.474</title></rect>
<rect x="300" y="123.3" width="8.633" height="8.633" fill="#3e89c2"><title>HF_severity 0, eGFR 67.5: 0.718</title></rect>
<rect x="308.3" y="123.3" width="8.633" height="8.633" fill="#438dc4"><title>HF_severity 0.174, eGFR 67.5: 0.702</title></rect>
<rect x="316.7" y="123.3" width="8.633" height="8.633" fill="#4992c7"><title>HF_severity 0.348, eGFR 67.5: 0.687</title></rect>
<rect x="325" y="123.3" width="8.633" height="8.633" fill="#4e96c9"><title>HF_severity 0.522, eGFR 67.5: 0.672</title></rect>
<rect x="333.3" y="123.3" width="8.633" height="8.633" fill="#549bcc"><title>HF_severity 0.696, eGFR 67.5: 0.657</title></rect>
<rect x="341.7" y="123.3" width="8.633" height="8.633" fill="#599fce"><title>HF_severity 0.87, eGFR 67.5: 0.642</title></rect>
<rect x="350" y="123.3" width="8.633" height="8.633" fill="#5fa4d1"><title>HF_severity 1.043, eGFR 67.5: 0.627</title></rect>
<rect x="358.3" y="123.3" width="8.633" height="8.633" fill="#64a8d3"><title>HF_severity 1.217, eGFR 67.5: 0.612</title></rect>
<rect x="366.7" y="123.3" width="8.633" height="8.633" fill="#6aadd5"><title>HF_severity 1.391, eGFR 67.5: 0.597</title></rect>
<rect x="375" y="123.3" width="8.633" height="8.633" fill="#70b1d7"><title>HF_severity 1.565, eGFR 67.5: 0.582</title></rect>
<rect x="383.3" y="123.3" width="8.633" height="8.633" fill="#77b4d9"><title>HF_severity 1.739, eGFR 67.5: 0.567</title></rect>
<rect x="391.7" y="123.3" width="8.633" height="8.633" fill="#7eb7db"><title>HF_severity 1.913, eGFR 67.5: 0.552</title></rect>
<rect x="400" y="123.3" width="8.633" height="8.633" fill="#87bcde"><title>HF_severity 2.087, eGFR 67.5: 0.531</title></rect>
<rect x="408.3" y="123.3" width="8.633" height="8.633" fill="#93c2e1"><title>HF_severity 2.261, eGFR 67.5: 0.505</title></rect>
<rect x="416.7" y="123.3" width="8.633" height="8.633" fill="#9ec7e4"><title>HF_severity 2.435, eGFR 67.5: 0.479</title></rect>
<rect x="425" y="123.3" width="8.633" height="8.633" fill="#aacde7"><title>HF_severity 2.609, eGFR 67.5: 0.453</title></rect>
<rect x="433.3" y="123.3" width="8.633" height="8.633" fill="#b6d3ea"><title>HF_severity 2.783, eGFR 67.5: 0.427</title></rect>
<rect x="441.7" y="123.3" width="8.633" height="8.633" fill="#c1d9ee"><title>HF_severity 2.957, eGFR 67.5: 0.401</title></rect>
<rect x="450" y="123.3" width="8.633" height="8.633" fill="#c1d8ee"><title>HF_severity 3.13, eGFR 67.5: 0.402</title></rect>
<rect x="458.3" y="123.3" width="8.633" height="8.633" fill="#bcd6ec"><title>HF_severity 3.304, eGFR 67.5: 0.412</title></rect>
<rect x="466.7" y="123.3" width="8.633" height="8.633" fill="#b8d4eb"><title>HF_severity 3.478, eGFR 67.5: 0.422</title></rect>
<rect x="475" y="123.3" width="8.633" height="8.633" fill="#b3d2ea"><title>HF_severity 3.652, eGFR 67.5: 0.432</title></rect>
<rect x="483.3" y="123.3" width="8.633" height="8.633" fill="#afd0e9"><title>HF_severity 3.826, eGFR 67.5: 0.442</title></rect>
<rect x="491.7" y="123.3" width="8.633" height="8.633" fill="#aacde7"><title>HF_severity 4, eGFR 67.5: 0.452</title></rect>
<rect x="300" y="131.7" width="8.633" height="8.633" fill="#4b93c8"><title>HF_severity 0, eGFR 62.5: 0.682</title></rect>
<rect x="308.3" y="131.7" width="8.633" height="8.633" fill="#5098ca"><title>HF_severity 0.174, eGFR 62.5: 0.668</title></rect>
<rect x="316.7" y="131.7" width="8.633" height="8.633" fill="#559ccc"><title>HF_severity 0.348, eGFR 62.5: 0.653</title></rect>
<rect x="325" y="131.7" width="8.633" height="8.633" fill="#5aa0cf"><title>HF_severity 0.522, eGFR 62.5: 0.639</title></rect>
<rect x="333.3" y="131.7" width="8.633" height="8.633" fill="#60a5d1"><title>HF_severity 0.696, eGFR 62.5: 0.625</title></rect>
<rect x="341.7" y="131.7" width="8.633" height="8.633" fill="#65a9d3"><title>HF_severity 0.87, eGFR 62.5: 0.61</title></rect>
<rect x="350" y="131.7" width="8.633" height="8.633" fill="#6aadd6"><title>HF_severity 1.043, eGFR 62.5: 0.596</title></rect>
<rect x="358.3" y="131.7" width="8.633" height="8.633" fill="#70b1d7"><title>HF_severity 1.217, eGFR 62.5: 0.582</title></rect>
<rect x="366.7" y="131.7" width="8.633" height="8.633" fill="#77b4d9"><title>HF_severity 1.391, eGFR 62.5: 0.567</title></rect>
<rect x="375" y="131.7" width="8.633" height="8.633" fill="#7db7db"><title>HF_severity 1.565, eGFR 62.5: 0.553</title></rect>
<rect x="383.3" y="131.7" width="8.633" height="8.633" fill="#84badd"><title>HF_severity 1.739, eGFR 62.5: 0.539</title></rect>
<rect x="391.7" y="131.7" width="8.633" height="8.633" fill="#8abdde"><title>HF_severity 1.913, eGFR 62.5: 0.524</title></rect>
<rect x="400" y="131.7" width="8.633" height="8.633" fill="#93c2e1"><title>HF_severity 2.087, eGFR 62.5: 0.505</title></rect>
<rect x="408.3" y="131.7" width="8.633" height="8.633" fill="#9ec7e4"><title>HF_severity 2.261, eGFR 62.5: 0.48</title></rect>
<rect x="416.7" y="131.7" width="8.633" height="8.633" fill="#a9cde7"><title>HF_severity 2.435, eGFR 62.5: 0.455</title></rect>
<rect x="425" y="131.7" width="8.633" height="8.633" fill="#b4d2ea"><title>HF_severity 2.609, eGFR 62.5: 0.43</title></rect>
<rect x="433.3" y="131.7" width="8.633" height="8.633" fill="#bfd8ed"><title>HF_severity 2.783, eGFR 62.5: 0.405</title></rect>
<rect x="441.7" y="131.7" width="8.633" height="8.633" fill="#c8dcf0"><title>HF_severity 2.957, eGFR 62.5: 0.381</title></rect>
<rect x="450" y="131.7" width="8.633" height="8.633" fill="#c8dcf0"><title>HF_severity 3.13, eGFR 62.5: 0.382</title></rect>
<rect x="458.3" y="131.7" width="8.633" height="8.633" fill="#c5dbef"><title>HF_severity 3.304, eGFR 62.5: 0.391</title></rect>
<rect x="466.7" y="131.7" width="8.633" height="8.633" fill="#c1d9ee"><title>HF_severity 3.478, eGFR 62.5: 0.401</title></rect>
<rect x="475" y="131.7" width="8.633" height="8.633" fill="#bdd7ed"><title>HF_severity 3.652, eGFR 62.5: 0.41</title></rect>
<rect x="483.3" y="131.7" width="8.633" height="8.633" fill="#b9d4eb"><title>HF_severity 3.826, eGFR 62.5: 0.42</title></rect>
<rect x="491.7" y="131.7" width="8.633" height="8.633" fill="#b4d2ea"><title>HF_severity 4, eGFR 62.5: 0.429</title></rect>
<rect x="300" y="140" width="8.633" height="8.633" fill="#589ecd"><title>HF_severity 0, eGFR 57.5: 0.646</title></rect>
<rect x="308.3" y="140" width="8.633" height="8.633" fill="#5da2d0"><title>HF_severity 0.174, eGFR 57.5: 0.633</title></rect>
<rect x="316.7" y="140" width="8.633" height="8.633" fill="#62a6d2"><title>HF_severity 0.348, eGFR 57.5: 0.619</title></rect>
<rect x="325" y="140" width="8.633" height="8.633" fill="#67aad4"><title>HF_severity 0.522, eGFR 57.5: 0.605</title></rect>
<rect x="333.3" y="140" width="8.633" height="8.633" fill="#6caed6"><title>HF_severity 0.696, eGFR 57.5: 0.592</title></rect>
<rect x="341.7" y="140" width="8.633" height="8.633" fill="#72b1d8"><title>HF_severity 0.87, eGFR 57.5: 0.578</title></rect>
<rect x="350" y="140" width="8.633" height="8.633" fill="#78b4da"><title>HF_severity 1.043, eGFR 57.5: 0.565</title></rect>
<rect x="358.3" y="140" width="8.633" height="8.633" fill="#7eb7db"><title>HF_severity 1.217, eGFR 57.5: 0.551</title></rect>
<rect x="366.7" y="140" width="8.633" height="8.633" fill="#84badd"><title>HF_severity 1.391, eGFR 57.5: 0.537</title></rect>
<rect x="375" y="140" width="8.633" height="8.633" fill="#8abddf"><title>HF_severity 1.565, eGFR 57.5: 0.524</title></rect>
<rect x="383.3" y="140" width="8.633" height="8.633" fill="#90c0e0"><title>HF_severity 1.739, eGFR 57.5: 0.51</title></rect>
<rect x="391.7" y="140" width="8.633" height="8.633" fill="#96c3e2"><title>HF_severity 1.913, eGFR 57.5: 0.497</title></rect>
<rect x="400" y="140" width="8.633" height="8.633" fill="#9fc7e4"><title>HF_severity 2.087, eGFR 57.5: 0.478</title></rect>
<rect x="408.3" y="140" width="8.633" height="8.633" fill="#a9cde7"><title>HF_severity 2.261, eGFR 57.5: 0.455</title></rect>
<rect x="416.7" y="140" width="8.633" height="8.633" fill="#b4d2ea"><title>HF_severity 2.435, eGFR 57.5: 0.431</title></rect>
<rect x="425" y="140" width="8.633" height="8.633" fill="#bed7ed"><title>HF_severity 2.609, eGFR 57.5: 0.408</title></rect>
<rect x="433.3" y="140" width="8.633" height="8.633" fill="#c7dcef"><title>HF_severity 2.783, eGFR 57.5: 0.384</title></rect>
<rect x="441.7" y="140" width="8.633" height="8.633" fill="#cde0f1"><title>HF_severity 2.957, eGFR 57.5: 0.361</title></rect>
<rect x="450" y="140" width="8.633" height="8.633" fill="#cddff1"><title>HF_severity 3.13, eGFR 57.5: 0.362</title></rect>
<rect x="458.3" y="140" width="8.633" height="8.633" fill="#cbdef1"><title>HF_severity 3.304, eGFR 57.5: 0.371</title></rect>
<rect x="466.7" y="140" width="8.633" height="8.633" fill="#c8ddf0"><title>HF_severity 3.478, eGFR 57.5: 0.38</title></rect>
<rect x="475" y="140" width="8.633" height="8.633" fill="#c6dbef"><title>HF_severity 3.652, eGFR 57.5: 0.389</title></rect>
<rect x="483.3" y="140" width="8.633" height="8.633" fill="#c2d9ee"><title>HF_severity 3.826, eGFR 57.5: 0.398</title></rect>
<rect x="491.7" y="140" width="8.633" height="8.633" fill="#bed7ed"><title>HF_severity 4, eGFR 57.5: 0.407</title></rect>
<rect x="300" y="148.3" width="8.633" height="8.633" fill="#65a9d3"><title>HF_severity 0, eGFR 52.5: 0.61</title></rect>
<rect x="308.3" y="148.3" width="8.633" height="8.633" fill="#69add5"><title>HF_severity 0.174, eGFR 52.5: 0.598</title></rect>
<rect x="316.7" y="148.3" width="8.633" height="8.633" fill="#6fb0d7"><title>HF_severity 0.348, eGFR 52.5: 0.585</title></rect>
<rect x="325" y="148.3" width="8.633" height="8.633" fill="#75b3d9"><title>HF_severity 0.522, eGFR 52.5: 0.572</title></rect>
<rect x="333.3" y="148.3" width="8.633" height="8.633" fill="#7ab6da"><title>HF_severity 0.696, eGFR 52.5: 0.559</title></rect>
<rect x="341.7" y="148.3" width="8.633" height="8.633" fill="#80b8dc"><title>HF_severity 0.87, eGFR 52.5: 0.546</title></rect>
<rect x="350" y="148.3" width="8.633" height="8.633" fill="#86bbdd"><title>HF_severity 1.043, eGFR 52.5: 0.533</title></rect>
<rect x="358.3" y="148.3" width="8.633" height="8.633" fill="#8cbedf"><title>HF_severity 1.217, eGFR 52.5: 0.521</title></rect>
<rect x="366.7" y="148.3" width="8.633" height="8.633" fill="#91c1e1"><title>HF_severity 1.391, eGFR 52.5: 0.508</title></rect>
<rect x="375" y="148.3" width="8.633" height="8.633" fill="#97c4e2"><title>HF_severity 1.565, eGFR 52.5: 0.495</title></rect>
<rect x="383.3" y="148.3" width="8.633" height="8.633" fill="#9dc7e4"><title>HF_severity 1.739, eGFR 52.5: 0.482</title></rect>
<rect x="391.7" y="148.3" width="8.633" height="8.633" fill="#a3c9e5"><title>HF_severity 1.913, eGFR 52.5: 0.469</title></rect>
<rect x="400" y="148.3" width="8.633" height="8.633" fill="#aacde7"><title>HF_severity 2.087, eGFR 52.5: 0.452</title></rect>
<rect x="408.3" y="148.3" width="8.633" height="8.633" fill="#b4d2ea"><title>HF_severity 2.261, eGFR 52.5: 0.43</title></rect>
<rect x="416.7" y="148.3" width="8.633" height="8.633" fill="#bed7ed"><title>HF_severity 2.435, eGFR 52.5: 0.407</title></rect>
<rect x="425" y="148.3" width="8.633" height="8.633" fill="#c7dcef"><title>HF_severity 2.609, eGFR 52.5: 0.385</title></rect>
<rect x="433.3" y="148.3" width="8.633" height="8.633" fill="#cddff1"><title>HF_severity 2.783, eGFR 52.5: 0.363</title></rect>
<rect x="441.7" y="148.3" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 2.957, eGFR 52.5: 0.341</title></rect>
<rect x="450" y="148.3" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 3.13, eGFR 52.5: 0.342</title></rect>
<rect x="458.3" y="148.3" width="8.633" height="8.633" fill="#d0e1f2"><title>HF_severity 3.304, eGFR 52.5: 0.35</title></rect>
<rect x="466.7" y="148.3" width="8.633" height="8.633" fill="#cee0f1"><title>HF_severity 3.478, eGFR 52.5: 0.359</title></rect>
<rect x="475" y="148.3" width="8.633" height="8.633" fill="#cbdff1"><title>HF_severity 3.652, eGFR 52.5: 0.367</title></rect>
<rect x="483.3" y="148.3" width="8.633" height="8.633" fill="#c9ddf0"><title>HF_severity 3.826, eGFR 52.5: 0.376</title></rect>
<rect x="491.7" y="148.3" width="8.633" height="8.633" fill="#c7dcef"><title>HF_severity 4, eGFR 52.5: 0.384</title></rect>
<rect x="300" y="156.7" width="8.633" height="8.633" fill="#73b2d8"><title>HF_severity 0, eGFR 47.5: 0.575</title></rect>
<rect x="308.3" y="156.7" width="8.633" height="8.633" fill="#79b5da"><title>HF_severity 0.174, eGFR 47.5: 0.563</title></rect>
<rect x="316.7" y="156.7" width="8.633" height="8.633" fill="#7eb7db"><title>HF_severity 0.348, eGFR 47.5: 0.55</title></rect>
<rect x="325" y="156.7" width="8.633" height="8.633" fill="#84badd"><title>HF_severity 0.522, eGFR 47.5: 0.538</title></rect>
<rect x="333.3" y="156.7" width="8.633" height="8.633" fill="#89bdde"><title>HF_severity 0.696, eGFR 47.5: 0.526</title></rect>
<rect x="341.7" y="156.7" width="8.633" height="8.633" fill="#8ec0e0"><title>HF_severity 0.87, eGFR 47.5: 0.514</title></rect>
<rect x="350" y="156.7" width="8.633" height="8.633" fill="#94c2e1"><title>HF_severity 1.043, eGFR 47.5: 0.502</title></rect>
<rect x="358.3" y="156.7" width="8.633" height="8.633" fill="#99c5e3"><title>HF_severity 1.217, eGFR 47.5: 0.49</title></rect>
<rect x="366.7" y="156.7" width="8.633" height="8.633" fill="#9fc8e4"><title>HF_severity 1.391, eGFR 47.5: 0.478</title></rect>
<rect x="375" y="156.7" width="8.633" height="8.633" fill="#a4cae6"><title>HF_severity 1.565, eGFR 47.5: 0.466</title></rect>
<rect x="383.3" y="156.7" width="8.633" height="8.633" fill="#a9cde7"><title>HF_severity 1.739, eGFR 47.5: 0.454</title></rect>
<rect x="391.7" y="156.7" width="8.633" height="8.633" fill="#afd0e9"><title>HF_severity 1.913, eGFR 47.5: 0.442</title></rect>
<rect x="400" y="156.7" width="8.633" height="8.633" fill="#b6d3eb"><title>HF_severity 2.087, eGFR 47.5: 0.425</title></rect>
<rect x="408.3" y="156.7" width="8.633" height="8.633" fill="#c0d8ed"><title>HF_severity 2.261, eGFR 47.5: 0.404</title></rect>
<rect x="416.7" y="156.7" width="8.633" height="8.633" fill="#c8dcf0"><title>HF_severity 2.435, eGFR 47.5: 0.383</title></rect>
<rect x="425" y="156.7" width="8.633" height="8.633" fill="#cddff1"><title>HF_severity 2.609, eGFR 47.5: 0.363</title></rect>
<rect x="433.3" y="156.7" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 2.783, eGFR 47.5: 0.342</title></rect>
<rect x="441.7" y="156.7" width="8.633" height="8.633" fill="#d7e6f4"><title>HF_severity 2.957, eGFR 47.5: 0.321</title></rect>
<rect x="450" y="156.7" width="8.633" height="8.633" fill="#d6e6f4"><title>HF_severity 3.13, eGFR 47.5: 0.322</title></rect>
<rect x="458.3" y="156.7" width="8.633" height="8.633" fill="#d5e5f4"><title>HF_severity 3.304, eGFR 47.5: 0.33</title></rect>
<rect x="466.7" y="156.7" width="8.633" height="8.633" fill="#d3e3f3"><title>HF_severity 3.478, eGFR 47.5: 0.338</title></rect>
<rect x="475" y="156.7" width="8.633" height="8.633" fill="#d1e2f2"><title>HF_severity 3.652, eGFR 47.5: 0.346</title></rect>
<rect x="483.3" y="156.7" width="8.633" height="8.633" fill="#cfe1f2"><title>HF_severity 3.826, eGFR 47.5: 0.354</title></rect>
<rect x="491.7" y="156.7" width="8.633" height="8.633" fill="#cddff1"><title>HF_severity 4, eGFR 47.5: 0.362</title></rect>
<rect x="300" y="165" width="8.633" height="8.633" fill="#82badc"><title>HF_severity 0, eGFR 42.5: 0.541</title></rect>
<rect x="308.3" y="165" width="8.633" height="8.633" fill="#87bcde"><title>HF_severity 0.174, eGFR 42.5: 0.53</title></rect>
<rect x="316.7" y="165" width="8.633" height="8.633" fill="#8cbfdf"><title>HF_severity 0.348, eGFR 42.5: 0.519</title></rect>
<rect x="325" y="165" width="8.633" height="8.633" fill="#92c1e1"><title>HF_severity 0.522, eGFR 42.5: 0.507</title></rect>
<rect x="333.3" y="165" width="8.633" height="8.633" fill="#97c4e2"><title>HF_severity 0.696, eGFR 42.5: 0.496</title></rect>
<rect x="341.7" y="165" width="8.633" height="8.633" fill="#9cc6e3"><title>HF_severity 0.87, eGFR 42.5: 0.484</title></rect>
<rect x="350" y="165" width="8.633" height="8.633" fill="#a1c9e5"><title>HF_severity 1.043, eGFR 42.5: 0.473</title></rect>
<rect x="358.3" y="165" width="8.633" height="8.633" fill="#a6cbe6"><title>HF_severity 1.217, eGFR 42.5: 0.462</title></rect>
<rect x="366.7" y="165" width="8.633" height="8.633" fill="#abcee8"><title>HF_severity 1.391, eGFR 42.5: 0.45</title></rect>
<rect x="375" y="165" width="8.633" height="8.633" fill="#b0d0e9"><title>HF_severity 1.565, eGFR 42.5: 0.439</title></rect>
<rect x="383.3" y="165" width="8.633" height="8.633" fill="#b5d3ea"><title>HF_severity 1.739, eGFR 42.5: 0.428</title></rect>
<rect x="391.7" y="165" width="8.633" height="8.633" fill="#bad5ec"><title>HF_severity 1.913, eGFR 42.5: 0.416</title></rect>
<rect x="400" y="165" width="8.633" height="8.633" fill="#c1d9ee"><title>HF_severity 2.087, eGFR 42.5: 0.401</title></rect>
<rect x="408.3" y="165" width="8.633" height="8.633" fill="#c8dcf0"><title>HF_severity 2.261, eGFR 42.5: 0.381</title></rect>
<rect x="416.7" y="165" width="8.633" height="8.633" fill="#cde0f1"><title>HF_severity 2.435, eGFR 42.5: 0.361</title></rect>
<rect x="425" y="165" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 2.609, eGFR 42.5: 0.342</title></rect>
<rect x="433.3" y="165" width="8.633" height="8.633" fill="#d6e6f4"><title>HF_severity 2.783, eGFR 42.5: 0.322</title></rect>
<rect x="441.7" y="165" width="8.633" height="8.633" fill="#dbe9f6"><title>HF_severity 2.957, eGFR 42.5: 0.302</title></rect>
<rect x="450" y="165" width="8.633" height="8.633" fill="#dbe9f6"><title>HF_severity 3.13, eGFR 42.5: 0.303</title></rect>
<rect x="458.3" y="165" width="8.633" height="8.633" fill="#d9e8f5"><title>HF_severity 3.304, eGFR 42.5: 0.311</title></rect>
<rect x="466.7" y="165" width="8.633" height="8.633" fill="#d7e6f5"><title>HF_severity 3.478, eGFR 42.5: 0.318</title></rect>
<rect x="475" y="165" width="8.633" height="8.633" fill="#d5e5f4"><title>HF_severity 3.652, eGFR 42.5: 0.326</title></rect>
<rect x="483.3" y="165" width="8.633" height="8.633" fill="#d4e4f3"><title>HF_severity 3.826, eGFR 42.5: 0.333</title></rect>
<rect x="491.7" y="165" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 4, eGFR 42.5: 0.341</title></rect>
<rect x="300" y="173.3" width="8.633" height="8.633" fill="#90c0e0"><title>HF_severity 0, eGFR 37.5: 0.511</title></rect>
<rect x="308.3" y="173.3" width="8.633" height="8.633" fill="#95c3e1"><title>HF_severity 0.174, eGFR 37.5: 0.5</title></rect>
<rect x="316.7" y="173.3" width="8.633" height="8.633" fill="#9ac5e3"><title>HF_severity 0.348, eGFR 37.5: 0.489</title></rect>
<rect x="325" y="173.3" width="8.633" height="8.633" fill="#9ec7e4"><title>HF_severity 0.522, eGFR 37.5: 0.478</title></rect>
<rect x="333.3" y="173.3" width="8.633" height="8.633" fill="#a3cae5"><title>HF_severity 0.696, eGFR 37.5: 0.468</title></rect>
<rect x="341.7" y="173.3" width="8.633" height="8.633" fill="#a8cce7"><title>HF_severity 0.87, eGFR 37.5: 0.457</title></rect>
<rect x="350" y="173.3" width="8.633" height="8.633" fill="#adcfe8"><title>HF_severity 1.043, eGFR 37.5: 0.446</title></rect>
<rect x="358.3" y="173.3" width="8.633" height="8.633" fill="#b2d1e9"><title>HF_severity 1.217, eGFR 37.5: 0.435</title></rect>
<rect x="366.7" y="173.3" width="8.633" height="8.633" fill="#b6d3eb"><title>HF_severity 1.391, eGFR 37.5: 0.425</title></rect>
<rect x="375" y="173.3" width="8.633" height="8.633" fill="#bbd6ec"><title>HF_severity 1.565, eGFR 37.5: 0.414</title></rect>
<rect x="383.3" y="173.3" width="8.633" height="8.633" fill="#c0d8ed"><title>HF_severity 1.739, eGFR 37.5: 0.403</title></rect>
<rect x="391.7" y="173.3" width="8.633" height="8.633" fill="#c5daef"><title>HF_severity 1.913, eGFR 37.5: 0.393</title></rect>
<rect x="400" y="173.3" width="8.633" height="8.633" fill="#c9ddf0"><title>HF_severity 2.087, eGFR 37.5: 0.378</title></rect>
<rect x="408.3" y="173.3" width="8.633" height="8.633" fill="#cde0f1"><title>HF_severity 2.261, eGFR 37.5: 0.359</title></rect>
<rect x="416.7" y="173.3" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 2.435, eGFR 37.5: 0.341</title></rect>
<rect x="425" y="173.3" width="8.633" height="8.633" fill="#d6e6f4"><title>HF_severity 2.609, eGFR 37.5: 0.322</title></rect>
<rect x="433.3" y="173.3" width="8.633" height="8.633" fill="#dbe9f6"><title>HF_severity 2.783, eGFR 37.5: 0.304</title></rect>
<rect x="441.7" y="173.3" width="8.633" height="8.633" fill="#dfecf7"><title>HF_severity 2.957, eGFR 37.5: 0.285</title></rect>
<rect x="450" y="173.3" width="8.633" height="8.633" fill="#dfebf7"><title>HF_severity 3.13, eGFR 37.5: 0.286</title></rect>
<rect x="458.3" y="173.3" width="8.633" height="8.633" fill="#ddeaf7"><title>HF_severity 3.304, eGFR 37.5: 0.293</title></rect>
<rect x="466.7" y="173.3" width="8.633" height="8.633" fill="#dce9f6"><title>HF_severity 3.478, eGFR 37.5: 0.3</title></rect>
<rect x="475" y="173.3" width="8.633" height="8.633" fill="#dae8f6"><title>HF_severity 3.652, eGFR 37.5: 0.307</title></rect>
<rect x="483.3" y="173.3" width="8.633" height="8.633" fill="#d8e7f5"><title>HF_severity 3.826, eGFR 37.5: 0.314</title></rect>
<rect x="491.7" y="173.3" width="8.633" height="8.633" fill="#d6e6f4"><title>HF_severity 4, eGFR 37.5: 0.322</title></rect>
<rect x="300" y="181.7" width="8.633" height="8.633" fill="#9ec7e4"><title>HF_severity 0, eGFR 32.5: 0.48</title></rect>
<rect x="308.3" y="181.7" width="8.633" height="8.633" fill="#a2c9e5"><title>HF_severity 0.174, eGFR 32.5: 0.47</title></rect>
<rect x="316.7" y="181.7" width="8.633" height="8.633" fill="#a7cce6"><title>HF_severity 0.348, eGFR 32.5: 0.46</title></rect>
<rect x="325" y="181.7" width="8.633" height="8.633" fill="#abcee8"><title>HF_severity 0.522, eGFR 32.5: 0.45</title></rect>
<rect x="333.3" y="181.7" width="8.633" height="8.633" fill="#b0d0e9"><title>HF_severity 0.696, eGFR 32.5: 0.439</title></rect>
<rect x="341.7" y="181.7" width="8.633" height="8.633" fill="#b4d2ea"><title>HF_severity 0.87, eGFR 32.5: 0.429</title></rect>
<rect x="350" y="181.7" width="8.633" height="8.633" fill="#b9d5eb"><title>HF_severity 1.043, eGFR 32.5: 0.419</title></rect>
<rect x="358.3" y="181.7" width="8.633" height="8.633" fill="#bdd7ed"><title>HF_severity 1.217, eGFR 32.5: 0.409</title></rect>
<rect x="366.7" y="181.7" width="8.633" height="8.633" fill="#c2d9ee"><title>HF_severity 1.391, eGFR 32.5: 0.399</title></rect>
<rect x="375" y="181.7" width="8.633" height="8.633" fill="#c6dbef"><title>HF_severity 1.565, eGFR 32.5: 0.389</title></rect>
<rect x="383.3" y="181.7" width="8.633" height="8.633" fill="#c9ddf0"><title>HF_severity 1.739, eGFR 32.5: 0.379</title></rect>
<rect x="391.7" y="181.7" width="8.633" height="8.633" fill="#cbdef1"><title>HF_severity 1.913, eGFR 32.5: 0.369</title></rect>
<rect x="400" y="181.7" width="8.633" height="8.633" fill="#cee0f2"><title>HF_severity 2.087, eGFR 32.5: 0.355</title></rect>
<rect x="408.3" y="181.7" width="8.633" height="8.633" fill="#d3e3f3"><title>HF_severity 2.261, eGFR 32.5: 0.338</title></rect>
<rect x="416.7" y="181.7" width="8.633" height="8.633" fill="#d7e6f4"><title>HF_severity 2.435, eGFR 32.5: 0.32</title></rect>
<rect x="425" y="181.7" width="8.633" height="8.633" fill="#dbe9f6"><title>HF_severity 2.609, eGFR 32.5: 0.303</title></rect>
<rect x="433.3" y="181.7" width="8.633" height="8.633" fill="#dfebf7"><title>HF_severity 2.783, eGFR 32.5: 0.285</title></rect>
<rect x="441.7" y="181.7" width="8.633" height="8.633" fill="#e3eef9"><title>HF_severity 2.957, eGFR 32.5: 0.268</title></rect>
<rect x="450" y="181.7" width="8.633" height="8.633" fill="#e3eef9"><title>HF_severity 3.13, eGFR 32.5: 0.269</title></rect>
<rect x="458.3" y="181.7" width="8.633" height="8.633" fill="#e2edf8"><title>HF_severity 3.304, eGFR 32.5: 0.275</title></rect>
<rect x="466.7" y="181.7" width="8.633" height="8.633" fill="#e0ecf7"><title>HF_severity 3.478, eGFR 32.5: 0.282</title></rect>
<rect x="475" y="181.7" width="8.633" height="8.633" fill="#deebf7"><title>HF_severity 3.652, eGFR 32.5: 0.289</title></rect>
<rect x="483.3" y="181.7" width="8.633" height="8.633" fill="#ddeaf6"><title>HF_severity 3.826, eGFR 32.5: 0.295</title></rect>
<rect x="491.7" y="181.7" width="8.633" height="8.633" fill="#dbe9f6"><title>HF_severity 4, eGFR 32.5: 0.302</title></rect>
<rect x="300" y="190" width="8.633" height="8.633" fill="#accee8"><title>HF_severity 0, eGFR 27.5: 0.449</title></rect>
<rect x="308.3" y="190" width="8.633" height="8.633" fill="#b0d0e9"><title>HF_severity 0.174, eGFR 27.5: 0.44</title></rect>
<rect x="316.7" y="190" width="8.633" height="8.633" fill="#b4d2ea"><title>HF_severity 0.348, eGFR 27.5: 0.43</title></rect>
<rect x="325" y="190" width="8.633" height="8.633" fill="#b8d4eb"><title>HF_severity 0.522, eGFR 27.5: 0.421</title></rect>
<rect x="333.3" y="190" width="8.633" height="8.633" fill="#bcd6ec"><title>HF_severity 0.696, eGFR 27.5: 0.411</title></rect>
<rect x="341.7" y="190" width="8.633" height="8.633" fill="#c1d8ee"><title>HF_severity 0.87, eGFR 27.5: 0.402</title></rect>
<rect x="350" y="190" width="8.633" height="8.633" fill="#c5daef"><title>HF_severity 1.043, eGFR 27.5: 0.392</title></rect>
<rect x="358.3" y="190" width="8.633" height="8.633" fill="#c8dcf0"><title>HF_severity 1.217, eGFR 27.5: 0.383</title></rect>
<rect x="366.7" y="190" width="8.633" height="8.633" fill="#cadef0"><title>HF_severity 1.391, eGFR 27.5: 0.374</title></rect>
<rect x="375" y="190" width="8.633" height="8.633" fill="#ccdff1"><title>HF_severity 1.565, eGFR 27.5: 0.364</title></rect>
<rect x="383.3" y="190" width="8.633" height="8.633" fill="#cfe1f2"><title>HF_severity 1.739, eGFR 27.5: 0.355</title></rect>
<rect x="391.7" y="190" width="8.633" height="8.633" fill="#d1e2f3"><title>HF_severity 1.913, eGFR 27.5: 0.345</title></rect>
<rect x="400" y="190" width="8.633" height="8.633" fill="#d4e4f4"><title>HF_severity 2.087, eGFR 27.5: 0.332</title></rect>
<rect x="408.3" y="190" width="8.633" height="8.633" fill="#d8e7f5"><title>HF_severity 2.261, eGFR 27.5: 0.316</title></rect>
<rect x="416.7" y="190" width="8.633" height="8.633" fill="#dce9f6"><title>HF_severity 2.435, eGFR 27.5: 0.3</title></rect>
<rect x="425" y="190" width="8.633" height="8.633" fill="#e0ecf7"><title>HF_severity 2.609, eGFR 27.5: 0.283</title></rect>
<rect x="433.3" y="190" width="8.633" height="8.633" fill="#e4eef9"><title>HF_severity 2.783, eGFR 27.5: 0.267</title></rect>
<rect x="441.7" y="190" width="8.633" height="8.633" fill="#e8f1fa"><title>HF_severity 2.957, eGFR 27.5: 0.251</title></rect>
<rect x="450" y="190" width="8.633" height="8.633" fill="#e7f1fa"><title>HF_severity 3.13, eGFR 27.5: 0.251</title></rect>
<rect x="458.3" y="190" width="8.633" height="8.633" fill="#e6f0f9"><title>HF_severity 3.304, eGFR 27.5: 0.258</title></rect>
<rect x="466.7" y="190" width="8.633" height="8.633" fill="#e4eff9"><title>HF_severity 3.478, eGFR 27.5: 0.264</title></rect>
<rect x="475" y="190" width="8.633" height="8.633" fill="#e3eef8"><title>HF_severity 3.652, eGFR 27.5: 0.27</title></rect>
<rect x="483.3" y="190" width="8.633" height="8.633" fill="#e1edf8"><title>HF_severity 3.826, eGFR 27.5: 0.277</title></rect>
<rect x="491.7" y="190" width="8.633" height="8.633" fill="#e0ecf7"><title>HF_severity 4, eGFR 27.5: 0.283</title></rect>
<rect x="300" y="198.3" width="8.633" height="8.633" fill="#b9d5ec"><title>HF_severity 0, eGFR 22.5: 0.418</title></rect>
<rect x="308.3" y="198.3" width="8.633" height="8.633" fill="#bdd7ed"><title>HF_severity 0.174, eGFR 22.5: 0.41</title></rect>
<rect x="316.7" y="198.3" width="8.633" height="8.633" fill="#c1d9ee"><title>HF_severity 0.348, eGFR 22.5: 0.401</title></rect>
<rect x="325" y="198.3" width="8.633" height="8.633" fill="#c5dbef"><title>HF_severity 0.522, eGFR 22.5: 0.392</title></rect>
<rect x="333.3" y="198.3" width="8.633" height="8.633" fill="#c8dcf0"><title>HF_severity 0.696, eGFR 22.5: 0.383</title></rect>
<rect x="341.7" y="198.3" width="8.633" height="8.633" fill="#caddf0"><title>HF_severity 0.87, eGFR 22.5: 0.374</title></rect>
<rect x="350" y="198.3" width="8.633" height="8.633" fill="#ccdff1"><title>HF_severity 1.043, eGFR 22.5: 0.366</title></rect>
<rect x="358.3" y="198.3" width="8.633" height="8.633" fill="#cee0f2"><title>HF_severity 1.217, eGFR 22.5: 0.357</title></rect>
<rect x="366.7" y="198.3" width="8.633" height="8.633" fill="#d0e2f2"><title>HF_severity 1.391, eGFR 22.5: 0.348</title></rect>
<rect x="375" y="198.3" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 1.565, eGFR 22.5: 0.339</title></rect>
<rect x="383.3" y="198.3" width="8.633" height="8.633" fill="#d4e4f4"><title>HF_severity 1.739, eGFR 22.5: 0.33</title></rect>
<rect x="391.7" y="198.3" width="8.633" height="8.633" fill="#d6e6f4"><title>HF_severity 1.913, eGFR 22.5: 0.322</title></rect>
<rect x="400" y="198.3" width="8.633" height="8.633" fill="#d9e8f5"><title>HF_severity 2.087, eGFR 22.5: 0.31</title></rect>
<rect x="408.3" y="198.3" width="8.633" height="8.633" fill="#ddeaf7"><title>HF_severity 2.261, eGFR 22.5: 0.294</title></rect>
<rect x="416.7" y="198.3" width="8.633" height="8.633" fill="#e1ecf8"><title>HF_severity 2.435, eGFR 22.5: 0.279</title></rect>
<rect x="425" y="198.3" width="8.633" height="8.633" fill="#e4eff9"><title>HF_severity 2.609, eGFR 22.5: 0.264</title></rect>
<rect x="433.3" y="198.3" width="8.633" height="8.633" fill="#e8f1fa"><title>HF_severity 2.783, eGFR 22.5: 0.249</title></rect>
<rect x="441.7" y="198.3" width="8.633" height="8.633" fill="#ecf4fb"><title>HF_severity 2.957, eGFR 22.5: 0.234</title></rect>
<rect x="450" y="198.3" width="8.633" height="8.633" fill="#ecf4fb"><title>HF_severity 3.13, eGFR 22.5: 0.234</title></rect>
<rect x="458.3" y="198.3" width="8.633" height="8.633" fill="#eaf3fb"><title>HF_severity 3.304, eGFR 22.5: 0.24</title></rect>
<rect x="466.7" y="198.3" width="8.633" height="8.633" fill="#e9f2fa"><title>HF_severity 3.478, eGFR 22.5: 0.246</title></rect>
<rect x="475" y="198.3" width="8.633" height="8.633" fill="#e7f1fa"><title>HF_severity 3.652, eGFR 22.5: 0.252</title></rect>
<rect x="483.3" y="198.3" width="8.633" height="8.633" fill="#e6f0f9"><title>HF_severity 3.826, eGFR 22.5: 0.258</title></rect>
<rect x="491.7" y="198.3" width="8.633" height="8.633" fill="#e5eff9"><title>HF_severity 4, eGFR 22.5: 0.263</title></rect>
<rect x="300" y="206.7" width="8.633" height="8.633" fill="#c6dbef"><title>HF_severity 0, eGFR 17.5: 0.39</title></rect>
<rect x="308.3" y="206.7" width="8.633" height="8.633" fill="#c8dcf0"><title>HF_severity 0.174, eGFR 17.5: 0.382</title></rect>
<rect x="316.7" y="206.7" width="8.633" height="8.633" fill="#cadef0"><title>HF_severity 0.348, eGFR 17.5: 0.374</title></rect>
<rect x="325" y="206.7" width="8.633" height="8.633" fill="#ccdff1"><title>HF_severity 0.522, eGFR 17.5: 0.366</title></rect>
<rect x="333.3" y="206.7" width="8.633" height="8.633" fill="#cee0f2"><title>HF_severity 0.696, eGFR 17.5: 0.358</title></rect>
<rect x="341.7" y="206.7" width="8.633" height="8.633" fill="#d0e1f2"><title>HF_severity 0.87, eGFR 17.5: 0.349</title></rect>
<rect x="350" y="206.7" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 1.043, eGFR 17.5: 0.341</title></rect>
<rect x="358.3" y="206.7" width="8.633" height="8.633" fill="#d4e4f3"><title>HF_severity 1.217, eGFR 17.5: 0.333</title></rect>
<rect x="366.7" y="206.7" width="8.633" height="8.633" fill="#d6e5f4"><title>HF_severity 1.391, eGFR 17.5: 0.325</title></rect>
<rect x="375" y="206.7" width="8.633" height="8.633" fill="#d8e7f5"><title>HF_severity 1.565, eGFR 17.5: 0.317</title></rect>
<rect x="383.3" y="206.7" width="8.633" height="8.633" fill="#dae8f5"><title>HF_severity 1.739, eGFR 17.5: 0.308</title></rect>
<rect x="391.7" y="206.7" width="8.633" height="8.633" fill="#dce9f6"><title>HF_severity 1.913, eGFR 17.5: 0.3</title></rect>
<rect x="400" y="206.7" width="8.633" height="8.633" fill="#deebf7"><title>HF_severity 2.087, eGFR 17.5: 0.289</title></rect>
<rect x="408.3" y="206.7" width="8.633" height="8.633" fill="#e2edf8"><title>HF_severity 2.261, eGFR 17.5: 0.275</title></rect>
<rect x="416.7" y="206.7" width="8.633" height="8.633" fill="#e5eff9"><title>HF_severity 2.435, eGFR 17.5: 0.261</title></rect>
<rect x="425" y="206.7" width="8.633" height="8.633" fill="#e9f2fa"><title>HF_severity 2.609, eGFR 17.5: 0.246</title></rect>
<rect x="433.3" y="206.7" width="8.633" height="8.633" fill="#ecf4fb"><title>HF_severity 2.783, eGFR 17.5: 0.232</title></rect>
<rect x="441.7" y="206.7" width="8.633" height="8.633" fill="#eff6fd"><title>HF_severity 2.957, eGFR 17.5: 0.218</title></rect>
<rect x="450" y="206.7" width="8.633" height="8.633" fill="#eff6fc"><title>HF_severity 3.13, eGFR 17.5: 0.218</title></rect>
<rect x="458.3" y="206.7" width="8.633" height="8.633" fill="#eef5fc"><title>HF_severity 3.304, eGFR 17.5: 0.224</title></rect>
<rect x="466.7" y="206.7" width="8.633" height="8.633" fill="#edf4fc"><title>HF_severity 3.478, eGFR 17.5: 0.229</title></rect>
<rect x="475" y="206.7" width="8.633" height="8.633" fill="#ebf3fb"><title>HF_severity 3.652, eGFR 17.5: 0.235</title></rect>
<rect x="483.3" y="206.7" width="8.633" height="8.633" fill="#eaf3fb"><title>HF_severity 3.826, eGFR 17.5: 0.24</title></rect>
<rect x="491.7" y="206.7" width="8.633" height="8.633" fill="#e9f2fa"><title>HF_severity 4, eGFR 17.5: 0.246</title></rect>
<rect x="300" y="215" width="8.633" height="8.633" fill="#cddff1"><title>HF_severity 0, eGFR 12.5: 0.362</title></rect>
<rect x="308.3" y="215" width="8.633" height="8.633" fill="#cee1f2"><title>HF_severity 0.174, eGFR 12.5: 0.355</title></rect>
<rect x="316.7" y="215" width="8.633" height="8.633" fill="#d0e2f2"><title>HF_severity 0.348, eGFR 12.5: 0.347</title></rect>
<rect x="325" y="215" width="8.633" height="8.633" fill="#d2e3f3"><title>HF_severity 0.522, eGFR 12.5: 0.34</title></rect>
<rect x="333.3" y="215" width="8.633" height="8.633" fill="#d4e4f4"><title>HF_severity 0.696, eGFR 12.5: 0.332</title></rect>
<rect x="341.7" y="215" width="8.633" height="8.633" fill="#d6e5f4"><title>HF_severity 0.87, eGFR 12.5: 0.324</title></rect>
<rect x="350" y="215" width="8.633" height="8.633" fill="#d8e7f5"><title>HF_severity 1.043, eGFR 12.5: 0.317</title></rect>
<rect x="358.3" y="215" width="8.633" height="8.633" fill="#dae8f5"><title>HF_severity 1.217, eGFR 12.5: 0.309</title></rect>
<rect x="366.7" y="215" width="8.633" height="8.633" fill="#dbe9f6"><title>HF_severity 1.391, eGFR 12.5: 0.301</title></rect>
<rect x="375" y="215" width="8.633" height="8.633" fill="#ddeaf7"><title>HF_severity 1.565, eGFR 12.5: 0.294</title></rect>
<rect x="383.3" y="215" width="8.633" height="8.633" fill="#dfebf7"><title>HF_severity 1.739, eGFR 12.5: 0.286</title></rect>
<rect x="391.7" y="215" width="8.633" height="8.633" fill="#e1edf8"><title>HF_severity 1.913, eGFR 12.5: 0.279</title></rect>
<rect x="400" y="215" width="8.633" height="8.633" fill="#e3eef9"><title>HF_severity 2.087, eGFR 12.5: 0.268</title></rect>
<rect x="408.3" y="215" width="8.633" height="8.633" fill="#e7f0fa"><title>HF_severity 2.261, eGFR 12.5: 0.255</title></rect>
<rect x="416.7" y="215" width="8.633" height="8.633" fill="#eaf2fb"><title>HF_severity 2.435, eGFR 12.5: 0.242</title></rect>
<rect x="425" y="215" width="8.633" height="8.633" fill="#edf4fc"><title>HF_severity 2.609, eGFR 12.5: 0.229</title></rect>
<rect x="433.3" y="215" width="8.633" height="8.633" fill="#f0f6fd"><title>HF_severity 2.783, eGFR 12.5: 0.215</title></rect>
<rect x="441.7" y="215" width="8.633" height="8.633" fill="#f3f9fe"><title>HF_severity 2.957, eGFR 12.5: 0.202</title></rect>
<rect x="450" y="215" width="8.633" height="8.633" fill="#f3f8fe"><title>HF_severity 3.13, eGFR 12.5: 0.203</title></rect>
<rect x="458.3" y="215" width="8.633" height="8.633" fill="#f2f8fd"><title>HF_severity 3.304, eGFR 12.5: 0.208</title></rect>
<rect x="466.7" y="215" width="8.633" height="8.633" fill="#f1f7fd"><title>HF_severity 3.478, eGFR 12.5: 0.213</title></rect>
<rect x="475" y="215" width="8.633" height="8.633" fill="#eff6fd"><title>HF_severity 3.652, eGFR 12.5: 0.218</title></rect>
<rect x="483.3" y="215" width="8.633" height="8.633" fill="#eef5fc"><title>HF_severity 3.826, eGFR 12.5: 0.223</title></rect>
<rect x="491.7" y="215" width="8.633" height="8.633" fill="#edf4fc"><title>HF_severity 4, eGFR 12.5: 0.228</title></rect>
<rect x="300" y="223.3" width="8.633" height="8.633" fill="#d3e4f3"><title>HF_severity 0, eGFR 7.5: 0.334</title></rect>
<rect x="308.3" y="223.3" width="8.633" height="8.633" fill="#d5e5f4"><title>HF_severity 0.174, eGFR 7.5: 0.327</title></rect>
<rect x="316.7" y="223.3" width="8.633" height="8.633" fill="#d7e6f4"><title>HF_severity 0.348, eGFR 7.5: 0.32</title></rect>
<rect x="325" y="223.3" width="8.633" height="8.633" fill="#d8e7f5"><title>HF_severity 0.522, eGFR 7.5: 0.313</title></rect>
<rect x="333.3" y="223.3" width="8.633" height="8.633" fill="#dae8f6"><title>HF_severity 0.696, eGFR 7.5: 0.306</title></rect>
<rect x="341.7" y="223.3" width="8.633" height="8.633" fill="#dce9f6"><title>HF_severity 0.87, eGFR 7.5: 0.299</title></rect>
<rect x="350" y="223.3" width="8.633" height="8.633" fill="#deeaf7"><title>HF_severity 1.043, eGFR 7.5: 0.292</title></rect>
<rect x="358.3" y="223.3" width="8.633" height="8.633" fill="#dfebf7"><title>HF_severity 1.217, eGFR 7.5: 0.285</title></rect>
<rect x="366.7" y="223.3" width="8.633" height="8.633" fill="#e1edf8"><title>HF_severity 1.391, eGFR 7.5: 0.278</title></rect>
<rect x="375" y="223.3" width="8.633" height="8.633" fill="#e3eef8"><title>HF_severity 1.565, eGFR 7.5: 0.271</title></rect>
<rect x="383.3" y="223.3" width="8.633" height="8.633" fill="#e4eff9"><title>HF_severity 1.739, eGFR 7.5: 0.264</title></rect>
<rect x="391.7" y="223.3" width="8.633" height="8.633" fill="#e6f0f9"><title>HF_severity 1.913, eGFR 7.5: 0.257</title></rect>
<rect x="400" y="223.3" width="8.633" height="8.633" fill="#e8f1fa"><title>HF_severity 2.087, eGFR 7.5: 0.248</title></rect>
<rect x="408.3" y="223.3" width="8.633" height="8.633" fill="#ebf3fb"><title>HF_severity 2.261, eGFR 7.5: 0.235</title></rect>
<rect x="416.7" y="223.3" width="8.633" height="8.633" fill="#eef5fc"><title>HF_severity 2.435, eGFR 7.5: 0.223</title></rect>
<rect x="425" y="223.3" width="8.633" height="8.633" fill="#f1f7fd"><title>HF_severity 2.609, eGFR 7.5: 0.211</title></rect>
<rect x="433.3" y="223.3" width="8.633" height="8.633" fill="#f4f9fe"><title>HF_severity 2.783, eGFR 7.5: 0.199</title></rect>
<rect x="441.7" y="223.3" width="8.633" height="8.633" fill="#f7fbff"><title>HF_severity 2.957, eGFR 7.5: 0.187</title></rect>
<rect x="450" y="223.3" width="8.633" height="8.633" fill="#f7fbff"><title>HF_severity 3.13, eGFR 7.5: 0.187</title></rect>
<rect x="458.3" y="223.3" width="8.633" height="8.633" fill="#f6faff"><title>HF_severity 3.304, eGFR 7.5: 0.192</title></rect>
<rect x="466.7" y="223.3" width="8.633" height="8.633" fill="#f5f9fe"><title>HF_severity 3.478, eGFR 7.5: 0.197</title></rect>
<rect x="475" y="223.3" width="8.633" height="8.633" fill="#f3f9fe"><title>HF_severity 3.652, eGFR 7.5: 0.201</title></rect>
<rect x="483.3" y="223.3" width="8.633" height="8.633" fill="#f2f8fd"><title>HF_severity 3.826, eGFR 7.5: 0.206</title></rect>
<rect x="491.7" y="223.3" width="8.633" height="8.633" fill="#f1f7fd"><title>HF_severity 4, eGFR 7.5: 0.211</title></rect>
<rect x="300" y="231.7" width="8.633" height="8.633" fill="#d3e4f3"><title>HF_severity 0, eGFR 2.5: 0.334</title></rect>
<rect x="308.3" y="231.7" width="8.633" height="8.633" fill="#d5e5f4"><title>HF_severity 0.174, eGFR 2.5: 0.327</title></rect>
<rect x="316.7" y="231.7" width="8.633" height="8.633" fill="#d7e6f4"><title>HF_severity 0.348, eGFR 2.5: 0.32</title></rect>
<rect x="325" y="231.7" width="8.633" height="8.633" fill="#d8e7f5"><title>HF_severity 0.522, eGFR 2.5: 0.313</title></rect>
<rect x="333.3" y="231.7" width="8.633" height="8.633" fill="#dae8f6"><title>HF_severity 0.696, eGFR 2.5: 0.306</title></rect>
<rect x="341.7" y="231.7" width="8.633" height="8.633" fill="#dce9f6"><title>HF_severity 0.87, eGFR 2.5: 0.299</title></rect>
<rect x="350" y="231.7" width="8.633" height="8.633" fill="#deeaf7"><title>HF_severity 1.043, eGFR 2.5: 0.292</title></rect>
<rect x="358.3" y="231.7" width="8.633" height="8.633" fill="#dfebf7"><title>HF_severity 1.217, eGFR 2.5: 0.285</title></rect>
<rect x="366.7" y="231.7" width="8.633" height="8.633" fill="#e1edf8"><title>HF_severity 1.391, eGFR 2.5: 0.278</title></rect>
<rect x="375" y="231.7" width="8.633" height="8.633" fill="#e3eef8"><title>HF_severity 1.565, eGFR 2.5: 0.271</title></rect>
<rect x="383.3" y="231.7" width="8.633" height="8.633" fill="#e4eff9"><title>HF_severity 1.739, eGFR 2.5: 0.264</title></rect>
<rect x="391.7" y="231.7" width="8.633" height="8.633" fill="#e6f0f9"><title>HF_severity 1.913, eGFR 2.5: 0.257</title></rect>
<rect x="400" y="231.7" width="8.633" height="8.633" fill="#e8f1fa"><title>HF_severity 2.087, eGFR 2.5: 0.248</title></rect>
<rect x="408.3" y="231.7" width="8.633" height="8.633" fill="#ebf3fb"><title>HF_severity 2.261, eGFR 2.5: 0.235</title></rect>
<rect x="416.7" y="231.7" width="8.633" height="8.633" fill="#eef5fc"><title>HF_severity 2.435, eGFR 2.5: 0.223</title></rect>
<rect x="425" y="231.7" width="8.633" height="8.633" fill="#f1f7fd"><title>HF_severity 2.609, eGFR 2.5: 0.211</title></rect>
<rect x="433.3" y="231.7" width="8.633" height="8.633" fill="#f4f9fe"><title>HF_severity 2.783, eGFR 2.5: 0.199</title></rect>
<rect x="441.7" y="231.7" width="8.633" height="8.633" fill="#f7fbff"><title>HF_severity 2.957, eGFR 2.5: 0.187</title></rect>
<rect x="450" y="231.7" width="8.633" height="8.633" fill="#f7fbff"><title>HF_severity 3.13, eGFR 2.5: 0.187</title></rect>
<rect x="458.3" y="231.7" width="8.633" height="8.633" fill="#f6faff"><title>HF_severity 3.304, eGFR 2.5: 0.192</title></rect>
<rect x="466.7" y="231.7" width="8.633" height="8.633" fill="#f5f9fe"><title>HF_severity 3.478, eGFR 2.5: 0.197</title></rect>
<rect x="475" y="231.7" width="8.633" height="8.633" fill="#f3f9fe"><title>HF_severity 3.652, eGFR 2.5: 0.201</title></rect>
<rect x="483.3" y="231.7" width="8.633" height="8.633" fill="#f2f8fd"><title>HF_severity 3.826, eGFR 2.5: 0.206</title></rect>
<rect x="491.7" y="231.7" width="8.633" height="8.633" fill="#f1f7fd"><title>HF_severity 4, eGFR 2.5: 0.211</title></rect>
<line x1="296" x2="300" y1="90" y2="90" stroke="currentColor"/>
<text x="294" y="94" text-anchor="end" font-size="11">90</text>
<line x1="296" x2="300" y1="140" y2="140" stroke="currentColor"/>
<text x="294" y="144" text-anchor="end" font-size="11">60</text>
<line x1="296" x2="300" y1="190" y2="190" stroke="currentColor"/>
<text x="294" y="194" text-anchor="end" font-size="11">30</text>
<line x1="296" x2="300" y1="215" y2="215" stroke="currentColor"/>
<text x="294" y="219" text-anchor="end" font-size="11">15</text>
<line x1="304.2" x2="304.2" y1="240" y2="244" stroke="currentColor"/>
<text x="304.2" y="257" text-anchor="middle" font-size="11">0</text>
<line x1="352.1" x2="352.1" y1="240" y2="244" stroke="currentColor"/>
<text x="352.1" y="257" text-anchor="middle" font-size="11">I</text>
<line x1="400" x2="400" y1="240" y2="244" stroke="currentColor"/>
<text x="400" y="257" text-anchor="middle" font-size="11">II</text>
<line x1="447.9" x2="447.9" y1="240" y2="244" stroke="currentColor"/>
<text x="447.9" y="257" text-anchor="middle" font-size="11">III</text>
<line x1="495.8" x2="495.8" y1="240" y2="244" stroke="currentColor"/>
<text x="495.8" y="257" text-anchor="middle" font-size="11">IV</text>
<text x="400" y="273" text-anchor="middle" font-size="12">HF_severity</text>
<rect x="300" y="282" width="10.3" height="8" fill="#f7fbff"/>
<rect x="310" y="282" width="10.3" height="8" fill="#edf4fc"/>
<rect x="320" y="282" width="10.3" height="8" fill="#e2eef8"/>
<rect x="330" y="282" width="10.3" height="8" fill="#d8e7f5"/>
<rect x="340" y="282" width="10.3" height="8" fill="#cee0f2"/>
<rect x="350" y="282" width="10.3" height="8" fill="#c1d9ee"/>
<rect x="360" y="282" width="10.3" height="8" fill="#aecfe8"/>
<rect x="370" y="282" width="10.3" height="8" fill="#9bc6e3"/>
<rect x="380" y="282" width="10.3" height="8" fill="#88bcde"/>
<rect x="390" y="282" width="10.3" height="8" fill="#75b3d9"/>
<rect x="400" y="282" width="10.3" height="8" fill="#63a8d3"/>
<rect x="410" y="282" width="10.3" height="8" fill="#549bcc"/>
<rect x="420" y="282" width="10.3" height="8" fill="#448ec5"/>
<rect x="430" y="282" width="10.3" height="8" fill="#3481be"/>
<rect x="440" y="282" width="10.3" height="8" fill="#2574b7"/>
<rect x="450" y="282" width="10.3" height="8" fill="#1d67a9"/>
<rect x="460" y="282" width="10.3" height="8" fill="#18599a"/>
<rect x="470" y="282" width="10.3" height="8" fill="#134b8a"/>
<rect x="480" y="282" width="10.3" height="8" fill="#0d3e7b"/>
<rect x="490" y="282" width="10.3" height="8" fill="#08306b"/>
<text x="300" y="304" font-size="11">0.19</text>
<text x="500" y="304" font-size="11" text-anchor="end">1</text>
<text x="550" y="18" font-weight="600">fu of an albumin-bound drug</text>
<text transform="translate(522 140) rotate(-90)" text-anchor="middle" font-size="12">eGFR (ml/min/1.73 m2)</text>
<rect x="550" y="40" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 0, eGFR 117.5: 1</title></rect>
<rect x="558.3" y="40" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 0.13, eGFR 117.5: 1.02</title></rect>
<rect x="566.7" y="40" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.261, eGFR 117.5: 1.03</title></rect>
<rect x="575" y="40" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0.391, eGFR 117.5: 1.05</title></rect>
<rect x="583.3" y="40" width="8.633" height="8.633" fill="#f0f7fd"><title>CP_severity 0.522, eGFR 117.5: 1.07</title></rect>
<rect x="591.7" y="40" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0.652, eGFR 117.5: 1.09</title></rect>
<rect x="600" y="40" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 0.783, eGFR 117.5: 1.11</title></rect>
<rect x="608.3" y="40" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.913, eGFR 117.5: 1.13</title></rect>
<rect x="616.7" y="40" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 1.043, eGFR 117.5: 1.16</title></rect>
<rect x="625" y="40" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 1.174, eGFR 117.5: 1.19</title></rect>
<rect x="633.3" y="40" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 1.304, eGFR 117.5: 1.22</title></rect>
<rect x="641.7" y="40" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 1.435, eGFR 117.5: 1.26</title></rect>
<rect x="650" y="40" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 1.565, eGFR 117.5: 1.3</title></rect>
<rect x="658.3" y="40" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.696, eGFR 117.5: 1.33</title></rect>
<rect x="666.7" y="40" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.826, eGFR 117.5: 1.38</title></rect>
<rect x="675" y="40" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 1.957, eGFR 117.5: 1.42</title></rect>
<rect x="683.3" y="40" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 2.087, eGFR 117.5: 1.47</title></rect>
<rect x="691.7" y="40" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 2.217, eGFR 117.5: 1.52</title></rect>
<rect x="700" y="40" width="8.633" height="8.633" fill="#bfd8ed"><title>CP_severity 2.348, eGFR 117.5: 1.58</title></rect>
<rect x="708.3" y="40" width="8.633" height="8.633" fill="#b5d3ea"><title>CP_severity 2.478, eGFR 117.5: 1.64</title></rect>
<rect x="716.7" y="40" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.609, eGFR 117.5: 1.71</title></rect>
<rect x="725" y="40" width="8.633" height="8.633" fill="#9ec7e4"><title>CP_severity 2.739, eGFR 117.5: 1.78</title></rect>
<rect x="733.3" y="40" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.87, eGFR 117.5: 1.86</title></rect>
<rect x="741.7" y="40" width="8.633" height="8.633" fill="#83badc"><title>CP_severity 3, eGFR 117.5: 1.94</title></rect>
<rect x="550" y="48.33" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 0, eGFR 112.5: 1</title></rect>
<rect x="558.3" y="48.33" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 0.13, eGFR 112.5: 1.02</title></rect>
<rect x="566.7" y="48.33" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.261, eGFR 112.5: 1.03</title></rect>
<rect x="575" y="48.33" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0.391, eGFR 112.5: 1.05</title></rect>
<rect x="583.3" y="48.33" width="8.633" height="8.633" fill="#f0f7fd"><title>CP_severity 0.522, eGFR 112.5: 1.07</title></rect>
<rect x="591.7" y="48.33" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0.652, eGFR 112.5: 1.09</title></rect>
<rect x="600" y="48.33" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 0.783, eGFR 112.5: 1.11</title></rect>
<rect x="608.3" y="48.33" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.913, eGFR 112.5: 1.13</title></rect>
<rect x="616.7" y="48.33" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 1.043, eGFR 112.5: 1.16</title></rect>
<rect x="625" y="48.33" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 1.174, eGFR 112.5: 1.19</title></rect>
<rect x="633.3" y="48.33" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 1.304, eGFR 112.5: 1.22</title></rect>
<rect x="641.7" y="48.33" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 1.435, eGFR 112.5: 1.26</title></rect>
<rect x="650" y="48.33" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 1.565, eGFR 112.5: 1.3</title></rect>
<rect x="658.3" y="48.33" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.696, eGFR 112.5: 1.33</title></rect>
<rect x="666.7" y="48.33" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.826, eGFR 112.5: 1.38</title></rect>
<rect x="675" y="48.33" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 1.957, eGFR 112.5: 1.42</title></rect>
<rect x="683.3" y="48.33" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 2.087, eGFR 112.5: 1.47</title></rect>
<rect x="691.7" y="48.33" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 2.217, eGFR 112.5: 1.52</title></rect>
<rect x="700" y="48.33" width="8.633" height="8.633" fill="#bfd8ed"><title>CP_severity 2.348, eGFR 112.5: 1.58</title></rect>
<rect x="708.3" y="48.33" width="8.633" height="8.633" fill="#b5d3ea"><title>CP_severity 2.478, eGFR 112.5: 1.64</title></rect>
<rect x="716.7" y="48.33" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.609, eGFR 112.5: 1.71</title></rect>
<rect x="725" y="48.33" width="8.633" height="8.633" fill="#9ec7e4"><title>CP_severity 2.739, eGFR 112.5: 1.78</title></rect>
<rect x="733.3" y="48.33" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.87, eGFR 112.5: 1.86</title></rect>
<rect x="741.7" y="48.33" width="8.633" height="8.633" fill="#83badc"><title>CP_severity 3, eGFR 112.5: 1.94</title></rect>
<rect x="550" y="56.67" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 0, eGFR 107.5: 1</title></rect>
<rect x="558.3" y="56.67" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 0.13, eGFR 107.5: 1.02</title></rect>
<rect x="566.7" y="56.67" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.261, eGFR 107.5: 1.03</title></rect>
<rect x="575" y="56.67" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0.391, eGFR 107.5: 1.05</title></rect>
<rect x="583.3" y="56.67" width="8.633" height="8.633" fill="#f0f7fd"><title>CP_severity 0.522, eGFR 107.5: 1.07</title></rect>
<rect x="591.7" y="56.67" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0.652, eGFR 107.5: 1.09</title></rect>
<rect x="600" y="56.67" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 0.783, eGFR 107.5: 1.11</title></rect>
<rect x="608.3" y="56.67" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.913, eGFR 107.5: 1.13</title></rect>
<rect x="616.7" y="56.67" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 1.043, eGFR 107.5: 1.16</title></rect>
<rect x="625" y="56.67" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 1.174, eGFR 107.5: 1.19</title></rect>
<rect x="633.3" y="56.67" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 1.304, eGFR 107.5: 1.22</title></rect>
<rect x="641.7" y="56.67" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 1.435, eGFR 107.5: 1.26</title></rect>
<rect x="650" y="56.67" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 1.565, eGFR 107.5: 1.3</title></rect>
<rect x="658.3" y="56.67" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.696, eGFR 107.5: 1.33</title></rect>
<rect x="666.7" y="56.67" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.826, eGFR 107.5: 1.38</title></rect>
<rect x="675" y="56.67" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 1.957, eGFR 107.5: 1.42</title></rect>
<rect x="683.3" y="56.67" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 2.087, eGFR 107.5: 1.47</title></rect>
<rect x="691.7" y="56.67" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 2.217, eGFR 107.5: 1.52</title></rect>
<rect x="700" y="56.67" width="8.633" height="8.633" fill="#bfd8ed"><title>CP_severity 2.348, eGFR 107.5: 1.58</title></rect>
<rect x="708.3" y="56.67" width="8.633" height="8.633" fill="#b5d3ea"><title>CP_severity 2.478, eGFR 107.5: 1.64</title></rect>
<rect x="716.7" y="56.67" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.609, eGFR 107.5: 1.71</title></rect>
<rect x="725" y="56.67" width="8.633" height="8.633" fill="#9ec7e4"><title>CP_severity 2.739, eGFR 107.5: 1.78</title></rect>
<rect x="733.3" y="56.67" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.87, eGFR 107.5: 1.86</title></rect>
<rect x="741.7" y="56.67" width="8.633" height="8.633" fill="#83badc"><title>CP_severity 3, eGFR 107.5: 1.94</title></rect>
<rect x="550" y="65" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 0, eGFR 102.5: 1</title></rect>
<rect x="558.3" y="65" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 0.13, eGFR 102.5: 1.02</title></rect>
<rect x="566.7" y="65" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.261, eGFR 102.5: 1.03</title></rect>
<rect x="575" y="65" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0.391, eGFR 102.5: 1.05</title></rect>
<rect x="583.3" y="65" width="8.633" height="8.633" fill="#f0f7fd"><title>CP_severity 0.522, eGFR 102.5: 1.07</title></rect>
<rect x="591.7" y="65" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0.652, eGFR 102.5: 1.09</title></rect>
<rect x="600" y="65" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 0.783, eGFR 102.5: 1.11</title></rect>
<rect x="608.3" y="65" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.913, eGFR 102.5: 1.13</title></rect>
<rect x="616.7" y="65" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 1.043, eGFR 102.5: 1.16</title></rect>
<rect x="625" y="65" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 1.174, eGFR 102.5: 1.19</title></rect>
<rect x="633.3" y="65" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 1.304, eGFR 102.5: 1.22</title></rect>
<rect x="641.7" y="65" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 1.435, eGFR 102.5: 1.26</title></rect>
<rect x="650" y="65" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 1.565, eGFR 102.5: 1.3</title></rect>
<rect x="658.3" y="65" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.696, eGFR 102.5: 1.33</title></rect>
<rect x="666.7" y="65" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.826, eGFR 102.5: 1.38</title></rect>
<rect x="675" y="65" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 1.957, eGFR 102.5: 1.42</title></rect>
<rect x="683.3" y="65" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 2.087, eGFR 102.5: 1.47</title></rect>
<rect x="691.7" y="65" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 2.217, eGFR 102.5: 1.52</title></rect>
<rect x="700" y="65" width="8.633" height="8.633" fill="#bfd8ed"><title>CP_severity 2.348, eGFR 102.5: 1.58</title></rect>
<rect x="708.3" y="65" width="8.633" height="8.633" fill="#b5d3ea"><title>CP_severity 2.478, eGFR 102.5: 1.64</title></rect>
<rect x="716.7" y="65" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.609, eGFR 102.5: 1.71</title></rect>
<rect x="725" y="65" width="8.633" height="8.633" fill="#9ec7e4"><title>CP_severity 2.739, eGFR 102.5: 1.78</title></rect>
<rect x="733.3" y="65" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.87, eGFR 102.5: 1.86</title></rect>
<rect x="741.7" y="65" width="8.633" height="8.633" fill="#83badc"><title>CP_severity 3, eGFR 102.5: 1.94</title></rect>
<rect x="550" y="73.33" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 0, eGFR 97.5: 1</title></rect>
<rect x="558.3" y="73.33" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 0.13, eGFR 97.5: 1.02</title></rect>
<rect x="566.7" y="73.33" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.261, eGFR 97.5: 1.03</title></rect>
<rect x="575" y="73.33" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0.391, eGFR 97.5: 1.05</title></rect>
<rect x="583.3" y="73.33" width="8.633" height="8.633" fill="#f0f7fd"><title>CP_severity 0.522, eGFR 97.5: 1.07</title></rect>
<rect x="591.7" y="73.33" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0.652, eGFR 97.5: 1.09</title></rect>
<rect x="600" y="73.33" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 0.783, eGFR 97.5: 1.11</title></rect>
<rect x="608.3" y="73.33" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.913, eGFR 97.5: 1.13</title></rect>
<rect x="616.7" y="73.33" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 1.043, eGFR 97.5: 1.16</title></rect>
<rect x="625" y="73.33" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 1.174, eGFR 97.5: 1.19</title></rect>
<rect x="633.3" y="73.33" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 1.304, eGFR 97.5: 1.22</title></rect>
<rect x="641.7" y="73.33" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 1.435, eGFR 97.5: 1.26</title></rect>
<rect x="650" y="73.33" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 1.565, eGFR 97.5: 1.3</title></rect>
<rect x="658.3" y="73.33" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.696, eGFR 97.5: 1.33</title></rect>
<rect x="666.7" y="73.33" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.826, eGFR 97.5: 1.38</title></rect>
<rect x="675" y="73.33" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 1.957, eGFR 97.5: 1.42</title></rect>
<rect x="683.3" y="73.33" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 2.087, eGFR 97.5: 1.47</title></rect>
<rect x="691.7" y="73.33" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 2.217, eGFR 97.5: 1.52</title></rect>
<rect x="700" y="73.33" width="8.633" height="8.633" fill="#bfd8ed"><title>CP_severity 2.348, eGFR 97.5: 1.58</title></rect>
<rect x="708.3" y="73.33" width="8.633" height="8.633" fill="#b5d3ea"><title>CP_severity 2.478, eGFR 97.5: 1.64</title></rect>
<rect x="716.7" y="73.33" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.609, eGFR 97.5: 1.71</title></rect>
<rect x="725" y="73.33" width="8.633" height="8.633" fill="#9ec7e4"><title>CP_severity 2.739, eGFR 97.5: 1.78</title></rect>
<rect x="733.3" y="73.33" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.87, eGFR 97.5: 1.86</title></rect>
<rect x="741.7" y="73.33" width="8.633" height="8.633" fill="#83badc"><title>CP_severity 3, eGFR 97.5: 1.94</title></rect>
<rect x="550" y="81.67" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 0, eGFR 92.5: 1</title></rect>
<rect x="558.3" y="81.67" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 0.13, eGFR 92.5: 1.02</title></rect>
<rect x="566.7" y="81.67" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.261, eGFR 92.5: 1.03</title></rect>
<rect x="575" y="81.67" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0.391, eGFR 92.5: 1.05</title></rect>
<rect x="583.3" y="81.67" width="8.633" height="8.633" fill="#f0f7fd"><title>CP_severity 0.522, eGFR 92.5: 1.07</title></rect>
<rect x="591.7" y="81.67" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0.652, eGFR 92.5: 1.09</title></rect>
<rect x="600" y="81.67" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 0.783, eGFR 92.5: 1.11</title></rect>
<rect x="608.3" y="81.67" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.913, eGFR 92.5: 1.13</title></rect>
<rect x="616.7" y="81.67" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 1.043, eGFR 92.5: 1.16</title></rect>
<rect x="625" y="81.67" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 1.174, eGFR 92.5: 1.19</title></rect>
<rect x="633.3" y="81.67" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 1.304, eGFR 92.5: 1.22</title></rect>
<rect x="641.7" y="81.67" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 1.435, eGFR 92.5: 1.26</title></rect>
<rect x="650" y="81.67" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 1.565, eGFR 92.5: 1.3</title></rect>
<rect x="658.3" y="81.67" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.696, eGFR 92.5: 1.33</title></rect>
<rect x="666.7" y="81.67" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.826, eGFR 92.5: 1.38</title></rect>
<rect x="675" y="81.67" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 1.957, eGFR 92.5: 1.42</title></rect>
<rect x="683.3" y="81.67" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 2.087, eGFR 92.5: 1.47</title></rect>
<rect x="691.7" y="81.67" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 2.217, eGFR 92.5: 1.52</title></rect>
<rect x="700" y="81.67" width="8.633" height="8.633" fill="#bfd8ed"><title>CP_severity 2.348, eGFR 92.5: 1.58</title></rect>
<rect x="708.3" y="81.67" width="8.633" height="8.633" fill="#b5d3ea"><title>CP_severity 2.478, eGFR 92.5: 1.64</title></rect>
<rect x="716.7" y="81.67" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.609, eGFR 92.5: 1.71</title></rect>
<rect x="725" y="81.67" width="8.633" height="8.633" fill="#9ec7e4"><title>CP_severity 2.739, eGFR 92.5: 1.78</title></rect>
<rect x="733.3" y="81.67" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.87, eGFR 92.5: 1.86</title></rect>
<rect x="741.7" y="81.67" width="8.633" height="8.633" fill="#83badc"><title>CP_severity 3, eGFR 92.5: 1.94</title></rect>
<rect x="550" y="90" width="8.633" height="8.633" fill="#f7fbff"><title>CP_severity 0, eGFR 87.5: 1.01</title></rect>
<rect x="558.3" y="90" width="8.633" height="8.633" fill="#f5fafe"><title>CP_severity 0.13, eGFR 87.5: 1.02</title></rect>
<rect x="566.7" y="90" width="8.633" height="8.633" fill="#f3f9fe"><title>CP_severity 0.261, eGFR 87.5: 1.04</title></rect>
<rect x="575" y="90" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0.391, eGFR 87.5: 1.06</title></rect>
<rect x="583.3" y="90" width="8.633" height="8.633" fill="#f0f6fd"><title>CP_severity 0.522, eGFR 87.5: 1.08</title></rect>
<rect x="591.7" y="90" width="8.633" height="8.633" fill="#eef5fc"><title>CP_severity 0.652, eGFR 87.5: 1.1</title></rect>
<rect x="600" y="90" width="8.633" height="8.633" fill="#ecf4fc"><title>CP_severity 0.783, eGFR 87.5: 1.12</title></rect>
<rect x="608.3" y="90" width="8.633" height="8.633" fill="#eaf3fb"><title>CP_severity 0.913, eGFR 87.5: 1.14</title></rect>
<rect x="616.7" y="90" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 1.043, eGFR 87.5: 1.16</title></rect>
<rect x="625" y="90" width="8.633" height="8.633" fill="#e5eff9"><title>CP_severity 1.174, eGFR 87.5: 1.2</title></rect>
<rect x="633.3" y="90" width="8.633" height="8.633" fill="#e2edf8"><title>CP_severity 1.304, eGFR 87.5: 1.23</title></rect>
<rect x="641.7" y="90" width="8.633" height="8.633" fill="#dfebf7"><title>CP_severity 1.435, eGFR 87.5: 1.26</title></rect>
<rect x="650" y="90" width="8.633" height="8.633" fill="#dce9f6"><title>CP_severity 1.565, eGFR 87.5: 1.3</title></rect>
<rect x="658.3" y="90" width="8.633" height="8.633" fill="#d8e7f5"><title>CP_severity 1.696, eGFR 87.5: 1.34</title></rect>
<rect x="666.7" y="90" width="8.633" height="8.633" fill="#d4e4f4"><title>CP_severity 1.826, eGFR 87.5: 1.38</title></rect>
<rect x="675" y="90" width="8.633" height="8.633" fill="#d0e2f2"><title>CP_severity 1.957, eGFR 87.5: 1.43</title></rect>
<rect x="683.3" y="90" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 2.087, eGFR 87.5: 1.48</title></rect>
<rect x="691.7" y="90" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 2.217, eGFR 87.5: 1.53</title></rect>
<rect x="700" y="90" width="8.633" height="8.633" fill="#bed7ed"><title>CP_severity 2.348, eGFR 87.5: 1.59</title></rect>
<rect x="708.3" y="90" width="8.633" height="8.633" fill="#b4d2ea"><title>CP_severity 2.478, eGFR 87.5: 1.65</title></rect>
<rect x="716.7" y="90" width="8.633" height="8.633" fill="#a9cce7"><title>CP_severity 2.609, eGFR 87.5: 1.72</title></rect>
<rect x="725" y="90" width="8.633" height="8.633" fill="#9cc6e4"><title>CP_severity 2.739, eGFR 87.5: 1.79</title></rect>
<rect x="733.3" y="90" width="8.633" height="8.633" fill="#8fc0e0"><title>CP_severity 2.87, eGFR 87.5: 1.87</title></rect>
<rect x="741.7" y="90" width="8.633" height="8.633" fill="#81b9dc"><title>CP_severity 3, eGFR 87.5: 1.95</title></rect>
<rect x="550" y="98.33" width="8.633" height="8.633" fill="#f5faff"><title>CP_severity 0, eGFR 82.5: 1.02</title></rect>
<rect x="558.3" y="98.33" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0.13, eGFR 82.5: 1.03</title></rect>
<rect x="566.7" y="98.33" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0.261, eGFR 82.5: 1.05</title></rect>
<rect x="575" y="98.33" width="8.633" height="8.633" fill="#f1f7fd"><title>CP_severity 0.391, eGFR 82.5: 1.07</title></rect>
<rect x="583.3" y="98.33" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0.522, eGFR 82.5: 1.09</title></rect>
<rect x="591.7" y="98.33" width="8.633" height="8.633" fill="#edf5fc"><title>CP_severity 0.652, eGFR 82.5: 1.11</title></rect>
<rect x="600" y="98.33" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.783, eGFR 82.5: 1.13</title></rect>
<rect x="608.3" y="98.33" width="8.633" height="8.633" fill="#e9f2fb"><title>CP_severity 0.913, eGFR 82.5: 1.15</title></rect>
<rect x="616.7" y="98.33" width="8.633" height="8.633" fill="#e7f1fa"><title>CP_severity 1.043, eGFR 82.5: 1.18</title></rect>
<rect x="625" y="98.33" width="8.633" height="8.633" fill="#e4eff9"><title>CP_severity 1.174, eGFR 82.5: 1.21</title></rect>
<rect x="633.3" y="98.33" width="8.633" height="8.633" fill="#e1edf8"><title>CP_severity 1.304, eGFR 82.5: 1.24</title></rect>
<rect x="641.7" y="98.33" width="8.633" height="8.633" fill="#deebf7"><title>CP_severity 1.435, eGFR 82.5: 1.28</title></rect>
<rect x="650" y="98.33" width="8.633" height="8.633" fill="#dae8f6"><title>CP_severity 1.565, eGFR 82.5: 1.32</title></rect>
<rect x="658.3" y="98.33" width="8.633" height="8.633" fill="#d7e6f4"><title>CP_severity 1.696, eGFR 82.5: 1.36</title></rect>
<rect x="666.7" y="98.33" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 1.826, eGFR 82.5: 1.4</title></rect>
<rect x="675" y="98.33" width="8.633" height="8.633" fill="#cfe1f2"><title>CP_severity 1.957, eGFR 82.5: 1.44</title></rect>
<rect x="683.3" y="98.33" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 2.087, eGFR 82.5: 1.49</title></rect>
<rect x="691.7" y="98.33" width="8.633" height="8.633" fill="#c5daef"><title>CP_severity 2.217, eGFR 82.5: 1.55</title></rect>
<rect x="700" y="98.33" width="8.633" height="8.633" fill="#bbd6ec"><title>CP_severity 2.348, eGFR 82.5: 1.61</title></rect>
<rect x="708.3" y="98.33" width="8.633" height="8.633" fill="#b1d0e9"><title>CP_severity 2.478, eGFR 82.5: 1.67</title></rect>
<rect x="716.7" y="98.33" width="8.633" height="8.633" fill="#a5cbe6"><title>CP_severity 2.609, eGFR 82.5: 1.74</title></rect>
<rect x="725" y="98.33" width="8.633" height="8.633" fill="#99c5e3"><title>CP_severity 2.739, eGFR 82.5: 1.81</title></rect>
<rect x="733.3" y="98.33" width="8.633" height="8.633" fill="#8cbedf"><title>CP_severity 2.87, eGFR 82.5: 1.89</title></rect>
<rect x="741.7" y="98.33" width="8.633" height="8.633" fill="#7db7db"><title>CP_severity 3, eGFR 82.5: 1.97</title></rect>
<rect x="550" y="106.7" width="8.633" height="8.633" fill="#f4f9fe"><title>CP_severity 0, eGFR 77.5: 1.03</title></rect>
<rect x="558.3" y="106.7" width="8.633" height="8.633" fill="#f3f8fe"><title>CP_severity 0.13, eGFR 77.5: 1.05</title></rect>
<rect x="566.7" y="106.7" width="8.633" height="8.633" fill="#f1f7fd"><title>CP_severity 0.261, eGFR 77.5: 1.06</title></rect>
<rect x="575" y="106.7" width="8.633" height="8.633" fill="#f0f6fd"><title>CP_severity 0.391, eGFR 77.5: 1.08</title></rect>
<rect x="583.3" y="106.7" width="8.633" height="8.633" fill="#eef5fc"><title>CP_severity 0.522, eGFR 77.5: 1.1</title></rect>
<rect x="591.7" y="106.7" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.652, eGFR 77.5: 1.12</title></rect>
<rect x="600" y="106.7" width="8.633" height="8.633" fill="#eaf3fb"><title>CP_severity 0.783, eGFR 77.5: 1.14</title></rect>
<rect x="608.3" y="106.7" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 0.913, eGFR 77.5: 1.16</title></rect>
<rect x="616.7" y="106.7" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 1.043, eGFR 77.5: 1.19</title></rect>
<rect x="625" y="106.7" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 1.174, eGFR 77.5: 1.22</title></rect>
<rect x="633.3" y="106.7" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 1.304, eGFR 77.5: 1.26</title></rect>
<rect x="641.7" y="106.7" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 1.435, eGFR 77.5: 1.29</title></rect>
<rect x="650" y="106.7" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.565, eGFR 77.5: 1.33</title></rect>
<rect x="658.3" y="106.7" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.696, eGFR 77.5: 1.37</title></rect>
<rect x="666.7" y="106.7" width="8.633" height="8.633" fill="#d1e2f3"><title>CP_severity 1.826, eGFR 77.5: 1.41</title></rect>
<rect x="675" y="106.7" width="8.633" height="8.633" fill="#cde0f1"><title>CP_severity 1.957, eGFR 77.5: 1.46</title></rect>
<rect x="683.3" y="106.7" width="8.633" height="8.633" fill="#c9ddf0"><title>CP_severity 2.087, eGFR 77.5: 1.51</title></rect>
<rect x="691.7" y="106.7" width="8.633" height="8.633" fill="#c2d9ee"><title>CP_severity 2.217, eGFR 77.5: 1.57</title></rect>
<rect x="700" y="106.7" width="8.633" height="8.633" fill="#b8d4eb"><title>CP_severity 2.348, eGFR 77.5: 1.62</title></rect>
<rect x="708.3" y="106.7" width="8.633" height="8.633" fill="#aecfe8"><title>CP_severity 2.478, eGFR 77.5: 1.69</title></rect>
<rect x="716.7" y="106.7" width="8.633" height="8.633" fill="#a2c9e5"><title>CP_severity 2.609, eGFR 77.5: 1.75</title></rect>
<rect x="725" y="106.7" width="8.633" height="8.633" fill="#96c3e2"><title>CP_severity 2.739, eGFR 77.5: 1.83</title></rect>
<rect x="733.3" y="106.7" width="8.633" height="8.633" fill="#88bdde"><title>CP_severity 2.87, eGFR 77.5: 1.91</title></rect>
<rect x="741.7" y="106.7" width="8.633" height="8.633" fill="#7ab5da"><title>CP_severity 3, eGFR 77.5: 1.99</title></rect>
<rect x="550" y="115" width="8.633" height="8.633" fill="#f3f9fe"><title>CP_severity 0, eGFR 72.5: 1.04</title></rect>
<rect x="558.3" y="115" width="8.633" height="8.633" fill="#f2f7fd"><title>CP_severity 0.13, eGFR 72.5: 1.06</title></rect>
<rect x="566.7" y="115" width="8.633" height="8.633" fill="#f0f6fd"><title>CP_severity 0.261, eGFR 72.5: 1.08</title></rect>
<rect x="575" y="115" width="8.633" height="8.633" fill="#eef5fc"><title>CP_severity 0.391, eGFR 72.5: 1.1</title></rect>
<rect x="583.3" y="115" width="8.633" height="8.633" fill="#ecf4fc"><title>CP_severity 0.522, eGFR 72.5: 1.12</title></rect>
<rect x="591.7" y="115" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.652, eGFR 72.5: 1.14</title></rect>
<rect x="600" y="115" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 0.783, eGFR 72.5: 1.16</title></rect>
<rect x="608.3" y="115" width="8.633" height="8.633" fill="#e7f0fa"><title>CP_severity 0.913, eGFR 72.5: 1.18</title></rect>
<rect x="616.7" y="115" width="8.633" height="8.633" fill="#e4eff9"><title>CP_severity 1.043, eGFR 72.5: 1.21</title></rect>
<rect x="625" y="115" width="8.633" height="8.633" fill="#e1edf8"><title>CP_severity 1.174, eGFR 72.5: 1.24</title></rect>
<rect x="633.3" y="115" width="8.633" height="8.633" fill="#deebf7"><title>CP_severity 1.304, eGFR 72.5: 1.27</title></rect>
<rect x="641.7" y="115" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 1.435, eGFR 72.5: 1.31</title></rect>
<rect x="650" y="115" width="8.633" height="8.633" fill="#d7e6f5"><title>CP_severity 1.565, eGFR 72.5: 1.35</title></rect>
<rect x="658.3" y="115" width="8.633" height="8.633" fill="#d4e4f3"><title>CP_severity 1.696, eGFR 72.5: 1.39</title></rect>
<rect x="666.7" y="115" width="8.633" height="8.633" fill="#d0e1f2"><title>CP_severity 1.826, eGFR 72.5: 1.43</title></rect>
<rect x="675" y="115" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 1.957, eGFR 72.5: 1.48</title></rect>
<rect x="683.3" y="115" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 2.087, eGFR 72.5: 1.53</title></rect>
<rect x="691.7" y="115" width="8.633" height="8.633" fill="#bfd7ed"><title>CP_severity 2.217, eGFR 72.5: 1.59</title></rect>
<rect x="700" y="115" width="8.633" height="8.633" fill="#b5d2ea"><title>CP_severity 2.348, eGFR 72.5: 1.64</title></rect>
<rect x="708.3" y="115" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.478, eGFR 72.5: 1.71</title></rect>
<rect x="716.7" y="115" width="8.633" height="8.633" fill="#9fc7e4"><title>CP_severity 2.609, eGFR 72.5: 1.78</title></rect>
<rect x="725" y="115" width="8.633" height="8.633" fill="#92c1e1"><title>CP_severity 2.739, eGFR 72.5: 1.85</title></rect>
<rect x="733.3" y="115" width="8.633" height="8.633" fill="#84bbdd"><title>CP_severity 2.87, eGFR 72.5: 1.93</title></rect>
<rect x="741.7" y="115" width="8.633" height="8.633" fill="#76b3d9"><title>CP_severity 3, eGFR 72.5: 2.02</title></rect>
<rect x="550" y="123.3" width="8.633" height="8.633" fill="#f2f8fd"><title>CP_severity 0, eGFR 67.5: 1.06</title></rect>
<rect x="558.3" y="123.3" width="8.633" height="8.633" fill="#f0f7fd"><title>CP_severity 0.13, eGFR 67.5: 1.08</title></rect>
<rect x="566.7" y="123.3" width="8.633" height="8.633" fill="#eef5fc"><title>CP_severity 0.261, eGFR 67.5: 1.09</title></rect>
<rect x="575" y="123.3" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 0.391, eGFR 67.5: 1.11</title></rect>
<rect x="583.3" y="123.3" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.522, eGFR 67.5: 1.13</title></rect>
<rect x="591.7" y="123.3" width="8.633" height="8.633" fill="#e9f2fa"><title>CP_severity 0.652, eGFR 67.5: 1.15</title></rect>
<rect x="600" y="123.3" width="8.633" height="8.633" fill="#e7f1fa"><title>CP_severity 0.783, eGFR 67.5: 1.18</title></rect>
<rect x="608.3" y="123.3" width="8.633" height="8.633" fill="#e5eff9"><title>CP_severity 0.913, eGFR 67.5: 1.2</title></rect>
<rect x="616.7" y="123.3" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 1.043, eGFR 67.5: 1.22</title></rect>
<rect x="625" y="123.3" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 1.174, eGFR 67.5: 1.26</title></rect>
<rect x="633.3" y="123.3" width="8.633" height="8.633" fill="#ddeaf6"><title>CP_severity 1.304, eGFR 67.5: 1.29</title></rect>
<rect x="641.7" y="123.3" width="8.633" height="8.633" fill="#d9e8f5"><title>CP_severity 1.435, eGFR 67.5: 1.33</title></rect>
<rect x="650" y="123.3" width="8.633" height="8.633" fill="#d6e5f4"><title>CP_severity 1.565, eGFR 67.5: 1.37</title></rect>
<rect x="658.3" y="123.3" width="8.633" height="8.633" fill="#d2e3f3"><title>CP_severity 1.696, eGFR 67.5: 1.41</title></rect>
<rect x="666.7" y="123.3" width="8.633" height="8.633" fill="#cee0f2"><title>CP_severity 1.826, eGFR 67.5: 1.45</title></rect>
<rect x="675" y="123.3" width="8.633" height="8.633" fill="#caddf0"><title>CP_severity 1.957, eGFR 67.5: 1.5</title></rect>
<rect x="683.3" y="123.3" width="8.633" height="8.633" fill="#c4daee"><title>CP_severity 2.087, eGFR 67.5: 1.55</title></rect>
<rect x="691.7" y="123.3" width="8.633" height="8.633" fill="#bbd5ec"><title>CP_severity 2.217, eGFR 67.5: 1.61</title></rect>
<rect x="700" y="123.3" width="8.633" height="8.633" fill="#b1d0e9"><title>CP_severity 2.348, eGFR 67.5: 1.67</title></rect>
<rect x="708.3" y="123.3" width="8.633" height="8.633" fill="#a6cbe6"><title>CP_severity 2.478, eGFR 67.5: 1.73</title></rect>
<rect x="716.7" y="123.3" width="8.633" height="8.633" fill="#9ac5e3"><title>CP_severity 2.609, eGFR 67.5: 1.8</title></rect>
<rect x="725" y="123.3" width="8.633" height="8.633" fill="#8dbfdf"><title>CP_severity 2.739, eGFR 67.5: 1.88</title></rect>
<rect x="733.3" y="123.3" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 2.87, eGFR 67.5: 1.96</title></rect>
<rect x="741.7" y="123.3" width="8.633" height="8.633" fill="#71b1d8"><title>CP_severity 3, eGFR 67.5: 2.05</title></rect>
<rect x="550" y="131.7" width="8.633" height="8.633" fill="#f0f7fd"><title>CP_severity 0, eGFR 62.5: 1.07</title></rect>
<rect x="558.3" y="131.7" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0.13, eGFR 62.5: 1.09</title></rect>
<rect x="566.7" y="131.7" width="8.633" height="8.633" fill="#edf4fc"><title>CP_severity 0.261, eGFR 62.5: 1.11</title></rect>
<rect x="575" y="131.7" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.391, eGFR 62.5: 1.13</title></rect>
<rect x="583.3" y="131.7" width="8.633" height="8.633" fill="#e9f2fb"><title>CP_severity 0.522, eGFR 62.5: 1.15</title></rect>
<rect x="591.7" y="131.7" width="8.633" height="8.633" fill="#e7f1fa"><title>CP_severity 0.652, eGFR 62.5: 1.17</title></rect>
<rect x="600" y="131.7" width="8.633" height="8.633" fill="#e5f0f9"><title>CP_severity 0.783, eGFR 62.5: 1.19</title></rect>
<rect x="608.3" y="131.7" width="8.633" height="8.633" fill="#e3eef9"><title>CP_severity 0.913, eGFR 62.5: 1.22</title></rect>
<rect x="616.7" y="131.7" width="8.633" height="8.633" fill="#e1edf8"><title>CP_severity 1.043, eGFR 62.5: 1.24</title></rect>
<rect x="625" y="131.7" width="8.633" height="8.633" fill="#deebf7"><title>CP_severity 1.174, eGFR 62.5: 1.28</title></rect>
<rect x="633.3" y="131.7" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 1.304, eGFR 62.5: 1.31</title></rect>
<rect x="641.7" y="131.7" width="8.633" height="8.633" fill="#d7e6f5"><title>CP_severity 1.435, eGFR 62.5: 1.35</title></rect>
<rect x="650" y="131.7" width="8.633" height="8.633" fill="#d4e4f3"><title>CP_severity 1.565, eGFR 62.5: 1.39</title></rect>
<rect x="658.3" y="131.7" width="8.633" height="8.633" fill="#d0e1f2"><title>CP_severity 1.696, eGFR 62.5: 1.43</title></rect>
<rect x="666.7" y="131.7" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 1.826, eGFR 62.5: 1.48</title></rect>
<rect x="675" y="131.7" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 1.957, eGFR 62.5: 1.52</title></rect>
<rect x="683.3" y="131.7" width="8.633" height="8.633" fill="#c0d8ed"><title>CP_severity 2.087, eGFR 62.5: 1.58</title></rect>
<rect x="691.7" y="131.7" width="8.633" height="8.633" fill="#b7d3eb"><title>CP_severity 2.217, eGFR 62.5: 1.63</title></rect>
<rect x="700" y="131.7" width="8.633" height="8.633" fill="#accee8"><title>CP_severity 2.348, eGFR 62.5: 1.69</title></rect>
<rect x="708.3" y="131.7" width="8.633" height="8.633" fill="#a2c9e5"><title>CP_severity 2.478, eGFR 62.5: 1.76</title></rect>
<rect x="716.7" y="131.7" width="8.633" height="8.633" fill="#96c3e2"><title>CP_severity 2.609, eGFR 62.5: 1.83</title></rect>
<rect x="725" y="131.7" width="8.633" height="8.633" fill="#89bdde"><title>CP_severity 2.739, eGFR 62.5: 1.91</title></rect>
<rect x="733.3" y="131.7" width="8.633" height="8.633" fill="#7bb6da"><title>CP_severity 2.87, eGFR 62.5: 1.99</title></rect>
<rect x="741.7" y="131.7" width="8.633" height="8.633" fill="#6caed6"><title>CP_severity 3, eGFR 62.5: 2.08</title></rect>
<rect x="550" y="140" width="8.633" height="8.633" fill="#eff6fc"><title>CP_severity 0, eGFR 57.5: 1.09</title></rect>
<rect x="558.3" y="140" width="8.633" height="8.633" fill="#edf5fc"><title>CP_severity 0.13, eGFR 57.5: 1.11</title></rect>
<rect x="566.7" y="140" width="8.633" height="8.633" fill="#ebf3fb"><title>CP_severity 0.261, eGFR 57.5: 1.13</title></rect>
<rect x="575" y="140" width="8.633" height="8.633" fill="#eaf2fb"><title>CP_severity 0.391, eGFR 57.5: 1.15</title></rect>
<rect x="583.3" y="140" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 0.522, eGFR 57.5: 1.17</title></rect>
<rect x="591.7" y="140" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 0.652, eGFR 57.5: 1.19</title></rect>
<rect x="600" y="140" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0.783, eGFR 57.5: 1.21</title></rect>
<rect x="608.3" y="140" width="8.633" height="8.633" fill="#e2edf8"><title>CP_severity 0.913, eGFR 57.5: 1.23</title></rect>
<rect x="616.7" y="140" width="8.633" height="8.633" fill="#dfecf7"><title>CP_severity 1.043, eGFR 57.5: 1.26</title></rect>
<rect x="625" y="140" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 1.174, eGFR 57.5: 1.3</title></rect>
<rect x="633.3" y="140" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.304, eGFR 57.5: 1.33</title></rect>
<rect x="641.7" y="140" width="8.633" height="8.633" fill="#d6e5f4"><title>CP_severity 1.435, eGFR 57.5: 1.37</title></rect>
<rect x="650" y="140" width="8.633" height="8.633" fill="#d2e3f3"><title>CP_severity 1.565, eGFR 57.5: 1.41</title></rect>
<rect x="658.3" y="140" width="8.633" height="8.633" fill="#cee0f2"><title>CP_severity 1.696, eGFR 57.5: 1.45</title></rect>
<rect x="666.7" y="140" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 1.826, eGFR 57.5: 1.5</title></rect>
<rect x="675" y="140" width="8.633" height="8.633" fill="#c5dbef"><title>CP_severity 1.957, eGFR 57.5: 1.55</title></rect>
<rect x="683.3" y="140" width="8.633" height="8.633" fill="#bcd6ec"><title>CP_severity 2.087, eGFR 57.5: 1.6</title></rect>
<rect x="691.7" y="140" width="8.633" height="8.633" fill="#b3d1ea"><title>CP_severity 2.217, eGFR 57.5: 1.66</title></rect>
<rect x="700" y="140" width="8.633" height="8.633" fill="#a8cce7"><title>CP_severity 2.348, eGFR 57.5: 1.72</title></rect>
<rect x="708.3" y="140" width="8.633" height="8.633" fill="#9dc7e4"><title>CP_severity 2.478, eGFR 57.5: 1.78</title></rect>
<rect x="716.7" y="140" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.609, eGFR 57.5: 1.86</title></rect>
<rect x="725" y="140" width="8.633" height="8.633" fill="#84badd"><title>CP_severity 2.739, eGFR 57.5: 1.93</title></rect>
<rect x="733.3" y="140" width="8.633" height="8.633" fill="#76b3d9"><title>CP_severity 2.87, eGFR 57.5: 2.02</title></rect>
<rect x="741.7" y="140" width="8.633" height="8.633" fill="#67abd4"><title>CP_severity 3, eGFR 57.5: 2.11</title></rect>
<rect x="550" y="148.3" width="8.633" height="8.633" fill="#edf5fc"><title>CP_severity 0, eGFR 52.5: 1.11</title></rect>
<rect x="558.3" y="148.3" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0.13, eGFR 52.5: 1.13</title></rect>
<rect x="566.7" y="148.3" width="8.633" height="8.633" fill="#eaf2fb"><title>CP_severity 0.261, eGFR 52.5: 1.15</title></rect>
<rect x="575" y="148.3" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 0.391, eGFR 52.5: 1.17</title></rect>
<rect x="583.3" y="148.3" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 0.522, eGFR 52.5: 1.19</title></rect>
<rect x="591.7" y="148.3" width="8.633" height="8.633" fill="#e4eff9"><title>CP_severity 0.652, eGFR 52.5: 1.21</title></rect>
<rect x="600" y="148.3" width="8.633" height="8.633" fill="#e2edf8"><title>CP_severity 0.783, eGFR 52.5: 1.23</title></rect>
<rect x="608.3" y="148.3" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 0.913, eGFR 52.5: 1.25</title></rect>
<rect x="616.7" y="148.3" width="8.633" height="8.633" fill="#deeaf7"><title>CP_severity 1.043, eGFR 52.5: 1.28</title></rect>
<rect x="625" y="148.3" width="8.633" height="8.633" fill="#dae8f6"><title>CP_severity 1.174, eGFR 52.5: 1.32</title></rect>
<rect x="633.3" y="148.3" width="8.633" height="8.633" fill="#d7e6f5"><title>CP_severity 1.304, eGFR 52.5: 1.35</title></rect>
<rect x="641.7" y="148.3" width="8.633" height="8.633" fill="#d4e4f3"><title>CP_severity 1.435, eGFR 52.5: 1.39</title></rect>
<rect x="650" y="148.3" width="8.633" height="8.633" fill="#d0e1f2"><title>CP_severity 1.565, eGFR 52.5: 1.43</title></rect>
<rect x="658.3" y="148.3" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 1.696, eGFR 52.5: 1.48</title></rect>
<rect x="666.7" y="148.3" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 1.826, eGFR 52.5: 1.52</title></rect>
<rect x="675" y="148.3" width="8.633" height="8.633" fill="#c1d9ee"><title>CP_severity 1.957, eGFR 52.5: 1.57</title></rect>
<rect x="683.3" y="148.3" width="8.633" height="8.633" fill="#b8d4eb"><title>CP_severity 2.087, eGFR 52.5: 1.62</title></rect>
<rect x="691.7" y="148.3" width="8.633" height="8.633" fill="#aecfe9"><title>CP_severity 2.217, eGFR 52.5: 1.68</title></rect>
<rect x="700" y="148.3" width="8.633" height="8.633" fill="#a4cae6"><title>CP_severity 2.348, eGFR 52.5: 1.74</title></rect>
<rect x="708.3" y="148.3" width="8.633" height="8.633" fill="#99c5e3"><title>CP_severity 2.478, eGFR 52.5: 1.81</title></rect>
<rect x="716.7" y="148.3" width="8.633" height="8.633" fill="#8cbfdf"><title>CP_severity 2.609, eGFR 52.5: 1.88</title></rect>
<rect x="725" y="148.3" width="8.633" height="8.633" fill="#7fb8dc"><title>CP_severity 2.739, eGFR 52.5: 1.96</title></rect>
<rect x="733.3" y="148.3" width="8.633" height="8.633" fill="#71b1d8"><title>CP_severity 2.87, eGFR 52.5: 2.05</title></rect>
<rect x="741.7" y="148.3" width="8.633" height="8.633" fill="#63a8d3"><title>CP_severity 3, eGFR 52.5: 2.14</title></rect>
<rect x="550" y="156.7" width="8.633" height="8.633" fill="#ecf4fb"><title>CP_severity 0, eGFR 47.5: 1.12</title></rect>
<rect x="558.3" y="156.7" width="8.633" height="8.633" fill="#eaf2fb"><title>CP_severity 0.13, eGFR 47.5: 1.14</title></rect>
<rect x="566.7" y="156.7" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 0.261, eGFR 47.5: 1.16</title></rect>
<rect x="575" y="156.7" width="8.633" height="8.633" fill="#e6f0fa"><title>CP_severity 0.391, eGFR 47.5: 1.18</title></rect>
<rect x="583.3" y="156.7" width="8.633" height="8.633" fill="#e4eff9"><title>CP_severity 0.522, eGFR 47.5: 1.21</title></rect>
<rect x="591.7" y="156.7" width="8.633" height="8.633" fill="#e2eef8"><title>CP_severity 0.652, eGFR 47.5: 1.23</title></rect>
<rect x="600" y="156.7" width="8.633" height="8.633" fill="#e0ecf8"><title>CP_severity 0.783, eGFR 47.5: 1.25</title></rect>
<rect x="608.3" y="156.7" width="8.633" height="8.633" fill="#deebf7"><title>CP_severity 0.913, eGFR 47.5: 1.27</title></rect>
<rect x="616.7" y="156.7" width="8.633" height="8.633" fill="#dce9f6"><title>CP_severity 1.043, eGFR 47.5: 1.3</title></rect>
<rect x="625" y="156.7" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 1.174, eGFR 47.5: 1.34</title></rect>
<rect x="633.3" y="156.7" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 1.304, eGFR 47.5: 1.37</title></rect>
<rect x="641.7" y="156.7" width="8.633" height="8.633" fill="#d2e3f3"><title>CP_severity 1.435, eGFR 47.5: 1.41</title></rect>
<rect x="650" y="156.7" width="8.633" height="8.633" fill="#cee0f2"><title>CP_severity 1.565, eGFR 47.5: 1.45</title></rect>
<rect x="658.3" y="156.7" width="8.633" height="8.633" fill="#cadef0"><title>CP_severity 1.696, eGFR 47.5: 1.5</title></rect>
<rect x="666.7" y="156.7" width="8.633" height="8.633" fill="#c5dbef"><title>CP_severity 1.826, eGFR 47.5: 1.54</title></rect>
<rect x="675" y="156.7" width="8.633" height="8.633" fill="#bdd7ed"><title>CP_severity 1.957, eGFR 47.5: 1.59</title></rect>
<rect x="683.3" y="156.7" width="8.633" height="8.633" fill="#b4d2ea"><title>CP_severity 2.087, eGFR 47.5: 1.65</title></rect>
<rect x="691.7" y="156.7" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 2.217, eGFR 47.5: 1.71</title></rect>
<rect x="700" y="156.7" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 2.348, eGFR 47.5: 1.77</title></rect>
<rect x="708.3" y="156.7" width="8.633" height="8.633" fill="#94c2e1"><title>CP_severity 2.478, eGFR 47.5: 1.84</title></rect>
<rect x="716.7" y="156.7" width="8.633" height="8.633" fill="#88bcde"><title>CP_severity 2.609, eGFR 47.5: 1.91</title></rect>
<rect x="725" y="156.7" width="8.633" height="8.633" fill="#7ab6da"><title>CP_severity 2.739, eGFR 47.5: 1.99</title></rect>
<rect x="733.3" y="156.7" width="8.633" height="8.633" fill="#6caed6"><title>CP_severity 2.87, eGFR 47.5: 2.08</title></rect>
<rect x="741.7" y="156.7" width="8.633" height="8.633" fill="#5fa4d1"><title>CP_severity 3, eGFR 47.5: 2.17</title></rect>
<rect x="550" y="165" width="8.633" height="8.633" fill="#eaf2fb"><title>CP_severity 0, eGFR 42.5: 1.15</title></rect>
<rect x="558.3" y="165" width="8.633" height="8.633" fill="#e8f1fa"><title>CP_severity 0.13, eGFR 42.5: 1.17</title></rect>
<rect x="566.7" y="165" width="8.633" height="8.633" fill="#e6f0f9"><title>CP_severity 0.261, eGFR 42.5: 1.19</title></rect>
<rect x="575" y="165" width="8.633" height="8.633" fill="#e4eff9"><title>CP_severity 0.391, eGFR 42.5: 1.21</title></rect>
<rect x="583.3" y="165" width="8.633" height="8.633" fill="#e2edf8"><title>CP_severity 0.522, eGFR 42.5: 1.23</title></rect>
<rect x="591.7" y="165" width="8.633" height="8.633" fill="#e0ecf8"><title>CP_severity 0.652, eGFR 42.5: 1.25</title></rect>
<rect x="600" y="165" width="8.633" height="8.633" fill="#deebf7"><title>CP_severity 0.783, eGFR 42.5: 1.28</title></rect>
<rect x="608.3" y="165" width="8.633" height="8.633" fill="#dce9f6"><title>CP_severity 0.913, eGFR 42.5: 1.3</title></rect>
<rect x="616.7" y="165" width="8.633" height="8.633" fill="#d9e8f5"><title>CP_severity 1.043, eGFR 42.5: 1.33</title></rect>
<rect x="625" y="165" width="8.633" height="8.633" fill="#d6e5f4"><title>CP_severity 1.174, eGFR 42.5: 1.36</title></rect>
<rect x="633.3" y="165" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 1.304, eGFR 42.5: 1.4</title></rect>
<rect x="641.7" y="165" width="8.633" height="8.633" fill="#cfe1f2"><title>CP_severity 1.435, eGFR 42.5: 1.44</title></rect>
<rect x="650" y="165" width="8.633" height="8.633" fill="#cbdef1"><title>CP_severity 1.565, eGFR 42.5: 1.48</title></rect>
<rect x="658.3" y="165" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 1.696, eGFR 42.5: 1.53</title></rect>
<rect x="666.7" y="165" width="8.633" height="8.633" fill="#c0d8ed"><title>CP_severity 1.826, eGFR 42.5: 1.58</title></rect>
<rect x="675" y="165" width="8.633" height="8.633" fill="#b7d4eb"><title>CP_severity 1.957, eGFR 42.5: 1.63</title></rect>
<rect x="683.3" y="165" width="8.633" height="8.633" fill="#aecfe8"><title>CP_severity 2.087, eGFR 42.5: 1.68</title></rect>
<rect x="691.7" y="165" width="8.633" height="8.633" fill="#a4cae6"><title>CP_severity 2.217, eGFR 42.5: 1.74</title></rect>
<rect x="700" y="165" width="8.633" height="8.633" fill="#99c5e3"><title>CP_severity 2.348, eGFR 42.5: 1.81</title></rect>
<rect x="708.3" y="165" width="8.633" height="8.633" fill="#8ebfe0"><title>CP_severity 2.478, eGFR 42.5: 1.88</title></rect>
<rect x="716.7" y="165" width="8.633" height="8.633" fill="#81b9dc"><title>CP_severity 2.609, eGFR 42.5: 1.95</title></rect>
<rect x="725" y="165" width="8.633" height="8.633" fill="#73b2d8"><title>CP_severity 2.739, eGFR 42.5: 2.03</title></rect>
<rect x="733.3" y="165" width="8.633" height="8.633" fill="#66aad4"><title>CP_severity 2.87, eGFR 42.5: 2.12</title></rect>
<rect x="741.7" y="165" width="8.633" height="8.633" fill="#599fce"><title>CP_severity 3, eGFR 42.5: 2.22</title></rect>
<rect x="550" y="173.3" width="8.633" height="8.633" fill="#e7f0fa"><title>CP_severity 0, eGFR 37.5: 1.18</title></rect>
<rect x="558.3" y="173.3" width="8.633" height="8.633" fill="#e5eff9"><title>CP_severity 0.13, eGFR 37.5: 1.2</title></rect>
<rect x="566.7" y="173.3" width="8.633" height="8.633" fill="#e3eef8"><title>CP_severity 0.261, eGFR 37.5: 1.22</title></rect>
<rect x="575" y="173.3" width="8.633" height="8.633" fill="#e1edf8"><title>CP_severity 0.391, eGFR 37.5: 1.24</title></rect>
<rect x="583.3" y="173.3" width="8.633" height="8.633" fill="#dfebf7"><title>CP_severity 0.522, eGFR 37.5: 1.27</title></rect>
<rect x="591.7" y="173.3" width="8.633" height="8.633" fill="#ddeaf6"><title>CP_severity 0.652, eGFR 37.5: 1.29</title></rect>
<rect x="600" y="173.3" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.783, eGFR 37.5: 1.31</title></rect>
<rect x="608.3" y="173.3" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 0.913, eGFR 37.5: 1.34</title></rect>
<rect x="616.7" y="173.3" width="8.633" height="8.633" fill="#d6e5f4"><title>CP_severity 1.043, eGFR 37.5: 1.37</title></rect>
<rect x="625" y="173.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 1.174, eGFR 37.5: 1.4</title></rect>
<rect x="633.3" y="173.3" width="8.633" height="8.633" fill="#cfe1f2"><title>CP_severity 1.304, eGFR 37.5: 1.44</title></rect>
<rect x="641.7" y="173.3" width="8.633" height="8.633" fill="#cbdff1"><title>CP_severity 1.435, eGFR 37.5: 1.48</title></rect>
<rect x="650" y="173.3" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 1.565, eGFR 37.5: 1.53</title></rect>
<rect x="658.3" y="173.3" width="8.633" height="8.633" fill="#c1d8ee"><title>CP_severity 1.696, eGFR 37.5: 1.57</title></rect>
<rect x="666.7" y="173.3" width="8.633" height="8.633" fill="#b9d4eb"><title>CP_severity 1.826, eGFR 37.5: 1.62</title></rect>
<rect x="675" y="173.3" width="8.633" height="8.633" fill="#b0d0e9"><title>CP_severity 1.957, eGFR 37.5: 1.67</title></rect>
<rect x="683.3" y="173.3" width="8.633" height="8.633" fill="#a6cbe6"><title>CP_severity 2.087, eGFR 37.5: 1.73</title></rect>
<rect x="691.7" y="173.3" width="8.633" height="8.633" fill="#9cc6e3"><title>CP_severity 2.217, eGFR 37.5: 1.79</title></rect>
<rect x="700" y="173.3" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 2.348, eGFR 37.5: 1.86</title></rect>
<rect x="708.3" y="173.3" width="8.633" height="8.633" fill="#85bbdd"><title>CP_severity 2.478, eGFR 37.5: 1.93</title></rect>
<rect x="716.7" y="173.3" width="8.633" height="8.633" fill="#78b5da"><title>CP_severity 2.609, eGFR 37.5: 2</title></rect>
<rect x="725" y="173.3" width="8.633" height="8.633" fill="#6baed6"><title>CP_severity 2.739, eGFR 37.5: 2.09</title></rect>
<rect x="733.3" y="173.3" width="8.633" height="8.633" fill="#5ea3d0"><title>CP_severity 2.87, eGFR 37.5: 2.18</title></rect>
<rect x="741.7" y="173.3" width="8.633" height="8.633" fill="#5198ca"><title>CP_severity 3, eGFR 37.5: 2.27</title></rect>
<rect x="550" y="181.7" width="8.633" height="8.633" fill="#e4eef9"><title>CP_severity 0, eGFR 32.5: 1.21</title></rect>
<rect x="558.3" y="181.7" width="8.633" height="8.633" fill="#e2edf8"><title>CP_severity 0.13, eGFR 32.5: 1.24</title></rect>
<rect x="566.7" y="181.7" width="8.633" height="8.633" fill="#e0ecf7"><title>CP_severity 0.261, eGFR 32.5: 1.26</title></rect>
<rect x="575" y="181.7" width="8.633" height="8.633" fill="#deebf7"><title>CP_severity 0.391, eGFR 32.5: 1.28</title></rect>
<rect x="583.3" y="181.7" width="8.633" height="8.633" fill="#dce9f6"><title>CP_severity 0.522, eGFR 32.5: 1.3</title></rect>
<rect x="591.7" y="181.7" width="8.633" height="8.633" fill="#dae8f5"><title>CP_severity 0.652, eGFR 32.5: 1.32</title></rect>
<rect x="600" y="181.7" width="8.633" height="8.633" fill="#d7e6f5"><title>CP_severity 0.783, eGFR 32.5: 1.35</title></rect>
<rect x="608.3" y="181.7" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 0.913, eGFR 32.5: 1.37</title></rect>
<rect x="616.7" y="181.7" width="8.633" height="8.633" fill="#d2e3f3"><title>CP_severity 1.043, eGFR 32.5: 1.4</title></rect>
<rect x="625" y="181.7" width="8.633" height="8.633" fill="#cfe1f2"><title>CP_severity 1.174, eGFR 32.5: 1.44</title></rect>
<rect x="633.3" y="181.7" width="8.633" height="8.633" fill="#cbdff1"><title>CP_severity 1.304, eGFR 32.5: 1.48</title></rect>
<rect x="641.7" y="181.7" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 1.435, eGFR 32.5: 1.52</title></rect>
<rect x="650" y="181.7" width="8.633" height="8.633" fill="#c2d9ee"><title>CP_severity 1.565, eGFR 32.5: 1.57</title></rect>
<rect x="658.3" y="181.7" width="8.633" height="8.633" fill="#bad5ec"><title>CP_severity 1.696, eGFR 32.5: 1.62</title></rect>
<rect x="666.7" y="181.7" width="8.633" height="8.633" fill="#b1d1e9"><title>CP_severity 1.826, eGFR 32.5: 1.67</title></rect>
<rect x="675" y="181.7" width="8.633" height="8.633" fill="#a8cce7"><title>CP_severity 1.957, eGFR 32.5: 1.72</title></rect>
<rect x="683.3" y="181.7" width="8.633" height="8.633" fill="#9ec7e4"><title>CP_severity 2.087, eGFR 32.5: 1.78</title></rect>
<rect x="691.7" y="181.7" width="8.633" height="8.633" fill="#94c2e1"><title>CP_severity 2.217, eGFR 32.5: 1.84</title></rect>
<rect x="700" y="181.7" width="8.633" height="8.633" fill="#89bdde"><title>CP_severity 2.348, eGFR 32.5: 1.91</title></rect>
<rect x="708.3" y="181.7" width="8.633" height="8.633" fill="#7cb7db"><title>CP_severity 2.478, eGFR 32.5: 1.98</title></rect>
<rect x="716.7" y="181.7" width="8.633" height="8.633" fill="#6fb0d7"><title>CP_severity 2.609, eGFR 32.5: 2.06</title></rect>
<rect x="725" y="181.7" width="8.633" height="8.633" fill="#63a7d2"><title>CP_severity 2.739, eGFR 32.5: 2.14</title></rect>
<rect x="733.3" y="181.7" width="8.633" height="8.633" fill="#569dcd"><title>CP_severity 2.87, eGFR 32.5: 2.23</title></rect>
<rect x="741.7" y="181.7" width="8.633" height="8.633" fill="#4892c7"><title>CP_severity 3, eGFR 32.5: 2.34</title></rect>
<rect x="550" y="190" width="8.633" height="8.633" fill="#e0ecf8"><title>CP_severity 0, eGFR 27.5: 1.25</title></rect>
<rect x="558.3" y="190" width="8.633" height="8.633" fill="#deebf7"><title>CP_severity 0.13, eGFR 27.5: 1.27</title></rect>
<rect x="566.7" y="190" width="8.633" height="8.633" fill="#dceaf6"><title>CP_severity 0.261, eGFR 27.5: 1.29</title></rect>
<rect x="575" y="190" width="8.633" height="8.633" fill="#dae8f6"><title>CP_severity 0.391, eGFR 27.5: 1.32</title></rect>
<rect x="583.3" y="190" width="8.633" height="8.633" fill="#d8e7f5"><title>CP_severity 0.522, eGFR 27.5: 1.34</title></rect>
<rect x="591.7" y="190" width="8.633" height="8.633" fill="#d6e6f4"><title>CP_severity 0.652, eGFR 27.5: 1.36</title></rect>
<rect x="600" y="190" width="8.633" height="8.633" fill="#d4e4f4"><title>CP_severity 0.783, eGFR 27.5: 1.39</title></rect>
<rect x="608.3" y="190" width="8.633" height="8.633" fill="#d2e3f3"><title>CP_severity 0.913, eGFR 27.5: 1.41</title></rect>
<rect x="616.7" y="190" width="8.633" height="8.633" fill="#cfe1f2"><title>CP_severity 1.043, eGFR 27.5: 1.44</title></rect>
<rect x="625" y="190" width="8.633" height="8.633" fill="#cbdef1"><title>CP_severity 1.174, eGFR 27.5: 1.48</title></rect>
<rect x="633.3" y="190" width="8.633" height="8.633" fill="#c8dcf0"><title>CP_severity 1.304, eGFR 27.5: 1.52</title></rect>
<rect x="641.7" y="190" width="8.633" height="8.633" fill="#c2d9ee"><title>CP_severity 1.435, eGFR 27.5: 1.57</title></rect>
<rect x="650" y="190" width="8.633" height="8.633" fill="#bad5ec"><title>CP_severity 1.565, eGFR 27.5: 1.61</title></rect>
<rect x="658.3" y="190" width="8.633" height="8.633" fill="#b2d1e9"><title>CP_severity 1.696, eGFR 27.5: 1.66</title></rect>
<rect x="666.7" y="190" width="8.633" height="8.633" fill="#a9cde7"><title>CP_severity 1.826, eGFR 27.5: 1.71</title></rect>
<rect x="675" y="190" width="8.633" height="8.633" fill="#a0c8e5"><title>CP_severity 1.957, eGFR 27.5: 1.77</title></rect>
<rect x="683.3" y="190" width="8.633" height="8.633" fill="#96c3e2"><title>CP_severity 2.087, eGFR 27.5: 1.83</title></rect>
<rect x="691.7" y="190" width="8.633" height="8.633" fill="#8bbedf"><title>CP_severity 2.217, eGFR 27.5: 1.89</title></rect>
<rect x="700" y="190" width="8.633" height="8.633" fill="#80b8dc"><title>CP_severity 2.348, eGFR 27.5: 1.96</title></rect>
<rect x="708.3" y="190" width="8.633" height="8.633" fill="#73b2d8"><title>CP_severity 2.478, eGFR 27.5: 2.03</title></rect>
<rect x="716.7" y="190" width="8.633" height="8.633" fill="#67aad4"><title>CP_severity 2.609, eGFR 27.5: 2.11</title></rect>
<rect x="725" y="190" width="8.633" height="8.633" fill="#5ba1cf"><title>CP_severity 2.739, eGFR 27.5: 2.2</title></rect>
<rect x="733.3" y="190" width="8.633" height="8.633" fill="#4e96c9"><title>CP_severity 2.87, eGFR 27.5: 2.3</title></rect>
<rect x="741.7" y="190" width="8.633" height="8.633" fill="#408ac3"><title>CP_severity 3, eGFR 27.5: 2.4</title></rect>
<rect x="550" y="198.3" width="8.633" height="8.633" fill="#ddeaf7"><title>CP_severity 0, eGFR 22.5: 1.29</title></rect>
<rect x="558.3" y="198.3" width="8.633" height="8.633" fill="#dbe9f6"><title>CP_severity 0.13, eGFR 22.5: 1.31</title></rect>
<rect x="566.7" y="198.3" width="8.633" height="8.633" fill="#d9e7f5"><title>CP_severity 0.261, eGFR 22.5: 1.33</title></rect>
<rect x="575" y="198.3" width="8.633" height="8.633" fill="#d7e6f5"><title>CP_severity 0.391, eGFR 22.5: 1.35</title></rect>
<rect x="583.3" y="198.3" width="8.633" height="8.633" fill="#d5e5f4"><title>CP_severity 0.522, eGFR 22.5: 1.38</title></rect>
<rect x="591.7" y="198.3" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0.652, eGFR 22.5: 1.4</title></rect>
<rect x="600" y="198.3" width="8.633" height="8.633" fill="#d0e2f2"><title>CP_severity 0.783, eGFR 22.5: 1.43</title></rect>
<rect x="608.3" y="198.3" width="8.633" height="8.633" fill="#cee0f2"><title>CP_severity 0.913, eGFR 22.5: 1.46</title></rect>
<rect x="616.7" y="198.3" width="8.633" height="8.633" fill="#cbdef1"><title>CP_severity 1.043, eGFR 22.5: 1.49</title></rect>
<rect x="625" y="198.3" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 1.174, eGFR 22.5: 1.53</title></rect>
<rect x="633.3" y="198.3" width="8.633" height="8.633" fill="#c1d9ee"><title>CP_severity 1.304, eGFR 22.5: 1.57</title></rect>
<rect x="641.7" y="198.3" width="8.633" height="8.633" fill="#bad5ec"><title>CP_severity 1.435, eGFR 22.5: 1.61</title></rect>
<rect x="650" y="198.3" width="8.633" height="8.633" fill="#b2d1ea"><title>CP_severity 1.565, eGFR 22.5: 1.66</title></rect>
<rect x="658.3" y="198.3" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 1.696, eGFR 22.5: 1.71</title></rect>
<rect x="666.7" y="198.3" width="8.633" height="8.633" fill="#a1c9e5"><title>CP_severity 1.826, eGFR 22.5: 1.76</title></rect>
<rect x="675" y="198.3" width="8.633" height="8.633" fill="#97c4e2"><title>CP_severity 1.957, eGFR 22.5: 1.82</title></rect>
<rect x="683.3" y="198.3" width="8.633" height="8.633" fill="#8dbfdf"><title>CP_severity 2.087, eGFR 22.5: 1.88</title></rect>
<rect x="691.7" y="198.3" width="8.633" height="8.633" fill="#82b9dc"><title>CP_severity 2.217, eGFR 22.5: 1.95</title></rect>
<rect x="700" y="198.3" width="8.633" height="8.633" fill="#76b4d9"><title>CP_severity 2.348, eGFR 22.5: 2.02</title></rect>
<rect x="708.3" y="198.3" width="8.633" height="8.633" fill="#6aadd5"><title>CP_severity 2.478, eGFR 22.5: 2.09</title></rect>
<rect x="716.7" y="198.3" width="8.633" height="8.633" fill="#5ea4d0"><title>CP_severity 2.609, eGFR 22.5: 2.17</title></rect>
<rect x="725" y="198.3" width="8.633" height="8.633" fill="#529acb"><title>CP_severity 2.739, eGFR 22.5: 2.26</title></rect>
<rect x="733.3" y="198.3" width="8.633" height="8.633" fill="#458fc5"><title>CP_severity 2.87, eGFR 22.5: 2.36</title></rect>
<rect x="741.7" y="198.3" width="8.633" height="8.633" fill="#3783bf"><title>CP_severity 3, eGFR 22.5: 2.47</title></rect>
<rect x="550" y="206.7" width="8.633" height="8.633" fill="#d3e3f3"><title>CP_severity 0, eGFR 17.5: 1.4</title></rect>
<rect x="558.3" y="206.7" width="8.633" height="8.633" fill="#d1e2f2"><title>CP_severity 0.13, eGFR 17.5: 1.42</title></rect>
<rect x="566.7" y="206.7" width="8.633" height="8.633" fill="#cee0f2"><title>CP_severity 0.261, eGFR 17.5: 1.45</title></rect>
<rect x="575" y="206.7" width="8.633" height="8.633" fill="#ccdff1"><title>CP_severity 0.391, eGFR 17.5: 1.47</title></rect>
<rect x="583.3" y="206.7" width="8.633" height="8.633" fill="#caddf0"><title>CP_severity 0.522, eGFR 17.5: 1.5</title></rect>
<rect x="591.7" y="206.7" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 0.652, eGFR 17.5: 1.53</title></rect>
<rect x="600" y="206.7" width="8.633" height="8.633" fill="#c4daee"><title>CP_severity 0.783, eGFR 17.5: 1.55</title></rect>
<rect x="608.3" y="206.7" width="8.633" height="8.633" fill="#bfd8ed"><title>CP_severity 0.913, eGFR 17.5: 1.58</title></rect>
<rect x="616.7" y="206.7" width="8.633" height="8.633" fill="#b9d5ec"><title>CP_severity 1.043, eGFR 17.5: 1.62</title></rect>
<rect x="625" y="206.7" width="8.633" height="8.633" fill="#b2d1ea"><title>CP_severity 1.174, eGFR 17.5: 1.66</title></rect>
<rect x="633.3" y="206.7" width="8.633" height="8.633" fill="#aacde7"><title>CP_severity 1.304, eGFR 17.5: 1.7</title></rect>
<rect x="641.7" y="206.7" width="8.633" height="8.633" fill="#a2c9e5"><title>CP_severity 1.435, eGFR 17.5: 1.75</title></rect>
<rect x="650" y="206.7" width="8.633" height="8.633" fill="#9ac5e3"><title>CP_severity 1.565, eGFR 17.5: 1.8</title></rect>
<rect x="658.3" y="206.7" width="8.633" height="8.633" fill="#91c1e0"><title>CP_severity 1.696, eGFR 17.5: 1.86</title></rect>
<rect x="666.7" y="206.7" width="8.633" height="8.633" fill="#87bcde"><title>CP_severity 1.826, eGFR 17.5: 1.91</title></rect>
<rect x="675" y="206.7" width="8.633" height="8.633" fill="#7db7db"><title>CP_severity 1.957, eGFR 17.5: 1.97</title></rect>
<rect x="683.3" y="206.7" width="8.633" height="8.633" fill="#72b2d8"><title>CP_severity 2.087, eGFR 17.5: 2.04</title></rect>
<rect x="691.7" y="206.7" width="8.633" height="8.633" fill="#67abd4"><title>CP_severity 2.217, eGFR 17.5: 2.11</title></rect>
<rect x="700" y="206.7" width="8.633" height="8.633" fill="#5da2d0"><title>CP_severity 2.348, eGFR 17.5: 2.19</title></rect>
<rect x="708.3" y="206.7" width="8.633" height="8.633" fill="#5299cb"><title>CP_severity 2.478, eGFR 17.5: 2.27</title></rect>
<rect x="716.7" y="206.7" width="8.633" height="8.633" fill="#458fc5"><title>CP_severity 2.609, eGFR 17.5: 2.36</title></rect>
<rect x="725" y="206.7" width="8.633" height="8.633" fill="#3884bf"><title>CP_severity 2.739, eGFR 17.5: 2.45</title></rect>
<rect x="733.3" y="206.7" width="8.633" height="8.633" fill="#2a79b9"><title>CP_severity 2.87, eGFR 17.5: 2.56</title></rect>
<rect x="741.7" y="206.7" width="8.633" height="8.633" fill="#1f6caf"><title>CP_severity 3, eGFR 17.5: 2.67</title></rect>
<rect x="550" y="215" width="8.633" height="8.633" fill="#c7dcef"><title>CP_severity 0, eGFR 12.5: 1.53</title></rect>
<rect x="558.3" y="215" width="8.633" height="8.633" fill="#c3daee"><title>CP_severity 0.13, eGFR 12.5: 1.56</title></rect>
<rect x="566.7" y="215" width="8.633" height="8.633" fill="#bfd7ed"><title>CP_severity 0.261, eGFR 12.5: 1.58</title></rect>
<rect x="575" y="215" width="8.633" height="8.633" fill="#bad5ec"><title>CP_severity 0.391, eGFR 12.5: 1.61</title></rect>
<rect x="583.3" y="215" width="8.633" height="8.633" fill="#b6d3ea"><title>CP_severity 0.522, eGFR 12.5: 1.64</title></rect>
<rect x="591.7" y="215" width="8.633" height="8.633" fill="#b1d0e9"><title>CP_severity 0.652, eGFR 12.5: 1.67</title></rect>
<rect x="600" y="215" width="8.633" height="8.633" fill="#accee8"><title>CP_severity 0.783, eGFR 12.5: 1.7</title></rect>
<rect x="608.3" y="215" width="8.633" height="8.633" fill="#a6cbe6"><title>CP_severity 0.913, eGFR 12.5: 1.73</title></rect>
<rect x="616.7" y="215" width="8.633" height="8.633" fill="#a0c8e5"><title>CP_severity 1.043, eGFR 12.5: 1.77</title></rect>
<rect x="625" y="215" width="8.633" height="8.633" fill="#98c4e2"><title>CP_severity 1.174, eGFR 12.5: 1.81</title></rect>
<rect x="633.3" y="215" width="8.633" height="8.633" fill="#90c0e0"><title>CP_severity 1.304, eGFR 12.5: 1.86</title></rect>
<rect x="641.7" y="215" width="8.633" height="8.633" fill="#87bcde"><title>CP_severity 1.435, eGFR 12.5: 1.91</title></rect>
<rect x="650" y="215" width="8.633" height="8.633" fill="#7eb8db"><title>CP_severity 1.565, eGFR 12.5: 1.97</title></rect>
<rect x="658.3" y="215" width="8.633" height="8.633" fill="#75b3d9"><title>CP_severity 1.696, eGFR 12.5: 2.03</title></rect>
<rect x="666.7" y="215" width="8.633" height="8.633" fill="#6aadd6"><title>CP_severity 1.826, eGFR 12.5: 2.09</title></rect>
<rect x="675" y="215" width="8.633" height="8.633" fill="#61a6d2"><title>CP_severity 1.957, eGFR 12.5: 2.15</title></rect>
<rect x="683.3" y="215" width="8.633" height="8.633" fill="#589ecd"><title>CP_severity 2.087, eGFR 12.5: 2.22</title></rect>
<rect x="691.7" y="215" width="8.633" height="8.633" fill="#4d95c9"><title>CP_severity 2.217, eGFR 12.5: 2.3</title></rect>
<rect x="700" y="215" width="8.633" height="8.633" fill="#428cc4"><title>CP_severity 2.348, eGFR 12.5: 2.38</title></rect>
<rect x="708.3" y="215" width="8.633" height="8.633" fill="#3682be"><title>CP_severity 2.478, eGFR 12.5: 2.47</title></rect>
<rect x="716.7" y="215" width="8.633" height="8.633" fill="#2978b9"><title>CP_severity 2.609, eGFR 12.5: 2.57</title></rect>
<rect x="725" y="215" width="8.633" height="8.633" fill="#1f6caf"><title>CP_severity 2.739, eGFR 12.5: 2.67</title></rect>
<rect x="733.3" y="215" width="8.633" height="8.633" fill="#1a5ea0"><title>CP_severity 2.87, eGFR 12.5: 2.78</title></rect>
<rect x="741.7" y="215" width="8.633" height="8.633" fill="#14508f"><title>CP_severity 3, eGFR 12.5: 2.9</title></rect>
<rect x="550" y="223.3" width="8.633" height="8.633" fill="#aecfe8"><title>CP_severity 0, eGFR 7.5: 1.68</title></rect>
<rect x="558.3" y="223.3" width="8.633" height="8.633" fill="#a9cde7"><title>CP_severity 0.13, eGFR 7.5: 1.71</title></rect>
<rect x="566.7" y="223.3" width="8.633" height="8.633" fill="#a4cae6"><title>CP_severity 0.261, eGFR 7.5: 1.74</title></rect>
<rect x="575" y="223.3" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.391, eGFR 7.5: 1.77</title></rect>
<rect x="583.3" y="223.3" width="8.633" height="8.633" fill="#9ac5e3"><title>CP_severity 0.522, eGFR 7.5: 1.8</title></rect>
<rect x="591.7" y="223.3" width="8.633" height="8.633" fill="#95c3e2"><title>CP_severity 0.652, eGFR 7.5: 1.83</title></rect>
<rect x="600" y="223.3" width="8.633" height="8.633" fill="#8fc0e0"><title>CP_severity 0.783, eGFR 7.5: 1.87</title></rect>
<rect x="608.3" y="223.3" width="8.633" height="8.633" fill="#8abdde"><title>CP_severity 0.913, eGFR 7.5: 1.9</title></rect>
<rect x="616.7" y="223.3" width="8.633" height="8.633" fill="#83badd"><title>CP_severity 1.043, eGFR 7.5: 1.94</title></rect>
<rect x="625" y="223.3" width="8.633" height="8.633" fill="#7bb6da"><title>CP_severity 1.174, eGFR 7.5: 1.99</title></rect>
<rect x="633.3" y="223.3" width="8.633" height="8.633" fill="#72b1d8"><title>CP_severity 1.304, eGFR 7.5: 2.04</title></rect>
<rect x="641.7" y="223.3" width="8.633" height="8.633" fill="#69acd5"><title>CP_severity 1.435, eGFR 7.5: 2.1</title></rect>
<rect x="650" y="223.3" width="8.633" height="8.633" fill="#61a5d1"><title>CP_severity 1.565, eGFR 7.5: 2.16</title></rect>
<rect x="658.3" y="223.3" width="8.633" height="8.633" fill="#589ece"><title>CP_severity 1.696, eGFR 7.5: 2.22</title></rect>
<rect x="666.7" y="223.3" width="8.633" height="8.633" fill="#4f97c9"><title>CP_severity 1.826, eGFR 7.5: 2.29</title></rect>
<rect x="675" y="223.3" width="8.633" height="8.633" fill="#458fc5"><title>CP_severity 1.957, eGFR 7.5: 2.36</title></rect>
<rect x="683.3" y="223.3" width="8.633" height="8.633" fill="#3b86c0"><title>CP_severity 2.087, eGFR 7.5: 2.44</title></rect>
<rect x="691.7" y="223.3" width="8.633" height="8.633" fill="#2f7dbb"><title>CP_severity 2.217, eGFR 7.5: 2.52</title></rect>
<rect x="700" y="223.3" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 2.348, eGFR 7.5: 2.61</title></rect>
<rect x="708.3" y="223.3" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 2.478, eGFR 7.5: 2.7</title></rect>
<rect x="716.7" y="223.3" width="8.633" height="8.633" fill="#195b9c"><title>CP_severity 2.609, eGFR 7.5: 2.81</title></rect>
<rect x="725" y="223.3" width="8.633" height="8.633" fill="#144e8d"><title>CP_severity 2.739, eGFR 7.5: 2.92</title></rect>
<rect x="733.3" y="223.3" width="8.633" height="8.633" fill="#0e407d"><title>CP_severity 2.87, eGFR 7.5: 3.04</title></rect>
<rect x="741.7" y="223.3" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 3, eGFR 7.5: 3.17</title></rect>
<rect x="550" y="231.7" width="8.633" height="8.633" fill="#aecfe8"><title>CP_severity 0, eGFR 2.5: 1.68</title></rect>
<rect x="558.3" y="231.7" width="8.633" height="8.633" fill="#a9cde7"><title>CP_severity 0.13, eGFR 2.5: 1.71</title></rect>
<rect x="566.7" y="231.7" width="8.633" height="8.633" fill="#a4cae6"><title>CP_severity 0.261, eGFR 2.5: 1.74</title></rect>
<rect x="575" y="231.7" width="8.633" height="8.633" fill="#9fc8e4"><title>CP_severity 0.391, eGFR 2.5: 1.77</title></rect>
<rect x="583.3" y="231.7" width="8.633" height="8.633" fill="#9ac5e3"><title>CP_severity 0.522, eGFR 2.5: 1.8</title></rect>
<rect x="591.7" y="231.7" width="8.633" height="8.633" fill="#95c3e2"><title>CP_severity 0.652, eGFR 2.5: 1.83</title></rect>
<rect x="600" y="231.7" width="8.633" height="8.633" fill="#8fc0e0"><title>CP_severity 0.783, eGFR 2.5: 1.87</title></rect>
<rect x="608.3" y="231.7" width="8.633" height="8.633" fill="#8abdde"><title>CP_severity 0.913, eGFR 2.5: 1.9</title></rect>
<rect x="616.7" y="231.7" width="8.633" height="8.633" fill="#83badd"><title>CP_severity 1.043, eGFR 2.5: 1.94</title></rect>
<rect x="625" y="231.7" width="8.633" height="8.633" fill="#7bb6da"><title>CP_severity 1.174, eGFR 2.5: 1.99</title></rect>
<rect x="633.3" y="231.7" width="8.633" height="8.633" fill="#72b1d8"><title>CP_severity 1.304, eGFR 2.5: 2.04</title></rect>
<rect x="641.7" y="231.7" width="8.633" height="8.633" fill="#69acd5"><title>CP_severity 1.435, eGFR 2.5: 2.1</title></rect>
<rect x="650" y="231.7" width="8.633" height="8.633" fill="#61a5d1"><title>CP_severity 1.565, eGFR 2.5: 2.16</title></rect>
<rect x="658.3" y="231.7" width="8.633" height="8.633" fill="#589ece"><title>CP_severity 1.696, eGFR 2.5: 2.22</title></rect>
<rect x="666.7" y="231.7" width="8.633" height="8.633" fill="#4f97c9"><title>CP_severity 1.826, eGFR 2.5: 2.29</title></rect>
<rect x="675" y="231.7" width="8.633" height="8.633" fill="#458fc5"><title>CP_severity 1.957, eGFR 2.5: 2.36</title></rect>
<rect x="683.3" y="231.7" width="8.633" height="8.633" fill="#3b86c0"><title>CP_severity 2.087, eGFR 2.5: 2.44</title></rect>
<rect x="691.7" y="231.7" width="8.633" height="8.633" fill="#2f7dbb"><title>CP_severity 2.217, eGFR 2.5: 2.52</title></rect>
<rect x="700" y="231.7" width="8.633" height="8.633" fill="#2373b6"><title>CP_severity 2.348, eGFR 2.5: 2.61</title></rect>
<rect x="708.3" y="231.7" width="8.633" height="8.633" fill="#1d68aa"><title>CP_severity 2.478, eGFR 2.5: 2.7</title></rect>
<rect x="716.7" y="231.7" width="8.633" height="8.633" fill="#195b9c"><title>CP_severity 2.609, eGFR 2.5: 2.81</title></rect>
<rect x="725" y="231.7" width="8.633" height="8.633" fill="#144e8d"><title>CP_severity 2.739, eGFR 2.5: 2.92</title></rect>
<rect x="733.3" y="231.7" width="8.633" height="8.633" fill="#0e407d"><title>CP_severity 2.87, eGFR 2.5: 3.04</title></rect>
<rect x="741.7" y="231.7" width="8.633" height="8.633" fill="#08306b"><title>CP_severity 3, eGFR 2.5: 3.17</title></rect>
<line x1="546" x2="550" y1="90" y2="90" stroke="currentColor"/>
<text x="544" y="94" text-anchor="end" font-size="11">90</text>
<line x1="546" x2="550" y1="140" y2="140" stroke="currentColor"/>
<text x="544" y="144" text-anchor="end" font-size="11">60</text>
<line x1="546" x2="550" y1="190" y2="190" stroke="currentColor"/>
<text x="544" y="194" text-anchor="end" font-size="11">30</text>
<line x1="546" x2="550" y1="215" y2="215" stroke="currentColor"/>
<text x="544" y="219" text-anchor="end" font-size="11">15</text>
<line x1="554.2" x2="554.2" y1="240" y2="244" stroke="currentColor"/>
<text x="554.2" y="257" text-anchor="middle" font-size="11">0</text>
<line x1="618.1" x2="618.1" y1="240" y2="244" stroke="currentColor"/>
<text x="618.1" y="257" text-anchor="middle" font-size="11">A</text>
<line x1="681.9" x2="681.9" y1="240" y2="244" stroke="currentColor"/>
<text x="681.9" y="257" text-anchor="middle" font-size="11">B</text>
<line x1="745.8" x2="745.8" y1="240" y2="244" stroke="currentColor"/>
<text x="745.8" y="257" text-anchor="middle" font-size="11">C</text>
<text x="650" y="273" text-anchor="middle" font-size="12">CP_severity</text>
<rect x="550" y="282" width="10.3" height="8" fill="#f7fbff"/>
<rect x="560" y="282" width="10.3" height="8" fill="#edf4fc"/>
<rect x="570" y="282" width="10.3" height="8" fill="#e2eef8"/>
<rect x="580" y="282" width="10.3" height="8" fill="#d8e7f5"/>
<rect x="590" y="282" width="10.3" height="8" fill="#cee0f2"/>
<rect x="600" y="282" width="10.3" height="8" fill="#c1d9ee"/>
<rect x="610" y="282" width="10.3" height="8" fill="#aecfe8"/>
<rect x="620" y="282" width="10.3" height="8" fill="#9bc6e3"/>
<rect x="630" y="282" width="10.3" height="8" fill="#88bcde"/>
<rect x="640" y="282" width="10.3" height="8" fill="#75b3d9"/>
<rect x="650" y="282" width="10.3" height="8" fill="#63a8d3"/>
<rect x="660" y="282" width="10.3" height="8" fill="#549bcc"/>
<rect x="670" y="282" width="10.3" height="8" fill="#448ec5"/>
<rect x="680" y="282" width="10.3" height="8" fill="#3481be"/>
<rect x="690" y="282" width="10.3" height="8" fill="#2574b7"/>
<rect x="700" y="282" width="10.3" height="8" fill="#1d67a9"/>
<rect x="710" y="282" width="10.3" height="8" fill="#18599a"/>
<rect x="720" y="282" width="10.3" height="8" fill="#134b8a"/>
<rect x="730" y="282" width="10.3" height="8" fill="#0d3e7b"/>
<rect x="740" y="282" width="10.3" height="8" fill="#08306b"/>
<text x="550" y="304" font-size="11">1</text>
<text x="750" y="304" font-size="11" text-anchor="end">3.2</text>
</svg>
</div>
<!-- end generated -->

The heat maps show three combinations relative to healthy, every cell a point of a continuous scan (hover a cell for its value): the GFR over `CP_severity` and the eGFR of the renal function (hepatorenal), the renal blood flow over `HF_severity` and the eGFR (cardiorenal) and the fu of a drug bound to albumin with a reference fu of 0.05 over `CP_severity` and the eGFR.

A continuous scan sets the severities before the initialization of every point. `severity_values` and `severity_changes` of `PBPKSimulationExperiment` give the values of a point by the severities or their clinical measures (`cp_severity` or `meld`, `egfr` or `renal_function`, `hf_severity`), `impairment_dimension(id, ...)` a dimension of a sbmlsim scan whose severities change together (a course, e.g. a cirrhosis and its hepatorenal decline); a scan of two dimensions is the grid of both. `color_of(value)` of a classification gives the color of a continuous value, interpolated between the colors of its classes:

```python
import numpy as np

from pkpdsim.impairment import CHILD_PUGH_CLASSES
from pkpdsim.pbpk import PBPKSimulationExperiment

values = PBPKSimulationExperiment.severity_values(
    meld=[6.43, 10.0, 15.0, 25.0], egfr=[100.0, 70.0, 45.0, 20.0]
)
for cp, f in zip(values["CP_severity"], values["f_renal_function"]):
    color = CHILD_PUGH_CLASSES.color_of(cp)
    print(f"CP_severity {cp:.2f} {color}, f_renal_function {f:.3f}")
grid = np.linspace(0.0, 3.0, 7)
print([CHILD_PUGH_CLASSES.color_of(x) for x in grid])
```

```python
# not executed
from sbmlsim.simulation import Scan


class HepatorenalScan(HctzSimulationExperiment):
    """A dose of 25 mg over the grid of cirrhosis and renal function."""

    def simulations(self) -> dict[str, Scan]:
        dose = Simulation(end=48 * 60, preinit_changes={"PODOSE_hctz": Q(25, "mg")})
        cirrhosis = self.impairment_dimension(
            "cirrhosis", cp_severity=np.linspace(0, 3, 13)
        )
        kidney = self.impairment_dimension("kidney", egfr=np.linspace(120, 5, 24))
        return {"scan_hepatorenal": Scan(dose, [cirrhosis, kidney])}
```

The PK of every point is `scan_nca` of `pkpdsim.pk` ([Scans and sensitivity](../scans.md)), drawn over the severities in the colors of `color_of`.

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

## References

<!-- generated: references (scripts/pbpk_tables.py), do not edit -->
- <span id="ref-Bozkurt2021"></span>Bozkurt B, Coats AJS, Tsutsui H, Abdelhamid CM, Adamopoulos S, Albert N, Anker SD, Atherton J, et al. Universal definition and classification of heart failure: a report of the Heart Failure Society of America, Heart Failure Association of the European Society of Cardiology, Japanese Heart Failure Society and Writing Committee of the Universal Definition of Heart Failure. Eur J Heart Fail 23(3):352-380, 2021. [doi:10.1002/ejhf.2115](https://doi.org/10.1002/ejhf.2115)
- <span id="ref-Bricker1960"></span>Bricker NS, Morrin PA, Kime SW Jr. The pathologic physiology of chronic Bright's disease. An exposition of the intact nephron hypothesis. Am J Med 28:77-98, 1960. [doi:10.1016/0002-9343(60)90225-4](https://doi.org/10.1016/0002-9343(60)90225-4)
- <span id="ref-ChildTurcotte1964"></span>Child CG, Turcotte JG. Surgery and portal hypertension. In: Child CG (ed.), The Liver and Portal Hypertension. Major Probl Clin Surg 1:1-85, 1964.
- <span id="ref-Edginton2008"></span>Edginton AN, Willmann S. Physiology-based simulations of a pathological condition: prediction of pharmacokinetics in patients with liver cirrhosis. Clin Pharmacokinet 47(11):743-752, 2008. [doi:10.2165/00003088-200847110-00005](https://doi.org/10.2165/00003088-200847110-00005)
- <span id="ref-EMA2005Hepatic"></span>European Medicines Agency. Guideline on the evaluation of the pharmacokinetics of medicinal products in patients with impaired hepatic function. CPMP/EWP/2339/02, 2005.
- <span id="ref-EMA2015Renal"></span>European Medicines Agency. Guideline on the evaluation of the pharmacokinetics of medicinal products in patients with decreased renal function. EMA/CHMP/83874/2014, 2015.
- <span id="ref-FDA2003Hepatic"></span>U.S. Food and Drug Administration. Pharmacokinetics in Patients with Impaired Hepatic Function: Study Design, Data Analysis, and Impact on Dosing and Labeling. Guidance for Industry, 2003. https://www.fda.gov/media/71311/download
- <span id="ref-FDA2024Renal"></span>U.S. Food and Drug Administration. Pharmacokinetics in Patients with Impaired Renal Function - Study Design, Data Analysis, and Impact on Dosing. Guidance for Industry, 2024. https://www.fda.gov/media/78573/download
- <span id="ref-Gines2018"></span>Ginès P, Solà E, Angeli P, Wong F, Nadim MK, Kamath PS. Hepatorenal syndrome. Nat Rev Dis Primers 4(1):23, 2018. [doi:10.1038/s41572-018-0022-7](https://doi.org/10.1038/s41572-018-0022-7)
- <span id="ref-Heidenreich2022"></span>Heidenreich PA, Bozkurt B, Aguilar D, Allen LA, Byun JJ, Colvin MM, Deswal A, Drazner MH, et al. 2022 AHA/ACC/HFSA Guideline for the Management of Heart Failure. Circulation 145(18):e895-e1032, 2022. [doi:10.1161/CIR.0000000000001063](https://doi.org/10.1161/CIR.0000000000001063)
- <span id="ref-ICHE11"></span>ICH Topic E11. Clinical investigation of medicinal products in the paediatric population, Section 2.5 Age classification of paediatric patients. CPMP/ICH/2711/99, European Medicines Agency, Step 5, January 2001 (ICH Step 4, July 2000).
- <span id="ref-ICHE7"></span>ICH Harmonised Tripartite Guideline E7. Studies in support of special populations: geriatrics, Section IV. Step 4, 24 June 1993; and ICH E7 Questions and Answers, Q2, 2010 (the subgroups 65-74, 75-84 and 85 years and older).
- <span id="ref-Iwakiri2014"></span>Iwakiri Y. Pathophysiology of portal hypertension. Clin Liver Dis 18(2):281-291, 2014. [doi:10.1016/j.cld.2013.12.001](https://doi.org/10.1016/j.cld.2013.12.001)
- <span id="ref-Johnson2010"></span>Johnson TN, Boussery K, Rowland-Yeo K, Tucker GT, Rostami-Hodjegan A. A semi-mechanistic model to predict the effects of liver cirrhosis on drug clearance. Clin Pharmacokinet 49(3):189-206, 2010. [doi:10.2165/11318160-000000000-00000](https://doi.org/10.2165/11318160-000000000-00000)
- <span id="ref-Kamath2001"></span>Kamath PS, Wiesner RH, Malinchoc M, Kremers W, Therneau TM, Kosberg CL, D'Amico G, Dickson ER, Kim WR. A model to predict survival in patients with end-stage liver disease. Hepatology 33(2):464-470, 2001. [doi:10.1053/jhep.2001.22172](https://doi.org/10.1053/jhep.2001.22172)
- <span id="ref-KDIGO2012"></span>Kidney Disease: Improving Global Outcomes (KDIGO) CKD Work Group. KDIGO 2012 clinical practice guideline for the evaluation and management of chronic kidney disease. Chapter 1: Definition and classification of CKD. Kidney Int Suppl 3(1):19-62, 2013. [doi:10.1038/kisup.2012.64](https://doi.org/10.1038/kisup.2012.64)
- <span id="ref-Kim2008"></span>Kim WR, Biggins SW, Kremers WK, Wiesner RH, Kamath PS, Benson JT, Edwards E, Therneau TM. Hyponatremia and mortality among patients on the liver-transplant waiting list. N Engl J Med 359(10):1018-1026, 2008. [doi:10.1056/NEJMoa0801209](https://doi.org/10.1056/NEJMoa0801209)
- <span id="ref-Leithe1984"></span>Leithe ME, Margorien RD, Hermiller JB, Unverferth DV, Leier CV. Relationship between central hemodynamics and regional blood flow in normal subjects and in patients with congestive heart failure. Circulation 69(1):57-64, 1984. [doi:10.1161/01.CIR.69.1.57](https://doi.org/10.1161/01.CIR.69.1.57)
- <span id="ref-Levey2005"></span>Levey AS, Eckardt KU, Tsukamoto Y, Levin A, Coresh J, Rossert J, de Zeeuw D, Hostetter TH, Lameire N, Eknoyan G. Definition and classification of chronic kidney disease: a position statement from Kidney Disease: Improving Global Outcomes (KDIGO). Kidney Int 67(6):2089-2100, 2005. [doi:10.1111/j.1523-1755.2005.00365.x](https://doi.org/10.1111/j.1523-1755.2005.00365.x)
- <span id="ref-Malinchoc2000"></span>Malinchoc M, Kamath PS, Gordon FD, Peine CJ, Rank J, ter Borg PC. A model to predict poor survival in patients undergoing transjugular intrahepatic portosystemic shunts. Hepatology 31(4):864-871, 2000. [doi:10.1053/he.2000.5852](https://doi.org/10.1053/he.2000.5852)
- <span id="ref-Nolin2008"></span>Nolin TD, Naud J, Leblond FA, Pichette V. Emerging evidence of the impact of kidney disease on drug metabolism and transport. Clin Pharmacol Ther 83(6):898-903, 2008. [doi:10.1038/clpt.2008.59](https://doi.org/10.1038/clpt.2008.59)
- <span id="ref-NYHA1994"></span>The Criteria Committee of the New York Heart Association. Nomenclature and Criteria for Diagnosis of Diseases of the Heart and Great Vessels. 9th ed. Boston, Little, Brown & Co, 253-256, 1994.
- <span id="ref-Pugh1973"></span>Pugh RNH, Murray-Lyon IM, Dawson JL, Pietroni MC, Williams R. Transection of the oesophagus for bleeding oesophageal varices. Br J Surg 60(8):646-649, 1973. [doi:10.1002/bjs.1800600817](https://doi.org/10.1002/bjs.1800600817)
- <span id="ref-Ronco2008"></span>Ronco C, Haapio M, House AA, Anavekar N, Bellomo R. Cardiorenal syndrome. J Am Coll Cardiol 52(19):1527-1539, 2008. [doi:10.1016/j.jacc.2008.07.051](https://doi.org/10.1016/j.jacc.2008.07.051)
- <span id="ref-Shammas1988"></span>Shammas FV, Dickstein K. Clinical pharmacokinetics in heart failure. An updated review. Clin Pharmacokinet 15(2):94-113, 1988. [doi:10.2165/00003088-198815020-00002](https://doi.org/10.2165/00003088-198815020-00002)
- <span id="ref-Tan2018"></span>Tan ML, Yoshida K, Zhao P, Zhang L, Nolin TD, Piquette-Miller M, Galetin A, Huang SM. Effect of chronic kidney disease on nonrenal elimination pathways: a systematic assessment of CYP1A2, CYP2C8, CYP2C9, CYP2C19, and OATP. Clin Pharmacol Ther 103(5):854-867, 2018. [doi:10.1002/cpt.807](https://doi.org/10.1002/cpt.807)
- <span id="ref-Verbeeck2008"></span>Verbeeck RK. Pharmacokinetics and dosage adjustment in patients with hepatic dysfunction. Eur J Clin Pharmacol 64(12):1147-1161, 2008. [doi:10.1007/s00228-008-0553-z](https://doi.org/10.1007/s00228-008-0553-z)
- <span id="ref-Verbeeck2009"></span>Verbeeck RK, Musuamba FT. Pharmacokinetics and dosage adjustment in patients with renal dysfunction. Eur J Clin Pharmacol 65(8):757-773, 2009. [doi:10.1007/s00228-009-0678-8](https://doi.org/10.1007/s00228-009-0678-8)
- <span id="ref-WHO2000"></span>World Health Organization. Obesity: preventing and managing the global epidemic. Report of a WHO Consultation. WHO Technical Report Series 894, Geneva, 2000.
- <span id="ref-WHOAsian2004"></span>WHO Expert Consultation. Appropriate body-mass index for Asian populations and its implications for policy and intervention strategies. Lancet 363(9403):157-163, 2004. [doi:10.1016/S0140-6736(03)15268-3](https://doi.org/10.1016/S0140-6736(03)15268-3)
- <span id="ref-Wiesner2003"></span>Wiesner R, Edwards E, Freeman R, Harper A, Kim R, Kamath P, Kremers W, Lake J, Howard T, Merion RM, Wolfe RA, Krom R. Model for end-stage liver disease (MELD) and allocation of donor livers. Gastroenterology 124(1):91-96, 2003. [doi:10.1053/gast.2003.50016](https://doi.org/10.1053/gast.2003.50016)
- <span id="ref-Wilkinson1975"></span>Wilkinson GR, Shand DG. A physiological approach to hepatic drug clearance. Clin Pharmacol Ther 18(4):377-390, 1975. [doi:10.1002/cpt1975184377](https://doi.org/10.1002/cpt1975184377)
- <span id="ref-Yoshida2016"></span>Yoshida K, Sun B, Zhang L, Zhao P, Abernethy DR, Nolin TD, Rostami-Hodjegan A, Zineh I, Huang SM. Systematic and quantitative assessment of the effect of chronic kidney disease on CYP2D6 and CYP3A4/5. Clin Pharmacol Ther 100(1):75-87, 2016. [doi:10.1002/cpt.337](https://doi.org/10.1002/cpt.337)
<!-- end generated -->
