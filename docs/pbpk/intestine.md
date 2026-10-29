# Intestine and routes

A route applies a dose of a substance: intravenous (bolus and infusion), oral tablet (with dissolution) and oral solution. A dose is a dose state in mg which an experiment sets or adds to (`PKPDSimulationExperiment.add_doses`, `"PODOSE_api + 10.0"`), applied by a first order process, the convention of the hand-written models; the dose in mg becomes mmole with the molecular weight `Mr_<s>`. The oral routes need an intestine with a lumen: the segmented intestine (the default) or the single lumen, which replace small and large intestine with one submodel `GU`, own the lumen, apply the oral doses and take the bile.

## Routes

| route | declared | dose state | parameters | applied into |
|---|---|---|---|---|
| intravenous | `Route.iv("api")` | `IVDOSE_api` (mg) | `ti_api` (half-life of the injection, s), `Ri_api` (infusion, mg/min) | the venous plasma `Cve_api` |
| oral tablet | `Route.oral("api", formulation)` | `PODOSE_api` (mg) | the dissolution of the formulation (`Ka_dis_api`, 1/hr, by default) | the stomach, undissolved |
| oral solution | `Route.oral_solution("api")` | `SOLDOSE_api` (mg) | `Ka_sol_api` (1/hr) | the stomach lumen, dissolved |

The dose states and the parameters of the routes are parameters of the top model; the intestine reads the oral ones through its ports and owns their reactions and rate rules. Subcutaneous, intramuscular, inhalation, sublingual and buccal, rectal, transdermal and vaginal routes and the glucose and insulin clamps are designed for (a route is a component with dose parameters) and implemented when the models which need them are migrated.

```equations title="intravenous (top model)"
J_iv_<s>: -> Cve_<s> = Kiv_<s> * IVDOSE_<s> / Mr_<s>
Kiv_<s> = ln(2) / ti_<s> * 60 s_per_min
d(IVDOSE_<s>)/dt = -J_iv_<s> * Mr_<s> + Ri_<s>
```

```equations title="oral solution and tablet (intestine)"
J_sol_<s>: -> Clumen_st_<s> = Ka_sol_<s> / 60 min_per_hr * SOLDOSE_<s> / Mr_<s>
d(SOLDOSE_<s>)/dt = -J_sol_<s> * Mr_<s>
J_po_<s>: -> Clumen_st_<s> = Ka_dis_<s> / 60 min_per_hr * PODOSE_<s> / Mr_<s>
J_trs_st_<s>: -> Asolid_du_<s> = PODOSE_<s> / Mr_<s> / ((T_st + FED * (T_st_fed - T_st)) * fimp_t_stomach)
d(PODOSE_<s>)/dt = -(J_po_<s> + J_trs_st_<s>) * Mr_<s>
```

`PODOSE_<s>` is the undissolved tablet in the stomach: it dissolves into the stomach lumen and empties undissolved with the stomach into the duodenum. In the next segments the undissolved drug is a species of its own, `Asolid_<seg>_<s>`, which dissolves there and moves on undissolved into the feces.

### Dissolution

The dissolution rate `kdis` (1/min) acts on the undissolved amount of every segment; the formulation of the tablet decides it:

- `FirstOrderDissolution(k, fed, release)`: `kdis = k` (`Ka_dis_<s>`, 1/hr), with two optional factors (see below); the default (2/hr, the default of the hand-written models, a choice of the library) when the route states no formulation;
- `WeibullDissolution(t50, shape, lag)`: the fraction dissolved `F = 1 - exp(-ln2*((t - lag)/t50)^shape)` after the lag time, applied as its hazard on the undissolved amount, with the time since the dose the mean age of the undissolved drug (`age_dis_<s>`); exact for one dose, doses which overlap share one clock (the amount-weighted mean age);
- `DoseDependentDissolution(k, dose50)`: `kdis = k*D50/(D50 + D)` with the undissolved amount `D` (mg) of all segments, a dissolution which saturates with the dose (solubility-limited).

The first order dissolution has two optional factors, each creating a dimensionless parameter of the body only if it is set (without them the model has neither parameter and the rate is `Ka_dis_<s>`): `fed` creates `f_dis_fed_<s>`, the factor of the rate in the fed state (`FED = 1`; food slows the dissolution), and `release` creates `f_dis_release_<s>`, the factor of a modified-release formulation (the given value, 1 for the immediate-release tablet; the experiments of a prolonged-release tablet set it, e.g. 0.2). The single lumen has no fed state, a `fed` factor there raises.

```equations title="first order dissolution with the factors of the formulation and the fed state"
J_po_<s>: -> Clumen_st_<s> = Ka_dis_<s> / 60 min_per_hr * f_dis_release_<s> * (1 + FED * (f_dis_fed_<s> - 1)) * PODOSE_<s> / Mr_<s>
```

```equations title="Weibull dissolution (t50, shape, lag of the tablet)"
kdis_<s> = piecewise(0, age_dis_<s> <= lag_dis_<s>, 0.693147180559945 * shape_dis_<s> / t50_dis_<s> * ((age_dis_<s> - lag_dis_<s>) / t50_dis_<s> + 1e-9)^(shape_dis_<s> - 1))
J_dis_du_<s>: Asolid_du_<s> -> Clumen_du_<s> = kdis_<s> * Asolid_du_<s>
```

## Segmented intestine

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/intestine.svg"
</div>

The lumen segments follow the anatomical segmentation of Thelen 2011, as PK-Sim: stomach `st`, duodenum `du`, upper and lower jejunum `jej1`, `jej2`, upper and lower ileum `ile1`, `ile2`, caecum `cae`, colon ascendens `asc`, transversum `tra`, descendens `des`, sigmoid `sig` and rectum `rec`, then the feces; the tokens in the order of the transit are `ids.GUT_SEGMENTS`, which a model reads instead of a list of its own. A substance which reaches the lumen (oral route, biliary excretion, transport into the lumen) has the dissolved concentration `Clumen_<seg>_<s>` in every segment, which moves on with the transit `V*C/T` of the segment and from the rectum into the feces `Afeces_<s>`. The fed state `FED = 1` slows the gastric emptying to the stomach transit time of the fed state.

The absorption rate constant of the small intestine is `ka_<s>` (1/hr, 1 by default for a substance with an oral route, a choice of the library; 0 for any other); a segment absorbs with `ka_<s>*f_area_<seg>`, the mucosal surface per lumen volume relative to the small intestine (the stomach does not absorb); two optional factors change it (see [Fraction absorbed and colon absorption](#fraction-absorbed-and-colon-absorption)). Without chemistry of the intestine the substance is absorbed into the plasma of the intestinal wall:

```equations title="transit and absorption"
J_tr_st_<s>: Clumen_st_<s> -> Clumen_du_<s> = Vlumen_st * Clumen_st_<s> / ((T_st + FED * (T_st_fed - T_st)) * fimp_t_stomach)
J_tr_du_<s>: Clumen_du_<s> -> Clumen_jej1_<s> = Vlumen_du * Clumen_du_<s> / (T_du * fimp_t_si)
J_tr_rec_<s>: Clumen_rec_<s> -> Afeces_<s> = Vlumen_rec * Clumen_rec_<s> / T_rec
J_abs_du_<s>: Clumen_du_<s> -> Cpl_<s> = ka_<s> / 60 min_per_hr * f_area_du * Vlumen_du * Clumen_du_<s>
J_trs_du_<s>: Asolid_du_<s> -> Asolid_jej1_<s> = Asolid_du_<s> / (T_du * fimp_t_si)
J_trs_rec_<s>: Asolid_rec_<s> -> Afeces_<s> = Asolid_rec_<s> / T_rec
```

The factors `fimp_t_stomach` and `fimp_t_si` are the factors of chronic kidney disease on the gastric emptying time and the small intestinal transit time: they multiply the transit time, a factor above 1 slows the transit. They are 1 until their source is read (see [Impairment](impairment.md#limitations)). The data of the segments:

<!-- generated: intestine-segments (scripts/pbpk_tables.py), do not edit -->
*Generated from `data/intestine_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the data file.*

| segment | volume [l] | transit [min] | transit fed [min] | pH [-] | area factor [-] | CYP3A4 [-] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| stomach | 0.028 | 22.79 | 63.48 | 1.7 | 0 |  |
| duodenum | 0.007819 | 14.3 |  | 6.2 | 1 | 31 |
| upper jejunum | 0.01955 | 35.74 |  | 6.6 | 1 | 31 |
| lower jejunum | 0.01955 | 35.74 |  | 6.6 | 1 | 24 |
| upper ileum | 0.02904 | 53.11 |  | 7.2 | 1 | 17 |
| lower ileum | 0.02904 | 53.11 |  | 7.5 | 1 | 17 |
| caecum | 0.002093 | 339 |  | 6.4 | 0.0362 | 0 |
| colon ascendens | 0.002093 | 339 |  | 6.6 | 0.0362 | 0 |
| colon transversum | 0.002111 | 342 |  | 6.8 | 0.0362 | 0 |
| colon descendens | 0.002111 | 342 |  | 7 | 0.0362 | 0 |
| colon sigmoid | 0.002296 | 372 |  | 7 | 0.0362 | 0 |
| rectum | 0.002296 | 372 |  | 7 | 0.0362 | 0 |

Sources of the columns: volume: `Mudie2010`, `Schiller2005`; transit: `Metcalf1987`, `Mudie2010`; transit fed: `Mudie2010`; pH: `Evans1988`, `Mudie2010`; area factor: `Helander2014`, `Mudie2010`, `Yu1999`; CYP3A4: `ModelChoice`, `Paine1997`.

??? note "Sources"

    - `Evans1988`: Evans DF, Pye G, Bramley R, Clark AG, Dyson TJ, Hardcastle JD. Measurement of gastrointestinal pH profiles in normal ambulant human subjects. Gut 29(8):1035-1041, 1988. [doi:10.1136/gut.29.8.1035](https://doi.org/10.1136/gut.29.8.1035)
    - `Helander2014`: Helander HF, Fändriks L. Surface area of the digestive tract - revisited. Scand J Gastroenterol 49(6):681-689, 2014. [doi:10.3109/00365521.2014.898326](https://doi.org/10.3109/00365521.2014.898326)
    - `Metcalf1987`: Metcalf AM, Phillips SF, Zinsmeister AR, MacCarty RL, Beart RW, Wolff BG. Simplified assessment of segmental colonic transit. Gastroenterology 92(1):40-47, 1987. [doi:10.1016/0016-5085(87)90837-7](https://doi.org/10.1016/0016-5085(87)90837-7)
    - `ModelChoice`: Choice of the pkpdsim PBPK library, not a value of the literature; the notes of the value give the reasoning.
    - `Mudie2010`: Mudie DM, Amidon GL, Amidon GE. Physiological parameters for oral delivery and in vitro testing. Mol Pharm 7(5):1388-1405, 2010. [doi:10.1021/mp100149j](https://doi.org/10.1021/mp100149j)
    - `Paine1997`: Paine MF, Khalighi M, Fisher JM, Shen DD, Kunze KL, Marsh CL, Perkins JD, Thummel KE. Characterization of interintestinal and intraintestinal variations in human CYP3A-dependent metabolism. J Pharmacol Exp Ther 283(3):1552-1562, 1997.
    - `Schiller2005`: Schiller C, Fröhlich CP, Giessmann T, Siegmund W, Mönnikes H, Hosten N, Weitschies W. Intestinal fluid volumes and transit of dosage forms as assessed by magnetic resonance imaging. Aliment Pharmacol Ther 22(10):971-979, 2005. [doi:10.1111/j.1365-2036.2005.02683.x](https://doi.org/10.1111/j.1365-2036.2005.02683.x)
    - `Yu1999`: Yu LX, Amidon GL. A compartmental absorption and transit model for estimating oral drug absorption. Int J Pharm 186(2):119-125, 1999. [doi:10.1016/S0378-5173(99)00147-7](https://doi.org/10.1016/S0378-5173(99)00147-7)
<!-- end generated -->

The small intestinal transit time is split onto duodenum, jejunum and ileum by their lengths (Mudie 2010), the colonic transit onto its segments by the segmental transits of Metcalf 1987, and the fluid volume of a region onto its segments in proportion to their transit times. The pH is data only, no process reads it yet.

### Fraction absorbed and colon absorption

Two optional attributes of `SegmentedIntestine` change the absorption of the segments, both `None` by default. Without them the model has neither parameter nor reaction, and the absorption is the one above.

- `fabs`: the fraction of the absorbed flux of every segment which enters the body, the parameter `fabs_<s>` of the top model (as in the single lumen). The rest, `1 - fabs_<s>`, goes from the lumen of the segment into the feces (`J_feces_<seg>_<s>`). The lumen empties as fast as without it, so the onset of the absorption stays and less of the drug enters the body: for a substance with incomplete absorption and a fast onset, which a smaller `ka_<s>` cannot give (it delays the onset and moves the absorption down the intestine).
- `f_ka_colon`: the factor `f_ka_colon_<s>` of the absorption of the segments of the large intestine (caecum to rectum, the segments of the region `Organ.LARGE_INTESTINE`), relative to `ka_<s>*f_area_<seg>`; the area factors stay. For a substance which the colon absorbs less than its area factor gives, or not at all (0): its slow transit otherwise adds a late absorption of what the small intestine left.

A set attribute is one value for every substance with an oral route (the substances whose `ka_<s>` has the default of the data), or a mapping from the substance id to a value (`ka={"cr": 2.0, "crtn": 0.5}`, `fabs={"cr": 0.6}`) for a model with several oral substances. In a mapping `ka` of a substance which is missing is the default of the data, `fabs_<s>` and `f_ka_colon_<s>` exist only for the substances of the mapping, and a key which is no substance of the definition (or, for `fabs` and `f_ka_colon`, a substance without an oral route) raises a `ValueError`; a `ka` for a substance without an oral route is allowed (its reabsorption after biliary excretion). The same holds for `ka` and `fabs` of the single lumen, where a substance missing from `fabs` has `fabs_<s> = 1`. A model fits the parameter of its substance by id (`fabs_api`, `f_ka_colon_api`). A substance without an oral route (a metabolite which reaches the lumen by the bile, `ka_<s> = 0`) gets neither parameter and no `J_feces_<seg>_<s>`, its absorption stays `ka_<s>*f_area_<seg>`. A value of 1 gives the model without the attribute numerically. `SingleLumen.fabs` differs: its `None` is `fabs_<s> = 1`, a parameter of every substance.

```equations title="absorption with fabs and f_ka_colon"
J_abs_du_<s>: Clumen_du_<s> -> Cpl_<s> = fabs_<s> * ka_<s> / 60 min_per_hr * f_area_du * Vlumen_du * Clumen_du_<s>
J_feces_du_<s>: Clumen_du_<s> -> Afeces_<s> = (1 - fabs_<s>) * ka_<s> / 60 min_per_hr * f_area_du * Vlumen_du * Clumen_du_<s>
J_abs_asc_<s>: Clumen_asc_<s> -> Cpl_<s> = fabs_<s> * ka_<s> / 60 min_per_hr * f_area_asc * f_ka_colon_<s> * Vlumen_asc * Clumen_asc_<s>
J_feces_asc_<s>: Clumen_asc_<s> -> Afeces_<s> = (1 - fabs_<s>) * ka_<s> / 60 min_per_hr * f_area_asc * f_ka_colon_<s> * Vlumen_asc * Clumen_asc_<s>
```

### Intestinal wall

The wall is the standard organ of the volume `Vgu = Vsi + Vla`, with the volume-weighted vascular and interstitial fractions and the distribution of the small intestine (its Kp; a Kp override of the large intestine is refused). It receives the arterial inflow of both organs, `Qgu = Qsi + Qla`, and drains into the portal vein.

```equations title="intestinal wall without mucosa"
Vgu = Vsi + Vla
f_vas_gu = (Vsi * f_vas_si + Vla * f_vas_la) / Vgu
Qgu = Qsi + Qla
J_in_ar_pl_<s>: Car_<s> -> Cpl_<s> = Qgu * (1 - HCT) * Car_<s>
J_out_pl_<s>: Cpl_<s> -> Cpo_<s> = Qgu * (1 - HCT) * Cout_<s>
```

### Mucosa

A substance with chemistry of the intestine (substrate, product, transported substance, inhibitor or inducer of a protein of a process of `Organ.SMALL_INTESTINE`) has a mucosa (enterocyte) species per segment, `Cmu_<seg>_<s>`: flow-limited with the mucosa to plasma coefficient `K_mu_<s>` (the Kp of the small intestine, 1 without), perfused by arterial blood with the mucosa share of the blood flow of its region times the transit share of the segment in the region, and draining into the portal vein; the absorbed drug enters it. Its volume is carved out of the tissue of the wall, so wall and mucosa together have the capacity of the wall of a substance without mucosa. The mucosa is a tenth of the tissue of the intestine and receives 70 % of its blood flow, choices of the library (`data/intestine_human.toml`). The chemistry of the small intestine is placed into the mucosa segments by the relative expression of its protein (CYP3A4 of Paine 1997) times the mucosa volume, the sum kept: the Vmax stays per litre of small intestine. A protein without expression data, and a process without protein, is placed by mucosa volume over the segments of the small intestine only, the organ of the process; the colon mucosa carries none of it, so the gut first pass of such a process equals that of the single lumen. An apical efflux moves the drug from the mucosa into the lumen of its segment, an apical influx back.

```equations title="mucosa of the duodenum (a gut metabolism of drg into drg_m)"
Vmu_du = f_mucosa_volume * Vsi * (1 - f_vas_si) * 0.0744792
Qmu_du = f_mucosa_flow * Qsi * 0.0744792
Qgu_wall = Qgu - Qmu_total
K_mu_<s> = Kp_<s>
J_mu_in_du_<s>: Car_<s> -> Cmu_du_<s> = Qmu_du * (1 - HCT) * Car_<s>
J_mu_out_du_<s>: Cmu_du_<s> -> Cpo_<s> = Qmu_du * (1 - HCT) * Cmu_du_<s> / K_mu_<s>
J_abs_du_<s>: Clumen_du_<s> -> Cmu_du_<s> = ka_<s> / 60 min_per_hr * f_area_du * Vlumen_du * Clumen_du_<s>
w_cyp3a4_du = 31 * Vmu_du / (31 * Vmu_du + 31 * Vmu_jej1 + 24 * Vmu_jej2 + 17 * Vmu_ile1 + 17 * Vmu_ile2)
DRGGU_du: Cmu_du_drg -> Cmu_du_drg_m = f_cyp3a4 * w_cyp3a4_du * DRGGU_Vmax * Vsi * (fu_drg / K_mu_drg * Cmu_du_drg) / (DRGGU_Km + fu_drg / K_mu_drg * Cmu_du_drg)
```

The number `0.0744792` is the transit share of the duodenum in the small intestine. A protein which a substance inactivates or induces has one turnover state for the intestine, driven by the unbound mucosa concentration weighted by the placement of the protein. ABCB1, ABCG2 and every protein but CYP3A4 are placed by mucosa volume over the small intestine until their abundances along the intestine (Drozdzik 2014) are transcribed and verified.

## Single lumen

`Variants(organs={Organ.SMALL_INTESTINE: "single_lumen"})` replaces the segments with one lumen, the fluid of the small intestine (105 ml, transit time 192 min): the first order absorption of most hand-written models, from the same interface. The lumen absorbs with `ka_<s>`, of which the fraction `fabs_<s>` enters the body and the rest the feces, and it keeps a physiological transit into the feces, so the fraction absorbed is `Fa = fabs*ka/(ka + 1/T)`: the fitted absorption constants of the hand-written models (which had no transit) do not carry over unchanged, a refit is expected anyway.

## Readouts

The lumen and the fluxes of the intestine are spread over its segments; a model which reads their sums names them in `PBPKDefinition.readouts`, and the intestine adds each as an output, as the kidney adds `Rurine_<s>`: the dissolved amount in the lumen `Alumen_<s>` and the summed fluxes of the absorption `Rabs_<s>`, of the efflux into the lumen `Reff_<s>` and of the loss to the feces `Rfeces_<s>`, the same ids for the segmented intestine and the single lumen (see [Observers](observers.md#readouts-of-the-intestine)). A model without readouts has none of them.

## Bile and enterohepatic circulation

The bile `Cbile_<s>` (compartment `Vbile`) collects the apical efflux of the liver and drains with the first order `k_bile` (0.01 1/min, a choice of the library) into the duodenum, of a substance without lumen species into the feces:

```equations title="bile"
J_bile_<s>: Cbile_<s> -> Clumen_du_<s> = k_bile * Vbile * Cbile_<s>
```

With an absorption of the substance from the lumen (`ka_<s>`, which is 0 for a substance without an oral route unless the model sets or fits it) the excreted drug is reabsorbed: the enterohepatic circulation.

## Sinks

Every substance has the amounts in urine `Aurine_<s>` and feces `Afeces_<s>` (compartments of amounts, `Vurine` and `Vfeces`) and the bile `Cbile_<s>`. With the standard organs as intestine, which have no lumen, the bile drains into the feces. The breath is a sink of a later project.

## Limitations

- The deconjugation in the colon lumen (glucuronides back to the parent) is not implemented: a custom reaction of the intestine is refused, a hook for the chemistry of the lumen is open.
- Only CYP3A4 has an expression profile along the intestine; every other protein is placed by mucosa volume over the small intestine, so a process expressed in the colon needs expression data (or a hook for the large intestine). The pH of the segments is data no process reads.
- The gastric emptying and small intestinal transit factors of chronic kidney disease are 1 (see [Impairment](impairment.md#limitations)).
