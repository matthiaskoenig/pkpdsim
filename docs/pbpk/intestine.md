# Intestine and routes

A route applies a dose of a substance: intravenous (bolus and infusion), oral tablet (with dissolution), oral solution, intramuscular, subcutaneous and rectal. A dose is a dose state in mg which an experiment sets or adds to (`PKPDSimulationExperiment.add_doses`, `"PODOSE_api + 10.0"`), applied by a first order process, the convention of the hand-written models; the dose in mg becomes mmole with the molecular weight `Mr_<s>`. The oral routes need an intestine with a lumen: the segmented intestine (the default) or the single lumen, which replace small and large intestine with one submodel `GU`, own the lumen, apply the oral doses and take the bile.

## Routes

| route | declared | dose state | parameters | applied into |
|---|---|---|---|---|
| intravenous | `Route.iv("api")` | `IVDOSE_api` (mg) | `ti_api` (half-life of the injection, s), `Ri_api` (infusion, mg/min) | the venous plasma `Cve_api` |
| oral tablet | `Route.oral("api", formulation)` | `PODOSE_api` (mg) | the dissolution of the formulation (`Ka_dis_api`, 1/hr, by default) | the stomach, undissolved |
| oral solution | `Route.oral_solution("api")` | `SOLDOSE_api` (mg) | `Ka_sol_api` (1/hr) | the stomach lumen, dissolved |
| intramuscular | `Route.im("api")` | `IMDOSE_api` (mg) | `ti_im_api` (s), `ka_im_api` (1/hr), `f_im_api` | the depot `Aim_api`, absorbed into the venous plasma `Cve_api` |
| subcutaneous | `Route.sc("api")` | `SCDOSE_api` (mg) | `ti_sc_api` (s), `ka_sc_api` (1/hr), `f_sc_api` | the depot `Asc_api`, absorbed into the venous plasma `Cve_api` |
| rectal | `Route.rectal("api")` | `RECDOSE_api` (mg) | `ti_rec_api` (s), `ka_rec_api` (1/hr), `f_rec_api`, `f_rec_portal_api` | the depot `Arec_api`, absorbed into the portal plasma `Cpo_api` and the venous plasma `Cve_api` |

The dose states and the parameters of the routes are parameters of the top model; the intestine reads the oral ones through its ports and owns their reactions and rate rules. Inhalation, sublingual and buccal, transdermal and vaginal routes and the glucose and insulin clamps are designed for (a route is a component with dose parameters) and implemented when the models which need them are migrated.

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

`PODOSE_<s>` is the undissolved tablet in the stomach: it dissolves into the stomach lumen and empties undissolved with the stomach into the duodenum. In the next segments the undissolved drug is a species of its own, `Asolid_<seg>_<s>`, which dissolves there and moves on undissolved into the feces. An enteric-coated tablet does not dissolve in the stomach (see [Dissolution](#dissolution)).

### Depot routes

The intramuscular, subcutaneous and rectal routes give their dose into a depot at the site of the administration, an amount of the top model in a compartment of the route (`Aim_<s>` in `Vim`, `Asc_<s>` in `Vsc`, `Arec_<s>` in `Vrec`, mmole; the compartments are compartments of amounts with the nominal volume 1 l, as the urine). The dose state drains into the depot with the half-life of the application `ti_<k>_<s>` (10 s by default, as the intravenous bolus, a choice of the library), so a dose which an experiment adds to the dose state (`add_doses`, `"IMDOSE_api + 10.0"`) adds to the depot. The depot is absorbed with first order: the fraction `f_<k>_<s>` (1 by default) into the venous plasma, the rest into a sink, `Aloss_im_<s>` and `Aloss_sc_<s>` (degraded at the site of the injection, amounts in the depot compartment) or the feces `Afeces_<s>` (the rectal dose which is not absorbed). Of the absorbed rectal dose the fraction `f_rec_portal_<s>` enters the portal plasma and passes the liver (the superior rectal vein drains into the portal vein), the rest enters the venous plasma (the middle and inferior rectal veins drain into the vena cava).

```equations title="intramuscular (subcutaneous: sc) and rectal (top model)"
J_im_<s>: -> Aim_<s> = ln(2) / ti_im_<s> * 60 s_per_min * IMDOSE_<s> / Mr_<s>
d(IMDOSE_<s>)/dt = -J_im_<s> * Mr_<s>
J_im_abs_<s>: Aim_<s> -> Cve_<s> = f_im_<s> * ka_im_<s> / 60 min_per_hr * Aim_<s>
J_im_loss_<s>: Aim_<s> -> Aloss_im_<s> = (1 - f_im_<s>) * ka_im_<s> / 60 min_per_hr * Aim_<s>
J_rec_<s>: -> Arec_<s> = ln(2) / ti_rec_<s> * 60 s_per_min * RECDOSE_<s> / Mr_<s>
d(RECDOSE_<s>)/dt = -J_rec_<s> * Mr_<s>
J_rec_abs_<s>: Arec_<s> -> Cve_<s> = f_rec_<s> * (1 - f_rec_portal_<s>) * ka_rec_<s> / 60 min_per_hr * Arec_<s>
J_rec_abs_po_<s>: Arec_<s> -> Cpo_<s> = f_rec_<s> * f_rec_portal_<s> * ka_rec_<s> / 60 min_per_hr * Arec_<s>
J_rec_loss_<s>: Arec_<s> -> Afeces_<s> = (1 - f_rec_<s>) * ka_rec_<s> / 60 min_per_hr * Arec_<s>
```

The defaults are those of the data (`database/tables/routes.toml`): the absorption rates `ka_<k>_<s>` are 1/hr, a choice of the library which is not verified, since the absorption from a depot depends on the substance, its formulation and the site (the fitted rates of the hand-written models: diazepam intramuscularly 1.3/hr, morphine intramuscularly 4.4/hr and subcutaneously 4.0/hr); the portal fraction of the rectal route is 0.45, from the systemic availability of lidocaine rectally (0.71, repeated 0.67) and orally (0.34, 0.27) in de Boer 1979 (`(1 - F_rectal)/(1 - F_oral)` with complete absorption, the abstract: slightly more than half of the rectal dose bypasses the liver), one drug, and the fraction depends on the site in the rectum. A model sets or fits its own values with the arguments of the route, numbers, `Value`s with a source or fitted values, which are the fields of `Depot`:

```python
# not executed
Route.im("mor", ka=fitted(FITS.parameter("ka_im_mor", unit=U.per_hr)))
Route.rectal("mor", ka=1.2, fraction=0.6, portal=0.5, application=600.0)
```

A depot dose is a dose of the equilibration (`DOSE_PREFIXES`): given before the initialization it moves to the start; the depot is a drug-free zero and the sink of an injected depot is set back at the start as the urine.

### Dissolution

The dissolution rate `kdis` (1/min) acts on the undissolved amount of every segment; the formulation of the tablet decides it:

- `FirstOrderDissolution(k, fed, release)`: `kdis = k` (`Ka_dis_<s>`, 1/hr), with two optional factors (see below); the default (2/hr, the default of the hand-written models, a choice of the library) when the route states no formulation;
- `WeibullDissolution(t50, shape, lag)`: the fraction dissolved `F = 1 - exp(-ln2*((t - lag)/t50)^shape)` after the lag time, applied as its hazard on the undissolved amount, with the time since the dose the mean age of the undissolved drug (`age_dis_<s>`); exact for one dose, doses which overlap share one clock (the amount-weighted mean age);
- `DoseDependentDissolution(k, dose50, fed)`: `kdis = k*D50/(D50 + D)` with the undissolved amount `D` (mg) of all segments, a dissolution which saturates with the dose (solubility-limited), with the optional factor of the fed state.

The first order dissolution has two optional factors, each creating a dimensionless parameter of the body only if it is set (without them the model has neither parameter and the rate is `Ka_dis_<s>`): `fed` creates `f_dis_fed_<s>`, the factor of the rate in the fed state (`FED = 1`; food slows the dissolution), and `release` creates `f_dis_release_<s>`, the factor of a modified-release formulation (the given value, 1 for the immediate-release tablet; the experiments of a prolonged-release tablet set it, e.g. 0.2). The dose-dependent dissolution has the factor of the fed state as well, `kdis = Ka_dis*(1 + FED*(f_dis_fed - 1))*D50/(D50 + D)`. The single lumen has no fed state, a `fed` factor there raises.

```equations title="first order dissolution with the factors of the formulation and the fed state"
J_po_<s>: -> Clumen_st_<s> = Ka_dis_<s> / 60 min_per_hr * f_dis_release_<s> * (1 + FED * (f_dis_fed_<s> - 1)) * PODOSE_<s> / Mr_<s>
```

```equations title="Weibull dissolution (t50, shape, lag of the tablet)"
kdis_<s> = piecewise(0, age_dis_<s> <= lag_dis_<s>, 0.693147180559945 * shape_dis_<s> / t50_dis_<s> * ((age_dis_<s> - lag_dis_<s>) / t50_dis_<s> + 1e-9)^(shape_dis_<s> - 1))
J_dis_du_<s>: Asolid_du_<s> -> Clumen_du_<s> = kdis_<s> * Asolid_du_<s>
```

An enteric-coated tablet does not dissolve in the stomach: `FirstOrderDissolution(..., enteric=True)` and `DoseDependentDissolution(..., enteric=True)` leave out the dissolution in the stomach, `J_po_<s>`, so the stomach has no dissolution rate of the tablet (it is 0, not a reaction with the rate 0). The undissolved tablet leaves the stomach with the gastric emptying of the condition, the transit time `(T_st + FED*(T_st_fed - T_st))*fimp_t_stomach` (`T_st_eff`), as the undissolved drug of every tablet does, and dissolves from the duodenum on with the dissolution of its formulation. The coated tablet in the stomach is not exposed to the intestinal fluid, so it does not saturate a dose-dependent dissolution: `D` is the undissolved amount of all segments without the coated tablet in the stomach, `Mr_<s>*(Asolid_du_<s> + ... + Asolid_rec_<s>)`, and the food effect on the emptying does not change the dissolution of the drug already emptied. The coat creates no parameter; without it the model is unchanged. The emptying is first order with the mean transit time of the stomach (of liquids, see the data of the segments); a tablet which does not disintegrate empties with the migrating motor complex, which the model does not have. The single lumen has no stomach, an enteric tablet there raises.

```equations title="enteric-coated tablet: the stomach only empties"
d(PODOSE_<s>)/dt = -J_trs_st_<s> * Mr_<s>
J_trs_st_<s>: -> Asolid_du_<s> = PODOSE_<s> / Mr_<s> / ((T_st + FED * (T_st_fed - T_st)) * fimp_t_stomach)
```

## Segmented intestine

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/intestine.svg"
</div>

The lumen segments follow the anatomical segmentation of Thelen 2011, as PK-Sim: stomach `st`, duodenum `du`, upper and lower jejunum `jej1`, `jej2`, upper and lower ileum `ile1`, `ile2`, caecum `cae`, colon ascendens `asc`, transversum `tra`, descendens `des`, sigmoid `sig` and rectum `rec`, then the feces; the tokens in the order of the transit are `ids.GUT_SEGMENTS`, which a model reads instead of a list of its own. A substance which reaches the lumen (oral route, biliary excretion, transport into the lumen) has the dissolved concentration `Clumen_<seg>_<s>` in every segment, which moves on with the transit `V*C/T` of the segment and from the rectum into the feces `Afeces_<s>`. The fed state `FED = 1` slows the gastric emptying to the stomach transit time of the fed state.

The absorption rate constant of the small intestine is `ka_<s>` (1/hr, 1 by default for a substance with an oral route, a choice of the library; 0 for any other); a segment absorbs with `ka_<s>*f_area_<seg>`, the mucosal surface per lumen volume relative to the small intestine (the stomach does not absorb); optional factors change it (see [Fraction absorbed and colon absorption](#fraction-absorbed-and-colon-absorption) and [Regional absorption](#regional-absorption)). Without chemistry of the intestine the substance is absorbed into the plasma of the intestinal wall:

```equations title="transit and absorption"
J_tr_st_<s>: Clumen_st_<s> -> Clumen_du_<s> = Vlumen_st * Clumen_st_<s> / ((T_st + FED * (T_st_fed - T_st)) * fimp_t_stomach)
J_tr_du_<s>: Clumen_du_<s> -> Clumen_jej1_<s> = Vlumen_du * Clumen_du_<s> / (T_du * fimp_t_si)
J_tr_rec_<s>: Clumen_rec_<s> -> Afeces_<s> = Vlumen_rec * Clumen_rec_<s> / T_rec
J_abs_du_<s>: Clumen_du_<s> -> Cpl_<s> = ka_<s> / 60 min_per_hr * f_area_du * Vlumen_du * Clumen_du_<s>
J_trs_du_<s>: Asolid_du_<s> -> Asolid_jej1_<s> = Asolid_du_<s> / (T_du * fimp_t_si)
J_trs_rec_<s>: Asolid_rec_<s> -> Afeces_<s> = Asolid_rec_<s> / T_rec
```

The factors `fimp_t_stomach` and `fimp_t_si` are the factors of chronic kidney disease on the gastric emptying time and the small intestinal transit time: they multiply the transit time, a factor above 1 slows the transit. They are 1 until their source is read (see [Impairment](impairment.md#limitations)). The female build multiplies the transit time of the stomach in the fed state by `fsex_t_stomach_fed` (1.2 from puberty to the menopause, the slower emptying of a meal of women; the fasted emptying of water does not differ) and the transit times of the segments of the large intestine by `fsex_t_colon` (1.333 of the adults, ICRP 89), see [Sex](physiology.md#sex). The transit time of the stomach of the condition, `(T_st + FED*(T_st_fed - T_st))*fimp_t_stomach`, is a constant of the top model as well, `T_st_eff` (`ids.T_STOMACH`, an output of the intestine evaluated after the changes before the simulation): a formula of the top model reads it, e.g. the rate of a conversion in the stomach which converts a given fraction of a solution, `kconv = 60*f/((1 - f)*T_st_eff)`, fasted and fed, instead of a copy of the fasted 22.79 min. The single lumen has no stomach and no `T_st_eff`. The data of the segments:

<!-- generated: intestine-segments (scripts/pbpk_tables.py), do not edit -->
*Generated from `database/tables/intestine_human.toml` by `scripts/pbpk_tables.py`; the notes of every number are in the database.*

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
    - `ModelChoice`: Choice of a library of pkpdsim (PBPK or machine perfusion), not a value of the literature; the notes of the value give the reasoning.
    - `Mudie2010`: Mudie DM, Amidon GL, Amidon GE. Physiological parameters for oral delivery and in vitro testing. Mol Pharm 7(5):1388-1405, 2010. [doi:10.1021/mp100149j](https://doi.org/10.1021/mp100149j)
    - `Paine1997`: Paine MF, Khalighi M, Fisher JM, Shen DD, Kunze KL, Marsh CL, Perkins JD, Thummel KE. Characterization of interintestinal and intraintestinal variations in human CYP3A-dependent metabolism. J Pharmacol Exp Ther 283(3):1552-1562, 1997.
    - `Schiller2005`: Schiller C, Fröhlich CP, Giessmann T, Siegmund W, Mönnikes H, Hosten N, Weitschies W. Intestinal fluid volumes and transit of dosage forms as assessed by magnetic resonance imaging. Aliment Pharmacol Ther 22(10):971-979, 2005. [doi:10.1111/j.1365-2036.2005.02683.x](https://doi.org/10.1111/j.1365-2036.2005.02683.x)
    - `Yu1999`: Yu LX, Amidon GL. A compartmental absorption and transit model for estimating oral drug absorption. Int J Pharm 186(2):119-125, 1999. [doi:10.1016/S0378-5173(99)00147-7](https://doi.org/10.1016/S0378-5173(99)00147-7)
<!-- end generated -->

The small intestinal transit time is split onto duodenum, jejunum and ileum by their lengths (Mudie 2010), the colonic transit onto its segments by the segmental transits of Metcalf 1987, and the fluid volume of a region onto its segments in proportion to their transit times. The pH is data only, no process reads it yet.

### Fraction absorbed and colon absorption

Two optional attributes of `SegmentedIntestine` change the absorption of the segments, both `None` by default. Without them the model has neither parameter nor reaction, and the absorption is the one above.

- `fabs`: the fraction of the absorbed flux of every segment which enters the body, the parameter `fabs_<s>` of the top model (as in the single lumen). The rest, `1 - fabs_<s>`, goes from the lumen of the segment into the feces (`J_feces_<seg>_<s>`). The lumen empties as fast as without it, so the onset of the absorption stays and less of the drug enters the body: for a substance with incomplete absorption and a fast onset, which a smaller `ka_<s>` cannot give (it delays the onset and moves the absorption down the intestine).
- `f_ka_colon`: the factor `f_ka_colon_<s>` of the absorption of the segments of the large intestine (caecum to rectum, the segments of the region `Organ.LARGE_INTESTINE`; the shorthand of the group `colon` of a [regional profile](#regional-absorption)), relative to `ka_<s>*f_area_<seg>`; the area factors stay. For a substance which the colon absorbs less than its area factor gives, or not at all (0): its slow transit otherwise adds a late absorption of what the small intestine left.

A set attribute is one value for every substance with an oral route (the substances whose `ka_<s>` has the default of the data), or a mapping from the substance id to a value (`ka={"cr": 2.0, "crtn": 0.5}`, `fabs={"cr": 0.6}`) for a model with several oral substances. In a mapping `ka` of a substance which is missing is the default of the data, `fabs_<s>` and `f_ka_colon_<s>` exist only for the substances of the mapping, and a key which is no substance of the definition (or, for `fabs` and `f_ka_colon`, a substance without an oral route) raises a `ValueError`; a `ka` for a substance without an oral route is allowed (its reabsorption after biliary excretion). The same holds for `ka` and `fabs` of the single lumen, where a substance missing from `fabs` has `fabs_<s> = 1`. A model fits the parameter of its substance by id (`fabs_api`, `f_ka_colon_api`). A substance without an oral route (a metabolite which reaches the lumen by the bile, `ka_<s> = 0`) gets neither parameter and no `J_feces_<seg>_<s>`, its absorption stays `ka_<s>*f_area_<seg>`. A value of 1 gives the model without the attribute numerically. `SingleLumen.fabs` differs: its `None` is `fabs_<s> = 1`, a parameter of every substance.

```equations title="absorption with fabs and f_ka_colon"
J_abs_du_<s>: Clumen_du_<s> -> Cpl_<s> = fabs_<s> * ka_<s> / 60 min_per_hr * f_area_du * Vlumen_du * Clumen_du_<s>
J_feces_du_<s>: Clumen_du_<s> -> Afeces_<s> = (1 - fabs_<s>) * ka_<s> / 60 min_per_hr * f_area_du * Vlumen_du * Clumen_du_<s>
J_abs_asc_<s>: Clumen_asc_<s> -> Cpl_<s> = fabs_<s> * ka_<s> / 60 min_per_hr * f_area_asc * f_ka_colon_<s> * Vlumen_asc * Clumen_asc_<s>
J_feces_asc_<s>: Clumen_asc_<s> -> Afeces_<s> = (1 - fabs_<s>) * ka_<s> / 60 min_per_hr * f_area_asc * f_ka_colon_<s> * Vlumen_asc * Clumen_asc_<s>
```

### Regional absorption

`SegmentedIntestine(regional=...)` gives a substance with an oral route a regional absorption profile, `RegionalAbsorption({"ile2": 0.5, "asc": 0.1})`: a factor `f_ka_<seg>_<s>` per segment (`du` to `rec`, the stomach does not absorb) which multiplies the absorption of the segment, `ka_<s>*f_area_<seg>*f_ka_<seg>_<s>`. The area factor holds the mucosal surface per lumen volume, the profile the permeability of the segment relative to the small intestine; a segment missing from the profile has no factor (1). One profile is that of every substance with an oral route, a mapping from the substance id gives each substance its own (`regional={"api": RegionalAbsorption(...)}`). A factor is a number or a `Value` (dimensionless, with its source or fitted). The group `colon` gives the segments of the large intestine (caecum to rectum) one shared factor, the parameter `f_ka_colon_<s>`: `f_ka_colon=x` is the shorthand of `RegionalAbsorption({"colon": x})`, both make the same model. A segment of the group with a factor of its own as well, or a substance with `f_ka_colon` and a factor of the large intestine in its profile, is refused (two factors of one absorption). A measured colonic absorption goes into the group (`RegionalAbsorption({"ile2": ..., "colon": Value(..., source="Byon2018")})`), a segment factor where the study resolves the segment. The notes of `f_ka_colon_<s>` say where its value comes from: of a `Value` with a source the regional absorption profile, the value of that source, not fitted (and its citation); of a number or a fitted value that a model fits it.

There is no generic default profile, the default is no factor (1, the absorption above): the absorption in the colon relative to the small intestine is a property of the substance. Of 42 drugs released in the colon in humans the highly permeable ones (BCS class I) had a relative colonic bioavailability above 70 % and those of low permeability (BCS III and IV) below 50 % (Tannergren 2009); a single profile would be wrong for one of the two groups. With the area factors of the large intestine alone a solution released in the caecum is absorbed to 70 % (in the colon ascendens 64 %) of the fraction absorbed of an oral solution at the default `ka` (1/hr), the range of a highly permeable drug (simulated with the data of this version, `test_the_colon_absorption_of_the_docs`). A substance whose regional absorption is measured takes its profile from the study, e.g. apixaban released by a capsule in the distal small intestine and the ascending colon (Byon 2018: the AUC 60 % and 84 % lower than after the oral solution).

A release study gives the exposure of a solution released in a segment relative to an oral solution, not the factors: they depend on the rest of the absorption of the model (`ka_<s>`, `fabs_<s>`, an efflux, the first pass), so they are derived with the model. `release_auc_ratios(r, substance, segments, dose)` simulates the oral solution (`SOLDOSE_<s>`) and the same dose dissolved in the lumen of each release segment (`GU__Clumen_<seg>_<s>`) and gives the AUC of the arm vein relative to the oral solution (14 days at one point per minute by default); `derive_regional_absorption(r, substance, targets, dose)` finds the factors which give the reported ratios, one factor per release segment (`regional_key`: the segment in the small intestine, the group `colon` in the large one), each the root of its ratio with the others fixed, in turns until no factor moves. A study in another condition than the reference (fed, a body weight, an impairment) passes it as `changes`, initial values by id, which the initial assignments follow; afterwards the condition is taken back and every value of the roadrunner is reset to its initial value (the factors of the model, not the derived ones). The model needs an oral solution of the substance and a first value of each factor in its profile. Because the factors depend on the fitted absorption, a model derives them on its adopted fit, writes them into its profile (a `Value` with the source of the study) and fits again from the adopted fit, until the factors derived on the fit are those of the model (apixaban: Byon 2018, the lower ileum 0.40 and the colon ascendens 0.16 give `f_ka_ile2_api` 0.149 and `f_ka_colon_api` 0.0567).

```python
# not executed
from pkdb_models.models.apixaban import MODEL
from pkpdsim.pbpk import derive_regional_absorption

r = MODEL.load_roadrunner()
r.integrator.absolute_tolerance = 1e-14
r.integrator.relative_tolerance = 1e-10
factors = derive_regional_absorption(r, "api", {"ile2": 0.40, "asc": 0.16}, dose=2.5)
```

```equations title="regional absorption (a factor of the lower ileum and of the group colon, with fabs)"
J_abs_ile2_<s>: Clumen_ile2_<s> -> Cpl_<s> = fabs_<s> * ka_<s> / 60 min_per_hr * f_area_ile2 * f_ka_ile2_<s> * Vlumen_ile2 * Clumen_ile2_<s>
J_abs_asc_<s>: Clumen_asc_<s> -> Cpl_<s> = fabs_<s> * ka_<s> / 60 min_per_hr * f_area_asc * f_ka_colon_<s> * Vlumen_asc * Clumen_asc_<s>
```

??? note "Sources"

    - `Tannergren2009`: Tannergren C, Bergendal A, Lennernäs H, Abrahamsson B. Toward an increased understanding of the barriers to colonic drug absorption in humans: implications for early controlled release candidate assessment. Mol Pharm 6(1):60-73, 2009. [doi:10.1021/mp800261a](https://doi.org/10.1021/mp800261a)
    - `Byon2018`: Byon W, Nepal S, Schuster AE, Shenker A, Frost CE. Regional gastrointestinal absorption of apixaban in healthy subjects. J Clin Pharmacol 58(7):965-971, 2018. Abstract read (Europe PMC, 2026-10-29); the full text is not open access. [doi:10.1002/jcph.1097](https://doi.org/10.1002/jcph.1097)

### Conversions in the lumen

`SegmentedIntestine(conversions=[LumenConversion(substrate, product, k, segments)])` converts a dissolved substance into another in the chosen lumen segments, first order with the rate constant `kconv_<s>_<m>` (1/hr, the unit of `ka_<s>`, a parameter of the top model), one mole of the product per mole of the substrate: a chemical conversion (creatine to creatinine in the acid of the stomach), the hydrolysis of a prodrug, a bacterial conversion in the colon (the deconjugation of a glucuronide). The undissolved drug of a tablet does not convert. The substrate reaches the lumen (an oral route, biliary excretion, a transport of the intestine or another conversion); the product gets the lumen species, the transit and the absorption of the intestine (its `ka_<m>` is 0 without an oral route unless `ka` gives it one). Without conversions the model has no rate constant and no reaction.

```equations title="conversion of <s> into <m> in the stomach"
J_conv_st_<s>_<m>: Clumen_st_<s> -> Clumen_st_<m> = kconv_<s>_<m> / 60 min_per_hr * Vlumen_st * Clumen_st_<s>
```

The stomach does not absorb, so the converted fraction `f` of an oral solution by a conversion in the stomach follows from the gastric transit time `T` (`T_st`, 22.79 min fasted, `T_st_fed` 63.48 min fed) and `k = kconv/60` (1/min): the dissolved drug leaves the stomach with `1/T` and converts with `k`,

```text
f = k*T / (1 + k*T)            kconv = 60 * f / ((1 - f) * T)
```

An immediate-release tablet with the first order dissolution `Ka_dis_<s>` (1/hr, without the factors `f_dis_fed_<s>` and `f_dis_release_<s>`) converts only what dissolves in the stomach, the fraction `g = kd*T/(1 + kd*T)` with `kd = Ka_dis/60` (1/min; not the rate `kdis_<s>` of a Weibull or dose-dependent dissolution, which is in 1/min already), so `f = g*k*T/(1 + k*T)`. Over several segments without absorption the drug survives each with `1/(1 + k*T_seg)`. A converted fraction of 0.0069 of a solution at the fasted transit is `kconv = 0.0183` 1/hr.

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

A substance with chemistry of the intestine (substrate, product, transported substance, inhibitor or inducer of a protein of a process of `Organ.SMALL_INTESTINE`) has a mucosa (enterocyte) species per segment, `Cmu_<seg>_<s>`: flow-limited by default (`mucosa_ps` makes it permeability-limited, see below) with the mucosa to plasma coefficient `K_mu_<s>` (the Kp of the small intestine, 1 without), perfused by arterial blood with the mucosa share of the blood flow of its region times the transit share of the segment in the region, and draining into the portal vein; the absorbed drug enters it. Its volume is carved out of the tissue of the wall, so wall and mucosa together have the capacity of the wall of a substance without mucosa. The mucosa is a tenth of the tissue of the intestine and receives 70 % of its blood flow, choices of the library (`database/tables/intestine_human.toml`). The chemistry of the small intestine is placed into the mucosa segments by the relative expression of its protein (CYP3A4 of Paine 1997) times the mucosa volume, the sum kept: the Vmax stays per litre of small intestine. A protein without expression data, and a process without protein, is placed by mucosa volume over the segments of the small intestine only, the organ of the process; the colon mucosa carries none of it, so the gut first pass of such a process equals that of the single lumen. An apical efflux moves the drug from the mucosa into the lumen of its segment, an apical influx back.

```equations title="mucosa of the duodenum (a gut metabolism of drg into drg_m)"
Vmu_du = f_mucosa_volume * Vsi * (1 - f_vas_si) * 0.0744792
Qmu_du = f_mucosa_flow * Qsi * 0.0744792
Qgu_wall = Qgu - Qmu_total
K_mu_<s> = Kp_<s>
J_mu_in_du_<s>: Car_<s> -> Cmu_du_<s> = Qmu_du * (1 - HCT) * Car_<s>
J_mu_out_du_<s>: Cmu_du_<s> -> Cpo_<s> = Qmu_du * (1 - HCT) * Cmu_du_<s> / K_mu_<s>
J_abs_du_<s>: Clumen_du_<s> -> Cmu_du_<s> = ka_<s> / 60 min_per_hr * f_area_du * Vlumen_du * Clumen_du_<s>
w_cyp3a4_du = 31 * Vmu_du / (31 * Vmu_du + 31 * Vmu_jej1 + 24 * Vmu_jej2 + 17 * Vmu_ile1 + 17 * Vmu_ile2)
DRGGU_du: Cmu_du_drg -> Cmu_du_drg_m = f_cyp3a4 * fimp_cyp3a4 * w_cyp3a4_du * DRGGU_Vmax * Vsi * (fu_drg / K_mu_drg * Cmu_du_drg) / (DRGGU_Km + fu_drg / K_mu_drg * Cmu_du_drg)
```

The number `0.0744792` is the transit share of the duodenum in the small intestine. A protein which a substance inactivates or induces has one turnover state for the intestine, driven by the unbound mucosa concentration weighted by the placement of the protein. ABCB1, ABCG2, CYP3A5 and every protein but CYP3A4 are placed by mucosa volume over the small intestine until their abundances along the intestine (Drozdzik 2014) are transcribed and verified.

### Permeability-limited mucosa

The flow-limited mucosa is perfused with the arterial blood of the mucosa, so the enterocytes see the total concentration of that blood and any gut extraction `E_G` of the absorbed drug extracts systemic drug with about `Qmu*E_G`, the mucosal blood flow being 70 % of the intestinal flow (about 29 l/hr of blood). A drug whose total clearance is far below that flow cannot have a gut first pass this way: with the gut availability 0.3 of a drug with fu 0.08 and a blood to plasma ratio of 142 (sirolimus, blood clearance 1.5 to 4 l/hr) the gut alone clears 20.6 l/hr of blood (`tests/pbpk/test_mucosa.py`).

`SegmentedIntestine(mucosa_ps={"rap": 0.5})` (and `SingleLumen(mucosa_ps=...)`) makes the mucosa of the substances it names permeability-limited: the enterocytes exchange with the mucosal blood through the permeability-surface product of their basolateral membrane on the unbound concentrations, in the convention of `ps` of an organ (`fu*(C_blood - Cmu/K_mu)`), while the absorbed drug still enters them from the lumen and the gut chemistry and the apical efflux still act on them. The mucosal blood is a fraction of the mucosa, so it is in quasi steady state, which leaves the hybrid flow of the Qgut model (Yang 2007, Gertz 2010), `Qgut = Qmu*carrier*fu*PS/(Qmu*carrier + fu*PS)` with the carrier of the mucosal blood flow (`1 - HCT` for a substance whose flows carry plasma, `BP_<s>` for erythrocytes in equilibrium, as in the flow-limited mucosa above; `Qmu*BP` is 69.8 l/min for sirolimus and the fence below shows the plasma form of its toy): the enterocytes see `Qgut` in place of the mucosal blood flow in the same two flows `J_mu_in_<seg>_<s>` and `J_mu_out_<seg>_<s>`, and the rest of the mucosal blood passes them from the artery into the portal vein (`J_mu_blood_<seg>_<s>`). A linear gut metabolism then has the gut availability of the Qgut model, `F_G = Qgut/(Qgut + fu*CLint)`, and the gut clears systemic drug with `Qgut*(1 - F_G)` of the arterial plasma: a small `PS` separates the two (sirolimus with a `CLperm` of 100 l/hr: the same gut availability 0.3 costs 0.039 l/hr of blood clearance instead of 20.6), and a large `PS` is the flow-limited mucosa, whose flows the model returns to (`Qgut` approaches `Qmu*carrier` and the bypass vanishes). Because the enterocytes are a compartment behind a barrier, the absorbed drug reaches the blood with the time constant `Vmu_<seg>*K_mu_<s>/Qgut_<seg>_<s>`, the same in every segment and about 0.47 minutes per unit of `K_mu_<s>` at a `CLperm` of 100 l/hr (0.0008 minutes with the flow-limited mucosa): a drug with a large mucosa to plasma coefficient is held in the mucosa for minutes, which is the physiology of the barrier and part of the option.

```equations title="permeability-limited mucosa of the duodenum (the exchange per litre of mucosa, a substance whose flows carry plasma)"
PS_mu_du_<s> = kps_mu_<s> * Vmu_du
Qgut_du_<s> = Qmu_du * (1 - HCT) * fu_<s> * PS_mu_du_<s> / (Qmu_du * (1 - HCT) + fu_<s> * PS_mu_du_<s>)
J_mu_in_du_<s>: Car_<s> -> Cmu_du_<s> = Qgut_du_<s> * Car_<s>
J_mu_out_du_<s>: Cmu_du_<s> -> Cpo_<s> = Qgut_du_<s> * Cmu_du_<s> / K_mu_<s>
J_mu_blood_du_<s>: Car_<s> -> Cpo_<s> = (Qmu_du * (1 - HCT) - Qgut_du_<s>) * Car_<s>
```

```equations title="the other form: the exchange of the mucosa of the small intestine"
Vmu_si = Vmu_du + Vmu_jej1 + Vmu_jej2 + Vmu_ile1 + Vmu_ile2
PS_mu_du_<s> = PS_mu_ref_<s> * Vmu_du / Vmu_si
```

The values are per litre of mucosa (1/min, `kps_mu_<s>`), so the exchange of a segment follows the mucosa of the individual; `mucosa_ps_per_volume=False` takes them as the exchange of the mucosa of the small intestine (l/min, `PS_mu_ref_<s>`, the `CLperm` of the Qgut model, which the literature gives for the small intestine and which the library normalises over the same organ as the chemistry of a protein without an expression profile). A segment has its share of it by mucosa volume, `PS_mu_ref_<s>*Vmu_<seg>/Vmu_si`, so the five segments of the small intestine sum to the given value and the colon adds its own share at the same exchange per litre of mucosa (its mucosa is 56.9 % of the small intestinal mucosa, so its segments sum to 0.569 of the value): the colon mucosa is permeability-limited as well, which keeps the drug absorbed there able to reach the blood (an exchange of zero would trap it in the colon enterocytes, which carry no chemistry). The colonic share is an extrapolation of the small intestinal density, and the one assumption of the option which the Qgut literature does not cover: `CLperm` is measured for the small intestine, and the colon has no villi, so its surface per litre of mucosa is smaller and the library overestimates its permeability (the error moves the colon towards the flow-limited mucosa, which is what a model gets without the option). The single lumen has one mucosa segment, the mucosa of the small intestine, so it gets the whole value and both intestines agree on the small intestine (it normalises over `Vmu_total`, which is that volume, and has no `Vmu_si`). Both forms mean the same exchange per litre of mucosa, `PS_mu_ref_<s>/Vmu_si`, and either is a number, a `Value` (fitted, in 1/min or l/min of its form) or a `Formula`; a number or a `Value` which is not positive is refused when the intestine is created (a negative exchange would run the flows of the mucosa backwards, zero would cut the enterocytes off from the blood), as is a `Value` in the unit of the other form; a `Formula` is not checked, as for the `ps` of an organ, since its value is only known when the model is built.

A substance which is not named keeps the flow-limited mucosa and nothing is created for it; a `mucosa_ps` of a substance without mucosa species is refused, naming the substances which have one. Use it for a drug with a gut first pass whose systemic clearance is below the mucosal blood flow, and the flow-limited mucosa (the default) for a drug which the mucosa extracts from the blood as fast as the flow delivers it. The blood cells of the mucosal blood pass the mucosa either way, and a basolateral transporter of the enterocytes is still not available (a basolateral transport of the small intestine is refused).

## Single lumen

`Variants(organs={Organ.SMALL_INTESTINE: "single_lumen"})` replaces the segments with one lumen, the fluid of the small intestine (105 ml, transit time 192 min): the first order absorption of most hand-written models, from the same interface. The lumen absorbs with `ka_<s>`, of which the fraction `fabs_<s>` enters the body and the rest the feces, and it keeps a physiological transit into the feces, so the fraction absorbed is `Fa = fabs*ka/(ka + 1/T)`: the fitted absorption constants of the hand-written models (which had no transit) do not carry over unchanged, a refit is expected anyway.

## Readouts

The lumen and the fluxes of the intestine are spread over its segments; a model which reads their sums names them in `PBPKDefinition.readouts`, and the intestine adds each as an output, as the kidney adds `Rurine_<s>`: the dissolved amount in the lumen `Alumen_<s>`, the summed fluxes of the absorption `Rabs_<s>`, of the efflux into the lumen `Reff_<s>` and of the loss to the feces `Rfeces_<s>`, the amount absorbed `Aabs_<s>` and the gross absorbed fraction of the oral dose `Fa_<s>` (absorbed from the lumen over the oral dose given, before the first pass, above the bioavailability), the same ids for the segmented intestine and the single lumen (see [Observers](observers.md#readouts-of-the-intestine)). A model without readouts has none of them.

## Bile and enterohepatic circulation

The bile `Cbile_<s>` (compartment `Vbile`) collects the apical efflux of the liver and drains with the first order `k_bile` (0.01 1/min, a choice of the library) into the duodenum, of a substance without lumen species into the feces:

```equations title="bile"
J_bile_<s>: Cbile_<s> -> Clumen_du_<s> = k_bile * Vbile * Cbile_<s>
```

With an absorption of the substance from the lumen (`ka_<s>`, which is 0 for a substance without an oral route unless the model sets or fits it) the excreted drug is reabsorbed: the enterohepatic circulation.

## Sinks

Every substance has the amounts in urine `Aurine_<s>` and feces `Afeces_<s>` (compartments of amounts, `Vurine` and `Vfeces`) and the bile `Cbile_<s>`. With the standard organs as intestine, which have no lumen, the bile drains into the feces. The breath is a sink of a later project.

## Limitations

- A conversion in the lumen is first order (`LumenConversion`): a saturable or enzyme-limited conversion in the lumen is not available, a custom reaction of the intestine is refused.
- The single lumen has no segments: `SingleLumen(regional=...)` and `SingleLumen(conversions=...)` are refused with a message, a regional profile and the conversions in the lumen need the segmented intestine.
- A regional profile is given to a substance with an oral route (as `fabs` and `f_ka_colon`): the product of a conversion without an oral route absorbs with `ka_<m>*f_area_<seg>` only.
- Only CYP3A4 has an expression profile along the intestine; every other protein is placed by mucosa volume over the small intestine, so a process expressed in the colon needs expression data (or a hook for the large intestine). The pH of the segments is data no process reads.
- The mucosa has no basolateral transporter: a basolateral transport of the small intestine is refused, the enterocytes exchange with the mucosal blood passively (flow-limited, or through the permeability of `mucosa_ps`). A permeability-limited mucosa is a passive permeability; a drug taken up into the enterocytes by a basolateral transporter needs one.
- The gastric emptying and small intestinal transit factors of chronic kidney disease are 1 (see [Impairment](impairment.md#limitations)).
