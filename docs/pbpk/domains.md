# Domains

A model of the library states what every one of its quantities is (an SBO term, a PBPKO term, an UBERON term, see [Annotation](annotation.md)) and, for a fitted or a measured value, what is known about it (see [Uncertainty](../uncertainty.md)). It also states what values a quantity may take at all: its **domain**. A fraction absorbed is a number in [0, 1], an organ volume or a concentration is not negative, a severity of a disease is in [0, 3] or [0, 4], a switch is 0 or 1. The domain is written into the SBML next to the quantity, so that a tool which samples, scans, fits or checks the model asks the model instead of carrying a copy of that knowledge of its own.

## What a domain is not

A domain is the set of values the model accepts. It is not the uncertainty of a value, which is what is known about the value, and it is not the posterior of a fit, which is what the data say about it. The three are different statements and each has its own element in the SBML: a `distrib` uncertainty named `source ...` for the source of a value, the uncertainty of a fit with its distribution, and the uncertainty named `domain`. A value of 0.94 for a fraction absorbed may be known to within, say, 0.02, and the data of a fit may support values from 0.85 to 0.99, and the domain of the quantity is still [0, 1], because a fraction absorbed is a fraction.

Conflating them is the root of [pkpdsim#207](https://github.com/matthiaskoenig/pkpdsim/issues/207). `library_parameters` varied every parameter of a sensitivity analysis by the same relative span, which knows nothing of what a value means: with the default span of 0.15 a fraction absorbed of 0.94 received the bounds 0.799 and 1.081 (0.94 times 0.85 and 0.94 times 1.15), and the scan asked for more than the whole dose absorbed. The clip belongs where the span is built, and the knowledge it needs belongs in the model.

A second example is `fraction_range` of the chemistry, which reads the support of the draws of a fitted fraction to decide whether the product fractions of an enzyme can sum above 1. It intersects the support with the domain of a fraction, the unit interval, and it treats the two ends differently. An infinite end of the support comes from the family of the distribution (a normal has no end) and states no intent, so it is cut at the unit interval in silence. A finite end outside the unit interval is a statement of the author of the fit which contradicts the domain, and `check_products` refuses it. The domain is what the support is checked against, the support is never taken for the domain.

## The shape in the SBML

The domain is a `distrib` uncertainty named `domain` on the parameter, the compartment or the species, holding one span of type `range`. It sits next to the uncertainty of the source of the value, which it does not replace. The element below is a sketch of the shape, written by hand rather than taken from a model, and `Example2020` stands for the citation of whichever source a real value carries:

```xml
<parameter id="fabs_lvn" value="0.94" units="dimensionless" constant="true">
  <distrib:listOfUncertainties>
    <distrib:uncertainty distrib:name="source Example2020">
      <distrib:uncertParameter distrib:value="0.02" distrib:type="standardDeviation"/>
    </distrib:uncertainty>
    <distrib:uncertainty distrib:name="domain">
      <distrib:uncertSpan distrib:type="range" distrib:valueLower="0" distrib:valueUpper="1"/>
    </distrib:uncertainty>
  </distrib:listOfUncertainties>
</parameter>
```

The `range` span is meant for the observed range of a measurement, and `distrib` has no element for the domain of validity, so the distinction is carried by the name of the uncertainty. A reader finds the domain by its name and its span type and ignores every other uncertainty. A bound which is not set, or `INF`, is the infinity on its side.

The values inside the interval are a second statement. They are continuous by default, which writes nothing further, so a document written before the value sets keeps its meaning. A domain of the integers (an index, a count, a class) or of the two truth values (a flag) adds one child to the uncertainty, an `externalParameter` whose `definitionURL` is the XSD datatype:

```xml
<distrib:uncertainty distrib:name="domain">
  <distrib:uncertSpan distrib:type="range" distrib:valueLower="0" distrib:valueUpper="1"/>
  <distrib:uncertParameter distrib:value="1" distrib:type="externalParameter" distrib:definitionURL="http://www.w3.org/2001/XMLSchema#boolean"/>
</distrib:uncertainty>
```

The child carries `value="1"` only for the tools: libsbml validates a bare `definitionURL` clean, but libAntimony rejects it ([sys-bio/antimony#190](https://github.com/sys-bio/antimony/issues/190)), and the archive of a model converts its SBML to antimony. Two kinds are boolean, `ckd` (`CKD_hepatic`) and `lobe_present`. `fed` and `lobe_patency` are continuous on purpose: they enter the model as an interpolation between two states and not as a test, so a value between 0 and 1 means something. Nothing enforces the value set at run time, it is a contract for the tools which read it.

## Silence is not a default

An element without a domain means the model says nothing about that quantity. It never means that a default applies: `domain_of` gives `None` and `domains` leaves the element out. The library states the domain of every parameter, compartment and species it writes, including the unbounded ones (`REAL` is written as an interval from `-INF` to `INF`), so that silence is unambiguous. A reaction holds no number and states none.

## The domain of a kind

The domain of a quantity comes from its kind, the row of `KINDS` of `pkpdsim.pbpk.annotation` which already says what the element is: a kind states its domain once and every element of the kind gets it, and the id of a fraction appears nowhere. `annotate` writes the domain of the kind unless the element already carries one, so an element which states its own domain wins over its kind. Four named domains are in use:

| name | interval | what it is |
|---|---|---|
| `FRACTION` | [0, 1] | a fraction of a dose, a flux, a volume, a protein, an unbound fraction |
| `NONNEGATIVE` | [0, inf) | the ordinary physiological or biochemical quantity: a volume, a flow, a rate constant, a concentration, a scaling factor |
| `REAL` | (-inf, inf) | a quantity which may be negative: a difference, an offset |
| `integers(lower, upper)`, `BOOLEAN` | stated | a class or a count, a switch |

`domains(model_path)` reads only the top-level model, so hand it the flat file of a model: the elements of a submodel of a composed document are not in the result, and their absence says nothing about their domain.

A positive and a non-negative quantity have the same domain, because `distrib` cannot say whether a bound is attained; the strict checks of the definition (`check_positive`) keep their own wording. Every kind of the library, the number of kinds by domain first:

<!-- generated: domains-summary (scripts/pbpk_tables.py), do not edit -->
*Generated from `KINDS` of `pkpdsim.pbpk.annotation` by `scripts/pbpk_tables.py`; a domain is classified by hand, no database entry stands behind it.*

| domain | kinds |
| --- | ---: |
| `(-inf, inf)` | 5 |
| `[0, 1] boolean` | 2 |
| `[0, 1]` | 25 |
| `[0, 24]` | 1 |
| `[0, 3]` | 1 |
| `[0, 4]` | 1 |
| `[0, inf)` | 172 |
| none | 36 |
| all kinds | 243 |
<!-- end generated -->

The full list follows, one row per kind. Most kinds are `[0, inf)`, because most quantities of a physiological model are a volume, a flow, a rate or a concentration; the rows worth looking at are the fractions, the few which may be negative and the switches.

<!-- generated: domains-kinds (scripts/pbpk_tables.py), do not edit -->
*Generated from `KINDS` of `pkpdsim.pbpk.annotation` by `scripts/pbpk_tables.py`; a domain is classified by hand, no database entry stands behind it.*

| kind | domain | PBPKO term |
| --- | --- | --- |
| `plasma` | `[0, inf)` | `PLASMA_COMPARTMENT` |
| `blood_cells` | `[0, inf)` | `RED_BLOOD_CELL_COMPARTMENT` |
| `capacity` | `[0, inf)` | `COMPARTMENT` |
| `organ` | `[0, inf)` | `COMPARTMENT` |
| `tissue` | `[0, inf)` | `COMPARTMENT` |
| `interstitial` | `[0, inf)` | `COMPARTMENT` |
| `intracellular` | `[0, inf)` | `COMPARTMENT` |
| `urine` | `[0, inf)` | `URINE_COMPARTMENT` |
| `feces` | `[0, inf)` | `FECES_COMPARTMENT` |
| `remodelled_tissue` | `[0, inf)` | `COMPARTMENT` |
| `bile` | `[0, inf)` | `BILE_COMPARTMENT` |
| `lumen` | `[0, inf)` | `GUT_LUMEN_COMPARTMENT` |
| `mucosa` | `[0, inf)` | `COMPARTMENT` |
| `gut_wall` | `[0, inf)` | `COMPARTMENT` |
| `tubule` | `[0, inf)` | `FILTRATE_COMPARTMENT` |
| `body_weight` | `[0, inf)` | `BODYWEIGHT` |
| `body_weight_change` | `(-inf, inf)` | `BODYWEIGHT` |
| `height` | `[0, inf)` | `HEIGHT` |
| `age` | `[0, inf)` | `AGE_IN_PBPK_MODEL` |
| `bmi` | `[0, inf)` | `BODY_MASS_INDEX` |
| `bsa` | `[0, inf)` | `SURFACE_AREA` |
| `lean_mass` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `scaling_factor` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `organ_mass` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `organ_volume` | `[0, inf)` | `VOLUME_OF_COMPARTMENT` |
| `blood_volume` | `[0, inf)` | `VOLUME_OF_BLOOD` |
| `pool_volume` | `[0, inf)` | `VOLUME_OF_COMPARTMENT` |
| `pool_fraction` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `vascular_fraction` | `[0, 1]` | `VOLUME_FRACTION_OF_BLOOD` |
| `interstitial_fraction` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `tissue_volume` | `[0, inf)` | `VOLUME_OF_COMPARTMENT` |
| `zone_fraction` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `lobe_volume` | `[0, inf)` | `VOLUME_OF_COMPARTMENT` |
| `lobe_fraction` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `lobe_supply_fraction` | `[0, 1]` | `FRACTION_OF_BLOOD_FLOW_TO_A_COMPARTMENT` |
| `lobe_flow` | `[0, inf)` | `BLOOD_FLOW_RATE_TO_COMPARTMENT` |
| `lobe_perfusion` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `lobe_growth` | `(-inf, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `remodelling_time` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `lobe_patency` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `lobe_present` | `[0, 1] boolean` | `PHYSIOLOGICAL_PARAMETER` |
| `remnant_fraction` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `lobe_pressure` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `lobe_resistance` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `buffer_response` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `hematocrit` | `[0, 1]` | `VOLUME_FRACTION_OF_RED_BLOOD_CELLS` |
| `cardiac_output` | `[0, inf)` | `CARDIAC_OUTPUT_RATE` |
| `flow_fraction` | `[0, 1]` | `FRACTION_OF_BLOOD_FLOW_TO_A_COMPARTMENT` |
| `organ_flow` | `[0, inf)` | `BLOOD_FLOW_RATE_TO_COMPARTMENT` |
| `hepatic_artery_flow` | `[0, inf)` | `BLOOD_FLOW_RATE_TO_COMPARTMENT` |
| `portal_flow` | `[0, inf)` | `BLOOD_FLOW_RATE_TO_PORTAL_VEIN` |
| `hepatic_flow` | `[0, inf)` | `HEPATIC_VEIN_BLOOD_FLOW` |
| `gfr` | `[0, inf)` | `GLOMERULAR_FILTRATION_RATE` |
| `gfr_base` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `egfr` | `[0, inf)` | `GLOMERULAR_FILTRATION_RATE` |
| `organ_function` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `tubular_share` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `shunt` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `child_pugh` | `[0, 3]` | `PHYSIOLOGICAL_PARAMETER` |
| `heart_failure` | `[0, 4]` | `PHYSIOLOGICAL_PARAMETER` |
| `ckd` | `[0, 1] boolean` | `PHYSIOLOGICAL_PARAMETER` |
| `cirrhosis_factor` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `ckd_factor` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `heart_failure_factor` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `impairment_factor` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `age_factor` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `sex_factor` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `body_composition` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `body_composition_change` | `(-inf, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `fat_mass` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `adipose_depot` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `obesity_factor` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `weight_change` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `gestational_age` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `pregnancy` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `placenta` | `[0, inf)` | `COMPARTMENT` |
| `fetus` | `[0, inf)` | `COMPARTMENT` |
| `clock_time` | `[0, 24]` | `PHYSIOLOGICAL_PARAMETER` |
| `circadian` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `molecular_weight` | `[0, inf)` | `MOLECULAR_WEIGHT` |
| `fu` | `[0, 1]` | `FRACTION_UNBOUND_PLASMA` |
| `fraction_reabsorbed` | `[0, 1]` | `RESORPTION` |
| `reabsorption_threshold` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `reabsorption_n` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `baseline` | `[0, inf)` | `CONCENTRATION_IN_PLASMA` |
| `steady_clearance` | `[0, inf)` | `CLEARANCE_RATE` |
| `steady_rate` | `(-inf, inf)` | `BIOCHEMICAL_PARAMETER` |
| `steady_vmax` | `[0, inf)` | `MAXIMUM_RATE` |
| `steady_km` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `steady_plasma` | `[0, inf)` | `CONCENTRATION_IN_PLASMA` |
| `kp` | `[0, inf)` | `TISSUE_PLASMA_PARTITION_COEFFICIENT` |
| `kp_organ` | `[0, inf)` | `TISSUE_PLASMA_PARTITION_COEFFICIENT` |
| `kp_factor` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `partition` | `[0, inf)` | `PARTITION_COEFFICIENT` |
| `slow_pool_exchange` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `tissue_composition` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `concentration_ratio` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `permeability` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `pore_permeability` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `diffusion` | `[0, inf)` | `DIFFUSION_COEFFICIENT` |
| `sinusoid_geometry` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `bp` | `[0, inf)` | `BLOOD_PLASMA_RATIO` |
| `kbc` | `[0, inf)` | `RED_BLOOD_CELLS_PARTITION_COEFFICIENT` |
| `exchange_rate` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `endothelial_exchange` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `activity` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `enzyme_amount` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `kdeg` | `[0, inf)` | `DEGRADATION_RATE_CONSTANT` |
| `ki` | `[0, inf)` | `INHIBITION_CONSTANT` |
| `kinact` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `emax` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `ec50` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `placement_weight` | `[0, 1]` | `MISCELLANEOUS_PARAMETER` |
| `relative_abundance` | `[0, inf)` | `MISCELLANEOUS_PARAMETER` |
| `vmax` | `[0, inf)` | `MAXIMUM_RATE` |
| `vmax_basolateral` | `[0, inf)` | `MAXIMUM_RATE_OF_BASOLATERAL_TRANSPORTER` |
| `vmax_apical` | `[0, inf)` | `MAXIMUM_RATE_OF_APICAL_TRANSPORTER` |
| `km` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `km_basolateral` | `[0, inf)` | `MICHAELIS_CONSTANT_BASO` |
| `km_apical` | `[0, inf)` | `MICHAELIS_CONSTANT_APICAL` |
| `hill` | `[0, inf)` | `HILL_COEFFICIENT` |
| `k_metabolism` | `[0, inf)` | `INTRINSIC_CLEARANCE_RATE` |
| `k_influx` | `[0, inf)` | `CLEARANCE_RATE` |
| `k_efflux` | `[0, inf)` | `CLEARANCE_RATE` |
| `k_metabolism_volume` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `k_transport_volume` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `rate_ratio` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `product_fraction` | `[0, 1]` | `BIOCHEMICAL_PARAMETER` |
| `custom_parameter` | none | `MISCELLANEOUS_PARAMETER` |
| `definition_parameter` | none | `MISCELLANEOUS_PARAMETER` |
| `custom_species` | `[0, inf)` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `custom_amount` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `dose` | `[0, inf)` | `PARAMETER` |
| `iv_time` | `[0, inf)` | `PARAMETER` |
| `iv_rate` | `[0, inf)` | `PARAMETER` |
| `infusion_rate` | `[0, inf)` | `PARAMETER` |
| `dissolution` | `[0, inf)` | `MISCELLANEOUS_PARAMETER` |
| `application_rate` | `[0, inf)` | `ABSORPTION_RATE_CONSTANT` |
| `absorption_rate` | `[0, inf)` | `ABSORPTION_RATE_CONSTANT_INTESTINE` |
| `transit_time` | `[0, inf)` | `INTESTINAL_TRANSIT_TIME` |
| `gastric_transit_time` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `fed` | `[0, 1]` | `MISCELLANEOUS_PARAMETER` |
| `intestine_geometry` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `gut_flow` | `[0, inf)` | `BLOOD_FLOW_RATE_TO_GUT` |
| `gut_volume` | `[0, inf)` | `VOLUME_OF_GUT` |
| `mucosa_flow` | `[0, inf)` | `BLOOD_FLOW_RATE_TO_COMPARTMENT` |
| `mucosa_volume` | `[0, inf)` | `ENTEROCYTE_VOLUME` |
| `bile_transit` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `tubule_fraction` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `fraction_absorbed` | `[0, 1]` | `UPTAKE_FRACTION_TO_BLOOD` |
| `colon_absorption` | `[0, inf)` | `PHYSIOLOGICAL_PARAMETER` |
| `regional_absorption` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `lumen_conversion_rate` | `[0, inf)` | `BIOCHEMICAL_PARAMETER` |
| `arm_weight` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `outflow` | `[0, inf)` | `TOTAL_VENOUS_CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `outflow_blood_cells` | `[0, inf)` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `outflow_lobe` | `[0, inf)` | `TOTAL_VENOUS_CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `organ_concentration` | `[0, inf)` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `organ_amount` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `plasma_concentration` | `[0, inf)` | `CONCENTRATION_IN_PLASMA` |
| `free_concentration` | `[0, inf)` | `UNBOUND_CONCENTRATION_IN_PLASMA` |
| `blood_concentration` | `[0, inf)` | `CONCENTRATION_IN_BLOOD` |
| `unbound_concentration` | `[0, inf)` | `UNBOUND_CONCENTRATION_OF_COMPUND_IN_COMPARTMENT` |
| `urine_rate` | `[0, inf)` | `OUTPUT_PARAMETER` |
| `lumen_amount` | `[0, inf)` | `AMOUNT_IN_LUMEN` |
| `intestine_absorption` | `[0, inf)` | `OUTPUT_PARAMETER` |
| `intestine_efflux` | `[0, inf)` | `OUTPUT_PARAMETER` |
| `feces_rate` | `[0, inf)` | `OUTPUT_PARAMETER` |
| `absorbed_amount` | `[0, inf)` | `OUTPUT_PARAMETER` |
| `absorbed_fraction` | `[0, inf)` | `OUTPUT_PARAMETER` |
| `oral_dose_applied` | `[0, inf)` | `OUTPUT_PARAMETER` |
| `observer` | `(-inf, inf)` | `OUTPUT_PARAMETER` |
| `concentration` | `[0, inf)` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `concentration_plasma` | `[0, inf)` | `CONCENTRATION_IN_PLASMA` |
| `concentration_tissue` | `[0, inf)` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `concentration_lumen` | `[0, inf)` | `CONCENTRATION_IN_LUMEN` |
| `concentration_mucosa` | `[0, inf)` | `CONCENTRATION_IN_INTESTINAL_TISSUE` |
| `concentration_bile` | `[0, inf)` | `CONCENTRATION_IN_BILE` |
| `concentration_tubule` | `[0, inf)` | `CONCENTRATION_IN_FILTRATE` |
| `amount_urine` | `[0, inf)` | `CUMULATIVE_AMOUNT_IN_URINE` |
| `amount_feces` | `[0, inf)` | `CUMULATIVE_AMOUNT_IN_FECES` |
| `amount_undissolved` | `[0, inf)` | `AMOUNT_IN_LUMEN` |
| `amount_remodelled` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `blood_flow` | none | `DISTRIBUTION_PROCESS` |
| `exchange` | none | `DISTRIBUTION_PROCESS` |
| `remodelling` | none | `DISTRIBUTION_PROCESS` |
| `metabolism` | none | `LINEAR_ELIMINATION` |
| `metabolism_saturable` | none | `NON_LINEAR_ELIMINATION` |
| `gut_metabolism` | none | `GUT_METABOLISM` |
| `metabolism_product` | none |  |
| `transport` | none | `DISTRIBUTION_PROCESS` |
| `biliary_excretion` | none | `BILIARY_ELIMINATION` |
| `renal_secretion` | none | `RENAL_EXCRETION` |
| `reabsorption` | none | `RESORPTION` |
| `filtration` | none | `RENAL_EXCRETION` |
| `tubular_flow` | none |  |
| `absorption` | none | `ABSORPTION` |
| `dissolution_reaction` | none |  |
| `transit` | none |  |
| `fecal_excretion` | none | `FECAL_EXCRETION` |
| `lumen_conversion` | none |  |
| `bile_flow` | none |  |
| `iv_application` | none |  |
| `oral_application` | none | `ORAL_EXPOSURE` |
| `depot_site_im` | `[0, inf)` | `COMPARTMENT` |
| `depot_amount_im` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `depot_time_im` | `[0, inf)` | `PARAMETER` |
| `depot_rate_im` | `[0, inf)` | `ABSORPTION_RATE_CONSTANT` |
| `depot_fraction_im` | `[0, 1]` | `UPTAKE_FRACTION_TO_BLOOD` |
| `depot_application_im` | none |  |
| `depot_absorption_im` | none | `ABSORPTION` |
| `depot_loss_im` | none |  |
| `depot_loss_amount_im` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `depot_site_sc` | `[0, inf)` | `COMPARTMENT` |
| `depot_amount_sc` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `depot_time_sc` | `[0, inf)` | `PARAMETER` |
| `depot_rate_sc` | `[0, inf)` | `ABSORPTION_RATE_CONSTANT` |
| `depot_fraction_sc` | `[0, 1]` | `UPTAKE_FRACTION_TO_BLOOD` |
| `depot_application_sc` | none |  |
| `depot_absorption_sc` | none | `ABSORPTION` |
| `depot_loss_sc` | none |  |
| `depot_loss_amount_sc` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `depot_site_rectal` | `[0, inf)` | `COMPARTMENT` |
| `depot_amount_rectal` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `depot_time_rectal` | `[0, inf)` | `PARAMETER` |
| `depot_rate_rectal` | `[0, inf)` | `ABSORPTION_RATE_CONSTANT` |
| `depot_fraction_rectal` | `[0, 1]` | `UPTAKE_FRACTION_TO_BLOOD` |
| `depot_application_rectal` | none |  |
| `depot_absorption_rectal` | none | `ABSORPTION` |
| `depot_loss_rectal` | none | `FECAL_EXCRETION` |
| `depot_site_inhaled` | `[0, inf)` | `COMPARTMENT` |
| `depot_amount_inhaled` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `depot_time_inhaled` | `[0, inf)` | `PARAMETER` |
| `depot_rate_inhaled` | `[0, inf)` | `ABSORPTION_RATE_CONSTANT` |
| `depot_fraction_inhaled` | `[0, 1]` | `UPTAKE_FRACTION_TO_BLOOD` |
| `depot_application_inhaled` | none |  |
| `depot_absorption_inhaled` | none | `ABSORPTION` |
| `depot_loss_inhaled` | none |  |
| `depot_loss_amount_inhaled` | `[0, inf)` | `AMOUNT_OF_COMPOUND_IN_COMPARTMENT` |
| `rectal_portal_fraction` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `inhaled_swallowed_fraction` | `[0, 1]` | `PHYSIOLOGICAL_PARAMETER` |
| `inhaled_swallow_rate` | `[0, inf)` | `OUTPUT_PARAMETER` |
| `custom_reaction` | none |  |
<!-- end generated -->

A raw sbmlutils `Parameter` in a definition is the model writing the element itself, so it carries no domain unless the model writes one (it attaches `domain_uncertainty` of `pkpdsim.sbml.domains` to the parameter). The two kinds of a constant which the model states, `definition_parameter` and `custom_parameter`, state no domain: the `ModelParameter` of the definition or of a custom reaction states its own (`domain=`), and `validate` refuses one which states none or a value outside it.

## Read the domains of a model

`domains(model_path)` of `pkpdsim.sbml.domains` reads the domain of every parameter, compartment and species of an SBML file which states one. The module imports libsbml and sbmlutils only, no simulator, so a model package and the other libraries can read it. The model of the snippet is built with one `ModelParameter` of its own, a fraction absorbed:

```python
from pathlib import Path

from pkpdsim.pbpk import (
    Binding,
    Excretion,
    Kp,
    ModelParameter,
    PBPKDefinition,
    Route,
    Substance,
    U,
    build,
)
from pkpdsim.sbml.domains import FRACTION, domains

definition = PBPKDefinition(
    sid="api_body",
    substances=[
        Substance(
            "api",
            name="apixaban",
            mr=459.5,
            binding=Binding(fu=0.13, protein="albumin"),
            distribution=Kp(1.2, ftissue=5.0),
            excretion=Excretion(urine=True),
        )
    ],
    routes=[Route.iv("api"), Route.oral("api")],
    parameters=[
        ModelParameter(
            "fabs_api", 0.94, U.dimensionless, name="fraction absorbed", domain=FRACTION
        )
    ],
)
model_path = build(definition, Path("models")).flat_path
by_id = domains(model_path)
for sid in ["fu_api", "fabs_api", "BW", "Vli", "Cve_api", "CP_severity", "CKD_hepatic"]:
    print(sid, by_id[sid])
print(len(by_id), by_id.get("J_iv_api"))
```

```text
fu_api [0, 1]
fabs_api [0, 1]
BW [0, inf)
Vli [0, inf)
Cve_api [0, inf)
CP_severity [0, 3]
CKD_hepatic [0, 1] boolean
403 None
```

<!-- generated: figures-domains-1 (scripts/visual/docs_figures.py), do not edit -->
The content of the model `api_body`: the anatomical figure (A) and the schematic (B), hover them for details ([figures](figures.md)).

<iframe class="pkpd-figure" src="img/figures/models/domains-1-anatomy.html" title="The anatomical figure of the model api_body" style="aspect-ratio:1460 / 885" loading="lazy"></iframe>

<iframe class="pkpd-figure" src="img/figures/models/domains-1-schematic.html" title="The schematic of the model api_body" style="aspect-ratio:1314 / 882" loading="lazy"></iframe>
<!-- end generated -->

The model has 403 parameters, compartments and species, and all of them state a domain. The reaction `J_iv_api` states none and reads as `None`, the silence of the model.

## Clip a span

`library_parameters` of `pkpdsim.sensitivity` intersects the bounds of every parameter of a sensitivity analysis with its domain (`clip=True`, the default). Every bound it moves is logged as a warning, and the parameters which state no domain are counted in one further warning. The severities (`SEVERITY_IDS`, `CP_severity` and `HF_severity`) are varied over their domain, from 0 to `relative` times its upper bound, and not around their value, which is 0 for a healthy individual. A parameter whose value is outside its own domain is an error, the clip does not hide it. `clip=False` keeps the span for the whole analysis and `unclipped` for the named parameters.

The span of 0.15 around 0.94 gives 0.799 to 1.081 and the upper bound is moved to 1. The unbound fraction `fu_api` of 0.13 gives 0.1105 to 0.1495 and stays inside its domain, and the severity of the Child-Pugh class gets 0 to 0.15 times 3, which is 0.45:

```python
import logging
import sys

from pkpdsim.sensitivity import library_parameters

logging.basicConfig(
    level=logging.WARNING, format="%(levelname)s %(message)s", stream=sys.stdout
)
ids = ["fabs_api", "fu_api", "CP_severity"]
for p in library_parameters(model_path, ids, fit_ids=[], scaling_ids=[]):
    print(p.uid, p.value, round(p.lower_bound, 4), round(p.upper_bound, 4))
print("--")
unclipped = library_parameters(
    model_path, ids, fit_ids=[], scaling_ids=[], unclipped=["fabs_api"]
)
print(unclipped[0].uid, round(unclipped[0].upper_bound, 4))
```

```text
WARNING The bounds [0.799, 1.081] of the sensitivity parameter fabs_api leave its domain [0, 1] and are clipped to [0.799, 1].
fabs_api 0.94 0.799 1.0
fu_api 0.13 0.1105 0.1495
CP_severity 0.0 0.0 0.45
--
fabs_api 1.081
```

The clip reaches the sampled analyses (`SamplingAnalysis`, `FastAnalysis`, `SobolAnalysis`) through the bounds of the parameter table, and not the local analysis, which perturbs a parameter by `delta` relative to its own value: a fraction at 0.99 with a `delta` of 0.15 still leaves its domain there.

## The domain itself

A `Domain` is a closed interval of two floats, infinities allowed, with the value set. `clip` limits each bound of an interval on its own side, `(max(lower, lo), min(upper, hi))`, so an interval whose bounds are reversed is returned as it is, and a caller may rely on that: an interval which lies wholly above the domain comes back as `(2.0, 1.0)`. A single value is limited to the domain and, for an integer or boolean domain, rounded to an integer it contains. `contains` rejects a value which is not finite, so `[0, inf)` means what it renders: infinity is no value of it. A domain of the integers or of the truth values whose interval holds no integer is rejected when it is constructed.

```python
import math

from pkpdsim.sbml.domains import BOOLEAN, FRACTION, NONNEGATIVE, integers

print(FRACTION.clip(0.799, 1.081), FRACTION.clip(2.0, 3.0))
print(NONNEGATIVE.contains(3.0), NONNEGATIVE.contains(math.inf))
print(integers(0, 3), BOOLEAN)
try:
    integers(0.2, 0.8)
except ValueError as err:
    print(err)
```

```text
(0.799, 1.0) (2.0, 1.0)
True False
[0, 3] integer [0, 1] boolean
the interval [0.2, 0.8] holds no integer, so no integer domain has a value
```

## What is checked and what is not

The table is a classification by hand, read off the meaning of each quantity (a fraction is [0, 1] because it is a fraction), never fitted to the values which happen to be in the reference models. `tests/pbpk/test_domains_build.py` is the evidence for it. It builds 20 shapes of a model, reads every quantity which has a domain at t=0 and every 30 minutes along a time course, and requires it inside its domain; it found no wrong domain and no wrong value. The shapes reach 204 of the 207 kinds which state a domain. The three it does not reach are `buffer_response`, `lobe_pressure` and `lobe_resistance`, all of the `PressureFlow` network of the lobes, whose algebraic rules roadrunner cannot simulate, so no shape can read their values.

The check has two limits, and a green test is no proof of more. The sweep shows that no domain is too narrow for the values these models actually reach; it cannot show that an upper bound is too wide when no value ever approaches it. And it reads sample points and not an envelope, so a brief excursion between two points is not seen; a monotone quantity like the fraction absorbed of pkpdsim#207 is caught either way. The sweep covers the individuals these 20 shapes build, so a domain can still be wrong for a covariate combination no shape instantiates.

The domain is also not enforced while a model runs: roadrunner ignores SBML constraints, and nothing in the stack stops a simulation whose concentration goes negative. The domain says what the model accepts, it does not police a simulation.
