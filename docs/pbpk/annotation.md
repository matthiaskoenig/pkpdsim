# Annotation

Every compartment, species, parameter and reaction which the PBPK library creates carries an SBO term and ontology terms (#63). The generators annotate an element when they create it, with `annotate(element, kind, ...)` of `pkpdsim.pbpk.annotation`: the kind names what the element is (an organ volume, the plasma of an organ, a saturable metabolism) and selects one row of the table `KINDS`, which holds the SBO term and the ontology terms of the kind. The place of the element adds the terms of where it is:

- the organ: its UBERON, NCIT, FMA and BTO terms (the compartment terms of `pkpdsim.annotations`) and its PBPKO compartment (`LIVER_COMPARTMENT`);
- the blood pool: arterial or venous blood, the portal vein;
- the lumen segment of the intestine or the segment of the renal tubule;
- the sampling site of an observer (hepatic vein, arm vein);
- the protein of a process or of its parameters: its UniProt and NCIT terms (the protein table `data/proteins.toml`);
- the substance of a species: the annotations of its specification (`Substance.annotations`, ChEBI, PubChem, NCIT).

An input of an organ or of a PD component, which the element of the body replaces in the comp model, takes the annotation of that element, so the comp model and the flat model carry the same terms.

## Ontology and qualifiers

The PBPKO terms are those of the pinned pymetadata (`pymetadata.ontologies.pbpko`). A kind takes the organ-specific PBPKO term where PBPKO has one for the organ (`VOLUME_OF_LIVER`, `BLOOD_FLOW_RATE_TO_KIDNEY`, `LIVER_PLASMA_PARTITION_COEFFICIENT`, `MAXIMUM_RATE_IN_LIVER`); the small and the large intestine take the term of the intestine where PBPKO has none of their own (`MAXIMUM_RATE_IN_INTESTINE`), adipose tissue is `FAT` in PBPKO.

| qualifier | meaning |
|---|---|
| `BQB.IS` | what the element is: its PBPKO term, the terms of its kind, the protein of a process, the substance of a species |
| `BQB.IS_PART_OF` | the place of a compartment |
| `BQB.IS_PROPERTY_OF` | the place of a parameter, the protein of an activity or an inhibition constant |
| `BQB.OCCURS_IN` | the place of a reaction |

The model carries `PHYSIOLOGICALLY_BASED_PHARMACOKINETIC_MODEL` and `WHOLE_BODY_PBPK`, the human taxon, the ODE model of MAMO and pharmacokinetics of NCIT. It is a `PERFUSION_LIMITED_MODEL` when every substance which enters tissue does so flow-limited (`Mode.PERFUSION`), a `PERMEABILITY_LIMITED_PBPK` when one has a finite exchange, and neither when no substance enters tissue.

## Resource format

Every resource of a collection which identifiers.org registers is written as its identifiers.org `collection/term`, the format of `pkpdsim.annotations`: `uberon/UBERON:0002107`, `uniprot/P08684`; sbmlutils writes the identifiers.org url. identifiers.org does not register PBPKO, so pymetadata cannot parse `pbpko/PBPKO:00558` (it logs a MIRIAM error and writes an identifiers.org url which does not resolve). A PBPKO term is therefore written as its bioregistry url, `https://bioregistry.io/PBPKO:00558`, which resolves and which sbmlutils writes as given (as the apixaban model of pkdb_models does). A build logs no annotation error.

## Terms requested upstream

PBPKO has no term for some of the quantities, compartments and processes of the library. The library annotates such an element with the nearest PBPKO term which is true for it: the superclass (`PHYSIOLOGICAL_PARAMETER`, `BIOCHEMICAL_PARAMETER`, `MISCELLANEOUS_PARAMETER`, `COMPARTMENT`), never a specific term whose definition contradicts the element (a fraction per body weight for a fraction of the blood volume, a flow rate for a fraction of the GFR, a clearance (a volume per time) for a clearance per litre of organ). A more specific PBPKO term which the element is a part or a version of is added with `BQB.IS_PART_OF` or `BQB.IS_VERSION_OF` (the mucosa is part of the `INTESTINAL_TISSUE_COMPARTMENT`, the dissolution part of the `ABSORPTION`). The missing terms are requested from PBPKO:

| missing term | broader term used |
|---|---|
| plasma-equivalent capacity of the blood (volume times BP) | `COMPARTMENT` |
| tissue compartment (interstitial and intracellular space) | `COMPARTMENT` |
| interstitial compartment | `COMPARTMENT` |
| intracellular space compartment | `COMPARTMENT` |
| intestinal mucosa compartment (enterocytes) | `COMPARTMENT` |
| intestinal tissue compartment without the mucosa | `COMPARTMENT` |
| fat-free mass | `PHYSIOLOGICAL_PARAMETER` |
| covariate scaling factor | `PHYSIOLOGICAL_PARAMETER` |
| mass of organ | `PHYSIOLOGICAL_PARAMETER` |
| fraction of the blood volume in a blood pool | `PHYSIOLOGICAL_PARAMETER` |
| interstitial volume fraction of an organ | `PHYSIOLOGICAL_PARAMETER` |
| volume fraction of a zone of the liver | `PHYSIOLOGICAL_PARAMETER` |
| hepatic artery blood flow | `BLOOD_FLOW_RATE_TO_COMPARTMENT` |
| glomerular filtration rate at normal renal function | `PHYSIOLOGICAL_PARAMETER` |
| relative organ function | `PHYSIOLOGICAL_PARAMETER` |
| share of the tubular capacity which follows the renal function | `PHYSIOLOGICAL_PARAMETER` |
| portosystemic shunt fraction | `PHYSIOLOGICAL_PARAMETER` |
| severity of hepatic impairment | `PHYSIOLOGICAL_PARAMETER` |
| severity of heart failure | `PHYSIOLOGICAL_PARAMETER` |
| stage of chronic kidney disease | `PHYSIOLOGICAL_PARAMETER` |
| impairment factor | `PHYSIOLOGICAL_PARAMETER` |
| endogenous baseline plasma concentration | `CONCENTRATION_IN_PLASMA` |
| factor of the partition coefficients | `BIOCHEMICAL_PARAMETER` |
| exchange rate constant of a slow tissue pool per tissue volume | `BIOCHEMICAL_PARAMETER` |
| ratio of the capacity of a process to another process | `BIOCHEMICAL_PARAMETER` |
| fraction of a metabolism which forms a product | `BIOCHEMICAL_PARAMETER` |
| renal threshold of a saturable reabsorption | `BIOCHEMICAL_PARAMETER` |
| sharpness of the renal threshold of a saturable reabsorption | `BIOCHEMICAL_PARAMETER` |
| tissue composition (water, protein, lipid fractions) | `PHYSIOLOGICAL_PARAMETER` |
| permeability surface area product (exchange flow) | `BIOCHEMICAL_PARAMETER` |
| exchange rate constant of plasma and blood cells | `BIOCHEMICAL_PARAMETER` |
| endothelial exchange rate constant per organ volume | `PHYSIOLOGICAL_PARAMETER` |
| relative activity of enzyme or transporter | `BIOCHEMICAL_PARAMETER` |
| relative amount of enzyme or transporter | `BIOCHEMICAL_PARAMETER` |
| maximal inactivation rate constant | `BIOCHEMICAL_PARAMETER` |
| maximal induction | `BIOCHEMICAL_PARAMETER` |
| half-maximal effective concentration | `BIOCHEMICAL_PARAMETER` |
| abundance weight of a protein or process | `MISCELLANEOUS_PARAMETER` |
| Michaelis constant | `BIOCHEMICAL_PARAMETER` |
| intrinsic clearance per organ volume | `BIOCHEMICAL_PARAMETER` |
| transport clearance per organ volume | `BIOCHEMICAL_PARAMETER` |
| dose | `PARAMETER` |
| injection time of an intravenous dose | `PARAMETER` |
| injection rate constant of an intravenous dose | `PARAMETER` |
| infusion rate of an intravenous infusion | `PARAMETER` |
| dissolution of a formulation | `MISCELLANEOUS_PARAMETER` |
| gastric emptying time | `PHYSIOLOGICAL_PARAMETER` |
| prandial state (fed or fasted) | `MISCELLANEOUS_PARAMETER` |
| intestinal mucosa and absorption geometry | `PHYSIOLOGICAL_PARAMETER` |
| absorption factor of the large intestine | `PHYSIOLOGICAL_PARAMETER` |
| relative absorption of a substance in an intestinal segment | `BIOCHEMICAL_PARAMETER` |
| rate constant of a conversion in the intestinal lumen | `BIOCHEMICAL_PARAMETER` |
| blood flow rate to intestinal mucosa | `BLOOD_FLOW_RATE_TO_COMPARTMENT` |
| bile transit rate constant | `PHYSIOLOGICAL_PARAMETER` |
| tubular flow fraction of the glomerular filtration rate | `PHYSIOLOGICAL_PARAMETER` |
| weight of a tissue in the arm vein sampling site | `PHYSIOLOGICAL_PARAMETER` |
| urinary excretion rate (amount of compound per time) | `OUTPUT_PARAMETER` |
| intestinal absorption rate (amount of compound per time) | `OUTPUT_PARAMETER` |
| intestinal efflux rate into the lumen (amount of compound per time) | `OUTPUT_PARAMETER` |
| fecal excretion rate (amount of compound per time) | `OUTPUT_PARAMETER` |
| `CONCENTRATION_IN_HEART` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `CONCENTRATION_IN_PANCREAS` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `CONCENTRATION_IN_MUSCLE` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `CONCENTRATION_IN_BONE` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `CONCENTRATION_IN_GONAD` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| amount of undissolved drug | `AMOUNT_IN_LUMEN` |
| dissolution process | none, the terms of other ontologies only |
| gastrointestinal transit process | none, the terms of other ontologies only |
| conversion in the intestinal lumen | none, the terms of other ontologies only |
| bone compartment | none, the terms of other ontologies only |
| portal vein blood compartment | none, the terms of other ontologies only |

PBPKO is further not registered at identifiers.org, so its terms are written as bioregistry urls.
