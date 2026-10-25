# Annotation

Every compartment, species, parameter and reaction which the PBPK library
creates carries an SBO term and ontology terms (#63). The generators annotate
an element when they create it, with `annotate(element, kind, ...)` of
`pkpdsim.pbpk.annotation`: the kind names what the element is (an organ
volume, the plasma of an organ, a saturable metabolism) and selects one row of
the table `KINDS`, which holds the SBO term and the ontology terms of the kind.
The place of the element adds the terms of where it is:

- the organ: its UBERON, NCIT, FMA and BTO terms (the compartment terms of
  `pkpdsim.annotations`) and its PBPKO compartment (`LIVER_COMPARTMENT`);
- the blood pool: arterial or venous blood, the portal vein;
- the lumen segment of the intestine or the segment of the renal tubule;
- the sampling site of an observer (hepatic vein, arm vein);
- the protein of a process or of its parameters: its UniProt and NCIT terms
  (the protein table `data/proteins.toml`);
- the substance of a species: the annotations of its specification
  (`Substance.annotations`, ChEBI, PubChem, NCIT).

An input of an organ or of a PD component, which the element of the body
replaces in the comp model, takes the annotation of that element, so the comp
model and the flat model carry the same terms.

## Ontology and qualifiers

The PBPKO terms are those of the pinned pymetadata (`pymetadata.ontologies.pbpko`).
A kind takes the organ-specific PBPKO term where PBPKO has one for the organ
(`VOLUME_OF_LIVER`, `BLOOD_FLOW_RATE_TO_KIDNEY`,
`LIVER_PLASMA_PARTITION_COEFFICIENT`, `MAXIMUM_RATE_IN_LIVER`); the small and
the large intestine take the term of the intestine where PBPKO has none of
their own (`MAXIMUM_RATE_IN_INTESTINE`), adipose tissue is `FAT` in PBPKO.

| qualifier | meaning |
|---|---|
| `BQB.IS` | what the element is: its PBPKO term, the terms of its kind, the protein of a process, the substance of a species |
| `BQB.IS_PART_OF` | the place of a compartment |
| `BQB.IS_PROPERTY_OF` | the place of a parameter, the protein of an activity or an inhibition constant |
| `BQB.OCCURS_IN` | the place of a reaction |

The model carries `PHYSIOLOGICALLY_BASED_PHARMACOKINETIC_MODEL` and
`WHOLE_BODY_PBPK`, the human taxon, the ODE model of MAMO and
pharmacokinetics of NCIT. It is a `PERFUSION_LIMITED_MODEL` when every
substance which enters tissue does so flow-limited (`Mode.PERFUSION`), a
`PERMEABILITY_LIMITED_PBPK` when one has a finite exchange, and neither when
no substance enters tissue.

## Resource format

Every resource is written in one format, the identifiers.org
`collection/term` of `pkpdsim.annotations`: `pbpko/PBPKO:00558`,
`uberon/UBERON:0002107`, `uniprot/P08684`. sbmlutils writes them as the
identifiers.org urls of pymetadata, `https://identifiers.org/PBPKO:00558`;
identifiers.org does not register PBPKO yet (bioregistry does), the url is the
`url` of the pymetadata term.

## Terms requested upstream

PBPKO has no term for some of the quantities, compartments and processes of
the library. The library annotates such an element with the nearest broader
PBPKO term (and its SBO term and the terms of other ontologies) and requests
the missing term from PBPKO:

| missing term | broader term used |
|---|---|
| tissue compartment (interstitial and intracellular space) | `COMPARTMENT` |
| interstitial compartment | `COMPARTMENT` |
| intracellular space compartment | `COMPARTMENT` |
| intestinal mucosa compartment (enterocytes) | `INTESTINAL_TISSUE_COMPARTMENT` |
| fat-free mass | `PHYSIOLOGICAL_PARAMETER` |
| covariate scaling factor | `PHYSIOLOGICAL_PARAMETER` |
| mass of organ | `PHYSIOLOGICAL_PARAMETER` |
| fraction of venous and of portal blood volume | `FRACTION_OF_VOLUME_OF_COMPARTMENT` |
| interstitial volume fraction of organ | `FRACTION_OF_VOLUME_OF_COMPARTMENT` |
| volume fraction of a zone of the liver | `FRACTION_OF_VOLUME_OF_COMPARTMENT` |
| hepatic artery blood flow | `BLOOD_FLOW_RATE_TO_COMPARTMENT` |
| relative organ function | `PHYSIOLOGICAL_PARAMETER` |
| portosystemic shunt fraction | `PHYSIOLOGICAL_PARAMETER` |
| severity of hepatic impairment | `PHYSIOLOGICAL_PARAMETER` |
| severity of heart failure | `PHYSIOLOGICAL_PARAMETER` |
| stage of chronic kidney disease | `PHYSIOLOGICAL_PARAMETER` |
| impairment factor | `PHYSIOLOGICAL_PARAMETER` |
| factor of the partition coefficients | `PARTITION_COEFFICIENT` |
| tissue composition (water, protein, lipid fractions) | `PHYSIOLOGICAL_PARAMETER` |
| permeability surface area product | `APPARENT_PERMEABILITY` |
| exchange rate constant | `DIFFUSION_RATE` |
| relative activity of enzyme or transporter | `BIOCHEMICAL_PARAMETER` |
| relative amount of enzyme or transporter | `BIOCHEMICAL_PARAMETER` |
| maximal inactivation rate constant | `BIOCHEMICAL_PARAMETER` |
| maximal induction | `BIOCHEMICAL_PARAMETER` |
| half-maximal effective concentration | `BIOCHEMICAL_PARAMETER` |
| abundance weight of a protein or process | `MISCELLANEOUS_PARAMETER` |
| Michaelis constant | `BIOCHEMICAL_PARAMETER` |
| dose | `PARAMETER` |
| injection time of an intravenous dose | `PARAMETER` |
| injection rate constant of an intravenous dose | `PARAMETER` |
| infusion rate of an intravenous infusion | `PARAMETER` |
| dissolution of a formulation | `PHYSICOCHEMICAL_PARAMETER` |
| fed state | `MISCELLANEOUS_PARAMETER` |
| intestinal mucosa and absorption geometry | `FRACTION_OF_INTESTINAL_TISSUE` |
| blood flow rate to intestinal mucosa | `BLOOD_FLOW_RATE_TO_COMPARTMENT` |
| tubular flow fraction of the glomerular filtration rate | `BODILY_FLUID_FLOW_RATE_PARAMETER` |
| weight of a tissue in the arm vein sampling site | `PHYSIOLOGICAL_PARAMETER` |
| `CONCENTRATION_IN_HEART` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `CONCENTRATION_IN_PANCREAS` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `CONCENTRATION_IN_MUSCLE` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `CONCENTRATION_IN_BONE` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| `CONCENTRATION_IN_GONAD` | `CONCENTRATION_OF_COMPOUND_IN_COMPARTMENT` |
| amount of undissolved drug | `AMOUNT_IN_LUMEN` |
| dissolution process | none, the terms of other ontologies only |
| gastrointestinal transit process | none, the terms of other ontologies only |
| bone compartment | none, the terms of other ontologies only |
| portal vein blood compartment | none, the terms of other ontologies only |

PBPKO is further not registered at identifiers.org, so its urls do not
resolve there yet.
