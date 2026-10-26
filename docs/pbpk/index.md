# PBPK library

`pkpdsim.pbpk` builds the whole-body PBPK model of a substance from a short definition. A model package writes only what is specific to its substances: the substances with their mechanism flags, the chemistry of its organs (metabolism, transport), its choices (organ variants, routes, population, sex) and its PD components. The library adds a standard human physiology (ICRP 89), a standard organ for every organ, the alternative organ variants, the routes and sinks, the impaired organ function, the observers and the annotation of every element, and writes the SBML: a comp model with one submodel per organ and the flat model which is simulated.

The pages of this section are for modellers first: this page and the [tutorial](tutorial.md) show how a model is built and how the parts fit together; the pages of the parts say what each part does, which equations it generates and where its numbers come from.

## How a model is built

<div class="pbpk-figure" markdown="0">
--8<-- "pbpk/img/composition.svg"
</div>

1. The model package describes its substances (`Substance`) and their chemistry (`Metabolism`, `Transport`, `CustomReaction`), and makes its choices in one `PBPKDefinition`: the organ variants (`Variants`, the standard organ everywhere else), the routes (`Route`), the PD components (`PDComponent`), the sex and the population.
2. `build(definition, output_dir)` in the factory of the model validates the definition and composes the comp model: a top model with the physiology, the blood pools, the parameters of the substances, the sinks, the dose states and the observers (of the library and of the definition), the creators of the definition as model history, and one submodel per organ (`LI`, `KI`, `GU`, ...), the standard organ or a variant, which owns everything inside its boundary. Every port and every replacement is generated; the model package never writes a `ReplacedElement`.
3. `build` writes the comp model `<sid>.xml` and flattens it into `<sid>_flat.xml`, the model which `MODEL_PATH` points at and which every simulation runs. Simulation never resolves comp.
4. The experiments of the model set the covariates (`BW`, `HEIGHT`), the doses (`PODOSE_<s>`), the impairment (`CP_severity`, `f_renal_function`) and the activities of the proteins (`LI__f_cyp2c9`) as they do today, and read the observers (`Carm_<s>`, `Aurine_<s>`).

## A definition

A simple substance needs a short specification and its chemistry: an oral drug, bound to albumin, with one Kp, renal filtration and a hepatic CYP3A4 metabolism into a metabolite.

```python
from pkpdsim.pbpk import (
    MM,
    Binding,
    Excretion,
    Kp,
    Metabolism,
    Organ,
    PBPKDefinition,
    Route,
    Substance,
)

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
        ),
        Substance("api_m2", name="O-desmethyl apixaban", mr=445.5),
    ],
    chemistry=[
        Metabolism(
            "API2M2",
            Organ.LIVER,
            "api",
            products={"api_m2": 1.0},
            kinetics=MM(vmax=0.01, km=0.05),
            protein="CYP3A4",
        ),
    ],
    routes=[Route.iv("api"), Route.oral("api")],
)
definition.validate()
print([s.sid for s in definition.substances], definition.sex, definition.population)
```

```text
['api', 'api_m2'] male european_icrp
```

Everything else is a default: the standard organ in every organ (the segmented intestine for small and large intestine), the male ICRP 89 reference individual, the arm vein preset of PK-Sim. The [tutorial](tutorial.md) builds and simulates such a model step by step.

## How the parts fit

| part | what it is | page |
|---|---|---|
| physiology | the organs, their volumes, blood flows and blood (ICRP 89), the topology of the circulation, the scaling to the individual from `BW` and `HEIGHT`, the populations | [Physiology](physiology.md) |
| standard organ | plasma, blood cells, interstitial and intracellular space of every organ; only the parts a substance needs are built (collapsing) | [Standard organ](organ.md) |
| substances | the mechanism flags of a substance: protein binding, erythrocytes, tissue distribution, excretion, inhibition and induction | [Substances](substances.md) |
| chemistry | the processes of an organ, defined once per organ: metabolism, transport, custom reactions; the proteins with their activity factors | [Chemistry](chemistry.md) |
| organ variants | alternatives of an organ with the boundary of the standard organ: the zonated liver, the nephron kidney, the intestines | [Organ variants](variants.md) |
| intestine and routes | the lumen segments, dissolution, absorption, the mucosa, the bile and the enterohepatic circulation; the routes and the sinks | [Intestine and routes](intestine.md) |
| impairment | liver cirrhosis, chronic kidney disease and heart failure as parameter vectors on the physiology | [Impairment](impairment.md) |
| observers | the plasma, free and blood concentrations, the hepatic and the arm vein, the urinary excretion rate, the observers of a definition, the sampling sites of the fit mappings | [Observers](observers.md) |
| annotation | the SBO, PBPKO, UBERON, NCIT, UniProt and ChEBI terms of every element | [Annotation](annotation.md) |

The [migration guide](migration.md) explains the move of a hand-written model of pkdb_models onto the library, with the glimepiride model as example.

## What is generated

Every quantity which depends on parameters only (the physiology, the impairment, the compartments, `fu`, `BP`, the coefficients and exchanges of the organs) is a constant with an initial assignment: an experiment which changes `BW` before the simulation changes every organ volume and flow with it, because sbmlsim evaluates the initial assignments again after the changes of an experiment. A change during a simulation is not followed, unless the definition names the parameter in `PBPKDefinition.dynamic` (an impairment changed in the middle of a study, a PD output): then every quantity downstream of it is an assignment rule. Only what reads a species is always an assignment rule: the outflow concentrations of the organs, the urinary excretion rate and the observers; an observer of a definition must therefore read a state (see [Observers](observers.md#observers-of-a-definition)). A model loaded with plain roadrunner sets a covariate with `r.setValue("init(BW)", 100.0)` and resets, so that the initial assignments are evaluated again.

Sex is a choice of the build (`PBPKDefinition.sex`): it selects the male or female reference values when the model is built, and a study of both sexes uses one model per sex. A sex switch inside the SBML made the models load about ten times slower.

## Values and their sources

Every number which the library writes into a model is a `Value`: the number, its unit and its source, a key of `data/sources.toml` (citation and DOI), and the generated SBML parameter carries the citation in its notes. The data files of the library (`pkpdsim/pbpk/data/*.toml`) give every number with its source; a number which no publication gives is a choice of the library with the source `ModelChoice` and a note which gives the reasoning. The tables of these pages are generated from the data files.

A number of a definition is a plain float (the unit of the model is implied), a `Value`, or a `Formula` over ids of the top model, which ties it to other parameters (a metabolite with the Kp of its parent, see [Substances](substances.md#values-fits-and-ids)); a fitted value enters with its uncertainty from the record of its fit, `fitted(FITS.parameter("Kp_api", unit=U.dimensionless))`, and keeps the id which `FITS.parameter` and `check_model` resolve. The units are those of the models: time in min, amount in mmole, volume in litre, area in m2.

## The id scheme

One id scheme for every model, built by the functions of `pkpdsim.pbpk.ids`, so that organs, observers, experiments and fits agree on the ids. The ids of the hand-written models are kept where they fit the scheme (`Cve_<s>`, `Aurine_<s>`, `BW`, `f_renal_function`, the submodel prefixes `LI__`, `KI__`, `GU__`). A local id of an organ is prefixed by its submodel in the flat model (`LI__Cpl_api`), an id of the top model is not. A substance id must be unique over all models (a metabolite takes the prefix of its parent, `api_m2`), must not contain `__` and must not start with a token of the scheme (`pl_`, `bc_`, `free_`, `z1_`, ...).

<div class="pbpk-ids" markdown>

| id | what | made by |
|---|---|---|
| `<O>__Cpl_<s>` | a local id of an organ in the flat model | `ids.flat(o, ids.species(part, s))` |
| `Cpl_<s>` | the plasma concentration in an organ (local) | `ids.species(part, s)` |
| `Cti_<s>` | the merged tissue (interstitial and intracellular) of an organ | `ids.species(Part.TISSUE, s)` |
| `Cu_ti_<s>` | the unbound tissue concentration of an organ | `ids.unbound(Part.TISSUE, s)` |
| `C<pool>_<s>` | the plasma concentration of a blood pool (`ar`, `ve`, `po`) | `ids.pool_species(pool, s)` |
| `C<pool>_bc_<s>` | the blood cells of a pool (kinetic erythrocytes) | `ids.pool_species(pool, s, Part.BLOOD_CELLS)` |
| `C<pool>_free_<s>` | the unbound plasma concentration of a pool | `ids.free(pool, s)` |
| `C<pool>_blood_<s>` | the blood concentration of a pool | `ids.blood(pool, s)` |
| `Carm_<s>` | the arm vein plasma concentration, the default fit target | `ids.arm(s)` |
| `Carm_free_<s>` | the unbound arm vein concentration | `ids.arm_free(s)` |
| `Chv_<s>` | the hepatic vein plasma concentration | `ids.hepatic_vein(s)` |
| `Cout_<o>_<s>` | the plasma leaving an organ | `ids.organ_outflow(o, s)` |
| `Aurine_<s>` | the amount in urine | `ids.urine(s)` |
| `Rurine_<s>` | the urinary excretion rate (mmole/min), an output of the kidney | `ids.urine_rate(s)` |
| `Afeces_<s>` | the amount in feces | `ids.feces(s)` |
| `Cbile_<s>` | the concentration in bile | `ids.bile(s)` |
| `Mr_<s>` | the molecular weight | `ids.mr(s)` |
| `fu_ref_<s>` | the unbound fraction of the healthy individual (input) | `ids.fu_ref(s)` |
| `fu_<s>` | the unbound fraction, `fu_ref_<s>` changed by the impairment | `ids.fu(s)` |
| `Kp_<s>` | the tissue to plasma coefficient of the substance | `ids.kp(s)` |
| `Kp_<o>_<s>` | the coefficient of one organ (an override or a Kp method) | `ids.kp_organ(o, s)` |
| `ftissue_<s>` | the whole-body exchange between plasma and tissue | `ids.ftissue(s)` |
| `PS_<s>` | the exchange of an organ (l/min, local id) | `ids.ps(s)` |
| `kps_<s>` | the exchange of an organ per litre of organ (`ps`, 1/min, local id) | `ids.kps(s)` |
| `PS_ref_<s>` | the exchange of the whole healthy liver (`ps` with `ps_per_volume=False`, l/min, local id) | `ids.ps_ref(s)` |
| `BP_<s>` | the blood to plasma ratio | `ids.bp(s)` |
| `Kbc_<s>` | the blood cell to plasma coefficient | `ids.kbc(s)` |
| `f_reabsorbed_<s>` | the fraction of the filtrate which is reabsorbed | `ids.fraction_reabsorbed(s)` |
| `IVDOSE_<s>` | the intravenous dose state (mg) | `ids.iv_dose(s)` |
| `Ri_<s>` | the intravenous infusion (mg/min) | `ids.infusion_rate(s)` |
| `PODOSE_<s>` | the oral dose state of a tablet (mg) | `ids.po_dose(s)` |
| `SOLDOSE_<s>` | the oral dose state of a solution (mg) | `ids.sol_dose(s)` |
| `Ka_dis_<s>` | the first order dissolution of a tablet (1/hr) | `ids.dissolution(s)` |
| `ka_<s>` | the absorption rate constant of the small intestine (1/hr) | `ids.absorption(s)` |
| `Clumen_<seg>_<s>` | the dissolved substance in a lumen segment | `ids.lumen(s, seg)` |
| `Cmu_<seg>_<s>` | the mucosa (enterocytes) of a segment | `ids.mucosa(s, seg)` |
| `V<o>` | the volume of an organ (tissue and blood) | `ids.organ_volume(o)` |
| `Q<o>` | the blood flow of an organ (the hepatic artery for the liver) | `ids.organ_flow(o)` |
| `M<o>` | the tissue mass of an organ without blood | `ids.organ_mass(o)` |
| `f_vas_<o>` | the vascular fraction of the organ volume | `ids.vascular_fraction(o)` |
| `f_<protein>` | the activity of a protein in an organ (`LI__f_cyp3a4`) | `ids.activity(protein)` |
| `E_<protein>` | the turnover state of a protein in an organ | `ids.enzyme(protein)` |
| `Ki_<protein>_<s>` | the inhibition constant of a substance on a protein | `ids.ki(protein, s)` |
| `fimp_<o>_<protein>` | the factor of the impairment on a protein of an organ | `ids.impairment_activity(protein, o)` |
| `BW` | the body weight (kg) | `ids.BW` |
| `HEIGHT` | the height (cm) | `ids.HEIGHT` |
| `HCT` | the hematocrit | `ids.HCT` |
| `CO` | the cardiac output | `ids.CO` |
| `GFR` | the glomerular filtration rate | `ids.GFR` |
| `GFR_base` | the glomerular filtration rate at renal function 1 | `ids.GFR_BASE` |
| `f_renal_function` | the renal function, the GFR relative to healthy | `ids.F_RENAL_FUNCTION` |
| `f_cardiac_function` | the cardiac function, a factor on the cardiac output | `ids.F_CARDIAC_FUNCTION` |
| `CP_severity` | the severity of liver cirrhosis (0, Child-Pugh A 1, B 2, C 3) | `ids.CP_SEVERITY` |
| `HF_severity` | the severity of heart failure (0, NYHA I to IV 1 to 4) | `ids.HF_SEVERITY` |
| `f_shunt` | the portosystemic shunt of the liver | `ids.F_SHUNT` |
| `f_functional` | the functional fraction of the liver tissue | `ids.F_FUNCTIONAL` |
| `f_tubular` | the share of the tubular capacity of the kidney which follows the renal function (`KI__f_tubular`) | `ids.F_TUBULAR` |

</div>

A process of the chemistry keeps the id its definition gives (`API2M2`), its parameters are `<p>_Vmax`, `<p>_Km`, `<p>_k`, `<p>_n` in the organ (`LI__API2M2_Vmax` in the flat model). The table is checked against `pkpdsim.pbpk.ids` by `tests/docs/test_pbpk_ids.py`.
