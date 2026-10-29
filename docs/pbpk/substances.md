# Substances

A `Substance` is a short id, a name, a molecular weight and the mechanism flags which decide the equations built for it. A flag which is off builds nothing: a substance without binding has `fu = 1`, without erythrocytes no blood cell species, without distribution no tissue species. The flags are dataclasses of `pkpdsim.pbpk`.

| field | type | off (`None`) means | what it switches on |
|---|---|---|---|
| `sid` | `str` | | the id, unique over all models; a metabolite takes the prefix of its parent (`api_m2`) |
| `name`, `mr` | `str`, number | | the name and the molecular weight (g/mole), `Mr_<s>` |
| `annotations` | list of `(BQB, resource)` | | the ChEBI, PubChem and NCIT terms of the species of the substance |
| `binding` | `Binding(fu, protein, ckd=True)` | `fu = 1` | protein binding: `fu_ref_<s>`, and `fu_<s>` which follows the binding protein in impairment |
| `erythrocytes` | `Erythrocytes(mode, bp, kbc, k_exchange)` | no blood cells, `BP = 1 - HCT` | the distribution into the erythrocytes, in equilibrium or kinetic |
| `distribution` | `Kp(...)` or `KpMethod(...)` | vascular only | the distribution into the tissues: one Kp with overrides, or a Kp method |
| `excretion` | `Excretion(urine, fraction_reabsorbed)` | the default `Excretion()`: no filtration | glomerular filtration `fu*GFR` into the urine |
| `inhibits` | tuple of `Inhibition` | | the inhibition of proteins by the substance, see [Chemistry](chemistry.md) |
| `induces` | tuple of `Induction` | | the induction of proteins by the substance |
| `baseline` | number, `Value` or `Formula` (mM) | every species starts at 0 | an endogenous substance: every species starts at the equilibrium with the plasma concentration `baseline_<s>`, see [Endogenous baseline](#endogenous-baseline) |

The routes are declared on the definition (`PBPKDefinition.routes`), not on the substance; the biliary and the intestinal excretion are processes of the chemistry (an apical efflux of the liver into the bile, of the intestine into the lumen).

## Protein binding

`Binding(fu=0.13, protein="albumin")` gives the fraction unbound in plasma of the healthy individual, `fu_ref_<s>`, and the binding protein: `albumin`, `agp` (alpha-1-acid glycoprotein) or `unknown`. The unbound fraction which every rate reads, `fu_<s>`, is an initial assignment of `fu_ref_<s>` and the impairment: the binding protein falls in liver cirrhosis and the binding capacity in chronic kidney disease, so fu rises (see [Impairment](impairment.md)); `ckd=False` keeps the substance free of the CKD factor of the binding capacity, for a substance whose binding is measured unchanged in renal impairment (apixaban: 92-94 % bound in healthy and ESRD subjects), so that only the factor of its binding protein (cirrhosis) acts on fu; healthy, `fu_<s> = fu_ref_<s>`. Fit `fu_ref_<s>`, not `fu_<s>`.

## Erythrocytes

| `Erythrocytes(...)` | parts | blood flows carry | blood to plasma ratio |
|---|---|---|---|
| `None` | plasma only | `Q*(1 - HCT)*C_p` | `BP_<s> = 1 - HCT` |
| `Erythrocytes("equilibrium", bp=0.8)` | plasma only, capacity `V_vas*BP` | `Q*BP*C_p` | `BP_<s>` given |
| `Erythrocytes("equilibrium", kbc=1.5)` | plasma only, capacity `V_vas*BP` | `Q*BP*C_p` | `BP_<s> = 1 - HCT + HCT*Kbc_<s>` |
| `Erythrocytes("kinetic", kbc=1.5, k_exchange=0.1)` | plasma and blood cells `Cbc_<s>` | `Q*(1 - HCT)*C_p` and `Q*HCT*C_bc` | `BP_<s> = 1 - HCT + HCT*Kbc_<s>` at equilibrium |

Without erythrocytes the substance stays in the plasma: the blood concentration observers are `(1 - HCT)*C_p`, the plasma concentration is not diluted by a hand-written conversion. With erythrocytes in equilibrium every plasma species has the capacity of the blood (`V_vas*BP` instead of `V_vas*(1 - HCT)`), and with kinetic erythrocytes every vascular space has a blood cell species which exchanges with the plasma, `kex_bc_<s>*Vbc*(Cpl_<s> - Cbc_<s>/Kbc_<s>)`.

## Tissue distribution

`Kp(value, overrides={...}, interstitial=..., organs=..., mode=..., ftissue=..., ps={...})` is the whole-tissue (interstitial and intracellular) to plasma coefficient of the substance, `Kp_<s>`, with the coefficients of single organs, `Kp_<o>_<s>` (`overrides={Organ.ADIPOSE: 0.4}`). An organ of `overrides`, `organs`, `ps` and of `Variants(organs=...)` may also be given by its value or name (`"ad"`, `"adipose"`) and the sex of the definition by its value (`sex="female"`): the same holds for the `mode` of a Kp and for the `organ`, `direction`, `membrane` and `part` of a process of the chemistry. The definition holds the enum members, and anything else raises an error which names the field and the allowed values, as does a text which is none of the choices of `Binding.protein`, `Erythrocytes.mode`, `KpMethod.method`/`compound`, `Inhibition.kind` or `Route.kind`. The Kp means the same in every form of the organ, so a fitted Kp keeps its meaning when the mode changes. The keyword-only options:

- `organs`: the organs the substance enters (`frozenset({Organ.LIVER, Organ.KIDNEY})`), all if `None`; in the other organs it stays in the blood.
- `mode`: how plasma and tissue exchange (the collapsing of the [standard organ](organ.md#collapsing)):
    - `Mode.TISSUE` (the default): plasma and one tissue species with the finite exchange `PS*fu*(C_p - C_t/Kp)`;
    - `Mode.PERFUSION`: one species for the whole organ, flow-limited;
    - `Mode.PERMEABILITY`: interstitial and intracellular species with a permeability each.
- `ftissue`: the whole-body exchange `ftissue_<s>` (l/min) of the modes `TISSUE` and `PERMEABILITY`; an organ has the share of its tissue volume, `PS_<s> = ftissue_<s>*Vti_o/Vti_total_ref_<s>` (the reference tissue volume of the organs the substance enters), so the whole-body exchange scales with the tissue volume of the individual. It is needed only when an organ the substance enters has no `ps` of its own: a metabolite which enters only the liver, with the liver exchange `ps={Organ.LIVER: ...}`, has no `ftissue` (and no `Vti_total_ref_<s>`).
- `ps`: the exchange of single organs instead of the share of `ftissue`, e.g. a fast exchange of the liver for a highly bound drug, `ps={Organ.LIVER: 42.9}` (1/min per litre of liver, 42.9 * 2.33 l = 100 l/min for the whole reference liver, the value of the [tutorial](tutorial.md)). A value is per litre of organ (1/min, `kps_<s>`), so the organ has `PS_<s> = kps_<s>*V<o>` (of the liver times `f_functional`) and its exchange follows the organ of the individual, as the share of `ftissue` and a `FirstOrder` per volume do; `ps_per_volume=False` gives the exchange of the whole organ in l/min (`PS_<s>`), which does not scale with the organ volume; in the liver the value is that of the healthy liver (`PS_ref_<s>`) and the organ has `PS_<s> = PS_ref_<s>*f_functional`, as a `FirstOrder` of the whole organ. The validation refuses a `Value` of `ps` in the other unit.

!!! warning "ftissue is not the ftissue of the hand-written models"
    `ftissue_<s>` multiplies `fu*(C_p - C_t/Kp)`, the unbound plasma concentration and the tissue concentration as its plasma equivalent; in the hand-written models it multiplied the total concentrations. The id is the same, so the fit parameter set stays the same, but the fitted values differ and are refitted. For a highly bound substance (`fu` 0.005) the exchange of the liver from the share of `ftissue` can limit the metabolism; give the liver its own `ps`.

### Extracellular distribution

A substance which reaches the interstitial fluid of the tissues and not their cells (a contrast agent such as gadoxetate, a marker of the extracellular volume, a highly bound acid such as torasemide, which distributes with its albumin) has in every organ the Kp of the interstitial share of its whole tissue, `f_int*K_int_pls(fu)/(1 - f_vas)`: the interstitial fraction of the tissue without blood times the interstitial to plasma coefficient of its fu (Schmitt 2008). `Kp(interstitial=1.0, ftissue=...)` gives it to every organ the substance enters which has no override: `Kp_<o>_<s> = fKp_<s> * share`, the share a number of the build (`interstitial_share(organ, fu)`, the fractions of the reference individuals, the larger of both sexes, at the fu of the definition, which must be a number) and the scale `fKp_<s>` a number, a fitted `Value` (one fitted fraction of the interstitial fluid the substance reaches) or a `Formula` of a parameter of the definition. An override keeps its own value (an organ of the full form needs a Kp above its share, `check_kp`, so the scale 1 is refused there), and the value of the substance is not needed: a `Kp` without `value` has no `Kp_<s>`, which is also allowed when the overrides cover every organ the substance enters (else the validation names the organs without one). Without `interstitial` the model has no `fKp_<s>` and its SBML is unchanged.

```python
from pkpdsim.pbpk import Binding, Kp, Organ, Substance, interstitial_share

contrast = Substance(
    "ctr",
    name="an extracellular contrast agent",
    mr=725.7,
    binding=Binding(fu=0.9, protein="albumin"),
    distribution=Kp(interstitial=1.0, ftissue=1.0),
)
for organ in (Organ.BRAIN, Organ.MUSCLE, Organ.SKIN):
    print(organ.name, round(interstitial_share(organ, 0.9), 4))
```

```text
BRAIN 0.0045
MUSCLE 0.1347
SKIN 0.5386
```

`KpMethod("rodgers_rowland", logp=..., pka=(...), compound=...)` computes the coefficient of every organ when the model is built, from the lipophilicity, the ionisation and the tissue composition (Rodgers, Leahy & Rowland 2005, Rodgers & Rowland 2006, in the form of Yau 2020), and writes `Kp_<o>_<s> = fKp_<s>*value`: one factor `fKp_<s>` (`factor=`, default 1) scales all coefficients, so a fit can adjust them together. A base with a pKa of at least 7 binds to the acidic phospholipids, whose affinity follows from the erythrocytes, so it needs `Erythrocytes(bp=...)` or a `kbc`. The coefficients are numbers of the build: a change of `fu` by an impairment does not recompute them. The rat tissue composition is applied to the human tissues (Rodgers & Rowland 2007, Coutinho 2024).

```python
from pkpdsim.pbpk import Binding, Erythrocytes, KpMethod, Organ, Sex, Substance
from pkpdsim.pbpk.kp import kp_values

base = Substance(
    "bas",
    name="a lipophilic base",
    mr=350.0,
    binding=Binding(fu=0.2, protein="agp"),
    erythrocytes=Erythrocytes("equilibrium", bp=1.1),
    distribution=KpMethod(
        "rodgers_rowland", logp=3.0, pka=(9.0,), compound="base", ftissue=10.0
    ),
)
kp = kp_values(base, Sex.MALE)
for organ in (Organ.ADIPOSE, Organ.LIVER, Organ.MUSCLE, Organ.BRAIN):
    print(organ.name, round(kp[organ].value, 2))
```

```text
ADIPOSE 1.81
LIVER 15.94
MUSCLE 5.6
BRAIN 1.89
```

The full form (`Mode.PERMEABILITY`) derives the intracellular coefficient from the Kp and the interstitial one (Schmitt 2008). A Kp at or below the interstitial share has no positive intracellular coefficient: the build refuses it for the reference individual, and the floor of `K_cell_pls_<s>` keeps a fit or another individual from crossing it. A highly bound, hydrophilic acid can lie below the interstitial share in adipose; such a substance uses the merged tissue or the flow-limited organ.

## Excretion

`Excretion(urine=True)` filters the unbound plasma concentration in the kidney, `fu*GFR*C_p`, into `Aurine_<s>`; `fraction_reabsorbed` (`f_reabsorbed_<s>`) returns a fraction of the filtrate, so the renal clearance without transport is `fu*GFR*(1 - f_reabsorbed)`. Tubular secretion and transporter-mediated reabsorption are transport chemistry of the kidney (the standard kidney refuses an uptake from the cumulative urine; a reabsorption by a transporter needs the [nephron kidney](variants.md#nephron-kidney)). The fecal excretion is what leaves the last lumen segment of the intestine.

## Endogenous baseline

An endogenous substance (creatinine, creatine, glucose) is in the body before any dose. `baseline=0.05` gives its plasma concentration (mM) at the start, the parameter `baseline_<s>`, and every species of the substance starts at the equilibrium of its exchange with plasma: the plasma at the baseline, every other part of the blood and the tissue at the baseline times its coefficient to plasma, which sets the net exchange of the library to zero. Without a baseline every species starts at 0 and the model has no `baseline_<s>`.

| species | start | the exchange it balances |
|---|---|---|
| plasma of an organ `Cpl_<s>`, pools `Cve_<s>`, `Car_<s>`, `Cpo_<s>` | `baseline_<s>` | the blood flows carry the same concentration in and out |
| kinetic blood cells `Cbc_<s>`, `Cve_bc_<s>` (pools) | `Kbc_<s>*baseline_<s>` | `kex_bc*Vbc*(Cpl - Cbc/Kbc)` |
| merged tissue `Cti_<s>` | `Kp*baseline_<s>` (`Kp_<s>` or `Kp_<o>_<s>`) | `PS*fu*(Cpl - Cti/Kp)` |
| interstitium `Cis_<s>` | `K_int_pls_<s>*baseline_<s>` | `PS_vas*fu*(Cpl - Cis/K_int_pls)` |
| intracellular space `Cic_<s>` | `K_cell_pls_<s>*baseline_<s>` | `PS*(fu/K_int_pls*Cis - fu/K_cell_pls*Cic)` |
| flow-limited organ `Cor_<s>` | `K_or_<s>*baseline_<s>` | the organ releases `Cor/K_or` |
| intestinal mucosa `Cmu_<seg>_<s>` | `K_mu_<s>*baseline_<s>` | the mucosa releases `Cmu/K_mu` |
| lumen, tubule, urine, feces, bile, dose states | 0 | |

With erythrocytes in equilibrium the plasma species live in the capacity `V*BP` and start at the baseline as well; the zones of the zonated liver, the wall of the intestine and the nephron kidney start as the standard organ. The arm vein `Carm_<s>`, a mean of venous outflows, is then at the baseline. The starts are initial assignments of the species, which read the coefficients and `baseline_<s>`: an impairment or a covariate changed before the simulation (`preinit_changes`) moves `fu`, the coefficients and the volumes, and every species follows, so the body without production and elimination stays where it starts. A substance which is eliminated (filtration, metabolism) and has no production decays from its baseline; its production is chemistry of the model (a `CustomReaction`), which balances the elimination at the baseline.

```equations
# the start of an endogenous substance at the equilibrium with its baseline
Cpl_<s> = baseline_<s>
Cve_<s> = baseline_<s>
Cbc_<s> = Kbc_<s> * baseline_<s>
Cve_bc_<s> = Kbc_<s> * baseline_<s>
Cti_<s> = Kp_<s> * baseline_<s>
Cis_<s> = K_int_pls_<s> * baseline_<s>
Cic_<s> = K_cell_pls_<s> * baseline_<s>
Cor_<s> = K_or_<s> * baseline_<s>
```

A `Formula` baseline is the initial assignment of `baseline_<s>`, so a model can state it as a production over a clearance, e.g. `Formula("0.006 mmole_per_min / GFR", U.mM)` or, with the production as a [parameter of the definition](#parameters-of-the-definition) which a custom reaction reads as well, `Formula("CRTNPROD_k / GFR", U.mM)`: a change of the renal function before the simulation moves the baseline and with it every species. Besides the parameters and compartments of the top model it may read the observers of the definition, which are evaluated at the start (a clearance observer which reads a dynamic GFR); what it reads, directly or through the observers, must not change during the simulation (species, reactions, the observers of the library, the outputs of the organs, and a parameter which a PD output replaces with everything that reads it, such as `CO` under a PD of the cardiac function), since the species start at the baseline, and the build names every such id. With a dynamic parameter `baseline_<s>` may become an assignment rule (it reads the parameter); the species keep their initial assignments and start once. The validation refuses a baseline which is not positive and a `Value` in another unit than mM.

## Values, fits and ids

A number of a substance is a plain float in the units of the model or a `Value` with unit and source (a key of `data/sources.toml` of the library); a fitted number is `fitted(FITS.parameter(pid, unit=...))`, which carries the distrib uncertainty of the fit into the SBML parameter (`notes=` adds a remark to the notes of the parameter, e.g. the literature of a process):

```python
# not executed: FITS are the adopted fits of a model package
from pkpdsim.pbpk import Binding, Kp, Substance, U, fitted

glimepiride = Substance(
    "gli",
    name="glimepiride",
    mr=490.6,
    binding=Binding(fu=0.005, protein="albumin"),
    distribution=Kp(
        fitted(FITS.parameter("Kp_gli", unit=U.dimensionless)),
        ftissue=fitted(FITS.parameter("ftissue_gli", unit=U.l_per_min)),
    ),
)
```

A value which is tied to other parameters is a `Formula(expression, unit, notes=None)`: an SBML L3 formula over ids of the top model, accepted wherever a number or a `Value` is (the Kp and its overrides, `ftissue`, `ps`, the binding, the kinetic constants of the chemistry, the routes). The hand-written glimepiride model gave its metabolites the distribution of glimepiride by assignment rules; with the library the metabolite reads the parameters of the parent, so the fit keeps the parameters of the parent only:

```python
from pkpdsim.pbpk import Binding, Formula, Kp, Organ, Substance, U

m1 = Substance(
    "gli_m1",
    name="hydroxy glimepiride (M1)",
    mr=506.6,
    binding=Binding(fu=0.48),
    distribution=Kp(
        Formula("Kp_gli", U.dimensionless, notes="the Kp of glimepiride"),
        ftissue=Formula("ftissue_gli * fu_ref_gli / fu_ref_gli_m1", U.l_per_min),
        overrides={Organ.LIVER: 1.0},
        ps={Organ.LIVER: 0.08},
    ),
)
print(m1.distribution.value.expression)
```

```text
Kp_gli
```

The parameter of a formula (`Kp_gli_m1`) has the formula as its initial assignment and names it in its notes, as the covariates do: a change of a parameter it reads before the simulation (`preinit_changes`, `init(Kp_gli)`) is followed, and an id of `PBPKDefinition.dynamic` makes it an assignment rule. The ids of the formula are constant parameters and compartments of the top model (the parameters of the substances such as `Kp_<s>`, `fu_ref_<s>` and `ftissue_<s>`, the covariates, the physiology, the [parameters of the definition](#parameters-of-the-definition)), also for a value of an organ (`ps`, a kinetic constant): the organ reads each of them through an input port under the same id, so in the flat model the formula of `LI__kps_<s>` reads `BW`, not an id of the liver. A local id of an organ is no id of a formula, also not that of an input whose local id differs from the id of the top model (the liver reads `fimp_li_cyp3a4` of the top model as `fimp_cyp3a4`; a formula reads `fimp_li_cyp3a4`). An organ reads the parameters and compartments which the top model has before the organs are built. The build names every id of a formula which the top model does not have, and refuses species and observers (and every element with a rule, such as the outflow of an organ): they change during the simulation, and a formula, evaluated before it, would freeze them at their initial values. The unit is the declared unit of the parameter, it is not inferred from the formula. A formula has no number before the simulation: the validation skips the positivity and interval checks of such a value, the check that a Kp leaves a positive `K_cell_pls` is skipped (the model bounds `K_cell_pls` below by 1e-6), and the Kp method, which computes with `fu`, `bp` and `kbc`, refuses a formula there.

### Parameters of the definition

A constant of the model which no field of the library holds is a parameter of the definition, `ModelParameter(sid, value, unit, name=None, notes=None)` in `PBPKDefinition.parameters`: a constant parameter of the top model under its id, e.g. the production rate of creatinine, which its baseline (production over clearance) and the production in the muscle (a custom reaction) share. The value is a number in `unit`, a `Value` (a fitted value, `fitted(FITS.parameter("CRTNPROD_k", unit=U.mmole_per_min))`, gets the provenance in the notes and the uncertainty of the fit as the `distrib` of the parameter, as every fitted value of the library, so `check_model` checks it against the record) or a `Formula` over ids of the top model (its initial assignment, which may read another parameter of the definition); a `Value` or a `Formula` is in `unit`. The formulas of the definition (a baseline, a Kp, a kinetic constant) and its observers read the parameter like any parameter of the top model, and so does the formula of a custom reaction: the organ reads it through an input port under the same id (see [Chemistry](chemistry.md#custom-reactions)), so in the flat model the reaction reads `CRTNPROD_k`, not an id of the muscle. A constant which only one custom reaction reads (a fitted constant of its mechanism) is a parameter of the reaction instead, a local id of its organ (`MU__CRIM_Km`, see [Chemistry](chemistry.md#custom-reactions)).

```python
from pkpdsim.pbpk import (
    CustomReaction,
    Excretion,
    Formula,
    Kp,
    ModelParameter,
    Organ,
    Part,
    PBPKDefinition,
    Substance,
    U,
    compose,
)

creatinine = PBPKDefinition(
    sid="creatinine_body",
    substances=[
        Substance(
            "crtn",
            name="creatinine",
            mr=113.1,
            distribution=Kp(1.0, ftissue=1.0),
            excretion=Excretion(urine=True),
            # production over filtration: the steady state of the arterial plasma
            baseline=Formula("CRTNPROD_k / GFR", U.mM),
        )
    ],
    chemistry=[
        CustomReaction(
            "CRTNPROD",
            Organ.MUSCLE,
            Part.TISSUE,
            "-> Cti_crtn",
            "CRTNPROD_k",
            parameters=[],
        )
    ],
    parameters=[
        ModelParameter(
            "CRTNPROD_k",
            0.0074,
            U.mmole_per_min,
            name="creatinine production",
            notes="About 1.2 g/day.",
        )
    ],
)
# build composes the model and checks every id against it; validate checks
# the definition alone (no clash with an id of a substance, sink or organ)
model = compose(creatinine)
print([p.sid for p in creatinine.parameters])
```

```text
['CRTNPROD_k']
```

The example starts every plasma species at `CRTNPROD_k / GFR`, the arterial plasma concentration at the steady state (the filtration clears `GFR*Car`, see [Organ](organ.md#kidney)). The body is not stationary at this start: the filtration takes its amount out of the plasma at once, while the production of the muscle reaches the plasma only as the tissue of the muscle rises. In a build of the example (`tests/docs/test_substances_example.py` checks these numbers) the arterial plasma drops about 17 % below the baseline within the first 20 minutes and is back within 0.2 % after two days; the arm vein plasma, which mixes the venous outflow of the muscle, settles about 5.5 % above the baseline. A model which needs a stationary start states the steady state of the producing organ as the [starts of its species](#starts-of-the-species).

The parameters are constants: a change before the simulation (`preinit_changes`) moves every formula which reads them, and an id of `PBPKDefinition.dynamic` makes them variable. The validation refuses an id which is no SBML SId, contains `__`, is given twice or is an id of the physiology or an organ submodel, and a `Value` or `Formula` in another unit than `unit`; the build (`compose`) refuses an id which the body has (a parameter of a substance such as `fu_<s>` or `Kp_<s>`, a sink, an observer of the library, a PD component), the id of an observer of the definition and a local id of an organ which is no input under the same id (`Vti`, `<process>_Vmax`): a custom reaction there would read the element of the organ instead. Formulas which read each other in a cycle are refused, naming them. Without parameters the SBML is that of a definition without the field.

### Starts of the species

A production or an elimination of the chemistry moves the organs where it acts away from the equilibrium with the baseline: at the steady state the plasma leaving the muscle carries the production, and the tissue of the muscle holds more than `Kp` times its plasma, by the production over the exchange. `PBPKDefinition.starts` gives such a species its own start: the flat id of the species (`MU__Cpl_crtn`, a pool `Cve_crtn`, a species of a custom reaction such as a pool of the muscle) to its initial concentration in mM, a number, a `Value` or a `Formula` over the constants of the top model, which replaces the start at the equilibrium. The organ reads the ids of a formula through input ports under the same ids, as for a `Formula` of the organ, and the start is an initial assignment of the species, so it follows a change of the individual before the simulation. The definition keeps the starts as a read-only mapping, as it keeps the values per substance of the intestine.

The steady state of the organs of the library is a closed form at a given inflow, and `pkpdsim.pbpk.steady_state` returns it as formulas over the ids of the top model, the ids the equations of the organs read, so a model does not restate the equations of the library:

| function | the steady state |
|---|---|
| `plasma_flow(s, organ, variants=)` | the flow which carries the substance, `Q<o>*(1 - HCT)` (`Q<o>*BP_<s>` with erythrocytes in equilibrium; the liver `(1 - f_shunt)*Qh*...`) |
| `tissue_volume(organ, variants=)`, `tissue_exchange(s, organ, variants=)` | the merged tissue `Vti` and its exchange `PS*fu`, the share of `ftissue_<s>` of the organ |
| `organ_plasma(s, organ, arterial=, rate=, clearance=, variants=)` | the plasma leaving a standard organ with a net source `rate` (mmole/min) and a clearance of its plasma, `(F*Car + rate)/(F + clearance)` |
| `organ_tissue(s, organ, plasma, rate=, variants=)` | the merged tissue with a source of the tissue, `Kp*(Cpl + rate/(PS*fu))` |
| `renal_filtration(s, variants=)`, `renal_clearance(s, clearance=, variants=)` | the net filtrate `G = fu*GFR*(1 - f_reabsorbed)` and the renal clearance of the arterial plasma with a secretion (glomerulus and tubule in series) |
| `kidney_plasma(s, arterial=, eliminated=, rate=, clearance=, secretion=, variants=)` | the kidney plasma `C` of `R*(Car - C) + rate = clearance*C + S(C)` (`R = F - G`), at a given arterial plasma or at a given net excretion into the urine (`eliminated`, the production of an endogenous substance which only the kidney excretes); a saturable secretion `secretion=Saturable(vmax, km)`, `S(C) = vmax*C/(km + C)`, makes `C` the positive root of a quadratic |
| `kidney_arterial(s, plasma, eliminated=, rate=, clearance=, secretion=, variants=)` | the arterial plasma of a kidney plasma, `C + (clearance*C + S(C) - rate)/R`; with the renal elimination `(E + R*C)/F`, which reads the kidney plasma once (the baseline of the elimination, read by every start of the substance) |
| `liver_plasma(s, arterial=, portal=, rate=, clearance=, variants=)` | the liver plasma, fed by the hepatic artery and the portal vein past the shunt, with a source and a clearance of the liver |
| `arm_vein(s, preset, outflows, arterial=)` | the arm vein of a preset, the mix of the outflows at the steady state |

The arterial plasma is the baseline by default. The functions take the variants of the definition (`variants=`, required), so a model with a variant of an organ never gets the formula of the standard organ silently (the flows of the variants are those of the physiology). An endogenous substance produced at a rate `P` and excreted only by the kidney with a saturable secretion (creatinine) has the kidney plasma `kidney_plasma(s, eliminated="P", secretion=Saturable("V", "km"), variants=...)` (a parameter of the definition, e.g. `C_ki_0`) and the baseline `kidney_arterial(s, "C_ki_0", eliminated="P", variants=...)`. The tests of the library build every function into the starts of a model with a production of the muscle, a supply of the liver plasma and tissue, a source and a secretion of the kidney and check that every species stays at its start, for another individual, a portosystemic shunt with a reduced functional liver and with erythrocytes in equilibrium, and in the zonated liver (its zones at the inflow), and the kidney with a saturable secretion of a production of the muscle (every species within 1e-9 of its start, `tests/pbpk/test_steady_state_saturable.py`); where no closed form applies a function raises `NotImplementedError` with the reason: kinetic erythrocytes, the full form and the flow-limited organ, an exchange of the organ of its own (`ps`), the intestine (it refuses custom reactions, its parts start at the equilibrium), the nephron kidney (its tubule starts empty) and a source in the zonated liver. With the steady state of the muscle the example starts stationary:

```python
import dataclasses

from pkpdsim.pbpk import steady_state

crtn = creatinine.substances[0]
variants = creatinine.variants
muscle = steady_state.organ_plasma(
    crtn, Organ.MUSCLE, rate="CRTNPROD_k", variants=variants
)
stationary = dataclasses.replace(
    creatinine,
    starts={
        "MU__Cpl_crtn": Formula(muscle, U.mM),
        "MU__Cti_crtn": Formula(
            steady_state.organ_tissue(
                crtn, Organ.MUSCLE, muscle, rate="CRTNPROD_k", variants=variants
            ),
            U.mM,
        ),
    },
    deltas=("Carm_crtn",),
)
model = compose(stationary)
print(sorted(stationary.starts))
```

```text
['MU__Cpl_crtn', 'MU__Cti_crtn']
```

In the build of this definition every species stays within 1e-6 of its start, for the reference individual as for another body weight or renal function (`tests/docs/test_substances_example.py`); the arm vein starts where the example without starts settled, 5.5 % above the baseline. A start which the library cannot know (a pool of a custom reaction, the steady state of a nonlinear process) is written by the model; a stationarity test of the model (every species over a day, several individuals) checks the starts. The validation refuses an id which is no SBML SId, a negative number and a `Value` or `Formula` in another unit than mM; the build refuses a start which names no species of the model (all of them at once) and a formula which reads an id which is no constant of the top model (a species or an observer would be read before the species start). Without starts the SBML is that of a definition without the field.

`PBPKDefinition.deltas` names observers (`Carm_crtn`, an observer of the definition, or an output of an organ: `Rurine_<s>`, a readout of the intestine, an exposed element of an organ) whose change over their value at the start a model reads, the increase of an endogenous substance over its value before a dose: each gets the constant `<id>_0`, the initial assignment of its value at the start (after the changes of the individual before the simulation), and the observer `<id>_delta = <id> - <id>_0`, which is 0 at the start of every simulation for any arm vein preset. The build refuses a delta which is no observer of the top model and an id `<id>_0` or `<id>_delta` which the model has; without deltas the SBML is unchanged.

The validation (`PBPKDefinition.validate()`, called by `build`) refuses an id which is no SBML SId, contains `__`, starts with a token of the id scheme or is an id of the physiology, a duplicate id, a fraction outside its interval, a non-positive constant, a mode which needs `ftissue` without it (an entered organ without `ps`, named in the message), kinetic erythrocytes with a flow-limited distribution, and names every problem at once.
