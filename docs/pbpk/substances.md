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

`Kp(value, overrides={...}, organs=..., mode=..., ftissue=..., ps={...})` is the whole-tissue (interstitial and intracellular) to plasma coefficient of the substance, `Kp_<s>`, with the coefficients of single organs, `Kp_<o>_<s>` (`overrides={Organ.ADIPOSE: 0.4}`). An organ of `overrides`, `organs`, `ps` and of `Variants(organs=...)` may also be given by its value or name (`"ad"`, `"adipose"`) and the sex of the definition by its value (`sex="female"`): the same holds for the `mode` of a Kp and for the `organ`, `direction`, `membrane` and `part` of a process of the chemistry. The definition holds the enum members, and anything else raises an error which names the field and the allowed values, as does a text which is none of the choices of `Binding.protein`, `Erythrocytes.mode`, `KpMethod.method`/`compound`, `Inhibition.kind` or `Route.kind`. The Kp means the same in every form of the organ, so a fitted Kp keeps its meaning when the mode changes. The keyword-only options:

- `organs`: the organs the substance enters (`frozenset({Organ.LIVER, Organ.KIDNEY})`), all if `None`; in the other organs it stays in the blood.
- `mode`: how plasma and tissue exchange (the collapsing of the [standard organ](organ.md#collapsing)):
    - `Mode.TISSUE` (the default): plasma and one tissue species with the finite exchange `PS*fu*(C_p - C_t/Kp)`;
    - `Mode.PERFUSION`: one species for the whole organ, flow-limited;
    - `Mode.PERMEABILITY`: interstitial and intracellular species with a permeability each.
- `ftissue`: the whole-body exchange `ftissue_<s>` (l/min) of the modes `TISSUE` and `PERMEABILITY`; an organ has the share of its tissue volume, `PS_<s> = ftissue_<s>*Vti_o/Vti_total_ref_<s>` (the reference tissue volume of the organs the substance enters), so the whole-body exchange scales with the tissue volume of the individual. It is needed only when an organ the substance enters has no `ps` of its own: a metabolite which enters only the liver, with the liver exchange `ps={Organ.LIVER: ...}`, has no `ftissue` (and no `Vti_total_ref_<s>`).
- `ps`: the exchange of single organs instead of the share of `ftissue`, e.g. a fast exchange of the liver for a highly bound drug, `ps={Organ.LIVER: 100.0}`. A value is per litre of organ (1/min, `kps_<s>`), so the organ has `PS_<s> = kps_<s>*V<o>` (of the liver times `f_functional`) and its exchange follows the organ of the individual, as the share of `ftissue` and a `FirstOrder` per volume do; `ps_per_volume=False` gives the exchange of the whole organ in l/min (`PS_<s>`), which does not scale with the organ volume; in the liver the value is that of the healthy liver (`PS_ref_<s>`) and the organ has `PS_<s> = PS_ref_<s>*f_functional`, as a `FirstOrder` of the whole organ. The validation refuses a `Value` of `ps` in the other unit.

!!! warning "ftissue is not the ftissue of the hand-written models"
    `ftissue_<s>` multiplies `fu*(C_p - C_t/Kp)`, the unbound plasma concentration and the tissue concentration as its plasma equivalent; in the hand-written models it multiplied the total concentrations. The id is the same, so the fit parameter set stays the same, but the fitted values differ and are refitted. For a highly bound substance (`fu` 0.005) the exchange of the liver from the share of `ftissue` can limit the metabolism; give the liver its own `ps`.

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

A value which is tied to other parameters is a `Formula(expression, unit, notes="")`: an SBML L3 formula over ids of the top model, accepted wherever a number or a `Value` is (the Kp and its overrides, `ftissue`, `ps`, the binding, the kinetic constants of the chemistry, the routes). The hand-written glimepiride model gave its metabolites the distribution of glimepiride by assignment rules; with the library the metabolite reads the parameters of the parent, so the fit keeps the parameters of the parent only:

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

The parameter of a formula (`Kp_gli_m1`) has the formula as its initial assignment and names it in its notes, as the covariates do: a change of a parameter it reads before the simulation (`preinit_changes`, `init(Kp_gli)`) is followed, and an id of `PBPKDefinition.dynamic` makes it an assignment rule. The ids of the formula are constant parameters and compartments of the top model (the parameters of the substances such as `Kp_<s>`, `fu_ref_<s>` and `ftissue_<s>`, the covariates, the physiology), also for a value of an organ (`ps`, a kinetic constant): the organ reads each of them through an input port under the same id, so in the flat model the formula of `LI__kps_<s>` reads `BW`, not an id of the liver. A local id of an organ is no id of a formula, also not that of an input whose local id differs from the id of the top model (the liver reads `fimp_li_cyp3a4` of the top model as `fimp_cyp3a4`; a formula reads `fimp_li_cyp3a4`). An organ reads the parameters and compartments which the top model has before the organs are built. The build names every id of a formula which the top model does not have, and refuses species and observers (and every element with a rule, such as the outflow of an organ): they change during the simulation, and a formula, evaluated before it, would freeze them at their initial values. The unit is the declared unit of the parameter, it is not inferred from the formula. A formula has no number before the simulation: the validation skips the positivity and interval checks of such a value, the check that a Kp leaves a positive `K_cell_pls` is skipped (the model bounds `K_cell_pls` below by 1e-6), and the Kp method, which computes with `fu`, `bp` and `kbc`, refuses a formula there.

The validation (`PBPKDefinition.validate()`, called by `build`) refuses an id which is no SBML SId, contains `__`, starts with a token of the id scheme or is an id of the physiology, a duplicate id, a fraction outside its interval, a non-positive constant, a mode which needs `ftissue` without it (an entered organ without `ps`, named in the message), kinetic erythrocytes with a flow-limited distribution, and names every problem at once.
