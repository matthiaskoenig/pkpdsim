# Context of use and credibility

The fit metrics of a model pool its validation data: one GMFE and one percentage within 2-fold for doses between the fitted doses, interactions, held-out subjects and processes the model is known not to describe. A number like that answers no question. What the model is credible for is not written down anywhere. This page describes how a model declares what it is used for, which evidence supports each use and when the use is accepted, and how `run_credibility` assesses the declaration against the adopted fits of the model.

## The frame

The frame is that of ASME V&V 40 (2018), which the FDA guidance on the credibility of computational modeling (2023) and ICH M15 on model-informed drug development (Step 4, 2026) took over, and which the EMA guideline on the reporting of PBPK modelling and simulation (2018) uses for the qualification of a platform. For PBPK models it is worked out by Kuemmel et al. (2020) and Musuamba et al. (2021), and Erdemir et al. (2020) give the rules of credible practice. The default acceptance uses fold criteria on the predictions of the validation data; for ratios between arms, e.g. impaired over healthy, a fixed fold is too lenient near 1, and Guest et al. (2011) define the limits which step 3 uses.

- The **question of interest** is what a decision needs to know, e.g. "Does the model predict the plasma concentration of the drug in renal impairment?".
- The **context of use** (COU) states the role of the model in answering the question: its quantities of interest and the domain it is claimed for.
- The **model risk** follows from the influence of the model on the decision (is it one argument among others or the only one) and the consequence of a wrong decision.
- The **evidence** is the data the model is compared with, each part with a role: calibration (the data the model was fitted to), validation (data held out of the fit) or limitation (data the model is known not to describe).
- The **assessment** decides per COU whether the validation evidence meets the acceptance, which part of the claimed domain the evidence covers, and which credibility factors of V&V 40 are assessed.

Credibility is per context of use. There is no score of a model.

## The ladder of contexts of use

A COU has one level, `CouLevel`. A model typically has two to five COUs.

| level | name | meaning | typical evidence |
|---|---|---|---|
| `C1` | describe | the model describes the data it integrates | calibration evidence |
| `C2` | interpolate | doses, routes, regimens and formulations inside the range of the data | held-out doses, subjects, studies |
| `C3` | extrapolate with data | special populations or conditions for which some data of the drug exists (renal or hepatic impairment, genotypes, food, interactions) | validation data of the population |
| `C4` | extrapolate without data | populations or conditions without data of the drug | the qualification of the platform across drugs, uncertainty (later steps) |
| `C5` | clinical endpoint | pharmacodynamics and clinical outcomes (blood pressure, glucose, coagulation) | PD validation data |

`CouLevel.label` gives the name of a level (`describe`, `interpolate`, ...). C1 is the only level decided on calibration evidence, since its question is whether the model describes the data it was fitted to. Every other level is decided on validation evidence only.

## Model risk, the default acceptance and the minimum of points

`Influence` and `Consequence` of a COU are `Level.LOW`, `Level.MEDIUM` or `Level.HIGH`; `ContextOfUse.risk` is the lookup `RISK_TABLE`:

| influence \ consequence | low | medium | high |
|---|---|---|---|
| low | low | low | medium |
| medium | low | medium | high |
| high | medium | high | high |

Influence is low when the model is one argument among others (a paper which also shows the data), medium when the prediction of the model is the main argument, high when it stands alone (no data of the condition, C4). Consequence is low for an explanation, medium for a statement a clinician or another modeller may act on (a dose adjustment) and high only with a rationale.

A COU which declares no `acceptance` takes the default of its risk (`DEFAULT_ACCEPTANCE`, `default_acceptance(risk)`), all rules of which must be met on the validation evidence:

| risk | rules |
|---|---|
| low | `within_2fold` >= 70 % |
| medium | `within_2fold` >= 80 % and GMFE <= 1.5 |
| high | `within_2fold` >= 90 %, `within_1.5fold` >= 70 % and GMFE <= 1.3 |

The thresholds are a starting point of this repository and not taken from a guideline, which leaves them to the COU. A COU overrides them with its own rules: `Acceptance(criterion, minimum)` (the percentage of the points a criterion of `sbmlsim.fit.evaluation`, e.g. `Fold(2.0)`, accepts is at least `minimum`), `MaxGmfe(maximum)` and `AfeRange(low, high)` (the bias).

A percentage of a few points says little, so a verdict also needs a minimum of evaluable decisive points. A COU which declares no `min_points` takes the default of its risk (`DEFAULT_MIN_POINTS`, `ContextOfUse.minimum_points`):

| risk | minimum of evaluable points |
|---|---|
| low | 5 |
| medium | 10 |
| high | 20 |

A COU overrides it with `min_points` (at least 1), e.g. for a pharmacodynamic endpoint of which a study reports a few values only; the rationale of the COU says why.

## Declaring the contexts of use

A model package gets a module `credibility.py` with the object `CREDIBILITY`, a `Credibility`, and references it in its model object like the experiments and the fit definitions, `credibility=".credibility:CREDIBILITY"` (see [Models and their commands](models.md)); it is resolved when a command uses it. A declaration is plain python, a tree of frozen dataclasses. This is the declaration of the test model of pkpdsim, a decay with an intravenous dose, which has calibration data, held-out doses, a renal impairment study and data the model does not describe:

```python
# not executed (the example is the module of the test package tests/pkpdsim_credibility)
from pkpdsim.credibility import (
    ContextOfUse,
    CouLevel,
    Credibility,
    Domain,
    Evidence,
    EvidenceRole,
    Level,
    QuantityOfInterest,
    Unclaimed,
)
from pkpdsim.fit import metadata_filter, named_filter
from pkpdsim.metadata import Health, Route
from pkpdsim_credibility.fitting import DOSE, LIMITATION, POPULATION

_A = (QuantityOfInterest("[A]", "the concentration of A"),)

CREDIBILITY = Credibility(
    contexts=(
        ContextOfUse(
            id="describe",
            title="Describe the decay",
            question="Does the model describe the intravenous decay of A?",
            level=CouLevel.C1,
            influence=Level.LOW,
            consequence=Level.LOW,
            quantities=_A,
            domain=Domain(),
            evidence=(
                Evidence(
                    "intravenous data",
                    EvidenceRole.CALIBRATION,
                    metadata_filter(route={Route.IV}),
                ),
            ),
        ),
        ContextOfUse(
            id="dose",
            title="Decay after a dose",
            question="Does the model predict the decay after a dose of the data?",
            level=CouLevel.C2,
            influence=Level.MEDIUM,
            consequence=Level.MEDIUM,
            quantities=_A,
            domain=Domain(
                metadata={"route": (Route.IV,), "health": (Health.HEALTHY,)},
                statement="the single doses of the data",
            ),
            evidence=(
                Evidence("doses", EvidenceRole.VALIDATION, named_filter(DOSE)),
                Evidence(
                    "other healthy data",
                    EvidenceRole.VALIDATION,
                    Unclaimed(
                        within=metadata_filter(
                            route={Route.IV}, health={Health.HEALTHY}
                        )
                    ),
                ),
            ),
        ),
        ContextOfUse(
            id="renal",
            title="Decay in renal impairment",
            question="Does the model predict the decay in renal impairment?",
            level=CouLevel.C3,
            influence=Level.MEDIUM,
            consequence=Level.HIGH,
            quantities=(
                *_A,
                QuantityOfInterest(
                    "A_ratio",
                    "the concentration in renal impairment over the healthy one",
                    ratio=True,
                ),
            ),
            domain=Domain(
                metadata={"health": (Health.HEALTHY, Health.RENAL_IMPAIRMENT)}
            ),
            evidence=(
                Evidence(
                    "renal impairment",
                    EvidenceRole.VALIDATION,
                    named_filter(POPULATION),
                ),
            ),
        ),
    ),
    limitations=(
        Evidence(
            "scaled data",
            EvidenceRole.LIMITATION,
            named_filter(LIMITATION),
            "data at three times the decay which the model does not describe",
        ),
    ),
)
```

- `ContextOfUse`: `id` (a slug matching `^[a-z][a-z0-9_]*$`, unique in the model), `title`, `question` (one sentence ending with `?`), `level`, `influence`, `consequence`, `quantities`, `domain`, `evidence`, and optional `acceptance` (empty means the default of the risk), `limitations` (evidence of role `LIMITATION` of this COU, subtracted like the limitations of the model), `rationale` (why this level, influence and consequence), `extrapolation` (see the verdicts) and `min_points` (`None` means the default of the risk). `ContextOfUse.risk`, `ContextOfUse.minimum_points` and `ContextOfUse.decisive_role` (calibration for C1, else validation) are derived.
- `QuantityOfInterest(id, description)`: an observable of the fit mappings (`Cve_hctz`) or the observable of a value per simulation (`pk.cmax`). The observable column of the data holds the selection, `[Cve_hctz]` for a concentration, so a quantity is matched after the brackets are removed (`[A]` and `A` are the same quantity); a quantity which matches no observable of the fits is an error. A quantity restricts the evidence of its COU to its observables. `ratio` marks a ratio between arms (impaired over healthy): `ratio=True` declares it without a pairing, it is matched against no observable, restricts no evidence and is shown as "ratio, no pairing declared"; a `Ratio` pairs its arms and is assessed, see [Ratios](#ratios). The absolute quantities decide the verdict, a ratio adds to it, so a COU has at least one quantity which is no ratio.
- `Domain`: `metadata`, the claimed values by field of the metadata (`route`, `health`, `genotype`, `phenotype`, `coadministration`, `dosing`, `fasting`, `application_form`, `tissue`; members of the enums of `pkpdsim.metadata`), and `statement`, a sentence for what the metadata cannot say (dose range, duration). The dose is not a field of the metadata, so a dose range is a statement only. The coverage is assessed on the metadata; the statement is shown, marked "statement, not assessed". A COU above C1 claims at least one value of the metadata, without one its domain would be covered by any data.
- `Evidence(name, role, selection, description, fit=None)`: the selection is a filter of fit mappings (`named_filter` of a `NamedMappings`, the mapping keys by study id, or `metadata_filter`, both of `pkpdsim.fit`) or `Unclaimed`; `fit` restricts the evidence to one subset of the adopted fits (`"PK"` or `"PD"`), and a subset the adopted fits do not have is an error naming the subsets they have.
- `Credibility(contexts, limitations, notes)`: the limitations are evidence of the whole model. A limitation, of the model or of a COU, is data the model does not describe and is subtracted from the evidence of every COU.

## The checks of a declaration

A declaration is validated when it is built (`DeclarationError` lists every problem):

- the ids match `^[a-z][a-z0-9_]*$` and are unique;
- a question ends with `?`;
- a COU has at least one quantity of interest which is no ratio;
- a COU of level C2 to C5 claims at least one value of the metadata in its domain, a C1 COU may claim none;
- the evidence of `limitations`, of the model and of a COU, has the role `LIMITATION`;
- a C1 COU has calibration evidence, a C4 COU has an `extrapolation` text (it has no validation data of its condition by definition), every other COU has validation evidence or an `extrapolation` text;
- `Unclaimed` is validation evidence and at most once per COU;
- a declared `min_points` is at least 1;
- the claimed metadata values are members of the enums.

 The assessment checks the declaration against the data on top of that, and any failure is an error and not an empty result:

- A quantity which matches no observable.
- An evidence which selects no point, named with its COU, evidence and the fits searched.
- Validation evidence on a mapping which is training data in any adopted fit. Circular evidence is a hard error and never a warning: a mapping fitted in one fit counts as fitted in every fit. The same holds for a limitation on a fitted mapping.
- Kinds of mappings which differ from the recorded ones (see below).

`check_declaration(model)` returns every problem of a declaration and of the records as messages, from the fit definitions and the records, without a simulation; the tests of a model repository use it.

## The evidence

| role | what it selects | what the assessment does with it |
|---|---|---|
| `CALIBRATION` | the training points of the selection | decides a C1 COU; for every other level it is reported and never decides |
| `VALIDATION` | held-out points | decides the acceptance of every level but C1; a training point is an error, an outlier is dropped and counted |
| `LIMITATION` | data the model does not describe | subtracted from every other evidence of every COU, reported apart, never in an acceptance |

A limitation stays a validation mapping for the fit (the kinds of sbmlsim are not changed); for the credibility it is a limitation. The points a limitation removes from an evidence are counted per evidence and shown in the report; a calibration keeps only training points, so it counts only the removed points it would have kept.

`Unclaimed(within=...)` is for validation data which is computed instead of named, e.g. held-out subjects: it takes the validation points which no other evidence of the model selects (the selected points of the other evidence, of any COU) and no limitation, restricted to the mappings `within` accepts. It is allowed once per COU and only as validation evidence. Without `within` it would take every validation point nobody claims, so a model gives it the metadata of the data it means. Two `Unclaimed` of different COUs may share points. A mapping which is training data in any adopted fit counts as fitted in every fit, so `Unclaimed` and the unclaimed report exclude it.

## The kinds come from the records

The kind of a mapping (`training`, `validation`, `outlier`) is the one the record of an adopted fit stores (`metrics_mappings.tsv`), keyed by experiment and mapping, never by the mapping key alone (two studies can share a key). The assessment compares the recorded kinds with the kinds of the current selection of the fit definitions and fails on any difference, naming the mappings: a split changed after the adoption without a refit would let fitted data count as validation. The remedy is a refit and a new adoption, never an override. The data points are evaluated fresh on the current model and data, not read from the metric tables of the record.

## The verdicts

Per COU, on the union of its validation evidence (each point once; the calibration evidence for C1):

| verdict | meaning |
|---|---|
| `supported` | the acceptance is met on at least the minimum of evaluable points, every claimed value of the domain is covered by validation data (for C1 by calibration data) and every quantity of interest which is no ratio is covered |
| `supported_with_extrapolation` | the acceptance is met on at least the minimum of evaluable points, some claimed value or quantity is not covered, and the COU has an `extrapolation` text which says why the claim holds; a C4 COU, which extrapolates by definition, is never plainly `supported` but at best this |
| `not_supported` | an acceptance rule is not met |
| `insufficient_evidence` | no evidence that can be evaluated, fewer evaluable points than the minimum, or claimed values or quantities not covered and no `extrapolation` text |
| `not_verified` | a check of the verification failed: the decisive predictions change with tighter tolerances or another integration method, or the setup has a problem, so the comparison with the data is not to be trusted |

The verdict is decided in this order, the first step which applies decides:

1. no decisive point: `insufficient_evidence`;
2. a failed check of the verification: `not_verified` (see [Verification](#verification));
3. an acceptance rule not met: `not_supported` (a failure on few points is still a failure);
4. an acceptance rule not evaluable: `insufficient_evidence`;
5. fewer evaluable points than `minimum_points`: `insufficient_evidence`;
6. every claimed value and every quantity which is no ratio covered: `supported` (C4: `supported_with_extrapolation`);
7. an `extrapolation` text: `supported_with_extrapolation`;
8. else `insufficient_evidence`.

Coverage counts only the points the criteria can evaluate: a fold error needs a finite `DV` and `IPRED`, both greater than 0. The numbers are reported as `n_decisive` (the points of the deciding evidence), `n_evaluable` and `minimum_points`; a rule without an evaluable point has the value `NaN` and `met` of `None`, is not evaluable, and gives `insufficient_evidence`, never `supported`. A `not_supported` verdict is a finding to keep and explain, not a reason to loosen the acceptance. A `supported` C1 COU rests on the data the model was fitted to, so the console, the report and the release page show next to every verdict what decided it (decided on: calibration or validation).

Domain coverage is per claimed (field, value): `validated` when validation points with the value exist, `calibrated` when only calibration points exist (this counts as covered for a C1 COU), `uncovered` otherwise. The values of the points outside of the claimed domain are reported as well; evidence that does not fit the claim is a smell of the declaration.

Quantity coverage is per quantity of interest which is no ratio: `covered` when an evaluable decisive point of its observable lies inside the claimed domain, i.e. for every claimed field of the metadata the value of the point is one of the claimed values, `uncovered` otherwise. A quantity validated only outside the claim, e.g. the PT ratio (an observable, no ratio between arms) validated only by a fed study for a fasted claim, is uncovered, although every claimed value may be validated by the points of other quantities. An assessed ratio quantity is covered by an evaluable pair whose numerator has a decisive point inside the claimed domain (its points are these pairs), a ratio without a pairing has no coverage. The coverage is computed for every level, C4 included, and the coverage of combinations of values is not assessed (each claimed value on its own).

## Ratios

A question of a special population is a ratio: how much higher is the exposure in renal impairment than in the healthy subjects of the same study? A fold criterion on the absolute predictions is too lenient for it near a ratio of 1, two predictions within 2-fold each may give a ratio which is 4-fold off. A ratio quantity with a `Ratio` compares the ratio directly (`pkpdsim.credibility.ratios`):

```python
# not executed
from pkpdsim.credibility import QuantityOfInterest, Ratio, Summary
from pkpdsim.fit import metadata_filter
from pkpdsim.metadata import Health

QuantityOfInterest(
    "auc_ratio_renal",
    "the AUC of HCTZ in renal impairment over the healthy one",
    ratio=Ratio(
        observable="Carm_hctz",
        numerator=metadata_filter(health={Health.RENAL_IMPAIRMENT}),
        denominator=metadata_filter(health={Health.HEALTHY}),
        summary=Summary.AUC,
    ),
)
```

- The arms are paired in the decisive points of the COU: per fit and study (the experiment) every mapping of the observable which `numerator` selects with the one mapping which `denominator` selects. A study without its reference arm gives no pair; a study with more than one mapping of the denominator is an error of the selection, so the reference arm is part of the decisive evidence (e.g. a validation evidence of the healthy arms of the studies). The observable of a ratio restricts the evidence of its COU like an absolute quantity.
- The summary of an arm is the same for the data and the predictions: `Summary.AUC` the area under its points by the trapezoid (at least two points with a finite data value and prediction), `Summary.VALUE` the value of its single point (a PK parameter). The observed ratio is the one of the data, the predicted ratio the one of the predictions.
- The limit of Guest et al. (2011) of an observed ratio `R` (its inverse below 1) is `(delta + 2 (R - 1)) / R` with `delta` 1.25: 1.25 at a ratio of 1, 1.625 at 2 or 0.5, towards 2 for a large ratio; a pair is within when the predicted ratio lies between the observed ratio divided by and times the limit.
- The rule `GuestAcceptance(minimum)`, the percentage of the evaluable pairs within the limits, is `ContextOfUse.ratio_acceptance`, by default `DEFAULT_RATIO_ACCEPTANCE` of the risk: 50 % for a low, 67 % for a medium and 80 % for a high risk. It joins the verdict like a rule of the points: not met is `not_supported`, no evaluable pair `insufficient_evidence`. The verification compares it with the alternative predictions as well.
- The factor relevance of the quantities of interest is `assessed` for a COU with an assessed ratio, by its pairs and the coverage of its ratios.

`ratios.tsv` has one row per COU and pair (`context`, `quantity`, `fit`, `experiment`, `numerator`, `denominator`, `observed`, `predicted`, `limit`, `within`), the report and the release page a table of the pairs per COU and the rule of the pairs (`Guest`) in the acceptance, which `contexts.tsv` has as the columns `Guest` and `Guest_limit`. A record of 0.18 or earlier has no pairs.

## Verification

A verdict compares the predictions of the adopted fits (`IPRED`) with the data. Whether these numbers are the solution of the model, and not of the settings of the integrator or of a broken setup, are three factors of V&V 40: the discretization error, the numerical solver error and the use error. `run_credibility` verifies the decisive predictions of every COU with one check per factor (`pkpdsim.credibility.verification`):

| check | factor | what is compared |
|---|---|---|
| `tolerances` | discretization error | every adopted fit evaluated again at its reference parameters with the absolute and the relative tolerance of its settings divided by 100 (`TIGHTENING`; every kind and override of an `AbsoluteTolerance`) |
| `integrator` | numerical solver error | the same at the tightened tolerances with the BDF of CVODE of at most order 2 (`SECOND_METHOD_ORDER`) instead of order 5, a second integration method: another formula with its own local error, which stays stable on a stiff model |
| `setup` | use error | every build of the model is valid SBML with consistent units (the consistency checks of libsbml with the units, every issue of severity warning or above) and every decisive prediction is finite |

The difference of a prediction is relative, `|alt - ref| / max(|ref|, |alt|, floor)`, with the floor `1e-6` (`SCALE_FLOOR`) of the largest reference prediction of its mapping: a difference below a millionth of the peak of the curve is none, and a point at zero compares with the curve. A numerical check of a COU passes when the largest difference of its decisive predictions is at most `0.01` (`LIMIT`, more than an order below the narrowest default criterion, 1.5-fold and a GMFE of 1.3) and no acceptance rule changes its result with the alternative predictions; a prediction which is finite in the reference and not in the alternative is an infinite difference, and an evaluation of a fit which fails is a failed check with its message. A failed check makes the verdict `not_verified`; the three factors are `assessed` with the result as their statement.

The checks run with every assessment and evaluate every adopted fit twice more; `--no-verification` skips them and the battery below for a quick look, the three factors are then `not_assessed` and the verdicts are those of the comparison alone. The non-stiff method of CVODE (Adams-Moulton) is no second method for these models: on the stiff PBPK models it takes too many steps (hctz: `CV_TOO_MUCH_WORK` on its long simulations), while the BDF of order 2 is A-stable. What the verification of the predictions does not do: the experiments on a second simulator (they run on roadrunner through sbmlsim; the battery below compares the builds with AMICI), the explicit Runge-Kutta integrator of roadrunner (sbmlsim sets the settings of CVODE only, and an explicit method is no method for a stiff model) and a check of the doses of the simulations against the interventions of the data (the link of a simulation to an intervention is not declared).

## Battery of the builds

The checks above ask whether the predictions a verdict is decided on are the solution of the model. The battery (`pkpdsim.credibility.battery`, `run_battery`) asks whether the model itself is sound, without any data: every build (`model_paths`) is simulated with roadrunner at `rtol 1e-10` and `atol 1e-16` for two days on a 10 min grid, without a dose and with a dose of 10 mg of every route of every substance (an infusion over 60 min), and the results are checked against what every PBPK model must fulfil:

| check | factor | what is checked |
|---|---|---|
| `mass_balance` | numerical code verification | the amount of the substances in the model, in its sinks and in its dose states (in mmole, a dose state in mg divided by the molecular weight) changes by the infusion and the flows of the reactions which read no dose alone, at most `1e-6` of the largest term (`MASS_LIMIT`): a dose enters with its molecular weight, and nothing appears or disappears outside of the reactions (a rule or an event on a species of a substance is a problem of its own) |
| `non_negativity` | numerical code verification | no species becomes negative beyond `1e-6` of its largest value in any simulation of the battery |
| `zero_dose` | numerical code verification | without a dose the species of a substance which starts at zero (a drug, not an endogenous substance) stay zero |
| `steady_state` | numerical code verification | without a dose every state but the sinks (the states nothing drains) stays within `1e-6` of its start (`STEADY_LIMIT`): the physiology and the endogenous substances start at their steady state |
| `superposition` | numerical code verification | twice a dose of 0.01 mg gives twice the response and two doses at 0 and 12 h the sum of their responses (the response is the difference to the simulation without a dose), within `0.01` of the largest response: the model is linear and time invariant in the limit of a small dose |
| `cross_simulator` | numerical code verification | AMICI, a second simulator with its own compilation of the SBML, gives the results of roadrunner without a dose and with every bolus dose and constant infusion, within `0.01` |
| `plausibility` | use error | the default individual has finite values, positive compartments and species which are not negative, and no value of an adopted fit lies within `1e-3` of the range of a bound (in the scale of the fit, `BOUND_MARGIN`): a value at a bound is decided by the bound, not by the data |
| `rebuild` | software quality assurance | the factory of the model writes its SBML files again, equal but for the dates of the model history; with the current pkpdsim, so the builds are those of the current library |

The battery reads the conventions of the models of pkdb_models from the flat SBML: a substance has its molecular weight `Mr_<s>`, its doses are the states `IVDOSE_<s>`, `PODOSE_<s>` and `SOLDOSE_<s>` (mg) and the infusion rate `Ri_<s>` (mg/min), and a species belongs to the substance its id (the part after the last `__`) ends with, `_<s>`, or else starts with, `<s>_`. A parameter in mg whose rate rule reads a dose (the drug in the stomach of some models) is a state of the dose while a reaction reads it; a sum of the given doses is none. A reaction is an input when it reads a dose of a substance it forms; a secretion which a dose of another substance stimulates is a flow. A check without what it needs is not applicable and assesses nothing: the checks of the simulations without a substance with a molecular weight, the mass balance and the superposition without a dose, the zero dose without a substance which starts at zero, the second simulator without AMICI (the extra `amici`: `pip install pkpdsim[amici]`, which builds AMICI from source and compiles every build with a C++ compiler, cmake and ninja, about a minute and a half for a PBPK model, cached in the temporary directory by the SBML and the version of AMICI; a build which AMICI does not compile within 15 min, `AMICI_TIMEOUT`, e.g. the creatinine model of pkdb_models with its large nested formulas of the starts, which did not compile within 45 min, is not compared and named in the statement of the check).

The mass balance records the fate of every dose (`Fate`): at the end of the two days, relative to the simulation without a dose, the amount in the body, in the urine and the feces, not yet applied, and the amount every reaction which reads no dose moved (negative for an elimination, e.g. a metabolism to products the model does not describe, positive for a production, e.g. of a metabolite by the fraction of a branch).

A failed check of the battery makes every context of use `not_verified`. Its checks are rows of the verification of every context, and the factors numerical code verification and software quality assurance are `assessed` with their platform statement followed by the results of the battery. A model whose build fails a check for a reason it accepts declares it with the reason, `CheckSettings(battery={"steady_state": "..."})`: `pkpdsim.testing.check_battery` (the test of a repository of models over its models) then passes, fails when the known failure no longer applies, and an assessment shows the reason and still fails the check.

## Cross-validation

A fit describes its training data; how well it predicts a study it has not seen is a question of the form of the model. `run_fit <package> --cross-validation <fit id>` refits the problem of an adopted fit once per training study without that study (its training mappings become validation mappings, `held_out`), by least squares from the adopted parameters with the adopted fits of the other problems fixed, and predicts the held-out study with the refitted parameters (`pkpdsim.credibility.cross_validation`, `cross_validate`). `--experiments` limits the held-out studies, `--cores` and `--seed` are those of the refits, which are written below `cross_validation/` of the output directory. The result is stored in the record of the adopted fit, `cross_validation.json` (`write_cross_validation`): per study its points, the GMFE and the share within 2-fold of the held-out predictions and of the predictions of the adopted fit for the same points, and the cost of the refit, and the GMFE of all held-out points together.

`run_credibility` reads the cross-validations of the adopted fits (`read_cross_validations`): the factor model form of a COU is `assessed` when an adopted fit of its decisive evidence is cross-validated, with the statement of the cross-validation; the record has them (`CredibilityAssessment.cross_validations`), `cross_validation.tsv` one row per fit and study, the report and the release page a section Cross-validation. A held-out GMFE close to the fitted one is a model form which carries across studies, a much larger one a model whose parameters are fitted to every study. The cross-validation fits the problem once per training study, so it runs on demand and is stored with the fit; a record of 0.18 or earlier has none.

## Locked splits

A validation counts only if its data was set apart before the predictions were seen. The kinds of the mappings (training, validation) are fixed by the records of the adopted fits, a changed kind is refused (`KindMismatchError`); which validation data answers which question is fixed by the declaration, which can be changed after a verdict was seen. The lock fixes it as well (`pkpdsim.credibility.lock`): `run_credibility <package> --lock` writes `credibility.lock.json` into the folder of the model, the mappings (`fit/experiment/mapping`) of every evidence and every limitation by `context/name` (a limitation of the model `/name`), to be committed and reviewed with the declaration. Afterwards:

- `run_credibility` compares the evidence with the lock before it evaluates anything and refuses a difference (`LockMismatchError`), naming every evidence which is new, gone or selects other mappings with the added (`+`) and removed (`-`) mappings; `check_declaration`, the test of a repository of models, reports the same.
- A changed split is renewed explicitly with `--lock`, so the change is a reviewed diff of the lock in the history of the repository.
- The record says whether the evidence matched the lock (`Provenance.lock`: `matched`, `written` when the lock was written with the assessment, `None` without a lock), the report shows it in the provenance.

A model without a lock is assessed as before.

## Across models: scorecard and qualification

Every model is assessed on its own. The scorecard reads the records of many models (`credibility.json`) into one table, one row per model and COU with the level, the risk, the verdict, the points, whether it is verified, the failed checks and the value of every acceptance rule (`pkpdsim.credibility.scorecard`, `scorecard`). Stored as a TSV it is the baseline of a repository of models: a change of the shared library, of sbmlsim or of the data which makes a COU worse shows across all models at once (`compare`): a COU of the baseline which is missing, a worse verdict (supported, supported with extrapolation, insufficient evidence, not verified, not supported), a COU which is no longer verified and a rule whose value got worse by more than the tolerance (5 points of a percentage, 5 % of a GMFE, an AFE 5 % further from 1; a value which is no longer evaluable). A new model, COU or rule is no regression.

The qualification pools the ratio pairs of all models by a value of the metadata of their numerator arm (`qualification`, by default the field `health`): per value, e.g. `RENAL_IMPAIRMENT` or `CIRRHOSIS`, the evaluable pairs of every model and of all models, how many lie within the limits of Guest et al. and the GMFE of the predicted over the observed ratios. The renal and hepatic scaling of `pkpdsim.pbpk.impairment` is shared by the models; it is qualified for a population when the ratios of the drugs are predicted across drugs, not only for one. A pair counts once even when two COUs of a model share it, and a numerator whose value is not reported (`NR`) counts for no value.

```bash
pkpdsim-verify scorecard results/ --out scorecard.tsv --qualification qualification.tsv
pkpdsim-verify scorecard results/ --baseline scorecard.tsv
```

The command prints the verdicts, the pooled qualification and the regressions against `--baseline`; it exits with 1 for a regression and 2 for a record or baseline it cannot read. `--field` pools by another field of the metadata (`genotype`, `coadministration`), `--tolerance` sets the tolerance.

## Running and reading the report

```bash
run_credibility pkdb_models.models.hctz
run_credibility pkdb_models.models.hctz --no-verification
run_credibility pkdb_models.models.hctz --from-json results/credibility/credibility.json
```

The console shows a table of the COUs, their verdicts, what decided them and the evaluable points as `evaluable/minimum`, and below it how many COUs are verified; see [Commands](commands.md). `index.html` of the report (one page, no network) has:

- the **overview**: the COUs with level, risk, verdict (in words and in a colour), the role of the evidence which decided it and the decisive and evaluable points with the minimum, and the limitations of the model;
- a **section per COU**: question, rationale, the evaluable points with their minimum, quantities with their coverage and points (a ratio without a pairing marked as such), the ratio pairs with their limits, the claimed domain with its coverage and the statement (not assessed), the evidence table (role, points, mappings, studies, GMFE, AFE, the percentages of the criteria, the points removed by limitations and outliers), the acceptance table (value, limit, met), the checks of the verification (points, largest difference and its point, whether the acceptance changed, the result and the problems), the table of the credibility factors and the extrapolation text;
- the **limitations**: their metrics, explicitly outside of any acceptance;
- the **unclaimed validation data**: the validation points per fit which no evidence and no limitation selects (held-out data which answers no declared question);
- the **verification**: the settings of the integrator of every adopted fit and of its alternatives, the issues of the builds, the battery of the builds with the statement of every check and the fate of every dose;
- the **provenance**: the versions of pkpdsim, sbmlsim, sbmlutils and libroadrunner, the commit of the model repository and whether the folder of the model had changes (`git status` of the model folder), the adopted fits and the date. A release passes the commit of its source repository (`--source-commit`), the assessment then records it and leaves the dirty state unknown.

The factors of V&V 40 are listed per COU with a status: `assessed` (output comparison, test samples, test conditions and the relevance of the validation activities, from the metrics, the evidence and the coverage; the discretization error, the numerical solver error and the use error from the verification, the numerical code verification and the software quality assurance from the battery after their platform statement), `platform` (a statement which holds for the whole stack, e.g. the code verification by the SBML Test Suite and the PEtab suites of sbmlsim) or `not_assessed` with the step which will assess it, so a reader sees what is missing. Step 1 checks that every claimed quantity which is no ratio has decisive points in the claimed domain (the quantity coverage, under the relevance of the validation activities); whether the quantities answer the question, e.g. the ratios of a special population, is the relevance of the quantities of interest of step 3.

The record `credibility.json` holds everything, and the tables `contexts.tsv`, `evidence.tsv`, `domain.tsv` and `quantities.tsv` are one row per COU, per COU and evidence, per COU, field and value and per COU and quantity, `ratios.tsv` one row per COU and ratio pair (see [Ratios](#ratios)) and `verification.tsv` one row per COU and check (the applicable checks of the battery included, the same for every COU), `battery.tsv` one row per check of the battery (`check`, `factor`, `status`, `cases`, `largest`, `at`, `statement`, `problems`) and `fates.tsv` one row per build, dose and part of its fate (`build`, `dose`, `substance`, `amount`, `part`, `value` in mmole and its `share` of the dose). `contexts.tsv` has the decisive and the evaluable points, the minimum (`minimum_points`), the uncovered quantities (`uncovered_quantities`, the ids comma separated, empty without) and whether the COU is verified (`verified`, empty without verification); `quantities.tsv` has the columns `context`, `quantity`, `description`, `ratio`, `coverage` (empty for a ratio without a pairing) and `n_points`; `verification.tsv` the columns `context`, `check`, `factor`, `n`, `max_difference`, `at`, `acceptance_changed`, `passed` and `problems` (joined by `; `), no rows without verification. A record written by pkpdsim 0.17 is still read (`--from-json`): it has no minimum (read as 0), no quantity coverage and no verification; a record of 0.18 has no verification and no ratio pairs. A release of a model repository carries the assessment, see [Release of a model repository](release.md).

## What the later steps add

Steps 1 and 2 (the verification) are this page. The next steps plug into the same declarations and the same report: ratio metrics between arms (impaired over healthy), leave-one-study-out cross-validation and the qualification of the platform across models, which a C3 and a C4 COU need; and the uncertainty and sensitivity of the quantities of interest.

## References

- ASME V&V 40-2018, Assessing Credibility of Computational Modeling through Verification and Validation: Application to Medical Devices.
- FDA, Assessing the Credibility of Computational Modeling and Simulation in Medical Device Submissions (2023).
- ICH M15, General Principles for Model-Informed Drug Development (Step 4, 2026).
- EMA, Guideline on the reporting of physiologically based pharmacokinetic (PBPK) modelling and simulation (2018).
- Kuemmel C et al. (2020), Consideration of a credibility assessment framework in model-informed drug development: potential application to physiologically-based pharmacokinetic modeling and simulation, CPT Pharmacometrics Syst Pharmacol 9:21-28.
- Musuamba FT et al. (2021), Scientific and regulatory evaluation of mechanistic in silico drug and disease models in drug development: building model credibility, CPT Pharmacometrics Syst Pharmacol 10:804-825.
- Erdemir A et al. (2020), Credible practice of modeling and simulation in healthcare: ten rules from a multidisciplinary perspective, J Transl Med 18:369.
- Guest EJ, Aarons L, Houston JB, Rostami-Hodjegan A, Galetin A (2011), Critique of the two-fold measure of prediction success for ratios: application for the assessment of drug-drug interactions, Drug Metab Dispos 39:170-173.
