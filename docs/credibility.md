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
- `QuantityOfInterest(id, description)`: an observable of the fit mappings (`Cve_hctz`) or the observable of a value per simulation (`pk.cmax`). The observable column of the data holds the selection, `[Cve_hctz]` for a concentration, so a quantity is matched after the brackets are removed (`[A]` and `A` are the same quantity); a quantity which matches no observable of the fits is an error. A quantity restricts the evidence of its COU to its observables. `ratio=True` marks a ratio between arms (impaired over healthy), which is assessed in step 3: it is matched against no observable, restricts no evidence and is shown as "ratio, assessed in step 3". Step 1 assesses the absolute quantities, so a COU has at least one quantity which is no ratio.
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

The verdict is decided in this order, the first step which applies decides:

1. no decisive point: `insufficient_evidence`;
2. an acceptance rule not met: `not_supported` (a failure on few points is still a failure);
3. an acceptance rule not evaluable: `insufficient_evidence`;
4. fewer evaluable points than `minimum_points`: `insufficient_evidence`;
5. every claimed value and every quantity which is no ratio covered: `supported` (C4: `supported_with_extrapolation`);
6. an `extrapolation` text: `supported_with_extrapolation`;
7. else `insufficient_evidence`.

Coverage counts only the points the criteria can evaluate: a fold error needs a finite `DV` and `IPRED`, both greater than 0. The numbers are reported as `n_decisive` (the points of the deciding evidence), `n_evaluable` and `minimum_points`; a rule without an evaluable point has the value `NaN` and `met` of `None`, is not evaluable, and gives `insufficient_evidence`, never `supported`. A `not_supported` verdict is a finding to keep and explain, not a reason to loosen the acceptance. A `supported` C1 COU rests on the data the model was fitted to, so the console, the report and the release page show next to every verdict what decided it (decided on: calibration or validation).

Domain coverage is per claimed (field, value): `validated` when validation points with the value exist, `calibrated` when only calibration points exist (this counts as covered for a C1 COU), `uncovered` otherwise. The values of the points outside of the claimed domain are reported as well; evidence that does not fit the claim is a smell of the declaration.

Quantity coverage is per quantity of interest which is no ratio: `covered` when an evaluable decisive point of its observable lies inside the claimed domain, i.e. for every claimed field of the metadata the value of the point is one of the claimed values, `uncovered` otherwise. A quantity validated only outside the claim, e.g. a coagulation ratio validated by a fed study for a fasted claim, is uncovered, although every claimed value may be validated by the points of other quantities. A ratio quantity has no coverage, it is assessed in step 3. The coverage is computed for every level, C4 included, and the coverage of combinations of values is not assessed (each claimed value on its own).

## Running and reading the report

```bash
run_credibility pkdb_models.models.hctz
run_credibility pkdb_models.models.hctz --from-json results/credibility/credibility.json
```

The console shows a table of the COUs, their verdicts, what decided them and the evaluable points as `evaluable/minimum`; see [Commands](commands.md). `index.html` of the report (one page, no network) has:

- the **overview**: the COUs with level, risk, verdict (in words and in a colour), the role of the evidence which decided it and the decisive and evaluable points with the minimum, and the limitations of the model;
- a **section per COU**: question, rationale, the evaluable points with their minimum, quantities with their coverage and points (a ratio marked as assessed in step 3), the claimed domain with its coverage and the statement (not assessed), the evidence table (role, points, mappings, studies, GMFE, AFE, the percentages of the criteria, the points removed by limitations and outliers), the acceptance table (value, limit, met), the table of the credibility factors and the extrapolation text;
- the **limitations**: their metrics, explicitly outside of any acceptance;
- the **unclaimed validation data**: the validation points per fit which no evidence and no limitation selects (held-out data which answers no declared question);
- the **provenance**: the versions of pkpdsim, sbmlsim, sbmlutils and libroadrunner, the commit of the model repository and whether the folder of the model had changes (`git status` of the model folder), the adopted fits and the date. A release passes the commit of its source repository (`--source-commit`), the assessment then records it and leaves the dirty state unknown.

The factors of V&V 40 are listed per COU with a status: `assessed` (output comparison, test samples, test conditions and the relevance of the validation activities, from the metrics, the evidence and the coverage), `platform` (a statement which holds for the whole stack, e.g. the code verification by the SBML Test Suite and the PEtab suites of sbmlsim) or `not_assessed` with the step which will assess it, so a reader sees what is missing.

The record `credibility.json` holds everything, and the tables `contexts.tsv`, `evidence.tsv`, `domain.tsv` and `quantities.tsv` are one row per COU, per COU and evidence, per COU, field and value and per COU and quantity. `contexts.tsv` has the decisive and the evaluable points, the minimum (`minimum_points`) and the uncovered quantities (`uncovered_quantities`, the ids comma separated, empty without); `quantities.tsv` has the columns `context`, `quantity`, `description`, `ratio`, `coverage` (empty for a ratio) and `n_points`. A record written by pkpdsim 0.17 is still read (`--from-json`): it has no minimum (read as 0) and no quantity coverage. A release of a model repository carries the assessment, see [Release of a model repository](release.md).

## What the later steps add

Step 1 is this page. The next steps plug into the same declarations and the same report: verification checks (mass balance, convergence of the tolerances, cross-simulator comparison) for the factors shown as not assessed; ratio metrics between arms (impaired over healthy), leave-one-study-out cross-validation and the qualification of the platform across models, which a C3 and a C4 COU need; and the uncertainty and sensitivity of the quantities of interest.

## References

- ASME V&V 40-2018, Assessing Credibility of Computational Modeling through Verification and Validation: Application to Medical Devices.
- FDA, Assessing the Credibility of Computational Modeling and Simulation in Medical Device Submissions (2023).
- ICH M15, General Principles for Model-Informed Drug Development (Step 4, 2026).
- EMA, Guideline on the reporting of physiologically based pharmacokinetic (PBPK) modelling and simulation (2018).
- Kuemmel C et al. (2020), Consideration of a credibility assessment framework in model-informed drug development: potential application to physiologically-based pharmacokinetic modeling and simulation, CPT Pharmacometrics Syst Pharmacol 9:21-28.
- Musuamba FT et al. (2021), Scientific and regulatory evaluation of mechanistic in silico drug and disease models in drug development: building model credibility, CPT Pharmacometrics Syst Pharmacol 10:804-825.
- Erdemir A et al. (2020), Credible practice of modeling and simulation in healthcare: ten rules from a multidisciplinary perspective, J Transl Med 18:369.
- Guest EJ, Aarons L, Houston JB, Rostami-Hodjegan A, Galetin A (2011), Critique of the two-fold measure of prediction success for ratios: application for the assessment of drug-drug interactions, Drug Metab Dispos 39:170-173.
