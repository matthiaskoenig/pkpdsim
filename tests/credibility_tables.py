"""Hand-built tables of the data points of the credibility test package.

The tables are assessed by `sbmlsim.credibility.assess_table` with the evidence
extensions of pkpdsim (`assess`), as a `PKPDModel` registers them.
"""

from collections.abc import Mapping
from types import MappingProxyType

import numpy as np
import pandas as pd
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.experiments import CredA, CredB
from sbmlsim.credibility import (
    CHECKS,
    Battery,
    BatteryResult,
    BatteryStatus,
    Check,
    Credibility,
    CredibilityAssessment,
    FitPropagation,
    FitTable,
    PredictiveChecks,
    PropagationSettings,
    Provenance,
    VariabilitySettings,
    Verification,
    VerificationSummary,
    assess_table,
)
from sbmlsim.credibility.uncertainty import FitDraws
from sbmlsim.credibility.verification import (
    INTEGRATOR,
    LIMIT,
    SCALE_FLOOR,
    TIGHTENING,
    TOLERANCES,
)
from sbmlsim.fit import FitMapping
from sbmlsim.fit.helpers import FitMappings
from sbmlsim.fit.metrics import metadata_dict
from sbmlsim.fit.variability import Variability
from sbmlsim.simulation.sampling import LogNormal

from pkpdsim.credibility import EXTENSIONS, PK_BATTERY, Fate

#: the provenance of the assessments of the hand-built tables
PROVENANCE = Provenance(
    model="credibility",
    title="Toy model with every kind of evidence",
    fits=["DECAY_20261001_120000__abcde"],
    versions={"pkpdsim": "0.0.0"},
    commit=None,
    dirty=None,
    date="2026-11-04T12:00:00+01:00",
)

#: the kinds of the adopted fit of the test package, see `test_credibility_package`
KINDS = {
    ("CredA", "fm_train"): "training",
    ("CredA", "fm_dose"): "validation",
    ("CredA", "fm_renal"): "validation",
    ("CredA", "fm_limit"): "validation",
    ("CredA", "fm_subject"): "validation",
    ("CredA", "fm_outlier"): "outlier",
    ("CredB", "fm_dose"): "training",
    ("CredB", "fm_zero"): "validation",
}
#: the number of data points of a mapping of the hand-built table, else 5
POINTS = {("CredA", "fm_limit"): 2, ("CredA", "fm_outlier"): 3}
#: the value of the data of a mapping, else 1: three times the decay and zeros
VALUES = {("CredA", "fm_limit"): 3.0, ("CredB", "fm_zero"): 0.0}


def fit_table(
    factor: float = 1.0,
    factors: Mapping[tuple[str, str], float] | None = None,
    *,
    times: bool = False,
) -> FitTable:
    """Get a table of the data points of the kinds of the test package.

    Every prediction is 1, the decay the model predicts, times a factor.

    Args:
        factor: the factor of the predictions.
        factors: the factor of the predictions of a mapping instead of `factor`.
        times: give the points of a mapping the times `x` 0, 1, 2, ...

    Returns:
        The table with the columns `DV`, `IPRED` and `DV_SD`.
    """
    factors = factors or {}
    mappings = FitMappings(
        experiment_classes=[CredA, CredB],
        base_path=MODEL.path,
        data_path=list(MODEL.data_paths),
    )

    def lookup(experiment: str, mapping: str) -> FitMapping:
        return mappings.experiments[experiment]._fit_mappings[mapping]

    rows = [
        {
            "fit": "DECAY",
            "experiment": experiment,
            "mapping": mapping,
            "kind": kind,
            "observable": "[A]",
            **metadata_dict(lookup(experiment, mapping).metadata, mapping),
            "DV": VALUES.get((experiment, mapping), 1.0),
            "IPRED": factors.get((experiment, mapping), factor),
            "DV_SD": 0.1,
            **({"x": float(k)} if times else {}),
        }
        for (experiment, mapping), kind in KINDS.items()
        for k in range(POINTS.get((experiment, mapping), 5))
    ]
    return FitTable(
        table=pd.DataFrame(rows),
        lookup=lookup,
        recorded_kinds={("DECAY", e, m): k for (e, m), k in KINDS.items()},
    )


def verification_of(
    fits: FitTable,
    factor: float = 1.0,
    *,
    factors: dict[Check, float] | None = None,
    failures: dict[Check, dict[str, str]] | None = None,
    problems: tuple[str, ...] = (),
    limit: float = LIMIT,
    battery: Battery | None = None,
) -> Verification:
    """Get a verification whose alternatives are the reference times a factor."""
    factors = factors or {}
    failures = failures or {}
    ipred = fits.table["IPRED"].astype(float)
    return Verification(
        predictions=MappingProxyType(
            {
                check: ipred * factors.get(check, factor)
                for check in (TOLERANCES, INTEGRATOR)
            }
        ),
        failures=MappingProxyType(
            {
                check: MappingProxyType(failures.get(check, {}))
                for check in (TOLERANCES, INTEGRATOR)
            }
        ),
        summary=VerificationSummary(
            tightening=TIGHTENING,
            limit=limit,
            scale_floor=SCALE_FLOOR,
            fits=[
                {
                    "fit": "DECAY",
                    "reference": "CVODE BDF, rtol 1e-06",
                    "tolerances": "CVODE BDF, rtol 1e-08",
                    "integrator": "CVODE BDF of order at most 2, rtol 1e-08",
                }
            ],
            problems=list(problems),
            battery=battery,
        ),
    )


def battery_of(failed: Check | None = None) -> Battery:
    """Get a battery of the PK hook without a simulation.

    The second simulator is not applicable, the fates are one dose.

    Args:
        failed: a check which fails with the problem `broken`, every other
            check passes.
    """
    results = []
    for check in (CHECKS[c] for c in PK_BATTERY.order):
        if check.id == "cross_simulator":
            status = BatteryStatus.NOT_APPLICABLE
        elif check == failed:
            status = BatteryStatus.FAILED
        else:
            status = BatteryStatus.PASSED
        results.append(
            BatteryResult(
                check=check,
                status=status,
                cases=0 if status is BatteryStatus.NOT_APPLICABLE else 4,
                largest=float("nan")
                if status is BatteryStatus.NOT_APPLICABLE
                else 1e-9,
                at="" if status is BatteryStatus.NOT_APPLICABLE else "male IVDOSE_d",
                problems=("broken",) if status is BatteryStatus.FAILED else (),
                statement=f"The check {check.id}.",
            )
        )
    fate = Fate(
        build="male",
        dose="IVDOSE_d",
        substance="d",
        amount=0.1,
        parts={"body": 0.01, "urine": 0.04, "ELIM": -0.05},
        imbalance=1e-15,
    )
    return Battery(
        results=tuple(results),
        sections={"fates": [fate.to_dict()]},
        settings="roadrunner CVODE BDF.",
    )


def propagation_of(spread: float, failed: int = 0) -> FitPropagation:
    """Get predictions per draw of the table, the reference times a factor of the spread."""
    settings = PropagationSettings(draws=256, seed=1)
    table = fit_table().table
    rng = np.random.default_rng(1)
    factors = np.exp(rng.normal(0.0, spread, settings.draws))
    values = np.outer(table["IPRED"].to_numpy(dtype=float), factors)
    values[:, :failed] = np.nan
    draws = FitDraws(
        fit_id="DECAY_20261001_120000__abcde",
        pids=("k",),
        values=factors[:, None],
        units={"k": "1/min"},
        correlated=False,
        note="nothing to correlate",
    )
    return FitPropagation(
        settings=settings,
        draws={"DECAY": draws},
        predictions=pd.DataFrame(values, index=table.index),
        failed={"DECAY": failed},
    )


def predictive_checks_of(width: float) -> PredictiveChecks:
    """Get bands of the points of the table, `width` around the observations.

    A width below 1 puts every observation outside its band.
    """
    table = fit_table().table
    dv = table["DV"].to_numpy(dtype=float)
    low, high = dv / (1 + width), dv * (1 + width)
    if width < 1:
        low, high = dv * 2.0, dv * 3.0
    points = pd.DataFrame(
        {
            "N": 4,
            "statistic": "mean",
            "PRED": table["IPRED"].to_numpy(dtype=float),
            "low": low,
            "median": dv,
            "high": high,
            "inside": [
                bool(lo <= v <= hi) for v, lo, hi in zip(dv, low, high, strict=True)
            ],
        },
        index=table.index,
    )
    return PredictiveChecks(
        settings=VariabilitySettings(subjects=16, groups=10),
        variability=Variability({"k": LogNormal(cv=0.2)}),
        points=points,
        vectors={"DECAY": 1},
        studies={"DECAY": 1},
        failed={"DECAY": 0},
    )


def assess(
    credibility: Credibility,
    fits: FitTable,
    verification: Verification | None = None,
    *,
    propagation: FitPropagation | None = None,
    variability: PredictiveChecks | None = None,
    prepared: Mapping[str, object] | None = None,
) -> CredibilityAssessment:
    """Assess a declaration on a table with the evidence extensions of pkpdsim.

    Args:
        credibility: the declaration.
        fits: the table.
        verification: the verification, none without.
        propagation: the predictions of the draws of the parameters, none
            without.
        variability: the predictive checks, none without.
        prepared: what an extension prepared by its name, e.g. the chronology
            `{"temporal": ...}`; an extension without assesses nothing but the
            ratio pairs, which need nothing prepared.
    """
    return assess_table(
        credibility,
        fits,
        PROVENANCE,
        verification,
        propagation=propagation,
        variability=variability,
        extensions=EXTENSIONS,
        prepared=prepared,
    )
