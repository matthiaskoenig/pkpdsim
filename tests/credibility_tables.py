"""Hand-built tables of the data points of the credibility test package."""

from collections.abc import Mapping
from types import MappingProxyType

import pandas as pd
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.experiments import CredA, CredB
from sbmlsim.fit import FitMapping
from sbmlsim.fit.helpers import FitMappings
from sbmlsim.fit.metrics import metadata_dict

from pkpdsim.credibility import (
    BATTERY_CHECKS,
    Battery,
    BatteryResult,
    BatteryStatus,
    Check,
    Fate,
    FitTable,
    Verification,
    VerificationSummary,
)
from pkpdsim.credibility.verification import LIMIT, SCALE_FLOOR, TIGHTENING

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
    factor: float = 1.0, factors: Mapping[tuple[str, str], float] | None = None
) -> FitTable:
    """Get a table of the data points of the kinds of the test package.

    Every prediction is 1, the decay the model predicts, times a factor.

    Args:
        factor: the factor of the predictions.
        factors: the factor of the predictions of a mapping instead of `factor`.

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
        }
        for (experiment, mapping), kind in KINDS.items()
        for _ in range(POINTS.get((experiment, mapping), 5))
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
                for check in (Check.TOLERANCES, Check.INTEGRATOR)
            }
        ),
        failures=MappingProxyType(
            {
                check: MappingProxyType(failures.get(check, {}))
                for check in (Check.TOLERANCES, Check.INTEGRATOR)
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
    """Get a battery without a simulation: the second simulator is not applicable.

    Args:
        failed: a check which fails with the problem `broken`, every other
            check passes.
    """
    results = []
    for check in BATTERY_CHECKS:
        if check is Check.CROSS_SIMULATOR:
            status = BatteryStatus.NOT_APPLICABLE
        elif check is failed:
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
                statement=f"The check {check.value}.",
            )
        )
    return Battery(
        results=tuple(results),
        fates=(
            Fate(
                build="male",
                dose="IVDOSE_d",
                substance="d",
                amount=0.1,
                parts={"body": 0.01, "urine": 0.04, "ELIM": -0.05},
                imbalance=1e-15,
            ),
        ),
        settings="roadrunner CVODE BDF.",
    )
