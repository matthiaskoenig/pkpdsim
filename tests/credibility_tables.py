"""Hand-built tables of the data points of the credibility test package."""

from collections.abc import Mapping

import pandas as pd
from pkpdsim_credibility import MODEL
from pkpdsim_credibility.experiments import CredA, CredB
from sbmlsim.fit import FitMapping
from sbmlsim.fit.helpers import FitMappings
from sbmlsim.fit.metrics import metadata_dict

from pkpdsim.credibility import FitTable

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
