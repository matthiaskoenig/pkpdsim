"""The credibility of many models at once: the scorecard and the qualification.

A repository of models assesses every model on its own (`run_credibility`);
the scorecard reads their records (`credibility.json`) into one table, one row
per model and context of use, and compares it with a stored baseline, so a
change of the shared library or of the data which makes a verdict or a metric
worse shows across all models at once (`compare`, a `Regression` each). The
qualification pools the ratio pairs of all models by a value of the metadata of
their numerator arm, e.g. the renal impairment of `pkpdsim.pbpk.impairment`:
the scaling the models share is qualified across drugs only if the ratios of
the drugs are predicted within the limits of Guest et al.

`pkpdsim-verify scorecard` runs both on the records of a directory.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd

from pkpdsim.credibility.assessment import CredibilityAssessment, Verdict
from pkpdsim.credibility.report import RECORD

__all__ = [
    "QUALIFICATION_COLUMNS",
    "SCORECARD_COLUMNS",
    "TOLERANCE",
    "VERDICT_RANK",
    "Regression",
    "compare",
    "load_records",
    "qualification",
    "scorecard",
]

#: the columns of the scorecard before the values of the acceptance rules
SCORECARD_COLUMNS: Final = (
    "model",
    "context",
    "title",
    "level",
    "risk",
    "verdict",
    "n_decisive",
    "n_evaluable",
    "verified",
    "failed_checks",
)
#: the columns of the qualification
QUALIFICATION_COLUMNS: Final = (
    "field",
    "value",
    "model",
    "n_pairs",
    "n_within",
    "within",
    "gmfe",
)
#: the rank of a verdict, a higher rank is worse
VERDICT_RANK: Final[Mapping[Verdict, int]] = {
    Verdict.SUPPORTED: 0,
    Verdict.SUPPORTED_WITH_EXTRAPOLATION: 1,
    Verdict.INSUFFICIENT_EVIDENCE: 2,
    Verdict.NOT_VERIFIED: 3,
    Verdict.NOT_SUPPORTED: 4,
}
#: the change of a metric which is no regression: points of a percentage, the
#: share of a fold error
TOLERANCE: Final = 0.05
#: the rules whose value is a fold error: lower is better, the AFE towards 1
_FOLD_ERRORS: Final = ("GMFE",)
_BIAS: Final = ("AFE",)


def load_records(paths: Iterable[Path]) -> dict[str, CredibilityAssessment]:
    """Read the records of assessments, by the model of their provenance.

    Args:
        paths: records (`credibility.json`) or directories, which are searched
            for records.

    Returns:
        The assessments by model, sorted.

    Raises:
        FileNotFoundError: if a path does not exist or no record is found.
        ValueError: if two records are of one model.
    """
    files: list[Path] = []
    for path in paths:
        if path.is_dir():
            files.extend(sorted(path.rglob(RECORD)))
        elif path.is_file():
            files.append(path)
        else:
            raise FileNotFoundError(f"The record or directory '{path}' does not exist.")
    if not files:
        raise FileNotFoundError(f"No record {RECORD} was found.")
    records: dict[str, CredibilityAssessment] = {}
    sources: dict[str, Path] = {}
    for file in files:
        assessment = CredibilityAssessment.from_json(file)
        model = assessment.provenance.model
        if model in records:
            raise ValueError(
                f"The records '{sources[model]}' and '{file}' are both of the "
                f"model '{model}'."
            )
        records[model] = assessment
        sources[model] = file
    return dict(sorted(records.items()))


def scorecard(records: Mapping[str, CredibilityAssessment]) -> pd.DataFrame:
    """Get one row per model and context of use.

    Args:
        records: the assessments by model.

    Returns:
        The columns of `SCORECARD_COLUMNS` (`verified` empty without
        verification, `failed_checks` the failed checks of the verification
        comma separated) and per acceptance rule its value, the rule of the
        ratio pairs `Guest` included.
    """
    rows: list[dict[str, object]] = []
    for model, assessment in records.items():
        for c in assessment.contexts:
            row: dict[str, object] = {
                "model": model,
                "context": c.id,
                "title": c.title,
                "level": c.level.value,
                "risk": c.risk.value,
                "verdict": c.verdict.value,
                "n_decisive": c.n_decisive,
                "n_evaluable": c.n_evaluable,
                "verified": c.verified,
                "failed_checks": ",".join(
                    r["check"] for r in c.verification if not r["passed"]
                ),
            }
            rules = [
                *c.acceptance,
                *([] if c.ratio_acceptance is None else [c.ratio_acceptance]),
            ]
            for rule in rules:
                row[rule["name"]] = rule["value"]
            rows.append(row)
    columns = list(
        dict.fromkeys([*SCORECARD_COLUMNS, *(k for row in rows for k in row)])
    )
    return pd.DataFrame(rows, columns=columns)


@dataclass(frozen=True, kw_only=True)
class Regression:
    """A context of use of a model which got worse than its baseline.

    Attributes:
        model: the model.
        context: the context of use.
        what: what got worse: `missing`, `verdict`, `verified` or the name of
            an acceptance rule.
        baseline: the value of the baseline, as text.
        current: the current value, as text.
    """

    model: str
    context: str
    what: str
    baseline: str
    current: str

    def __str__(self) -> str:
        """Describe the regression in one line."""
        return (
            f"{self.model}/{self.context}: {self.what} {self.baseline} -> "
            f"{self.current}"
        )


def _worse(name: str, baseline: float, current: float, tolerance: float) -> bool:
    """Check whether the value of a rule got worse beyond the tolerance."""
    if math.isnan(baseline):
        return False
    if math.isnan(current):
        return True
    if name in _FOLD_ERRORS:
        return current > baseline * (1 + tolerance)
    if name in _BIAS:
        if baseline <= 0 or current <= 0:
            return current != baseline
        return abs(math.log(current)) > abs(math.log(baseline)) + math.log1p(tolerance)
    # a percentage of the points or pairs, higher is better
    return current < baseline - 100 * tolerance


def compare(
    baseline: pd.DataFrame, current: pd.DataFrame, tolerance: float = TOLERANCE
) -> list[Regression]:
    """Get the contexts of use which got worse than the baseline.

    A context of the baseline which is missing is a regression, and so is a
    worse verdict (`VERDICT_RANK`), a verified context which is no longer
    verified and the value of a rule of both which got worse by more than the
    tolerance: a percentage by `100 * tolerance` points, a GMFE by the share
    `tolerance`, an AFE by `tolerance` further from 1. A value which is no
    longer evaluable is worse, a new context or rule is none.

    Args:
        baseline: a scorecard, see `scorecard`, e.g. read from its TSV.
        current: the current scorecard.
        tolerance: the change of a value which is no regression.

    Returns:
        The regressions, in the order of the baseline.
    """
    found: list[Regression] = []
    rules = [
        c
        for c in baseline.columns
        if c not in SCORECARD_COLUMNS and c in current.columns
    ]
    now = {(str(r["model"]), str(r["context"])): r for _, r in current.iterrows()}
    for _, before in baseline.iterrows():
        key = (str(before["model"]), str(before["context"]))
        after = now.get(key)
        if after is None:
            found.append(
                Regression(
                    model=key[0],
                    context=key[1],
                    what="missing",
                    baseline="present",
                    current="missing",
                )
            )
            continue
        old, new = Verdict(before["verdict"]), Verdict(after["verdict"])
        if VERDICT_RANK[new] > VERDICT_RANK[old]:
            found.append(
                Regression(
                    model=key[0],
                    context=key[1],
                    what="verdict",
                    baseline=old.value,
                    current=new.value,
                )
            )
        if _true(before["verified"]) and not _true(after["verified"]):
            found.append(
                Regression(
                    model=key[0],
                    context=key[1],
                    what="verified",
                    baseline="True",
                    current=str(after["verified"]),
                )
            )
        for rule in rules:
            a, b = _number(before[rule]), _number(after[rule])
            if _worse(rule, a, b, tolerance):
                found.append(
                    Regression(
                        model=key[0],
                        context=key[1],
                        what=rule,
                        baseline=f"{a:.4g}",
                        current=f"{b:.4g}",
                    )
                )
    return found


def _true(value: object) -> bool:
    """Check whether a value of a table is true, also as the text of a TSV."""
    return value is True or value == "True" or value == np.True_


def _number(value: object) -> float:
    """Get a value of a table as a float, `NaN` for an empty one."""
    if value is None or value == "":
        return math.nan
    return float(value)  # ty: ignore[invalid-argument-type]


def qualification(
    records: Mapping[str, CredibilityAssessment], field: str = "health"
) -> pd.DataFrame:
    """Pool the ratio pairs of the models by a value of the metadata of their numerator.

    Every evaluable pair of every context of use counts once per study and
    arms: a pair which two contexts of a model share counts once; a pair
    whose numerator has the value `NR` (not reported) counts for no value.

    Args:
        records: the assessments by model.
        field: the field of the metadata of the numerator, e.g. `health`.

    Returns:
        The columns of `QUALIFICATION_COLUMNS`: per value of the field one row
        per model and one pooled row with the model `all`: the evaluable
        pairs, the pairs within the limits of Guest et al., their percentage
        and the geometric mean fold error of the predicted over the observed
        ratios; sorted by value with the pooled row first.
    """
    pairs: dict[tuple[str, str, str, str, str], tuple[str, float, float, bool]] = {}
    for model, assessment in records.items():
        for c in assessment.contexts:
            for r in c.ratios:
                value = r["metadata"].get(field)
                if r["within"] is None or value in (None, "NR"):
                    # a value which is not reported qualifies nothing
                    continue
                key = (
                    model,
                    r["fit"],
                    r["experiment"],
                    r["numerator"],
                    r["denominator"],
                )
                pairs[key] = (
                    value,
                    r["observed"],
                    r["predicted"],
                    bool(r["within"]),
                )
    table = pd.DataFrame(
        [
            {"model": key[0], "value": v, "observed": o, "predicted": p, "within": w}
            for key, (v, o, p, w) in pairs.items()
        ],
        columns=["model", "value", "observed", "predicted", "within"],
    )
    rows: list[dict[str, object]] = []
    for value, group in table.groupby("value", sort=True):
        for model, part in [("all", group), *group.groupby("model", sort=True)]:
            fold = np.abs(np.log10(part["predicted"] / part["observed"]))
            rows.append(
                {
                    "field": field,
                    "value": value,
                    "model": model,
                    "n_pairs": len(part),
                    "n_within": int(part["within"].sum()),
                    "within": 100.0 * float(part["within"].mean()),
                    "gmfe": float(10 ** fold.mean()),
                }
            )
    return pd.DataFrame(rows, columns=list(QUALIFICATION_COLUMNS))
