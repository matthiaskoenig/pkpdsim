"""The verification of the predictions a credibility assessment decides on.

A context of use is decided on the predictions of the adopted fits at the data
points of its evidence (`IPRED`). Whether these numbers are the solution of the
model, and not of the settings of the integrator or of a broken setup, are the
factors discretization error, numerical solver error and use error of V&V 40.
One check per factor (`Check`):

- `tolerances`: every adopted fit evaluated again at its reference parameters
  with the absolute and the relative tolerance of its settings divided by
  `TIGHTENING` (`tightened`).
- `integrator`: the same at the tightened tolerances with the BDF of CVODE of
  at most order 2 (`SECOND_METHOD_ORDER`) instead of order 5, a second
  integration method: another formula with its own local error, stable on a
  stiff model (the non-stiff Adams-Moulton of CVODE takes too many steps on
  the PBPK models, e.g. on the long simulations of hctz).
- `setup`: every build of the model is valid SBML with consistent units
  (`build_problems`) and every decisive prediction is finite.

A numerical check compares the alternative predictions with the reference
predictions of the assessment (`relative_differences`) and passes when the
largest relative difference of the decisive points of a context is at most
`LIMIT` and the acceptance of the context does not change, see
`pkpdsim.credibility.assessment`.
"""

from __future__ import annotations

import dataclasses
import logging
import math
import time
from dataclasses import dataclass
from enum import StrEnum, unique
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final, Self, TypedDict

import numpy as np
import pandas as pd
from sbmlsim.model.tolerances import AbsoluteTolerance

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path

    from sbmlsim.fit.options import FitSettings
    from sbmlsim.fit.report import FitReport

    from pkpdsim.credibility.battery import Battery
    from pkpdsim.credibility.records import FitTable
    from pkpdsim.model import PKPDModel

__all__ = [
    "BATTERY_CHECKS",
    "LIMIT",
    "LISTED_ISSUES",
    "NUMERICAL_CHECKS",
    "SCALE_FLOOR",
    "SECOND_METHOD_ORDER",
    "TIGHTENING",
    "Check",
    "CheckRecord",
    "FitSettingsRecord",
    "Verification",
    "VerificationSummary",
    "build_problems",
    "check_record_from_dict",
    "check_record_to_dict",
    "relative_differences",
    "tightened",
    "verify",
]

logger = logging.getLogger(__name__)

#: the factor by which the tolerances of the alternatives are divided
TIGHTENING: Final = 100.0
#: the largest relative difference of a decisive prediction a check accepts,
#: more than an order below the narrowest default criterion (1.5-fold, GMFE 1.3)
LIMIT: Final = 0.01
#: the floor of a relative difference, of the largest absolute reference
#: prediction of the mapping: a difference below it is none
SCALE_FLOOR: Final = 1e-6
#: the largest order of the BDF of the second integration method
SECOND_METHOD_ORDER: Final = 2
#: the issues of a build which are listed, the others are counted
LISTED_ISSUES: Final = 5

_LABELS: Final = {
    "tolerances": "tolerances",
    "integrator": "second integration method",
    "setup": "setup of the simulations",
    "mass_balance": "mass balance",
    "non_negativity": "non-negativity",
    "zero_dose": "zero dose",
    "steady_state": "physiological steady state",
    "superposition": "superposition",
    "cross_simulator": "second simulator",
    "plausibility": "plausibility of the values",
    "rebuild": "rebuild",
}
_FACTORS: Final = {
    "tolerances": "discretization_error",
    "integrator": "numerical_solver_error",
    "setup": "use_error",
    "mass_balance": "numerical_code_verification",
    "non_negativity": "numerical_code_verification",
    "zero_dose": "numerical_code_verification",
    "steady_state": "numerical_code_verification",
    "superposition": "numerical_code_verification",
    "cross_simulator": "numerical_code_verification",
    "plausibility": "use_error",
    "rebuild": "software_quality_assurance",
}


@unique
class Check(StrEnum):
    """A check of the verification and the factor of V&V 40 it assesses.

    The checks of the predictions (`tolerances`, `integrator`, `setup`)
    compare the decisive predictions of a context of use; the checks of the
    battery (`pkpdsim.credibility.battery`) check the builds of the model,
    for every context of use alike.
    """

    TOLERANCES = "tolerances"
    INTEGRATOR = "integrator"
    SETUP = "setup"
    MASS_BALANCE = "mass_balance"
    NON_NEGATIVITY = "non_negativity"
    ZERO_DOSE = "zero_dose"
    STEADY_STATE = "steady_state"
    SUPERPOSITION = "superposition"
    CROSS_SIMULATOR = "cross_simulator"
    PLAUSIBILITY = "plausibility"
    REBUILD = "rebuild"

    @property
    def label(self) -> str:
        """Get the check in words, e.g. `second integration method`."""
        return _LABELS[self.value]

    @property
    def factor(self) -> str:
        """Get the id of the factor of V&V 40 the check assesses."""
        return _FACTORS[self.value]


#: the checks which compare alternative predictions with the reference
NUMERICAL_CHECKS: Final = (Check.TOLERANCES, Check.INTEGRATOR)
#: the checks of the battery of the builds, in the order of the report, see
#: `pkpdsim.credibility.battery`
BATTERY_CHECKS: Final = (
    Check.MASS_BALANCE,
    Check.NON_NEGATIVITY,
    Check.ZERO_DOSE,
    Check.STEADY_STATE,
    Check.SUPERPOSITION,
    Check.CROSS_SIMULATOR,
    Check.PLAUSIBILITY,
    Check.REBUILD,
)


class CheckRecord(TypedDict):
    """The result of a check for a context of use.

    A check of the battery (`pkpdsim.credibility.battery`) is a check of the
    builds of the model, the same for every context.

    Attributes:
        check: the check, see `Check`.
        factor: the id of the factor of V&V 40 the check assesses.
        n_points: the decisive points compared, the cases of a check of the
            battery.
        max_difference: the largest relative difference of a decisive
            prediction, `NaN` for the setup and without a point; the largest
            deviation of a check of the battery.
        at: the point of the largest difference, `experiment.mapping at x`,
            empty without one; the case of the largest deviation of a check of
            the battery.
        acceptance_changed: whether a rule of the acceptance is met with the
            alternative predictions and not with the reference or the other
            way, `None` for the setup and the battery.
        passed: whether the check passed.
        problems: what failed beyond the numbers, e.g. a failed evaluation or
            an issue of a build.
    """

    check: str
    factor: str
    n_points: int
    max_difference: float
    at: str
    acceptance_changed: bool | None
    passed: bool
    problems: list[str]


class FitSettingsRecord(TypedDict):
    """The settings of the integrator of an adopted fit and its alternatives.

    Attributes:
        fit: the subset of the fit.
        reference: the settings of the reference predictions, the fit's own.
        tolerances: the settings of the check `tolerances`.
        integrator: the settings of the check `integrator`.
    """

    fit: str
    reference: str
    tolerances: str
    integrator: str


def check_record_from_dict(data: Mapping[str, Any]) -> CheckRecord:
    """Get a check from plain data, `None` as `NaN`.

    Args:
        data: the check as `CheckRecord` with `None` for `NaN`.

    Returns:
        The check.
    """
    value = data["max_difference"]
    return CheckRecord(
        check=data["check"],
        factor=data["factor"],
        n_points=data["n_points"],
        max_difference=math.nan if value is None else float(value),
        at=data["at"],
        acceptance_changed=data["acceptance_changed"],
        passed=data["passed"],
        problems=list(data["problems"]),
    )


def check_record_to_dict(record: CheckRecord) -> dict[str, object]:
    """Get a check as plain data, `NaN` as `None`."""
    value = record["max_difference"]
    return {**record, "max_difference": None if math.isnan(value) else value}


@dataclass(frozen=True, kw_only=True)
class VerificationSummary:
    """What the verification of an assessment compared, as plain data.

    Attributes:
        tightening: the factor the tolerances were divided by.
        limit: the largest relative difference a check accepts.
        scale_floor: the floor of a relative difference, of the largest
            reference prediction of the mapping.
        fits: the settings of every adopted fit and of its alternatives.
        problems: the issues of the builds of the model, see
            `build_problems`.
        battery: the checks of the builds of the model, see
            `pkpdsim.credibility.battery`, `None` if it did not run.
    """

    tightening: float
    limit: float
    scale_floor: float
    fits: list[FitSettingsRecord]
    problems: list[str]
    battery: Battery | None = None

    def to_dict(self) -> dict[str, object]:
        """Get the summary as plain data."""
        return {
            "tightening": self.tightening,
            "limit": self.limit,
            "scale_floor": self.scale_floor,
            "fits": [dict(f) for f in self.fits],
            "problems": list(self.problems),
            "battery": None if self.battery is None else self.battery.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Get the summary from plain data."""
        from pkpdsim.credibility.battery import Battery

        battery = data.get("battery")
        return cls(
            tightening=float(data["tightening"]),
            limit=float(data["limit"]),
            scale_floor=float(data["scale_floor"]),
            fits=[
                FitSettingsRecord(
                    fit=f["fit"],
                    reference=f["reference"],
                    tolerances=f["tolerances"],
                    integrator=f["integrator"],
                )
                for f in data["fits"]
            ],
            problems=list(data["problems"]),
            battery=None if battery is None else Battery.from_dict(battery),
        )


@dataclass(frozen=True, kw_only=True, eq=False)
class Verification:
    """The alternative predictions and the problems an assessment is verified with.

    Attributes:
        predictions: the predictions of every numerical check (`tolerances`,
            `integrator`) per row of the table of the fits, `NaN` where the
            evaluation of a fit failed.
        failures: the failed evaluations of every numerical check, the message
            by the subset of the fit.
        summary: what was compared, see `VerificationSummary`.
    """

    predictions: Mapping[Check, pd.Series]
    failures: Mapping[Check, Mapping[str, str]]
    summary: VerificationSummary

    @property
    def problems(self) -> list[str]:
        """Get the issues of the builds of the model."""
        return list(self.summary.problems)


def tightened(
    tolerance: float | AbsoluteTolerance, factor: float = TIGHTENING
) -> AbsoluteTolerance:
    """Divide an absolute tolerance by a factor, every kind and override.

    Args:
        tolerance: the tolerance, a float is the same for every kind.
        factor: the factor, greater than 1 tightens.

    Returns:
        The tightened tolerance.
    """
    tolerance = AbsoluteTolerance.of(tolerance)
    return AbsoluteTolerance(
        amount=tolerance.amount / factor,
        concentration=tolerance.concentration / factor,
        other=tolerance.other / factor,
        ids={sid: value / factor for sid, value in tolerance.overrides.items()},
    )


def _tightened_settings(settings: FitSettings, factor: float) -> FitSettings:
    """Get the settings of a fit with the tolerances divided by a factor."""
    return dataclasses.replace(
        settings,
        relative_tolerance=settings.relative_tolerance / factor,
        absolute_tolerance=tightened(settings.absolute_tolerance, factor),
    )


def _describe(settings: FitSettings, method: str = "BDF") -> str:
    """Describe the integration of settings, e.g. `CVODE BDF, rtol 1e-06, atol ...`.

    The method of the reference is the BDF of CVODE of variable order up to 5.
    """
    tolerance = AbsoluteTolerance.of(settings.absolute_tolerance)
    return f"CVODE {method}, rtol {settings.relative_tolerance:g}, atol {tolerance}"


def _second_method(report: FitReport) -> None:
    """Integrate the problem of a report with the second integration method.

    The BDF of CVODE of at most `SECOND_METHOD_ORDER`; the setting is given to
    the simulator of the problem and to its loaded models, as
    `OptimizationProblem.initialize` does with the settings of a fit.
    """
    from sbmlsim.simulator import Simulator

    problem = report.problem
    settings = report.settings
    problem.set_simulator(
        Simulator(
            n_workers=1,
            absolute_tolerance=settings.absolute_tolerance,
            relative_tolerance=settings.relative_tolerance,
            variable_step_size=settings.variable_step_size,
            initial_time_step=settings.initial_time_step,
            maximum_bdf_order=SECOND_METHOD_ORDER,
        )
    )
    for model in {id(m): m for m in problem.models}.values():
        model.set_integrator_settings(maximum_bdf_order=SECOND_METHOD_ORDER)


def _alternative_points(
    report: FitReport, check: Check, tightening: float
) -> pd.DataFrame:
    """Evaluate the reference parameters of a report with the settings of a check."""
    from sbmlsim.fit.report import FitReport

    alternative = FitReport(
        problem=report.problem,
        settings=_tightened_settings(report.settings, tightening),
        parameter_sets=report.reference_set,
        mapping_figures=False,
    )
    if check is Check.INTEGRATOR:
        _second_method(alternative)
    return alternative.metrics(alternative.reference_set).datapoints_df()


def _predictions(
    fits: FitTable, check: Check, tightening: float
) -> tuple[pd.Series, dict[str, str]]:
    """Get the predictions of a numerical check per row of the table of the fits.

    Returns:
        The predictions, `NaN` for the rows of a fit whose evaluation failed,
        and the message of every failed evaluation by the subset of its fit.

    Raises:
        ValueError: if the points of an alternative are not the points of the
            reference, in their order.
    """
    table = fits.table
    predictions = pd.Series(np.nan, index=table.index, dtype=float)
    failures: dict[str, str] = {}
    for subset, report in fits.reports.items():
        rows = table.index[table["fit"] == subset]
        start = time.perf_counter()
        try:
            points = _alternative_points(report, check, tightening)
        except (RuntimeError, ValueError) as err:
            # sbmlsim raises a `ValueError` for the data of a failed simulation
            logger.warning(
                "The %s of the fit '%s' failed: %s", check.label, subset, err
            )
            failures[subset] = f"the evaluation of the fit '{subset}' failed: {err}"
            continue
        logger.info(
            "Evaluate the fit '%s' with the %s in %.1f s",
            subset,
            check.label,
            time.perf_counter() - start,
        )
        reference = table.loc[rows]
        if (
            len(points) != len(rows)
            or not (
                (points["experiment"].to_numpy() == reference["experiment"].to_numpy())
                & (points["mapping"].to_numpy() == reference["mapping"].to_numpy())
                & (
                    points["x"].to_numpy(dtype=float)
                    == reference["x"].to_numpy(dtype=float)
                )
            ).all()
        ):
            raise ValueError(
                f"The points of the {check.label} of the fit '{subset}' are not "
                f"the points of its reference."
            )
        predictions.loc[rows] = points["IPRED"].to_numpy(dtype=float)
    return predictions, failures


def relative_differences(
    table: pd.DataFrame, predictions: pd.Series, scale_floor: float = SCALE_FLOOR
) -> pd.Series:
    """Get the relative differences of alternative predictions to the reference.

    The difference of a point is `|alt - ref| / max(|ref|, |alt|, floor)`; the
    floor is `scale_floor` of the largest absolute finite reference prediction
    of the mapping of the point (fit, experiment and key), so a difference below
    it is none and a point at zero compares with the curve. A point without a
    finite reference prediction has no difference (`NaN`, the setup reports
    it); a finite reference with an alternative which is not finite is an
    infinite difference.

    Args:
        table: the data points with `fit`, `experiment`, `mapping` and the
            reference predictions `IPRED`.
        predictions: the alternative predictions per row of the table.
        scale_floor: the floor relative to the largest reference prediction of
            the mapping.

    Returns:
        The relative difference per row of the table.
    """
    reference = table["IPRED"].to_numpy(dtype=float)
    alternative = predictions.reindex(table.index).to_numpy(dtype=float)
    finite = np.isfinite(reference)
    scale = (
        pd.Series(np.where(finite, np.abs(reference), np.nan), index=table.index)
        .groupby([table["fit"], table["experiment"], table["mapping"]], sort=False)
        .transform("max")
        .fillna(0.0)
        .to_numpy(dtype=float)
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        denominator = np.maximum.reduce(
            [np.abs(reference), np.abs(alternative), scale_floor * scale]
        )
        difference = np.abs(alternative - reference) / denominator
    # equal values, also two zeros of a mapping at zero, differ by nothing
    difference = np.where(alternative == reference, 0.0, difference)
    difference = np.where(finite & ~np.isfinite(alternative), np.inf, difference)
    difference = np.where(finite, difference, np.nan)
    return pd.Series(difference, index=table.index, dtype=float)


def build_problems(paths: Mapping[str, Path]) -> list[str]:
    """Get the issues of the SBML of the builds of a model.

    Every build is checked by libsbml with the consistency checks and the
    checks of the units; an issue of severity warning or above is a problem.

    Args:
        paths: the SBML files by the name of their build, e.g. `male`.

    Returns:
        Per build with issues its file, the number of issues and the first
        `LISTED_ISSUES` with their id and message; empty for valid builds.
    """
    import libsbml

    problems: list[str] = []
    for name, path in paths.items():
        doc = libsbml.readSBMLFromFile(str(path))
        doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, True)
        doc.checkConsistency()
        issues = [
            doc.getError(k)
            for k in range(doc.getNumErrors())
            if doc.getError(k).getSeverity() >= libsbml.LIBSBML_SEV_WARNING
        ]
        if not issues:
            continue
        listed = "; ".join(
            f"[{issue.getErrorId()}] {' '.join(issue.getMessage().split())}"
            for issue in issues[:LISTED_ISSUES]
        )
        more = len(issues) - LISTED_ISSUES
        problems.append(
            f"the {name} build '{path.name}' has {len(issues)} "
            f"{'issue' if len(issues) == 1 else 'issues'} of the SBML: {listed}"
            + (f" and {more} more" if more > 0 else "")
        )
    return problems


def _battery(model: PKPDModel) -> Battery:
    """Run the battery of the checks of the builds of a model.

    A failed check which the model knows (`CheckSettings.battery`) gets its
    reason as a further problem and stays failed.
    """
    from pkpdsim.credibility.battery import run_battery

    logger.info("Run the battery of the checks of the builds of '%s'", model.name)
    battery = run_battery(model)
    known = model.checks.battery
    return dataclasses.replace(
        battery,
        results=tuple(
            dataclasses.replace(
                result,
                problems=(*result.problems, f"known failure: {known[result.check]}"),
            )
            if result.passed is False and result.check.value in known
            else result
            for result in battery.results
        ),
    )


def verify(
    model: PKPDModel,
    fits: FitTable,
    tightening: float = TIGHTENING,
    *,
    battery: bool = True,
) -> Verification:
    """Evaluate the adopted fits of a model with the alternatives of the checks.

    Args:
        model: the model, whose builds are checked.
        fits: the data points of the adopted fits with their reports, see
            `pkpdsim.credibility.records.evaluate_fits`.
        tightening: the factor the tolerances are divided by.
        battery: run the battery of the checks of the builds, see
            `pkpdsim.credibility.battery.run_battery`.

    Returns:
        The alternative predictions, the issues of the builds and the battery.

    Raises:
        ValueError: if the table of the fits has no reports.
    """
    if not fits.reports:
        raise ValueError(
            "The table of the fits has no reports to verify, evaluate the fits "
            "with `evaluate_fits`."
        )
    predictions: dict[Check, pd.Series] = {}
    failures: dict[Check, Mapping[str, str]] = {}
    for check in NUMERICAL_CHECKS:
        logger.info("Verify the predictions of '%s': %s", model.name, check.label)
        predictions[check], failed = _predictions(fits, check, tightening)
        failures[check] = MappingProxyType(failed)
    summary = VerificationSummary(
        tightening=tightening,
        limit=LIMIT,
        scale_floor=SCALE_FLOOR,
        fits=[
            FitSettingsRecord(
                fit=subset,
                reference=_describe(report.settings),
                tolerances=_describe(_tightened_settings(report.settings, tightening)),
                integrator=_describe(
                    _tightened_settings(report.settings, tightening),
                    f"BDF of order at most {SECOND_METHOD_ORDER}",
                ),
            )
            for subset, report in fits.reports.items()
        ],
        problems=build_problems(
            {sex.value: path for sex, path in model.model_paths.items()}
        ),
        battery=_battery(model) if battery else None,
    )
    return Verification(
        predictions=MappingProxyType(predictions),
        failures=MappingProxyType(failures),
        summary=summary,
    )
