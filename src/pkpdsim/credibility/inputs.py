"""The equivalency of the inputs: the simulated subjects and doses against the reported ones.

The decisive evidence of a context of use supports it only for the inputs the
simulations of the evidence had. An experiment simulates the subjects of its
study (the body weight and the height its simulations set, else those of the
model) and its doses; the study reports its subjects (`TabGroups`,
`TabIndividuals`) and its interventions (`TabInterventions`), the tables of the
collection (`pkpdsim.data`). `input_checks` compares, for every mapping of the
adopted fits, the simulated inputs with the reported ones of the groups, the
individuals and the interventions of the rows of its data, without a
simulation:

- a covariate (`Credibility.covariates`, the target of the model and the
  measurement type of PK-DB, `BW` the `weight` and `HEIGHT` the `height` by
  default) is `equivalent` within `COVARIATE_TOLERANCE` of the reported mean,
  median or geometric mean (the value of an individual), `in_range` within the
  reported range (the minimum and the maximum, else the mean with `SPREAD`
  standard deviations), else it `differs`; a group without the measurement
  reads the one of its parent, an individual without it the one of its group;
- a dose (a target of `dose_prefixes`, in the pre-initialization changes or in
  a change) is `equivalent` within `DOSE_TOLERANCE` of a dose of the reported
  interventions of the data, a dose per body weight times the simulated body
  weight, else it `differs`;
- an input without a reported value is `not_reported`.

The values are compared in the unit of the target in the model; a target
which a fit fits is no input and not checked. A check
decides nothing about a verdict: `inputs_record` summarizes the checks of the
decisive mappings of a context of use, with the range of every input the
evidence covers, and the factor equivalency of inputs states it. The
simulations of a scan are checked with the inputs of their base simulation.

The checks are the evidence extension `INPUTS` (`InputsExtension`, name
`inputs`) of the credibility of sbmlsim: the checks of every mapping are the
section `inputs` of the record of an assessment (`inputs.tsv`), the summary of
the decisive mappings of a context is the section `inputs` of its record
(`input_ranges.tsv`, the block Inputs of the report after the parameter
uncertainty), its statement comes before the one of the variability in the
factor equivalency of inputs, and the scorecard has the inputs of the decisive
mappings which differ (`inputs_differ`). The covariates are the ones of the
declaration (`Credibility.covariates`), `DEFAULT_COVARIATES` when it declares
none.
"""

from __future__ import annotations

import json
import logging
import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from functools import cache
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final, TypedDict, override

import pandas as pd
import pint
from pint.errors import DimensionalityError, UndefinedUnitError
from sbmlsim.credibility.extensions import (
    ContextBlock,
    ContextContribution,
    EvidenceExtension,
    FactorStatement,
)

from pkpdsim.credibility.blocks import render_macro
from pkpdsim.credibility.studies import study_folders

if TYPE_CHECKING:
    from pathlib import Path

    from pint import UnitRegistry
    from sbmlsim.credibility.declaration import ContextOfUse, Credibility
    from sbmlsim.credibility.records import FitTable
    from sbmlsim.fit import FitMapping
    from sbmlsim.fit.optimization import OptimizationProblem
    from sbmlsim.project import Project
    from sbmlsim.units import Quantity

__all__ = [
    "COVARIATE_TOLERANCE",
    "DEFAULT_COVARIATES",
    "DOSE_TOLERANCE",
    "INPUTS",
    "INPUT_COLUMNS",
    "INPUT_RANGE_COLUMNS",
    "NAMED_DIFFERENCES",
    "SPREAD",
    "InputCheck",
    "InputStatus",
    "InputSummary",
    "InputsExtension",
    "InputsRecord",
    "Reported",
    "StudyTables",
    "check_covariate",
    "check_dose",
    "dose_prefixes",
    "input_checks",
    "inputs_record",
    "inputs_statement",
    "read_study_tables",
    "reported_doses",
    "reported_subjects",
]

logger = logging.getLogger(__name__)

#: the covariates of the subjects whose simulated values are checked against
#: the reported ones when the declaration declares none
#: (`Credibility.covariates`): the target of the model and the measurement
#: type of PK-DB
DEFAULT_COVARIATES: Final[Mapping[str, str]] = MappingProxyType(
    {"BW": "weight", "HEIGHT": "height"}
)

#: the relative difference of a covariate from the reported central value
#: within which it is equivalent
COVARIATE_TOLERANCE = 0.1

#: the relative difference of a dose from a reported dose within which it is
#: equivalent
DOSE_TOLERANCE = 0.02

#: the standard deviations around the mean of the reported range of a group
#: without its minimum and maximum
SPREAD = 2.0

#: the columns of the table of the checks, see `InputCheck`
INPUT_COLUMNS = (
    "fit",
    "experiment",
    "mapping",
    "target",
    "kind",
    "measurement",
    "unit",
    "simulated",
    "source",
    "reported",
    "low",
    "high",
    "reported_by",
    "status",
    "ratio",
)

#: the value of a change which adds a dose to a dose state, see
#: `pkpdsim.experiment.PKPDSimulationExperiment.add_doses`
_ADDED_DOSE = re.compile(r"^\s*(?P<target>\w+)\s*\+\s*(?P<value>[-+0-9.eE]+)\s*$")


class InputStatus(StrEnum):
    """How a simulated input compares with the reported one."""

    EQUIVALENT = "equivalent"
    IN_RANGE = "in_range"
    DIFFERS = "differs"
    NOT_REPORTED = "not_reported"


class InputCheck(TypedDict):
    """The check of a simulated input of a mapping against the reported one.

    Attributes:
        fit: the subset of the fit.
        experiment: the experiment, the study.
        mapping: the key of the mapping.
        target: the target of the model.
        kind: `covariate` or `dose`.
        measurement: the measurement type of the reported value (`weight`,
            `height`, ..., `dosing`).
        unit: the unit of the target in the model, of every value.
        simulated: the simulated value.
        source: `simulation` (set by the simulation) or `model` (the value of
            the model).
        reported: the reported central value (mean, median, geometric mean,
            the value of an individual, the dose) the simulated one is
            compared with, `None` without.
        low: the low end of the reported range, `None` without.
        high: the high end of the reported range, `None` without.
        reported_by: the groups, individuals or interventions of the reported
            value.
        status: the result, see `InputStatus`.
        ratio: the simulated over the reported central value, `None` without.
    """

    fit: str
    experiment: str
    mapping: str
    target: str
    kind: str
    measurement: str
    unit: str
    simulated: float
    source: str
    reported: float | None
    low: float | None
    high: float | None
    reported_by: list[str]
    status: str
    ratio: float | None


class InputSummary(TypedDict):
    """The checks of one input over the decisive mappings of a context.

    Attributes:
        target: the target of the model.
        kind: `covariate` or `dose`.
        unit: the unit of the target in the model.
        n_mappings: the decisive mappings with the input.
        low: the smallest simulated value, the low end of the inputs the
            evidence covers.
        high: the largest simulated value.
        n_equivalent: the checks which are equivalent.
        n_in_range: the checks in the reported range.
        n_differs: the checks which differ.
        n_not_reported: the checks without a reported value.
    """

    target: str
    kind: str
    unit: str
    n_mappings: int
    low: float
    high: float
    n_equivalent: int
    n_in_range: int
    n_differs: int
    n_not_reported: int


class InputsRecord(TypedDict):
    """The equivalency of the inputs of the decisive mappings of a context.

    Attributes:
        n_mappings: the decisive mappings.
        n_checked: the decisive mappings with at least one input checked.
        inputs: the summary of every input, sorted by kind and target.
        differs: the checks which differ.
    """

    n_mappings: int
    n_checked: int
    inputs: list[InputSummary]
    differs: list[InputCheck]


@dataclass(frozen=True)
class Reported:
    """A reported value of a group, an individual or an intervention.

    Attributes:
        by: the group, individual or intervention, e.g. `group all`.
        central: the central value, `None` without.
        low: the low end of the range, `None` without.
        high: the high end of the range, `None` without.
    """

    by: str
    central: float | None
    low: float | None = None
    high: float | None = None


@dataclass(frozen=True)
class StudyTables:
    """The tables of the subjects and the interventions of a study.

    Attributes:
        groups: the groups (`TabGroups`), empty without.
        individuals: the individuals (`TabIndividuals`), empty without.
        interventions: the interventions (`TabInterventions`), empty without.
    """

    groups: pd.DataFrame
    individuals: pd.DataFrame
    interventions: pd.DataFrame


@cache
def dose_prefixes() -> tuple[str, ...]:
    """Get the prefixes of the targets which are a dose.

    The dose states of the models, `pkpdsim.credibility.battery.dose_states`:
    the ones the PBPK library creates and the ones of the models which are
    not built with it. The infusion rate is none of them, it is a rate in
    mg/min and no dose. The battery imports the library inside the function,
    so that `pkpdsim.credibility` imports without it, and so does this.

    Returns:
        The prefixes, the library's first.
    """
    from pkpdsim.credibility.battery import dose_states

    return dose_states()


def _read(folder: Path, study: str, name: str) -> pd.DataFrame:
    """Get a table of a study, empty without the file."""
    path = folder / f".{study}_{name}.tsv"
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=True)


def read_study_tables(folder: Path | None, study: str) -> StudyTables:
    """Read the tables of the subjects and the interventions of a study.

    Args:
        folder: the folder of the study, `None` for a study without one.
        study: the name of the study.

    Returns:
        The tables, empty for a missing file.
    """
    if folder is None:
        empty = pd.DataFrame()
        return StudyTables(groups=empty, individuals=empty, interventions=empty)
    return StudyTables(
        groups=_read(folder, study, "TabGroups"),
        individuals=_read(folder, study, "TabIndividuals"),
        interventions=_read(folder, study, "TabInterventions"),
    )


def _number(value: object) -> float | None:
    """Get a finite number of a cell, `None` for an empty or other one."""
    try:
        number = float(value)  # ty: ignore[invalid-argument-type]
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _convert(
    value: float | None,
    unit: str,
    to: str,
    ureg: UnitRegistry,
    body_weight: Quantity | None = None,
) -> float | None:
    """Convert a reported value into the unit of the model.

    A value per body weight is multiplied by the body weight, a quantity.

    Raises:
        ValueError: if the unit cannot be converted.
    """
    if value is None:
        return None
    try:
        quantity = ureg.Quantity(value, unit)
        try:
            return float(quantity.to(to).magnitude)
        except DimensionalityError:
            if body_weight is None:
                raise
            return float((quantity * body_weight).to(to).magnitude)
    except (DimensionalityError, UndefinedUnitError, ValueError, TypeError) as err:
        raise ValueError(f"'{unit}' cannot be converted into '{to}': {err}") from err


def _group_values(row: pd.Series) -> tuple[float | None, float | None, float | None]:
    """Get the central value and the range of a row of a group."""
    mean, median, gmean = (_number(row.get(c)) for c in ("mean", "median", "gmean"))
    central = next((v for v in (mean, median, gmean) if v is not None), None)
    low, high = _number(row.get("min")), _number(row.get("max"))
    if low is not None and high is not None:
        return central, low, high
    sd = _number(row.get("sd"))
    se, count = _number(row.get("se")), _number(row.get("count"))
    if sd is None and se is not None and count:
        sd = se * math.sqrt(count)
    if mean is not None and sd is not None:
        return central, max(mean - SPREAD * sd, 0.0), mean + SPREAD * sd
    gsd = _number(row.get("gsd"))
    if gmean is not None and gsd is not None and gsd > 0:
        return central, gmean / gsd**SPREAD, gmean * gsd**SPREAD
    return central, None, None


def _rows(table: pd.DataFrame, name: str, measurement: str) -> pd.DataFrame:
    """Get the rows of a subject and a measurement type."""
    if table.empty or not {"name", "measurement_type"} <= set(table.columns):
        return table.iloc[0:0]
    return table[
        (table["name"].astype(str) == name)
        & (table["measurement_type"].astype(str) == measurement)
    ]


def _parent(groups: pd.DataFrame, name: str) -> str | None:
    """Get the parent of a group, `None` for a group without."""
    if groups.empty or "parent" not in groups.columns:
        return None
    parents = groups.loc[groups["name"].astype(str) == name, "parent"].dropna()
    return str(parents.iloc[0]) if len(parents) else None


def _group_reported(
    tables: StudyTables, group: str, measurement: str, unit: str, ureg: UnitRegistry
) -> Reported | None:
    """Get the reported value of a group or of the nearest of its parents."""
    seen: set[str] = set()
    name: str | None = group
    while name is not None and name not in seen:
        seen.add(name)
        rows = _rows(tables.groups, name, measurement)
        for _, row in rows.iterrows():
            central, low, high = _group_values(row)
            if central is None and low is None:
                continue
            row_unit = str(row.get("unit"))
            return Reported(
                by=f"group {name}",
                central=_convert(central, row_unit, unit, ureg),
                low=_convert(low, row_unit, unit, ureg),
                high=_convert(high, row_unit, unit, ureg),
            )
        name = _parent(tables.groups, name)
    return None


def reported_subjects(
    tables: StudyTables,
    groups: Iterable[str],
    individuals: Iterable[str],
    measurement: str,
    unit: str,
    ureg: UnitRegistry,
) -> list[Reported]:
    """Get the reported values of a covariate of the subjects of data.

    An individual with the measurement gives its value; one without it the
    value of its group. A group without the measurement gives the one of the
    nearest parent which has it.

    Args:
        tables: the tables of the study.
        groups: the groups of the rows of the data.
        individuals: the individuals of the rows of the data.
        measurement: the measurement type, e.g. `weight`.
        unit: the unit of the target in the model.
        ureg: the unit registry.

    Returns:
        The reported values, each subject once.

    Raises:
        ValueError: if a reported unit cannot be converted into `unit`.
    """
    found: dict[str, Reported] = {}
    for individual in individuals:
        rows = _rows(tables.individuals, individual, measurement)
        values = [
            (_number(row.get("value")), str(row.get("unit")))
            for _, row in rows.iterrows()
        ]
        values = [(v, u) for v, u in values if v is not None]
        if values:
            value, row_unit = values[0]
            found[f"individual {individual}"] = Reported(
                by=f"individual {individual}",
                central=_convert(value, row_unit, unit, ureg),
            )
            continue
        own = (
            tables.individuals.loc[
                tables.individuals["name"].astype(str) == individual, "group"
            ].dropna()
            if not tables.individuals.empty
            and {"name", "group"} <= set(tables.individuals.columns)
            else pd.Series(dtype=object)
        )
        if len(own):
            reported = _group_reported(
                tables, str(own.iloc[0]), measurement, unit, ureg
            )
            if reported is not None:
                found.setdefault(reported.by, reported)
    for group in groups:
        reported = _group_reported(tables, group, measurement, unit, ureg)
        if reported is not None:
            found.setdefault(reported.by, reported)
    return list(found.values())


def reported_doses(
    tables: StudyTables,
    interventions: Iterable[str],
    unit: str,
    ureg: UnitRegistry,
    body_weight: Quantity | None = None,
) -> list[Reported]:
    """Get the reported doses of the interventions of data.

    Args:
        tables: the tables of the study.
        interventions: the interventions of the rows of the data.
        unit: the unit of the dose target in the model.
        ureg: the unit registry.
        body_weight: the simulated body weight, a quantity, for a dose per
            body weight; such a dose is left out without.

    Returns:
        The dose of every intervention with one, which `unit` reads.
    """
    table = tables.interventions
    if table.empty or not {"name", "measurement_type", "value", "unit"} <= set(
        table.columns
    ):
        return []
    doses: list[Reported] = []
    for name in dict.fromkeys(interventions):
        rows = table[
            (table["name"].astype(str) == name)
            & (table["measurement_type"].astype(str) == "dosing")
        ]
        for _, row in rows.iterrows():
            value = _number(row["value"])
            try:
                dose = _convert(value, str(row["unit"]), unit, ureg, body_weight)
            except ValueError as err:
                logger.debug("The dose of '%s' is left out: %s", name, err)
                continue
            if dose is not None:
                doses.append(Reported(by=f"intervention {name}", central=dose))
    return doses


def _ratio(simulated: float, reported: float | None) -> float | None:
    """Get the ratio of the simulated to the reported value."""
    if reported is None or reported == 0:
        return None
    return simulated / reported


def _distance(simulated: float, reported: Reported) -> float:
    """Get the distance of a simulated value from a reported one in log."""
    ratio = _ratio(simulated, reported.central)
    if ratio is None or ratio <= 0:
        return math.inf
    return abs(math.log(ratio))


def check_covariate(
    simulated: float, reported: Sequence[Reported]
) -> tuple[InputStatus, Reported | None]:
    """Check a simulated covariate against the reported values of the subjects.

    Args:
        simulated: the simulated value.
        reported: the reported values of the subjects of the data.

    Returns:
        The status and the reported value it was decided by: the first
        equivalent one, else the first range which holds it, else the nearest
        central value; `None` without a reported value.
    """
    comparable = [r for r in reported if r.central is not None or r.low is not None]
    if not comparable:
        return InputStatus.NOT_REPORTED, None
    for r in comparable:
        ratio = _ratio(simulated, r.central)
        if ratio is not None and abs(ratio - 1.0) <= COVARIATE_TOLERANCE:
            return InputStatus.EQUIVALENT, r
    for r in comparable:
        if r.low is not None and r.high is not None and r.low <= simulated <= r.high:
            return InputStatus.IN_RANGE, r
    return InputStatus.DIFFERS, min(comparable, key=lambda r: _distance(simulated, r))


def check_dose(
    simulated: float, doses: Sequence[Reported]
) -> tuple[InputStatus, Reported | None]:
    """Check a simulated dose against the reported doses of the interventions.

    Args:
        simulated: the simulated dose.
        doses: the reported doses of the interventions of the data.

    Returns:
        The status and the nearest reported dose, `None` without one.
    """
    if not doses:
        return InputStatus.NOT_REPORTED, None
    nearest = min(doses, key=lambda d: _distance(simulated, d))
    ratio = _ratio(simulated, nearest.central)
    if ratio is not None and abs(ratio - 1.0) <= DOSE_TOLERANCE:
        return InputStatus.EQUIVALENT, nearest
    return InputStatus.DIFFERS, nearest


def _nominal(value: object) -> object:
    """Get the nominal value of an uncertain value, else the value."""
    return getattr(value, "nominal", value)


def _magnitude(value: object, unit: str) -> float | None:
    """Get a value of a change in the unit of its target, `None` for a formula.

    A number is in the unit of the target; a dose added to its state (`<state>
    + <number>`) is the number.
    """
    value = _nominal(value)
    if isinstance(value, str):
        match = _ADDED_DOSE.match(value)
        return _number(match["value"]) if match else None
    if isinstance(value, pint.Quantity):
        try:
            return float(value.to(unit).magnitude)
        except (DimensionalityError, UndefinedUnitError):
            return None
    return _number(value)


def _mapping_rows(fit_mapping: FitMapping) -> pd.DataFrame:
    """Get the rows of the dataset of the reference of a mapping."""
    reference = fit_mapping.reference
    experiment = reference.experiment
    datasets = getattr(experiment, "_datasets", None) or experiment.datasets()
    dataset = datasets.get(reference.dset_id) if reference.dset_id else None
    if dataset is None:
        return pd.DataFrame()
    rows = pd.DataFrame(dataset)
    for column, value in (reference.sel or {}).items():
        if column not in rows.columns:
            return rows.iloc[0:0]
        values = value if isinstance(value, list | tuple | set) else [value]
        rows = rows[rows[column].isin(values)]
    return rows


def _names(rows: pd.DataFrame, column: str, split: bool = False) -> list[str]:
    """Get the names of a column of the rows, each once."""
    if column not in rows.columns:
        return []
    names: list[str] = []
    for value in rows[column].dropna():
        text = str(value).strip()
        if not text or text == "NA":
            continue
        parts = text.split(",") if split else [text]
        names.extend(p.strip() for p in parts if p.strip())
    return list(dict.fromkeys(names))


def _check(
    key: tuple[str, str, str],
    target: str,
    kind: str,
    measurement: str,
    unit: str,
    simulated: float,
    source: str,
    status: InputStatus,
    reported: Reported | None,
    ureg: UnitRegistry,
) -> InputCheck:
    """Get the check of an input, its unit as a symbol."""
    fit, experiment, mapping = key
    try:
        unit = f"{ureg.Unit(unit):~}"
    except (UndefinedUnitError, ValueError, TypeError):
        logger.debug("The unit '%s' has no symbol", unit)
    return InputCheck(
        fit=fit,
        experiment=experiment,
        mapping=mapping,
        target=target,
        kind=kind,
        measurement=measurement,
        unit=unit,
        simulated=simulated,
        source=source,
        reported=None if reported is None else reported.central,
        low=None if reported is None else reported.low,
        high=None if reported is None else reported.high,
        reported_by=[] if reported is None else [reported.by],
        status=status.value,
        ratio=None if reported is None else _ratio(simulated, reported.central),
    )


def _mapping_checks(
    key: tuple[str, str, str],
    problem: OptimizationProblem,
    k: int,
    rows: pd.DataFrame,
    tables: StudyTables,
    covariates: Mapping[str, str],
) -> list[InputCheck]:
    """Get the checks of the inputs of the mapping `k` of a problem."""
    model = problem.models[k]
    uinfo = model.uinfo
    ureg = uinfo.ureg
    defined = problem.defined_changes[k]
    groups = _names(rows, "group")
    individuals = _names(rows, "individual")
    checks: list[InputCheck] = []
    body_weight = None
    # a fitted parameter is no input of the evidence
    fitted = {parameter.pid for parameter in problem.parameters}
    for target, measurement in covariates.items():
        unit = uinfo.get(target)
        if unit is None or target in fitted:
            continue
        if target in defined:
            simulated, source = _magnitude(defined[target], unit), "simulation"
        else:
            try:
                simulated, source = model.loaded_value(target), "model"
            except (KeyError, AttributeError):
                continue
        if simulated is None:
            continue
        if measurement == "weight":
            body_weight = ureg.Quantity(simulated, unit)
        try:
            reported = reported_subjects(
                tables, groups, individuals, measurement, unit, ureg
            )
        except ValueError as err:
            logger.warning("%s/%s %s: %s", key[1], key[2], target, err)
            reported = []
        status, by = check_covariate(simulated, reported)
        checks.append(
            _check(
                key,
                target,
                "covariate",
                measurement,
                unit,
                simulated,
                source,
                status,
                by,
                ureg,
            )
        )
    interventions = _names(rows, "intervention", split=True)
    simulation = problem.simulations[k]
    values: list[tuple[str, Any]] = list(defined.items())
    for change in simulation.changes or ():
        values.extend(change.values.items())
    doses: dict[str, list[float]] = {}
    for target, value in values:
        if not target.startswith(dose_prefixes()):
            continue
        unit = uinfo.get(target)
        if unit is None:
            continue
        amount = _magnitude(value, unit)
        if amount is None or amount == 0:
            continue
        if amount not in doses.setdefault(target, []):
            doses[target].append(amount)
    for target, amounts in doses.items():
        unit = str(uinfo.get(target))
        reported_dose = reported_doses(tables, interventions, unit, ureg, body_weight)
        for amount in amounts:
            status, by = check_dose(amount, reported_dose)
            checks.append(
                _check(
                    key,
                    target,
                    "dose",
                    "dosing",
                    unit,
                    amount,
                    "simulation",
                    status,
                    by,
                    ureg,
                )
            )
    return checks


def input_checks(
    model: Project, fits: FitTable, covariates: Mapping[str, str]
) -> list[InputCheck]:
    """Check the simulated inputs of every mapping of the adopted fits.

    Args:
        model: the model, the folders of its studies hold the reported values.
        fits: the evaluated adopted fits with their reports, see
            `sbmlsim.credibility.records.evaluate_fits`.
        covariates: the target of every covariate and the measurement type of
            its reported value, see `Credibility.covariates`.

    Returns:
        The checks of every mapping of the points, in the order of the fits and
        their mappings.
    """
    folders = study_folders(model)
    tables: dict[str, StudyTables] = {}
    checks: list[InputCheck] = []
    keys = fits.table[["fit", "experiment", "mapping"]].drop_duplicates()
    for fit, experiment, mapping in keys.itertuples(index=False):
        report = fits.reports.get(str(fit))
        if report is None:
            continue
        problem = report.problem
        indices = [
            k
            for k, (e, m) in enumerate(
                zip(problem.experiment_keys, problem.mapping_keys, strict=True)
            )
            if e == experiment and m == mapping
        ]
        if not indices:
            continue
        study = str(experiment)
        if study not in tables:
            folder = folders.get(study, [None])[0]
            tables[study] = read_study_tables(folder, study)
        rows = _mapping_rows(fits.lookup(study, str(mapping)))
        checks.extend(
            _mapping_checks(
                (str(fit), study, str(mapping)),
                problem,
                indices[0],
                rows,
                tables[study],
                covariates,
            )
        )
    return checks


def inputs_record(checks: Sequence[InputCheck], decisive: pd.DataFrame) -> InputsRecord:
    """Get the equivalency of the inputs of the decisive mappings of a context.

    Args:
        checks: the checks of the mappings of the model, see `input_checks`.
        decisive: the decisive points of the context, with the columns `fit`,
            `experiment` and `mapping`.

    Returns:
        The summary of the checks of the decisive mappings.
    """
    keys = (
        set(
            decisive[["fit", "experiment", "mapping"]]
            .astype(str)
            .itertuples(index=False, name=None)
        )
        if len(decisive)
        else set()
    )
    own = [c for c in checks if (c["fit"], c["experiment"], c["mapping"]) in keys]
    by_input: dict[tuple[str, str], list[InputCheck]] = {}
    for check in own:
        by_input.setdefault((check["kind"], check["target"]), []).append(check)
    summaries: list[InputSummary] = []
    for (kind, target), group in sorted(by_input.items()):
        values = [c["simulated"] for c in group]
        status = [c["status"] for c in group]
        summaries.append(
            InputSummary(
                target=target,
                kind=kind,
                unit=group[0]["unit"],
                n_mappings=len(
                    {(c["fit"], c["experiment"], c["mapping"]) for c in group}
                ),
                low=min(values),
                high=max(values),
                n_equivalent=status.count(InputStatus.EQUIVALENT.value),
                n_in_range=status.count(InputStatus.IN_RANGE.value),
                n_differs=status.count(InputStatus.DIFFERS.value),
                n_not_reported=status.count(InputStatus.NOT_REPORTED.value),
            )
        )
    return InputsRecord(
        n_mappings=len(keys),
        n_checked=len({(c["fit"], c["experiment"], c["mapping"]) for c in own}),
        inputs=summaries,
        differs=[c for c in own if c["status"] == InputStatus.DIFFERS.value],
    )


#: the checks which differ a statement names
NAMED_DIFFERENCES = 5


def _value(value: float | None, unit: str) -> str:
    """Get a value with its unit in words."""
    if value is None:
        return "n/a"
    return f"{value:.3g} {unit}".strip()


def _range(summary: InputSummary) -> str:
    """Get the simulated range of an input in words."""
    if summary["low"] == summary["high"]:
        return _value(summary["low"], summary["unit"])
    return f"{_value(summary['low'], '')} to {_value(summary['high'], summary['unit'])}"


def _counts(summary: InputSummary) -> str:
    """Get the statuses of the checks of an input in words."""
    parts = [
        f"{summary['n_equivalent']} equivalent",
        f"{summary['n_in_range']} in the reported range",
        f"{summary['n_differs']} differ" if summary["n_differs"] != 1 else "1 differs",
        f"{summary['n_not_reported']} not reported",
    ]
    return ", ".join(p for p in parts if not p.startswith("0 "))


def inputs_statement(record: InputsRecord) -> str:
    """Get the statement of the factor equivalency of inputs.

    Args:
        record: the equivalency of the inputs of a context.

    Returns:
        The statement: the checked mappings, per input the statuses and the
        values the evidence covers, and the first checks which differ.
    """
    inputs = "; ".join(
        f"{s['target']} ({_counts(s)}; simulated {_range(s)})" for s in record["inputs"]
    )
    parts = [
        f"The simulated inputs of {record['n_checked']} of the "
        f"{record['n_mappings']} decisive mappings against the subjects and the "
        f"doses their studies report: {inputs}."
    ]
    differs = record["differs"]
    if differs:
        named = "; ".join(
            f"{c['experiment']}/{c['mapping']} {c['target']} "
            f"{_value(c['simulated'], c['unit'])} against "
            f"{_value(c['reported'], c['unit'])} ({', '.join(c['reported_by'])})"
            for c in differs[:NAMED_DIFFERENCES]
        )
        more = len(differs) - NAMED_DIFFERENCES
        parts.append(f"Differ: {named}" + (f"; and {more} more." if more > 0 else "."))
    return " ".join(parts)


#: the columns of `input_ranges.tsv`, see `InputSummary`
INPUT_RANGE_COLUMNS: Final = (
    "context",
    "target",
    "kind",
    "unit",
    "n_mappings",
    "low",
    "high",
    "n_equivalent",
    "n_in_range",
    "n_differs",
    "n_not_reported",
)


def _plain(value: object) -> Any:  # noqa: ANN401  # plain data of JSON
    """Get a record as plain data of JSON, a copy."""
    return json.loads(json.dumps(value))


def section_record(section: Mapping[str, Any] | None) -> InputsRecord | None:
    """Get the equivalency of the inputs of the section `inputs` of the record of a context.

    Args:
        section: the section, `{"inputs": {...}}`, also of a record of pkpdsim
            0.19; `None` without.

    Returns:
        The record, `None` without.
    """
    data = None if section is None else section.get("inputs")
    if data is None:
        return None
    return InputsRecord(
        n_mappings=int(data["n_mappings"]),
        n_checked=int(data["n_checked"]),
        inputs=[InputSummary(**i) for i in data["inputs"]],
        differs=[InputCheck(**c) for c in data["differs"]],
    )


def section_checks(section: Mapping[str, Any] | None) -> list[InputCheck]:
    """Get the checks of the section `inputs` of the record of an assessment.

    Args:
        section: the section, `{"inputs": [...]}`, also of a record of pkpdsim
            0.19; `None` without.

    Returns:
        The checks, empty without.
    """
    if section is None:
        return []
    return [InputCheck(**check) for check in section.get("inputs", ())]


class InputsExtension(EvidenceExtension):
    """The equivalency of the inputs of the contexts of use, see the module.

    Prepared: the checks of the inputs of every mapping of the adopted fits.
    """

    name = "inputs"
    scorecard_names = ("inputs_differ",)

    @override
    def prepare(
        self, project: Project, credibility: Credibility, fits: FitTable
    ) -> object:
        """Check the inputs of the mappings with the covariates of the declaration."""
        covariates = (
            DEFAULT_COVARIATES
            if credibility.covariates is None
            else credibility.covariates
        )
        return input_checks(project, fits, covariates)

    @override
    def assess_context(
        self,
        cou: ContextOfUse,
        decisive: pd.DataFrame,
        fits: FitTable,
        prepared: object,
    ) -> ContextContribution | None:
        """Summarize the checks of the decisive mappings, none without the checks."""
        if not isinstance(prepared, list):
            return None
        record = inputs_record(prepared, decisive)
        return ContextContribution(
            section={"inputs": _plain(record)},
            factors=(
                {
                    "equivalency_of_inputs": FactorStatement(
                        inputs_statement(record), before_core=True
                    )
                }
                if record["n_checked"]
                else {}
            ),
        )

    @override
    def assessment_section(self, prepared: object) -> Mapping[str, Any] | None:
        """Get the checks of every mapping as the section of the assessment."""
        if not isinstance(prepared, list):
            return None
        return {"inputs": _plain(prepared)}

    @override
    def context_report_blocks(
        self, context: Mapping[str, Any], section: Mapping[str, Any] | None
    ) -> Sequence[ContextBlock]:
        """Show the inputs of the decisive mappings after the parameter uncertainty."""
        record = section_record(section)
        if record is None or not record["n_checked"]:
            return ()
        return (ContextBlock(render_macro("inputs", record), after="uncertainty"),)

    @override
    def report_tables(self, record: Mapping[str, Any]) -> Mapping[str, pd.DataFrame]:
        """Write `inputs.tsv` and `input_ranges.tsv`.

        `inputs.tsv`: one row per mapping of the adopted fits and simulated
        input, the reporting groups, individuals and interventions joined by
        `;`; `input_ranges.tsv`: one row per context of use and input of its
        decisive mappings.
        """
        checks = section_checks(record.get("extensions", {}).get(self.name))
        rows = [
            {**check, "reported_by": ";".join(check["reported_by"])} for check in checks
        ]
        ranges = [
            {"context": c["id"], **summary}
            for c in record["contexts"]
            if (r := section_record(c.get("extensions", {}).get(self.name))) is not None
            for summary in r["inputs"]
        ]
        return {
            "inputs": pd.DataFrame(rows, columns=list(INPUT_COLUMNS)),
            "input_ranges": pd.DataFrame(ranges, columns=list(INPUT_RANGE_COLUMNS)),
        }

    @override
    def scorecard_columns(
        self, context: Mapping[str, Any] | None
    ) -> Mapping[str, object]:
        """Get the inputs of the decisive mappings which differ, empty without the checks."""
        record = section_record(context)
        return {"inputs_differ": None if record is None else len(record["differs"])}


#: the evidence extension of the equivalency of the inputs, registered by
#: `pkpdsim.model.PKPDModel`
INPUTS: Final = InputsExtension()
